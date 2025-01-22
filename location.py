from goatools import obo_parser
import requests
import pandas as pd
import wget
from collections import defaultdict
import time

def get_protein_cellular_locations(uniprot_id, go_obo):
    """
    Retrieve cellular locations for a protein using UniProt and GO annotations
    
    Args:
        uniprot_id (str): UniProt ID of the protein
        go_obo: GO database object
        
    Returns:
        list: List of cellular locations
    """
    try:
        # Get GO annotations from UniProt
        go_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}?format=json"
        response = requests.get(go_url)
        protein_data = response.json()
        
        cellular_locations = []
        
        # Extract GO terms related to cellular component
        for annotation in protein_data.get('uniProtKBCrossReferences', []):
            if annotation.get('database') == 'GO':
                go_id = annotation.get('id')
                if go_id in go_obo and go_obo[go_id].namespace == 'cellular_component':
                    cellular_locations.append(go_obo[go_id].name)
        
        return "; ".join(cellular_locations) if cellular_locations else "No location data found"
    
    except Exception as e:
        return f"Error retrieving location data: {str(e)}"

def process_mdm2_interactions(input_file, go_obo):
    """
    Process MDM2 interactions file and add cellular locations
    
    Args:
        input_file (str): Path to input file
        go_obo: GO database object
        
    Returns:
        pandas.DataFrame: Processed data with locations
    """
    # Read the interaction data
    df = pd.read_csv(input_file, sep='\t', header=None)
    
    # Assign column names
    df.columns = ['ID', 'UniprotID1', 'UniprotName1', 'GeneName1', 
                 'UniprotID2', 'UniprotName2', 'GeneName2',
                 'Score1', 'Score2', 'Score3', 'Score4', 'FinalScore']
    
    # Create a dictionary to store locations (for caching)
    location_cache = {}
    
    # Process each unique protein that interacts with MDM2
    unique_proteins = df[df['UniprotID2'] == 'Q00987']['UniprotID1'].unique()
    
    print(f"Processing {len(unique_proteins)} unique proteins...")
    
    # Process proteins in batches to avoid overwhelming the API
    for i, uniprot_id in enumerate(unique_proteins):
        if i % 10 == 0:
            print(f"Processing protein {i+1} of {len(unique_proteins)}")
            
        if uniprot_id not in location_cache:
            # Add small delay to avoid overwhelming the UniProt API
            time.sleep(0.5)
            location_cache[uniprot_id] = get_protein_cellular_locations(uniprot_id, go_obo)
    
    # Add locations to the dataframe
    df['Location'] = df.apply(lambda row: 
        location_cache.get(row['UniprotID1']) if row['UniprotID2'] == 'Q00987' 
        else location_cache.get(row['UniprotID2']) if row['UniprotID1'] == 'Q00987'
        else "Not MDM2 interaction", axis=1)
    
    return df

def main():
    # Download GO database if not already present
    try:
        go_obo = obo_parser.GODag("go-basic.obo")
    except:
        print("Downloading GO database...")
        wget.download("http://purl.obolibrary.org/obo/go/go-basic.obo")
        go_obo = obo_parser.GODag("go-basic.obo")
    
    # Process the interactions file
    input_file = "../apid_humanMDM2.txt"
    print("Processing MDM2 interactions file...")
    
    result_df = process_mdm2_interactions(input_file, go_obo)
    
    # Save the results
    output_file = "mdm2_interactions_with_locations.tsv"
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults saved to {output_file}")
    
    # Print some statistics
    print("\nLocation Statistics:")
    location_counts = result_df[result_df['Location'] != "Not MDM2 interaction"]['Location'].value_counts()
    print(f"Total interactions with locations: {len(location_counts)}")
    print("\nMost common locations:")
    print(location_counts.head())

if __name__ == "__main__":
    main()