import pandas as pd
import numpy as np

def transform_fragment_name(name):
    # Transform fragment name to match the sequencing format
    parts = name.split("_")
    position = int(int(parts[2])/3) + 1
    return f"{parts[0]}_{parts[1]}:{position}"



def sort_names(input_list):
    # Group items by their base name (everything before the first ":")
    groups = {}
    for fragment in input_list:
        base_name = fragment.split(':')[0]
        if base_name not in groups:
            groups[base_name] = []
        groups[base_name].append(fragment)
    
    # Sort items within each group based on the numerical values
    for base_name in groups:
        groups[base_name].sort(key=lambda x: [int(n) for n in x.split(':')[1:]])
    
    # Maintain original order of base names and combine sorted groups
    sorted_list = []
    seen_base_names = []
    for fragment in input_list:
        base_name = fragment.split(':')[0]
        if base_name not in seen_base_names:
            seen_base_names.append(base_name)
            sorted_list.extend(groups[base_name])
    
    return sorted_list


#HPylori genome-preserving-order synthesized fragments
fragments = pd.read_excel("../FragmentList.xlsx", sheet_name="Sheet")
fragments["Name"] = fragments["Name"].apply(lambda x: transform_fragment_name(x))
sorted_fragments = sort_names(fragments["Name"].tolist())
print(sorted_fragments)

output_file = open("sorted_genomeFrags.txt", "w")
for fragment in sorted_fragments:
	output_file.write(fragment + "\n")



