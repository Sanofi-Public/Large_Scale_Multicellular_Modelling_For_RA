"""
This script calculates the mean vector of all the stable states that have the highest similarity score
it takes as input the CSV file containing all the SMC scores associated to all the stable states. It filters this file to keep the ones with the highest score
The average vector is then calculated over those stable states
"""
import pandas as pd
import os

my_path = os.path.join("a-computational-framework-to-build-and-calibrate-large-scale-boolean-models-main","Attractors_filtering_based_on_similarity_scores")
os.chdir(my_path)

#filter the fixpoints and keep the ones with max similarity:
fixpoints_max_similarity = pd.read_csv("M1_similarity_score_on_nodes.csv", index_col=None, header=0,sep=";")
fixpoints_max_similarity=fixpoints_max_similarity[fixpoints_max_similarity['similarity_score']>0.98]

#read the list of nodes that are upstream of the phenotypes of interest
json_names_filtered=pd.read_csv("M1_macrophage_nodes_upstream_phenotypes_in_json.csv", index_col=None, header=0)


#calculate the number of occurences
my_path_2 = os.path.join("a-computational-framework-to-build-and-calibrate-large-scale-boolean-models-main","Attractors_search","BioCheckConsoleMulti","BioCheckConsoleMulti","bin","Debug","netcoreapp3.1")
os.chdir(my_path_2)
store=[0] * len(json_names_filtered)
for i in fixpoints_max_similarity['combination_index']:
    file_exists = os.path.exists("stable_nodes_"+str(i)+".csv")
    if file_exists:
        data = pd.read_csv("stable_nodes_"+str(i)+".csv", index_col=None, header=0)
        data=data[data["stable_nodes_name"].isin(list(json_names_filtered))]
        data=data.sort_values(by=['stable_nodes_name'])
        store += data["stable_nodes_val"]
    
    

occurences={'stable_nodes_name':data["stable_nodes_name"],'occurence':store/len(fixpoints_max_similarity)}    
occurences=pd.DataFrame(occurences)
os.chdir(my_path)  
occurences.to_csv('M1_occurence_per_node_in_fixpoints_with_max_similarity_on_nodes.csv')
