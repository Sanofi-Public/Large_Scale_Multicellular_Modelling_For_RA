"""
This script calculates the similarity score for each stable state generated with BMA
It takes as inputs the csv files corresponding to the stable states of the model and 
the Boolean values associated to some nodes in the model that were observed in the literature or in the RNA-seq dataset 
"""


import pandas as pd
import os
import glob

#define the function that calculates the SMC score
def SMC(x,y):
    shared=0
    for i in range(len(x)):
        if x[i] == y[i]:
            shared=shared+1
    return (shared/len(x))

#read the csv file containing the observed Boolean values
my_path = os.path.join("a-computational-framework-to-build-and-calibrate-large-scale-boolean-models-main","Attractors_filtering_based_on_similarity_scores")
os.chdir(my_path)
inputs_with_expected_values= pd.read_csv("M1_nodes_after_upstream_export_with_known_value_literature.csv",index_col=None,header=0,sep=";")
inputs_with_expected_values=inputs_with_expected_values.sort_values(by=['node'])
inputs_with_expected_values=inputs_with_expected_values.reset_index(drop=True)
inputs_val=inputs_with_expected_values["value"]

#create en empty list to store the similarity scores
i=0
score_list=[]
combi_index=[] 
my_path_2 = os.path.join("a-computational-framework-to-build-and-calibrate-large-scale-boolean-models-main","Attractors_search","BioCheckConsoleMulti","BioCheckConsoleMulti","bin","Debug","netcoreapp3.1")
os.chdir(my_path_2)
maximum=glob.glob('*.csv')

#for each stable state, the SMC score is calculated based on the nodes having a known Boolean value
while i < maximum :
    file_exists = os.path.exists("stable_nodes_"+str(i)+".csv")
    if file_exists:
        df_1 = pd.read_csv("stable_nodes_"+str(i)+".csv", index_col=None, header=0)
        df_1=df_1[df_1["stable_nodes_name"].isin(list(inputs_with_expected_values["node"]))]
        df_1=df_1.sort_values(by=['stable_nodes_name'])
        df_1=df_1.reset_index(drop=True)
        score=SMC(df_1["stable_nodes_val"],inputs_val) 
        if score > 0.90:
            score_list.append(score)
            combi_index.append(i)
        
    i=i+1   
os.chdir(my_path)   
score_similarity_on_inputs={'combination_index':combi_index,'similarity_score':score_list}
score_similarity_on_inputs=pd.DataFrame(score_similarity_on_inputs)
score_similarity_on_inputs.to_csv('M1_similarity_score_on_nodes.csv')
