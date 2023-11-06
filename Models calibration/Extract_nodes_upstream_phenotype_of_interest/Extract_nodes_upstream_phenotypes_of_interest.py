# -*- coding: utf-8 -*-
"""
This script extracts the nodes present upstream of the phenotype of interest from the RA macrophage models in a JSON format

The inputs files are:
-M1_macrophage_all_nodes_in_json.csv that contains all the nodes present in the entire model 
-M1_model_export_upstream_phenotype_of_interest.sbml : the submodule containing the nodes upstream of the phenotypes of interest in SBML format.
This file was obtained using the export option of CaSQ tool

The output file is:
-M1_macrophage_nodes_upstream_phenotypes_in_json.csv that contains the nodes upstream of the phenotypes of interest in the JSON format

"""

import pandas as pd
import biolqm
import os


my_path = os.path.join("a-computational-framework-to-build-and-calibrate-large-scale-boolean-models-main","Extract_nodes_upstream_phenotype_of_interest")
os.chdir(my_path)

#load sbml file containing the nodes upstream of the phenotypes of interest and generate the list of nodes upstream of the phenotypes of interest 
model_sbml = biolqm.load("M1_model_export_upstream_phenotype_of_interest.sbml")
model_sbml = biolqm.sanitize(model_sbml)
nodes_in_sbml=list(c.getNodeID() for c in model_sbml.getComponents())
dic = {"nodes_in_sbml": nodes_in_sbml}
nodes_in_sbml= pd.DataFrame(dic)

#process the nodes name in the SBML format by removing the underscores
sbml_nodes_without_underscore=list(s.replace('_', '') for s in nodes_in_sbml["nodes_in_sbml"])
sbml_df={'nodes_in_sbml':nodes_in_sbml["nodes_in_sbml"],'sbml_nodes_without_underscore': sbml_nodes_without_underscore}
sbml_df=pd.DataFrame(sbml_df)
sbml_df=sbml_df.sort_values('sbml_nodes_without_underscore')
sbml_df=sbml_df.reset_index(drop=True)


#read the list of nodes present in the entire model in JSON format
json_names = pd.read_csv("M1_macrophage_all_nodes_in_json.csv", index_col=None, header=0)
#process the nodes name by removing the underscores
json_names["json_nodes_without_underscore"]=list(s.replace('_', '') for s in json_names['Name'])
json_names=json_names.sort_values('json_nodes_without_underscore') 
json_names=json_names.reset_index(drop=True)

#filter out the nodes in the JSON file that are not upstream of the phenotypes of interest
json_names_filtered=json_names[json_names['json_nodes_without_underscore'].isin(list(sbml_df["sbml_nodes_without_underscore"]))]
json_names_filtered=json_names_filtered["Name"]
json_names_filtered.to_csv("M1_macrophage_nodes_upstream_phenotypes_in_json.csv")
