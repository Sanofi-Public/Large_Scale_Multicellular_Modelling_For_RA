# -*- coding: utf-8 -*-
"""
This script extracts the nodes present upstream of the phenotype of interest from the models in a JSON format
"""

import pandas as pd
import biolqm
import os


my_path = os.path.join("Large-scale multicellular modeling of the arthritic joint","Extract_nodes_upstream_phenotype_of_interest")
os.chdir(my_path)

#load sbml file containing the nodes upstream of the phenotypes of interest and generate the list of nodes upstream of the phenotypes of interest 
model_sbml = biolqm.load("RA_fibroblast_export_uptream_phenotype.sbml")
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
json_names = pd.read_csv("fibroblast_all_nodes_in_json.csv", index_col=None, header=0)
#process the nodes name by removing the underscores
json_names["json_nodes_without_underscore"]=list(s.replace('_', '') for s in json_names['Name'])
json_names=json_names.sort_values('json_nodes_without_underscore') 
json_names=json_names.reset_index(drop=True)

#filter out the nodes in the JSON file that are not upstream of the phenotypes of interest
json_names_filtered=json_names[json_names['json_nodes_without_underscore'].isin(list(sbml_df["sbml_nodes_without_underscore"]))]
json_names_filtered=json_names_filtered["Name"]
json_names_filtered.to_csv("fibroblast_nodes_upstream_phenotypes_in_json.csv")
