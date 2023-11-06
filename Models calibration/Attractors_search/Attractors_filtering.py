"""
This script filters the attractors obtained via BMA to keep the stable states only
it takes as inputs the csv files containing the list of stable nodes for each inputs combination (stable_nodes.csv)
it filters out the csv files where not all the nodes are stable
"""
import os
import pandas as pd
import glob
my_path = os.path.join("a-computational-framework-to-build-and-calibrate-large-scale-boolean-models-main","Attractors_search","BioCheckConsoleMulti","BioCheckConsoleMulti","bin","Debug","netcoreapp3.1")
os.chdir(my_path)
co=0
maximum=glob.glob('*.csv')
while co < maximum :
    file_exists = os.path.exists("stable_nodes_"+str(co)+".csv")
    if file_exists:
        df = pd.read_csv("stable_nodes_"+str(co)+".csv", index_col=None, header=0)
        #check if all the nodes are stable
        if len(df)!= 309:
            os.remove("stable_nodes_"+str(co)+".csv")
    co=co+1