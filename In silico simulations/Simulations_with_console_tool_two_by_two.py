# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 13:52:29 2023

@author: I0471594
"""



import pandas as pd
import json
import subprocess
import itertools
from itertools import combinations
from joblib import Parallel, delayed
import os
my_path = os.path.join("Large-scale multicellular modeling of the arthritic joint","In_silico_simulations","BMA")
os.chdir(my_path)


#file=open("multicellular_model.json.json")
#global model
#model = json.load(file)
#file.close()
global model
with open(''multicellular_model.json", 'r') as file:
    model = json.load(file)

df = pd.read_csv('multicell_targets_nodes.csv',sep=';', header=0)
ko=df['Target_nodes']
ko_sep=ko.str.split(',', expand=True)
all_ko=ko_sep.values.tolist()
all_ko=list(itertools.chain.from_iterable(all_ko))
all_ko=list(filter(None,all_ko))


ko_id=[]
ko_name=[]
for i in model['Model']["Variables"]:
    if i['Name'] in all_ko:
        ko_id.append(i['Id'])
        ko_name.append(i["Name"])
        
ko_id_name= {'name':ko_name,'id':ko_id}
ko_id_name=pd.DataFrame(ko_id_name)
ko_id_name=ko_id_name.drop_duplicates(subset=['id'])

for col in ko_sep.columns:
    for i in ko_sep.index:
        if str(ko_sep[col][i]) in set(ko_id_name["name"]):
            n=ko_sep[col][i]
            index=ko_id_name[ko_id_name["name"]==n].index.values
            new_val=ko_id_name['id'][index]
            ko_sep.at[i,col] = new_val.values[0]
            
M1_targets=df[['Target_nodes']]
global ko_sep_with_targets
ko_sep_with_targets=pd.concat([M1_targets,ko_sep],axis=1)
ko_sep_with_targets=pd.concat([df[['Names']],ko_sep_with_targets],axis=1)


combi=[",".join(map(str, comb)) for comb in combinations(ko_sep_with_targets['Names'],2)]
combi = [(i, combi[i]) for i in range(len(combi))]


def drug_combination(item):
    co,pair=item
    args='BioCheckConsole.exe -engine VMCAI -model multicellular_model.json -prove stability_analysis_'+str(co)+'.json'
    pair=pair.split(',')
    data=ko_sep_with_targets[ko_sep_with_targets['Names'].isin(pair)]
    data=data.drop(columns=['Names', 'Target_nodes'])
    ko_list=[]
    for col in data.columns:
       ko_list.append(data[col].values.tolist())
    ko_list=list(itertools.chain.from_iterable(ko_list))
    ko_list=[x for x in ko_list if str(x) != 'nan']
    ko_list=[x for x in ko_list if str(x) != ' ']
    ko_list=list(filter(None,ko_list))
    ko_list=[int(item) for item in ko_list]
    
    for ko in ko_list:
      args=args +' '+ '-ko '+ str(ko) +' '+ str(0)

    subprocess.run(args,shell=True)



    stable_nodes_id=[]
    stable_nodes_val=[]
    stable_nodes_name=[]
      
    
    #file=open('stability_analysis_'+str(co)+'.json')
    #result=json.load(file)
    #file.close()
    with open('stability_analysis_'+str(co)+'.json', 'r') as file:
    result = json.load(file)

      
    for i in result['Ticks'][0]['Variables']:
        if i["Lo"]==i["Hi"]:
            stable_nodes_id.append(i['Id'])
            stable_nodes_val.append(i['Lo'])
    
    
    for i in stable_nodes_id:
        for k in model['Model']["Variables"]:
            if i == k['Id']:
                stable_nodes_name.append(k['Name'])
                
    
   
        
    df_stable={'stable_nodes_id':stable_nodes_id,'stable_nodes_val':stable_nodes_val,'stable_nodes_name':stable_nodes_name}
    df_stable=pd.DataFrame(df_stable)
    df_stable.to_csv('Large-scale multicellular modeling of the arthritic joint\\In_silico_simulations\\Two_by_Two_simulation_results_multicell\\Simulation'+'_'+str(co)+'.csv')
    os.remove('stability_analysis_'+str(co)+'.json')

Parallel(n_jobs=7)(delayed(drug_combination)(item) for item in combi)

  
combi={'drug_pairs':combi}
combi=pd.DataFrame(combi)
combi.to_csv('Large-scale multicellular modeling of the arthritic joint\\In_silico_simulations\\Two_by_Two_simulation_results_multicell\\two_by_two_combination_KO_multicell.csv')
  
