"""
This script runs the parallelized attractors search for the macrophage models in JSON format.
This script runs inside the BMA Linux version built with dotnet
It takes as inputs : the model in a JSON format and a csv file containing the list of inputs with their observed Boolean values ()

"""
#!/usr/bin/python 
import os
from itertools import product
import subprocess
import json
import csv
import pandas as pd
from joblib import Parallel, delayed

my_path = os.path.join("Large-scale multicellular modeling of the arthritic joint","Attractors_search","BioCheckConsoleMulti","BioCheckConsoleMulti","bin","Debug","netcoreapp3.1")
os.chdir(my_path)

#read the model JSON file
file=open("RA_fibroblast.json.json")
global model
model = json.load(file)
file.close()


#read the list of all the inputs in the chosen cell-specific model
global inputs_name
inputs_name=[]
global fixed_inputs_val
fixed_inputs_val=[]
global combi_len

with open('fibroblast_all_inputs_with_fixed_ones.csv', newline = '') as file:                                                                                          
    file.reader = csv.reader(file, delimiter=';')
    for f in file.reader:
        inputs_name.append(f[0])
        fixed_inputs_val.append(f[1])

fixed_inputs_val=list(filter(None,fixed_inputs_val))

#the combination length includes the inputs with no associated values
combi_len=len(inputs_name)-len(fixed_inputs_val)

#retrieve the IDs of the corresponding inputs nodes
global inputs_id
inputs_id=[]
global inputs
inputs=[]
for i in model['Model']["Variables"]:
    if i['Name'] in inputs_name:
        inputs_id.append(i['Id'])
        inputs.append(i["Name"])

#generate all the possible combinations of zeros and ones 
l=list(product(range(2), repeat=combi_len))
l = [(i, l[i]) for i in range(len(l))]

#for each inputs combinations, attractors function calculate the corresponding attractor and filter them depending on the phenotypes values
def attractors(item):
    ind,combi=item
      
    args='./BioCheckConsoleMulti -engine VMCAI -model RA_fibroblast.json -prove stability_analysis_'+str(ind)+'.json'
    d= {'inputs_name':inputs,'inputs_id':inputs_id}
    d=pd.DataFrame(d)
    d['inputs_name_cat'] = pd.Categorical(d['inputs_name'],categories=list(inputs_name), ordered=True)
    d=d.sort_values('inputs_name_cat')       
      
    val=list(fixed_inputs_val)+list(combi)
    df= {'inputs_id':d['inputs_id'],'val':val,'inputs_name':d['inputs_name']}
    df=pd.DataFrame(df)
    df.to_csv('inputs_combination_'+str(ind)+'.csv')
    for index, row in df.iterrows():
          args=args +' '+ '-ko '+ str(row["inputs_id"]) +' '+ str(row['val'])
    subprocess.run(args,shell=True)
      
    #read the results
    stable_nodes_id=[]
    stable_nodes_val=[]
    stable_nodes_name=[]
      
    file_exists = os.path.exists('stability_analysis_'+str(ind)+'.json')
    if file_exists:
        file=open('stability_analysis_'+str(ind)+'.json')
        result=json.load(file)
        file.close()
          
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
        df_stable.to_csv('stable_nodes_'+str(ind)+'.csv')


        #phenotype are stables?
        proliferation=df_stable.index[df_stable['stable_nodes_name']=="proliferation_survival_Fibroblast_phenotype"].tolist()
        apoptosis=df_stable.index[df_stable['stable_nodes_name']=="apoptosis_Fibroblast_phenotype"].tolist()
        migration=df_stable.index[df_stable['stable_nodes_name']=="migration_Fibroblast_phenotype"].tolist()
        
        #phenotype stable and with the right value?
        if proliferation and apoptosis and osteo: 

                proliferation=df_stable["stable_nodes_val"].iloc[proliferation]
                proliferation=proliferation.reset_index(drop=True)
                proliferation=proliferation[[0]]==1
                proliferation=proliferation.bool()


                apoptosis=df_stable["stable_nodes_val"].iloc[apoptosis] 
                apoptosis=apoptosis.reset_index(drop=True)
                apoptosis=apoptosis[[0]]==0
                apoptosis=apoptosis.bool()
                
                
                migration=df_stable["stable_nodes_val"].iloc[migration] 
                migration=osteo.reset_index(drop=True)
                migration=migration[[0]]==1
                migration=migration.bool()
                
                if (proliferation and apoptosis and migration):

                    try:
                           os.remove('stability_analysis_'+str(ind)+'.json')
                    except Exception:
                           pass

                    try:
                           os.remove('inputs_combination_'+str(ind)+'.csv')
                    except Exception:
                           pass

                    file_exists = os.path.exists('stability_analysis_'+str(ind)+'_cex.json')
                    if file_exists:
                        try:
                               os.remove('stability_analysis_'+str(ind)+'_cex.json')
                        except Exception:
                               pass

                else:
                    try:
                        os.remove('stable_nodes_'+str(ind)+'.csv')
                    except Exception:
                        pass
                    try:
                        os.remove('stability_analysis_'+str(ind)+'.json')
                    except Exception:
                        pass
                    try:
                        os.remove('inputs_combination_'+str(ind)+'.csv')
                    except Exception:
                        pass
                    file_exists = os.path.exists('stability_analysis_'+str(ind)+'_cex.json')
                    if file_exists:
                        try:
                            os.remove('stability_analysis_'+str(ind)+'_cex.json')
                        except Exception:
                            pass


        else:
                    try:
                        os.remove('stable_nodes_'+str(ind)+'.csv')
                    except Exception:
                        pass
                    try:
                        os.remove('stability_analysis_'+str(ind)+'.json')
                    except Exception:
                        pass
                    try:
                        os.remove('inputs_combination_'+str(ind)+'.csv')
                    except Exception:
                        pass
                    file_exists = os.path.exists('stability_analysis_'+str(ind)+'_cex.json')
                    if file_exists:
                        try:
                            os.remove('stability_analysis_'+str(ind)+'_cex.json')
                        except Exception:
                            pass
                    
    else:
        try:
            os.remove('inputs_combination_'+str(ind)+'.csv')
        except Exception:
              pass
        file_exists = os.path.exists('stability_analysis_'+str(ind)+'_cex.json')
        if file_exists:
            try:
                os.remove('stability_analysis_'+str(ind)+'_cex.json')
            except Exception:
                  pass

Parallel(n_jobs=96)(delayed(attractors)(item) for item in l)
