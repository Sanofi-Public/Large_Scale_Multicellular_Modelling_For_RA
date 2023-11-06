library(readxl)
library(xlsx)
library(readxl)
library(xlsx)
library(sys)
library(plyr)
library(stringi)
setwd("C:/Users/I0471594/PhD_2/In_silico_simulations_multicellular_model/One_by_One_simulation_results_multicell")
calibrated <- read.csv("multicell_calibrated_state.csv", sep=";")
vec=calibrated$stable_nodes_name
drug_combi<- read.csv("Multicell_drug_targets_in_model.csv",sep=";")
files <- list.files(pattern = "Simulation*")
interesting=c()
for (i in files){
  print(i)
  modif=read.csv(i)
  modif <- modif[order(modif$stable_nodes_name),]
  calibrated=calibrated[order(calibrated$stable_nodes_name),]
  modif=cbind(modif,calibrated$stable_nodes_val)
  modif=modif[modif$stable_nodes_val!=modif$`calibrated$stable_nodes_val`,]
  if ("angiogenesis_signal_phenotype" %in% modif$stable_nodes_name){ 
    res=stri_replace_all_regex(i,pattern=c('Simulation_', '.csv'),replacement=c('', ''),vectorize=FALSE)
    interesting=append(strtoi(res)+1,interesting)}
}

drug_interesting=as.data.frame(drug_combi[interesting,])
