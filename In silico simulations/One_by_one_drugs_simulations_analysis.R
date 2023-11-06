library(readxl)
library(xlsx)
library(readxl)
library(xlsx)
library(sys)
library(plyr)
library(stringi)
setwd("Large-scale multicellular modeling of the arthritic joint//In silico simulations")
calibrated <- read.csv("Multicell_model_calibrated_state_with_no_oscillations.csv", sep=";")
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
