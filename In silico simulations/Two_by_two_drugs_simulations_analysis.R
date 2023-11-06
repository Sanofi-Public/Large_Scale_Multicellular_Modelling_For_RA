library(readxl)
library(xlsx)
library(sys)
library(plyr)
library(stringi)
setwd("large-scale-multicellular-modeling-of-the-arthritic-joint//In silico simulations")
calibrated <- read.csv("Multicell_model_calibrated_state_with_no_oscillations.csv", sep=";")
vec=calibrated$stable_nodes_name
drug_combi<- read.csv("two_by_two_combination_KO_multicell.csv")

files <- list.files(pattern = "Simulation*")
interesting=c()
for (i in files){
  modif=read.csv(i)
  modif=modif[modif$stable_nodes_name %in% calibrated$stable_nodes_name,]
  calibrated=calibrated[calibrated$stable_nodes_name %in% modif$stable_nodes_name,]
  modif <- modif[order(modif$stable_nodes_name),]
  calibrated=calibrated[order(calibrated$stable_nodes_name),]
  modif=cbind(modif,calibrated$stable_nodes_val)
  modif=modif[modif$stable_nodes_val!=modif$`calibrated$stable_nodes_val`,]
  if ("inflammation_signal_phenotype" %in% modif$stable_nodes_name){ 
    res=stri_replace_all_regex(i,pattern=c('Simulation_', '.csv'),replacement=c('', ''),vectorize=FALSE)
    interesting=append(strtoi(res),interesting)}
}

drug_combi_interesting=subset(drug_combi, X %in% interesting)
drug_combi_interesting =drug_combi_interesting[- grep("| ",drug_combi_interesting$drug_pairs),]

