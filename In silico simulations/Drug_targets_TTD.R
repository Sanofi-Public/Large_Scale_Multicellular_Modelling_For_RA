library(plyr)
library(stringr)
library(tidyr)
library(dplyr)
library(readxl)
library(fuzzyjoin)
TTD <- read_excel("large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\P1-07-Drug-TargetMapping.xlsx")

Filtered_TTD=subset(TTD,MOA != "."& MOA != "Agonist" & MOA != "Modulator"& MOA != "Activator"
                    & MOA != "Stimulator"& MOA != "Stabilizer"& MOA != "Agonis; Inverse agonist"
                    & MOA != "Binder"& MOA != "Replacement"& MOA != "Reactivator"&
                    MOA != "Regulator (upregulator)"& MOA != "Regulator"
                    & MOA != "Partial agonist"& MOA != "Opener"& MOA != "Modulator (allosteric modulator)"
                    & MOA != "Ligand"& MOA != "Inducer"& MOA != "CAR-T-Cell-Therapy"
                    & MOA != "CAR-T-Cell-Therapy(Dual specific)" & MOA != "Stablizer"
                    & MOA != "Chelator" & MOA != "Co-agonist" & MOA != "Cofactor" 
                    & MOA != "Enhancer" & MOA != "Immune response agent" 
                    & MOA != "Immunomodulator" & MOA != "Immunostimulant"
                    & MOA != "Immunomodulator (Immunostimulant)")


TTD_Target_GeneName <- read.delim("large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\P1-01-TTD_target_download.txt", header=FALSE, comment.char="#")
TTD_Target_GeneName_f=subset(TTD_Target_GeneName,V2=="GENENAME")
TTD_Target_GeneName_f=TTD_Target_GeneName_f[,c(1,3)]
colnames(TTD_Target_GeneName_f)=c("TargetID","GeneName")
TTD_Target_GeneName_f= TTD_Target_GeneName_f%>% 
  mutate(GeneName = strsplit(as.character(GeneName), "; ")) %>% unnest(GeneName)

Filtered_TTD_with_genenames=join(Filtered_TTD,TTD_Target_GeneName_f,type="inner", by="TargetID",match="all")
length(unique(Filtered_TTD_with_genenames$TargetID)) #1643 unique targets


Multicellular_model_occurence<- read.csv("large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\Multicell_model_calibrated_state_with_no_oscillations.csv", sep=";")
Multicellular_nodes_1=as.data.frame(Multicellular_model_occurence$stable_nodes_name)
Multicellular_nodes=as.data.frame(str_replace_all(Multicellular_nodes_1$`Multicellular_model_occurence$stable_nodes_name`, c("_complex" = "", "_rna" = "",
                                                                                                   "_M2_macrophage_nucleus"="",
                                                                                                   "_phosphorylated"="",
                                                                                                   "_signal_phenotype"="",
                                                                                                   "_fibroblast_phenotype"="",
                                                                                                   "_M1_macrophage_phenotype"="",
                                                                                                   "_M2_macrophage_phenotype"="",
                                                                                                   "_TH1_phenotype"="",
                                                                                                   "_Fibroblast___Extracellular_Space"="",
                                                                                                   "_M1_macrophage___nucleus"="",
                                                                                                   "_M2_macrophage___Secreted_components"="",
                                                                                                   "_Fibroblast__Cytoplasm"="",
                                                                                                   "_Fibroblast___Mitochondrion"="",
                                                                                                   "_M1_macrophage__Mitochondria_membrane"="",
                                                                                                   "_M2_macrophage_mitochondrion_membrane"="",
                                                                                                   "_active"="",
                                                                                                   "_Fibroblast___nucleus"="",
                                                                                                   "_M2_macrophage__secreted_components"="",
                                                                                                   "_simple_molecule"="",
                                                                                                   "_M1_macrophage__Cytoplasmic_membrane_up"="",
                                                                                                   "_M2_macrophage__Cytoplasmic_membrane_up"="",
                                                                                                   "_M1_macrophage___Extracellular_Space"="",
                                                                                                   "_M2_macrophage__extracellular_space"="",
                                                                                                   "TH1___cytoplasmic_membrane_up"="",
                                                                                                   "_Fibroblast___secreted_components"="",
                                                                                                   "_M1_macrophage___Secreted_components"="",
                                                                                                   "_TH1___extracellular_space"="",
                                                                                                   "_TH1___secreted_components"="",
                                                                                                   "_M1_macrophage__Cytoplasmic_membrane_down"="",
                                                                                                   "_M2_macrophage__cytoplasmic_membrane_down"="",
                                                                                                   "_TH1___cytoplasmic_membrane_down"="",
                                                                                                   "_Fibroblast__Cytoplasmic_membrane_down"="",
                                                                                                   "_TH1___cytoplasmic_membrane_up"="",
                                                                                                   "_TH1___cytoplasmic_membrane_up"="",
                                                                                                   "_M1_macrophage___Cytoplasm"="",
                                                                                                   "_TH1___cytoplasm"="",
                                                                                                   "_TH1___nucleus"="",
                                                                                                   "_Fibroblast___Cytoplasm"="",
                                                                                                   "_M2_macrophage__cytoplasm"="",
                                                                                                   "_M1_macrophage__Mitochondria"="",
                                                                                                   "_M2_macrophage__mitochondria"="",
                                                                                                   "ic_membrane_up"="",
                                                                                                   "ic_membrane_down"="",
                                                                                                   "_ubiquitinated"="",
                                                                                                   "_empty"="",
                                                                                                   "_INFLAMMASOME"="")))
                                                                                                   
colnames(Multicellular_nodes)="Name"
colnames(Multicellular_nodes_1)="Name"
Multicellular_nodes_1$Name<-as.character(Multicellular_nodes_1$Name)
Multicellular_nodes = separate_rows(Multicellular_nodes,1,sep = "_")
multicell_Gene_names_Node_names=regex_left_join(Multicellular_nodes_1,Multicellular_nodes,by = c("Name" = "Name"))
write.csv(multicell_Gene_names_Node_names,"large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\multicell_Gene_names_Node_names.csv")
multicell_Gene_names_Node_names <- read.csv("large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\multicell_Gene_names_Node_names.csv", sep=";")
multicell_Gene_names_Node_names=multicell_Gene_names_Node_names %>% group_by(Gene) %>% 
  mutate(Node = paste(Node, collapse=","))
multicell_Gene_names_Node_names=multicell_Gene_names_Node_names[!duplicated(multicell_Gene_names_Node_names$Gene),]
write.csv(multicell_Gene_names_Node_names,"large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\multicell_Gene_names_Node_names.csv")




multicell_Gene_names_Node_names <- read.csv("large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\multicell_Gene_names_Node_names.csv", sep=";")
multicell_drug_targets_in_model= subset(Filtered_TTD_with_genenames, GeneName %in% multicell_Gene_names_Node_names$Name)
multicell_drug_targets_in_model_name=as.data.frame(unique(multicell_drug_targets_in_model$GeneName))
colnames(multicell_drug_targets_in_model_name)="Name"

multicell_drug_targets_in_model_target_nodes=join(multicell_drug_targets_in_model_name,multicell_Gene_names_Node_names,type="left", by="Name",match="all")

multicell_drug_targets_in_model_target_nodes=multicell_drug_targets_in_model_target_nodes[!duplicated(multicell_drug_targets_in_model_target_nodes$Name), ]

write.csv(multicell_drug_targets_in_model_target_nodes,"large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\multicell_targets_nodes.csv")
write.csv(multicell_drug_targets_in_model_name,"large-scale-multicellular-modeling-of-the-arthritic-joint\\In silico simulations\\Multicell_drug_targets_in_model.csv")
