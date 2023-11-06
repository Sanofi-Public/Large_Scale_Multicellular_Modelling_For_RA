library(dplyr)
library(fuzzyjoin)
library(xlsx)
setwd("C:/Users/I0471594/OneDrive - Sanofi/PhD/Macrophage_fibro_interactions")
Cell_phone_Db_interaction_curated <- read.csv("C:/Users/I0471594/OneDrive - Sanofi/PhD/Macrophage_fibro_interactions/Cell_phone_Db_interaction_curated.csv", sep=";")
M1_fibroblast_list <- read.table("C:/Users/I0471594/OneDrive - Sanofi/PhD/Macrophage_fibro_interactions/M2_fibroblast_list.txt", quote="\"", comment.char="")
Cell_phone_Db_M1_fibro= data.frame(ncol(4))

n_ligne=1
l=1
for (i in Cell_phone_Db_interaction_curated$protein_name_a){
for (j in M1_fibroblast_list$V1)
{if (grepl(j,i, fixed = TRUE)==TRUE){Cell_phone_Db_M1_fibro[l,1]=Cell_phone_Db_interaction_curated[n_ligne,1];
Cell_phone_Db_M1_fibro[l,2]=Cell_phone_Db_interaction_curated[n_ligne,2];
Cell_phone_Db_M1_fibro[l,3]=Cell_phone_Db_interaction_curated[n_ligne,3];
Cell_phone_Db_M1_fibro[l,4]=Cell_phone_Db_interaction_curated[n_ligne,4];
l=l+1}
}
n_ligne=n_ligne+1}

colnames(Cell_phone_Db_M1_fibro)=colnames(Cell_phone_Db_interaction_curated)


#filter with protein_name_b
CPDB_M1_fibro=data.frame()
n_ligne=1
l=1
for (i in Cell_phone_Db_M1_fibro$protein_name_b){
  for (j in M1_fibroblast_list$V1)
  {if (grepl(j,i, fixed = TRUE)==TRUE){CPDB_M1_fibro[l,1]=Cell_phone_Db_M1_fibro[n_ligne,1];
  CPDB_M1_fibro[l,2]=Cell_phone_Db_M1_fibro[n_ligne,2];
  CPDB_M1_fibro[l,3]=Cell_phone_Db_M1_fibro[n_ligne,3];
  CPDB_M1_fibro[l,4]=Cell_phone_Db_M1_fibro[n_ligne,4];
  l=l+1}
  }
  n_ligne=n_ligne+1}

colnames(CPDB_M1_fibro)=colnames(Cell_phone_Db_interaction_curated)

write.xlsx(CPDB_M1_fibro,"CellphoneDb_M2_fibroblast_interaction.xlsx")

