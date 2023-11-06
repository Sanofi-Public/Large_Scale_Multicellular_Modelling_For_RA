library(Seurat)
library(plyr)
library(sva)
library(xlsx)

setwd("C:/Users/I0471594/OneDrive - Sanofi/PhD/Macrophage_fibro_interactions/SANOFI_datasets")


#read fibroblast dataset
exp_matrix_filtered_fibro=read_excel("C:/Users/I0471594/OneDrive - Sanofi (1)/Desktop/multicellular_map_model/datasets/GSE109449/GSE109449_RA_fibroblast_SC_matrix.xlsx")

colnames(exp_matrix_filtered_fibro)[1]="ID"
exp_matrix_filtered_fibro=exp_matrix_filtered_fibro %>% select(matches("ID|RA"))


#read Th1 gene expression matrix
TH1_gene_expression_matrix <- read_excel("C:/Users/I0471594/OneDrive - Sanofi (1)/Desktop/multicellular_map_model/datasets/SDY998/TH1_gene_expression_matrix.xlsx")
colnames(TH1_gene_expression_matrix)[1]="ID"

#join the two datasets

total_ra_data=join(as.data.frame(exp_matrix_filtered_fibro),
                   as.data.frame(TH1_gene_expression_matrix),
                   by="ID",
                   type="inner")

names=total_ra_data[,1]
total_ra_data=total_ra_data[,-c(1)]
rownames(total_ra_data)=names
colnames(total_ra_data)=c(rep("fibroblast",times = 192),rep("TH1",times=12))


# batch effect correction
names=total_ra_data[,1]
total_ra_data=total_ra_data[,-c(1)]
rownames(total_ra_data)=names
colnames(total_ra_data)=c(rep("fibroblast",times = 192),rep("TH1",times=12))
batch=c(rep("fibroblast",times = 192),rep("TH1",times=12))


total_ra_dataset_batch_removed_para= sva::ComBat(dat=as.matrix(total_ra_data),
                                                batch=batch, mod=NULL,
                                                par.prior=T, mean.only=FALSE )
rownames(total_ra_dataset_batch_removed_para)=names

openxlsx::write.xlsx(total_ra_dataset_batch_removed_para,
                     "total_RA_data_batch_effect_free.xlsx",
                     row.names = TRUE)
