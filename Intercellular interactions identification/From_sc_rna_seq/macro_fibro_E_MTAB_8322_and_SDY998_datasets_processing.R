library(Seurat)
library(plyr)
library(sva)
library(xlsx)

setwd("C:/Users/I0471594/OneDrive - Sanofi/PhD/Macrophage_fibro_interactions/SANOFI_datasets")

#read macrophage dataset
macro_metadata <- read.delim("Bioturing_E_MTAB_8322_macrophage_SC/STM/E_MTAB_8322_metadata.tsv")


#keep only treatment naive RA
macro_metadata=macro_metadata[macro_metadata$Group=="Naive rheumatoid arthritis",c(1)]
exp_matrix=read.csv("Bioturing_E_MTAB_8322_macrophage_SC/STM/Merged.matrix.csv")
exp_matrix_filtered_macro=exp_matrix[,macro_metadata]
rownames(exp_matrix_filtered_macro)=exp_matrix[,1]
#remove genes with var equal zero 
exp_matrix_filtered_macro=as.data.frame(exp_matrix_filtered_macro[apply(as.data.frame(exp_matrix_filtered_macro),1, var) != 0,])
names_macro=rownames(exp_matrix_filtered_macro)
exp_matrix_filtered_macro=cbind(names_macro,exp_matrix_filtered_macro) #normalized dataset with NormalizeData
colnames(exp_matrix_filtered_macro)[1]="ID"
exp_matrix_filtered_macro$ID=as.character(exp_matrix_filtered_macro$ID)
rownames(exp_matrix_filtered_macro)=NULL



#read fibroblast dataset
fibro_metadata <- read.delim("BioTuring_SDY998_fibroblast_SC/SDY998_metadata.tsv")

#filter dataset
fibro_metadata=fibro_metadata[fibro_metadata$Cell.type=="Fibroblast" & fibro_metadata$Diagnosis=="Rheumatoid arthritis",c(1)]
seurat <- readRDS(file = "BioTuring_SDY998_fibroblast_SC/SDY998_seurat_1635258168512.Rds")
exp_matrix_fibro=as.matrix(GetAssayData(object = seurat, slot = "counts"))
exp_matrix_filtered_fibro=as.data.frame(exp_matrix_fibro[,fibro_metadata]) # normalized with log2(CPM)
#remove genes with null variance
exp_matrix_filtered_fibro=as.data.frame(exp_matrix_filtered_fibro[apply(as.data.frame(exp_matrix_filtered_fibro),1, var) != 0,])

name_fibro=rownames(exp_matrix_filtered_fibro)
exp_matrix_filtered_fibro=cbind(name_fibro,exp_matrix_filtered_fibro)
colnames(exp_matrix_filtered_fibro)[1]="ID"
exp_matrix_filtered_fibro$ID=as.character(exp_matrix_filtered_fibro$ID)
rownames(exp_matrix_filtered_fibro)=NULL

#join the two datasets

total_ra_data=join(as.data.frame(exp_matrix_filtered_fibro),
                   as.data.frame(exp_matrix_filtered_macro),
                   by="ID",
                   type="inner")


openxlsx::write.xlsx(total_ra_data,
                     "total_RA_data.xlsx",
                     row.names = TRUE)


# batch effect correction
names=total_ra_data[,1]
total_ra_data=total_ra_data[,-c(1)]
rownames(total_ra_data)=names
colnames(total_ra_data)=c(rep("fibroblast",times = 1654),rep("macrophage",times=5815))
batch=c(rep("fibroblast",times = 1654),rep("macrophage",times=5815))


total_ra_dataset_batch_removed_para= sva::ComBat(dat=as.matrix(total_ra_data),
                                                batch=batch, mod=NULL,
                                                par.prior=T, mean.only=FALSE )
rownames(total_ra_dataset_batch_removed_para)=names

openxlsx::write.xlsx(total_ra_dataset_batch_removed_para,
                     "total_RA_data_batch_effect_free.xlsx",
                     row.names = TRUE)
