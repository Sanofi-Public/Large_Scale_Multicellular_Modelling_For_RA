library(Seurat)
library(plyr)
library(sva)
library(xlsx)
library(dplyr)

setwd("Large-scale multicellular modeling of the arthritic joint\\From_sc_rna_seq")

#read macrophage dataset
macro_metadata <- read.delim("Large-scale multicellular modeling of the arthritic joint\\Datasets\\E_MTAB_8322\\E_MTAB_8322_metadata.tsv")


#keep only treatment naive RA
macro_metadata=macro_metadata[macro_metadata$Group=="Naive rheumatoid arthritis",c(1)]
exp_matrix=read.csv("Large-scale multicellular modeling of the arthritic joint\\Datasets\\E_MTAB_8322\\Merged.matrix.csv")
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
exp_matrix_filtered_fibro=read_excel("Large-scale multicellular modeling of the arthritic joint\\Datasets\\GSE109449\\GSE109449_RA_fibroblast_SC_matrix.xlsx")

colnames(exp_matrix_filtered_fibro)[1]="ID"
exp_matrix_filtered_fibro=exp_matrix_filtered_fibro %>% select(matches("ID|RA"))

#join the two datasets

total_ra_data=join(as.data.frame(exp_matrix_filtered_fibro),
                   as.data.frame(exp_matrix_filtered_macro),
                   by="ID",
                   type="inner")


# batch effect correction
names=total_ra_data[,1]
total_ra_data=total_ra_data[,-c(1)]
rownames(total_ra_data)=names
colnames(total_ra_data)=c(rep("fibroblast",times = 192),rep("macrophage",times=5815))
batch=c(rep("fibroblast",times = 192),rep("macrophage",times=5815))


total_ra_dataset_batch_removed_para= sva::ComBat(dat=as.matrix(total_ra_data),
                                                batch=batch, mod=NULL,
                                                par.prior=T, mean.only=FALSE )
rownames(total_ra_dataset_batch_removed_para)=names

openxlsx::write.xlsx(total_ra_dataset_batch_removed_para,
                     "total_RA_data_batch_effect_free.xlsx",
                     row.names = TRUE)
