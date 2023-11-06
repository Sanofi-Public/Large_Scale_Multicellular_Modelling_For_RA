library(Seurat)
library(plyr)
library(sva)
library(xlsx)

setwd("Large-scale multicellular modeling of the arthritic joint\\From_sc_rna_seq")


#read macrophage dataset
macro_metadata <- read.delim("Large-scale multicellular modeling of the arthritic joint\\Datasets\\SDY998\\SDY998_metadata.tsv")

#filter dataset
macro_metadata=macro_metadata[macro_metadata$Cell.type=="Monocyte" & macro_metadata$Diagnosis=="Rheumatoid arthritis",c(1)]
seurat <- readRDS(file = "Large-scale multicellular modeling of the arthritic joint\\Datasets\\SDY998\\SDY998_seurat_1635258168512.Rds")
exp_matrix_macro=as.matrix(GetAssayData(object = seurat, slot = "counts"))
exp_matrix_filtered_macro=as.data.frame(exp_matrix_macro[,macro_metadata]) # normalized with log2(CPM)
#remove genes with null variance
exp_matrix_filtered_macro=as.data.frame(exp_matrix_filtered_macro[apply(as.data.frame(exp_matrix_filtered_macro),1, var) != 0,])

name_acro=rownames(exp_matrix_filtered_macro)
exp_matrix_filtered_macro=cbind(name_macro,exp_matrix_filtered_macro)
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


openxlsx::write.xlsx(total_ra_data,
                     "total_RA_data.xlsx",
                     row.names = TRUE)


# batch effect correction
names=total_ra_data[,1]
total_ra_data=total_ra_data[,-c(1)]
rownames(total_ra_data)=names
colnames(total_ra_data)=c(rep("fibroblast",times = 192),rep("macrophage",times=602))
batch=c(rep("fibroblast",times = 192),rep("macrophage",times=5815))


total_ra_dataset_batch_removed_para= sva::ComBat(dat=as.matrix(total_ra_data),
                                                batch=batch, mod=NULL,
                                                par.prior=T, mean.only=FALSE )
rownames(total_ra_dataset_batch_removed_para)=names

openxlsx::write.xlsx(total_ra_dataset_batch_removed_para,
                     "total_RA_data_batch_effect_free.xlsx",
                     row.names = TRUE)
