library(Seurat)
library(plyr)
library(sva)
library(xlsx)

setwd("C:/Users/I0471594/OneDrive - Sanofi/PhD/Macrophage_fibro_interactions/SANOFI_datasets")

#read macrophage dataset
metadata <- read.delim("BioTuring_SDY998_fibroblast_SC/SDY998_metadata.tsv")

#filter dataset
metadata=metadata[(metadata$Cell.type=="Monocyte")& metadata$Diagnosis=="Rheumatoid arthritis" ,]
seurat <- readRDS(file = "BioTuring_SDY998_fibroblast_SC/SDY998_seurat_1635258168512.Rds")
exp_matrix=as.matrix(GetAssayData(object = seurat, slot = "counts"))
exp_matrix_filtered=as.data.frame(exp_matrix[,metadata]) # normalized with log2(CPM)
#remove genes with null variance
exp_matrix_filtered=as.data.frame(exp_matrix_filtered[apply(as.data.frame(exp_matrix_filtered),1, var) != 0,])

name=rownames(exp_matrix_filtered)
exp_matrix_filtered=cbind(name,exp_matrix_filtered)
colnames(exp_matrix_filtered)[1]="ID"
exp_matrix_filtered$ID=as.character(exp_matrix_filtered$ID)
rownames(exp_matrix_filtered)=NULL

#read Th1 gene expression matrix
TH1_gene_expression_matrix <- read_excel("C:/Users/I0471594/OneDrive - Sanofi (1)/Desktop/multicellular_map_model/datasets/SDY998/TH1_gene_expression_matrix.xlsx")
colnames(TH1_gene_expression_matrix)[1]="ID"

#join the two datasets

total_ra_data=join(as.data.frame(exp_matrix_filtered),
                   as.data.frame(TH1_gene_expression_matrix),
                   by="ID",
                   type="inner")

names=total_ra_data[,1]
total_ra_data=total_ra_data[,-c(1)]
rownames(total_ra_data)=names
colnames(total_ra_data)=c(rep("macrophage",times = 602),rep("TH1",times=12))

openxlsx::write.xlsx(exp_matrix_filtered_cells,
                     "total_RA_data_monocytes_Th1.xlsx",
                     rowNames = F,colNames=T)



