library(Seurat)
library(plyr)
library(sva)
library(xlsx)

setwd("C:/Users/I0471594/OneDrive - Sanofi/PhD/Macrophage_fibro_interactions/SANOFI_datasets")


#read fibroblast dataset
metadata <- read.delim("BioTuring_SDY998_fibroblast_SC/SDY998_metadata.tsv")

#filter dataset
metadata=metadata[(metadata$Cell.type=="Fibroblast" | metadata$Cell.type=="Monocyte")& metadata$Diagnosis=="Rheumatoid arthritis" ,1]
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


metadata <- read.delim("BioTuring_SDY998_fibroblast_SC/SDY998_metadata.tsv")
metadata=metadata[(metadata$Cell.type=="Fibroblast" | metadata$Cell.type=="Monocyte")& metadata$Diagnosis=="Rheumatoid arthritis" ,c(1,7)]


exp_matrix_filtered_colnames=as.data.frame(colnames(exp_matrix_filtered))
colnames(exp_matrix_filtered_colnames)="Barcodes"
colnames=join(exp_matrix_filtered_colnames,metadata,type="left",by="Barcodes")

exp_matrix_filtered_cells=exp_matrix_filtered

colnames(exp_matrix_filtered_cells)=colnames$Cell.type
colnames(exp_matrix_filtered_cells)[1]="ID"
openxlsx::write.xlsx(exp_matrix_filtered_cells,
                     "total_RA_data_fibro_monocytes.xlsx",
                     rowNames = F,colNames=T)



