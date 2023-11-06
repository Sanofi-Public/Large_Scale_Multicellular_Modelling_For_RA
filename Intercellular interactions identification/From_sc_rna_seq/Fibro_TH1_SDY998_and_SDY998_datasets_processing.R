library(Seurat)
library(plyr)
library(sva)
library(xlsx)

setwd("Large-scale multicellular modeling of the arthritic joint\\From_sc_rna_seq")


#read fibroblast dataset
metadata <- read.delim("Large-scale multicellular modeling of the arthritic joint\\Datasets\\SDY998\\SDY998_metadata.tsv")

#filter dataset
metadata=metadata[metadata$Cell.type=="Fibroblast"& metadata$Diagnosis=="Rheumatoid arthritis" ,1]
seurat <- readRDS(file = "Large-scale multicellular modeling of the arthritic joint\\Datasets\\SDY998\\SDY998_seurat_1635258168512.Rds")
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
TH1_gene_expression_matrix <- read_excel("Large-scale multicellular modeling of the arthritic joint\\Datasets\\SDY998\\TH1_gene_expression_matrix.xlsx")
colnames(TH1_gene_expression_matrix)[1]="ID"

#join the two datasets

total_ra_data=join(as.data.frame(exp_matrix_filtered),
                   as.data.frame(TH1_gene_expression_matrix),
                   by="ID",
                   type="inner")

names=total_ra_data[,1]
total_ra_data=total_ra_data[,-c(1)]
rownames(total_ra_data)=names
colnames(total_ra_data)=c(rep("fibroblast",times = 1654),rep("TH1",times=12))


openxlsx::write.xlsx(exp_matrix_filtered_cells,
                     "total_RA_data_fibro_Th1.xlsx",
                     rowNames = F,colNames=T)



