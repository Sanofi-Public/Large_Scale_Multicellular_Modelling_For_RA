setwd("a-computational-framework-to-build-and-calibrate-large-scale-boolean-models-main/GSE97779_dataset_analysis")
library(openxlsx)
library(preprocessCore)
library(biomaRt)
library(AnnotationHub)
library(plyr)
library(dplyr)
library(limma)
library(nortest)
#read the expression matrix
expression_matrix= read.delim("GSE97779_series_matrix.txt", comment.char="!",row.names=1)

#quantile normalization 
norm_expression_matrix=normalize.quantiles(as.matrix(expression_matrix))
norm_expression_matrix=log2(norm_expression_matrix)
#restore rows' and columns' names after normalization 
colnames(norm_expression_matrix)=c("Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","RA","RA","RA","RA","RA","RA","RA","RA","RA")
rownames(norm_expression_matrix) = rownames(expression_matrix)

#test the normality and variance of Ctrl and RA populations: 
Ctrl_samples= norm_expression_matrix[,c(1:5)]
RA_samples = norm_expression_matrix[,c(6:ncol(norm_expression_matrix))]
hist(unlist(log2(Ctrl_samples)), breaks=100,main="Gene expression distribution in control samples")
hist(unlist(log2(RA_samples)), breaks=100,main="Gene expression distribution in RA samples")
# with Anderson-Darling test for normality :P-value < 0.05 = not normal 
ad.test(Ctrl_samples)
ad.test(RA_samples)


# assign samples to groups and set up design matrix
gsms <- "11111000000000"
sml <- strsplit(gsms, split="")[[1]]
gs <- factor(sml)
groups <- make.names(c("RA","control"))
levels(gs) <- groups

design <- model.matrix( ~ 0 + gs )
colnames(design) <- levels(gs)

fit <- lmFit(norm_expression_matrix, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000000)

tT <- subset(tT, select=c("adj.P.Val","logFC"))
tT$ID=row.names(tT)

GPL570 <- read.delim2("GPL570-55999.txt", comment.char="#")

GPL570 = GPL570[,c("ID","Gene.Symbol")]
#join annotation file to the expression matrix

annotated_DE_genes=join(GPL570,tT,type="right",by="ID") 

annotated_DE_genes=annotated_DE_genes[!(annotated_DE_genes$Gene.Symbol==""),]
annotated_DE_genes=annotated_DE_genes[annotated_DE_genes$adj.P.Val<=0.05,]
write.xlsx(annotated_DE_genes,'DEG_macrophage.xlsx')


