library(devtools)
library("dplyr")
install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")
setwd("Large-scale multicellular modeling of the arthritic joint\\From_sc_rna_seq")
library(readxl)
library(xlsx)
library(icellnet)
library(Biobase)
library(Seurat)

db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"),
                          sep="\t",header = T,
                          check.names=FALSE, stringsAsFactors = FALSE,
                          na.strings = ""))
db.name.couple=name.lr.couple(db, type="Family")


#read data
total_RA_data_fibro_monocytes <- read_excel("total_RA_data_fibro_monocytes.xlsx")
fibro=as.data.frame(total_RA_data_fibro_monocytes[,grepl("Fibroblast", colnames(total_RA_data_fibro_monocytes))])
macro=as.data.frame(total_RA_data_fibro_monocytes[,grepl("Monocyte", colnames(total_RA_data_fibro_monocytes))])


fibro_average=as.data.frame(rowMeans(fibro,na.rm = TRUE))                                                          
symbol_fibro=total_RA_data_fibro_monocytes$ID
fibro_average_IDs=cbind(fibro_average,symbol_fibro)
colnames(fibro_average_IDs)=c("fibroblast","Symbol")

macro_average=as.data.frame(rowMeans(macro,na.rm = TRUE))                                                 
symbol_macro=total_RA_data_fibro_monocytes$ID
macro_average_IDs=cbind(macro_average,symbol_macro)
colnames(macro_average_IDs)=c("macrophage","Symbol")


TH1_gene_expression_matrix <- read_excel("TH1_gene_expression_matrix.xlsx")
TH1=TH1_gene_expression_matrix[,-c(1)]
TH1_average=as.data.frame(rowMeans(TH1,na.rm = TRUE))  
symbol_TH1=TH1_gene_expression_matrix[,1]
TH1_average_IDs=cbind(TH1_average,symbol_TH1)
colnames(TH1_average_IDs)=c("TH1","Symbol")




#compute icellnet scores
PC.target=data.frame("Class"=c("fibroblasts","monocytes","TH1"),
                     "ID"= c("fibroblasts","monocytes","TH1"),
                     "Cell_type"=c("fibroblasts","monocytes","TH1"))

rownames(PC.target)=c("fibroblasts","monocytes","TH1")
my.selection=c("fibroblasts","monocytes")

PC_data=cbind(macro_average,fibro_average)
colnames(PC_data)=c("monocytes","fibroblasts")
rownames(PC_data)=symbol_fibro

CC_data=TH1_average
colnames(CC_data)="TH1"
rownames(CC_data)=pull(TH1_gene_expression_matrix, ...1) 


score_TH1_to_fibro_and_mono= icellnet.score(direction="out",
                                                  PC.data=PC_data,#partner cell
                                                  CC.data=TH1_average_IDs, #central cell
                                                  CC.type = "RNAseq", 
                                                  PC.type = "RNAseq",
                                                  PC.target = PC.target,
                                                  PC=my.selection,
                                                  db = db)
score_TH1_to_fibro_and_mono_1=as.data.frame( score_TH1_to_fibro_and_mono[[1]])
score_TH1_to_fibro_and_mono_lr1=score_TH1_to_fibro_and_mono[[2]] 
                                               

write.xlsx(score_TH1_to_fibro_and_mono_lr1,"ICELLNET_scores_from_TH1_to_fibro_and_mono.xlsx")



score_fibro_and_mono_to_TH1= icellnet.score(direction="in",
                                  PC.data=PC_data,#partner cell
                                  CC.data=CC_data, #central cell
                                  CC.type = "RNAseq", 
                                  PC.type = "RNAseq",
                                  PC.target = PC.target,
                                  PC=my.selection,
                                  db = db)
score_fibro_and_mono_to_TH1_1=as.data.frame(score_fibro_and_mono_to_TH1[[1]])
score_fibro_and_mono_to_TH1_lr1=score_fibro_and_mono_to_TH1[[2]] 


write.xlsx(score_fibro_and_mono_to_TH1_lr1,"ICELLNET_scores_from_fibroblast_and_monocyte_to_TH1.xlsx")



#############################Plots##################################################
my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding", "ECM")

contrib_TH1_to_fibro_and_mono=LR.family.score(lr=score_TH1_to_fibro_and_mono_lr1,
                                              my.family=my.family, 
                                              db.couple=db.name.couple, plot=F) # table of contribution of each family of molecule to the scores



library(ggplot2)
ICELLNET_scores_TH1_to_fibro= as.numeric(contrib_TH1_to_fibro_and_mono[,1])
ICELLNET_scores_TH1_to_mono= as.numeric(contrib_TH1_to_fibro_and_mono[,2])
molecules_families = rownames(contrib_TH1_to_fibro_and_mono)


plot_TH1_to_fibro= ggplot(as.data.frame(contrib_TH1_to_fibro_and_mono),
                          aes(x="",y=ICELLNET_scores_TH1_to_fibro,fill=molecules_families))+
  
  geom_bar(width = 1, stat = "identity")

pie_TH1_to_fibro <-plot_TH1_to_fibro + coord_polar("y", start=0)+
  ggtitle("molecules families contributions to TH1/fibroblast interactions")

pie_TH1_to_fibro=pie_TH1_to_fibro +scale_fill_brewer(palette="Dark2")

pie_TH1_to_fibro

#########
plot_TH1_to_mono= ggplot(as.data.frame(contrib_TH1_to_fibro_and_mono),
                          aes(x="",y=ICELLNET_scores_TH1_to_mono,fill=molecules_families))+
  
  geom_bar(width = 1, stat = "identity")

pie_TH1_to_mono<-plot_TH1_to_mono + coord_polar("y", start=0)+
  ggtitle("molecules families contributions to TH1/monocytes interactions")

pie_TH1_to_mono=pie_TH1_to_mono +scale_fill_brewer(palette="Dark2")

pie_TH1_to_mono


#######################
my.family=c("Growth factor","Chemokine","Checkpoint","Cytokine","Notch family","Antigen binding", "ECM")

contrib_fibro_and_mono_to_TH1=LR.family.score(lr=score_fibro_and_mono_to_TH1_lr1,
                                              my.family=my.family, 
                                              db.couple=db.name.couple, plot=F) # table of contribution of each family of molecule to the scores


ICELLNET_scores_fibro_to_TH1= as.numeric(contrib_fibro_and_mono_to_TH1[,1])
ICELLNET_scores_mono_to_TH1= as.numeric(contrib_fibro_and_mono_to_TH1[,2])
molecules_families = rownames(contrib_fibro_and_mono_to_TH1)


plot_fibro_to_TH1= ggplot(as.data.frame(contrib_fibro_and_mono_to_TH1),
                          aes(x="",y=ICELLNET_scores_fibro_to_TH1,fill=molecules_families))+
  
  geom_bar(width = 1, stat = "identity")

pie_fibro_to_TH1 <-plot_fibro_to_TH1 + coord_polar("y", start=0)+
  ggtitle("molecules families contributions to fibroblasts/TH1 interactions")

pie_fibro_to_TH1=pie_fibro_to_TH1 +scale_fill_brewer(palette="Dark2")

pie_fibro_to_TH1

########
plot_mono_to_TH1= ggplot(as.data.frame(contrib_fibro_and_mono_to_TH1),
                          aes(x="",y=ICELLNET_scores_mono_to_TH1,fill=molecules_families))+
  
  geom_bar(width = 1, stat = "identity")

pie_mono_to_TH1 <-plot_mono_to_TH1 + coord_polar("y", start=0)+
  ggtitle("molecules families contributions to monocytes/TH1 interactions")

pie_mono_to_TH1=pie_mono_to_TH1 +scale_fill_brewer(palette="Dark2")

pie_mono_to_TH1
