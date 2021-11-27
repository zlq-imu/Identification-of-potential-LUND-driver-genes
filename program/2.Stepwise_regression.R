##stepwise_regression.R

###Packages and parameters============================================

##Required packages====================================================

rm(list = ls())
#if (!requireNamespace("survival", quietly = TRUE))
#  install.packages("survival",repos = "http://cran.rstudio.com/")
library(survival)
#if (!requireNamespace("argparser", quietly = TRUE))
#  install.packages("argparser",repos = "http://cran.rstudio.com/")
library(argparser)

##Necessary parameters================================================

parser <- arg_parser("multi cox regression plot")
parser <- add_argument(parser, "--expr_file", help="expr data file")  
parser <- add_argument(parser, "--clinical_file", help="clinical data file") 
parser <- add_argument(parser, "--Siggene_file", help="significant genes from univariable cox regression analysis")
parser <- add_argument(parser, "--Outfile_stepwise_Sig", help="output file with significant genes through stepwise regression analysis")
parameters <- parse_args(parser)

##Input parameters====================================================

rm(parameters)
parameters <- parse_args(parser,c("--expr_file","E:/桌面/Github/raw_data/LUAD_exp_file(TPM).txt",
                            "--clinical_file","E:/桌面/Github/raw_data/LUAD_clinical_file.txt", 
                            "--Siggene_file","E:/桌面/Github/results/Univariate_cox_egression_signifcant_gene.txt", 
                            "--Outfile_stepwise_Sig","E:/桌面/Github/results/Stepwise_regression_signifcant_gene.txt" )
)

##Gene expression file================================================

expr <- read.table(parameters$expr_file,header = TRUE, row.names=1,stringsAsFactors=F,sep="\t",check.names = FALSE)

##Clinical file=======================================================

clinical <- read.table(parameters$clinical_file,header = TRUE,sep="\t")

##Significant genes===================================================

Siggene <- read.table(parameters$Siggene_file,header = T,sep="\t")

##Data preprocessing===================================================

if(sum(which(clinical$OS.time==0))!=0){clinical<-clinical[-which(clinical$OS.time==0),]}#delete samples with OS.time=0
if(sum(which(is.na(clinical$OS.time)))!=0){clinical<-clinical[-which(is.na(clinical$OS.time)),]} #delete samples with OS.time=NA
row.names(clinical)<-clinical$Sample
sample_overlap<-as.character(intersect(as.character(clinical[,1]),as.character(names(expr))))#obtain samples containing both expression data and clinical data  
clinical<-clinical[sample_overlap,]
expr<-expr[,sample_overlap]
expr<-expr[as.character(Siggene[,1]),]
gene_exp_t<-as.data.frame(t(expr))
gene_exp_t$Sample<-row.names(gene_exp_t)
clinical_expr<-merge(clinical[,c("Sample","OS","OS.time")],gene_exp_t,all=F)
surv<-Surv(clinical_expr$OS.time, clinical_expr$OS)
gene_mutlicox_all<-coxph(surv~., data=clinical_expr[,-(1:3)])

##Stepwise regression analysis================================================

gene_mutlicox_all<-step(gene_mutlicox_all,direction = "both") 
gene_summary<-summary(gene_mutlicox_all)
Stepwise_sig_gene<-rownames(gene_summary[[7]])
rownames(Siggene)<-Siggene$gene
Siggene<-Siggene[Stepwise_sig_gene,]

##Output file with Significant genes, HR, Pvalue through stepwise regression analysis===================

write.table(Siggene,parameters$Outfile_stepwise_Sig,quote=F,row.names = F,sep="\t")
print(sprintf('Complete'))
