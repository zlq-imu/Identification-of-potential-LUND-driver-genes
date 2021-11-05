##univariate_cox.R

###Packages and parameters============================================

##Required packages====================================================

rm(list = ls())
#if (!requireNamespace("survival", quietly = TRUE))
#  install.packages("survival", repo = "http://cran.rstudio.com/")
library(survival)
#if (!requireNamespace("argparser", quietly = TRUE))
#  install.packages("argparser", repo = "http://cran.rstudio.com/")
library(argparser)

##Necessary parameters================================================

parser <- arg_parser("Single cox regression")
parser <- add_argument(parser, "--p", help="p.value", type="numeric", default=0.05)
parser <- add_argument(parser, "--fdr", help="FDR", type="numeric", default=0.05) 
parser <- add_argument(parser, "--q", help="percent 1/2 or 1/4", type="numeric", default=0.5)
parser <- add_argument(parser, "--diffFilter", help="Filter value for the difference between the two survival curves", type="numeric", default=0.1)
parser <- add_argument(parser, "--expr_file", help="expr data file") 
parser <- add_argument(parser, "--clinical_file", help="clinical data file") 
parser <- add_argument(parser, "--Siggene_file", help="Siggene data file") 
parser <- add_argument(parser, "--outfilesinglecox", help="output file with significant genes through univariate cox analysis")
parameters <- parse_args(parser)

##Input parameters=====================================================

rm(parameters)
parameters <- parse_args(parser,c("-f", "1",
                             "-p", "0.05",
                             "--diffFilter","0.1",
                             "--expr_file","E:/桌面/Github/raw_data/LUAD_exp_file(TPM).txt",
                             "--clinical_file","E:/桌面/Github/raw_data/LUAD_clinical_file.txt",
                             "--Siggene_file","E:/桌面/Github/raw_data/DK79Gs_only_name.txt",
                             "--outfilesinglecox","E:/桌面/Github/results/Univariate_cox_egression_signifcant_gene.txt"
                             )
)


##Gene expression file=================================================

expr <- read.table(parameters$expr_file,header = TRUE, row.names=1,stringsAsFactors=F,sep="\t",check.names = FALSE)

##Clinical file========================================================

clinical <- read.table(parameters$clinical_file,header = TRUE,sep="\t")

##Significant genes====================================================

Siggene <- read.table(parameters$Siggene_file,header = F,sep="\t")

##Data preprocessing===================================================

if(sum(which(clinical$OS.time==0))!=0){clinical<-clinical[-which(clinical$OS.time==0),]}   #delete samples with OS.time=0
if(sum(which(is.na(clinical$OS.time)))!=0){clinical<-clinical[-which(is.na(clinical$OS.time)),]}  #delete samples with OS.time=NA
row.names(clinical)<-clinical$Sample
Sample_overlap<-as.character(intersect(as.character(clinical[,1]),as.character(names(expr))))
clinical<-clinical[Sample_overlap,]
expr<-expr[,Sample_overlap]
expr<-expr[as.character(Siggene[,1]),]#obtain expression file for siggnificant genes
Sur_Obj<-as.data.frame(Surv(as.numeric(clinical$OS.time),as.numeric(clinical$OS)))  #create a survival object
row.names(Sur_Obj)<-row.names(clinical)
expr<-expr[,row.names(Sur_Obj)]

##Univariate cox analysis================================================

gene_p<-data.frame()
i=1
for(x in unique(rownames(expr))){
  gene_expr<-expr[x,]
  gene_expr_sort<-sort(gene_expr) 
  gene_expr_sort_zhuanzhi<-as.data.frame(t(gene_expr_sort))
  if(parameters$q!=0.5){
    low<-floor(nrow(clinical)*parameters$q)
    high<-floor(nrow(clinical)*parameters$q)
    gene_expr_sort_zhuanzhi[,2]<-rep(c("low","NA","high"),c(low,nrow(clinical)-low-high,high))
  }else{
    gene_expr_sort_zhuanzhi[,2]<-rep(c("low","high"),c(floor(nrow(clinical)*parameters$q),ceiling(nrow(clinical)*parameters$q)))
  }
  gene_expr_sort_zhuanzhi[,3]<-names(gene_expr_sort) 
  names(gene_expr_sort_zhuanzhi)<-c("expr","class","Sample")
  Sur_Obj[,2]<-row.names(Sur_Obj)
  names(Sur_Obj)<-c("time","Sample")
  combine<-merge(gene_expr_sort_zhuanzhi,Sur_Obj)
  if(parameters$q!=0.5){
  combine<-combine[-which(combine$class=="NA"),]}
  combine1<-merge(clinical,combine,by="Sample")
  combine1<-combine1[,-c(4:11,14)]
  combine1 <- combine1[,c(1,3,2,4,5)]
  a<-combine1[which(combine$class=="low"),]
  a<-a[,-c(5)]
  b<-combine1[which(combine$class=="high"),]
  b<-b[,-c(5)]
  surTab1=summary(survfit(Surv(OS.time, OS) ~ 1, data = a))
  surTab2=summary(survfit(Surv(OS.time, OS) ~ 1, data = b))
  survivalTab1=cbind(time=surTab1$time, surv=surTab1$surv,lower=surTab1$lower,upper=surTab1$upper)
	   survivalTab1=survivalTab1[survivalTab1[,"time"]<5,]
	   if(class(survivalTab1)=="matrix"){
	     survivalTab1=survivalTab1[nrow(survivalTab1),]
	   }
	   survivalTab2=cbind(time=surTab2$time, surv=surTab2$surv,lower=surTab2$lower,upper=surTab2$upper)
	   survivalTab2=survivalTab2[survivalTab2[,"time"]<5,]
	   if(class(survivalTab2)=="matrix"){
	     survivalTab2=survivalTab2[nrow(survivalTab2),]
	   }
	fiveYearsDiff=abs(survivalTab1["surv"]-survivalTab2["surv"])
  sdf= survdiff(combine$time~combine$class) 
  p.val<-1-pchisq(sdf$chisq,length(sdf$n)-1)
  gene_p[i,1]<-x
  gene_p[i,2]<-p.val
  gene_p[i,3]<-fiveYearsDiff
  i=i+1
}
gene_p[,4]=p.adjust(gene_p[,2], method ="fdr",n=length(gene_p[,2]))
names(gene_p)<-c("gene","p","curve_diff","adj.p")
gene_p_Sig<-gene_p[which(gene_p$p<parameters$p),]
gene_curve_diff_Sig<-gene_p[which(gene_p$curve_diff>parameters$diffFilter),]
gene_fdr_Sig<-gene_curve_diff_Sig[which(gene_curve_diff_Sig$adj.p<parameters$fdr),]
gene_marker <- c() 
gene_pvalue<-c() 
gene_HR<-c()
gene_low<-c() #Low 95%CI
gene_high<-c() #High 95%CI
curve_diff<-c() 
for(x in gene_fdr_Sig[,1]){
  cox.fit_dan <- coxph(Sur_Obj$time~as.numeric(expr[x,]))
  coxresult<-summary(cox.fit_dan)
  pvalue=coxresult$coefficients[5]
  gene_marker<-c(gene_marker,x)
  gene_pvalue<-c(gene_pvalue,pvalue)
  gene_HR<-c(gene_HR,coxresult$conf.int[1])
  gene_low<-c(gene_low,coxresult$conf.int[3])
  gene_high<-c(gene_high,coxresult$conf.int[4])
}
p.adj=p.adjust(gene_pvalue, method ="fdr",n=length(gene_pvalue))
gene_singlecox_all<-cbind(gene_marker,gene_pvalue,gene_HR,gene_low,gene_high,p.adj)
gene_singlecox_all<-data.frame(gene_singlecox_all)
names(gene_singlecox_all)<-c("gene","P_value","HR","Low 95%CI","High 95%CI","fdr")
gene_singlecox_sig<-gene_singlecox_all[which(as.numeric(as.character(gene_singlecox_all$P_value))<parameters$p),]
gene_singlecox_sig<-gene_singlecox_sig[which(as.numeric(as.character(gene_singlecox_sig$fdr ))<parameters$fdr),]

##Output file with Significant genes, HR, Pvalue through univariate cox analysis===================

write.table(gene_singlecox_sig,parameters$outfilesinglecox,quote=F,row.names = F,sep="\t")
print(sprintf('Complete'))
