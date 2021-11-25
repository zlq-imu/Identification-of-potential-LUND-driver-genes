##K_Iteration.R


##Required pPackages====================================================

rm(list = ls())
#if (!requireNamespace("survival", quietly = TRUE))
#  install.packages("survival",repos = "http://cran.rstudio.com/")
library(survival)
#if (!requireNamespace("survminer", quietly = TRUE))
#  install.packages("survminer",repos = "http://cran.rstudio.com/")
library(survminer)
#if (!requireNamespace("argparser", quietly = TRUE))
#  install.packages("argparser",repos = "http://cran.rstudio.com/")
library(argparser)
#if (!requireNamespace("timeROC", quietly = TRUE))
# install.packages("timeROC",repos = "http://cran.rstudio.com/")
library(timeROC)
#if (!requireNamespace("glmnet", quietly = TRUE))
#  install.packages("glmnet",repos = "http://cran.rstudio.com/")
library(glmnet)

##Necessary parameters================================================

parser <- arg_parser("nfold lasso regression plot")
parser <- add_argument(parser, "--expr_file", help="expr data file") 
parser <- add_argument(parser, "--clinical_file", help="clinical data file") 
parser <- add_argument(parser, "--Siggene_file", help="significant genes from stepwise regression analysis")
parser <- add_argument(parser, "--outfile_coef", help="output file with coefficient for each significant gene")
parser <- add_argument(parser, "--outfile_clinical", help="output file with clinical data, exp data and Risk_Score data")
parser <- add_argument(parser, "--Outfile_ite_Sig", help="output file with significant genes through K-times LASSO iteration regression analysis")

##Adjustable parameters===============================================

parser <- add_argument(parser, "--Ktimes",help="times of LASSO iteration regression", type="numeric",default=1000)
parser <- add_argument(parser, "--nfold",help="n-fold cross-validation", type="numeric",default=10)
parser <- add_argument(parser, "--ROCcutyear1",help="survivalROCplot cut year",default=3)
parser <- add_argument(parser, "--ROCcutyear2",help="survivalROCplot cut year",default=5)
parser <- add_argument(parser, "--pdf", help="path for output pictures")
parser <- add_argument(parser, "--survplot", help="name for survival curve")
parser <- add_argument(parser, "--ROCplot",help="survivalROCplot")
parser <- add_argument(parser, "--PHtest",help="cox ph test")
parser <- add_argument(parser, "--coxdiagnostics",help="Identifying influential cases")
parameters <- parse_args(parser)

##Input parameters====================================================

rm(parameters)
parameters <- parse_args(parser,c("--expr_file","E:/桌面/Github/raw_data/LUAD_exp_file(TPM).txt",
                            "--clinical_file","E:/桌面/Github/raw_data/LUAD_clinical_file.txt",
                            "--Siggene_file","E:/桌面/Github/results/Stepwise_regression_signifcant_gene.txt", 
                            "--outfile_coef","E:/桌面/Github/results/Cofficients_of_K_iterations_for_sig_genes.txt",                                                                                                 
                            "--outfile_clinical","E:/桌面/Github/results/Clinical_expr.txt", 
                            "--Outfile_ite_Sig","E:/桌面/Github/results/K_Iteration_signifcant_gene.txt",                                                                               
                            "--survplot","survival_regression",
                            "--ROCplot","ROCplot",
                            "--PHtest","PHtest",
                            "--coxdiagnostics","coxdiagnostics",
                            "--Ktimes","1000","--nfold","10",
                            "--ROCcutyear1","3","--ROCcutyear2","5",
                            "--pdf","E:/桌面/Github/results/image/")                             
)
                                                                                                      

##Gene expression file================================================

expr <- read.table(parameters$expr_file,header = TRUE, row.names=1,stringsAsFactors=F,sep="\t",check.names = FALSE)

##Clinical file=======================================================

clinical <- read.table(parameters$clinical_file,header = TRUE,sep="\t")
for (i in 1:nrow(clinical)){if (clinical$OS.time[i]>=10){clinical$OS.time[i]=10}} #OS.time that is over 10 years is replaced by 10 years

##Significant genes===================================================

Siggene <- read.table(parameters$Siggene_file,header = T,sep="\t")

##Perform K-times LASSO iterations of n-fold cross-validation=========

k.nfold.cv.glmnet=function(dset,           # a data.frame
                           yvar,           # name or numeric index of response variable in dset
                           xvars,          # names or numeric indices of predictor variables in dset
                           fam,            # indicates model family to fit, see help(glmnet) for more details
                           alpha=1,        # elastic net mixing parameter, alpha = 1 is lasso, alpha=0 is ridge, see help(glmnet) for more info
                           nfolds=parameters$nfold, 
                           k=parameters$Ktimes) 
{
  # Check input data
  if (!is.data.frame(dset))
    stop("dset must be a data.frame")
  if (is.character(yvar))
  {
    if (!is.element(yvar,colnames(dset)))
      stop("yvar must be the numeric index or character name of a column of the response variable in dset.")
  }
  if (is.character(xvars))
  {
    if (any(!is.element(xvars,colnames(dset))))
      stop("xvars must be a vector of numeric indices or character names of the candidate predictor variables in dset.")
  }
  y=dset[,yvar]                   # extract response variable
  X=dset[,xvars]                  # extract candidate predictor variables
  X=as.matrix(X)                  # represent candidate predictors as a matrix
  B=matrix(NA,ncol(X),k)          # initialize matrix of coefficient estimates
  rownames(B)=colnames(X)         # label rows of coefficient estimates with variable names
  colnames(B)=paste0("CV_",1:k)   # label columns of coefficient estimates with 
  for (i in 1:k)                  # loop over k iterations of nfold cross-validation
  {
    message(paste0("Performing iteration ",i," of ",k, 
                   " iterations of ",nfolds,
                   "-fold cross-validation: ",
                   date()))
    fit.res=cv.glmnet(y=y,x=X,                  # perform the n-fold cross-validation of glmnet fit
                      nfolds=nfolds,            
                      family=fam)
    fit.est=as.numeric(coef(fit.res,            # extract the estimates from this iteration
                            s="lambda.min"))
    B[,i]=fit.est                               # place those estimates in the coefficient estimate matrix
  }
  
  class(B)="k.nfold.cv.glmnet"                  # assign a class to the result matrix for printing
  
  return(B)                                     # return the result matrix
}  # end of function k.nfold.cv.glmnet

##Display results of K-times LASSO iterations of n-fold cross-validation============

print.k.nfold.cv.glmnet=function(B,alpha=0.05,
                                 max.rows=24)# B is the result of k.nfold.cv.glmnet

{
  non.zero=rowMeans(B!=0)
  mean.coef=rowMeans(B)
  median.coef=apply(B,1,median)
  q1.coef=apply(B,1,quantile,0.25)
  q3.coef=apply(B,1,quantile,0.75)
  lower.bound=apply(B,1,quantile,alpha/2)
  upper.bound=apply(B,1,quantile,1-alpha/2)
  
  
  res=cbind(non.zero=non.zero,
            mean.coef=mean.coef,
            median.coef=median.coef,
            q1.coef=q1.coef,
            q3.coef=q3.coef,
            lower.bound=lower.bound,
            upper.bound=upper.bound)
  colnames(res)[6:7]=paste0(c("lower.","upper."),
                            1-alpha,".bound")
            
  ord=rev(order(res[,"non.zero"],abs(res[,"mean.coef"])))
  res=res[ord,]
  n.non.zero=sum(non.zero!=0)
  if (n.non.zero>max.rows)
    warning(paste0("There are ",n.non.zero,
                   " variables with a non-zero ",
                   " coefficient in at least one iteration.  ",
                   "Only the top ",max.rows," are shown here."))
  print(head(res,max.rows))
}

##Data preprocessing===================================================

if(sum(which(clinical$OS.time==0))!=0){clinical<-clinical[-which(clinical$OS.time==0),]}#delete samples with OS.time=0
if(sum(which(is.na(clinical$OS.time)))!=0){  clinical<-clinical[-which(is.na(clinical$OS.time)),]} #delete samples with OS.time=NA
row.names(clinical)<-clinical$Sample
sample_overlap<-as.character(intersect(as.character(clinical[,1]),as.character(names(expr))))#obtain samples containing both expression data and clinical data 
clinical<-clinical[sample_overlap,]
expr<-expr[,sample_overlap]
expr<-expr[as.character(Siggene[,1]),]
gene_exp_t<-as.data.frame(t(expr))
gene_exp_t$Sample<-row.names(gene_exp_t)
clinical_expr<-merge(clinical[,c("Sample","OS","OS.time")],gene_exp_t,all=F)
X=as.matrix(clinical_expr[,4:length(clinical_expr)])
obs.surv<-Surv(clinical_expr$OS.time, clinical_expr$OS)
dset=cbind.data.frame(obs.surv=obs.surv,X)
xvars<-2:ncol(dset)
sim.res=k.nfold.cv.glmnet(dset,
                          yvar="obs.surv",
                          xvars=xvars,
                          fam="cox",k=parameters$Ktimes)

###Retain genes with non-zero coefficient in at least 95% of K-times LASSO iterations=========
c=c() 
for(i in 1:nrow(sim.res) ){ 
if(sum(sim.res[i,]==0)>= ncol(sim.res)*0.05){ 
c=append(c,i) 
} } 
sim.res=sim.res[-c,]

##Output file with Significant genes and coefficient for each significant gene================

write.table(t(sim.res),parameters$outfile_coef,quote=F,row.names = T,sep="\t")
rownames(Siggene)<-Siggene$gene
Siggene<-Siggene[rownames(sim.res),]
write.table(Siggene,parameters$Outfile_ite_Sig,quote=F,row.names = T,sep="\t")

##Risk assessment model construction====================================

gene_coef=rowMeans(sim.res)
clinical_expr<-clinical_expr[,c("Sample","OS","OS.time",rownames(sim.res))]
clinical_expr[,ncol(clinical_expr)+1]<-as.matrix(clinical_expr[,rownames(sim.res)])%*%gene_coef
colnames(clinical_expr)[ncol(clinical_expr)]=c("Risk_Score")

##Output file with clinical data, exp data and Risk_Score data==========

write.table(clinical_expr,parameters$outfile_clinical,quote=F,row.names = F,sep="\t")

##Draw survival curve and obtain the cutoff for Risk Score==============

clinical_expr$Risk_Score<-as.numeric(clinical_expr$Risk_Score)
surv.cut<-surv_cutpoint(clinical_expr, time="OS.time", event="OS",variables="Risk_Score")
surv.cat<-surv_categorize(surv.cut)
surv.fit<-survfit(Surv(OS.time,OS)~Risk_Score, data=surv.cat)
surv.dif<-survdiff(Surv(OS.time,OS)~Risk_Score, data=surv.cat)
pValue=paste0("=",format(1-pchisq(surv.dif$chisq,df=1),scientific=TRUE,digit=4))
surPlot<-ggsurvplot(surv.fit, risk.table=TRUE, pval.method = TRUE, pval=paste0("p",pValue), pval.size=5,add.all=TRUE)
pdf(paste(parameters$pdf,parameters$survplot,".pdf",sep = ""),onefile = FALSE, width=6,height=6)
print(surPlot)
dev.off()

##Draw the ROC cruve====================================================

cutoff1 <- parameters$ROCcutyear1
cutoff2 <- parameters$ROCcutyear2
pdf(paste(parameters$pdf,parameters$ROCplot,".pdf",sep = ""),width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
aucText=c()
rocCol=c("#377EB8","#4DAF4A")

##Draw the ROC cruve at 3 year========================================== 

mRNA_ROC_3<-timeROC(T=clinical_expr$OS.time,delta=clinical_expr$OS,marker=clinical_expr$Risk_Score,cause=1, weighting="cox",times=c(cutoff1),iid=F) 
plot(mRNA_ROC_3$FP, mRNA_ROC_3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1],
		 xlab="1-Specificity", ylab="Sensitivity",lwd = 2, cex.main=1.2, 
		 cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("three year"," (AUC=",sprintf("%.3f",mRNA_ROC_3$AUC[2]),")"))
abline(0,1)

##Draw the ROC cruve at 5 year========================================== 

mRNA_ROC_5<-timeROC(T=clinical_expr$OS.time,delta=clinical_expr$OS,marker=clinical_expr$Risk_Score,cause=1, weighting="cox",times=c(cutoff2),iid=F) 
aucText=c(aucText,paste0("five year"," (AUC=",sprintf("%.3f",mRNA_ROC_5$AUC[2]),")"))
lines(mRNA_ROC_5$FP, mRNA_ROC_5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[2],lwd = 2)
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()

##PH test for the proportional hazards assumption of a Cox regression===================================

gene_cox_all<-coxph(obs.surv~., data=clinical_expr[,-c(1:3,ncol(clinical_expr))])
test.ph <- cox.zph(gene_cox_all)
ggph<-ggcoxzph(test.ph,point.col = "red", point.size = 1)
pdf(paste(parameters$pdf,parameters$PHtest,".pdf",sep = ""),onefile = FALSE, width=14,height=12)
print(ggph)
dev.off()

##Diagnostic Plots for Cox Proportional Hazards Model===================================================

gene_cox_all<-coxph(obs.surv~., data=clinical_expr[,-c(1:3,ncol(clinical_expr))])
ggcoxdiag<-ggcoxdiagnostics(gene_cox_all,type =  "dfbeta",linear.predictions = FALSE,ggtheme = theme_bw())
pdf(paste(parameters$pdf,parameters$coxdiagnostics,".pdf",sep = ""),onefile = FALSE, width=14,height=12)
print(ggcoxdiag)
dev.off()

##=====================================================================
print(sprintf('Complete'))


