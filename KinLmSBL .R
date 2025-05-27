#The code below written by Zhenghan Wu,Gang Han,and Mingzhi Cai
#Dependent R/Rcpp Packages

# File:KinLmSBL.R
# ------------------
# 
#
#
#
#
#
#
#


# ---------------------------------------------Initialization----------------------------------------------------------------------------------------------------- 

library(rMVP)
library(data.table)
library(doParallel)
library(RcppArmadillo)
library(Rcpp)
library(RcppEigen)

Rcpp::sourceCpp("sbl3.cpp")# form Zhenghan Wu
Rcpp::sourceCpp("lmm_diago_fit.cpp")# from Package gaston
Rcpp::sourceCpp("kinCorrect.cpp")# from FASTmrEMMA (Wen et al.,2018,2020)
source("emma_REMLE.r", encoding = "UTF-8")# from FASTmrEMMA (Wen et al.,2018,2020)

# fileGen <- "snp.all.csv"
fileGen <- "snp.all_plant_height.csv"
filePhe <- "rice01.raw"
fileKin <- "kk.csv"
filePs <- "fixedpc.csv"
fileMaf <- "MAF_check.frq"
dir = getwd()

KinLmSBL <- function(dir = NULL, fileGen = NULL, filePhe = NULL, fileKin = NULL, filePS = NULL, fileMAF = NULL) { 
  normalized_Kin <- normalizePath(fileKin, mustWork = FALSE)
  if (!file.exists(normalized_Kin)) {
    stop("fileKin dose not exist: ", normalized_Kin)
  }
  
  k<-read.table(file = normalized_Kin,sep = ",",header = TRUE)
  K<-matrix(unlist(k),nrow = dim(k)[1],ncol = dim(k)[2])
  
  #phenotype
  normalized_Phe <- normalizePath(filePhe, mustWork = FALSE)
  if (!file.exists(normalized_Phe)) {
    stop("filePhe dose not exist: ", normalized_Phe)
  }
  allData3 <- fread(file = normalized_Phe, showProgress = FALSE)
  dim(allData3)
  snpName<-colnames(allData3[1,7:dim(allData3)[2]])
  snpName2<-substr(snpName,1,nchar(snpName)-2)
  pheno<-allData3[,6]
  yo<-pheno
  if(is.numeric(yo)==FALSE){
    yo<-as.numeric(unlist(yo)) 
  }
  
  #Genotype  
  normalized_Gen <- normalizePath(fileGen, mustWork = FALSE)
  if (!file.exists(normalized_Gen)) {
    stop("fileGen dose not exist: ", normalized_Gen)
  }
  x<-read.table(file = normalized_Gen,sep = ",",header = TRUE)
  #SnpMaf
  normalized_MAF <- normalizePath(fileMaf, mustWork = FALSE)
  if (!file.exists(normalized_MAF)) {
    stop("fileMAF dose not exist: ", normalized_MAF)
  }
  snpmaf<-read.table(file = normalized_MAF,quote = "",header = TRUE)
  
  gen<-t(x[,-c(1,2,3,4)])
  pos<-x[,3]
  n1 <- length(yo)
  x1<-rep(1,times=n1)
  w0<-as.matrix(x1)
  z<-gen
  
  #If No Q matirx
    Y <- yo
  #Q matirx
    pcs<-read.table(file = normalized_Ps,sep = ",",header = TRUE)
    cvs<-matrix(unlist(pcs[,-1]),nrow = dim(pcs)[1],ncol = dim(pcs)[2]-1)
    
    system.time(pKinLm0<-RcppEigen::fastLmPure(y=yo,X=cvs)) #Check residual output
    Y=yo-cvs %*% pKinLm0$coefficients
  
  # -------------------------------------------- first stage: correct for polygenic background ----------------------------------

  qq<-eigen(K,symmetric=T)
  
  timeEmRemle<-system.time(remle2<-emma_REMLE(Y=Y, X=w0,qq)) #lmm_diago_fit.cpp and emma_REMLE.r
  timeRemle1<-system.time(remle1_B1<-emma_maineffects_B1(Z=matrix(),K,remle2$delta,complete=TRUE)) # kinCorrect.cpp
  
  C2<-remle1_B1$mC
  Y_c <- C2%*%Y
  t01<-proc.time()
  G_c <- C2%*%z
  t02<-proc.time()
  tcorrecttime<-(t02-t01)[3]
  # -------------------------------------------- first stage: linear regression scan ----------------------------------
  
  yName<-matrix(1:nrow(Y_c),ncol = 1)
  Y2<-cbind(yName,Y_c)
  timeBig<-system.time(G_c01<-as.big.matrix(t(G_c))) 
  timeMvp<-system.time(pKinLm01<-MVP.GLM(Y2, G_c01, cpu=1, mrk_bycol=FALSE))
  timerMVP<-timeBig+timeMvp
  head(pKinLm01)
  
  
  pKinLm<-pKinLm01[,3]
  names(pKinLm)<-snpName2  
  tax<-!is.na(pKinLm)
  pKinLm<-pKinLm[tax]
  m<-dim(gen)[2]
  seqloc<-c(1:m)
  selectLoc0.01<-seqloc[pKinLm<0.01]
  snpSelect0.01<-snpName2[selectLoc0.01]
  ps<-pos
  snpPsSelect0.01<-ps[selectLoc0.01]
  
  gen0.01<-gen[,selectLoc0.01]
  # ----------------------------------------- second stage: sparse Bayesian learning -------------------------------------
  sbltime0.01<-system.time(fit0.01<-sblgwas(x=w0, y=Y, z=gen0.01, t=-1)) 
  my.blup0.01<-fit0.01$blup
  my.parm0.01 <- fit0.01$parm ##
  sbl_pvalue <- my.blup0.01[,4]
  pvalue_new<-p.adjust(sbl_pvalue, method="BH", n=length(sbl_pvalue)) #False Discovery Rate correction
  ids_marked <- which(pvalue_new<0.05)
  
  var_Y <- var(Y)
  m2<-my.blup0.01[pvalue_new<0.05,]
  qtnIdKinlmsbl0.01<-snpSelect0.01[pvalue_new<0.05] 
  qtnPsKinlmsbl0.01<-snpPsSelect0.01[pvalue_new<0.05]
  beta <- my.parm0.01[ids_marked, 4]
  gamma <- my.blup0.01[ids_marked, 1]
  MAF <- snpmaf[ids_marked, 5]
  pvalue_final <- pvalue_new[pvalue_new<0.05]
  # ------------------------------------------------------ output -------------------------------------------------
  output_file <- file.path(dir, "KinLmSBL.result.csv")
  if (length(qtnIdKinlmsbl0.01) == 0) {
    print("Find no snp!")
    result <- data.frame(Chr = c(), Snp =  c(), Position =  c(), MAF =  c(), r2 =  c(), beta =  c(), pvalue =  c())
    write.csv(file = output_file, result)
  } else {
    Chr = x$chr[which(pvalue_new<0.05)]
    Snp = qtnIdKinlmsbl0.01
    Position = qtnPsKinlmsbl0.01
    MAF = round(MAF, 4)
    beta = round(gamma, 4)
    pvalue = format(pvalue_final, scientific = TRUE, digits = 4)
    r2 = round(fit0.01$blup$vg[pvalue_new<0.05] / rep(var_Y, length(fit0.01$blup$vg[pvalue_new<0.05])) * 100, 4)
    result <- data.frame(Chr = Chr, Snp = Snp, Position = Position, MAF = MAF, 'r2(%)' = r2, beta = beta, pvalue = pvalue, check.names = FALSE) #p_value,effect,MAF?,pve,MrMLM,3v?
    write.csv(file = output_file, result)
    return(result)
  }  
}
KinLmSBL(dir = dir, fileGen = fileGen, filePhe = filePhe, fileKin = fileKin, filePS = filePs, fileMAF = fileMaf  )

