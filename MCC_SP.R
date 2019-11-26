library(dplyr)
library(energy)
library(acepack)
#---------Before running the following code, please run the following 3 source files and ensure that the required variables have been arranged as required.
source("KDE.r")
source("dijkstra.r")
source("kshortestpath.r")
#------------------The following comment lines indicate that the variable can be obtained by reading the input variables in the format of csv
#------------------Readers can also get the input variables through R.data for this function
# ADs_gene <- read.csv("ADnetwork.csv",header = T,stringsAsFactors=F)
# sig_gene <- read.csv("sig_gene.csv",stringsAsFactors=F)
# sigma <- read.csv("network.csv",header=TRUE)
# family_xuhao <- read.csv("family_xuhao.csv",header = T,stringsAsFactors=F)
MCC_SP <- function(gene,cg_table1,family_xuhao,sigma,sig_gene, path_test=10000){

  all_gene <- gene
  apoe_type <- ifelse(cg_table1$apoe_genotype==24|cg_table1$apoe_genotype==34|cg_table1$apoe_genotype==44,1,0)
  all_gene <- cbind(all_gene,cg_table1$ending)
  all_gene <- cbind(all_gene,apoe_type)
  colnames(all_gene) <- c(colnames(gene),'ending','apoe_type')
  all_gene <- as.data.frame(all_gene)

  #-------------------------------------------------------------
  mean_data <-gene[,sig_gene$gene]
  mean_data2 <- as.data.frame(mean_data)
  xuhao <- family_xuhao$gene
  library(dplyr)
  aa <- mean_data2 %>% select(xuhao)
  sort_data <- as.matrix(aa)
  #The variable 'ex' indicates whether the APOE genotype is mutated
  ex <- ifelse(cg_table1$apoe_genotype==24|cg_table1$apoe_genotype==34|cg_table1$apoe_genotype==44,1,0)
  mydata <- cbind(ex,sort_data,cg_table1$ending)
  colnames(mydata)[1] <- 'ex'
  colnames(mydata)[ncol(mydata)] <- 'Y'
  #-------------------------------------------------
  sigma <- as.matrix(sigma)
  #The variable path_test is a temporary variable set to determine how many paths there are in the network. Generally, it is considered that it is sufficient to set 10000 in a small network (about 40 points or less).
  start <- 1
  destination <- ncol(mydata)
  #---------------------------------------------------------
  corpearson<-matrix(0,nrow(sigma),ncol(sigma))
  for(i in 1:nrow(sigma)){
    for(j in 1:nrow(sigma)){
      if(sigma[i,j]!=0){
        data.x<-mydata[,i]
        data.y<-mydata[,j]
        corpearson[i,j]<-cor(data.x,data.y,method = "pearson")
      }
    }
  }
  Ccorpearson<-abs(corpearson)
  for(i in 1:nrow(corpearson)){
    for(j in 1:nrow(corpearson)){
      Ccorpearson[i,j]<-log(1/Ccorpearson[i,j])
    }
  }
  
  for(i in 1:nrow(corpearson)){
    for(j in 1:nrow(corpearson)){
      if(Ccorpearson[i,j]==0){
        Ccorpearson[i,j]<-Inf
      }
    }
  }
  #The list 'resultpe' is the result of pearson_SP, and uses this result as a reference to sort the paths found by other methods. Therefore, the order of the resultspe's path is a reference, which is 1-33.
  resultpe<-kshortestpath( Ccorpearson,start,destination,path_test)
  paths <- length(resultpe[[1]])
  #spearman
  corspearman<-matrix(0,nrow(sigma),ncol(sigma))
  for(i in 1:nrow(sigma)){
    for(j in 1:nrow(sigma)){
      if(sigma[i,j]!=0){
        data.x<-mydata[,i]
        data.y<-mydata[,j]
        corspearman[i,j]<-cor(data.x,data.y,method = "spearman")
      }
    }
  }
  Ccorspearman<-abs(corspearman)
  for(i in 1:nrow(corspearman)){
    for(j in 1:nrow(corspearman)){
      Ccorspearman[i,j]<-log(1/Ccorspearman[i,j])
    }
  }
  
  for(i in 1:nrow(corspearman)){
    for(j in 1:nrow(corspearman)){
      if(Ccorspearman[i,j]==0){
        Ccorspearman[i,j]<-Inf
      }
    }
  }
  rresultsp<-kshortestpath( Ccorspearman,start,destination,paths)
  
  
  library(energy)
  cordistance<-matrix(0,nrow(sigma),ncol(sigma))
  for(i in 1:nrow(sigma)){
    for(j in 1:nrow(sigma)){
      if(sigma[i,j]!=0){
        data.x<-mydata[,i]
        data.y<-mydata[,j]
        cordistance[i,j]<-dcor(data.x,data.y)
      }
    }
  }
  Ccordistance<-abs(cordistance)
  for(i in 1:nrow(cordistance)){
    for(j in 1:nrow(cordistance)){
      Ccordistance[i,j]<-log(1/Ccordistance[i,j])
    }
  }
  
  for(i in 1:nrow(cordistance)){
    for(j in 1:nrow(cordistance)){
      if(Ccordistance[i,j]==0){
        Ccordistance[i,j]<-Inf
      }
    }
  }
  rresultdi<-kshortestpath( Ccordistance,start,destination,paths)
  

  cormcc<-matrix(0,nrow(sigma),ncol(sigma))
  for(i in 1:nrow(sigma)){
    for(j in 1:nrow(sigma)){
      if(sigma[i,j]!=0){
        data.x<-mydata[,i]
        data.y<-mydata[,j]
        argmax = ace(data.x, data.y)
        cormcc[i,j]<-cor(argmax$tx, argmax$ty)
      }
    }
  }
  Ccormcc<-abs(cormcc)
  for(i in 1:nrow(cormcc)){
    for(j in 1:nrow(cormcc)){
      Ccormcc[i,j]<-log(1/Ccormcc[i,j])
    }
  }
  
  for(i in 1:nrow(cormcc)){
    for(j in 1:nrow(cormcc)){
      if(Ccormcc[i,j]==0){
        Ccormcc[i,j]<-Inf
      }
    }
  }
  rresultmcc<-kshortestpath( Ccormcc,start,destination,paths)
  library(minerva)
  
  cormic<-matrix(0,nrow(sigma),ncol(sigma))
  for(i in 1:nrow(sigma)){
    for(j in 1:nrow(sigma)){
      if(sigma[i,j]!=0){
        data.x<-mydata[,i]
        data.y<-mydata[,j]
        cormic[i,j]<-mine(data.x,data.y,alpha=0.3)$MIC
      }
    }
  }
  Ccormic<-abs(cormic)
  for(i in 1:nrow(cormic)){
    for(j in 1:nrow(cormic)){
      Ccormic[i,j]<-log(1/Ccormic[i,j])
    }
  }
  
  for(i in 1:nrow(cormic)){
    for(j in 1:nrow(cormic)){
      if(Ccormic[i,j]==0){
        Ccormic[i,j]<-Inf
      }
    }
  }
  rresultmic<-kshortestpath( Ccormic,start,destination,paths)
  #----------------------
  jieguo <- matrix()
  sigmaI<-sigma
  sigmaI[lower.tri(sigmaI,diag=T)]=0
  edge <- which(sigmaI!=0,arr.ind=T)
  
  cormi<-diag(rep(1,nrow(sigma)))
  jieguo<-apply(edge,DataFit=mydata,1,denPre2D)
  for(i in 1:nrow(edge)){
    cormi[edge[i,1],edge[i,2]]<-mean(jieguo[,i])
  }
  
  Ccormi<-abs(cormi)
  for(i in 1:nrow(cormi)){
    for(j in 1:nrow(cormi)){
      Ccormi[i,j]<-log(1/Ccormi[i,j])
    }
  }
  
  for(i in 1:nrow(cormi)){
    for(j in 1:nrow(cormi)){
      if(Ccormi[i,j]==0){
        Ccormi[i,j]<-Inf
      }
    }
  }
  rresultmi<-kshortestpath( Ccormi,start,destination,paths)
  #----------------------------------------------------The following code gets path sorting for various methods
  allpaths <- resultpe[[1]]
  names(allpaths) <- paste0('path_',1:paths)
  #pearson
  pathorderpk <- c()
  for(k in 1:paths){
    for(p in 1:paths){
      if(all(allpaths[[k]]==rresultpearson$shortestpaths[[p]])){
        order <- p
      }
    }
    pathorderpk <- c(pathorderpk,order)
  }
  #spearman
  pathordersk <- c()
  for(k in 1:paths){
    for(p in 1:paths){
      if(all(allpaths[[k]]==rresultsp$shortestpaths[[p]])){
        order <- p
      }
    }
    pathordersk <- c(pathordersk,order)
  }
  #di
  pathorderdk <- c()
  for(k in 1:paths){
    for(p in 1:paths){
      if(all(allpaths[[k]]==rresultdi$shortestpaths[[p]])){
        order <- p
      }
    }
    pathorderdk <- c(pathorderdk,order)
  }
  #mic
  pathordermick <- c()
  for(k in 1:paths){
    for(p in 1:paths){
      if(all(allpaths[[k]]==rresultmic$shortestpaths[[p]])){
        order <- p
      }
    }
    pathordermick <- c(pathordermick,order)
  }
  #mi
  pathordermik <- c()
  for(k in 1:paths){
    for(p in 1:paths){
      if(all(allpaths[[k]]==rresultmi$shortestpaths[[p]])){
        order <- p
      }
    }
    pathordermik <- c(pathordermik,order)
  }
  #mcc
  pathordermcck <- c()
  for(k in 1:paths){
    for(p in 1:paths){
      if(all(allpaths[[k]]==rresultmcc$shortestpaths[[p]])){
        order <- p
      }
    }
    pathordermcck <- c(pathordermcck,order)
    
  }
  
  Spearman_SP <- order(pathordersk)
  DC_SP <- order(pathorderdk)
  MIC_SP <- order(pathordermick)
  MI_SP <- order(pathordermik)
  MCC_SP <- order(pathordermcck)
  result_end <- list(Spearman_SP,DC_SP,MIC_SP,MI_SP,MCC_SP)
  names(result_end) <- c('Spearman_SP','DC_SP','MIC_SP','MI_SP','MCC_SP')
  return(result_end)
}

