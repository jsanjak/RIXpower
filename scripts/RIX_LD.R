##### haploDMR.R simulates QTL mapping and haplotype effect estimation in a multiparental cross
##### and uses DMR method to collapse levels
# last edited: 02/07/2017 
# Jaleal Sanjak
#####
#####
rm(list=ls(all=TRUE))
require(plyr)
require(MASS)
require(haplo.stats)
library(magic)
library(gdata)
library(lme4)
library(ggplot2)
source("genomat.utils.R")
source("myhaploLMM.R")

#####
##### Variable arguments
h <- 2 # Number of alleles
Nind <-c(24,48,72,96,192,384,768) # Number of individuals
H2 <- 0.1 # Heritability explained by the locus
reps <-3 # Replicate samples
domi <- "add" #Model for genetic dominance

set.seed(7)
#####
##### Hard coded parameters
f <- 18# Number of founders
p <- 1/f# Equal expected allele frequencies
r <- 10 # Number of loci analyzed
ntrials <- 1 # Number of simulation trials
ksites <- 1 # Number of sites per locus
gamma <- 0.1
eps <- 10^(-6) # The gamma squared parameter from Burger defining the mutational variance
#####
output <- NULL
for(N in Nind){
  #RIL genotypes as founder haplotype number N by r matrix
  #We have A's and Alphas
  HAP_RIL_A = as.matrix(ldply(1:N,.fun=function(x) sample(seq(1,f),size=r,prob=rep(p,f),replace=TRUE)))
  HAP_RIL_ALPHA =  as.matrix(ldply(1:N,.fun=function(x) sample(seq(1,f),size=r,prob=rep(p,f),replace=TRUE)))
  #Model matrix N by nalpha
  HAP_RIL_A = t( ((1:r -1)*(f))+ t(HAP_RIL_A))
  HAP_RIL_ALPHA = t(((1:r -1)*(f))+ t(HAP_RIL_ALPHA))
  
  #RIL genotypes
  GENO_RIL_A = matrix(0,nrow=N,ncol=r*f)
  GENO_RIL_ALPHA = matrix(0,nrow=N,ncol=r*f)
  
  for (i in 1:N ){
    
    GENO_RIL_A[i,HAP_RIL_A[i,]] = GENO_RIL_A[i,HAP_RIL_A[i,]] + 1
    GENO_RIL_ALPHA[i,HAP_RIL_ALPHA[i,]] = GENO_RIL_ALPHA[i,HAP_RIL_ALPHA[i,]] + 1
    
  }
  sum(rowSums(GENO_RIL_A)==r)==N
  sum(rowSums(GENO_RIL_ALPHA)==r)==N
  
  #RIX genotypes will depend on the design
  #model 1  96 different 8 by 8 RIL intercross  96*(8*8) = 6144 RIX's
  #model 2  64 different 16 RIL intercross  64*(12*12) = 9216 RIX's
  base_1 <-  matrix(1,8,8)
  base_2 <-  matrix(1,12,12)
  #base_1[upper.tri(base_1)] <- 1
  #base_2[upper.tri(base_2)] <- 1
  
  #Design matrices are block upper triangular
  Design_1 <- base_1
  Design_2 <- base_2
  
  for ( i in 1:(N/nrow(base_1)-1)){
    if(i<= (N/nrow(base_2)-1)) {
      Design_2 <- adiag(Design_2,base_2)
      
    }
    Design_1 <- adiag(Design_1,base_1)
    
  }
  
  #RIX genotypes in condensed format
  RIX_1 <- which(Design_1 == 1, arr.ind = T)
  RIX_2 <- which(Design_2 == 1, arr.ind = T)
  
  #Z is a list matrices that are M by N, where M is the number of RIX and N is the number of RILs
  #There will be an A and Alpha matrix to represent the A and Alpha RIL of each RIX
  fill_Z_RIX = function(RIX,N_RIL,P){
    Z_RIX <- rep(0,N_RIL)
    if(P>length(RIX)){
      break
    }
    Z_RIX[RIX[P]] <- 1
    return(Z_RIX)
  }
  
  Z_RIX_1 <- list(A= t(apply(RIX_1,1,function(x){fill_Z_RIX(x,N,1)})),
                  ALPHA = t(apply(RIX_1,1,function(x){fill_Z_RIX(x,N,2)})))
  
  Z_RIX_2 <- list(A= t(apply(RIX_2,1,function(x){fill_Z_RIX(x,N,1)})),
                  ALPHA = t(apply(RIX_2,1,function(x){fill_Z_RIX(x,N,2)})))
  
  #make sure each RIX has exactly one A and one Alpha RIL parent
  sum(rowSums(Z_RIX_1$A)==1) == nrow(RIX_1)
  sum(rowSums(Z_RIX_2$A)==1) == nrow(RIX_2)
  
  #Impute the multi locus genotypes of the RIX's in terms of the founder haplotypes based on the RIL genotypes
  fill_G_RIX = function(Z_RIX,GENO_RIL){
    RILS <- which(Z_RIX==1)
    G_RIX <- GENO_RIL[RILS[1],]
    return(G_RIX)
  }
  
  GENO_RIX_A_1 <- t(apply(Z_RIX_1$A,1,function(x){fill_G_RIX(x,GENO_RIL_A)}))
  GENO_RIX_ALPHA_1 <- t(apply(Z_RIX_1$ALPHA,1,function(x){fill_G_RIX(x,GENO_RIL_ALPHA)}))
  
  GENO_RIX_1 <- GENO_RIX_A_1 + GENO_RIX_ALPHA_1
  
  GENO_RIX_A_2 <- t(apply(Z_RIX_2$A,1,function(x){fill_G_RIX(x,GENO_RIL_A)}))
  GENO_RIX_ALPHA_2 <- t(apply(Z_RIX_2$ALPHA,1,function(x){fill_G_RIX(x,GENO_RIL_ALPHA)}))
  
  GENO_RIX_2 <- GENO_RIX_A_2 + GENO_RIX_ALPHA_2
  
  NRIX_1<-nrow(GENO_RIX_1)
  NRIX_2<-nrow(GENO_RIX_2)
  
  allele_freq_1 <- colSums(GENO_RIX_1)/(2*NRIX_1)
  allele_freq_2 <- colSums(GENO_RIX_2)/(2*NRIX_2)
  
  #R squared calculation
  r_sq_hap_1 <- R_sq_hap(r,f,allele_freq_1,GENO_RIX_A_1,GENO_RIX_ALPHA_1,NRIX_1)
  r_sq_hap_2 <- R_sq_hap(r,f,allele_freq_2,GENO_RIX_A_2,GENO_RIX_ALPHA_2,NRIX_2)
  
  r_sq_data <- rbind(cbind(rep(N,r*(r-1)/2), rep(1,r*(r-1)/2),r_sq_hap_1),
                     cbind(rep(N,r*(r-1)/2), rep(2,r*(r-1)/2),r_sq_hap_2))
  
  
  output <- rbind(output,r_sq_data)
}
output <- data.table(setNames(output,c("N","Design","rsq")))

pdf("../plots/RIX_LD_DECAY.pdf",width=9,height=8,pointsize =18 )
p<- ggplot(filter(output,Design==1),aes(N,rsq,group=N)) + geom_boxplot() +theme_bw()
p <- p + xlab("Number of RILS in the RIX") + ylab(expression(paste(r^2)))
p + theme(axis.text = element_text(size=20),
          axis.title= element_text(size = 22),
          legend.title= element_text(size = 20),
          legend.text= element_text(size = 18),
          plot.title = element_text(size=24, hjust = 0),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"))
dev.off()

