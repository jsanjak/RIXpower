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
source("genomat.utils.R")
source("myhaploLMM.R")

#####
##### Variable arguments
h <- 2 # Number of alleles
N <-768 # Number of individuals
H2 <- 0.1 # Heritability explained by the locus
reps <-3 # Replicate samples
domi <- "add" #Model for genetic dominance

set.seed(7)
#####
##### Hard coded parameters
f <- 3 # Number of founders
p <- 1/f# Equal expected allele frequencies
r <- 3 # Number of loci analyzed
ntrials <- 1 # Number of simulation trials
ksites <- 1 # Number of sites per locus
gamma <- 0.1
eps <- 10^(-6) # The gamma squared parameter from Burger defining the mutational variance
#####
output <- data.frame()

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

sum(rowSums(GENO_RIX_1)==(2*r))==nrow(GENO_RIX_1)
sum(rowSums(GENO_RIX_2)==(2*r))==nrow(GENO_RIX_2)

#Calculate phenotypes
#The approach here is to look at a focal locus on an otherwise additive polygenic background
#The genetic architechture of the trait will be controlle by the following parameters
#The total heritability of the trait
#The heritability due to the focal locus
#The proportion of the non-focal loci which are causal

H2_trait <- 0.8
H2_focus <- 0.1
prop_causal <- 0.5
poly_causal <- r*prop_causal - 1
H2_polygenic <- H2_trait - H2_focus
VE <- 1 - H2_trait

#focal locus is a middle one
focus <- r/2
causal <- sample(sample(c(rep(1,poly_causal),rep(0,r-poly_causal-1))))
causal <- c(causal[1:(focus-1)],0,causal[(focus):(r-1)])

a_focus_uniq <- c(0,rgamma(h-1,shape=ksites,scale=gamma)) # Unique allelic effects
a_focus <- c(a_focus_uniq,sample(a_focus_uniq,f-h,replace=T))
#Genetic value at focal locus
G_FOCUS_RIX_1 <- GBR_Gvalue(GENO_RIX_1[,(f*(focus-1)+1):(f*focus)],f,1,a_focus,domi)
G_FOCUS_RIX_2 <- GBR_Gvalue(GENO_RIX_2[,(f*(focus-1)+1):(f*focus)],f,1,a_focus,domi)


#algebra to figure our what the polygenic VG should be
poly_VG_1 <-var(G_FOCUS_RIX_1)*(H2_trait +  (H2_trait*(1-H2_focus)/(H2_focus)) - 1)
poly_VG_2 <-var(G_FOCUS_RIX_2)*(H2_trait +  (H2_trait*(1-H2_focus)/(H2_focus)) - 1)


get_poly_effects = function(causal,f,poly_VG){
  
  if(causal==1) {
    effects <- rnorm(f,0,sqrt(poly_VG))
  } else {
    effects <- rep(0,f) 
  }
  return(effects)
  
}


a_poly_1 <- unmatrix(sapply(causal,function(x) get_poly_effects(x,f,poly_VG_1/(2*poly_causal))))
a_poly_2 <- unmatrix(sapply(causal,function(x) get_poly_effects(x,f,poly_VG_2/(2*poly_causal))))

G_POLY_RIX_1 <- rowSums(t(a_poly_1 * t(GENO_RIX_1)))
G_POLY_RIX_2 <- rowSums(t(a_poly_2 * t(GENO_RIX_2)))

G_RIX_1 <- G_POLY_RIX_1 + G_FOCUS_RIX_1
G_RIX_2 <- G_POLY_RIX_2 + G_FOCUS_RIX_2

VE_1 <- var(G_RIX_1)*(1-H2_trait)/H2_trait
VE_2 <- var(G_RIX_2)*(1-H2_trait)/H2_trait

E_1 <- sapply(rep(0,reps),function(x)rnorm(length(G_RIX_1),0,sqrt(VE_1)))
E_2 <- sapply(rep(0,reps),function(x)rnorm(length(G_RIX_2),0,sqrt(VE_2)))

P_RIX_1 <- G_RIX_1 + rowMeans(E_1)
P_RIX_2 <- G_RIX_2 + rowMeans(E_2)



#Fit LMM
# Y = X * a + Z_A * A + Z_ALPHA * ALPHA + e
# X is the genotype matrix and a is a vector of fixed genotype effects
# Z_ALPHA is the ALPHA RIL type and ALPHA is the randome effect of the ALPHA RIL
# Z_A is the A RIL type and A is the random effect of the A RIL
library(lme4)

#one choice is to model each parent as a dummy variable
#MODEL_MAT_1 <- data.frame(P_RIX_1,GENO_RIX_1[,(f*(focus-1)+1):(f*focus)],Z_RIX_1$A,Z_RIX_1$ALPHA)
#MODEL_MAT_2 <- data.frame(P_RIX_2,GENO_RIX_2[,(f*(focus-1)+1):(f*focus)],Z_RIX_2$A,Z_RIX_2$ALPHA)


#the other approach is to model the A and Alpha as the two random effects with a lot of levels
MODEL_MAT_1 <- data.frame(P_RIX_1,GENO_RIX_1[,(f*(focus-1)+1):(f*focus)],RIX_1)
MODEL_MAT_2 <- data.frame(P_RIX_2,GENO_RIX_2[,(f*(focus-1)+1):(f*focus)],RIX_2)

haplofreq_focus_1 <- colSums(GENO_RIX_1[,(f*(focus-1)+1):(f*focus)])/(2*nrow(GENO_RIX_1))
haplofreq_focus_2 <- colSums(GENO_RIX_1[,(f*(focus-1)+1):(f*focus)])/(2*nrow(GENO_RIX_2))

colnames(MODEL_MAT_1) <- c("Y",sapply(seq(1:(ncol(MODEL_MAT_1)-1)),function(x) paste("X",x,sep="")))
colnames(MODEL_MAT_2) <- c("Y",sapply(seq(1:(ncol(MODEL_MAT_2)-1)),function(x) paste("X",x,sep="")))

names_1 <- colnames(MODEL_MAT_1) 
names_2 <- colnames(MODEL_MAT_2)

LMM_formula_1 <- as.formula(paste(names_1[1]," ~ ",
                              paste(c(names_1[3:(f+1)],
                                      sapply(names_1[(f+2):length(names_1)],
                                             function(x) paste("(1 | ",x,")",sep=""))),
                                      collapse=" + "),sep=""))
LMM_null_formula_1 <- as.formula(paste(names_1[1]," ~ ",
                                  paste(sapply(names_1[(f+2):length(names_1)],
                                                 function(x) paste("(1 | ",x,")",sep="")),
                                        collapse=" + "),sep=""))

simple_formula <- as.formula(paste(names_1[1]," ~ ",
                                   paste(names_1[3:(f+1)],
                                         collapse=" + "),sep=""))


LMM_1 <- lmer(formula=LMM_formula_1,data=MODEL_MAT_1)
LMM_1_null <- lmer(formula=LMM_null_formula_1,data=MODEL_MAT_1)

LOD <- (-2*(logLik(LMM_1_null) - logLik(LMM_1)))/(2*log(10))

GLM_1 <- glm(formula=simple_formula,data=MODEL_MAT_1,family="gaussian")
GLM_NULL_1  <-glm(formula="Y ~ 1",data=MODEL_MAT_1,family="gaussian")
source("~/Documents/Research/yeast_cross/Power/RIX/myhaploLMM.R")
source("~/Documents/Research/yeast_cross/Power/RIX/myhaplopen.R")

LMM_1.summ <- summary(LMM_1)
LMM_1.summ$fitted.values <- fitted(LMM_1)
LMM_1.summ$family <- "gaussian"
LMM_1.summ$y <- MODEL_MAT_1$Y
hapdata <- cbind(HAP_RIL_A[,focus],HAP_RIL_ALPHA[,focus]) - ((focus - 1) * (f))
lambdas <- seq(0.01,5.01,0.1)
eps=1e-6

CAS_LMM_1 <- my.haplo.pen.lmm(LMM_1.summ, data=MODEL_MAT_1[,c(1,3:(f+1))] ,
                              haplofreq=haplofreq_focus_1, hapdata=hapdata,hap.base = 1,
                              lambdas = lambdas, eps = eps)

CAS_GLM <- my.haplo.pen.glm.1(formula = as.formula(Y~.), data=MODEL_MAT_1[,c(1,3:(f+1))], family = "gaussian" , 
                              haplofreq = haplofreq_focus_1, hap.base = 1, hapdata = hapdata,
                              lambdas = lambdas, eps = eps)


beta.t <- round(CAS_LMM_1$coefficients[-1],max(round(-log10(CAS_LMM_1$init.fit$coefficients[,2]))))
alleles <- 1+ length(unique(beta.t[beta.t!=0])) #round(beta[-1],-log10(eps)-1)))

par(mfrow=c(1,2))
plot(LMM_1.summ$coefficients[,1],col="red",pch=19,main="LMM")
points(CAS_LMM_1$coefficients,col="black",pch=19)

plot(GLM_1$coefficients,col="red",pch=19,main="GLM")
points(CAS_GLM$coefficients,col="black",pch=19)



P_dat_null <- prod(dnorm(MODEL_MAT_1$Y, mean = mean(MODEL_MAT_1$Y), sd = sqrt(var(MODEL_MAT_1$Y))))
P_dat_LMM <- prod(haplo.stats:::dglm.fit(CAS_LMM_1))
LOD <- log10(P_dat_null/P_dat_LMM)
#TODO
#Impute the RIX dominance matrices for each locus

