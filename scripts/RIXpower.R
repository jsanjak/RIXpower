##### haploDMR.R simulates QTL mapping and haplotype effect estimation in a multiparental cross
##### and uses DMR method to collapse levels
# last edited: 02/07/2017
# Jaleal Sanjak
#####
#####
rm(list=ls(all=TRUE))
require(plyr)
require(MASS)
library(magic)
library(gdata)
library(lme4)
source("genomat.utils.R")
source("helpers.R")
source("casanova.R")
#####
##### Variable arguments
h <- 2 # Number of alleles
N <- 768 # Number of effective diploid individuals to be sequences
H2 <- 0.1 # Heritability explained by the locus
reps <-10 # Replicate samples
domi <- "add" #Model for genetic dominance
f <- 3 
r <- 10
ofile="test.txt"
ntrials <- 1
set.seed(9)

H2_trait <- 0.8
H2_focus <- 0.1
prop_causal <- 0.5
poly_causal <- r*prop_causal - 1
H2_polygenic <- H2_trait - H2_focus
VE <- 1 - H2_trait

#####
##### Hard coded parameters
# Number of founders
p <- 1/f# Equal expected allele frequencies
 # Number of loci analyzed
 # Number of simulation trials
ksites <- 1 # Number of sites per locus
gamma <- 1
eps <- 10^(-6) # The gamma squared parameter from Burger defining the mutational variance
lambda_pt <- seq(0.0,10,0.01)

#Information for Calculating phenotypes
#The approach here is to look at a focal locus on an otherwise additive polygenic background
#The genetic architechture of the trait will be controlle by the following parameters
#The total heritability of the trait
#The heritability due to the focal locus
#The proportion of the non-focal loci which are causal
#focal locus is a middle one
focus <- r/2
causal <- sample(sample(c(rep(1,poly_causal),rep(0,r-poly_causal-1))))
causal <- c(causal[1:(focus-1)],0,causal[(focus):(r-1)])
a_focus_uniq <- c(0,rgamma(h-1,shape=ksites,scale=gamma)) # Unique allelic effects
a_focus <- c(a_focus_uniq,sample(a_focus_uniq,f-h,replace=T))
#####
#prepare output file
write("nfounders reps H2_trait H2_locus nloci prop_causal design LOD_add LOD_dom LOD_full alleles lambda sigma_sq",ofile,ncol=13)

#
#RIL genotypes as founder haplotype number N by r matrix
#We have A's and Alphas
HAP_RIL_A = as.matrix(plyr:::ldply(1:N,.fun=function(x) sample(seq(1,f),size=r,prob=rep(p,f),replace=TRUE)))
HAP_RIL_ALPHA =  as.matrix(plyr:::ldply(1:N,.fun=function(x) sample(seq(1,f),size=r,prob=rep(p,f),replace=TRUE)))
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
#model 2  64 different 12 by 12 RIL intercross  64*(12*12) = 9216 RIX's
#model 3  76 different 8 by 12 RIL intercross 76*(8*12) = 7296 RIX's

base_1 <-  matrix(1,8,8)
base_2 <-  matrix(1,12,12)
base_3 <-  matrix(1,8,12)

Designs <- list(D1=base_1,
                D2=base_2,
                D3=base_3)


for ( i in 1:(N/nrow(base_1)-1)){
  if(i<= (N/nrow(base_2)-1)) {
    Designs$D2 <- adiag(Designs$D2,base_2)

  }
  if(i<= floor(N/mean(c(nrow(base_3),ncol(base_3))))-1) {
    base_3 <- t(base_3)
    Designs$D3 <- adiag(Designs$D3,base_3)
  }
  Designs$D1 <- adiag(Designs$D1,base_1)

}

#pdf("../plots/Design1.pdf",height=8,width=8)
#image(Design_1,asp=1,xlim=c(0,(8*4)/768),ylim=c(0,(8*4)/768),col=c("white","black"))
#dev.off()

#pdf("../plots/Design2.pdf",height=8,width=8)
#image(Design_2,asp=1,xlim=c(0,(12*3)/768),ylim=c(0,(12*3)/768),col=c("white","black"))
#dev.off()

#pdf("../plots/Design3.pdf",height=8,width=8)
#image(Design_3,asp=1,xlim=c(0,40/760),ylim=c(0,40/760),col=c("white","black"))
#dev.off()

#pdf("../plots/Design3full.pdf",height=8,width=8)
#image(Design_3,asp=1,col=c("white","black"))
#dev.off()


#A for loop over designs
#there are just a few so this isnt a big deal
for(di in 1:length(Designs)){
#di<-1
Design <- Designs[[di]]

RIX <- RIX_geno(Design)

#Z is a list matrices that are M by N, where M is the number of RIX and N is the number of RILs
#There will be an A and Alpha matrix to represent the A and Alpha RIL of each RIX

Z_RIX <- list(A= t(apply(RIX,1,function(x){fill_Z_RIX(x,N,1)})),
                ALPHA = t(apply(RIX,1,function(x){fill_Z_RIX(x,N,2)})))

#make sure each RIX has exactly one A and one Alpha RIL parent
sum(rowSums(Z_RIX$A)==1) == nrow(RIX)

GENO_RIX_A <- t(apply(Z_RIX$A,1,function(x){fill_G_RIX(x,GENO_RIL_A)}))
GENO_RIX_ALPHA <- t(apply(Z_RIX$ALPHA,1,function(x){fill_G_RIX(x,GENO_RIL_ALPHA)}))

GENO_RIX <- GENO_RIX_A + GENO_RIX_ALPHA

sum(rowSums(GENO_RIX)==(2*r))==nrow(GENO_RIX)

#Genetic value at focal locus
G_FOCUS_RIX <- GBR_Gvalue(GENO_RIX[,(f*(focus-1)+1):(f*focus)],f,1,a_focus,domi)

#algebra to figure our what the polygenic VG should be
poly_VG <-var(G_FOCUS_RIX)*(H2_trait +  (H2_trait*(1-H2_focus)/(H2_focus)) - 1)

a_poly <- unmatrix(sapply(causal,function(x) get_poly_effects(x,f,poly_VG/(2*poly_causal))))

G_POLY_RIX <- rowSums(t(a_poly * t(GENO_RIX)))

G_RIX <- G_POLY_RIX + G_FOCUS_RIX

VE <- var(G_RIX)*(1-H2_trait)/H2_trait

#Here we cheat in the model a bit by ignoring the repeated
#measurements as a random effect and just
#take the line mean.
E <- matrix(rnorm(length(G_RIX)*reps,0,sqrt(VE)),nrow=length(G_RIX),ncol=reps)

#P_RIX <- rep(G_RIX,reps)  + E
P_RIX <- G_RIX + rowMeans(E)

#CROSS <- rep(1:length(G_RIX),reps)
CROSS <- 1:length(G_RIX)
#Fit LMM
# Y = X_a * a + X_d * d + Z_A * A + Z_ALPHA * ALPHA  Z_CROSS * CROSS + e
# X_a is the genotype matrix and a is a vector of fixed additive genotype effects
# X_d is the dominance encoding matrix and d is a vector of fixed  dominance genotype effects
# Z_ALPHA is the ALPHA RIL type and ALPHA is the random effect of the ALPHA RIL
# Z_A is the A RIL type and A is the random effect of the A RIL
# Z_CROSS is the RIX type and CROSS is the random effect of the RIX -- the combine AxAlpha

DOM <- fill_d_zero(GENO_RIX[,(f*(focus-1)+1):(f*focus)],f,1)

MODEL_MAT <- data.frame(P_RIX,CROSS,
                        do.call("rbind", 
                                rep(list(data.frame(GENO_RIX[,(f*(focus-1)+1):(f*focus)],
                                                    DOM,
                                                    RIX)), 
                                    1)))

colnames(MODEL_MAT) <- c("Y","Cross",sapply(seq(1:(ncol(MODEL_MAT)-4)),function(x) paste("X",x,sep="")),"MAT_A","MAT_ALPHA")

modnames <- colnames(MODEL_MAT)

add_names <- modnames[c(3:(f+2))]
dom_names <- modnames[(f+3):(length(modnames)-2)]
RIX_names <- modnames[c((length(modnames)-1):length(modnames),2)]


LMM_formula_null <-as.formula( paste(modnames[1]," ~ 0 + ",
                                     paste(sapply(RIX_names[-3],
                                                  function(x) paste("(1 | ",x,")",sep="") ),
                                           collapse=" + ") ))

LMM_formula_add <- as.formula( paste(modnames[1]," ~ 0 + ",
                               paste(c(add_names,
                                      sapply(RIX_names[-3],
                                             function(x) paste("(1 | ",x,")",sep="") )),
                                      collapse=" + ") ))

LMM_formula_full <- as.formula( paste(modnames[1]," ~ 0 + ",
                                      paste(c(add_names,dom_names,
                                              sapply(RIX_names[-3],
                                                     function(x) paste("(1 | ",x,")",sep="") )),
                                            collapse=" + ") ))



LMM_null <- lmer(formula=LMM_formula_null,data=MODEL_MAT,REML=TRUE)
LMM_add <- lmer(formula=LMM_formula_add,data=MODEL_MAT,REML=TRUE)
LMM_full <- lmer(formula=LMM_formula_full,data=MODEL_MAT,REML=TRUE)


##################
#LOD scores
LOD_add <- (-2*(logLik(LMM_null) - logLik(LMM_add)))/(2*log(10))
LOD_dom  <- (-2*(logLik(LMM_add) - logLik(LMM_full)))/(2*log(10))
LOD_full <- (-2*(logLik(LMM_null) - logLik(LMM_full)))/(2*log(10))

#####
#CASANOVA 
x <- as.matrix(MODEL_MAT[add_names]) 
y <- as.matrix(MODEL_MAT[,1])

LAMBDA  <- getME(LMM_null,"Lambda")
sigma_sq <- getME(LMM_null,"sigma")^2
Z <- t(do.call(Matrix::rBind,getME(LMM_null,"Ztlist")))
bigSigma <- Z%*%tcrossprod(LAMBDA)%*%t(Z)

ZId <- sparseMatrix(1:nrow(Z),1:nrow(Z),x=rep(1,nrow(Z)))
beta.start <- fixef(LMM_add)

eig <- eigen(bigSigma)

U <- eig$vectors
d <- eig$values

DIinv <- diag(1/sqrt(d+1))

UTy <- crossprod(U,y)
UTx <- crossprod(U,x)

y_tilde <- DIinv%*%UTy
x_tilde <- DIinv%*%UTx


#lambdas <- seq(lambda_pt[1],lambda_pt[2],0.1*sigma_sq)/sigma_sq

CAS_LMM <- CasANOVA(y_tilde,x_tilde,beta.start,sigma_sq,lambdas = lambda_pt*sigma_sq, eps = eps)

beta.t <- round(CAS_LMM$coefficients,-log10(eps)-1)
alleles <- length(unique(beta.t))
lambdafit <- CAS_LMM$lambda/sigma_sq
###
#Write output
output <- c(f,reps,H2_trait,H2_focus,r,prop_causal,di,LOD_add,LOD_dom,LOD_full,alleles,lambdafit,sigma_sq)
write(output,ofile,ncol=length(output),append=T)
}


