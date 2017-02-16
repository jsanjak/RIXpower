load("RIXtest.RData")
#library(lmmlasso)

#X must not be overparamerized, but instead must include the intercerpt term
X <- as.matrix(cbind(rep(1,nrow(MODEL_MAT_1)),MODEL_MAT_1[,c(4,5)]))
Y <- as.matrix(MODEL_MAT_1[,1])
Z <- as.matrix(MODEL_MAT_1[,c(6,7)]) # do not include the cross
GRP <- as.matrix(MODEL_MAT_1[,2])
LMM_lasso <- lmmlasso(x=X,y=Y,z=Z,grp=GRP,lambda=0.1,pdMat="pdDiag")
