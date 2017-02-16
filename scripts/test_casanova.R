##################################################################
#
#
#   Function to fit CAS-ANOVA method of Bondell and Reich (2009)
#
#
#
#####################################################################


CasANOVA=function(init.fit,data = sys.parent(), lambdas = seq(0.1,5,0.1), eps = 10^{-6}){
  
  #simple model 4 levels, two effects
  require(MASS)
  N=nrow(data)
  #dat <- data.frame(Y,X)
  lambdas =sort(lambdas)
  eps=1e-6
  
  
  
  #init.fit = lm('Y ~ 0 + .',data=data,x=T,y=T)
  #plot(init.fit$coefficients)
  
  #n.terms = length(attr(terms(init.fit),"dataClasses")[-1])
  #n.levels = as.vector(2*(attr(terms(init.fit),"dataClasses")[-1]=="numeric"))
  #if (min(n.levels)==0) {n.levels[n.levels==0] = as.vector(sapply(init.fit$xlevels,length))}
  #n.factors = length(n.levels)
  
  
  total.p = length(init.fit$coefficients[,1])
  full.x = as.matrix(data[,-1])#init.fit$
  if (ncol(full.x)>nrow(full.x)) {return(cat("ERROR: Number of predictors (including dummy variables created for factors) \n cannot exceed the number of observations. \n"))}
  
  y = data[,1]
  n = length(y)
  
  # Create Difference Matrix
  
  
  p = length(init.fit$coefficients[,1])
  qp = p*(p-1)/2
  
  if (qp==0) temp.dm = matrix(0,nrow=1,ncol=1)
  if (qp>0)
  {
    v1 = rep(0,qp)
    c1 = qp
    c2 = 0
    for(w1 in 1:(p-1))
    {
      c1 = c1-w1
      c2 = c2+(w1-1)
      v1 = cbind(c(rep(0,c1),rep(1,w1),rep(0,c2)),v1)
    }
    v2 = matrix(0,qp,p)
    c1 = 1
    c2 = p-1
    for (w1 in 1:(p-1))
    {
      v2[c1:c2,(w1+1):p] = diag(p-w1)
      c1 = c2+1
      c2 = c2+p-w1-1
    }
    temp.dm = v1 - v2
    
  }
  
  dm = temp.dm
  
  dm = dm[apply(abs(dm),1,sum)>0,]
  if (length(dm) == 0) {dm = matrix(0,nrow=1,ncol=total.p)}
  if (length(dm) == total.p) {dm = matrix(dm,nrow=1)}
  
  beta.mle = init.fit$coefficients[,1]
  
  beta.wt = beta.mle#[-1]
  abs.wt.beta = abs(beta.wt)
  abs.wt.diff = abs(dm %*% beta.mle)#[-1])		
  
  # the design matrix without the intercept times the ginve of the gamma matix
  # this is "Z" in the Bondell and Reich manuscript
  x.sc.mat = full.x %*% t(ginv(cbind(diag(total.p), t(dm))))
  
  #The L2 norm of each column of Z -- this is the proposed non-adaptive weights
  sc.vec = sqrt(apply(x.sc.mat^2,2,sum))
  
  #by dividing the weights by the absolute value of either the mle beta or difference
  #the weights become "adaptive" in that estimates that are small under OLS are more likely
  #to shrink and differences that are small are more like to move towards each other
  a.sc.vec = as.vector(sc.vec / c(abs.wt.beta, abs.wt.diff))
  
  full.mat = NULL
  beta.t = beta.mle
  
  for(lambda in lambdas)
  {		
    beta.eps = rep(1,length(beta.wt)) 
    while(max(beta.eps)>eps)
    {
      beta.tilde = beta.t#[-1]
      abs.tilde = abs(beta.tilde)
      abs.diff = as.vector(abs(dm %*% beta.tilde))
      if (length(abs.diff)<2) {quad.mat.2 = a.sc.vec[-(1:total.p)]/(abs.diff+eps)} else {quad.mat.2 = diag(a.sc.vec[-(1:total.p)]/(abs.diff+eps))}
      if (length(abs.tilde)<2) {quad.mat.1 = a.sc.vec[1:total.p]/(abs.tilde+eps)} else {quad.mat.1 = diag(a.sc.vec[1:total.p]/(abs.tilde+eps))}
      
      if (length(abs.diff)==1) {quad.mat.2 = 0}
      penalty.mat = lambda*(quad.mat.1 + t(dm) %*% quad.mat.2 %*% dm)
      Sigma.Mat = t(full.x) %*%ginv(V_reml_est)%*% full.x + penalty.mat#rbind(0,cbind(0,penalty.mat))
      beta.t = solve(Sigma.Mat) %*% t(full.x) %*% y
      beta.eps = abs(beta.tilde-beta.t)/(abs(beta.tilde)+eps)
    }     
    df = sum((unique(abs(round(beta.t,6)))!=0))
    bic = n*log(sum((y-full.x %*% beta.t)^2))+log(n)*df #this maybe should be mean(RSS), but it doesn't matter for choosing the smallest BIC
    full.mat = rbind(full.mat,c(beta.t,bic,df))
  }
  
  BIC.loc = which.min(full.mat[,(total.p+2)])
  if (lambdas[BIC.loc] == max(lambdas)) {cat("Note: The chosen lambda is the largest of the input grid. \n You should try larger values of lambda \n to continue the search for the minimum BIC. \n")}
  if (lambdas[BIC.loc] == min(lambdas)) {cat("Note: The chosen lambda is the smallest of the input grid. \n You should try smaller values of lambda \n to continue the search for the minimum BIC. \n")}
  BIC.coefs = round(full.mat[BIC.loc, 1:(total.p)],6)
  names(BIC.coefs) = colnames(full.x)
  fit = NULL
  fit$full.mat = full.mat
  fit$coefficients = BIC.coefs
  fit$lambda.grid = lambdas
  fit$BIC = full.mat[, total.p+1]
  fit$df = full.mat[, total.p+2]
  fit$lambda = lambdas[BIC.loc]
  fit$fitted = full.x %*% BIC.coefs
  fit$init.coefs = init.fit$coefficients[,1]
  return(fit)
}
