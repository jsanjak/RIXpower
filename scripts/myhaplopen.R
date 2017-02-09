

my.haplo.pen.glm.1 = function (formula = formula(data), family = gaussian, data = sys.parent(), haplofreq = haplofreq, hap.base = hap.base,
                            hapdata = hapdata, na.action = "na.geno.keep", lambdas = seq(0.1,5,0.1), eps = 10^(-6), miss.val = c(0,NA),
                            locus.label = NA, allele.lev = NULL, control = haplo.glm.control(), 
                            method = "glm.fit", model = FALSE, x = FALSE, y = TRUE, contrasts = NULL,  
                            ...)
{
  require(haplo.stats)
  require(MASS)
  lambdas = sort(lambdas)
  
  if (min(lambdas) <=0) {return(cat("ERROR: All values for lambda must be > 0. \n"))}
  if (eps <=0) {return(cat("ERROR: Eps must be > 0. \n"))}
  
  init.fit <- glm(formula =formula,data=data,family = family )
  #init.fit <- glm(formula =as.formula(Y~.), data=mydat, family = "gaussian" )
  init.fit$haplo.freq <- haplofreq
  
  init.est = init.fit$coef #just the coefficients
  p = length(init.est)-1
  cov.est = vcov(init.fit)# vcov for regression coefficients
  inv.cov.est = ginv(cov.est)#[1:(p+1),1:(p+1)] 
  
  X.star = chol(inv.cov.est)
  Y.star = X.star%*%init.est 	
  
  #none of this may be needed
  N=nrow(data)
  X = as.matrix(cbind(rep(1,N),data[,-1])) #the genotype encodings with a column of 1's for the intercept
  y = data[,1] #just the phenotype data
  w = rep(1,N)#init.fit$haplo.post.info$post #posterior probabilities of the haplotypes
  #haplo.names = init.fit$haplo.names #variable names
  
  #init some stuff to be used later
  y1 = y	
  X1 = X
  
  #if we assume we know the haplotypes, then this step is unnecessary
  y = sqrt(w)*y
  X = diag(sqrt(w))%*%X
  ###
  #why do this?
  #this actually gets rid of the column of 1's
  haplo.x.mat = X[,-1]
  
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
  if (length(dm) == p) {dm = matrix(dm,nrow=1)}
  
  #coefs no intercept
  abs.wt.beta = abs(init.est[-1])
  abs.wt.diff = abs(dm %*% init.est[-1])
  
  a.sc.vec = as.vector(1 / c(abs.wt.beta, abs.wt.diff))
  
  full.mat = NULL
  beta.t = as.matrix(init.est)
  
  fit = init.fit
  
  #simply the haplotype assignments
  #without the base removed!!
  gdat1 = hapdata[,1]
  gdat2 = hapdata[,2]
  
  #
  prior.coef = ifelse(gdat1 != gdat2, 2, 1)
  
  #subject IDs
  subj.indx = seq(1,N)#as.numeric(factor(fit$haplo.post.info$indx))
  len.subj.indx = length(subj.indx)
  n.subj = length(unique(subj.indx))
  prior.tot = rep(0, n.subj)
  fitted_mat = NULL
  for(lambda in lambdas)
  {		
    beta.eps = rep(1,length(abs.wt.beta)) 
    while(max(beta.eps)>eps)
    {
      beta.tilde = beta.t[-1]
      abs.tilde = abs(beta.tilde)
      abs.diff = as.vector(abs(dm %*% beta.tilde))
      if (length(abs.diff)<2) {quad.mat.2 = a.sc.vec[-(1:p)]/(abs.diff+eps)} else {quad.mat.2 = diag(a.sc.vec[-(1:p)]/(abs.diff+eps))}
      if (length(abs.tilde)<2) {quad.mat.1 = a.sc.vec[1:p]/(abs.tilde+eps)} else {quad.mat.1 = diag(a.sc.vec[1:p]/(abs.tilde+eps))}
      
      if (length(abs.diff)==1) {quad.mat.2 = 0}
      penalty.mat = lambda*(quad.mat.1 + t(dm) %*% quad.mat.2 %*% dm)
      Sigma.Mat = t(X.star) %*% X.star + rbind(0,cbind(0,penalty.mat))
      beta.t = solve(Sigma.Mat) %*% t(X.star) %*% Y.star
      beta.eps = abs(beta.tilde-beta.t[-1])/(abs(beta.tilde)+eps)
    }     
    df = sum((unique(abs(round(beta.t,-log10(eps)-1)))!=0))
    
    #this is present in the glm object
    fit$fitted.values = X1%*%beta.t
    
    switch(as.character(fit$family[1]), Binomial = , binomial = {
      fit$fitted.values = 1/(1+exp(-fit$fitted.values))
    }, Poisson = , poisson = {
      fit$linear.predictors = fit$fitted.values
    })
    fitted_mat = cbind(fitted_mat,fit$fitted.values)
    #works just fine
    dfit = haplo.stats:::dglm.fit(fit)
    
    #need a haplotype frequency object based on haplotypes withouth base type removed
    prior = prior.coef*fit$haplo.freq[gdat1]*fit$haplo.freq[gdat2]*dfit
    
    tmp.sum = .C("groupsum", x = as.double(prior), indx = as.integer(subj.indx), 
                 n = as.integer(len.subj.indx), grouptot = as.double(prior.tot), 
                 ngroup = as.integer(n.subj), PACKAGE = "haplo.stats")
    pr.pheno = tmp.sum$grouptot
    lnlike.new = sum(log(pr.pheno))
    
    bic = -2*lnlike.new+log(sum(w))*df
    full.mat = rbind(full.mat,c(beta.t,bic,df))
  }
  
  
  BIC.loc = which.min(full.mat[,(p+2)])
  if (lambdas[BIC.loc] == max(lambdas)) {cat("Note: The chosen lambda is the largest of the input grid. \n You should try larger values of lambda \n to continue the search for the minimum BIC. \n")}
  if (lambdas[BIC.loc] == min(lambdas)) {cat("Note: The chosen lambda is the smallest of the input grid. \n You should try smaller values of lambda \n to continue the search for the minimum BIC. \n")}
  BIC.coefs = round(full.mat[BIC.loc, 1:(p+1)],6)
  names(BIC.coefs) = c("(Intercept)",colnames(haplo.x.mat))
  
  fit$fitted_mat = fitted_mat
  fit$haplo.freq = round(fit$haplo.freq,4)
  fit$coefficients = BIC.coefs
  fit$lambda.grid = lambdas
  fit$BIC = full.mat[, p+2]
  fit$df = full.mat[, p+3]
  fit$lambda = lambdas[BIC.loc]
  fit$init.fit = init.fit
  fit$CoefMat = full.mat[, 1:(p+1)]
  fit$BIC.loc = BIC.loc
  return(fit)
}