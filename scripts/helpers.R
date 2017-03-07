RIX_geno = function(design){
  rg <- which(design == 1, arr.ind = T)
  return(rg)
}

fill_Z_RIX = function(RIX,N_RIL,P){
  Z_RIX <- rep(0,N_RIL)
  if(P>length(RIX)){
    break
  }
  Z_RIX[RIX[P]] <- 1
  return(Z_RIX)
}

get_poly_effects = function(causal,f,poly_VG){
  
  if(causal==1) {
    effects <- rnorm(f,0,sqrt(poly_VG))
  } else {
    effects <- rep(0,f)
  }
  return(effects)
  
}

#Impute the multi locus genotypes of the RIX's in terms of the founder haplotypes based on the RIL genotypes
fill_G_RIX = function(Z_RIX,GENO_RIL){
  RILS <- which(Z_RIX==1)
  G_RIX <- GENO_RIL[RILS[1],]
  return(G_RIX)
}