GBR_Gvalue = function(W,h,r,a,domi="geom") {#
  N <- nrow(W)
  G <- rep(0,N)
  for (ind in 1:N ){
    G_l <- rep(0,r)
    for ( l in 0:(r-1) ){
      W_l <- W[ind,(l*(h-1)+1):((l+1)*(h-1))]
      a_l <- a[(l*(h)+1):((l+1)*(h))]
      i <- h
      j <- h
      if (sum(W_l!=0) == 2){ # if there are two non zeros then individual is a het
        i <- which(W_l!=0)[1]
        j <- which(W_l!=0)[2]
      }
      if (sum(W_l!=0) == 1){ # if there is one non zero, then  the individual is either a het
        i <- which(W_l!=0)
        if (W_l[i] == 2) { # or homozygous for not the "last" allele
          j <- i
        }
      }
      
      a_i <- a_l[i]
      a_j <- a_l[j]
      #Now asign the individual value based on the haplotypes and the dominance model
      if (domi == "min"){
        G_l[l+1] <- min(c(a_i,a_j))
      } else if ( domi == "harm" ) {
        if (a_i == 0. | a_j == 0.) {
          G_l[l+1] <- 0.
        } else{
          G_l[l+1] <- 1./(((1./a_i+1./a_j)*0.5))
        }
      } else if ( domi == "geom") {
        G_l[l+1] <- sqrt(a_i*a_j)
      } else if ( domi == "add" ) {
        G_l[l+1] <- (a_i+a_j)*0.5
      } else if ( domi == "max" ) {
        G_l[l+1] <- max(c(a_i,a_j))
      } else {
        print("Must provide a valid dominance model")
        break
      }
    }
    G[ind] = sum(G_l)
  }
  return(G)
}

#Assign dominance encodings as 0,1,0 as often done in GWAS.
fill_d_zero = function(W,h,r=1) {
  N = nrow(W)
  d = matrix(0,nrow=N,ncol=(h*(h-1)/2)*r)
  for ( l in 0:(r-1)){
    W_l <-  W[,(l*(h-1)+1):((l+1)*(h-1))]
    for (ind in 1:N) {
      i = h
      j = h
      if (sum(W_l[ind,]!=0) == 2){ # if there are two non zeros then individual is a het
        i <- which(W_l[ind,]!=0)[1]
        j <- which(W_l[ind,]!=0)[2]
      }
      if (sum(W_l[ind,]!=0) == 1){ # if there is one non zero, then  the individual is either a het
        i <- which(W_l[ind,]!=0)
        if (W_l[ind,i] == 2) { # or homozygous for not the "last" allele
          j <- i
        }
      }
      di <- 0
      for (  k in 1:(h-1)){
        for ( f in (k+1):(h)){
          di <- di + 1
          Wd_ij_kf <- 0
          if ( (sum( c(i,j) == c(k,f))==2) | (sum(c(j,i) == c(k,f))==2 ) ){
            Wd_ij_kf <- 1
          }
          d[ind,((h*(h-1)/2)*l)+di] <-  Wd_ij_kf
        }
      }
    }
  }
  return(d)
}

#TDOD
#don't leave this as stupid for loops
R_sq_hap <- function(r,f,allele_freq,GENO_MAT,GENO_PAT,N){
  
  r_sq_hap <- rep(0,r*(r-1)/2)
  for ( K in 1:(r-1)){
    for( L in (K+1):r){
      r_sq_KL <- 0
      for (i in 1:f){
        pA_i <- allele_freq[f*(K-1)+i]
        if( pA_i == 0){next}
        for (j in 1:f){
          pB_j <- allele_freq[f*(L-1)+j]
          if( pB_j == 0){next}
          pAB_ij <- (sum(rowSums(GENO_MAT[,c(f*(K-1)+i,f*(L-1)+j)])==2) +  
                       sum(rowSums(GENO_PAT[,c(f*(K-1)+i,f*(L-1)+j)])==2))/(2*N)
          
          D_KL_ij <- pAB_ij - (pA_i*pB_j)
          r_sq_KL_ij <- (D_KL_ij^2)/(pA_i*(1-pA_i)*pB_j*(1-pB_j))
          
          r_sq_KL <- r_sq_KL + r_sq_KL_ij*pAB_ij
        }
      }
      r_sq_hap[choose(r,2)-choose(r-K+1,2)+L-K] <- r_sq_KL
    }
  }
  return(r_sq_hap)
}