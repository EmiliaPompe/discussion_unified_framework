## pcluster_M[z,l] <- marginal of observed cluster
for (j in 1:n){
  if (clsize[lambda[j] ] == 1) ksize <- ksize - 1;
  for (q in 1:n){
    psamp[q] <- 1.0
    old <- F
    old <- any(lambda[-j] == q) ## if old == TRUE, then q labels an observed cluster
    if (old){
      for (h in 1:H){
        # calculate calculate pcluster_M : probability of the cluster q with the xclision of record j  
        pcluster_M[q,h] <- 1.0
        first_element <- TRUE
        for ( i in 1:n){
          if (i != j & lambda[i] == q){
            # print(c(i, first_element))
            if (first_element){
              first_element <- FALSE
              pcluster_M[q,h] <- p[cumdime[h] + V[i,h] ]
            }else{
              cumproduct <- 1.0
              for (k in 1:(i-1)){
                if (k != j & lambda[k] == q) cumproduct <- cumproduct * (( 1- a[q,h]) * (V[k,h] == V[i,h]) + a[q,h] * p[cumdime[h] +  V[k,h]])
              }
              pcluster_M[q,h] <- pcluster_M[q,h] * a[q,h] * p[cumdime[h] + V[i,h]] + (1 - a[q,h]) * p[cumdime[h] + V[i,h]] * cumproduct 
            }
          }
        }
        # calculate numerator 
        prob_M[q, h] <- 1.0
        for (i in 1 : n){
          if (i != j & lambda[i] == q){
            prob_M[q,h] <- prob_M[q,h]*( (V[i,h] == V[j,h])*(1 - a[q,h]) + a[q,h] * p[cumdime[h]+V[i,h]] );
          }
        }
        prob_M[q,h] <- a[q,h]*p[cumdime[h] + V[j,h]] + (1 - a[q,h]) * p[cumdime[h] + V[j,h]] * prob_M[q,h] / pcluster_M[q,h];
        psamp[q] <- psamp[q]*prob_M[q,h];
        if(is.na(psamp[q])) print(c(j,q,h,k,i))
      }
    }else{
      for (h in 1:H){
        psamp[q] = psamp[q] * p[cumdime[h] + V[j,h]]
      }
    }
    if( clsize[q] == 0 | (clsize[q] ==1 & q == lambda[j])) psamp[q] = psamp[q] * (N - ksize ) / (n - ksize)
  }
  ## process psamp 
  psamp <- psamp / sum(psamp)
  clsize[lambda[j]] <- clsize[lambda[j]] - 1 
  lambda[j] <- sample.int(n = n, size = 1, prob = psamp)
  clsize[lambda[j]] <- clsize[lambda[j] ] + 1 
  if (clsize[lambda[j] ] == 1) ksize = ksize + 1
}
