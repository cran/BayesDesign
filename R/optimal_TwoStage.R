optimal_TwoStage <- function(alphacutoff, powercutoff, S0, x, ta,
                             tf, a = 2, delta, frac = .5, ntrial,
                             complete = "partial", seed = 8232){
  shape <- 1
  dist <- "WB"
  eta <- seq(0.8, 0.95, by = 0.05)
  xi <- seq(0.01, 0.15, by = 0.01)
  oc_mat <- NULL
  for (i in 1:length(eta)){
    for (j in 1:length(xi)){
      bsize <- B2Size(shape=shape,S0=S0,x=x,ta=ta,tf=tf,a = a,delta=delta,
                      eta=eta[i],frac=frac,xi=xi[j],emax=100,dist = dist)
      h0 <- pow_TwoStage(S0=S0,x=x,delta=1,ta=ta,tf=tf,m1=bsize[6],m2=bsize[9],
                         k1=bsize[8],k2=bsize[11],N=500,M=ntrial, seed = seed)
      h1 <- pow_TwoStage(S0=S0,x=x,delta=delta,ta=ta,tf=tf,m1=bsize[6],m2=bsize[9],
                         k1=bsize[8],k2=bsize[11],N=500,M=ntrial, seed = seed)
      ## type I error
      typeI <- h0[1]
      ## power
      power <- h1[1]
      if((typeI <= alphacutoff) & (powercutoff <= power)){
        oc_mat <- rbind(oc_mat, c(eta[i], xi[j], bsize[6], bsize[7], bsize[8],
                                  bsize[9], bsize[10], bsize[11], typeI, power,
                                  h1[2], h1[3], h0[2], h0[3]))
      }
    }
  }
  if (complete == "partial"){
    oc_df <- data.frame(oc_mat)[, -c(1, 2)]
    colnames(oc_df) <- c("m1", "n1", "k1", "m", "n", "k", "typeI", "power", "PET1", "ES1", "PET0", "ES0")
    out <- oc_df[oc_df$n == min(oc_df$n), ]
  }
  else if (complete == "complete"){
    oc_df <- data.frame(oc_mat)
    colnames(oc_df) <- c("eta", "xi", "m1", "n1", "k1", "m", "n", "k", "typeI", "power", "PET1", "ES1", "PET0", "ES0")
    out <- oc_df[oc_df$n == min(oc_df$n), ]
  }
  return(out)
}
