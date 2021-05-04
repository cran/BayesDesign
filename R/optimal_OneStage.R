optimal_OneStage <- function(alphacutoff, powercutoff, S0,
                             x, ta, tf, a = 2, delta, ntrial,
                             complete = "partial", seed = 8232){
  shape <- 1
  dist <- "WB"
  eta <- seq(0.8, 0.95, by = 0.05)
  zeta <- seq(0.8, 0.9, by = 0.01)
  oc_mat <- NULL
  for (i in 1:length(eta)){
    for (j in 1:length(zeta)){
      bsize <- BSize(shape=shape,S0=S0,x=x,ta=ta,tf=tf,a = a,delta=delta,
                      eta=eta[i],zeta=zeta[j],emax=100,dist = dist)
      h0 <- pow_OneStage(S0=S0,x=x,delta=1,ta=ta,tf=tf,m1=bsize[3],k1=bsize[5],
                         N=500,M=ntrial, seed = seed)
      h1 <- pow_OneStage(S0=S0,x=x,delta=delta,ta=ta,tf=tf,m1=bsize[3],k1=bsize[5],
                         N=500,M=ntrial, seed = seed)
      ## type I error
      typeI <- h0[1]
      ## power
      power <- h1[1]
      if((typeI <= alphacutoff) & (powercutoff <= power)){
        oc_mat <- rbind(oc_mat, c(eta[i], zeta[j], bsize[3], bsize[4], bsize[5],
                                  typeI, power, h1[2], h0[2]))
      }
    }
  }
  if (complete == "partial"){
    oc_df <- data.frame(oc_mat)[, -c(1, 2)]
    colnames(oc_df) <- c("m", "n", "k", "typeI", "power", "ES1", "ES0")
    out <- oc_df[oc_df$n == min(oc_df$n), ]
  }
  else if (complete == "complete"){
    oc_df <- data.frame(oc_mat)
    colnames(oc_df) <- c("eta", "zeta", "m", "n", "k", "typeI", "power", "ES1", "ES0")
    out <- oc_df[oc_df$n == min(oc_df$n), ]
  }
  return(out)
}
