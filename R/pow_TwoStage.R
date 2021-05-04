pow_TwoStage=function(S0,x,delta,ta,tf,m1,m2,k1,k2,N,M,seed)
{
  lambda0=-log(S0)/x
  lambda=delta*lambda0
  s0=function(t){exp(-lambda0*t)}
  tau=ta+tf
  set.seed(seed)
  s=0
  earlystopping <- 0
  pts <- 0
  for (j in 1:M){
    time=rexp(N, rate=lambda)
    A=runif(N, 0, ta)
    data=data.frame(Entry=A, time=time)
    data0=data[order(data$Entry),]
    time0=data0$time
    A0=data0$Entry
    cens=as.numeric(time0<tau-A0)
    x=pmin(time0,tau-A0)
    d = cumsum(cens)
    s1 = min(which(d == m1))
    s2 = min(which(d == m2))
    x1=x[1:s1]
    cens1=cens[1:s1]
    w1=-log(s0(x1))
    U1=sum(w1)
    x2=x[1:s2]
    w2=-log(s0(x2))
    U2=sum(w2)
    if (U1>k1 & U2>k2) {
      s=s+1
      pts <- pts + s2
    }
    else if (U1 <= k1) {
      earlystopping <- earlystopping + 1
      pts <- pts + s1
    }
    else {
      pts <- pts + s2
    }
  }
  pow=round(s/M,3)
  es_rate <- round(earlystopping/M, 3)
  pts_rate <- round(pts/M, 3)
  return(c(pow, es_rate, pts_rate))
}
