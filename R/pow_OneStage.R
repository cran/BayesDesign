pow_OneStage=function(S0,x,delta,ta,tf,m1,k1,N,M,seed)
{
  lambda0=-log(S0)/x
  lambda=delta*lambda0
  s0=function(t){exp(-lambda0*t)}
  tau=ta+tf
  set.seed(seed)
  s=0
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
    x1=x[1:s1]
    cens1=cens[1:s1]
    w1=-log(s0(x1))
    U1=sum(w1)
    if (U1>k1) {
      s=s+1
      pts <- pts + s1
    }
    else{
      pts <- pts + s1
    }
  }
  pow=round(s/M,3)
  pts_rate <- round(pts/M, 3)
  return(c(pow, pts_rate))
}
