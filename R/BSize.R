BSize=function(shape,S0,x,ta,tf,a,delta,eta,zeta,emax,dist)
{
  S1=S0^delta
  tau=ta+tf
  if (dist=="WB"){
    f=function(a,b,u){dweibull(u,a,b)}
    scale1=x/(-log(S1))^(1/shape)
  }
  if (dist=="LN"){
    f=function(a,b,u){dlnorm(u,b,a)}
    scale1=log(x)-shape*qnorm(1-S1)
  }
  if (dist=="LG"){
    f=function(a,b,u){(a/b)*(u/b)^(a-1)/(1+(u/b)^a)^2}
    scale1=x/(1/S1-1)^(1/shape)
  }
  if (dist=="GM"){
    s=function(a,b,u){1-pgamma(u,a,b)}
    f=function(a,b,u){dgamma(u,a,b)}
    root1=function(t){s(shape,t,x)-S1}
    scale1=uniroot(root1,c(0,10))$root
  }
  G=function(u){1-punif(u, tf, tau)}
  g=function(u){f(shape,scale1,u)*G(u)}
  p=integrate(g, 0, tau)$value
  b=2*a/(1+delta)
  for (s in 1:emax){
    ff=function(k){
      pgamma(1,shape=a+s,rate=b+k)-eta
    }
    k=uniroot(ff,c(0,100))$root
    m=s; k=round(k,2)
    exp=pgamma(delta,shape=a+s,rate=b+k)
    if (exp<=1-zeta) break
  }
  n=ceiling(m/p)
  return(c(eta,zeta,m,n,k))
}
BSize(shape=1,S0=0.53,x=3,ta=4,tf=2,a=2,delta=0.6,eta=0.95,
      zeta=0.85,emax=100,dist="WB")
