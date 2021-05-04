B2Size=function(shape,S0,x,ta,tf,a,delta,eta,zeta,frac,xi,emax,dist)
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
  zeta=eta*(1-xi)
  for (i in 1:emax){
    ff=function(k){
      eta-pgamma(1,shape=a+i,rate=b+k)
    }
    k=uniroot(ff,c(0, 999))$root
    m=i; k=round(k,2)
    exp=pgamma(delta,shape=a+i,rate=b+k)
    if (exp<=1-zeta) break
  }
  m1=ceiling(m*frac)
  zeta1=eta*(1-xi*frac)
  f1=function(k){
    zeta1-(1-pgamma(delta,shape=a+m1,rate=b+k))
  }
  k1=uniroot(f1,c(0, 999))$root
  n1=ceiling(m1/p)
  n=ceiling(m/p)
  ans=round(c(a, frac, eta, zeta1, zeta, m1,n1,k1,m,n,k),3)
  return(ans)
}
