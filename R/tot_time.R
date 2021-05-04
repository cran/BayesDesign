tot_time <- function(obs_time, S0, x){
  lambda0 <- -log(S0)/x
  s0 <- function(t){exp(-lambda0*t)}
  w <- -log(s0(obs_time))
  U <- sum(w)
  return(U)
}
