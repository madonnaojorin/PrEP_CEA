rpert <- function(n,min=-1,mode=0,max=1,shape=4)
  #ISALIAS dpert
  #--------------------------------------------
{
  if (length(n) > 1) 
    n <- length(n)
  if (length(n) == 0 || as.integer(n) == 0) 
    return(numeric(0))
  n <- as.integer(n)
  
  min <- rep(as.vector(min),length.out=n)
  mode <- rep(as.vector(mode),length.out=n)
  max <- rep(as.vector(max),length.out=n)
  shape <- rep(as.vector(shape),length.out=n)
  
  a1 <- 1 + shape * (mode - min)/(max - min)
  a2 <- 1 + shape * (max - mode)/(max - min)
  oldw <- options(warn = -1)
  r <- rbeta(n, shape1 = a1, shape2 = a2) * (max - min) + min
  options(warn = oldw$warn)
  minmodemax <- (abs(min - max) < (.Machine$double.eps^0.5))
  r <- ifelse(minmodemax, min, r)
  if (any(is.na(r))) 
    warning("NaN in rpert")
  return(r)
}