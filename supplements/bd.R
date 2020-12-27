BDN_sapply <- function(z,n1,n2,p,q) {
  sum <- 0
  
  
  if(is.na(z) || is.nan(z) || is.infinite(z)) {
    return(0)
  }
  
  if (z>=0){  
    domain1 <- 0:n1
    domain2 <- domain1 + z
    sum<-sum(suppressWarnings(sapply(domain2, dbinom, size=n1, prob=p))* suppressWarnings(sapply(domain1, dbinom, size=n2, prob=q)))
  } else {
    domain1 <- 0:n2
    domain2 <- domain1 - z
    sum<-sum(suppressWarnings(sapply(domain1, dbinom, size=n1, prob=p))* suppressWarnings(sapply(domain2, dbinom, size=n2, prob=q)))
  }
  

  return(sum)
}


BDN_vec <- Vectorize(BDN_sapply)


binomprop_estim <- function(x, x0, n0, x_, n_, delta) {
  A <- n_
  B <- (n_ + n0)*delta - (x_ + n_)
  C <- x_ - (2*x0+n_)*delta + n0*delta^2
  D <- delta*(1-delta)*x0
  A*x^3 + B*x^2 + C*x + D
}

binomprop_diffscore <- function(x1, n1, pi, delta) {
  (x1-n1*(pi+delta))/((pi+delta)*(1-(pi+delta)))
}

binomprop_var <- function(n0, n1, pi, delta) {
  1/((pi*(1-pi)/n0) + (pi+delta)*(1-(pi+delta))/n1)
}


binomprop_test <- function(N1,N2,X1,X2,delta=0) {
  
  a <- N1
  b <- N2
  e <- X1
  f <- X2
  
  pi1_mle <- sum(e)/sum(a)
  pi <- pi1_mle
  
    tryCatch(
      pi <- uniroot(binomprop_estim,interval=c(0.0001,0.9999-delta), x0=sum(e), n0=sum(a), x_ = sum(e+f), 
                    n_ = sum(a+b), delta=delta)$root, 
      error=function(err) {
        return(pi1_mle)
      }
    )
  
  g = binomprop_diffscore(sum(f),sum(b),pi, delta)
  gvar = binomprop_var(sum(a),sum(b),pi, delta)
  z = g/sqrt(gvar)
  
  #this might happen in extreme cases
  if(is.nan(z)) {
    print(paste("setting to zero at pi:",pi, " delta:", delta, "g:", g, "gvar:", gvar, "mle", pi1_mle))
    z <- 0
  }
  
  return(z)
}


binomll <- function(x,n,pi) {
  a <- x*log(pi)
  #taking care of all instances where we have no counts
  a <- replace(a,x==0,0)
  b <- (n-x)*log(1-pi)
  #taking care of all instances where we have no counts
  b <- replace(b,(n-x)==0,0)
  
  a+b
}

binomll_logodds <- function(x,n,phi) {
  x*phi + n*log((1/(1+exp(phi))))
}

invlogodds <- function(phi) {
  exp(phi)/(1+exp(phi))
}

binomll_deriv <- function(x,n,pi) {
  -x/(pi) + (n-x)/(1-pi)
}

binomll_deriv_logodds <- function(x,n,phi) {
  x-((n*exp(phi))/(1+exp(phi)))
}


binomdiffll_deriv_delta <- function(x,n,pi,delta) {
  -x/(pi-delta) + (n-x)/(1-pi+delta)
}


binomdiffll_deriv_pi <- function(x1, x2, n1, n2, pi1, delta) {
  x1/(pi1) - (n1-x1)/(1-pi1) + x2/(pi1-delta) - (n2-x2)/(1-pi1+delta) 
}

binomdiffll2_deriv_delta <- function(x1,x2,n1,n2,pi1,pi2,delta1,delta2) {
  t1 <- binomdiffll_deriv_delta(x1,n2,pi1,delta1)
  t2 <- binomdiffll_deriv_delta(x2,n2,pi2,delta2)
  t1+t2
}

binomdiffll2_deriv_pi <- function(x1,x2,x3,x4,n1,n2,n3,n4,pi1,pi2,delta1,delta2) {
  t1 <- binomdiffll_deriv_pi(x1, x2, n1, n2, pi1, delta1)
  t2 <- binomdiffll_deriv_pi(x3, x4, n3, n4, pi3, delta2)
  t1+t2
}

binomdiffll_pi_find <- function(x1,x2,n1,n2,delta,eps=10^-6, maxiter=100) {
  
  k = 0
  grad = 1
  
  l <- delta+eps
  r <- 1-eps
  
  if(l > 1-eps) return(1-eps);
  
  #changed here for groups of samples
  mid <- sum(x1)/sum(n1)
  if(mid < l || mid >r) {
    mid = (l+r)/2
  }
  
  x <- mid
  while (k == 0 || (k < maxiter && abs(grad) > eps)) {
    
    # gradient of objective function
    x <- mid
    #changed here for groups of samples
    grad <- sum(binomdiffll_deriv_pi(x1, x2, n1, n2, x, delta))
    if(is.infinite(grad) || is.nan(grad)) {
      print("exiting irregularily")
      return(x)
    }
    # search direction
    dir <- sign(grad)
    #go right
    if(dir == 1) {
      l <- mid
      mid <- max((mid+r)/2,delta+eps)
      #left boundary reached
      if(l == mid && mid == delta+eps) {
        return(mid)
      }
    }
    #go left
    if(dir == -1) {
      r <- mid
      mid <- max((mid+l)/2,delta+eps)
      #right boundary reached
      if(mid == r && mid == delta+eps) {
        return(mid)
      }
    }
    # iterations counter
    k <- k + 1
  }
  
  return(x)
}

optim_pi1 <- function(x1,x2,n1,n2,delta1,eps=10^-6, maxiter=100) {
  if(delta1 >= 0) {
    pi1 <- binomdiffll_pi_find(x1,x2,n1,n2,delta=delta1,eps=eps,maxiter=maxiter) 
    pi2 <- max(eps,pi1-delta1)
  } else {
    pi2 <- binomdiffll_pi_find(x2,x1,n2,n1,delta=-delta1,eps=eps, maxiter=maxiter) 
    pi1 <- max(eps,pi2+delta1)
  }
  
  c(sum(binomll(x1,n1,pi1)) + sum(binomll(x2,n2,pi2)), pi1, pi2)
}

optim_pi2 <- function(x1,x2,x3,x4,n1,n2,n3,n4,delta1,delta2,eps=10^-6, maxiter=100) {
  
  sum1 <- optim_pi1(x1,x2,n1,n2,delta1,eps=eps,maxiter=maxiter)
  sum2 <- optim_pi1(x3,x4,n3,n4,delta2,eps=eps,maxiter=maxiter)
  c(sum1[1]+sum2[1], sum1[2], sum2[2], sum1[3], sum2[3])  
}

ternary_H0 <- function(x1,x2,x3,x4,n1,n2,n3,n4,eps=10^-6) {
  
  l <- -1+eps
  r <- 1-eps
  
  while (abs(r - l) >= eps) {
    lt = l + (r - l) / 3
    rt = r - (r - l) / 3
    
    v_l <- optim_pi2(x1,x2,x3,x4,n1,n2,n3,n4,lt,lt)
    v_r <- optim_pi2(x1,x2,x3,x4,n1,n2,n3,n4,rt,rt)
    
    if(v_l[1] < v_r[1]) {
      l <- lt 
    } else {
      r <- rt      
    }
  }
  c(v_l[1],v_l[2],v_l[4],l)
} 

ternary_H1 <- function(x1,x2,n1,n2,eps=10^-6) {
  
  l <- -1+eps
  r <- 1-eps
  
  while (abs(r - l) >= eps) {
    lt = l + (r - l) / 3
    rt = r - (r - l) / 3
    
    v_l <- optim_pi1(x1,x2,n1,n2,lt)
    v_r <- optim_pi1(x1,x2,n1,n2,rt)
    
    if(v_l[1] < v_r[1]) {
      l <- lt 
    } else {
      r <- rt      
    }
  }
  c(v_l[1],v_l[2],v_l[3],l)
} 


bdnci_search <- function(a,b,N1,N2,X1,X2,pi1,off,eps=10^-6,searchdir=1,maxiter=100) {
  
  k = 0
  grad = 100
  
  #left, right, mid
  l = a
  r = b
  mid <- (a+b)/2
  
  while (k == 0 || (k < maxiter && abs(grad) > eps)) {
    # gradient of objective function
    grad<- optim_pi1(X1,X2,N1,N2,mid)[1] + off #binomll(X2,N2,pi1+mid) + binomll(X1,N1,pi1) + off
    # search direction
    dir <- sign(grad)
    #go right
    if(dir == searchdir) {
      l <- mid
      mid <- (mid+r)/2
    }
    #go left
    if(dir == -1*searchdir) {
      r <- mid
      mid <- (l+mid)/2
    }
    k <- k + 1
  }
  
  return(mid)
}

bdnci <- function(a,b,e,f,pi1,delta,epsd=10^-6,eps=10^-6,maxiter=100,chialphahalf=1.920729) {
  
  a1 <- -1+eps 
  b1 <- 1-eps 
  
  cil <- delta
  cir <- delta
  
  if(a1 < b1) {
    off <- binomll(f,b,pi1-delta) + binomll(e,a,pi1)
    if(!is.nan(off) && !is.infinite(off)) {
      if(a1 < delta) {
        #shift upwards and find left root
        cil <- bdnci_search(a1,delta,a,b,e,f,pi1,-off+chialphahalf,searchdir=-1,eps=eps,maxiter=maxiter)
      } else {
        cil <- delta
      }
      #shift upwards and find right root
      if(delta < b1) {
        cir <- bdnci_search(delta,b1,a,b,e,f,pi1,-off+chialphahalf,eps=eps,maxiter=maxiter)
      } else {
        cir <- delta
      }
      
      if(abs(cil) == epsd) {
        cil = 0
      }
      
      if(abs(cir) == epsd) {
        cir = 0
      }
    }
  }
  return(tibble(cil,cir))
}



bdnci2_search <- function(a,b,N1,N2,N3,N4,X1,X2,X3,X4,delta1,delta2,off,eps=10^-6,searchdir=1,maxiter=100) {
  
  k = 0
  grad = 100
  
  #left, right, mid
  l = a
  r = b
  mid <- (a+b)/2
  
  while (k == 0 || (k < maxiter && abs(grad) > eps)) {
    # gradient of objective function
    grad <- optim_pi2(X1,X2,X3,X4,N1,N2,N3,N4,
                      delta1-mid,delta2-mid)[1]+off
    # search direction
    dir <- sign(grad)
    #go right
    if(dir == searchdir) {
      l <- mid
      mid <- (mid+r)/2
    }
    #go left
    if(dir == -1*searchdir) {
      r <- mid
      mid <- (l+mid)/2
    }
    k <- k + 1
  }
  #print(k)
  return(mid)
}


bdnci2 <- function(a,b,c,d,e,f,g,h,pi1,pi2,delta1,delta2, eps=10^-6,maxiter=100, chialphahalf=1.920729) {
  
  a1 <- -1+(max(abs(delta2),abs(delta1)))
  b1 <- 1-(max(abs(delta1),abs(delta2)))
  dd <- 0
  
  cil <- dd
  cir <- dd
  
  if(a1 < b1) {
    off <- binomll(f,b,pi1-delta1) + binomll(e,a,pi1) + binomll(h,d,pi2-delta2) + binomll(g,c,pi2)
    if(!is.nan(off) && !is.infinite(off)) {
      if(a1 < dd) {
        #shift upwards and find left root
        cil <- bdnci2_search(a1,dd,a,b,c,d,e,f,g,h, delta1, delta2, -off+chialphahalf,searchdir=-1,eps=eps, maxiter=maxiter)
      } else {
        cil <- dd
      }
      #shift upwards and find right root
      if(dd < b1) {
        cir <- bdnci2_search(dd,b1,a,b,c,d,e,f,g,h, delta1, delta2, -off+chialphahalf, eps=eps,maxiter=maxiter)
      } else {
        cir <- dd
      }
      
      if(abs(cil) == eps) {
        cil = 0
      }
      
      if(abs(cir) == eps) {
        cir = 0
      }
    }
  }
  cil <- max(-2,min(2,cil))
  cir <- max(-2,min(2,cir))
  return(tibble(cil,cir))
}

