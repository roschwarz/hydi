library(tidyverse)
library(multidplyr)
source("bd.R")


simulate_grid_script <- function(p1, q1, p2, q2, replicates=5, 
                     cov.mean.modC = 40, cov.sd.modC = 10, 
                     cov.mean.mC = cov.mean.modC, cov.sd.mC = cov.sd.modC) {
  
  a <- pmax(1,round(rnbinom(replicates,mu=cov.mean.modC,size=cov.mean.modC/2)))
  b <- pmax(1,round(rnbinom(replicates,mu=cov.mean.mC,size=cov.mean.mC/2)))
  c <- pmax(1,round(rnbinom(replicates,mu=cov.mean.modC,size=cov.mean.modC/2)))
  d <- pmax(1,round(rnbinom(replicates,mu=cov.mean.mC,size=cov.mean.mC/2)))
  
  e <- rbinom(replicates,size=a,prob=p1)
  f <- rbinom(replicates,size=b,prob=q1)

  g <- rbinom(replicates,size=c,prob=p2)
  h <- rbinom(replicates,size=d,prob=q2)

  X <- data.frame(cov.modC.1 = a, cov.mC.1 = b, cov.modC.2 = c, cov.mC.2 = d,
               p1=p1,q1=q1,p2=p2,q2=q2,modC.1=e,mC.1=f,modC.2=g,mC.2=h)
   
  X
}


simulate_random_script <- function(replicates=5, 
                                 cov.mean.modC = 40, cov.sd.modC = 10, 
                                 cov.mean.mC = cov.mean.modC, cov.sd.mC = cov.sd.modC, min=0, max=1) {
  
  
  a <- pmax(1,round(rnbinom(replicates,mu=cov.mean.modC,size=cov.mean.modC/2)))
  b <- pmax(1,round(rnbinom(replicates,mu=cov.mean.mC,size=cov.mean.mC/2)))
  c <- pmax(1,round(rnbinom(replicates,mu=cov.mean.modC,size=cov.mean.modC/2)))
  d <- pmax(1,round(rnbinom(replicates,mu=cov.mean.mC,size=cov.mean.mC/2)))
  
  p1 <- signif(runif(1,min=min,max=max),digits=3)
  p2 <- signif(runif(1,min=min,max=max),digits=3)
  
  delta1_ <- signif(runif(1,min=0,max=p1),digits=3)
  delta2_ <- signif(runif(1,min=0,max=p2),digits=3)
  
  q1 <- p1 - delta1_
  q2 <- p2 - delta2_
  
  
  e <- rbinom(replicates,size=a,prob=p1)
  f <- rbinom(replicates,size=b,prob=q1)
  g <- rbinom(replicates,size=c,prob=p2)
  h <- rbinom(replicates,size=d,prob=q2)
  

  X <- data.frame(cov.modC.1 = a, cov.mC.1 = b, cov.modC.2 = c, cov.mC.2 = d,
                  p1=p1,q1=q1,p2=p2,q2=q2,delta1=delta1_,delta2=delta2_,modC.1=e,mC.1=f,modC.2=g,mC.2=h)
  
  X
}


simulate_random_commondelta_script <- function(min=0, max=1, replicates=5, 
                                cov.mean.modC = 40, cov.sd.modC = 10, 
                                cov.mean.mC = cov.mean.modC, cov.sd.mC = cov.sd.modC, reverse=F) {
  
  a <- pmax(1,round(rnbinom(replicates,mu=cov.mean.modC,size=cov.mean.modC/2)))
  b <- pmax(1,round(rnbinom(replicates,mu=cov.mean.mC,size=cov.mean.mC/2)))
  c <- pmax(1,round(rnbinom(replicates,mu=cov.mean.modC,size=cov.mean.modC/2)))
  d <- pmax(1,round(rnbinom(replicates,mu=cov.mean.mC,size=cov.mean.mC/2)))
  
  p <- runif(2,min=min,max=max)
  delta <- runif(1,min=0,max=max(0,min(p)-min))
  
  if(reverse) {
    e <- rbinom(replicates,size=a,prob=p[1]-delta)
    f <- rbinom(replicates,size=b,prob=p[1])
    
    
    g <- rbinom(replicates,size=c,prob=p[2]-delta)
    h <- rbinom(replicates,size=d,prob=p[2])
    
  } else {
    e <- rbinom(replicates,size=a,prob=p[1])
    f <- rbinom(replicates,size=b,prob=p[1]-delta)
    
    g <- rbinom(replicates,size=c,prob=p[2])
    h <- rbinom(replicates,size=d,prob=p[2]-delta)
    
  }

  
  
  X <- data.frame(cov.modC.1 = a, cov.mC.1 = b, cov.modC.2 = c, cov.mC.2 = d,
                  p1=p[1],p2=p[2],delta=delta,modC.1=e,mC.1=f,modC.2=g,mC.2=h)
  
  X
}

simulate_random_deltanull_script <- function(min=0.0, max=1.0, replicates=5, 
                                   cov.mean.modC = 40, cov.sd.modC = 10, 
                                   cov.mean.mC = cov.mean.modC, cov.sd.mC = cov.sd.modC) {
  
  a <- pmax(1,round(rnbinom(replicates,mu=cov.mean.modC,size=cov.mean.modC/2)))
  b <- pmax(1,round(rnbinom(replicates,mu=cov.mean.mC,size=cov.mean.mC/2)))
  c <- pmax(1,round(rnbinom(replicates,mu=cov.mean.modC,size=cov.mean.modC/2)))
  d <- pmax(1,round(rnbinom(replicates,mu=cov.mean.mC,size=cov.mean.mC/2)))
  
  p <- runif(2,min=min,max=max)
  delta <- 0
  
  
  e <- rbinom(replicates,size=a,prob=p[1])
  f <- rbinom(replicates,size=b,prob=p[1]-delta)
  
  g <- rbinom(replicates,size=c,prob=p[2])
  h <- rbinom(replicates,size=d,prob=p[2]-delta)
  
  X <- data.frame(cov.modC.1 = a, cov.mC.1 = b, cov.modC.2 = c, cov.mC.2 = d,
                  p1=p[1],p2=p[2],delta=delta,modC.1=e,mC.1=f,modC.2=g,mC.2=h)
  
  X
}


hydi_test <- function(X,flip=F,verbose=F) {

  ###
  #  data
  ##
  
  a <- X$cov.modC.1
  b <- X$cov.mC.1
  c <- X$cov.modC.2
  d <- X$cov.mC.2
  
  e <- X$modC.1
  f  <- X$mC.1
  g <- X$modC.2
  h  <- X$mC.2
  
  a_pooled <- sum(a)
  b_pooled <- sum(b)
  c_pooled <- sum(c)
  d_pooled <- sum(d)
  
  e_pooled <- sum(e)
  f_pooled <- sum(f)
  g_pooled <- sum(g)
  h_pooled <- sum(h)
  
  ###
  #  mle
  ##
  
  #group 1
  p1_mle <- sum(e)/sum(a)
  q1_mle <- sum(f)/sum(b) 
  d1_mle <- p1_mle-q1_mle
  
  #group 2
  p2_mle <- sum(g)/sum(c)
  q2_mle <- sum(h)/sum(d)
  d2_mle <- p2_mle-q2_mle
  
  d0_mle <- (p1_mle-q1_mle) - (p2_mle-q2_mle)
  
  ###
  #  testing absence of hydroxymethylation
  ##
  
  # testing absence of hydroxymethylation or overshoots in G1
  z_C1 <- binomprop_test(a,b,e,f)
  p_C1 <- 2*pnorm(abs(z_C1),lower.tail=F)

  # testing absence of hydroxymethylation or overshoots in G2
  z_C2 <- binomprop_test(c,d,g,h)
  p_C2 <- 2*pnorm(abs(z_C2),lower.tail=F)

  ###
  #  testing equality of hydroxymethylation (single)
  ##
  
  H00 <- ternary_H0(e,f,g,h,a,b,c,d)
  H11 <- ternary_H1(e,f,a,b)
  H12 <- ternary_H1(g,h,c,d)
  lr <- -2*(H00[1] - (H11[1]+H12[1]))
  p <- pchisq(lr, df=1, lower.tail=F)
  
  ###
  #  testing equality of hydroxymethylation (pooled)
  ##
  
  H00_pooled <- ternary_H0(e_pooled,f_pooled,g_pooled,h_pooled,a_pooled,b_pooled,c_pooled,d_pooled)
  H11_pooled <- ternary_H1(e_pooled,f_pooled,a_pooled,b_pooled)
  H12_pooled <- ternary_H1(g_pooled,h_pooled,c_pooled,d_pooled)
  lr_pooled <- -2*(H00_pooled[1] - (H11_pooled[1]+H12_pooled[1]))
  p_pooled <- pchisq(lr_pooled, df=1, lower.tail=F)
  
  ###
  #  confidence intervals for each group
  ##
  
  #first binomial proportion
  d12_pll<- H11_pooled[2] - H11_pooled[3]
  pi121 <- H11_pooled[2]
  ci1 <- bdnci(a_pooled,b_pooled,e_pooled,f_pooled,pi121,d12_pll)
  
  #second binomial proportion
  d22_pll<- H12_pooled[2] - H12_pooled[3]
  pi221 <- H12_pooled[2]
  ci2 <- bdnci(c_pooled,d_pooled,g_pooled,h_pooled,pi221,d22_pll)
  
  #difference of differences
  d0_pll <- d12_pll - d22_pll
  ci3 <- bdnci2(a_pooled,b_pooled,c_pooled,d_pooled,
                       e_pooled,f_pooled,g_pooled,h_pooled, pi121,pi221,d12_pll,d22_pll)
  
  ci3 <- ci3 %>% mutate(lo=(d12_pll+cil)-(d22_pll+cir), hi=(d12_pll+cir)-(d22_pll+cil))
  
  
  return(tibble(p1_mle, q1_mle, p2_mle, q2_mle,
    d1_mle,  d1_cil=ci1$cil, d12_pll, d1_cir=ci1$cir, z_C1, p_C1,
    d2_mle,  d2_cil=ci2$cil, d22_pll, d2_cir=ci2$cir, z_C2, p_C2,
    d0_mle,  d0_lo=ci3$lo, d0_pll, d0_hi=ci3$hi, lr, lr_pooled, p, p_pooled))
}

simulate_random_worker <- function(coverage, replicates) {
  
  S <- simulate_random_script(replicates=replicates, cov.mean.modC = coverage, cov.sd.modC=coverage*0.1, 
                              cov.sd.mC=coverage*0.1)
  p1=S$p1[1]
  p2=S$p2[1]
  q1=S$q1[1]
  q2=S$q2[1]
  delta1 = S$delta1[1]
  delta2 = S$delta2[1]
  
  tab = hydi_test(S)
  return(tibble(p1,q1,p2,q2,delta1,delta2,tab))
}


simulate_random_commondelta_worker <- function(coverage, replicates, reverse=F) {
  
  S <- simulate_random_commondelta_script(replicates=replicates, cov.mean.modC = coverage, cov.sd.modC=coverage*0.1, 
                          cov.sd.mC=coverage*0.1, reverse = reverse)
  p1=S$p1[1]
  p2=S$p2[1]
  delta=S$delta[1]

  tab = hydi_test(S)
  return(tibble(p1,p2,q1=p1-delta,q2=p2-delta,delta,tab))
}


simulate_random_deltanull_worker <- function(coverage, replicates) {
  
  S <- simulate_random_deltanull_script(replicates=replicates, cov.mean.modC = coverage, cov.sd.modC=coverage*0.1, 
                              cov.sd.mC=coverage*0.1)
  p1=S$p1[1]
  p2=S$p2[1]
  delta=S$delta[1]

  tab = hydi_test(S)
  return(tibble(p1,p2,q1=p1-delta,q2=p2-delta,delta,tab))
}



simulate_grid_worker <- function(coverage, replicates, p1, q1, p2, q2,devfrac=0.1,flip=T) {
  
  S <- simulate_grid_script(p1, q1, p2, q2, cov.mean.modC = coverage, cov.sd.modC=coverage*devfrac, 
                           cov.sd.mC=coverage*devfrac, replicates=replicates)
  tab = hydi_test(S,flip)
  return(tab)
}

simulate_grid_parallel <- function(N=5000,nclust=50) {
  
  cl <- new_cluster(nclust)
  cluster_library(cl, "tidyverse")
  cluster_library(cl, "dplyr")
  cluster_library(cl, "purrr")
  cluster_copy(cl, "simulate_grid_worker")
  cluster_copy(cl, "simulate_grid_script")
  
  cluster_copy(cl, "hydi_test")
  cluster_copy(cl, "binomll")
  cluster_copy(cl, "binomdiffll_deriv_pi")
  cluster_copy(cl, "binomdiffll_pi_find")
  cluster_copy(cl, "optim_pi1")
  cluster_copy(cl, "optim_pi2")
  cluster_copy(cl, "ternary_H1")
  cluster_copy(cl, "ternary_H0")
  
  
  cluster_copy(cl,"binomprop_estim")
  cluster_copy(cl,"binomprop_diffscore")
  cluster_copy(cl,"binomprop_var")
  cluster_copy(cl,"binomprop_test")
  cluster_copy(cl, "bdnci_search")
  cluster_copy(cl, "bdnci")
  cluster_copy(cl, "bdnci2_search")
  cluster_copy(cl, "bdnci2")

  
  cov <- c(30,40,50,100,200)
  reps <- c(5,10,20)
  
  p1s <- c(0.10,0.50,0.90, 0.10, 0.50,0.50,0.50,  0.90,0.90,0.90, 0.10,0.10,0.50,0.50, 0.90,0.90,0.50,0.90,1.00,0.80)
  q1s <- c(0.10,0.50,0.90, 0.00, 0.45,0.10,0.00,  0.50,0.10,0.00, 0.10,0.10,0.50,0.50, 0.50,0.50,0.40,0.80,0.75,0.40)
  
  p2s <- c(0.10,0.50,0.90, 0.10, 0.50,0.50,0.50,  0.90,0.90,0.90, 0.50,0.90,0.90,0.00, 0.50,0.40,0.20,0.40,0.25,0.60)
  q2s <- c(0.10,0.50,0.90, 0.00, 0.45,0.10,0.00,  0.50,0.10,0.00, 0.50,0.90,0.90,0.00, 0.10,0.00,0.10,0.30,0.00,0.20)

  
  Q <- crossing(coverage=cov, replicates=reps, combination=1:length(p1s), sampleno=1:N) %>% 
    mutate(p1=p1s[combination],q1=q1s[combination],p2=p2s[combination],q2=q2s[combination])#,
          #  d0_mle=NA, d0_lin=NA, d0_par=NA, ll_mle=NA, ll_lin=NA, ll_par=NA, p_mle=NA, p_lin=NA, p_par=NA)
  
  Q_part <- Q %>% partition(cl)
  tmp <- Q_part %>% mutate(lrt=pmap(list(coverage, replicates, p1, q1, p2, q2), simulate_grid_worker)) 
  Q <- tmp %>% collect() %>% unnest(lrt)
  
  save(Q,file="nullmodel_grid_19.Rda")
}


simulate_random_commondelta_parallel <- function(N=50000,nclust=50) {
  
  cl <- new_cluster(nclust)
  cluster_library(cl, "tidyverse")
  cluster_library(cl, "dplyr")
  cluster_library(cl, "purrr")
  cluster_copy(cl, "simulate_random_commondelta_worker")
  cluster_copy(cl, "simulate_random_commondelta_script")
  
  cluster_copy(cl, "hydi_test")
  cluster_copy(cl, "binomll")
  cluster_copy(cl, "binomdiffll_deriv_pi")
  cluster_copy(cl, "binomdiffll_pi_find")
  cluster_copy(cl, "optim_pi1")
  cluster_copy(cl, "optim_pi2")
  cluster_copy(cl, "ternary_H1")
  cluster_copy(cl, "ternary_H0")
  
  
  cluster_copy(cl,"binomprop_estim")
  cluster_copy(cl,"binomprop_diffscore")
  cluster_copy(cl,"binomprop_var")
  cluster_copy(cl,"binomprop_test")
  cluster_copy(cl, "bdnci_search")
  cluster_copy(cl, "bdnci")
  cluster_copy(cl, "bdnci2_search")
  cluster_copy(cl, "bdnci2")
  
  
  cov <- c(30,40,50,100,200)
  reps <- c(5,10,20)

  Q <- crossing(coverage=cov, replicates=reps, combination=99, sampleno=1:N, reverse=F) 
  
  Q_part <- Q %>% partition(cl)
  tmp <- Q_part %>% mutate(lrt=pmap(list(coverage, replicates, reverse), simulate_random_commondelta_worker)) 
  Q <- tmp %>% collect() %>% unnest(lrt)
  
  save(Q,file="nullmodel_random_commondelta_19.Rda")
}

simulate_random_deltanull_parallel <- function(N=50000,nclust=50) {
  
  cl <- new_cluster(nclust)
  cluster_library(cl, "tidyverse")
  cluster_library(cl, "dplyr")
  cluster_library(cl, "purrr")
  cluster_copy(cl, "simulate_grid_worker")
  cluster_copy(cl, "simulate_grid_script")
  cluster_copy(cl, "simulate_random_deltanull_worker")
  cluster_copy(cl, "simulate_random_deltanull_script")
  
  cluster_copy(cl, "hydi_test")
  cluster_copy(cl, "binomll")
  cluster_copy(cl, "binomdiffll_deriv_pi")
  cluster_copy(cl, "binomdiffll_pi_find")
  cluster_copy(cl, "optim_pi1")
  cluster_copy(cl, "optim_pi2")
  cluster_copy(cl, "ternary_H1")
  cluster_copy(cl, "ternary_H0")
  
  
  cluster_copy(cl,"binomprop_estim")
  cluster_copy(cl,"binomprop_diffscore")
  cluster_copy(cl,"binomprop_var")
  cluster_copy(cl,"binomprop_test")
  cluster_copy(cl, "bdnci_search")
  cluster_copy(cl, "bdnci")
  cluster_copy(cl, "bdnci2_search")
  cluster_copy(cl, "bdnci2")
  
  cov <- c(30,40,50,100,200)
  reps <- c(5,10,20)
  
  Q <- crossing(coverage=cov, replicates=reps, combination=99, sampleno=1:N) 
  
  Q_part <- Q %>% partition(cl)
  tmp <- Q_part %>% mutate(lrt=pmap(list(coverage, replicates), simulate_random_deltanull_worker)) 
  Q <- tmp %>% collect() %>% unnest(lrt)
  
  save(Q,file="nullmodel_random_deltanull_18.Rda")
}


simulate_random_parallel <- function(N=50000,nclust=50) {
  
  cl <- new_cluster(nclust)
  cluster_library(cl, "tidyverse")
  cluster_library(cl, "dplyr")
  cluster_library(cl, "purrr")
  cluster_copy(cl, "simulate_grid_worker")
  cluster_copy(cl, "simulate_grid_script")
  cluster_copy(cl, "simulate_random_worker")
  cluster_copy(cl, "simulate_random_script")
  
  cluster_copy(cl, "hydi_test")
  cluster_copy(cl, "binomll")
  cluster_copy(cl, "binomdiffll_deriv_pi")
  cluster_copy(cl, "binomdiffll_pi_find")
  cluster_copy(cl, "optim_pi1")
  cluster_copy(cl, "optim_pi2")
  cluster_copy(cl, "ternary_H1")
  cluster_copy(cl, "ternary_H0")
  
  
  cluster_copy(cl,"binomprop_estim")
  cluster_copy(cl,"binomprop_diffscore")
  cluster_copy(cl,"binomprop_var")
  cluster_copy(cl,"binomprop_test")
  cluster_copy(cl, "bdnci_search")
  cluster_copy(cl, "bdnci")
  cluster_copy(cl, "bdnci2_search")
  cluster_copy(cl, "bdnci2")
  
  cov <- c(30,40,50,100,200)
  reps <- c(5,10,20)
  
  Q <- crossing(coverage=cov, replicates=reps, combination=99, sampleno=1:N) 
  
  Q_part <- Q %>% partition(cl)
  tmp <- Q_part %>% mutate(lrt=pmap(list(coverage, replicates), simulate_random_worker)) 
  Q <- tmp %>% collect() %>% unnest(lrt)
  
  save(Q,file="nullmodel_random_19.Rda")
}


simulate_grid_ll_parallel <- function(nclust=50) {
  
  cl <- new_cluster(nclust)
  cluster_library(cl, "tidyverse")
  cluster_library(cl, "dplyr")
  cluster_library(cl, "purrr")
  cluster_copy(cl, "simulate_grid_worker")
  cluster_copy(cl, "simulate_grid_script")
  
  cluster_copy(cl, "hydi_test")
  cluster_copy(cl, "binomll")
  cluster_copy(cl, "binomdiffll_deriv_pi")
  cluster_copy(cl, "binomdiffll_pi_find")
  cluster_copy(cl, "optim_pi1")
  cluster_copy(cl, "optim_pi2")
  cluster_copy(cl, "ternary_H1")
  cluster_copy(cl, "ternary_H0")
  
  
  cluster_copy(cl,"binomprop_estim")
  cluster_copy(cl,"binomprop_diffscore")
  cluster_copy(cl,"binomprop_var")
  cluster_copy(cl,"binomprop_test")
  cluster_copy(cl, "bdnci_search")
  cluster_copy(cl, "bdnci")
  cluster_copy(cl, "bdnci2_search")
  cluster_copy(cl, "bdnci2")
  
  cov <- c(200)
  reps <- c(20)
  devfrac <- c(0.0,0.1)
  flip <- c(T,F)
  
  p1s <- c(0.10,0.25,0.50,0.90,0.25,0.5,0.25)
  q1s <- c(0.10,0.25,0.50,0.90,0.5,0.25,0)
  
  Q <- crossing(coverage=cov, devfrac=devfrac, replicates=reps, p1=p1s, q1=q1s, p2=seq(0.05,0.95,by=0.025),q2=seq(0.05,0.95,by=0.025))
  
  Q_part <- Q %>% partition(cl)
  tmp <- Q_part %>% mutate(lrt=pmap(list(coverage, replicates, p1, q1, p2, q2, devfrac), simulate_grid_worker)) 
  Q <- tmp %>% collect() %>% unnest(lrt)
  
  save(Q,file="nullmodel_grid_ll_19.Rda")
}



simulate_hydi_data <- function(N=10000, replicates=5, coverage=30, analyse=T, dir=".") {
  Q1 <- tibble()
  Q2 <- tibble()
  
  pval1 <- c()
  pval2 <- c()
  pval3 <- c()
  
  P <- tibble()
  for(i in 1:N) {
    S <- simulate_random_commondelta_script(min=0.01, replicates=replicates, cov.mean.modC = coverage, cov.sd.modC=coverage*0.1, 
                                cov.sd.mC=coverage*0.1)
    p1=S$p1[1]
    p2=S$p2[1]
    delta=S$delta[1]


    selectstr1 <- rep(c("cov.modC.1_","modC.1_","cov.mC.1_","mC.1_"), times=nrow(S))
    selectstr1 <- paste0(selectstr1, rep(1:nrow(S), each=4))
    selectstr2 <- rep(c("cov.modC.2_","modC.2_","cov.mC.2_","mC.2_"), times=nrow(S))
    selectstr2 <- paste0(selectstr2, rep(1:nrow(S), each=4))

    V1 <- S %>% select(cov.modC.1, modC.1, cov.mC.1, mC.1) %>% mutate(id= 1:n())
    V2 <- S %>% select(cov.modC.2, modC.2, cov.mC.2, mC.2) %>% mutate(id= 1:n())
    
    V1 <- V1 %>% pivot_wider(id_cols=id,names_from=c(id),values_from=c(cov.modC.1,modC.1,cov.mC.1,mC.1), names_glue = '{.value}_{id}')
    V1 <- V1 %>% select(selectstr1)

    V2 <- V2 %>% pivot_wider(id_cols=id,names_from=c(id),values_from=c(cov.modC.2,modC.2,cov.mC.2,mC.2), names_glue = '{.value}_{id}')
    V2 <- V2 %>% select(selectstr2)
    
    V1 <- bind_cols(chrom="chr1",pos=i,strand="+",V1)
    V2 <- bind_cols(chrom="chr1",pos=i,strand="+",V2)
    
    if(analyse) {
      tab = hydi_test(S)
      
      P <- bind_rows(P,tibble(chrom="chr1", pos=i, strand="+",
                              cil1=sprintf("%f",tab$d1_cil), mid1=sprintf("%f",tab$d12_pll), cir1=sprintf("%f",tab$d1_cir), p1=sprintf("%f",tab$p_C1), fdr1=sprintf("%f",tab$p_C1),
                              cil2=sprintf("%f",tab$d2_cil), mid2=sprintf("%f",tab$d22_pll), cir2=sprintf("%f",tab$d2_cir), p2=sprintf("%f",tab$p_C2), fdr2=sprintf("%f",tab$p_C2),
                              cil3=sprintf("%f",tab$d0_lo),  mid3=sprintf("%f",tab$d0_pll),  cir3=sprintf("%f",tab$d0_hi),  p3=sprintf("%f",tab$p), fdr3=sprintf("%f",tab$p)))
      
      
      pval1 <- c(pval1,tab$p_C1)
      pval2 <- c(pval2,tab$p_C2)
      pval3 <- c(pval3,tab$p)
    }
    
    Q1 <- bind_rows(Q1,V1)
    Q2 <- bind_rows(Q2,V2)
  }
  
  if(analyse) {
    f1 <- p.adjust(pval1, method="fdr")
    f2 <- p.adjust(pval2, method="fdr")  
    f3 <- p.adjust(pval3, method="fdr")
    
    
    P <- P %>% mutate(fdr1=sprintf("%f",f1))
    P <- P %>% mutate(fdr2=sprintf("%f",f2))
    P <- P %>% mutate(fdr3=sprintf("%f",f3))
  }
  
  name1 <- paste0(dir,"/simulate_hydi_data_counts_rep",replicates,"_cov",coverage,"_group1.txt")
  name2 <- paste0(dir,"/simulate_hydi_data_counts_rep",replicates,"_cov",coverage,"_group2.txt")
  name3 <- paste0(dir,"/simulate_hydi_data_counts_rep",replicates,"_cov",coverage,"_stats.txt")
  print(name1)
  print(name2)
  
  write_delim(Q1,name1, delim="\t")
  write_delim(Q2,name2, delim="\t")
  
  if(analyse) {
    write_delim(P,name3, delim="\t", col_names=F)
  }
  return(S)
}

simulate_foreground <- function() {
  load("nullmodel_grid_19.Rda")
  
  silent = F
  M <- 500
  set.seed(789)
  
  p1s <- c(0.5,0.5,0.5,0.5,0.5,   0.9,0.9,0.9,0.9,0.9,    0.1,0.1,0.1,0.1,0.1)
  q1s <- c(0.5,0.5,0.5,0.5,0.5,   0.9,0.9,0.9,0.9,0.9,    0.1,0.1,0.1,0.1,0.1)
  p2s <- c(0.5,0.5,0.5,0.5,0.5,   0.9,0.9,0.9,0.9,0.9,    0.1,0.1,0.1,0.1,0.1)
  q2s <- c(0.5,0.5,0.5,0.5,0.5,   0.9,0.9,0.9,0.9,0.9,    0.1,0.1,0.1,0.1,0.1)
  
  di <- c(0.1,0.15,0.2,0.25,0.3,  0.1,0.15,0.2,0.25,0.3, 0.1,0.15,0.2,0.25,0.3)
  
  R1 <- tibble()
  R_ <- tibble()
  
  for(i in 1:length(p1s)) {
    temp <-  Q %>% filter(p1==p1s[i], q1==q1s[i], p2==p2s[i], q2==p2s[i], coverage %in% c(30,40,100)) %>% 
      mutate(diff=0) %>% mutate(combination=i)
    R_ <- bind_rows(R_, temp)
  }
  
  #simulate Diff
  for(c in unique(R_$coverage)) {
    for(k in unique(R_$replicates)) {
      for(i in 1:length(unique(R_$combination))) {
        for(j in 1:M) {
          if(!silent && j %% 1000 == 1) {
            print(paste("c:",c,"k:",k,"i:",i,"j:",j))
          }
          if(p2s[i] + di[i] < 1 && p2s[i] + di[i] > 0) {
            S <- simulate_grid_script(p1s[i], q1s[i], p2s[i]+di[i], q2s[i], replicates=k, cov.mean.modC = c)
          } else if(q2s[i] - di[i] > 0 && q2s[i] - di[i] < 1) {
            S <- simulate_grid_script(p1s[i], q1s[i], p2s[i], q2s[i]-di[i], replicates=k, cov.mean.modC = c)
          } else {
            print("not allowed")
          }
          tab = hydi_test(S)
          R1 <- bind_rows(R1, tibble(coverage=c, replicates=k,combination=i,diff=1, sampleno=j, 
                                     p1=p1s[i], q1=q1s[i], p2=p2s[i], q2=q2s[i]-di[i], tab))
        }
      }
    }
  }
  R <- bind_rows(R_,R1)
  save(R,file="simulation_foreground_4.Rda")
}

runtime_benchmark <- function(N=10000) {
  for(reps in c(5,10,20)) {
    for(cov in c(30,40,100,200)) {
      print(paste("cov:",cov,"reps:",reps))
      simulate_hydi_data(N=N, replicates=reps, coverage=cov, analyse=F, dir="./runtime_memory_benchmarks")
    }
  }
}
