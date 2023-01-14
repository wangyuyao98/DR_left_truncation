#functions for analysis of censored data with P(C<T)>0

get.S <- function(t, time, surv)	# assumes that t is a scalar
{
  n <- length(time)
  if (t<time[1] && abs(t-time[1])>100*.Machine$double.eps) return(1)
  if (t>time[n] && abs(t-time[n])>100*.Machine$double.eps) return(surv[n])
  return( surv[sum(time<=t | abs(time-t)<100*.Machine$double.eps)] )
}

# inner function for CW, allows ties:
get.CW <- function(X, T, delta, Z)
{
  n <- length(T)
  m <- max(c(X, T))
  d.t <- rep(1,n)
  sol <- try(coxph(Surv(time=m-X, time2=m-T, event=d.t) ~ Z) ) 
  if (class(sol)!="coxph") 
  {
    cat("coxph: Error occurred in coxph(srv.trunc ~ Z).","\n")
  }
  alpha <- matrix(sol$coef, length(sol$coef), 1)
# cat("CW alpha=", alpha, "\n") 
# cat("SE alpha=", sqrt(diag(sol$var)), "\n") 

  bh <- basehaz(sol, centered=FALSE)
  H <- bh$hazard # ordered cumulative hazard

#---------------  
#   # diagnostic plots for GOF of the Cox model vs nonparametric estimates:
#   # Z_1:
#   fit.0 <- summary(survfit(Surv(time=m-X[Z[,1]==0], time2=m-T[Z[,1]==0],
#                                 event=delta[Z[,1]==0])~1))
#   fit.1 <- summary(survfit(Surv(time=m-X[Z[,1]==1], time2=m-T[Z[,1]==1],
#                                 event=delta[Z[,1]==1])~1))
#   t1.r <- c(0, fit.1$time)
#   S1.r <- c(1, fit.1$surv)
#   t0.r <- c(0, fit.0$time)
#   S0.r <- c(1, fit.0$surv)
#   # diagnostic plots:
#   plot(t1.r, log(-log(S1.r)), type="s", main="check PH for reverse time-Z_1",
#        xlim=c( min(c(t1.r,t0.r)), max(c(t1.r,t0.r)) )  )  
#   #  ylim=c( min(c(log(-log(S1.r[S1.r!=0])), log(-log(S0.r[S0.r!=0])))), max(c(log(-log(S1.r[S1.r!=0])), log(-log(S0.r[S0.r!=0])))) ) ) 
#   lines(t0.r, log(-log(S0.r)), type="s", col="purple")
#   
#   # Z_2:
#   fit.0 <- summary(survfit(Surv(time=m-X[Z[,2]==0], time2=m-T[Z[,2]==0],
#                                 event=delta[Z[,2]==0])~1))
#   fit.1 <- summary(survfit(Surv(time=m-X[Z[,2]==1], time2=m-T[Z[,2]==1],
#                                 event=delta[Z[,2]==1])~1))
#   
#   t1.r <- c(0, fit.1$time)
#   cat(t1.r,"\n")
#   S1.r <- c(1, fit.1$surv)
#   t0.r <- c(0, fit.0$time)
#   S0.r <- c(1, fit.0$surv)
#   
#   plot(t1.r, log(-log(S1.r)), type="s", main="check PH for reverse time-Z_2",
#        xlim=c( min(c(t1.r,t0.r)), max(c(t1.r,t0.r)) )  )  
# #  ylim=c( min(c(log(-log(S1.r[S1.r!=0])), log(-log(S0.r[S0.r!=0])))), max(c(log(-log(S1.r[S1.r!=0])), log(-log(S0.r[S0.r!=0])))) ) ) 
#   lines(t0.r, log(-log(S0.r)), type="s", col="purple")
#   
  #---------------------
  
  K.S.NA <- exp(-H) # in the order of m-T[n], m-T[n-1], ...
  n1 <- length(bh$hazard) 

  K.CDF.forw <- K.S.NA[n1:1] # this is only baseline pr-ties for z=0
  base <- 1/(1-sapply(X, get.S, time=m-bh$time[n1:1], surv=1-K.CDF.forw, simplify=TRUE, USE.NAMES=FALSE))

  s <- as.vector(exp(Z %*% alpha ) )
  w <- base^s #

  f <- survfit(Surv(X, delta) ~ 1, weights=w) # survfit treats the ties in X correctly, but smth was wrong there in the length times
  list(time=summary(f)$time, surv=summary(f)$surv)
}

# CW estimator - the function allows for having several covariates in Z:
CW.deptr <- function(X, T, delta, Z, bs=FALSE, nbs.rep=200) 
{
# Z is a matrix (n X p) or a data frame
  q <- dim(Z)
  n <- length(T)
  Z <- as.matrix(Z, q[1], q[2])
  # removing missing observations:
  na.i <- is.na(X) | is.na(T) | is.na(delta) | apply(is.na(Z), 1, any)
  if (sum(na.i)>0)
  {
    cat(sum(na.i), " missing observations were omitted.\n")
    X <- X[!na.i]
    T <- T[!na.i]
    delta <- delta[!na.i]
    Z <- Z[!na.i,]
  }
  est.CW <- get.CW(X, T, delta, Z)

# bootstrap:
# estimate SE at these time points:
  x=est.CW$time
  
  se.bs=cil=ciu=NULL
  
  if (bs) {
    b.b <-  NULL
    for (b in 1:nbs.rep) 
    {
      samp.b <- sample(n, size=n, replace=TRUE)
      bs.est <- get.CW(X[samp.b], T[samp.b], delta[samp.b], Z[samp.b,])
      res.bs <- sapply(x, get.S, time=bs.est$time, surv=bs.est$surv, simplify=TRUE, USE.NAMES=FALSE)
      b.b <- rbind(b.b, res.bs)
    }
    se.bs <- apply(b.b, 2, sd, na.rm=TRUE)
    cil <- apply(b.b, 2, quantile, prob=0.025, na.rm=TRUE)
    ciu <- apply(b.b, 2, quantile, prob=0.975, na.rm=TRUE)
    avg <- apply(b.b, 2, mean, na.rm=TRUE)
  }
  list(time=est.CW$time, surv=est.CW$surv, SE=se.bs, CI.L=cil, CI.U=ciu)
}




### the TVW estimator, the function that allows for several covariates in Z:
### function for the size of at-risk set
m.fun<-function(surv, trun, cens)
{
  a0=surv[cens>0]
  a1=sort(a0)
  n0=length(surv)
  m0=length(a1)
  xx=rep(surv, m0)
  tt=rep(trun, m0)
  xx2=rep(a1, each=n0)
  yy=(xx>=xx2)*(tt<=xx2)
  yy.mat=matrix(yy, nrow=n0, ncol=m0, byrow=F)
  risk.set=apply(yy.mat, 2, sum)
  return(risk.set)  #output in order of survival times
}

###survival function under indep truncation
surv.indep.trun<-function(surv, trun, cens)
{
  a0=surv[cens>0] #number of failure times
  a1=sort(a0)
  n0=length(surv)
  m0=length(a1)
  vec.pl= 1-1/m.fun(surv, trun, cens)
  surv.est=cumprod(vec.pl)
  return(surv.est)
}

### cox model (regression with reversed truncation time)
# no censoring here!
cox.reg.coef.fun=function(surv, trun, cens, zz)
{  ############### input ################################
  #surv: original observed survival time Y
  #trun: original truncation time T
  #cens: censoring indicator, 1-failure; 0-censored.
  #zz: covariate matrix
  #####################################################
  #
  #reverse T and Y
  tau.y=max(surv, trun)+0.1
  surv.rev=tau.y-surv
  trun.rev=tau.y-trun
  cens2=rep(1, length(trun.rev))
  cox.fit=coxph(Surv(time=surv.rev, time2=trun.rev, event=cens2) ~ zz)
  reg.coef=as.vector(as.numeric(cox.fit$coef)) #estimate of reg coef
  return(reg.coef)
}

###Breslow estimator of baseline hazard function
# output is a vector with length n0 (# of distinctive truncation times)
breslow.fun=function(surv, trun, cens, zz)
{
  #reverse T and Y
  tau.y=max(surv, trun)+0.1
  surv.rev=tau.y-surv
  trun.rev=tau.y-trun
  exp.z=exp(zz%*%cox.reg.coef.fun(surv, trun, cens, zz))
  exp.z.neg=exp(-zz%*%cox.reg.coef.fun(surv, trun, cens, zz))
  #
  ### calculate at risk set
  a1=sort(trun.rev) #sort T_i reverse
  n0=length(trun.rev)
  xx=rep(trun.rev, n0)
  tt=rep(surv.rev, n0)
  xx2=rep(a1, each=n0)
  yy=(xx>=xx2)*(tt<=xx2)
  yy.exp.z=yy*rep(exp.z, n0)
  yy.exp.z.mat=matrix(yy.exp.z, nrow=n0, ncol=n0, byrow=F)
  risk.set=apply(yy.exp.z.mat, 2, sum)
  #
  exp.z.sort=exp.z[order(trun.rev)] #revised on 04-09-2012!!!! (should have the same order as of the failure time)
  exp.z.sort.neg=exp.z.neg[order(trun.rev)] 
  vec.bres=1 - (1-exp.z.sort/risk.set)^(exp.z.sort.neg)
  breslow.culm=cumsum(vec.bres) #length of n0
  ###jump of breslow
  breslow.diff= breslow.culm - c(0, breslow.culm[-length(breslow.culm)])
  breslow.est=list(culm=breslow.culm, jump=breslow.diff)
  return(breslow.est)
}

###survival function under dep truncation
surv.dep.trun<-function(surv, trun, cens, zz, bs=FALSE, nbs.rep=200)
{
  ### estimation of K_i(t)
  #output: a matrix, the row is given covariate, 
  #the column is given time point t
  #note that the model is for truncation time
  zz=as.matrix(zz)
  #tau.y=max(surv, trun)+0.1
  tau.y=0
  trun.rev=tau.y-trun
  a0=surv[cens>0] #number of failure times
  a1=sort(a0)
  t.rev=tau.y-a1
  n0=length(surv)
  #n1=length(trun.rev) n1==n0
  m0=length(a1)
  
  xx=rep(surv, m0)
  tt=rep(trun, m0)
  tt.rev=rep(trun.rev, m0)
  xx2=rep(a1, each=n0)
  xx2.rev=rep(t.rev, each=n0)
  
  yy=(xx>=xx2)*(tt<=xx2)
  yy.rev=(tt.rev<=xx2.rev)  #at risk set for ipw calculation
  yy.rev.mat=matrix(yy.rev, nrow=n0, ncol=m0, byrow=F)
  risk.set.ipw=apply(yy.rev.mat, 2, sum) #(r_1, r_2, ..., R_m), length of m0
  #dealing with risk.set.ipw==0
  risk.set.ipw.2=risk.set.ipw+(risk.set.ipw==0)*(n0+1)
  
  alpha=cox.reg.coef.fun(surv, trun, cens, zz)
  cat("TVW alpha=", alpha,"\n")    
  
  exp.z=exp(zz%*%as.matrix(alpha)) # n0*1 matrix
  breslow.jump=breslow.fun(surv, trun, cens, zz)$jump  #length of n0    ############### revised on 04-09-2012
  A.mat=1- exp.z%*%t(as.matrix(breslow.jump))
  A.mat.2=((A.mat>0)*1)*A.mat + ((A.mat<=0)*1) # deal with negative values
  B.mat=t(apply(A.mat.2, 1, cumprod)) #n0*n0 matrix, each row corresponds to each z_i
  B.mat.2=cbind(B.mat, rep(1,n0)) #n0*(n0+1) matrix, dealing with risk.set.ipw==0
  #K_j(Y_(i)) is the (j, r_i) element of B.mat, if (risk.set.ipw==0) then K(.)=1
  #
  #now, the estimator of S_X(x)
  #first, calculate element of the product for each ordered failure time (totally m0)
  ##denominator
  ipw.mat=B.mat.2[ , risk.set.ipw.2] #n0*m0 matrix
  yy.mat=matrix(yy, nrow=n0, ncol=m0, byrow=F) #n0*m0 matrix
  yy.mat.2=yy.mat/ipw.mat #elementwise dividing
  risk.set=apply(yy.mat.2, 2, sum)  #m0*1 vector (n_i in asymptotic varaice formula)
  ##numerator
  #index of uncensored obs
  index.mat=cbind(surv, cens, (1:n0))
  index.mat.2=index.mat[cens>0,]
  surv.uncen=index.mat.2[,1]
  index.vec=index.mat.2[order(surv.uncen),3] #original position of Y_(1),...,Y_(m0)
  ipw.mat.num=B.mat.2[index.vec, risk.set.ipw.2] #m0*m0 matrix for numerator
  ipw.vec.num=diag(ipw.mat.num) #m0*1 vector for numerator, K_(i)(Y_(i)), sorted by Y_(i)
  ##element of the product in eqn (2) [ref: 05-12-2011 manuscript]
  vec.pl= 1-(1/ipw.vec.num)/risk.set
  risk.set.2= risk.set^2/apply(yy.mat.2^2, 2, sum) #computing m_i in asymptotic variance formula
  var.est.comp= (1/ipw.vec.num)/(risk.set.2*(risk.set - 1/ipw.vec.num))
  surv.est=cumprod(vec.pl)
  var.est=(surv.est)^2*cumsum(var.est.comp)
  se.est=var.est^0.5

  # bootstrap..........................
  x = sort(surv[cens>0]) # estimate SE at these time points:
  cil.bs=ciu.bs=se.bs=NULL
  n <- length(surv)
  if (bs) 
  {
    b.b <-  NULL
    for (b in 1:nbs.rep) 
    {
      samp.b <- sample(n, size=n, replace=TRUE)
      zz.bs=as.matrix(zz[samp.b,])
      surv.bs =surv[samp.b]
      cens.bs <- cens[samp.b]
      trun.bs <- trun[samp.b]
      #tau.y=max(surv, trun)+0.1
      tau.y.bs=0
      trun.rev.bs=tau.y.bs-trun.bs
      a0.bs=surv.bs[cens.bs>0] #number of failure times
      a1.bs=sort(a0.bs)
      t.rev.bs=tau.y.bs-a1.bs
      n0.bs=length(surv.bs)
      #n1=length(trun.rev) n1==n0
      m0.bs=length(a1.bs)
      
      xx.bs=rep(surv.bs, m0.bs)
      tt.bs=rep(trun.bs, m0.bs)
      tt.rev.bs=rep(trun.rev.bs, m0.bs)
      xx2.bs=rep(a1.bs, each=n0.bs)
      xx2.rev.bs=rep(t.rev.bs, each=n0.bs)
      
      yy.bs=(xx.bs>=xx2.bs)*(tt.bs<=xx2.bs)
      yy.rev.bs=(tt.rev.bs<=xx2.rev.bs)  #at risk set for ipw calculation
      yy.rev.mat.bs=matrix(yy.rev.bs, nrow=n0.bs, ncol=m0.bs, byrow=F)
      risk.set.ipw.bs=apply(yy.rev.mat.bs, 2, sum) #(r_1, r_2, ..., R_m), length of m0
      #dealing with risk.set.ipw==0
      risk.set.ipw.2.bs=risk.set.ipw.bs+(risk.set.ipw.bs==0)*(n0.bs+1)
      
      exp.z.bs=exp(zz.bs%*%as.matrix(cox.reg.coef.fun(surv.bs, trun.bs, cens.bs, zz.bs))) # n0*1 matrix
      breslow.jump.bs=breslow.fun(surv.bs, trun.bs, cens.bs, zz.bs)$jump  #length of n0    ############### revised on 04-09-2012
      A.mat.bs=1- exp.z.bs%*%t(as.matrix(breslow.jump.bs))
      A.mat.2.bs=((A.mat.bs>0)*1)*A.mat.bs + ((A.mat.bs<=0)*1) # deal with negative values
      B.mat.bs=t(apply(A.mat.2.bs, 1, cumprod)) #n0*n0 matrix, each row corresponds to each z_i
      B.mat.2.bs=cbind(B.mat.bs, rep(1,n0.bs)) #n0*(n0+1) matrix, dealing with risk.set.ipw==0
      #K_j(Y_(i)) is the (j, r_i) element of B.mat, if (risk.set.ipw==0) then K(.)=1
      #
      #now, the estimator of S_X(x)
      #first, calculate element of the product for each ordered failure time (totally m0)
      ##denominator
      ipw.mat.bs=B.mat.2.bs[ , risk.set.ipw.2.bs] #n0*m0 matrix
      yy.mat.bs=matrix(yy.bs, nrow=n0.bs, ncol=m0.bs, byrow=F) #n0*m0 matrix
      yy.mat.2.bs=yy.mat.bs/ipw.mat.bs #elementwise dividing
      risk.set.bs=apply(yy.mat.2.bs, 2, sum)  #m0*1 vector (n_i in asymptotic variance formula)
      ##numerator
      #index of uncensored obs
      index.mat.bs=cbind(surv.bs, cens.bs, (1:n0.bs))
      index.mat.2.bs=index.mat.bs[cens.bs>0,]
      surv.uncen.bs=index.mat.2.bs[,1]
      index.vec.bs=index.mat.2.bs[order(surv.uncen.bs),3] #original position of Y_(1),...,Y_(m0)
      ipw.mat.num.bs=B.mat.2.bs[index.vec.bs, risk.set.ipw.2.bs] #m0*m0 matrix for numerator
      ipw.vec.num.bs=diag(ipw.mat.num.bs) #m0*1 vector for numerator, K_(i)(Y_(i)), sorted by Y_(i)
      ##element of the product in eqn (2) [ref: 05-12-2011 manuscript]
      vec.pl.bs= 1-(1/ipw.vec.num.bs)/risk.set.bs
      risk.set.2.bs= risk.set.bs^2/apply(yy.mat.2.bs^2, 2, sum) #computing m_i in asymptotic varaice formula
      var.est.comp.bs= (1/ipw.vec.num.bs)/(risk.set.2.bs*(risk.set.bs - 1/ipw.vec.num.bs))
      surv.est.bs=cumprod(vec.pl.bs)
      
      res.bs <- sapply(x, get.S, time=sort(surv.bs[cens.bs>0]), surv=surv.est.bs, simplify=TRUE, USE.NAMES=FALSE)
      b.b <- rbind(b.b, res.bs)
    }
    
    se.bs <- apply(b.b, 2, sd, na.rm=TRUE)
    cil.bs <- apply(b.b, 2, quantile, prob=0.025, na.rm=TRUE)
    ciu.bs <- apply(b.b, 2, quantile, prob=0.975, na.rm=TRUE)
  }
  
  results=list(times=x, surv.est = surv.est, se.est=se.est, risk.set=risk.set, di=1/ipw.vec.num,
               SE.bs=se.bs, CI.L.bs=cil.bs, CI.U.bs=ciu.bs)
}


