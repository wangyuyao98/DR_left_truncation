# function for simulating data from different scenarios

gen.mixture <- function(n, multi = 100, Tmodel, Qmodel, T.min, Q.max, tau,
                        beta0.T, beta0.Q, beta.T = NA, beta.Q = NA, beta.T2 = NA, beta.Q2 = NA,
                        shape.T = NA, shape.Q = NA, epsT.sd = NA, epsQ.sd = NA){
    Q2.min = tau - Q.max
    
    Z1 = runif(multi*n, -1, 1)
    Z2 = rbinom(multi*n, size = 1, prob = 0.5) - 0.5
    Z = cbind(Z1 = Z1, Z2 = Z2)
    Zsq = cbind(Z1 = Z1, Z2 = Z2, Z1_sq = Z1^2-1/3, Z1Z2 = Z1*Z2)
    id = which(Z%*%beta.T<0)
    nn1 = length(id)
    nn2 = multi*n - nn1
    TT = rep(NA, multi*n)
    Q2 = rep(NA, multi*n)
    
    if(Tmodel == 'weibull1'){
        TT = T.min + rweibull(multi*n, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + Z%*%beta.T)))
    }else if(Tmodel == 'weibull2'){
        TT = T.min + rweibull(multi*n, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + Zsq%*%beta.T2)))
    }else if(Tmodel == 'AFT2_weibull2'){
        TT[id] = T.min + exp(beta0.T + (Zsq%*%beta.T2)[id] + rnorm(nn1, mean = 0, sd = epsT.sd))
        TT[-id] = T.min + rweibull(nn2, shape = shape.T, scale = exp(-1/shape.T * (beta0.T + (Zsq%*%beta.T2)[-id])))
    }
    
    if(Qmodel =='weibull1'){
        Q2 = Q2.min + rweibull(multi*n, shape = shape.Q, scale = exp(-1/shape.Q * (beta0.Q + Z%*%beta.Q)))
    }else if(Qmodel == 'weibull2'){
        Q2 = Q2.min + rweibull(multi*n, shape = shape.Q, scale = exp(-1/shape.Q * (beta0.Q + Zsq%*%beta.Q2)))
    }else if(Qmodel == 'weibull2_AFT2'){
        Q2[id] = Q2.min + rweibull(nn1, shape = shape.Q, scale = exp(-1/shape.Q * (beta0.Q + (Zsq%*%beta.Q2)[id])))
        # Q2[-id] = Q2.min + exp(beta0.Q + (Zsq%*%beta.Q2)[-id] + rnorm(nn2, mean = 0, sd = epsQ.sd) )   # this will lead to more Q less than 0
        Q2[-id] = Q2.min + exp(beta0.Q + (Zsq%*%beta.Q2)[-id] + 1 * rweibull(nn2,shape=1.5,scale=1)-gamma(1+1/1.5))
    }
    
    QQ = tau - Q2
    
    dat.full = data.frame(time = TT, Q = QQ, Z)
    dat.all = dat.full
    obs.id = which(dat.full$Q < dat.full$time)
    dat.obs = dat.full[obs.id,]
    if(length(obs.id)<n){
        stop('Truncation rate is high. Need to increase the mutiplier.')
    }
    dat = dat.obs[1:n,]
    
    # bound the left truncation and event times at 0 and tau respectively
    dat$time = pmin(tau, dat$time)
    dat$Q = pmax(0, dat$Q)
    # indicators for whether the times are in [0,tau] or not (1 or 0)
    delta.T = as.numeric(dat$time < tau)
    delta.Q = as.numeric(dat$Q > 0)
    
    dat = cbind(delta.T = delta.T, delta.Q = delta.Q, dat)
    dat.full = dat.full[(1:obs.id[n]), ]
    
    return(list(dat = dat, dat.full = dat.full, dat.all = dat.all))
}
