# function for simulating data from different scenarios

gen <- function(n, multi = 100, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                        beta0.T, beta0.Q, beta.T = NA, beta.Q = NA, beta.T2 = NA, 
                        beta.Q2 = NA,
                        shape.T = NA, shape.Q = NA, epsT.sd = NA, epsQ.sd = NA){
    
    Z1 = runif(multi*n, -1, 1)
    Z2 = rbinom(multi*n, size = 1, prob = 0.5) - 0.5
    Z = cbind(Z1 = Z1, Z2 = Z2)
    Zsq = cbind(Z1 = Z1, Z2 = Z2, Z1_sq = Z1^2-1/3, Z1Z2 = Z1*Z2)
    # Zsq = cbind(Z1 = Z1, Z2 = Z2, Z1_sq = sqrt(abs(Z1))-2/3, Z1Z2 = Z1*Z2)
    
    id = which(Z2<0)
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
    
    
    Q2.max = tau - Q.min
    Q2.min = tau - Q.max
    U2 = runif(multi*n, min = 0, max = 1)
    if(Qmodel == "Cox1"){
        Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-Z %*% beta.Q))))
    }else if(Qmodel == "Cox2"){
        Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-Zsq %*% beta.Q2))))
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
    dat.full = dat.full[(1:obs.id[n]), ]
    
    return(list(dat = dat, dat.full = dat.full, dat.all = dat.all))
}



# function for simulating data from mixture models

gen.mixture <- function(n, multi = 100, Tmodel, Qmodel, T.min, Q.min, Q.max, tau,
                beta0.T_1, beta0.T_2, beta0.Q_1, beta0.Q_2,
                beta.T_1, beta.T_2, beta.Q_1, beta.Q_2, 
                shape.T = NA, shape.Q = NA, epsQ.min, epsQ.max){
    
    Z1 = runif(multi*n, -1, 1)
    Z2 = rbinom(multi*n, size = 1, prob = 0.5) - 0.5
    Z = cbind(Z1 = Z1, Z2 = Z2)
    Zsq = cbind(Z1 = Z1, Z1_sq = Z1^2)
    
    id = which(Z2<0)
    nn1 = length(id)
    nn2 = multi*n - nn1
    TT = rep(NA, multi*n)
    Q2 = rep(NA, multi*n)
    
    if(Tmodel == 'AFT2_weibull2'){
        TT[id] = T.min + exp(beta0.T_1 + Zsq[id,]%*%beta.T_1 + runif(nn1,-1,1))
        TT[-id] = T.min + rweibull(nn2, shape = shape.T, scale = exp(-1/shape.T * (beta0.T_2 + (Zsq[-id,]%*%beta.T_2))))
    }else if(Tmodel == 'weibull2_AFT2'){
        TT[id] = T.min + rweibull(nn1, shape = shape.T, scale = exp(-1/shape.T * (beta0.T_1 + (Zsq[id,]%*%beta.T_1))))
        TT[-id] = T.min + exp(beta0.T_2 + Zsq[-id,]%*%beta.T_2 + runif(nn2,-1,1))
    }
    
    
    Q2.max = tau - Q.min
    Q2.min = tau - Q.max
    U2 = runif(multi*n, min = 0, max = 1)
    z1 = seq(-1,1, by = 0.01)
    z1_sq = cbind(z1, z1^2)
    if(Qmodel == "Cox2_AFT2"){
        Q2[id] = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2[id]^(exp(-(beta0.Q_1+Zsq[id,] %*% beta.Q_1)))))
        Q2[-id] = Q2.min + exp(beta0.Q_2 + Zsq[-id,]%*%beta.Q_2 + runif(nn2, epsQ.min, epsQ.max))
        
        # check if the simulated Q2 is less that Q2.max
        if(Q2.min + exp(beta0.Q_2 + max(z1_sq%*%beta.Q_2) + epsQ.max) > Q2.max){
            stop("Q2.max exceed tau_2")
        }
        
    }else if(Qmodel == "AFT2_Cox2"){
        Q2[id] = Q2.min + exp(beta0.Q_1 + Zsq[-id,]%*%beta.T_1 + runif(nn1, epsQ.min, epsQ.max))
        Q2[-id] = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2[-id]^(exp(-(beta0.Q_2+Zsq[id,] %*% beta.Q_2)))))
        
        # check if the simulated Q2 is less that Q2.max
        if(Q2.min + exp(beta0.Q_1 + max(z1_sq%*%beta.Q_1) + epsQ.max) > Q2.max){
            stop("Q2.max exceed tau_2")
        }
    }else if(Qmodel == "Cox2_step2"){
        Q2[id] = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2[id]^(exp(-(beta0.Q_1+Zsq[id,] %*% beta.Q_1)))))
        Q2[-id] = tau - (Q.min + (Q.max - Q.min)*(Z1[-id]^2)*rbeta(nn2, 2, 1))
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
    dat.full = dat.full[(1:obs.id[n]), ]
    
    return(list(dat = dat, dat.full = dat.full, dat.all = dat.all))
}

