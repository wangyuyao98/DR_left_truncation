## functions used in applying random forest for LTRC data 
## from LTRCforests package


# return a matrix m(v,Z_i):= \int_0^v\nu(u)dF(u)
# row - subjects; col - each v value
int_nu_dF.mx.LTRCforests <- function(T.pred, v, nu, delta, Sc){
    time = T.pred$survival.times
    Ftz.mx = 1 - t(T.pred$survival.probs)  # row - subject; col - different times points
    dFtz.mx = Ftz.mx - cbind(0, Ftz.mx[,-ncol(Ftz.mx)]) 
    
    inte = sapply(v, int_nu_dF.mx.LTRCforests.single, 
                  nu = nu, time = time, dFtz.mx = dFtz.mx, delta = delta, Sc = Sc)

    return(inte)
}


# return a vector of integrals for each subject from 0 to v (a single time point)
int_nu_dF.mx.LTRCforests.single <- function(v, nu, time, dFtz.mx, delta, Sc){
    n.est = nrow(dFtz.mx)
    delta.mx = matrix( rep(delta, length(time)), nrow = n.est, byrow = FALSE)   # each row represent a subject
    nutc.mx = matrix(rep(nu(time)/Sc(time), n.est), nrow = n.est, byrow = TRUE)
    nut = delta.mx * nutc.mx
    vals = cbind(0, nut*dFtz.mx)
    
    id = findInterval(v, c(0, time, Inf))
    if(id == 1){
        inte = vals[,1]
    }else{
        inte = rowSums(vals[,(1:id)])
    }

    return(inte)
}




# return a vector m(v_i,Z_i):= \int_0^{v_i}\nu(u)dF(u|Z_i)
# row - subjects; col - each v value
int_nu_dF.LTRCforests <- function(T.pred, v, nu, delta, Sc){
    time = T.pred$survival.times
    Ftz.mx = 1 - t(T.pred$survival.probs)  # row - subject; col - different times points
    dFtz.mx = Ftz.mx - cbind(0, Ftz.mx[,-ncol(Ftz.mx)]) 
    n.est = nrow(Ftz.mx)
    m = length(v)
    if(m!=n.est){
        stop("The number of time points does not equal to the number of subjects!")
    }
    if(n.est!=length(delta)){
        stop("The length of delta does not equal to the number of subjects!")
    }
    
    inte = rep(NA, n.est)
    for(i in 1:n.est){
        dF = dFtz.mx[i,] 
        vals = c(0,nu(time)/Sc(time)*dF)
        
        id = findInterval(v[i], c(0, time, Inf))
        inte[i] = delta[i]*sum(vals[(1:id)])
    }
    
    return(inte)
}



# ## compute the integral \int F(v|Z)/((G(v|Z))^2*(1-F(v|Z))) J(v) dG(v|Z)
# # G.jumps: jump times of G
# int_FJ_dG.LTRCforests <- function(Fvz.mx, Gvz.mx, atrisk, tau, nu, trim){
#     
#     dGvz.mx = Gvz.mx - cbind(0, Gvz.mx[,-ncol(Gvz.mx)])
#     
#     result = Fvz.mx/((pmax(Gvz.mx,trim))^2*pmax(1-Fvz.mx, trim)) * atrisk * dGvz.mx
#     inte = rowSums(result)
#     
#     return(inte)
# }

## compute the integral \int m(v,Z;F)/((G(v|Z))^2*(1-F(v|Z))) J(v) dG(v|Z)
int_JG_m.LTRCforests <- function(T.pred.dF,jumps.Q, Fvz.mx, Gvz.mx, atrisk, 
                                 tau, nu, trim, delta, Sc){
    
    mvz.mx = int_nu_dF.mx.LTRCforests(T.pred.dF, jumps.Q, nu, delta, Sc)
    dGvz.mx = Gvz.mx - cbind(0, Gvz.mx[,-ncol(Gvz.mx)])

    result = mvz.mx/((pmax(Gvz.mx,trim))^2*pmax(1-Fvz.mx,trim)) * atrisk * dGvz.mx
    inte = rowSums(result)
    
    return(inte)
}



# function that return the predicted the survival probabilities for all subjects
# at one given time points
predict.survProb.mx.single <- function(newtime, predict.obj){
    time = predict.obj$survival.times
    surv = t(predict.obj$survival.probs)   # row - subject, col - time points
    
    if(newtime < min(time)){
        surv.predict = rep(1, nrow(surv))
    }else{
        id = max(which(time <= newtime))
        surv.predict = surv[,id]
    }
    
    return(surv.predict)
}



# function that predict the survival probabilities at a vector of given time points
predict.survProb.mx <- function(predict.obj, newtime){
    surv.predict = sapply(newtime, predict.survProb.mx.single, 
                          predict.obj = predict.obj)
    
    return(surv.predict)
}



predictProb.vec <- function(newdata.teval, object){
    teval = newdata.teval$time.eval
    surv = predictProb(object = object, newdata = newdat.teval,
                       time.eval = teval)$survival.probs
    
    return(surv)
}


