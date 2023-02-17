
# package
library(glmnet)
library(mvtnorm)
library(MASS)
library(rmutil)
library(nleqslv)
library(doParallel)
library(rootSolve)
library(parallel) #
library(data.table)
library(plyr)
library(SGL)
library(graper)
library(gcdnet)
library(fwelnet)
library(gtools)
library(grpreg)
# library(gtools)
# library(ipflasso)

mycombine_rep_delta = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  length.list.out = length(list.combined[[1]])
  list.out = vector('list', length.list.out)
  for(k in 1:length.list.combined){
    
    for(l in 1:(length.list.out-2)){
      
      list.out[[l]] = rbind(list.out[[l]], list.combined[[k]][[l]])
    }
    
    for(l in (length.list.out-1):length.list.out){
      
      list.out[[l]] = abind::abind(list.out[[l]], 
                                   list.combined[[k]][[l]],
                                   rev.along = 3)
    }
  }
  
  list.out
}

mycombine_rep = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  length.list.out = length(list.combined[[1]])
  list.out = vector('list', length.list.out)
  for(k in 1:length.list.combined){
    
    for(l in 1:(length.list.out-2)){
      
      if(k==1) list.out[[l]] = vector('list', 2)
      list.out[[l]][[1]] = rbind(list.out[[l]][[1]], list.combined[[k]][[l]][[1]])
      list.out[[l]][[2]] = rbind(list.out[[l]][[2]], list.combined[[k]][[l]][[2]])
    }
    
    for(l in (length.list.out-1):length.list.out){
      
      if(k==1) list.out[[l]] = vector('list', 2)
      list.out[[l]][[1]] = abind::abind(list.out[[l]][[1]],
                                        list.combined[[k]][[l]][[1]],
                                        rev.along = 3)
      list.out[[l]][[2]] = abind::abind(list.out[[l]][[2]],
                                        list.combined[[k]][[l]][[2]],
                                        rev.along = 3)
    }
  }
  
  list.out
}

# generating design matrix
design = function(type = "iid", n, p){
  if(type=="iid"){
    
    X = matrix(rnorm((n*p), sd = sqrt(1/n) ),n,p)
    
  }else if(type=="ar") {
    
    Sigma = matrix(NA,p,p)
    for(i in 1:p) Sigma[i,] = 0.5^(abs(i-(1:p)))
    X = matrix(rnorm(n*p),n,p)
    X = t(t(chol(Sigma))%*%t(X))
    # X = scale(X)
    
  }else if(type=="eq") {
    
    Sigma = matrix(0.5, p, p)
    diag(Sigma) = 1
    X = matrix(rnorm(n*p),n,p)
    X = t(t(chol(Sigma))%*%t(X))
    # X = scale(X)
    
  }else if(type=="binary"){
    
    X = matrix(sample(x=c(-1,1), size = n*p, replace = T, prob = c(.5,.5)),
               n,p)
    X = X/sqrt(n)
  }
  
  return(X)
}


# generating true \beta modeling through covariate
gen.coeff = function(n = 1, pi0, x, nSignal, l = -1, u = 1, a = 0, b = 1){
  
  # probability of being 0
  if(missing(pi0)&&(!missing(x))) pi0 = 1/(1 + exp(-(a + b*x)))
  
  if(missing(nSignal)){
    
    # indicators of being 0 or not. 1 indicates a signal
    ind = rbinom(n, 1, 1-pi0)
    
    # generating nonzero signals
    if(sum(ind)==0){
      
      return(numeric(n))
      
    }else{
      
      out = numeric(n)
      out[which(ind==1)] = runif(sum(ind), l, u)*(2*rbinom(sum(ind), 1, 0.5) -1)
      
      return(out)
    }
    
  }else if(nSignal==0){
    
    return(numeric(n))
    
  }else{
    
    # indicators of being 0 or not. 1 indicates a signal
    ind = rbinom(n, 1, 1-pi0)
    
    # ensuring nSignal many signals
    if(sum(ind)>nSignal){
      
      turnsignaloff = sample(x = which(ind==1), 
                             size = sum(ind)-nSignal)
      ind[turnsignaloff] = 0
      
    }else if(sum(ind)<nSignal){
      
      turnsignalon = sample(x = which(ind==0), 
                            size = nSignal-sum(ind))
      ind[turnsignalon] = 1
    }
    
    # generating nonzero signals
    out = numeric(n)
    out[which(ind==1)] = runif(nSignal, l, u)*(2*rbinom(nSignal, 1, 0.5) -1)
    
    return(out)
  }
}


# soft thresholding function
eta = function(x, theta) {
  
  t = sign(x)*max((abs(x)-theta),0)
  return(t)
}


# # objective function for covariate structure
# obj.func = function( b, beta.init, cov.vec, gam=1, root=0){
# 
#   w = exp(b*cov.vec)
#   if(gam==1){
#     sum((w*abs(beta.init)) - log(w)) -root
#   }else{
#     sum((w*abs(beta.init)) - w^(1- (1/gam))/(1- (1/gam)) ) -root
#   }
# }

# # solution of tau1 in cov structure for univariate u
# f = function(tau1, beta.v, cov.v, gam = 1){
#   
#   if(gam==1){
#     
#     return(sum(exp(tau1*cov.v[beta.v!=0])*cov.v[beta.v!=0]*abs(beta.v[beta.v!=0]))/
#              sum(exp(tau1*cov.v[beta.v!=0])*abs(beta.v[beta.v!=0])))
#     
#   }else{
#     
#     return((sum(exp(tau1*cov.v[beta.v!=0])*cov.v[beta.v!=0]*abs(beta.v[beta.v!=0]))/
#               sum(exp(tau1*cov.v[beta.v!=0])*abs(beta.v[beta.v!=0]))) -
#              (sum(exp((1-(1/gam))*tau1*cov.v)*cov.v)/
#                 sum(exp((1-(1/gam))*tau1*cov.v))))
#   }
# }


# Lambdas <- function(...) {
#   cv <- cv.glmnet(...)
#   return(data.table(cvm=cv$cvm, lambda=cv$lambda))
# }
# 
# Optimcvm <- function(k, nCore, ...) {
#   
#   MSEs <- data.table(rbind.fill(mclapply(seq(k), function(dummy) Lambdas(...),
#                                          mc.cores = nCore)))
#   return(MSEs[, list(mean.cvm=mean(cvm)), lambda][order(mean.cvm)][1]$mean.cvm)
# }
# 
# OptimLambda <- function(k, nCore, ...) {
#   
#   MSEs <- data.table(rbind.fill(mclapply(seq(k), function(dummy) Lambdas(...),
#                                          mc.cores = nCore)))
#   return(MSEs[, list(mean.cvm=mean(cvm)), lambda][order(mean.cvm)][1]$lambda)
# }


# SA.Lasso = function(resp, desg, struct, lambda.fit = .1, gamma.fit = 1,
#                     Intercept = F, std.desg = T, std.resp = F,
#                     iter.max = 100, eps = 1e-4, ...){
#   
#   if(missing(struct)==T){
#     
#     n.fit = length(resp)  # n
#     p.fit = ncol(desg)  # p
#     nlambda = length(lambda.fit)
#     iter.max = iter.max +1  # max iteration for storing obj val, 1st one for lasso
#     
#     resp.obj = matrix(resp, n.fit, nlambda)
#     
#     # rescaled weights in std.desg
#     if(std.desg == T){
#       sd.v = apply(desg, 2, sd)*sqrt((n.fit -1)/n.fit)
#     }else{sd.v = rep(1, p.fit)}
#     sd.v.obj = matrix(sd.v, p.fit, nlambda)
#     
#     
#     # creating storage for tracking the value of objective function
#     obj = matrix(NA, iter.max, nlambda)
#     
#     # 1st iteration
#     itr = 1
#     itr.v = rep(itr, nlambda)
#     
#     # new estimates
#     beta0 = matrix(NA, n.fit, nlambda)
#     beta.est = matrix(NA, p.fit, nlambda)
#     
#     # updating beta: Lasso step
#     for (l in seq(nlambda)) {
#       
#       fit.sal = glmnet(y = resp, x = desg, lambda = lambda.fit[l], 
#                        intercept = Intercept,
#                        standardize = std.desg, standardize.response = std.resp, ...)
#       beta0[,l] = as.numeric(fit.sal$a0)
#       beta.est[,l] = as.numeric(fit.sal$beta)
#     }
#     
#     
#     # evaluating the joint objective function
#     if(gamma.fit==1){
#       
#       obj[itr,] = colSums((resp.obj - beta0 - as.matrix(desg%*%beta.est))^2)/(2*n.fit) + 
#         lambda.fit*colSums(sd.v.obj*abs(beta.est))
#       
#     }else{
#       
#       obj[itr,] = colSums((resp.obj - beta0 - as.matrix(desg%*%beta.est))^2)/(2*n.fit) + 
#         lambda.fit*colSums(sd.v.obj*abs(beta.est)) -
#         lambda.fit*(colSums(sd.v.obj)/(1- (1/gamma.fit)))
#     }
#     converged.v = rep(F, nlambda)
#     converged = F
#     
#     # running while loop until convergence
#     while(!converged) {
#       
#       # starting new iteration
#       itr = itr +1
#       itr.v[!converged.v] = itr
#       
#       # sending current estimate to old
#       beta.old = beta.est
#       
#       # updating weights for SA-Lasso for group structure
#       if(itr==2) wts = matrix(NA, p.fit, nlambda)
#       
#       wts[,!converged.v] = (abs(as.matrix(beta.old[,!converged.v])))^(-gamma.fit)
#       
#       # substituting inf weights by a large number for evaluating the joint objective function
#       wts.obj = wts
#       wts.obj[wts==Inf] = 1e+100
#       
#       # updating beta: A-Lasso step
#       for(l in which(!converged.v)){
#         
#         if(sum(wts[,l]!=Inf)==0){
#           
#           # new estimates
#           beta0[,l] = mean(resp) # intercept
#           beta.est[,l] = numeric(p.fit) # reg coef est
#           
#         }else{
#           
#           wts.l = wts[,l]
#           fit.sal = glmnet(y = resp, x = desg, lambda = lambda.fit[l]*mean(wts.l[wts.l!=Inf]),
#                            intercept = Intercept,
#                            standardize = std.desg, standardize.response = std.resp,
#                            penalty.factor = wts.l/mean(wts.l[wts.l!=Inf]), ...)
#           
#           # new estimates
#           beta0[,l] = as.numeric(fit.sal$a0) # intercept
#           beta.est[,l] = as.numeric(fit.sal$beta) # reg coef
#         }
#       }
#       
#       obj[itr,converged.v] = obj[itr-1,converged.v]
#       if(gamma.fit==1){
#         
#         obj[itr,!converged.v] = colSums((as.matrix(resp.obj[,!converged.v]) - as.matrix(beta0[,!converged.v]) - 
#                                            as.matrix(desg%*%beta.est[,!converged.v]))^2)/(2*n.fit) + 
#           lambda.fit[!converged.v]*
#           colSums(as.matrix(wts.obj[,!converged.v])*as.matrix(sd.v.obj[,!converged.v])*abs(as.matrix(beta.est[,!converged.v]))) - 
#           lambda.fit[!converged.v]*colSums(as.matrix(sd.v.obj[,!converged.v])*log(as.matrix(wts.obj[,!converged.v])))
#         
#       }else{
#         
#         obj[itr,!converged.v] = colSums((as.matrix(resp.obj[,!converged.v]) - as.matrix(beta0[,!converged.v]) - 
#                                            as.matrix(desg%*%beta.est[,!converged.v]))^2)/(2*n.fit) + 
#           lambda.fit[!converged.v]*
#           colSums(as.matrix(wts.obj[,!converged.v])*as.matrix(sd.v.obj[,!converged.v])*abs(as.matrix(beta.est[,!converged.v]))) - 
#           lambda.fit[!converged.v]*colSums(as.matrix(sd.v.obj[,!converged.v])*((as.matrix(wts.obj[,!converged.v])^(1- (1/gamma.fit)))/(1- (1/gamma.fit))))
#       }
#       
#       converged.v = abs(obj[itr-1,]-obj[itr,])<=eps
#       converged = (itr==iter.max)||(sum(!converged.v)==0)
#     }
#     
#     if(nlambda>1){
#       
#       return(list('beta0' = beta0[1,], 'beta' = beta.est,
#                   'obj' = as.matrix(obj[1:itr,]), 'n.iter' = itr.v -1,
#                   'lambda' = lambda.fit,
#                   
#                   # tracking if objective func decrease at each iteration
#                   'sign.change' = sign(diff(as.matrix(obj[1:itr,]))),
#                   
#                   # how many steps had an increase
#                   'nIncrease' = apply(sign(diff(as.matrix(obj[1:itr,]))), 2, FUN = function(v){sum(v>0)})
#       ))
#       
#     }else{
#       
#       return(list('beta0' = beta0[1,], 'beta' = as.numeric(beta.est),
#                   'obj' = obj[1:itr,], 'n.iter' = itr.v -1,
#                   
#                   # tracking if objective func decrease at each iteration
#                   'sign.change' = sign(diff(obj[1:itr,])),
#                   
#                   # how many steps had an increase
#                   'nIncrease' = sum(sign(diff(obj[1:itr,]))>0)
#       ))
#     }
#     
#   }else{
#     
#     if(names(struct)=='grp'){
#       
#       n.fit = length(resp)  # n
#       p.fit = ncol(desg)  # p
#       nlambda = length(lambda.fit)
#       iter.max = iter.max +1  # max iteration for storing obj val, 1st one for lasso
#       
#       resp.obj = matrix(resp, n.fit, nlambda)
#       
#       # rescaled weights in std.desg
#       if(std.desg == T){
#         sd.v = apply(desg, 2, sd)*sqrt((n.fit -1)/n.fit)
#       }else{sd.v = rep(1, p.fit)}
#       sd.v.obj = matrix(sd.v, p.fit, nlambda)
#       
#       
#       # group structure information
#       grp.struct = struct$grp
#       indx.NA = which(is.na(grp.struct))
#       grp.id = unique(grp.struct[!is.na(grp.struct)])
#       D = length(grp.id)
#       
#       grp.id.indx = list()
#       for (d in 1:D) grp.id.indx = c(grp.id.indx, 
#                                      list(which(grp.struct==grp.id[d])))
#       
#       # creating storage for tracking the value of objective function
#       obj = matrix(NA, iter.max, nlambda)
#       
#       # 1st iteration
#       itr = 1
#       itr.v = rep(itr, nlambda)
#       
#       # new estimates
#       beta0 = matrix(NA, n.fit, nlambda)
#       beta.est = matrix(NA, p.fit, nlambda)
#       
#       # updating beta: Lasso step
#       for (l in seq(nlambda)) {
#         
#         fit.sal = glmnet(y = resp, x = desg, lambda = lambda.fit[l], 
#                          intercept = Intercept,
#                          standardize = std.desg, standardize.response = std.resp, ...)
#         beta0[,l] = as.numeric(fit.sal$a0)
#         beta.est[,l] = as.numeric(fit.sal$beta)
#       }
#       
#       
#       # evaluating the joint objective function
#       if(gamma.fit==1){
#         
#         obj[itr,] = colSums((resp.obj - beta0 - as.matrix(desg%*%beta.est))^2)/(2*n.fit) + 
#           lambda.fit*colSums(sd.v.obj*abs(beta.est))
#         
#       }else{
#         
#         obj[itr,] = colSums((resp.obj - beta0 - as.matrix(desg%*%beta.est))^2)/(2*n.fit) + 
#           lambda.fit*colSums(sd.v.obj*abs(beta.est)) -
#           lambda.fit*(colSums(sd.v.obj)/(1- (1/gamma.fit)))
#       }
#       converged.v = rep(F, nlambda)
#       converged = F #itr==iter.max
#       
#       # running while loop until convergence
#       while(!converged) {
#         
#         # starting new iteration
#         itr = itr +1
#         itr.v[!converged.v] = itr
#         
#         # sending current estimate to old
#         beta.old = beta.est
#         
#         # updating weights for SA-Lasso for group structure
#         if(itr==2) wts = matrix(NA, p.fit, nlambda)
#         
#         for (d in 1:D){
#           
#           new.wts = (colSums(as.matrix(sd.v.obj[grp.id.indx[[d]],!converged.v])*
#                                abs(as.matrix(beta.old[grp.id.indx[[d]],!converged.v])))/sum(sd.v[grp.id.indx[[d]]]))^(-gamma.fit)
#           wts[grp.id.indx[[d]],!converged.v] = matrix(new.wts, length(grp.id.indx[[d]]),
#                                                       sum(!converged.v), byrow = T)
#         }
#         
#         wts[indx.NA,!converged.v] = (abs(as.matrix(beta.old[indx.NA,!converged.v])))^(-gamma.fit)
#         
#         # substituting inf weights by a large number for evaluating the joint objective function
#         wts.obj = wts
#         wts.obj[wts==Inf] = 1e+100
#         
#         # updating beta: A-Lasso step
#         for(l in which(!converged.v)){
#           
#           if(sum(wts[,l]!=Inf)==0){
#             
#             # new estimates
#             beta0[,l] = mean(resp) # intercept
#             beta.est[,l] = numeric(p.fit) # reg coef est
#             
#           }else{
#             
#             wts.l = wts[,l]
#             fit.sal = glmnet(y = resp, x = desg, lambda = lambda.fit[l]*mean(wts.l[wts.l!=Inf]),
#                              intercept = Intercept,
#                              standardize = std.desg, standardize.response = std.resp,
#                              penalty.factor = wts.l/mean(wts.l[wts.l!=Inf]), ...)
#             
#             # new estimates
#             beta0[,l] = as.numeric(fit.sal$a0) # intercept
#             beta.est[,l] = as.numeric(fit.sal$beta) # reg coef
#           }
#         }
#         
#         obj[itr,converged.v] = obj[itr-1,converged.v]
#         if(gamma.fit==1){
#           
#           obj[itr,!converged.v] = colSums((as.matrix(resp.obj[,!converged.v]) - as.matrix(beta0[,!converged.v]) - 
#                                              as.matrix(desg%*%beta.est[,!converged.v]))^2)/(2*n.fit) + 
#             lambda.fit[!converged.v]*
#             colSums(as.matrix(wts.obj[,!converged.v])*as.matrix(sd.v.obj[,!converged.v])*abs(as.matrix(beta.est[,!converged.v]))) - 
#             lambda.fit[!converged.v]*colSums(as.matrix(sd.v.obj[,!converged.v])*log(as.matrix(wts.obj[,!converged.v])))
#           
#         }else{
#           
#           obj[itr,!converged.v] = colSums((as.matrix(resp.obj[,!converged.v]) - as.matrix(beta0[,!converged.v]) - 
#                                              as.matrix(desg%*%beta.est[,!converged.v]))^2)/(2*n.fit) + 
#             lambda.fit[!converged.v]*
#             colSums(as.matrix(wts.obj[,!converged.v])*as.matrix(sd.v.obj[,!converged.v])*abs(as.matrix(beta.est[,!converged.v]))) - 
#             lambda.fit[!converged.v]*colSums(as.matrix(sd.v.obj[,!converged.v])*((as.matrix(wts.obj[,!converged.v])^(1- (1/gamma.fit)))/(1- (1/gamma.fit))))
#         }
#         
#         converged.v = abs(obj[itr-1,]-obj[itr,])<=eps
#         converged = (itr==iter.max)||(sum(!converged.v)==0)
#       }
#       
#       if(nlambda>1){
#         
#         return(list('beta0' = beta0[1,], 'beta' = beta.est,
#                     'obj' = as.matrix(obj[1:itr,]), 'n.iter' = itr.v -1,
#                     'lambda' = lambda.fit,
#                     
#                     # tracking if objective func decrease at each iteration
#                     'sign.change' = sign(diff(as.matrix(obj[1:itr,]))),
#                     
#                     # how many steps had an increase
#                     'nIncrease' = apply(sign(diff(as.matrix(obj[1:itr,]))), 2, FUN = function(v){sum(v>0)})
#         ))
#         
#       }else{
#         
#         return(list('beta0' = beta0[1,], 'beta' = as.numeric(beta.est),
#                     'obj' = obj[1:itr,], 'n.iter' = itr.v -1,
#                     
#                     # tracking if objective func decrease at each iteration
#                     'sign.change' = sign(diff(obj[1:itr,])),
#                     
#                     # how many steps had an increase
#                     'nIncrease' = sum(sign(diff(obj[1:itr,]))>0)
#         ))
#       }
#       
#     }else if(names(struct)=='cov'){
#       
#       n.fit = length(resp)  # n
#       p.fit = ncol(desg)  # p
#       nlambda = length(lambda.fit)
#       iter.max = iter.max +1  # max iteration for storing obj val, 1st one for lasso
#       
#       resp.obj = matrix(resp, n.fit, nlambda)
#       
#       # rescaled weights in std.desg
#       if(std.desg == T){
#         sd.v = apply(desg, 2, sd)*sqrt((n.fit -1)/n.fit)
#       }else{sd.v = rep(1, p.fit)}
#       sd.v.obj = matrix(sd.v, p.fit, nlambda)
#       
#       
#       # covariate structure information
#       u.fit = struct$cov
#       indx.NA = which(is.na(u.fit))
#       indx.noNA = which(!is.na(u.fit))
#       p.new.fit = length(indx.noNA)
#       u.std = u.fit
#       a.center = sum(sd.v[indx.noNA]*u.fit[indx.noNA])/sum(sd.v[indx.noNA])
#       b.scale = sqrt(mean(sd.v[indx.noNA]*((u.fit[indx.noNA] - a.center)^2)))
#       u.std[indx.noNA] = (u.fit[indx.noNA] - a.center)/b.scale
#       
#       
#       # creating storage for tracking the value of objective function
#       obj = matrix(NA, iter.max, nlambda)
#       
#       # 1st iteration
#       itr = 1
#       itr.v = rep(itr, nlambda)
#       
#       # new estimates
#       beta0 = matrix(NA, n.fit, nlambda)
#       beta.est = matrix(NA, p.fit, nlambda)
#       
#       # updating beta: Lasso step
#       for (l in seq(nlambda)) {
#         
#         fit.sal = glmnet(y = resp, x = desg, lambda = lambda.fit[l], 
#                          intercept = Intercept,
#                          standardize = std.desg, standardize.response = std.resp, ...)
#         beta0[,l] = as.numeric(fit.sal$a0)
#         beta.est[,l] = as.numeric(fit.sal$beta)
#       }
#       
#       
#       # evaluating the joint objective function
#       if(gamma.fit==1){
#         
#         obj[itr,] = colSums((resp.obj - beta0 - as.matrix(desg%*%beta.est))^2)/(2*n.fit) + 
#           lambda.fit*colSums(sd.v.obj*abs(beta.est))
#         
#       }else{
#         
#         obj[itr,] = colSums((resp.obj - beta0 - as.matrix(desg%*%beta.est))^2)/(2*n.fit) + 
#           lambda.fit*colSums(sd.v.obj*abs(beta.est)) -
#           lambda.fit*(colSums(sd.v.obj)/(1- (1/gamma.fit)))
#       }
#       
#       # print(obj[1,])  ###
#       converged.v = rep(F, nlambda)
#       converged = F
#       
#       
#       # running while loop until convergence
#       while(!converged) {
#         
#         # starting new iteration
#         itr = itr +1
#         itr.v[!converged.v] = itr
#         
#         # sending current estimate to old
#         beta.old = beta.est
#         
#         # updating weights for SA-Lasso for group structure
#         if(itr==2) wts = matrix(NA, p.fit, nlambda)
#         
#         # print(apply(beta.old, 2, FUN = function(v){sum(v!=0)})) ###
#         
#         wts[indx.noNA,!converged.v] = apply(X = matrix(which(!converged.v), ncol = 1), 1, 
#                                             FUN = function(l){
#                                               
#                                               Sbeta = which(beta.old[,l]!=0)
#                                               if(gamma.fit==1){
#                                                 
#                                                 opt.out = nlminb(start = c(0,0),
#                                                                  lower = c(-30,-30), upper = c(30,30),
#                                                                  objective = function(x){
#                                                                    
#                                                                    sum(exp(x[1] + x[2]*u.std[intersect(Sbeta, indx.noNA)])*
#                                                                          sd.v[intersect(Sbeta, indx.noNA)]*abs(beta.old[intersect(Sbeta, indx.noNA),l])) -
#                                                                      x[1]*sum(sd.v[indx.noNA]) #- x[2]*sum(sd.v[indx.noNA]*u.std[indx.noNA])
#                                                                    
#                                                                  }, control = list(eval.max = 5000, iter.max = 5000))
#                                                 tau.opt = opt.out$par
#                                                 
#                                               }else{
#                                                 
#                                                 opt.out = nlminb(start = c(0,0),
#                                                                  lower = c(-30,-30), upper = c(30,30),
#                                                                  objective = function(x){
#                                                                    
#                                                                    sum(exp(x[1] + x[2]*u.std[intersect(Sbeta, indx.noNA)])*
#                                                                          sd.v[intersect(Sbeta, indx.noNA)]*abs(beta.old[intersect(Sbeta, indx.noNA),l])) -
#                                                                      sum(sd.v[indx.noNA]*exp((1- (1/gamma.fit))*(x[1] + x[2]*u.std[indx.noNA])))/
#                                                                      (1- (1/gamma.fit))
#                                                                    
#                                                                  }, control = list(eval.max = 5000, iter.max = 5000))
#                                                 tau.opt = opt.out$par
#                                               }
#                                               
#                                               exp(tau.opt[1] + tau.opt[2]*u.std[indx.noNA])
#                                             })
#         wts[indx.NA,!converged.v] = (abs(as.matrix(beta.old[indx.NA,!converged.v])))^(-gamma.fit)
#         
#         
#         # substituting inf weights by a large number for evaluating the joint objective function
#         wts.obj = wts
#         wts.obj[wts==Inf] = 1e+100
#         
#         
#         # updating beta: A-Lasso step
#         for(l in which(!converged.v)){
#           
#           if(sum(wts[,l]!=Inf)==0){
#             
#             # new estimates
#             beta0[,l] = mean(resp) # intercept
#             beta.est[,l] = numeric(p.fit) # reg coef est
#             
#           }else{
#             
#             wts.l = wts[,l]
#             fit.sal = glmnet(y = resp, x = desg, lambda = lambda.fit[l]*mean(wts.l[wts.l!=Inf]),
#                              intercept = Intercept,
#                              standardize = std.desg, standardize.response = std.resp,
#                              penalty.factor = wts.l/mean(wts.l[wts.l!=Inf]), ...)
#             
#             # new estimates
#             beta0[,l] = as.numeric(fit.sal$a0) # intercept
#             beta.est[,l] = as.numeric(fit.sal$beta) # reg coef
#           }
#         }
#         
#         obj[itr,converged.v] = obj[itr-1,converged.v]
#         if(gamma.fit==1){
#           
#           obj[itr,!converged.v] = colSums((as.matrix(resp.obj[,!converged.v]) - as.matrix(beta0[,!converged.v]) - 
#                                              as.matrix(desg%*%beta.est[,!converged.v]))^2)/(2*n.fit) + 
#             lambda.fit[!converged.v]*
#             colSums(as.matrix(wts.obj[,!converged.v])*as.matrix(sd.v.obj[,!converged.v])*abs(as.matrix(beta.est[,!converged.v]))) - 
#             lambda.fit[!converged.v]*colSums(as.matrix(sd.v.obj[,!converged.v])*log(as.matrix(wts.obj[,!converged.v])))
#           
#         }else{
#           
#           obj[itr,!converged.v] = colSums((as.matrix(resp.obj[,!converged.v]) - as.matrix(beta0[,!converged.v]) - 
#                                              as.matrix(desg%*%beta.est[,!converged.v]))^2)/(2*n.fit) + 
#             lambda.fit[!converged.v]*
#             colSums(as.matrix(wts.obj[,!converged.v])*as.matrix(sd.v.obj[,!converged.v])*abs(as.matrix(beta.est[,!converged.v]))) - 
#             lambda.fit[!converged.v]*colSums(as.matrix(sd.v.obj[,!converged.v])*((as.matrix(wts.obj[,!converged.v])^(1- (1/gamma.fit)))/(1- (1/gamma.fit))))
#         }
#         
#         converged.v = abs(obj[itr-1,]-obj[itr,])<=eps
#         converged = (itr==iter.max)||(sum(!converged.v)==0)
#         
#         # print(itr)
#         # print(obj[itr,])
#       }
#       
#       if(nlambda>1){
#         
#         return(list('beta0' = beta0[1,], 'beta' = beta.est,
#                     'obj' = as.matrix(obj[1:itr,]), 'n.iter' = itr.v -1,
#                     'lambda' = lambda.fit,
#                     
#                     # tracking if objective func decrease at each iteration
#                     'sign.change' = sign(diff(as.matrix(obj[1:itr,]))),
#                     
#                     # how many steps had an increase
#                     'nIncrease' = apply(sign(diff(as.matrix(obj[1:itr,]))), 2, FUN = function(v){sum(v>0)})
#         ))
#         
#       }else{
#         
#         return(list('beta0' = beta0[1,], 'beta' = as.numeric(beta.est),
#                     'obj' = obj[1:itr,], 'n.iter' = itr.v -1,
#                     
#                     # tracking if objective func decrease at each iteration
#                     'sign.change' = sign(diff(obj[1:itr,])),
#                     
#                     # how many steps had an increase
#                     'nIncrease' = sum(sign(diff(obj[1:itr,]))>0)
#         ))
#       }
#     }
#   }
# }
# 
# 
# 
# 
# cvSA.Lasso = function(resp, desg, struct, n.fold = 10,
#                       Intercept = F, std.desg = T, std.resp = F,
#                       opt.lambda.type = 'min',
#                       gam.seq = NULL, gam.low = .01, gam.up = 1, n.gam = 50,
#                       lambda.seq = NULL, lambda.min = 1e-8, lambda.max = NULL,
#                       n.lambda = 100, iter.max = 100, eps = 1e-4, #seed = 1,
#                       nCore = 3, ...){
#   
#   n.data = ifelse(length(resp)==nrow(desg), length(resp), 
#                   print('ERROR! length of `resp` and number of rows in `desg` should be same!'))
#   
#   # sequence of gamma
#   if(missing(gam.seq)) gam.seq = seq(gam.low, gam.up, length.out = n.gam)
#   
#   ## decreasing sequence of lambda
#   
#   # determining lambda.max
#   if(missing(lambda.max)){
#     
#     if(std.desg){
#       
#       A.l = scale(desg)*sqrt(n.data/(n.data -1))
#       
#     }else{A.l = desg}
#     
#     if(std.resp){
#       
#       a.l = resp/(sqrt((n.data -1)/n.data)*sd(resp))
#       
#     }else{a.l = resp}
#     
#     lambda.max = max(abs(as.numeric(t(A.l)%*%a.l)))/n.data
#   }
#   
#   # if(missing(lambda.seq)) lambda.seq = 10^(seq(log10(lambda.max), log10(lambda.min), by = -.01))
#   
#   if(missing(lambda.seq)) lambda.seq = 10^(seq(log10(lambda.max), log10(lambda.min),
#                                                length.out = n.lambda))
#   
#   # if(missing(lambda.seq)) lambda.seq = exp(seq(log(lambda.max), log(lambda.min), 
#   #                                              length.out = n.lambda))
#   # n.lambda = length(lambda.seq)
#   
#   
#   if(missing(struct)==T){
#     
#     # randomly partitioning the data
#     folds = cut(seq(length(resp)), breaks = n.fold, labels=FALSE)
#     # set.seed(seed)
#     folds = permute(folds)
#     
#     if(opt.lambda.type=='min'){
#       
#       registerDoParallel(cores = nCore)
#       cv.out = foreach(g = seq(n.gam), .combine = 'rbind') %dopar% {
#         
#         cvm.temp = t(apply(X = matrix(seq(n.fold), ncol = 1), 1, 
#                            FUN = function(v){
#                              
#                              test.indx = which(folds==v)
#                              out.sal = SA.Lasso(resp = resp[-test.indx], desg = desg[-test.indx,],
#                                                 lambda.fit = lambda.seq, gamma.fit = gam.seq[g],
#                                                 Intercept = Intercept, std.desg = std.desg, std.resp = std.resp,
#                                                 iter.max = iter.max, eps = eps, ...)
#                              
#                              # Mean Square Prediction Error
#                              colMeans((matrix(resp[test.indx], length(test.indx), n.lambda) -
#                                          matrix(out.sal$beta0, length(test.indx), n.lambda, byrow = T) -
#                                          as.matrix(desg[test.indx,]%*%out.sal$beta))^2)
#                            }))
#         
#         cvm.seq = colMeans(cvm.temp)
#         
#         indx.min.cvm = min(which(cvm.seq==min(cvm.seq)))
#         
#         c(gam.seq[g], lambda.seq[indx.min.cvm], cvm.seq[indx.min.cvm])
#       }
#       
#     }else if(opt.lambda.type=='1se'){
#       
#       registerDoParallel(cores = nCore)
#       cv.out = foreach(g = seq(n.gam), .combine = 'rbind') %dopar% {
#         
#         cvm.temp = t(apply(X = matrix(seq(n.fold), ncol = 1), 1, 
#                            FUN = function(v){
#                              
#                              test.indx = which(folds==v)
#                              out.sal = SA.Lasso(resp = resp[-test.indx], desg = desg[-test.indx,],
#                                                 lambda.fit = lambda.seq, gamma.fit = gam.seq[g],
#                                                 Intercept = Intercept, std.desg = std.desg, std.resp = std.resp,
#                                                 iter.max = iter.max, eps = eps, ...)
#                              
#                              # Mean Square Prediction Error
#                              colMeans((matrix(resp[test.indx], length(test.indx), n.lambda) -
#                                          matrix(out.sal$beta0, length(test.indx), n.lambda, byrow = T) -
#                                          as.matrix(desg[test.indx,]%*%out.sal$beta))^2)
#                            }))
#         
#         cvm.seq = colMeans(cvm.temp) + apply(cvm.temp, 2, sd)
#         
#         indx.min.cvm = min(which(cvm.seq==min(cvm.seq)))
#         
#         c(gam.seq[g], lambda.seq[indx.min.cvm], cvm.seq[indx.min.cvm])
#       }
#     }
#     
#     cv.out = as.data.frame(cv.out)
#     colnames(cv.out) = c('g', 'l', 'cvm')
#     
#     indx.opt = max(which(cv.out$cvm==min(cv.out$cvm)))
#     g.opt = cv.out$g[indx.opt]
#     l.opt = cv.out$l[indx.opt]
#     
#     sal.out = SA.Lasso(resp = resp, desg = desg, lambda.fit = l.opt, gamma.fit = g.opt,
#                        Intercept = Intercept, std.desg = std.desg, std.resp = std.resp,
#                        iter.max = iter.max, eps = eps, ...)
#     
#     return(list('beta0' = sal.out$beta0, 'beta' = sal.out$beta,
#                 'lambda.opt' = l.opt, 'gamma.opt' = g.opt,
#                 'cvm' = min(cv.out$cvm), 'cv.df' = cv.out))
#     
#   }else{
#     
#     # randomly partitioning the data
#     folds = cut(seq(length(resp)), breaks = n.fold, labels=FALSE)
#     # set.seed(seed)
#     folds = permute(folds)
#     
#     if(opt.lambda.type=='min'){
#       
#       registerDoParallel(cores = nCore)
#       cv.out = foreach(g = seq(n.gam), .combine = 'rbind') %dopar% {
#         
#         cvm.temp = t(apply(X = matrix(seq(n.fold), ncol = 1), 1, 
#                            FUN = function(v){
#                              
#                              test.indx = which(folds==v)
#                              out.sal = SA.Lasso(resp = resp[-test.indx], desg = desg[-test.indx,],
#                                                 lambda.fit = lambda.seq, gamma.fit = gam.seq[g],
#                                                 Intercept = Intercept, std.desg = std.desg, std.resp = std.resp,
#                                                 struct = struct,
#                                                 iter.max = iter.max, eps = eps, ...)
#                              
#                              # Mean Square Prediction Error
#                              colMeans((matrix(resp[test.indx], length(test.indx), n.lambda) -
#                                          matrix(out.sal$beta0, length(test.indx), n.lambda, byrow = T) -
#                                          as.matrix(desg[test.indx,]%*%out.sal$beta))^2)
#                            }))
#         
#         cvm.seq = colMeans(cvm.temp)
#         
#         indx.min.cvm = min(which(cvm.seq==min(cvm.seq)))
#         
#         c(gam.seq[g], lambda.seq[indx.min.cvm], cvm.seq[indx.min.cvm])
#       }
#       
#     }else if(opt.lambda.type=='1se'){
#       
#       registerDoParallel(cores = nCore)
#       cv.out = foreach(g = seq(n.gam), .combine = 'rbind') %dopar% {
#         
#         cvm.temp = t(apply(X = matrix(seq(n.fold), ncol = 1), 1, 
#                            FUN = function(v){
#                              
#                              test.indx = which(folds==v)
#                              out.sal = SA.Lasso(resp = resp[-test.indx], desg = desg[-test.indx,],
#                                                 lambda.fit = lambda.seq, gamma.fit = gam.seq[g],
#                                                 Intercept = Intercept, std.desg = std.desg, std.resp = std.resp,
#                                                 struct = struct,
#                                                 iter.max = iter.max, eps = eps, ...)
#                              
#                              # Mean Square Prediction Error
#                              colMeans((matrix(resp[test.indx], length(test.indx), n.lambda) -
#                                          matrix(out.sal$beta0, length(test.indx), n.lambda, byrow = T) -
#                                          as.matrix(desg[test.indx,]%*%out.sal$beta))^2)
#                            }))
#         
#         cvm.seq = colMeans(cvm.temp) + apply(cvm.temp, 2, sd)
#         
#         indx.min.cvm = min(which(cvm.seq==min(cvm.seq)))
#         
#         c(gam.seq[g], lambda.seq[indx.min.cvm], cvm.seq[indx.min.cvm])
#       }
#     }
#     
#     cv.out = as.data.frame(cv.out)
#     colnames(cv.out) = c('g', 'l', 'cvm')
#     
#     indx.opt = max(which(cv.out$cvm==min(cv.out$cvm)))
#     g.opt = cv.out$g[indx.opt]
#     l.opt = cv.out$l[indx.opt]
#     
#     sal.out = SA.Lasso(resp = resp, desg = desg, lambda.fit = l.opt, gamma.fit = g.opt,
#                        Intercept = Intercept, std.desg = std.desg, std.resp = std.resp,
#                        struct = struct,
#                        iter.max = iter.max, eps = eps, ...)
#     
#     return(list('beta0' = sal.out$beta0, 'beta' = sal.out$beta,
#                 'lambda.opt' = l.opt, 'gamma.opt' = g.opt,
#                 'cvm' = min(cv.out$cvm), 'cv.df' = cv.out))
#   }
# }



# asymptotic objective function for covariate structure
asymp.obj.func = function( b, alpha.l, tau.l, B.vec, Z.vec, U.vec, gam = 1){
  
  if(gam==1){
    
    # parameterization
    eta.vec = mapply( eta, x = B.vec + tau.l*Z.vec,
                      MoreArgs = list(theta = (alpha.l*tau.l)) )
    w.v = exp(b*U.vec)
    out = mean( w.v*abs(eta.vec) - log(w.v))
    
  }else if((gam>0) && (gam<1)){
    
    # parameterization
    eta.vec = mapply( eta, x = B.vec + tau.l*Z.vec,
                      MoreArgs = list(theta = (alpha.l*tau.l)) )
    w.v = exp(b*U.vec)
    out = mean( w.v*abs(eta.vec) - (w.v^(1- gam)) )
  }
  
  return(out)
}

asymp.f = function(tau, alpha.l, tau.l, B.vec, Z.vec, U.vec, gam=1) {
  
  eta.vec = mapply( eta, x = B.vec + tau.l*Z.vec,
                    MoreArgs = list(theta = (alpha.l*tau.l)) )
  
  mean(U.vec*exp(tau*U.vec)*abs(eta.vec) - U.vec)*as.numeric(gam==1) +
    mean(U.vec*exp(tau*U.vec)*(abs(eta.vec) - exp(-(tau*U.vec)/gam)))*as.numeric(gam!=1)
}

