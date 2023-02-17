
source('functions.R')

SAEnet = function(y, X, structure.info, nIteration = 5,
                  std.y = F, std.X = F,
                  lambda2.seq = exp(seq(0, -40, length.out = 50)),
                  gamma.seq = seq(.1, 1, by = .1),
                  nfolds = 10, foldid, nCore, verbose = T){
  
  # number of observations
  n.X = nrow(X)
  n.y = length(y)
  if(n.X!=n.y){
    
    return('Number of observations from y and X should match. The length of y and
           the number of rows in X do not match!')
    
  }else{n = n.y}
  
  is.foldid = !missing(foldid)
  
  # default no. of cores if parallel
  if(missing(nCore)) nCore = max(detectCores() - 1, 1)
  
  # number of predictor variables
  p = ncol(X)
  
  ## transforming to a "0 intercept" regression
  # standardizing y
  y.scale = ifelse(std.y, sd(y), 1)
  y.std = y/y.scale
  X.std = X
  
  if(missing(structure.info)){
    
    #### no structural info ####
    
    # storage for outputs
    betahat_std.out = w.out = matrix(nrow = p, ncol = 1 + nIteration)
    intercept_std.out = cvm.out = lambda1.out = lambda2.out = 
      gamma.out = rep(NA, 1 + nIteration)
    
    # Implementing Adaptive Enet (AEnet) with adaptive L1 and non-adaptive L2 penalties
    # the problem is rewritten by merging L2 norm with the loss function as per Zou and Hastie (2005)
    for(itr.SAEnet in 0:nIteration){
      
      if(itr.SAEnet==0){
        
        #### Iteration 0 ####
        doParallel::registerDoParallel(cores = nCore)
        cv.SAEnet = foreach::foreach(lambda2 = lambda2.seq, .combine = 'rbind', .multicombine = T, .export = 'foldid.cv') %dopar% {
          
          if(is.foldid){
            
            cv.SAEnet_lambda2 = 
              gcdnet::cv.gcdnet(x = X.std, y = y.std, 
                                method = 'ls', pred.loss = 'loss',
                                standardize = std.X,
                                foldid = foldid, lambda2 = lambda2)
            
          }else{
            
            cv.SAEnet_lambda2 = 
              gcdnet::cv.gcdnet(x = X.std, y = y.std, 
                                method = 'ls', pred.loss = 'loss',
                                standardize = std.X,
                                nfolds = nfolds, lambda2 = lambda2)
          }
          
          c(min(cv.SAEnet_lambda2$cvm), cv.SAEnet_lambda2$lambda.min, lambda2)
        }
        
        opt.tuning.id = max(which(cv.SAEnet[,1]==min(cv.SAEnet[,1])))
        fit.SAEnet = gcdnet::gcdnet(x = X.std, y = y.std, method = 'ls', 
                                    standardize = std.X,
                                    lambda = cv.SAEnet[opt.tuning.id,2],
                                    lambda2 = cv.SAEnet[opt.tuning.id,3])
        
        intercept_std.out[itr.SAEnet+1] = as.numeric(fit.SAEnet$b0)
        betahat_std.out[,itr.SAEnet+1] = as.numeric(fit.SAEnet$beta)
        cvm.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,1]
        lambda1.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,2]
        lambda2.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,3]
        
      }else{
        
        #### Iteration 1 to nIteration ####
        
        # base adaptive penalties
        wts.init.SAEnet = 1/abs(betahat_std.out[,itr.SAEnet])
        
        if(any(sum(betahat_std.out[,itr.SAEnet]!=0)==0, sum(wts.init.SAEnet!=Inf)==0)){
          
          betahat_std.out[,itr.SAEnet+1] = numeric(p)
          
        }else{
          
          doParallel::registerDoParallel(cores = nCore)
          cv.SAEnet = foreach::foreach(lambda2 = lambda2.seq, .combine = 'rbind', .multicombine = T) %dopar% {
            
            if(length(gamma.seq)==1){
              
              if(is.foldid){
                
                cv.SAEnet_lambda2.g = 
                  gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                    method = 'ls', pred.loss = 'loss',
                                    standardize = std.X, pf = pmin(wts.init.SAEnet^gamma.seq, 1e+100),
                                    foldid = foldid, lambda2 = lambda2)
                
              }else{
                
                cv.SAEnet_lambda2.g = 
                  gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                    method = 'ls', pred.loss = 'loss',
                                    standardize = std.X, pf = pmin(wts.init.SAEnet^gamma.seq, 1e+100),
                                    nfolds = nfolds, lambda2 = lambda2)
              }
              
              cv.SAEnet_lambda2 = matrix(c(min(cv.SAEnet_lambda2.g$cvm), 
                                           cv.SAEnet_lambda2.g$lambda.min, gamma.seq),
                                         nrow = 1, ncol = 3, byrow = T)
              opt.tuning.id_lambda2 = 1
              
            }else{
              
              cv.SAEnet_lambda2 = foreach(g = gamma.seq, .combine = 'rbind', .multicombine = T) %do% {
                
                wts.SAEnet_g = pmin(wts.init.SAEnet^g, 1e+100)
                
                if(is.foldid){
                  
                  cv.SAEnet_lambda2.g = 
                    gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                      method = 'ls', pred.loss = 'loss',
                                      standardize = std.X, pf = pmin(wts.init.SAEnet^gamma.seq, 1e+100),
                                      foldid = foldid, lambda2 = lambda2)
                  
                }else{
                  
                  cv.SAEnet_lambda2.g = 
                    gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                      method = 'ls', pred.loss = 'loss',
                                      standardize = std.X, pf = pmin(wts.init.SAEnet^gamma.seq, 1e+100),
                                      nfolds = nfolds, lambda2 = lambda2)
                }
                
                c(min(cv.SAEnet_lambda2.g$cvm), cv.SAEnet_lambda2.g$lambda.min, g)
              }
              
              opt.tuning.id_lambda2 = max(which(cv.SAEnet_lambda2[,1]==min(cv.SAEnet_lambda2[,1])))
            }
            
            c(cv.SAEnet_lambda2[opt.tuning.id_lambda2,], lambda2)
          }
          
          opt.tuning.id = max(which(cv.SAEnet[,1]==min(cv.SAEnet[,1])))
          opt.wts.SAEnet = pmin(wts.init.SAEnet^cv.SAEnet[opt.tuning.id,3], 1e+100)
          fit.SAEnet = gcdnet::gcdnet(x = X.std, y = y.std, method = 'ls', 
                                      standardize = std.X, pf = opt.wts.SAEnet,
                                      lambda = cv.SAEnet[opt.tuning.id,2],
                                      lambda2 = cv.SAEnet[opt.tuning.id,4])
          
          intercept_std.out[itr.SAEnet+1] = as.numeric(fit.SAEnet$b0)
          betahat_std.out[,itr.SAEnet+1] = as.numeric(fit.SAEnet$beta)
          w.out[,itr.SAEnet+1] = opt.wts.SAEnet
          cvm.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,1]
          lambda1.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,2]
          gamma.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,3]
          lambda2.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,4]
        }
      }
      
      if(verbose) print(paste('Iteration', itr.SAEnet))
    }
    
  }else if(names(structure.info)=='grp'){
    
    #### group structural info ####
    
    # provided group structure information
    grp.id = unique(structure.info$grp)
    D = length(grp.id)
    
    grp.id.indx = vector(mode = 'list', length = D)
    for (d in 1:D) grp.id.indx[[d]] = which(structure.info$grp==grp.id[d])
    
    # storage for outputs
    betahat_std.out = w.out = matrix(nrow = p, ncol = 1 + nIteration)
    intercept_std.out = cvm.out = lambda1.out = lambda2.out = 
      gamma.out = rep(NA, 1 + nIteration)
    
    # Implementing Adaptive Enet (AEnet) with adaptive L1 and non-adaptive L2 penalties
    # the problem is rewritten by merging L2 norm with the loss function as per Zou and Hastie (2005)
    for(itr.SAEnet in 0:nIteration){
      
      if(itr.SAEnet==0){
        
        #### Iteration 0 ####
        doParallel::registerDoParallel(cores = nCore)
        cv.SAEnet = foreach::foreach(lambda2 = lambda2.seq, .combine = 'rbind', .multicombine = T) %dopar% {
          
          if(is.foldid){
            
            cv.SAEnet_lambda2 = 
              gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                method = 'ls', pred.loss = 'loss',
                                standardize = std.X,
                                foldid = foldid, lambda2 = lambda2)
            
          }else{
            
            cv.SAEnet_lambda2 = 
              gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                method = 'ls', pred.loss = 'loss',
                                standardize = std.X,
                                nfolds = nfolds, lambda2 = lambda2)
          }
          
          c(min(cv.SAEnet_lambda2$cvm), cv.SAEnet_lambda2$lambda.min, lambda2)
        }
        
        opt.tuning.id = max(which(cv.SAEnet[,1]==min(cv.SAEnet[,1])))
        fit.SAEnet = gcdnet::gcdnet(x = X.std, y = y.std, method = 'ls',
                                    standardize = std.X,
                                    lambda = cv.SAEnet[opt.tuning.id,2],
                                    lambda2 = cv.SAEnet[opt.tuning.id,3])
        
        intercept_std.out[itr.SAEnet+1] = as.numeric(fit.SAEnet$b0)
        betahat_std.out[,itr.SAEnet+1] = as.numeric(fit.SAEnet$beta)
        cvm.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,1]
        lambda1.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,2]
        lambda2.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,3]
        
      }else{
        
        #### Iteration 1 to nIteration ####
        
        # base adaptive penalties
        wts.init.SAEnet = rep(NA, p)
        for (d in 1:D) wts.init.SAEnet[grp.id.indx[[d]]] = 1/mean(abs(betahat_std.out[grp.id.indx[[d]],
                                                                                      itr.SAEnet]))
        
        if(any(sum(betahat_std.out[,itr.SAEnet]!=0)==0, sum(wts.init.SAEnet!=Inf)==0)){
          
          betahat_std.out[,itr.SAEnet+1] = numeric(p)
          
        }else{
          
          doParallel::registerDoParallel(cores = nCore)
          cv.SAEnet = foreach::foreach(lambda2 = lambda2.seq, .combine = 'rbind', .multicombine = T) %dopar% {
            
            if(length(gamma.seq)==1){
              
              if(is.foldid){
                
                cv.SAEnet_lambda2.g = 
                  gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                    method = 'ls', pred.loss = 'loss',
                                    standardize = std.X, pf = pmin(wts.init.SAEnet^gamma.seq, 1e+100),
                                    foldid = foldid, lambda2 = lambda2)
                
              }else{
                
                cv.SAEnet_lambda2.g = 
                  gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                    method = 'ls', pred.loss = 'loss',
                                    standardize = std.X, pf = pmin(wts.init.SAEnet^gamma.seq, 1e+100),
                                    nfolds = nfolds, lambda2 = lambda2)
              }
              
              cv.SAEnet_lambda2 = matrix(c(min(cv.SAEnet_lambda2.g$cvm), 
                                           cv.SAEnet_lambda2.g$lambda.min, gamma.seq),
                                         nrow = 1, ncol = 3, byrow = T)
              opt.tuning.id_lambda2 = 1
              
            }else{
              
              cv.SAEnet_lambda2 = foreach(g = gamma.seq, .combine = 'rbind', .multicombine = T) %do% {
                
                wts.SAEnet_g = pmin(wts.init.SAEnet^g, 1e+100)
                
                if(is.foldid){
                  
                  cv.SAEnet_lambda2.g = 
                    gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                      method = 'ls', pred.loss = 'loss',
                                      standardize = std.X, pf = wts.SAEnet_g,
                                      foldid = foldid, lambda2 = lambda2)
                  
                }else{
                  
                  cv.SAEnet_lambda2.g = 
                    gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                      method = 'ls', pred.loss = 'loss',
                                      standardize = std.X, pf = wts.SAEnet_g,
                                      nfolds = nfolds, lambda2 = lambda2)
                }
                
                c(min(cv.SAEnet_lambda2.g$cvm), cv.SAEnet_lambda2.g$lambda.min, g)
              }
              
              opt.tuning.id_lambda2 = max(which(cv.SAEnet_lambda2[,1]==min(cv.SAEnet_lambda2[,1])))
            }
            
            c(cv.SAEnet_lambda2[opt.tuning.id_lambda2,], lambda2)
          }
          
          opt.tuning.id = max(which(cv.SAEnet[,1]==min(cv.SAEnet[,1])))
          opt.wts.SAEnet = pmin(wts.init.SAEnet^cv.SAEnet[opt.tuning.id,3], 1e+100)
          fit.SAEnet = gcdnet::gcdnet(x = X.std, y = y.std, method = 'ls', 
                                      standardize = std.X, pf = opt.wts.SAEnet,
                                      lambda = cv.SAEnet[opt.tuning.id,2],
                                      lambda2 = cv.SAEnet[opt.tuning.id,4])
          
          intercept_std.out[itr.SAEnet+1] = as.numeric(fit.SAEnet$b0)
          betahat_std.out[,itr.SAEnet+1] = as.numeric(fit.SAEnet$beta)
          w.out[,itr.SAEnet+1] = opt.wts.SAEnet
          cvm.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,1]
          lambda1.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,2]
          gamma.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,3]
          lambda2.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,4]
        }
      }
      
      if(verbose) print(paste('Iteration', itr.SAEnet))
    }
    
  }else if(names(structure.info)=='cov'){
    
    #### covariate-dependent structural info ####
    
    # provided group structure information
    u = structure.info$cov
    u.center = mean(u)
    u.scale = sd(u)
    u.std = (u - u.center)/u.scale
    
    # storage for outputs
    betahat_std.out = w.out = matrix(nrow = p, ncol = 1 + nIteration)
    intercept_std.out = cvm.out = lambda1.out = lambda2.out = 
      gamma.out = rep(NA, 1 + nIteration)
    
    # Implementing Adaptive Enet (AEnet) with adaptive L1 and non-adaptive L2 penalties
    # the problem is rewritten by merging L2 norm with the loss function as per Zou and Hastie (2005)
    for(itr.SAEnet in 0:nIteration){
      
      if(itr.SAEnet==0){
        
        #### Iteration 0 ####
        doParallel::registerDoParallel(cores = nCore)
        cv.SAEnet = foreach::foreach(lambda2 = lambda2.seq, .combine = 'rbind', .multicombine = T) %dopar% {
          
          if(is.foldid){
            
            cv.SAEnet_lambda2 = 
              gcdnet::cv.gcdnet(x = X.std, y = y.std, 
                                standardize = std.X,
                                method = 'ls', pred.loss = 'loss',
                                foldid = foldid, lambda2 = lambda2)
            
          }else{
            
            cv.SAEnet_lambda2 = 
              gcdnet::cv.gcdnet(x = X.std, y = y.std, 
                                standardize = std.X,
                                method = 'ls', pred.loss = 'loss',
                                nfolds = nfolds, lambda2 = lambda2)
          }
          
          c(min(cv.SAEnet_lambda2$cvm), cv.SAEnet_lambda2$lambda.min, lambda2)
        }
        
        opt.tuning.id = max(which(cv.SAEnet[,1]==min(cv.SAEnet[,1])))
        fit.SAEnet = gcdnet::gcdnet(x = X.std, y = y.std, method = 'ls',
                                    standardize = std.X,
                                    lambda = cv.SAEnet[opt.tuning.id,2],
                                    lambda2 = cv.SAEnet[opt.tuning.id,3])
        
        intercept_std.out[itr.SAEnet+1] = as.numeric(fit.SAEnet$b0)
        betahat_std.out[,itr.SAEnet+1] = as.numeric(fit.SAEnet$beta)
        cvm.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,1]
        lambda1.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,2]
        lambda2.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,3]
        
      }else{
        
        #### Iteration 1 to nIteration ####
        
        doParallel::registerDoParallel(cores = nCore)
        cv.SAEnet = foreach::foreach(lambda2 = lambda2.seq, .combine = 'rbind', .multicombine = T) %dopar% {
          
          if(length(gamma.seq)==1){
            
            # adaptive penalties
            if(sum(betahat_std.out[,itr.SAEnet]!=0)==0){
              
              if(gamma.seq==1){
                
                opt.out_gamma = nlminb(start = c(0,0),
                                       lower = c(-30,-30), upper = c(30,30),
                                       objective = function(x){
                                         
                                         -(p*x[1])
                                         
                                       }, control = list(eval.max = 5000, iter.max = 5000))
                tau.opt_gamma = opt.out_gamma$par
                
              }else if((gamma.seq>0)&&(gamma.seq<1)){
                
                opt.out_gamma = nlminb(start = c(0,0),
                                       lower = c(-30,-30), upper = c(30,30),
                                       objective = function(x){
                                         
                                         -(sum(exp((1-(1/gamma.seq))*(x[1] + x[2]*u.std))))/(1-(1/gamma.seq))
                                         
                                       }, control = list(eval.max = 5000, iter.max = 5000))
                tau.opt_gamma = opt.out_gamma$par
              }
              
            }else{
              
              Sbeta = which(betahat_std.out[,itr.SAEnet]!=0)
              
              if(gamma.seq==1){
                
                opt.out_gamma = nlminb(start = c(0,0),
                                       lower = c(-30,-30), upper = c(30,30),
                                       objective = function(x){
                                         
                                         sum(exp(x[1] + x[2]*u.std[Sbeta])*abs(betahat_std.out[Sbeta,itr.SAEnet])) - (p*x[1])
                                         
                                       }, control = list(eval.max = 5000, iter.max = 5000))
                tau.opt_gamma = opt.out_gamma$par
                
              }else if((gamma.seq>0)&&(gamma.seq<1)){
                
                opt.out_gamma = nlminb(start = c(0,0),
                                       lower = c(-30,-30), upper = c(30,30),
                                       objective = function(x){
                                         
                                         sum(exp(x[1] + x[2]*u.std[Sbeta])*abs(betahat_std.out[Sbeta,itr.SAEnet])) - 
                                           (sum(exp((1-(1/gamma.seq))*(x[1] + x[2]*u.std))))/(1-(1/gamma.seq))
                                         
                                       }, control = list(eval.max = 5000, iter.max = 5000))
                tau.opt_gamma = opt.out_gamma$par
              }
            }
            
            wts.SAEnet_g = pmin(exp(tau.opt_gamma[1] + tau.opt_gamma[2]*u.std), 1e+100)
            
            if(is.foldid){
              
              cv.SAEnet_lambda2.g = 
                gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                  method = 'ls', pred.loss = 'loss',
                                  standardize = std.X, pf = wts.SAEnet_g,
                                  foldid = foldid, lambda2 = lambda2)
              
            }else{
              
              cv.SAEnet_lambda2.g = 
                gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                  method = 'ls', pred.loss = 'loss',
                                  standardize = std.X, pf = wts.SAEnet_g,
                                  nfolds = nfolds, lambda2 = lambda2)
            }
            
            cv.SAEnet_lambda2 = matrix(c(min(cv.SAEnet_lambda2.g$cvm), 
                                         cv.SAEnet_lambda2.g$lambda.min, gamma.seq),
                                       nrow = 1, ncol = 3, byrow = T)
            opt.tuning.id_lambda2 = 1
            
          }else{
            
            cv.SAEnet_lambda2 = foreach(g = gamma.seq, .combine = 'rbind', .multicombine = T) %do% {
              
              # adaptive penalties
              if(sum(betahat_std.out[,itr.SAEnet]!=0)==0){
                
                if(g==1){
                  
                  opt.out_gamma = nlminb(start = c(0,0),
                                         lower = c(-30,-30), upper = c(30,30),
                                         objective = function(x){
                                           
                                           -(p*x[1])
                                           
                                         }, control = list(eval.max = 5000, iter.max = 5000))
                  tau.opt_gamma = opt.out_gamma$par
                  
                }else if((g>0)&&(g<1)){
                  
                  opt.out_gamma = nlminb(start = c(0,0),
                                         lower = c(-30,-30), upper = c(30,30),
                                         objective = function(x){
                                           
                                           -(sum(exp((1-(1/g))*(x[1] + x[2]*u.std))))/(1-(1/g))
                                           
                                         }, control = list(eval.max = 5000, iter.max = 5000))
                  tau.opt_gamma = opt.out_gamma$par
                }
                
              }else{
                
                Sbeta = which(betahat_std.out[,itr.SAEnet]!=0)
                
                if(g==1){
                  
                  opt.out_gamma = nlminb(start = c(0,0),
                                         lower = c(-30,-30), upper = c(30,30),
                                         objective = function(x){
                                           
                                           sum(exp(x[1] + x[2]*u.std[Sbeta])*abs(betahat_std.out[Sbeta,itr.SAEnet])) - (p*x[1])
                                           
                                         }, control = list(eval.max = 5000, iter.max = 5000))
                  tau.opt_gamma = opt.out_gamma$par
                  
                }else if((g>0)&&(g<1)){
                  
                  opt.out_gamma = nlminb(start = c(0,0),
                                         lower = c(-30,-30), upper = c(30,30),
                                         objective = function(x){
                                           
                                           sum(exp(x[1] + x[2]*u.std[Sbeta])*abs(betahat_std.out[Sbeta,itr.SAEnet])) - 
                                             (sum(exp((1-(1/g))*(x[1] + x[2]*u.std))))/(1-(1/g))
                                           
                                         }, control = list(eval.max = 5000, iter.max = 5000))
                  tau.opt_gamma = opt.out_gamma$par
                }
              }
              
              wts.SAEnet_g = pmin(exp(tau.opt_gamma[1] + tau.opt_gamma[2]*u.std), 1e+100)
              
              if(is.foldid){
                
                cv.SAEnet_lambda2.g = 
                  gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                    method = 'ls', pred.loss = 'loss',
                                    standardize = std.X, pf = wts.SAEnet_g,
                                    foldid = foldid, lambda2 = lambda2)
                
              }else{
                
                cv.SAEnet_lambda2.g = 
                  gcdnet::cv.gcdnet(x = X.std, y = y.std,
                                    method = 'ls', pred.loss = 'loss',
                                    standardize = std.X, pf = wts.SAEnet_g,
                                    nfolds = nfolds, lambda2 = lambda2)
              }
              
              c(min(cv.SAEnet_lambda2.g$cvm), cv.SAEnet_lambda2.g$lambda.min, g)
            }
            
            opt.tuning.id_lambda2 = max(which(cv.SAEnet_lambda2[,1]==min(cv.SAEnet_lambda2[,1])))
          }
          
          c(cv.SAEnet_lambda2[opt.tuning.id_lambda2,], lambda2)
        }
        
        opt.tuning.id = max(which(cv.SAEnet[,1]==min(cv.SAEnet[,1])))
        opt.wts.SAEnet = pmin(exp(cv.SAEnet[opt.tuning.id,4] + tau.opt_gamma[5]*u.std), 1e+100)
        fit.SAEnet = gcdnet::gcdnet(x = X.std, y = y.std, method = 'ls',
                                    standardize = std.X, pf = opt.wts.SAEnet,
                                    lambda = cv.SAEnet[opt.tuning.id,2],
                                    lambda2 = cv.SAEnet[opt.tuning.id,6])
        
        intercept_std.out[itr.SAEnet+1] = as.numeric(fit.SAEnet$b0)
        betahat_std.out[,itr.SAEnet+1] = as.numeric(fit.SAEnet$beta)
        w.out[,itr.SAEnet+1] = opt.wts.SAEnet
        cvm.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,1]
        lambda1.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,2]
        gamma.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,3]
        lambda2.out[itr.SAEnet+1] = cv.SAEnet[opt.tuning.id,6]
      }
      
      if(verbose) print(paste('Iteration', itr.SAEnet))
    }
  }
  
  if(std.y){
    
    intercept.out = y.scale*intercept_std.out
    betahat.out = y.scale*betahat_std.out
    
  }else{
    
    intercept.out = intercept_std.out
    betahat.out = betahat_std.out
  }
  
  # return
  return(list('intercept' = intercept.out, 'betahat' = as.matrix(betahat.out), 
              'w' = w.out, 'cvm' = cvm.out,
              'lambda1' = lambda1.out, 'lambda2' = lambda2.out,
              'gamma' = gamma.out, 'nIteration' = nIteration))
  
}


