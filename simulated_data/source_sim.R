
 

#######libraries

library(MASS)
library(brglm2)
library(optimx)
library(lme4)
library(numDeriv)
library(glmmTMB)
library(blme)
library(mvtnorm)


#######rok's functions and aux

expit<-function(x) 1/(1+exp(-x))
"%^%" <- function(x, n)   with(eigen(x), vectors %*% (values^n * t(vectors)))

make_pseudo_data_rand_eigen_general_psi_v3_glmm<-function(psi,nu,const=1e8,param="precision",link_fun=function(x) 1/(1+exp(-x))){
  if (is.null(match.arg(param,c("precision","variance")))) stop("param needs to be one of: precision,variance")
  
  q<-ncol(psi)
  if (param=="precision") cc<-(nu-q-1)/q
  if (param=="variance") cc<-(nu+q+1)/q
  
  cc<-max(c(floor(cc),1))
  true<-psi/cc
  ee<-eigen(true,TRUE)
  ui<-list()
  for (j in 1:q){
    ui[[j]]<-sqrt(ee$values[j])*ee$vectors[,j]
  }
  
  pi<-list()
  
  for (j in 1:length(ui)){
    I<-diag(rep(1,length(ui[[j]])))
    pi[[j]]<-link_fun(I%*%matrix(ui[[j]],ncol=1))
  }
  Y<-unlist(pi)
  
  id<-rep(1:q,each=q)
  n<-rep(const,length(id))
  
  Zi<-matrix(0,ncol=q,nrow=q)
  
  for (j in 1:q){
    Zi[j,j]<-1
  }
  for (j in 1:q){
    if (j==1) Z=Zi else Z<-rbind(Z,Zi)
  }
  
  fact<-cc
  if (fact>1){
    Y<-rep(Y,fact)
    n<-rep(n,fact)
    id<-rep(1:(q*fact),each=q)
    for (j in 1:(q*fact)){
      if (j==1) Z=Zi else Z<-rbind(Z,Zi)
    }
  }
  
  data0<-list(Y=Y,grouping=id,nn=n,Z=Z)
  
  
  list(data=data0)
}

rok.iter.onestep<-function(data,nu,psi,weight.c=1/2,penX=TRUE,penZ=FALSE, penTina=TRUE,fit.TMB=TRUE,REML=FALSE,plot.coef.path=FALSE,fit.frth=TRUE,nAGQ=1,optimizer=c("bobyqa", "Nelder_Mead"),typehat=2){
  
  if (penTina) pd1<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi,nu,const=1e8,param="variance",link_fun=function(x) 1/(1+exp(-x)))
  
  xdf=data
  N<-length(unique(xdf$grouping))
  ##here the iterative procedure starts:
  if(penX==TRUE){
    if (fit.frth==TRUE){
      firth<-glm(Y~-1+X,data=xdf,family=binomial,method="brglmFit",type="MPL_Jeffreys",a=mv_multiplier(xdf)/2)
      #hi<-firth$hat.diag
      pi0<-c(expit(xdf$X%*%coef(firth)))
      W<-diag(pi0*(1-pi0))
      #H<-(W%^%(1/2))%*%xdf$X%*%solve(t(xdf$X)%*%W%*%xdf$X)%*%t(xdf$X)%*%(W%^%(1/2))
      Xw <- (W %^% 0.5) %*% xdf$X
      qrXw <- qr(Xw)
      H <- tcrossprod(qr.Q(qrXw))
      hi<-diag(H)
      
    } else {
      pi0<-rep(0.5,nrow(xdf$X))
      W<-diag(pi0*(1-pi0))
      #H<-(W%^%(1/2))%*%xdf$X%*%solve(t(xdf$X)%*%W%*%xdf$X)%*%t(xdf$X)%*%(W%^%(1/2))
      Xw <- (W %^% 0.5) %*% xdf$X
      qrXw <- qr(Xw)
      H <- tcrossprod(qr.Q(qrXw))
      hi<-diag(H)
    }
  } else {
    hi<-rep(0,nrow(xdf$X))
  }
  
  ##add pseudoobsr:
  
  
  ##add pseudoobsr:
  
  
  if (penTina){
    X=rbind(xdf$X,xdf$X,xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
    if (penZ==FALSE){
      Z=rbind(xdf$Z,matrix(0,ncol=ncol(xdf$Z),nrow=nrow(xdf$Z)*2),pd1$data$Z) #no penalty on Z, note that this means that we have tina penalty on Z but we dont have extra penalty on Z through the fixef pen!
    } else {
      Z=rbind(xdf$Z,xdf$Z,xdf$Z,pd1$data$Z) #georg, does not really make sense
    }
    Y=c(xdf$Y,xdf$Y,1-xdf$Y,pd1$data$Y)
    M=c(xdf$M,hi*weight.c,hi*weight.c,pd1$data$nn)
    grouping=c(xdf$grouping,
               xdf$grouping,xdf$grouping,
               max(xdf$grouping)+pd1$data$grouping)
    #grouping=rep(xdf$grouping,3) #old version, the new makes more sense (although the res are similar)
  } else {
    X=rbind(xdf$X,xdf$X,xdf$X)
    if (penZ==FALSE){
      Z=rbind(xdf$Z,matrix(0,ncol=ncol(xdf$Z),nrow=nrow(xdf$Z)*2)) #no penalty on Z, note that this means that we dont have tina penalty on Z but we dont have extra penalty on Z through the fixef pen!
    } else {
      Z=rbind(xdf$Z,xdf$Z,xdf$Z) #georg, does not really make sense
    }
    Y=c(xdf$Y,xdf$Y,1-xdf$Y)
    M=c(xdf$M,hi*weight.c,hi*weight.c)
    grouping=c(xdf$grouping,
               xdf$grouping,xdf$grouping)
    
  }
  
  pseudo<-list(X=X,
               Z=Z,
               Y=Y,
               M=M,
               grouping=grouping
  )
  
  
  if (fit.TMB==FALSE){
    gmm.fit0<-glmer(Y~-1+X+(-1+Z|grouping), data = pseudo, family = binomial(link="logit"),weights=M ,nAGQ=nAGQ,control = glmerControl(optimizer=optimizer))
  } else {
    gmm.fit0<-glmmTMB(Y~-1+X+(-1+Z|grouping), data = pseudo, family = binomial(link="logit"),weights=M, se = TRUE, verbose = FALSE, doFit = TRUE, REML = REML)
  }
  
  
  
  
  gmm.fit<-gmm.fit0
  #gmm.fit
  if (fit.TMB){
    #refit to ogrig + pd Z
    
    if (penTina){
      X=rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
      
      Z=rbind(xdf$Z,pd1$data$Z)
      
      Y=c(xdf$Y,pd1$data$Y)
      M=c(xdf$M,pd1$data$nn)
      grouping=c(xdf$grouping,
                 max(xdf$grouping)+pd1$data$grouping)
      #grouping=rep(xdf$grouping,3) #old version, the new makes more sense (although the res are similar)
      
      
      xdfa<-list(X=X,
                 Z=Z,
                 Y=Y,
                 M=M,
                 grouping=grouping
      )
      
      gmm.fit2<-glmmTMB(gmm.fit$call$formula,data=xdfa,family="binomial",
                        start=list(beta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"beta"],theta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"theta"]),
                        control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
      )
    }
    
    X=rbind(xdf$X )
    
    Z=rbind(xdf$Z )
    
    Y=c(xdf$Y )
    M=c(xdf$M )
    grouping=c(xdf$grouping )
    #grouping=rep(xdf$grouping,3) #old version, the new makes more sense (although the res are similar)
    
    
    xdfa<-list(X=X,
               Z=Z,
               Y=Y,
               M=M,
               grouping=grouping
    )
    
    gmm.fit3<-glmmTMB(gmm.fit$call$formula,data=xdfa,family="binomial",
                      start=list(beta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"beta"],theta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"theta"]),
                      control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
    )
    
    
    #gmm.fit2$fit$objective<-gmm.fit3$fit$objective
    
    if (penTina)  list(fit=gmm.fit2,loglik=-gmm.fit3$fit$objective) else list(fit=gmm.fit3,loglik=-gmm.fit3$fit$objective) 
  } else list(fit=gmm.fit,loglik="warning: for lmer, the model was not refited, SEs are too small! Use fit.TMB=TRUE to get correct SEs!")
}




rok.iter<-function(data,nu,psi,weight.c=1/2,tol=1e-5,maxIter=15,penX=TRUE,penZ=FALSE, penTina=TRUE,fit.TMB=TRUE,REML=FALSE,plot.coef.path=FALSE,fit.frth=FALSE,nAGQ=1,optimizer=c("bobyqa", "Nelder_Mead"),typehat=2){
  
  if (penTina) pd1<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi,nu,const=1e8,param="variance",link_fun=function(x) 1/(1+exp(-x)))
  
  xdf=data
  N<-length(unique(xdf$grouping))
  ##here the iterative procedure starts:
  if(penX==TRUE){
    if (fit.frth==TRUE){
      #firth<-logistf(Y~-1+X,data=xdf,pl=FALSE,firth=TRUE)
      #hi<-firth$hat.diag
      firth<-glm(Y~-1+X,data=xdf,family=binomial,method="brglmFit",type="MPL_Jeffreys",a=mv_multiplier(xdf)/2)
      #hi<-firth$hat.diag
      pi0<-c(expit(xdf$X%*%coef(firth)))
      W<-diag(pi0*(1-pi0))
      #H<-(W%^%(1/2))%*%xdf$X%*%solve(t(xdf$X)%*%W%*%xdf$X)%*%t(xdf$X)%*%(W%^%(1/2))
      Xw <- (W %^% 0.5) %*% xdf$X
      qrXw <- qr(Xw)
      H <- tcrossprod(qr.Q(qrXw))
      hi<-diag(H)
      
    } else {
      pi0<-rep(0.5,nrow(xdf$X))
      W<-diag(pi0*(1-pi0))
      #H<-(W%^%(1/2))%*%xdf$X%*%solve(t(xdf$X)%*%W%*%xdf$X)%*%t(xdf$X)%*%(W%^%(1/2))
      Xw <- (W %^% 0.5) %*% xdf$X
      qrXw <- qr(Xw)
      H <- tcrossprod(qr.Q(qrXw))
      hi<-diag(H)
    }
  } else {
    hi<-rep(0,nrow(xdf$X))
  }
  
  ##add pseudoobsr:
  cfs<-list()
  
  
  flag=FALSE
  
  nIter<-0
  while(flag==FALSE&nIter<maxIter){
    nIter<-nIter+1
    
    ##add pseudoobsr:
    
    if (nIter==1) {
      
      if (penTina){
        X=rbind(xdf$X,xdf$X,xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
        if (penZ==FALSE){
          Z=rbind(xdf$Z,matrix(0,ncol=ncol(xdf$Z),nrow=nrow(xdf$Z)*2),pd1$data$Z) #no penalty on Z, note that this means that we have tina penalty on Z but we dont have extra penalty on Z through the fixef pen!
        } else {
          Z=rbind(xdf$Z,xdf$Z,xdf$Z,pd1$data$Z) #georg, does not really make sense
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y,pd1$data$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c,pd1$data$nn)
        
         
        
        grouping=c(xdf$grouping,
                   xdf$grouping,xdf$grouping,
                   max(xdf$grouping)+pd1$data$grouping)
        #grouping=rep(xdf$grouping,3) #old version, the new makes more sense (although the res are similar)
      } else {
        X=rbind(xdf$X,xdf$X,xdf$X)
        if (penZ==FALSE){
          Z=rbind(xdf$Z,matrix(0,ncol=ncol(xdf$Z),nrow=nrow(xdf$Z)*2)) #no penalty on Z, note that this means that we dont have tina penalty on Z but we dont have extra penalty on Z through the fixef pen!
        } else {
          Z=rbind(xdf$Z,xdf$Z,xdf$Z) #georg, does not really make sense
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c)
        
           
                grouping=c(xdf$grouping,
                   xdf$grouping,xdf$grouping)
        
      }
      
      pseudo<-list(X=X,
                   Z=Z,
                   Y=Y,
                   M=M,
                   grouping=grouping
      )
      
      
      if (fit.TMB==FALSE){
        gmm.fit0<-glmer(Y~-1+X+(-1+Z|grouping), data = pseudo, family = binomial(link="logit"),weights=M ,nAGQ=nAGQ,control = glmerControl(optimizer=optimizer))
      } else {
        gmm.fit0<-glmmTMB(Y~-1+X+(-1+Z|grouping), data = pseudo, family = binomial(link="logit"),control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)),weights=M, se = TRUE, verbose = FALSE, doFit = TRUE, REML = REML)
        }
    } else {
      gmm.fit0<-gmm.fit
    }
    if (fit.TMB==FALSE){
      cfs[[nIter]]<-fixef(gmm.fit0)
    } else {
      #cfs[[nIter]]<-fixef(gmm.fit0)$cond
      cfs[[nIter]]<-gmm.fit0$sdr$par.fixed
    }
    #xdf$lp<-predict(gmm.fit0,random.only=TRUE)[1:nrow(xdf$X)]
    #xdf$lp<-get.re(gmm.fit0,nrow(xdf$X),xdf$grouping,xdf$Z)
    
    
    if (fit.TMB==FALSE){
      pi1<-predict(gmm.fit0)[1:nrow(xdf$X)]
      pi0<-predict(gmm.fit0,re.form = NA)[1:nrow(xdf$X)]
    } else {
      pi1<-predict(gmm.fit0,newdata=xdf)
      pi0<-predict(gmm.fit0,newdata=xdf,re.form = NA)
    }
    
    xdf$lp<-pi1-pi0
    
    
    if(penX==TRUE){
      if (fit.frth==TRUE){
        #firth<-logistf(Y~-1+X+offset(lp),data=xdf,pl=FALSE,firth=TRUE)
        #hi<-firth$hat.diag
        firth<-glm(Y~-1+X+offset(lp),data=xdf,family=binomial,method="brglmFit",type="MPL_Jeffreys",a=mv_multiplier(xdf)/2)
        #hi<-firth$hat.diag
        pi0<-c(expit(xdf$X%*%coef(firth)))
        W<-diag(pi0*(1-pi0))
        #H<-(W%^%(1/2))%*%xdf$X%*%solve(t(xdf$X)%*%W%*%xdf$X)%*%t(xdf$X)%*%(W%^%(1/2))
        Xw <- (W %^% 0.5) %*% xdf$X
        qrXw <- qr(Xw)
        H <- tcrossprod(qr.Q(qrXw))
        hi<-diag(H)
        
      } else {
        if (fit.TMB==FALSE){
          if (typehat==1) pi0<-predict(gmm.fit0)[1:nrow(xdf$X)] else pi0<-predict(gmm.fit0,re.form = NA)[1:nrow(xdf$X)]
        } else {
          if (typehat==1) pi0<-predict(gmm.fit0,newdata=xdf) else pi0<-predict(gmm.fit0,newdata=xdf,re.form = NA)
        }
        pi0<-expit(pi0)
        W<-diag(pi0*(1-pi0))
        #H<-(W%^%(1/2))%*%xdf$X%*%solve(t(xdf$X)%*%W%*%xdf$X)%*%t(xdf$X)%*%(W%^%(1/2))
        Xw <- (W %^% 0.5) %*% xdf$X
        qrXw <- qr(Xw)
        H <- tcrossprod(qr.Q(qrXw))
        
        hi<-diag(H)
      }
    } else {
      hi<-rep(0,nrow(xdf$X))
    }
    
    if (penTina==TRUE){
      
       pseudo$M=c(xdf$M,hi*weight.c,hi*weight.c,pd1$data$nn)
     
      
    } else {
      
        pseudo$M=c(xdf$M,hi*weight.c,hi*weight.c)
     }
    #pseudo<-list(X=X,
    #            Z=Z,
    #             Y=Y,
    #             M=M,
    #             grouping=grouping
    #)
    
    
    if (fit.TMB==FALSE){
      gmm.fit<-glmer(Y~-1+X+(-1+Z|grouping), data = pseudo, family = binomial(link="logit"),weights=M,nAGQ=nAGQ,control = glmerControl(optimizer=optimizer) )
      tol2<-max(c(max(abs(fixef(gmm.fit0)-fixef(gmm.fit))),
                  max(abs(ranef(gmm.fit0)$grouping[1:N,]-ranef(gmm.fit)$grouping[1:N,])) ))
      
      
    } else {
       gmm.fit<-glmmTMB(Y~-1+X+(-1+Z|grouping), data = pseudo, family = binomial(link="logit"),control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)),weights=M, se = TRUE, verbose = FALSE, doFit = TRUE, REML = REML)
       
      tol2<-max(c(max(abs(fixef(gmm.fit0)$cond-fixef(gmm.fit)$cond)),
                  max(abs(ranef(gmm.fit0)[[1]]$grouping[1:N,]-ranef(gmm.fit)[[1]]$grouping[1:N,])) ))
      
    }
    
    
    if (tol2<tol) flag=TRUE
  }
  
  if (fit.TMB==FALSE){
    cfs[[nIter+1]]<-fixef(gmm.fit)
  } else {
    #cfs[[nIter+1]]<-fixef(gmm.fit)$cond
    cfs[[nIter+1]]<-gmm.fit$sdr$par.fixed
  }
  
  if (plot.coef.path==TRUE){
    pp<-length(cfs[[1]])
    ce<-matrix(unlist(cfs),byrow=TRUE,ncol=pp)
    par(mfrow=c(1,pp))
    for (i in 1:pp){
      plot(ce[,i],type="l")
    }
    
  }
  
  ###added so that we in fact have the same as what we claim in the paper!
  
   
  #gmm.fit
  if (fit.TMB){
    #refit to ogrig + pd Z
    
    if (penTina){
      X=rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
      
      Z=rbind(xdf$Z,pd1$data$Z)
      
      Y=c(xdf$Y,pd1$data$Y)
      M=c(xdf$M,pd1$data$nn)
      
      
        
      grouping=c(xdf$grouping,
                 max(xdf$grouping)+pd1$data$grouping)
      #grouping=rep(xdf$grouping,3) #old version, the new makes more sense (although the res are similar)
      
      
      xdfa<-list(X=X,
                 Z=Z,
                 Y=Y,
                 M=M,
                 grouping=grouping
      )
      
      gmm.fit2<-glmmTMB(gmm.fit$call$formula,data=xdfa,family="binomial",
                        start=list(beta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"beta"],theta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"theta"]),
                        control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
      )
      
      
    
      
      
    }
    
    X=rbind(xdf$X )
    
    Z=rbind(xdf$Z )
    
    Y=c(xdf$Y )
    M=c(xdf$M )
    
      
    grouping=c(xdf$grouping )
    #grouping=rep(xdf$grouping,3) #old version, the new makes more sense (although the res are similar)
    
    
    xdfa<-list(X=X,
               Z=Z,
               Y=Y,
               M=M,
               grouping=grouping
    )
    
    gmm.fit3<-glmmTMB(gmm.fit$call$formula,data=xdfa,family="binomial",
                      start=list(beta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"beta"],theta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"theta"]),
                      control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
    )
    
    
    #gmm.fit2$fit$objective<-gmm.fit3$fit$objective
    
    if (penTina)  list(fit=gmm.fit2,loglik=-gmm.fit3$fit$objective,cfs=cfs,hi=hi) else list(fit=gmm.fit3,loglik=-gmm.fit3$fit$objective,cfs=cfs,hi=hi) 
  } else list(fit=gmm.fit,loglik="warning: for lmer, the model was not refited, SEs are too small! Use fit.TMB=TRUE to get correct SEs!",cfs=cfs,hi=hi)
}


######KOS

source("mv_MSPAL.R") #need this file!



######## Supporting functions for data generation in simulation study #######

# Generate grouping variables and determine the number of observations in every cluster
generate_grouping <-
  function(
    # number of clusters for each level of random effects
    num_cluster,
    # a scalar or vector of length(num_cluster);
    # a fixed number of observations within clusters or the expected value for Poisson variable
    num_obs,
    # should the number of observations be sampled from Poisson distribution?
    # In that case, the num_obs is used as lambda in zero-truncated poisson
    pois = FALSE) {
    N <- num_cluster

    if (!(length(num_obs) == 1 | length(num_obs) == num_cluster)) {
      stop("Parameter num_obs should be of length 1 or num_cluster")
    }

    if (pois) {
      flag <- TRUE
      while (flag) {
        n <- rpois(N, num_obs)
        flag <- sum(n >= 1) != N
      }
    } else{
      n <- num_obs
    }

    gr <- cbind(cl = 1:N, n_j = n)
    return(gr)
  }


# Generate fixed effects model matrix according to Geroldinger et. al. (2022)

generate_fixed_effects <-
  function(
    # matrix of grouping factors and number of observations
    grouping,
    # integer scalar; how many columns in X will be considered cluster level
    cluster_level = NULL,
    # correlation matrix of covariates in X; should be identity in case of independent variables
    corr,
    # list of functions to transform to transform the variables
    trans) {
    # if (ncol(corr) < cluster_level)
    #   stop(
    #     "Number of variables according to the covariate matrix and cluster_level argument does not match."
    #   )

    chol_cov_mat <- chol(corr)
    N <- max(grouping[, 1])
    n <- grouping[, 2]

    # generate cluster level variables (they have constant values for observations within the same cluster)
    if (cluster_level > 0)
      z_cl <-
      matrix(rnorm(N * cluster_level), ncol = cluster_level)
    # generate within level variables (the values for observations differ also within cluster)
    if ((ncol(corr) - cluster_level) > 0)
      z_with <-
      matrix(rnorm(sum(n) * (ncol(corr) - cluster_level)), ncol = (ncol(corr) -
                                                                     cluster_level)) # different within cluster

    # cluster level variables are always at the beginning
    z_x <- cbind(z_cl[rep(1:N, n), ], z_with)
    # introduce the covariance between covariates
    Z_x <- z_x %*% chol_cov_mat

    # transform the standard normal variables according to defined transformations
    x <- lapply(1:ncol(Z_x), function(x)
      trans[[x]](Z_x[, x]))
    X <- cbind(1, do.call(cbind, x))

    return(X)
  }


### dimension of the random effects covariance matrix
d = 2

## define covariance matrix between standard normal distributions between fixed effects variables according to Geroldinger et. al. (2022)
q = 5
cov_mat <- matrix(0, ncol = q, nrow = q)
cov_mat[upper.tri(cov_mat)] <- c(0.5,
                                 0.5, 0,
                                 0, 0.5, 0,
                                 0, 0, -0.3, 0.5)
cov_mat <- cov_mat + t(cov_mat)
diag(cov_mat) <- 1

## define transformations of standard normal distributions for fixed effects variables according to Geroldinger et. al. (2022)
trans_x <- list(
  x1 = function(x) {
    as.integer(x < (-1))
  },
  x2 = function(x) {
    as.integer(x < (-1))
  },
  x3 = function(x) {
    as.integer(x < 0)
  },
  x4 = function(x) {
    as.integer(x >= 0.5) +
      as.integer(x >= 1.5)
  },
  x5 = function(x) {
    a <- 10 * x + 55
    q <- quantile(a, c(0.25, 0.75))
    iqr <- abs(diff(q))
    ifelse(a < q[1] - 3 * iqr, q[1] - 3 * iqr,
           ifelse(a > q[2] + 3 * iqr,  a > q[2] + 3 * iqr, a))
  }
)



get_param<-function(fit,b0,betas){

  if (class(fit)=="try-error"){rep(NA,6+3+1+6+6+6)} else {
    fixef<-fixef(fit)$cond
    vrcv<-VarCorr(fit)$cond$group[1:2,1:2]

    vars<-c(vrcv[1,],vrcv[2,2])
    rho<-cov2cor(vrcv)[1,2]

    cib0<-confint(fit)[1,1:2]

    if (sum(is.na(cib0))==0){
      hit_b0<-ifelse(cib0[1]<=b0&cib0[2]>=b0,1,0)
      w_b0<-cib0[2]-cib0[1]
    } else {
      hit_b0<- w_b0<-NA
    }

    hits<-wdths<-rep(NA,5)
    for (ii in 1:5){
      cib1<-confint(fit)[ii+1,1:2]
      if (sum(is.na(cib1))==0){
        hits[ii]<-ifelse(cib1[1]<=betas[ii]&cib1[2]>=betas[ii],1,0)
        wdths[ii]<-cib1[2]-cib1[1]
      }

    }

    c(fixef,vars,rho,hit_b0,hits,w_b0,wdths,sqrt(diag(vcov(fit)$cond)))
  }
}




get_param_bglmer<-function(fit,b0,betas){

  if (class(fit)=="try-error"){rep(NA,6+3+1+6+6+6)} else {
    fixef<-fixef(fit)
    vrcv<-VarCorr(fit)$group[1:2,1:2]

    vars<-c(vrcv[1,],vrcv[2,2])
    rho<-cov2cor(vrcv)[1,2]

    ses<-sqrt(diag(vcov(fit)))

    make_wald_ci<-function(est,se,alpha=0.05){
      cr<-qnorm(1-alpha/2)
      cbind(est-cr*se,  est+cr*se )
    }
    cis<-make_wald_ci(fixef,ses)

    cib0<-cis[1,1:2]

    if (sum(is.na(cib0))==0){
      hit_b0<-ifelse(cib0[1]<=b0&cib0[2]>=b0,1,0)
      w_b0<-cib0[2]-cib0[1]
    } else {
      hit_b0<- w_b0<-NA
    }

    hits<-wdths<-rep(NA,5)
    for (ii in 1:5){
      cib1<-cis[ii+1,1:2]
      if (sum(is.na(cib1))==0){
        hits[ii]<-ifelse(cib1[1]<=betas[ii]&cib1[2]>=betas[ii],1,0)
        wdths[ii]<-cib1[2]-cib1[1]
      }

    }

    c(fixef,vars,rho,hit_b0,hits,w_b0,wdths,ses)
  }
}


get_res_kos<-function(fit,b0,betas){
  if (class(fit)[1]=="try-error"){rep(NA,6+3+1+6+6+6)} else {

    fixef <- fit[1:6, 1]
    ses<-fit[1:6, 2]

    make_wald_ci<-function(est,se,alpha=0.05){
      cr<-qnorm(1-alpha/2)
      cbind(est-cr*se,  est+cr*se )
    }
    cis<-make_wald_ci(fixef,ses)

    cib0<-cis[1,1:2]

    if (sum(is.na(cib0))==0){
      hit_b0<-ifelse(cib0[1]<=b0&cib0[2]>=b0,1,0)
      w_b0<-cib0[2]-cib0[1]
    } else {
      hit_b0<- w_b0<-NA
    }

    hits<-wdths<-rep(NA,5)
    for (ii in 1:5){
      cib1<-cis[ii+1,1:2]
      if (sum(is.na(cib1))==0){
        hits[ii]<-ifelse(cib1[1]<=betas[ii]&cib1[2]>=betas[ii],1,0)
        wdths[ii]<-cib1[2]-cib1[1]
      }

    }



    l11 <- exp(fit[7, 1])

    l12 <- fit[8, 1]

    l22 <- exp(fit[9, 1])

    ch <- matrix(c(l11,0,l12,l22), ncol=2)



    vrcv<- try(t(ch)%*%ch,silent=TRUE)
    if (class(vrcv)[1]=="try-error") {
      vars<-rep(NA,3)
      rho<-NA
    } else {
      vars<-c(vrcv[1,],vrcv[2,2])
      rho<-cov2cor(vrcv)[1,2]
    }
    c(fixef,vars,rho,hit_b0,hits,w_b0,wdths,ses)


  }


}


get_res_kos2<-function(fit,b0,betas){
  if (class(fit)[1]=="try-error"){rep(NA,6+3+1+6+6+6)} else {

    fixef <- fit[1:6, 1]
    ses<-fit[1:6, 2]

    make_wald_ci<-function(est,se,alpha=0.05){
      cr<-qnorm(1-alpha/2)
      cbind(est-cr*se,  est+cr*se )
    }
    cis<-make_wald_ci(fixef,ses)

    cib0<-cis[1,1:2]

    if (sum(is.na(cib0))==0){
      hit_b0<-ifelse(cib0[1]<=b0&cib0[2]>=b0,1,0)
      w_b0<-cib0[2]-cib0[1]
    } else {
      hit_b0<- w_b0<-NA
    }

    hits<-wdths<-rep(NA,5)
    for (ii in 1:5){
      cib1<-cis[ii+1,1:2]
      if (sum(is.na(cib1))==0){
        hits[ii]<-ifelse(cib1[1]<=betas[ii]&cib1[2]>=betas[ii],1,0)
        wdths[ii]<-cib1[2]-cib1[1]
      }

    }



    l11 <- exp(fit[7, 1])

    l12 <- fit[8, 1]

    l22 <- exp(fit[9, 1])

    ch <- matrix(c(l11,0,l12,l22), ncol=2)



    vrcv<- try(t(ch)%*%ch,silent=TRUE)
    if (class(vrcv)[1]=="try-error") {
      vars<-rep(NA,3)
      rho<-NA
    } else {
      vars<-c(vrcv[2,2],vrcv[1,2],vrcv[1,1])
      rho<-cov2cor(vrcv)[1,2]
    }
    c(fixef,vars,rho,hit_b0,hits,w_b0,wdths,ses)


  }


}




funi<-function(ii,num_cluster,num_subj,cov_mat,trans_x,corr_mat,b0,betas) {

  gr <- generate_grouping(num_cluster,num_subj , pois = TRUE)

  flag=FALSE
  while(flag==FALSE){
    X <- generate_fixed_effects(
      grouping = gr,
      cluster_level = 2,
      corr = cov_mat,
      trans = trans_x
    )

    X[,6]<-(X[,6]-mean(X[,6]))/sd(X[,6])
    if (length(unique(X[,2]))==1|length(unique(X[,3]))==1|length(unique(X[,4]))==1|sum((X[,2]-X[,3])!=0)==0|sum((X[,2]-X[,4])!=0)==0|sum((X[,4]-X[,3])!=0)==0) flag=FALSE else flag=TRUE

  }


  n<-num_cluster

  for (i in 1:n){
    idi<-rep(gr[i,1],gr[i,2])
    if(i==1) id<-idi else id<-c(id,idi)
  }
  for (i in 1:n){
    ni<-gr[i,2]

    b<-rmvnorm(1, mean = rep(0, 2), sigma = corr_mat)


    lpfi<-X[which(id==i),]%*%matrix(c(b0,betas),ncol=1)
    lpri<-rep(b[1],ni)+X[which(id==i),6]*rep(b[2],ni)

    lp<-lpfi+lpri


    y<-rbinom(ni,size=1,prob=1/(1+exp(-lp)))


    dfi<-data.frame(y=y,x1=X[which(id==i),2],
                    x2=X[which(id==i),3],
                    x3=X[which(id==i),4],
                    x4=X[which(id==i),5],
                    x5=X[which(id==i),6],id=id[which(id==i)])
    if (i==1) df<-dfi else df<-rbind(df,dfi)
  }

  nms<-paste(num_cluster,num_subj,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")

  #saveRDS(df, file=paste0("data/",ii,nms,".Rda"))

  X<-model.matrix(~x1+x2+x3+x4+x5,data=df)
  Z<-model.matrix(~x5,data=df)

  xdf<-list(Y=df$y,X=X,Z=Z,grouping=df$id,M=rep(1,nrow(X)))


  fit_ml<-try(glmmTMB(Y~-1+X+(-1+Z|grouping),data=xdf,family="binomial"),silent=TRUE)
  fit_blm.t<-try(bglmer(Y~-1+X+(-1+Z|grouping),data=xdf,
                    family="binomial",fixef.prior = "t"),silent=TRUE)
  fit_blm.n<-try(bglmer(Y~-1+X+(-1+Z|grouping),data=xdf,
                    family="binomial",fixef.prior = "normal"),silent=TRUE)

  
 
  ####my method, single step estimator
  fit_rok0<-try(rok.iter.onestep(xdf,nu=NULL,psi=NULL,penTina = FALSE,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)
  
  if (class(fit_rok0)=="try-error") D_est<-diag(1,2,2) else D_est<-VarCorr(fit_rok0$fit)$cond$group[1:2,1:2]
  q<-nrow(D_est)
  tau=0
  ee<-eigen(D_est)

  lm<-1 #identity as the shrink target
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok1<-try(rok.iter.onestep(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)
  
  
  
  get_cond_lik<-function(model,Y){
      
    p_hat<-predict(model,type="response")[1:length(Y)]
    
    sum(  Y*log(p_hat)+(1-Y)*log(1-p_hat)  )
    
  }
  
  
  
  lik0<-get_cond_lik(fit_rok1$fit,xdf$Y)
    
  
  get_lim<-function(model,Y,X){
      
    
    p_hat<-predict(model,type="response")[1:length(Y)]
    
    dl_deta<-Y-p_hat
    
    
    grad_beta <- t(X) %*% dl_deta
    
    
    vcov_beta <- vcov(model)$cond
    if (any(is.na(vcov_beta))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) 
    if (any(diag(vcov_beta<0))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) #res from Kos
    sqrt(t(grad_beta)%*%vcov_beta%*%grad_beta)
    
    
  }  
   
  
  lim<-get_lim(fit_rok1$fit,xdf$Y,xdf$X) #1SE rule 
  
    
    
  
  tau_finder_maxi<-function(tau){
    li<-ee$values+tau*(lm-ee$values)
    psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
    nu0=2*q-1
    
    fit_rok<-try(rok.iter.onestep(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)
    -get_cond_lik(fit_rok$fit,xdf$Y) 
    
  }
  
  tau_finder<-function(tau){
    li<-ee$values+tau*(lm-ee$values)
    psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
    nu0=2*q-1
    
    fit_rok<-try(rok.iter.onestep(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)
    abs(get_cond_lik(fit_rok$fit,xdf$Y)-lik0)-lim
    
  }
  
  tol_min<-.Machine$double.eps^0.2
  tol_max<-.Machine$double.eps^0.2
  
  opti<-optimise(tau_finder_maxi,c(0,1))
  
  if (opti$minimum<tol_min){
    uni<-  try(uniroot(tau_finder,lower=0,upper = 1),silent=TRUE)
    if (class(uni)=="try-error") opt_tau<-1 else opt_tau<-uni$root
  } else {
    if ((-opti$objective-lik0)<lim) opt_tau<-opti$minimum else {
      uni<-  uniroot(tau_finder,lower=tol_min,upper = opti$minimum) 
      opt_tau<-uni$root
    }
  }
  if (opt_tau>(1-tol_max)) opt_tau<-1-tol_max
  if (opt_tau<tol_min) opt_tau<-tol_min
  
tau2<-opt_tau
  
  li<-ee$values+opt_tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok_1<-try(rok.iter.onestep(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)
  
  ##full iteration
  
  
  fit_rok0<-try(rok.iter(xdf,nu=NULL,psi=NULL,penTina = FALSE,weight.c = mv_multiplier(xdf)/2,
                         plot.coef.path=FALSE,fit.frth=FALSE),silent=TRUE)
  
  if (class(fit_rok0)=="try-error") D_est<-diag(1,2,2) else D_est<-VarCorr(fit_rok0$fit)$cond$group[1:2,1:2]
  q<-nrow(D_est)
  tau=0
  ee<-eigen(D_est)
  
  lm<-1 #identity as the shrink target
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok1<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,
                         plot.coef.path=FALSE,fit.frth=FALSE),silent=TRUE)
  
  
  
  
  
  
  
  lik0<-get_cond_lik(fit_rok1$fit,xdf$Y)
  
  
  
  
  
  lim<-get_lim(fit_rok1$fit,xdf$Y,xdf$X) #1SE rule 
  
  
  
  
  tau_finder_maxi<-function(tau){
    li<-ee$values+tau*(lm-ee$values)
    psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
    nu0=2*q-1
    
    fit_rok<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,
                          plot.coef.path=FALSE,fit.frth=FALSE),silent=TRUE)
    -get_cond_lik(fit_rok$fit,xdf$Y) 
    
  }
  
  tau_finder<-function(tau){
    li<-ee$values+tau*(lm-ee$values)
    psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
    nu0=2*q-1
    
    fit_rok<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,
                          plot.coef.path=FALSE,fit.frth=FALSE),silent=TRUE)
    abs(get_cond_lik(fit_rok$fit,xdf$Y)-lik0)-lim
    
  }
  
  tol_min<-.Machine$double.eps^0.2
  tol_max<-.Machine$double.eps^0.2
  
  opti<-optimise(tau_finder_maxi,c(0,1))
  
  if (opti$minimum<tol_min){
    uni<-  try(uniroot(tau_finder,lower=0,upper = 1),silent=TRUE)
    if (class(uni)=="try-error") opt_tau<-1 else opt_tau<-uni$root
  } else {
    if ((-opti$objective-lik0)<lim) opt_tau<-opti$minimum else {
      uni<-  uniroot(tau_finder,lower=tol_min,upper = opti$minimum) 
      opt_tau<-uni$root
    }
  }
  if (opt_tau>(1-tol_max)) opt_tau<-1-tol_max
  if (opt_tau<tol_min) opt_tau<-tol_min
  
  
  tau1<-opt_tau
  
  li<-ee$values+opt_tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,
                        plot.coef.path=FALSE,fit.frth=FALSE),silent=TRUE)
  
  
  
  

  #ORCALE
  psi0<-corr_mat*3*q
  nu0=2*q-1

  fit_rok_o<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=FALSE),silent=TRUE)
  fit_rok_o1<-try(rok.iter.onestep(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)

 
  

  #KOS
  np<-9

  fit_kos <- try(mv_get_MSPAL(start = rep(0, np),
                               data = xdf,
                               mult=c(mv_multiplier(xdf), mv_multiplier(xdf)),
                               nAGQ = 1,
                               pen_log_sigma = function(x) sum(nHuber(x, a=1)),
                               method_name = "MSPL[KosWeight, KosWeight]"  ),silent=TRUE)


#KOS with reversed Z
  xdf$Z<-xdf$Z[,c(2,1)]
  fit_kos2 <- try(mv_get_MSPAL(start = rep(0, np),
                              data = xdf,
                              mult=c(mv_multiplier(xdf), mv_multiplier(xdf)),
                              nAGQ = 1,
                              pen_log_sigma = function(x) sum(nHuber(x, a=1)),
                              method_name = "MSPL[KosWeight, KosWeight]"  ),silent=TRUE)

  #####

  #opt taus
  
  opt.taus<-c(tau1,tau2)
  
#ML
  resml<-get_param(fit_ml,b0=b0,betas=betas)
#BM
  resbglmern<-get_param_bglmer(fit_blm.n,b0=b0,betas=betas)
  resbglmert<-get_param_bglmer(fit_blm.t,b0=b0,betas=betas)


#OUR
  respml<-get_param(fit_rok$fit,b0=b0,betas=betas)
  resoracle<-get_param(fit_rok_o$fit,b0=b0,betas=betas)

  respml1<-get_param(fit_rok_1$fit,b0=b0,betas=betas)
  resoracle1<-get_param(fit_rok_o1$fit,b0=b0,betas=betas)


#KOS


  reskos<-get_res_kos(fit_kos,b0=b0,betas=betas)

  reskos2<-get_res_kos2(fit_kos2,b0=b0,betas=betas)


  write(opt.taus,paste0("results/taus",nms,".txt"),ncolumns=length(opt.taus),append = TRUE,sep="\t")
  
  
  write(reskos,paste0("results/kos",nms,".txt"),ncolumns=length(reskos),append = TRUE,sep="\t")
  write(reskos2,paste0("results/kos2",nms,".txt"),ncolumns=length(reskos2),append = TRUE,sep="\t")


  write(resml,paste0("results/ml",nms,".txt"),ncolumns=length(resml),append = TRUE,sep="\t")
   write(resbglmern,paste0("results/bmln",nms,".txt"),ncolumns=length(resbglmern),append = TRUE,sep="\t")
   write(resbglmert,paste0("results/bmlt",nms,".txt"),ncolumns=length(resbglmert),append = TRUE,sep="\t")


   write(respml,paste0("results/pml",nms,".txt"),ncolumns=length(respml),append = TRUE,sep="\t")
  write(resoracle,paste0("results/o",nms,".txt"),ncolumns=length(resoracle),append = TRUE,sep="\t")

  write(respml1,paste0("results/pml1",nms,".txt"),ncolumns=length(respml1),append = TRUE,sep="\t")
  write(resoracle1,paste0("results/o1",nms,".txt"),ncolumns=length(resoracle1),append = TRUE,sep="\t")


  print(ii)
}







####################
####################




##sim params, fixed!
sd_slope<-sd_int*multiplier_int

b0 <- -3.5

# regression coefficients of fixed effects
betas <- c(1, 1, 1, 0.35, -0.1)




# list of covariance matrices for random effects
corr_mat <-  matrix(c(sd_int**2, rho*sd_int*sd_slope, rho*sd_int*sd_slope, sd_slope**2), ncol = 2, nrow = 2)



B=1000




  m=0

  while(m<B){

    m=m+1

    dd<-try( funi(m,num_cluster,num_subj,cov_mat,trans_x,corr_mat,b0,betas)  ,silent=TRUE)

    if (class(dd)=="try-error") m<-m-1

  }




