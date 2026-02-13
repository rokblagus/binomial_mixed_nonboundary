

library(MASS)
library(brglm2)
library(optimx)
library(lme4)
library(numDeriv)
library(glmmTMB)
library(blme)
library(mvtnorm)
library(xtable)


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

rok.iter.onestep<-function(data,nu,psi,weight.c=1/2,penX=TRUE,penZ=FALSE, penTina=TRUE,fit.TMB=TRUE,REML=FALSE,plot.coef.path=FALSE,fit.frth=FALSE,nAGQ=1,optimizer=c("bobyqa", "Nelder_Mead"),typehat=2){
  
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




source("mv_MSPAL_illustration.R")  

dd<-read.table("dataBirds.txt",header=T,dec=",",sep="\t",fill=TRUE)


X <- model.matrix(parasites~migration+food, dd)
Z<-model.matrix(parasites~migration, dd)


 

 
xdf<-list(Y=dd$parasites,X=X,Z=Z,grouping=as.numeric(as.factor(dd$Phylogenetic.Tomi)),M=rep(1,nrow(X)))




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


##show the plot that illustrates how we set tau
taus<-seq(from=0,to=1,by=0.01)
liks<-rep(NA,length(taus))
zz=0
for (jjj in taus){
  zz=zz+1
liks[zz]<-tau_finder_maxi(jjj)
}

plot(taus,-liks,type="l")
abline(h=-liks[1]+lim)
abline(v=tau2)

df_taus_1<-data.frame(taus=taus,liks=liks,lim=lim,opttau=tau2)

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

li<-ee$values+tau1*(lm-ee$values)
psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
nu0=2*q-1

fit_rok<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,
                      plot.coef.path=TRUE,fit.frth=FALSE),silent=TRUE)


##show the plot that illustrates how we set tau
taus<-seq(from=0,to=1,by=0.01)
liks<-rep(NA,length(taus))
zz=0
for (jjj in taus){
  zz=zz+1
  liks[zz]<-tau_finder_maxi(jjj)
}

plot(taus,-liks,type="l")
abline(h=-liks[1]+lim)
abline(v=tau2)

df_taus<-data.frame(taus=taus,liks=liks,lim=lim,opttau=tau1)



##full iteration, firth


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
                       plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)







lik0<-get_cond_lik(fit_rok1$fit,xdf$Y)





lim<-get_lim(fit_rok1$fit,xdf$Y,xdf$X) #1SE rule 




tau_finder_maxi<-function(tau){
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,
                        plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)
  -get_cond_lik(fit_rok$fit,xdf$Y) 
  
}

tau_finder<-function(tau){
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,
                        plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)
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


tau1_f<-opt_tau

li<-ee$values+tau1_f*(lm-ee$values)

psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
nu0=2*q-1

fit_rok_frth<-try(rok.iter(xdf,nu=nu0,psi=psi0,weight.c = mv_multiplier(xdf)/2,
                      plot.coef.path=TRUE,fit.frth=TRUE),silent=TRUE)


##show the plot that illustrates how we set tau
taus<-seq(from=0,to=1,by=0.01)
liks<-rep(NA,length(taus))
zz=0
for (jjj in taus){
  zz=zz+1
  liks[zz]<-tau_finder_maxi(jjj)
}

plot(taus,-liks,type="l")
abline(h=-liks[1]+lim)
abline(v=tau2)

df_taus_f<-data.frame(taus=taus,liks=liks,lim=lim,opttau=tau1_f)

####plot cond liks by tau, with limits and optimal tau values in a single 3x1 plot

df_taus_1$algorithm<-"MPL-A2(1)"
df_taus$algorithm<-"MPL-A1"
df_taus_f$algorithm<-"MPL-A2"

df_taus_a<-rbind(df_taus_1,df_taus,df_taus_f)

library(dplyr)

df_ref <- df_taus_a %>%
  group_by(algorithm) %>%
  arrange(taus) %>%
  summarise(
    opttau = first(opttau),
    hline = -first(liks) + first(lim),
    .groups = "drop"
  )

library(ggplot2)

g_lim<-ggplot(df_taus_a, aes(x = taus, y = -liks)) +
  geom_line() +
  geom_vline(
    data = df_ref,
    aes(xintercept = opttau),
    linetype = "dashed",
    color = "blue"
  ) +
  geom_hline(
    data = df_ref,
    aes(yintercept = hline),
    linetype = "dashed",
    color = "red"
  ) +
  facet_wrap(~ algorithm, ncol = 3, scales = "free_y") +
  labs(
    x = expression(tau),
    y = "Conditional log-likelihood"
  ) +
  theme_minimal()

pdf("lims.pdf",height=4,width=9)
g_lim
dev.off()


####plot coef paths for all methods
 
m1<-matrix(unlist(fit_rok_frth$cfs),ncol=length(fit_rok_frth$cfs[[1]]),byrow = TRUE)
  
m2<-matrix(unlist(fit_rok$cfs),ncol=length(fit_rok$cfs[[1]]),byrow = TRUE)

m3<-matrix( fit_rok_1$fit$sdr$par.fixed,ncol=length(fit_rok_1$fit$sdr$par.fixed),nrow=max(c(nrow(m1),nrow(m2))),byrow = TRUE)



library(ggplot2)
library(tidyr)
library(dplyr)

# Convert matrices to data frames
df1 <- as.data.frame(m1)
df2 <- as.data.frame(m2)
df3 <- as.data.frame(m3)

# Add row index
df1$row <- 1:nrow(df1)
df2$row <- 1:nrow(df2)
df3$row <- 1:nrow(df3)

# Add source label
df1$matrix <- "MPL-A2"
df2$matrix <- "MPL-A1"
df3$matrix <- "MPL-A2(1)"

# Combine
df <- bind_rows(df1, df2,df3 )

# Pivot to long format
df_long <- pivot_longer(
  df,
  cols = -c(row, matrix),
  names_to = "column",
  values_to = "value"
)


df_long$column[df_long$column=="V1"]<-"b0"
df_long$column[df_long$column=="V2"]<-"b1"
df_long$column[df_long$column=="V3"]<-"b2"
df_long$column[df_long$column=="V4"]<-"b3"
df_long$column[df_long$column=="V5"]<-"b4"
df_long$column[df_long$column=="V6"]<-"theta1"

df_long$column[df_long$column=="V7"]<-"theta2"
df_long$column[df_long$column=="V8"]<-"theta3"
 

 
 
limits <- df_long %>%
  group_by(column) %>%
  summarise(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE)
  ) %>%
  mutate(
    center = (ymin + ymax) / 2,
    half_range = (ymax - ymin) / 2,
    ymin = center - 10 * half_range,
    ymax = center + 10* half_range
  )


#rather set to 0.01SE!

ses<-sqrt(diag(fit_rok$fit$sdr$cov.fixed))

limits$ymin<-limits$center-0.01*ses/2
limits$ymax<-limits$center+0.01*ses/2

dfpallete <- c(   "#CC6677","#AA4499","#88CCEE" )

g1<-ggplot(df_long, aes(x = row, y = value, color = matrix)) +
  geom_line() +
  geom_blank(
    data = limits,
    aes(y = ymin),
    inherit.aes = FALSE
  ) +
  geom_blank(
    data = limits,
    aes(y = ymax),
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values=dfpallete)+
  facet_wrap(~ column, scales = "free_y",nrow=2,ncol=4) +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)
  ) +
  labs(
    x = "Iteration",
    y = "",
    color = "Algorithm"
  ) +
  theme_minimal()


pdf("figure_real.pdf", width = 9/3*4, height = 6.5/3*2)
g1
dev.off()


#make table (for all results)

fit_ml<-try(glmmTMB(Y~-1+X+(-1+Z|grouping),data=xdf,family="binomial"),silent=TRUE)
fit_blm.t<-try(bglmer(Y~-1+X+(-1+Z|grouping),data=xdf,
                      family="binomial",fixef.prior = "t"),silent=TRUE)
fit_blm.n<-try(bglmer(Y~-1+X+(-1+Z|grouping),data=xdf,
                      family="binomial",fixef.prior = "normal"),silent=TRUE)

library(optimx)
library(lme4)
library(numDeriv)

np<-8

fit_kos <- try(mv_get_MSPAL(start = rep(0, np),
                            data = xdf,
                            mult=c(mv_multiplier(xdf), mv_multiplier(xdf)),
                            nAGQ = 1,
                            pen_log_sigma = function(x) sum(nHuber(x, a=1)),
                            method_name = "MSPL[KosWeight, KosWeight]"  ),silent=TRUE)



mod<-fit_ml

#sd for   RE

vr<-mod$sdr$cov.fixed[(7:8)-1,(7:8)-1]
est<-mod$sdr$par.fixed[(7:8)-1]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[9-1,9-1]
est<-mod$sdr$par.fixed[9-1]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:5,1:5]
est<-mod$sdr$par.fixed[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_ml<-c(fxf,
           paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
           paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
           
           paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
           
)


mod<-fit_rok_1$fit

#sd for   RE

vr<-mod$sdr$cov.fixed[(7:8)-1,(7:8)-1]
est<-mod$sdr$par.fixed[(7:8)-1]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[9-1,9-1]
est<-mod$sdr$par.fixed[9-1]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:5,1:5]
est<-mod$sdr$par.fixed[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_rok_1<-c(fxf,
            paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
            paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
            
            paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
            
)

 

mod<-fit_rok_frth$fit

#sd for   RE

vr<-mod$sdr$cov.fixed[(7:8)-1,(7:8)-1]
est<-mod$sdr$par.fixed[(7:8)-1]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[9-1,9-1]
est<-mod$sdr$par.fixed[9-1]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:5,1:5]
est<-mod$sdr$par.fixed[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_rok_f<-c(fxf,
              paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
              paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
              
              paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
              
)


mod<-fit_rok$fit

#sd for   RE

vr<-mod$sdr$cov.fixed[(7:8)-1,(7:8)-1]
est<-mod$sdr$par.fixed[(7:8)-1]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[9-1,9-1]
est<-mod$sdr$par.fixed[9-1]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:5,1:5]
est<-mod$sdr$par.fixed[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_rok<-c(fxf,
              paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
              paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
              
              paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
              
)



mod<-fit_blm.t

#sd for 1st RE

h <- mod@optinfo$derivs$Hessian
h <- solve(h)
v <- forceSymmetric(h + t(h))

mod@optinfo$val



#sd for  RE

vr<-v[(2:4)-1,(2:4)-1]
est<-mod@optinfo$val[(2:4)-1]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"theta1",vr)
s3.2<-deltaMethod(est,"sqrt(theta2^2+theta3^2)",vr)

s3.12<-deltaMethod(est,"theta2/sqrt(theta2^2+theta3^2)",vr)


#fixef
vr<-vcov(mod)[1:5,1:5]
est<-fixef(mod)[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_bglmert<-c(fxf,
                
                paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
                paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
                
                paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
                
)

mod<-fit_blm.n

#sd for 1st RE

h <- mod@optinfo$derivs$Hessian
h <- solve(h)
v <- forceSymmetric(h + t(h))

mod@optinfo$val



#sd for  RE

vr<-v[(2:4)-1,(2:4)-1]
est<-mod@optinfo$val[(2:4)-1]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"theta1",vr)
s3.2<-deltaMethod(est,"sqrt(theta2^2+theta3^2)",vr)

s3.12<-deltaMethod(est,"theta2/sqrt(theta2^2+theta3^2)",vr)


#fixef
vr<-vcov(mod)[1:5,1:5]
est<-fixef(mod)[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_bglmern<-c(fxf,
                
                paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
                paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
                
                paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
                
)

#kos
fit<-fit_kos$out   

fixef <- fit[1:5, 1]
ses<-fit[1:5, 2]




l11 <- exp(fit[6, 1])

l12 <- fit[7, 1]

l22 <- exp(fit[8, 1])

ch <- matrix(c(l11,0,l12,l22), ncol=2)



vrcv<- try(t(ch)%*%ch,silent=TRUE)

sds<-sqrt(diag(vrcv))
rho<-cov2cor(vrcv)[1,2]



#for sds 1, this should work (note that l11=sds[1])

est<-fit[6, 1]
vr<-fit[6, 2]**2
names(est)<-names(vr)<-"theta1"
s3.1<-deltaMethod(est,"exp(theta1)",vr)

#for sds2 and rho we would need the entire covariance matrix, we now have it!

est<-fit[7:8, 1]
vr<-fit_kos$vcv[7:8,7:8]
names(est)<-c("theta1","theta2")
colnames(vr)<-rownames(vr)<-c("theta1","theta2")
s3.2<-deltaMethod(est,"sqrt(theta1^2+exp(theta2)^2)",vr)

#rho:
s3.3<-deltaMethod(est,"theta1/sqrt(theta1^2+exp(theta2)^2)",vr)


resi_kos<-c(paste0(round(fixef,3)," (",round(ses,3),")"),
            
            paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
            paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
            
            paste0(round(s3.3$Estimate,3)," (",round(s3.3$SE,3),")")
            
)
 




#as in paper
res<- cbind(resi_ml,resi_rok,resi_rok_f,resi_rok_1,resi_bglmert,resi_bglmern,resi_kos)
colnames(res)<-c("ML","MPL-A1","MPL-A2","MPL-A2(1)","BMT","BMN","SPM")


rownames(res)<-c("Intercept","migration","food[granivourous]","food[all]","food[insectivourous]"
                 ,"plyo (Intercept)","plyo (migration)","rho")

resx<-xtable(res)

############################


#####model with 2 REs

rok.iter.2<-function(data,nu1,psi1,nu2,psi2,weight.c=1/2,tol=1e-6,maxIter=25,penX=TRUE,penZ=FALSE,fit.TMB=TRUE,REML=FALSE,plot.coef.path=FALSE,fit.frth=TRUE,nAGQ=1,optimizer=c("bobyqa", "Nelder_Mead"),typehat=2){
  if (!is.null(psi1))  pd1<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi1,nu1,const=1e8,param="variance",link_fun=function(x) 1/(1+exp(-x)))
  if (!is.null(psi2))  pd2<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi2,nu2,const=1e8,param="variance",link_fun=function(x) 1/(1+exp(-x)))
  
  xdf=data
  N1<-length(unique(xdf$grouping1))
  N2<-length(unique(xdf$grouping2))
  ##here the iterative procedure starts:
  if(penX==TRUE){
    if (fit.frth==TRUE){
      #firth<-logistf(Y~-1+X,data=xdf,pl=FALSE,firth=TRUE)
      firth<-glm(Y~-1+X,data=xdf,family=binomial,method="brglmFit",type="MPL_Jeffreys",a=mv_multiplier(xdf)/2)
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
      #hi<-diag(H)
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
      
      
      if (is.null(psi1)&!is.null(psi2)){
        X=rbind(xdf$X,xdf$X,xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)))
        if (penZ==FALSE){
          Z2=rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(xdf$Z2)*2),pd2$data$Z) #no penalty on Z
          Z1=rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(xdf$Z1)*2),matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)))
        } else {
          Z2=rbind(xdf$Z2,xdf$Z2,xdf$Z2,pd2$data$Z) #georg, does not really make sense
          Z1=rbind(xdf$Z1,xdf$Z1,xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)) )
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y,pd2$data$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c,pd2$data$nn)
        grouping2=c(xdf$grouping2,
                    xdf$grouping2,xdf$grouping2,
                    max(xdf$grouping2)+pd2$data$grouping)
        
        grouping1=c(xdf$grouping1,
                    xdf$grouping1,xdf$grouping1,
                    max(xdf$grouping1)+pd2$data$grouping)
      }
      
      if (!is.null(psi1)&is.null(psi2)){            
        X=rbind(xdf$X,xdf$X,xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
        if (penZ==FALSE){
          Z2=rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(xdf$Z2)*2),matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z))) #no penalty on Z
          Z1=rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(xdf$Z1)*2),pd1$data$Z)
        } else {
          Z2=rbind(xdf$Z2,xdf$Z2,xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd2$data$Z))) #georg, does not really make sense
          Z1=rbind(xdf$Z1,xdf$Z1,xdf$Z1,pd1$data$Z )
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y,pd1$data$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c,pd1$data$nn)
        grouping2=c(xdf$grouping2,
                    xdf$grouping2,xdf$grouping2,
                    max(xdf$grouping2)+pd1$data$grouping)
        
        grouping1=c(xdf$grouping1,
                    xdf$grouping1,xdf$grouping1,
                    max(xdf$grouping1)+pd1$data$grouping)
      }
      
      
      if (is.null(psi1)&is.null(psi2)){            
        X=rbind(xdf$X,xdf$X,xdf$X)
        if (penZ==FALSE){
          Z2=rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(xdf$Z2)*2)) #no penalty on Z
          Z1=rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(xdf$Z1)*2))
        } else {
          Z2=rbind(xdf$Z2,xdf$Z2,xdf$Z2) #georg, does not really make sense
          Z1=rbind(xdf$Z1,xdf$Z1,xdf$Z1 )
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c)
        grouping2=c(xdf$grouping2,
                    xdf$grouping2,xdf$grouping2)
        
        grouping1=c(xdf$grouping1,
                    xdf$grouping1,xdf$grouping1)
      }
      
      
      if (!is.null(psi1)&!is.null(psi2)){
        X=rbind(xdf$X,xdf$X,xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)),matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
        if (penZ==FALSE){
          Z2=rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(xdf$Z2)*2),pd2$data$Z,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z))) #no penalty on Z
          Z1=rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(xdf$Z1)*2),matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)),pd1$data$Z)
        } else {
          Z2=rbind(xdf$Z2,xdf$Z2,xdf$Z2,pd2$data$Z,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z))) #georg, does not really make sense
          Z1=rbind(xdf$Z1,xdf$Z1,xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)),pd1$data$Z )
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y,pd2$data$Y,pd1$data$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c,pd2$data$nn,pd1$data$nn)
        grouping2=c(xdf$grouping2,
                    xdf$grouping2,xdf$grouping2,
                    max(xdf$grouping2)+pd2$data$grouping,max(max(xdf$grouping2)+pd2$data$grouping)+pd1$data$grouping)
        
        grouping1=c(xdf$grouping1,
                    xdf$grouping1,xdf$grouping1,
                    max(xdf$grouping1)+pd2$data$grouping,max(max(xdf$grouping1)+pd2$data$grouping)+pd1$data$grouping)
      }
      
      pseudo<-list(X=X,
                   Z1=Z1,Z2=Z2,
                   Y=Y,
                   M=M,
                   grouping1=grouping1,grouping2=grouping2
      )
      
      
      if (fit.TMB==FALSE){
        gmm.fit0<-glmer(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), data = pseudo, family = binomial(link="logit"),weights=M ,nAGQ=nAGQ,control = glmerControl(optimizer=optimizer))
      } else {
        gmm.fit0<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), data = pseudo, family = binomial(link="logit"),weights=M, se = TRUE, verbose = FALSE, doFit = TRUE, REML = REML)
      }
    } else {
      gmm.fit0<-gmm.fit
    }
    if (fit.TMB==FALSE){
      cfs[[nIter]]<-fixef(gmm.fit0)
    } else {
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
        #hi<-diag(H)
        Xw <- (W %^% 0.5) %*% xdf$X
        qrXw <- qr(Xw)
        H <- tcrossprod(qr.Q(qrXw))
        hi<-diag(H)
      }
    } else {
      hi<-rep(0,nrow(xdf$X))
    }
    
    if (!is.null(psi1)&is.null(psi2))  pseudo$M=c(xdf$M,hi*weight.c,hi*weight.c,pd1$data$nn)
    if (is.null(psi1)&!is.null(psi2))  pseudo$M=c(xdf$M,hi*weight.c,hi*weight.c,pd2$data$nn)
    if (!is.null(psi1)&!is.null(psi2))  pseudo$M=c(xdf$M,hi*weight.c,hi*weight.c,pd2$data$nn,pd1$data$nn)
    if (is.null(psi1)&is.null(psi2))  pseudo$M=c(xdf$M,hi*weight.c,hi*weight.c)
    
    #pseudo<-list(X=X,
    #             Z1=Z1,Z2=Z2,
    #             Y=Y,
    #             M=M,
    #             grouping1=grouping1,grouping2=grouping2
    #)
    
    
    if (fit.TMB==FALSE){
      gmm.fit<-glmer(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), data = pseudo, family = binomial(link="logit"),weights=M,nAGQ=nAGQ,control = glmerControl(optimizer=optimizer) )
      tol2<-max(c(max(abs(fixef(gmm.fit0)-fixef(gmm.fit))),
                  max(abs(ranef(gmm.fit0)$grouping1[1:N1,]-ranef(gmm.fit)$grouping1[1:N1,])),
                  max(abs(ranef(gmm.fit0)$grouping2[1:N2,]-ranef(gmm.fit)$grouping2[1:N2,]))))
      
      
    } else {
      gmm.fit<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), data = pseudo, family = binomial(link="logit"),weights=M, se = TRUE, verbose = FALSE, doFit = TRUE, REML = REML)
      
      tol2<-max(c(max(abs(fixef(gmm.fit0)$cond-fixef(gmm.fit)$cond)),
                  max(abs(ranef(gmm.fit0)[[1]]$grouping1[1:N1,]-ranef(gmm.fit)[[1]]$grouping1[1:N1,])),
                  max(abs(ranef(gmm.fit0)[[1]]$grouping2[1:N2,]-ranef(gmm.fit)[[1]]$grouping2[1:N2,]))))
      
    }
    
    
    if (tol2<tol) flag=TRUE
  }
  
  if (fit.TMB==FALSE){
    cfs[[nIter+1]]<-fixef(gmm.fit)
  } else {
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
  
  
  #gmm.fit
  
  #refit to get correct SEs only works for TMB!
  
  if (is.null(psi1)&!is.null(psi2)){
    X=rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)))
    
    Z2=rbind(xdf$Z2,pd2$data$Z)
    Z1=rbind(xdf$Z1,matrix(0,ncol(xdf$Z1),nrow=nrow(pd2$data$Z)))
    Y=c(xdf$Y,pd2$data$Y)
    M=c(xdf$M,pd2$data$nn)
    grouping2=c(xdf$grouping2,
                max(xdf$grouping2)+pd2$data$grouping)
    grouping1=c(xdf$grouping1,
                max(xdf$grouping1)+pd2$data$grouping)
    
  }
  
  if (!is.null(psi1)&is.null(psi2)){
    X=rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
    
    Z2=rbind(xdf$Z2,matrix(0,ncol(xdf$Z2),nrow=nrow(pd1$data$Z)))
    Z1=rbind(xdf$Z1,pd1$data$Z)
    Y=c(xdf$Y,pd1$data$Y)
    M=c(xdf$M,pd1$data$nn)
    grouping2=c(xdf$grouping2,
                max(xdf$grouping2)+pd1$data$grouping)
    grouping1=c(xdf$grouping1,
                max(xdf$grouping1)+pd1$data$grouping)
    
  }
  
  if (is.null(psi1)&is.null(psi2)){
    X=rbind(xdf$X)
    
    Z2=rbind(xdf$Z2)
    Z1=rbind(xdf$Z1)
    Y=c(xdf$Y)
    M=c(xdf$M)
    grouping2=c(xdf$grouping2)
    grouping1=c(xdf$grouping1)
  }
  
  if (!is.null(psi1)&!is.null(psi2)){
    X=rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)),matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
    
    Z2=rbind(xdf$Z2,pd2$data$Z,matrix(0,ncol(xdf$Z2),nrow=nrow(pd1$data$Z)))
    Z1=rbind(xdf$Z1,matrix(0,ncol(xdf$Z1),nrow=nrow(pd2$data$Z)),pd1$data$Z)
    Y=c(xdf$Y,pd2$data$Y,pd1$data$Y)
    M=c(xdf$M,pd2$data$nn,pd1$data$nn)
    grouping2=c(xdf$grouping2,
                max(xdf$grouping2)+pd2$data$grouping,max(max(xdf$grouping2)+pd2$data$grouping)+pd1$data$grouping)
    grouping1=c(xdf$grouping1,
                max(xdf$grouping1)+pd2$data$grouping,max(max(xdf$grouping1)+pd2$data$grouping)+pd1$data$grouping)
    
  }
  
  xdfa<-list(X=X,
             Z1=Z1,Z2=Z2,
             Y=Y,
             M=M,
             grouping1=grouping1,grouping2=grouping2
  )
  
  gmm.fit2<-glmmTMB(gmm.fit$call$formula,data=xdfa,family="binomial",
                    start=list(beta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"beta"],theta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"theta"]),
                    control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
  )
  
  #get mar lik on orig data only works for TMB!
  
  if (!is.null(psi1)|!is.null(psi2)){
    X=rbind(xdf$X )
    
    Z2=rbind(xdf$Z2 )
    Z1=rbind(xdf$Z1 )
    Y=c(xdf$Y )
    M=c(xdf$M )
    grouping2=c(xdf$grouping2 )
    grouping1=c(xdf$grouping1 )
    
    
    xdfa<-list(X=X,
               Z1=Z1,Z2=Z2,
               Y=Y,
               M=M,
               grouping1=grouping1,grouping2=grouping2
    )
    
    gmm.fit3<-glmmTMB(gmm.fit$call$formula,data=xdfa,family="binomial",
                      start=list(beta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"beta"],theta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"theta"]),
                      control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
    ) } else gmm.fit3<-gmm.fit2
  
  
  #gmm.fit2$fit$objective<-gmm.fit3$fit$objective
  list(fit=gmm.fit2,loglik=gmm.fit3$fit$objective,cfs=cfs)
}




rok.iter.2.onestep<-function(data,nu1,psi1,nu2,psi2,weight.c=1/2,penX=TRUE,penZ=FALSE,fit.TMB=TRUE,REML=FALSE,plot.coef.path=FALSE,fit.frth=TRUE,nAGQ=1,optimizer=c("bobyqa", "Nelder_Mead"),typehat=2){
  if (!is.null(psi1))  pd1<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi1,nu1,const=1e8,param="variance",link_fun=function(x) 1/(1+exp(-x)))
  if (!is.null(psi2))  pd2<-make_pseudo_data_rand_eigen_general_psi_v3_glmm(psi2,nu2,const=1e8,param="variance",link_fun=function(x) 1/(1+exp(-x)))
  
  xdf=data
  N1<-length(unique(xdf$grouping1))
  N2<-length(unique(xdf$grouping2))
  ##here the iterative procedure starts:
  if(penX==TRUE){
    if (fit.frth==TRUE){
      #firth<-logistf(Y~-1+X,data=xdf,pl=FALSE,firth=TRUE)
      firth<-glm(Y~-1+X,data=xdf,family=binomial,method="brglmFit",type="MPL_Jeffreys",a=mv_multiplier(xdf)/2)
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
      #hi<-diag(H)
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
    
     
      
      
      if (is.null(psi1)&!is.null(psi2)){
        X=rbind(xdf$X,xdf$X,xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)))
        if (penZ==FALSE){
          Z2=rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(xdf$Z2)*2),pd2$data$Z) #no penalty on Z
          Z1=rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(xdf$Z1)*2),matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)))
        } else {
          Z2=rbind(xdf$Z2,xdf$Z2,xdf$Z2,pd2$data$Z) #georg, does not really make sense
          Z1=rbind(xdf$Z1,xdf$Z1,xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)) )
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y,pd2$data$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c,pd2$data$nn)
        grouping2=c(xdf$grouping2,
                    xdf$grouping2,xdf$grouping2,
                    max(xdf$grouping2)+pd2$data$grouping)
        
        grouping1=c(xdf$grouping1,
                    xdf$grouping1,xdf$grouping1,
                    max(xdf$grouping1)+pd2$data$grouping)
      }
      
      if (!is.null(psi1)&is.null(psi2)){            
        X=rbind(xdf$X,xdf$X,xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
        if (penZ==FALSE){
          Z2=rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(xdf$Z2)*2),matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z))) #no penalty on Z
          Z1=rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(xdf$Z1)*2),pd1$data$Z)
        } else {
          Z2=rbind(xdf$Z2,xdf$Z2,xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd2$data$Z))) #georg, does not really make sense
          Z1=rbind(xdf$Z1,xdf$Z1,xdf$Z1,pd1$data$Z )
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y,pd1$data$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c,pd1$data$nn)
        grouping2=c(xdf$grouping2,
                    xdf$grouping2,xdf$grouping2,
                    max(xdf$grouping2)+pd1$data$grouping)
        
        grouping1=c(xdf$grouping1,
                    xdf$grouping1,xdf$grouping1,
                    max(xdf$grouping1)+pd1$data$grouping)
      }
      
      
      if (is.null(psi1)&is.null(psi2)){            
        X=rbind(xdf$X,xdf$X,xdf$X)
        if (penZ==FALSE){
          Z2=rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(xdf$Z2)*2)) #no penalty on Z
          Z1=rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(xdf$Z1)*2))
        } else {
          Z2=rbind(xdf$Z2,xdf$Z2,xdf$Z2) #georg, does not really make sense
          Z1=rbind(xdf$Z1,xdf$Z1,xdf$Z1 )
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c)
        grouping2=c(xdf$grouping2,
                    xdf$grouping2,xdf$grouping2)
        
        grouping1=c(xdf$grouping1,
                    xdf$grouping1,xdf$grouping1)
      }
      
      
      if (!is.null(psi1)&!is.null(psi2)){
        X=rbind(xdf$X,xdf$X,xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)),matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
        if (penZ==FALSE){
          Z2=rbind(xdf$Z2,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(xdf$Z2)*2),pd2$data$Z,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z))) #no penalty on Z
          Z1=rbind(xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(xdf$Z1)*2),matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)),pd1$data$Z)
        } else {
          Z2=rbind(xdf$Z2,xdf$Z2,xdf$Z2,pd2$data$Z,matrix(0,ncol=ncol(xdf$Z2),nrow=nrow(pd1$data$Z))) #georg, does not really make sense
          Z1=rbind(xdf$Z1,xdf$Z1,xdf$Z1,matrix(0,ncol=ncol(xdf$Z1),nrow=nrow(pd2$data$Z)),pd1$data$Z )
        }
        Y=c(xdf$Y,xdf$Y,1-xdf$Y,pd2$data$Y,pd1$data$Y)
        M=c(xdf$M,hi*weight.c,hi*weight.c,pd2$data$nn,pd1$data$nn)
        grouping2=c(xdf$grouping2,
                    xdf$grouping2,xdf$grouping2,
                    max(xdf$grouping2)+pd2$data$grouping,max(max(xdf$grouping2)+pd2$data$grouping)+pd1$data$grouping)
        
        grouping1=c(xdf$grouping1,
                    xdf$grouping1,xdf$grouping1,
                    max(xdf$grouping1)+pd2$data$grouping,max(max(xdf$grouping1)+pd2$data$grouping)+pd1$data$grouping)
      }
      
      pseudo<-list(X=X,
                   Z1=Z1,Z2=Z2,
                   Y=Y,
                   M=M,
                   grouping1=grouping1,grouping2=grouping2
      )
      
      
      if (fit.TMB==FALSE){
        gmm.fit<-glmer(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), data = pseudo, family = binomial(link="logit"),weights=M ,nAGQ=nAGQ,control = glmerControl(optimizer=optimizer))
      } else {
        gmm.fit<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2), data = pseudo, family = binomial(link="logit"),weights=M, se = TRUE, verbose = FALSE, doFit = TRUE, REML = REML)
      }
     
  
  #gmm.fit
  
  #refit to get correct SEs only works for TMB!
  
  if (is.null(psi1)&!is.null(psi2)){
    X=rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)))
    
    Z2=rbind(xdf$Z2,pd2$data$Z)
    Z1=rbind(xdf$Z1,matrix(0,ncol(xdf$Z1),nrow=nrow(pd2$data$Z)))
    Y=c(xdf$Y,pd2$data$Y)
    M=c(xdf$M,pd2$data$nn)
    grouping2=c(xdf$grouping2,
                max(xdf$grouping2)+pd2$data$grouping)
    grouping1=c(xdf$grouping1,
                max(xdf$grouping1)+pd2$data$grouping)
    
  }
  
  if (!is.null(psi1)&is.null(psi2)){
    X=rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
    
    Z2=rbind(xdf$Z2,matrix(0,ncol(xdf$Z2),nrow=nrow(pd1$data$Z)))
    Z1=rbind(xdf$Z1,pd1$data$Z)
    Y=c(xdf$Y,pd1$data$Y)
    M=c(xdf$M,pd1$data$nn)
    grouping2=c(xdf$grouping2,
                max(xdf$grouping2)+pd1$data$grouping)
    grouping1=c(xdf$grouping1,
                max(xdf$grouping1)+pd1$data$grouping)
    
  }
  
  if (is.null(psi1)&is.null(psi2)){
    X=rbind(xdf$X)
    
    Z2=rbind(xdf$Z2)
    Z1=rbind(xdf$Z1)
    Y=c(xdf$Y)
    M=c(xdf$M)
    grouping2=c(xdf$grouping2)
    grouping1=c(xdf$grouping1)
  }
  
  if (!is.null(psi1)&!is.null(psi2)){
    X=rbind(xdf$X,matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd2$data$Z)),matrix(0,ncol=ncol(xdf$X),nrow=nrow(pd1$data$Z)))
    
    Z2=rbind(xdf$Z2,pd2$data$Z,matrix(0,ncol(xdf$Z2),nrow=nrow(pd1$data$Z)))
    Z1=rbind(xdf$Z1,matrix(0,ncol(xdf$Z1),nrow=nrow(pd2$data$Z)),pd1$data$Z)
    Y=c(xdf$Y,pd2$data$Y,pd1$data$Y)
    M=c(xdf$M,pd2$data$nn,pd1$data$nn)
    grouping2=c(xdf$grouping2,
                max(xdf$grouping2)+pd2$data$grouping,max(max(xdf$grouping2)+pd2$data$grouping)+pd1$data$grouping)
    grouping1=c(xdf$grouping1,
                max(xdf$grouping1)+pd2$data$grouping,max(max(xdf$grouping1)+pd2$data$grouping)+pd1$data$grouping)
    
  }
  
  xdfa<-list(X=X,
             Z1=Z1,Z2=Z2,
             Y=Y,
             M=M,
             grouping1=grouping1,grouping2=grouping2
  )
  
  gmm.fit2<-glmmTMB(gmm.fit$call$formula,data=xdfa,family="binomial",
                    start=list(beta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"beta"],theta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"theta"]),
                    control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
  )
  
  #get mar lik on orig data only works for TMB!
  
  if (!is.null(psi1)|!is.null(psi2)){
    X=rbind(xdf$X )
    
    Z2=rbind(xdf$Z2 )
    Z1=rbind(xdf$Z1 )
    Y=c(xdf$Y )
    M=c(xdf$M )
    grouping2=c(xdf$grouping2 )
    grouping1=c(xdf$grouping1 )
    
    
    xdfa<-list(X=X,
               Z1=Z1,Z2=Z2,
               Y=Y,
               M=M,
               grouping1=grouping1,grouping2=grouping2
    )
    
    gmm.fit3<-glmmTMB(gmm.fit$call$formula,data=xdfa,family="binomial",
                      start=list(beta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"beta"],theta=gmm.fit$sdr$par.fixed[names(gmm.fit$sdr$par.fixed)%in%"theta"]),
                      control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
    ) } else gmm.fit3<-gmm.fit2
  
  
  #gmm.fit2$fit$objective<-gmm.fit3$fit$objective
  list(fit=gmm.fit2,loglik=gmm.fit3$fit$objective)
}






dd<-read.table("dataBirds.txt",header=T,dec=",",sep="\t",fill=TRUE)

 

#put into our format
X <- model.matrix(parasites~migration+food, dd)
Z2<-model.matrix(parasites~migration, dd)
Z1<-model.matrix(parasites~1, dd)
#Z<-matrix(1, nrow=length(dd$parasites), ncol=1) #only R int
xdf <- list(Y=dd$parasites,
            nn=1-dd$parasites,
            X=as.matrix(X),
            M=rep(1, nrow(dd)),
            Z1=Z1,Z2=Z2,
            grouping1=as.numeric(as.factor(dd$species)),
            grouping2=as.numeric(as.factor(dd$Phylogenetic.Tomi))
)

fit_ml<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,family="binomial")

#summarize the data

mod<-fit_ml

#sd for 1st RE
vr<-mod$sdr$cov.fixed[6,6]
est<-mod$sdr$par.fixed[6]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 2nd RE

vr<-mod$sdr$cov.fixed[7:8,7:8]
est<-mod$sdr$par.fixed[7:8]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[9,9]
est<-mod$sdr$par.fixed[9]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:5,1:5]
est<-mod$sdr$par.fixed[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_tmb<-c(fxf,
            paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
            paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
            paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
            
            paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
            
)

##

fit_mle<-glmer(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,family="binomial")


mod<-fit_mle

#sd for 1st RE

h <- mod@optinfo$derivs$Hessian
h <- solve(h)
v <- forceSymmetric(h + t(h))

mod@optinfo$val

vr<-v[1,1]
est<-mod@optinfo$val[1]
names(est)<-"theta"

s1<-deltaMethod(est,"theta",vr)

#sd for 2nd RE

vr<-v[2:4,2:4]
est<-mod@optinfo$val[2:4]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"theta1",vr)
s3.2<-deltaMethod(est,"sqrt(theta2^2+theta3^2)",vr)

s3.12<-deltaMethod(est,"theta2/sqrt(theta2^2+theta3^2)",vr)


#fixef
vr<-vcov(fit_mle)[1:5,1:5]
est<-fixef(fit_mle)[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_glmer<-c(fxf,
              paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
              paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
              paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
              
              paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
              
)

####

res<-cbind(resi_glmer,resi_tmb)
colnames(res)<-c("lme4","glmmTMB")
rownames(res)<-c("Intercept","migration","food[granivourous]","food[all]","food[insectivourous]","species (Intercept)"
                 ,"plyo (Intercept)","plyo (migration)","rho")

resx<-xtable(res)
resx




#####fit penalized models

fit_ml<-glmmTMB(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,family="binomial")
fit_blm.t<-bglmer(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,
                  family="binomial",fixef.prior = "t")
fit_blm.n<-bglmer(Y~-1+X+(-1+Z1|grouping1)+(-1+Z2|grouping2),data=xdf,
                  family="binomial",fixef.prior = "normal")


fit_rok0<-rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=NULL,psi2=NULL,weight.c = mv_multiplier(xdf)/2,plot.coef.path=TRUE,fit.frth=FALSE,tol=1e-5,maxIter = 50)


if (class(fit_rok0)=="try-error") D_est<-diag(1,2,2) else D_est<-VarCorr(fit_rok0$fit)$cond$grouping2[1:2,1:2]
q<-nrow(D_est)
tau=0
ee<-eigen(D_est)

lm<-1 #identity as the shrink target
li<-ee$values+tau*(lm-ee$values)
psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
nu0=2*q-1

fit_rok1<-rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=FALSE,tol=1e-5,maxIter = 50)



get_cond_lik<-function(model,Y){
  
  p_hat<-predict(model,type="response")[1:length(Y)]
  
  sum(  Y*log(p_hat)+(1-Y)*log(1-p_hat)  )
  
}




get_lim<-function(model,Y,X){
  
  
  p_hat<-predict(model,type="response")[1:length(Y)]
  
  dl_deta<-Y-p_hat
  
  
  grad_beta <- t(X) %*% dl_deta
  
  
  vcov_beta <- vcov(model)$cond
  if (any(is.na(vcov_beta))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) 
  if (any(diag(vcov_beta<0))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) #res from Kos
  sqrt(t(grad_beta)%*%vcov_beta%*%grad_beta)
  
  
}  





lik0<-get_cond_lik(fit_rok1$fit,xdf$Y)





lim<-get_lim(fit_rok1$fit,xdf$Y,xdf$X) #1SE rule 




tau_finder_maxi<-function(tau){
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,
                          plot.coef.path=FALSE,fit.frth=FALSE,tol=1e-5,maxIter = 15)
               ,silent=TRUE)
  -get_cond_lik(fit_rok$fit,xdf$Y) 
  
}

tau_finder<-function(tau){
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,
                          plot.coef.path=FALSE,fit.frth=FALSE,tol=1e-5,maxIter = 15),silent=TRUE)
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

fit_rok<-try(rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=FALSE,tol=1e-5,maxIter = 50)
             ,silent=TRUE)


##show the plot that illustrates how we set tau
taus<-seq(from=0,to=1,by=0.01)
liks<-rep(NA,length(taus))
zz=0
for (jjj in taus){
  zz=zz+1
  liks[zz]<-tau_finder_maxi(jjj)
}

plot(taus,-liks,type="l")
abline(h=-liks[1]+lim)
abline(v=tau2)

df_taus<-data.frame(taus=taus,liks=liks,lim=lim,opttau=tau1)


##firth

fit_rok0<-rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=NULL,psi2=NULL,weight.c = mv_multiplier(xdf)/2,plot.coef.path=TRUE,fit.frth=TRUE,tol=1e-5,maxIter = 50)


if (class(fit_rok0)=="try-error") D_est<-diag(1,2,2) else D_est<-VarCorr(fit_rok0$fit)$cond$grouping2[1:2,1:2]
q<-nrow(D_est)
tau=0
ee<-eigen(D_est)

lm<-1 #identity as the shrink target
li<-ee$values+tau*(lm-ee$values)
psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
nu0=2*q-1

fit_rok1<-rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE,tol=1e-5,maxIter = 50)



get_cond_lik<-function(model,Y){
  
  p_hat<-predict(model,type="response")[1:length(Y)]
  
  sum(  Y*log(p_hat)+(1-Y)*log(1-p_hat)  )
  
}




get_lim<-function(model,Y,X){
  
  
  p_hat<-predict(model,type="response")[1:length(Y)]
  
  dl_deta<-Y-p_hat
  
  
  grad_beta <- t(X) %*% dl_deta
  
  
  vcov_beta <- vcov(model)$cond
  if (any(is.na(vcov_beta))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) 
  if (any(diag(vcov_beta<0))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) #res from Kos
  sqrt(t(grad_beta)%*%vcov_beta%*%grad_beta)
  
  
}  





lik0<-get_cond_lik(fit_rok1$fit,xdf$Y)





lim<-get_lim(fit_rok1$fit,xdf$Y,xdf$X) #1SE rule 




tau_finder_maxi<-function(tau){
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,
                          plot.coef.path=FALSE,fit.frth=TRUE,tol=1e-5,maxIter = 15)
               ,silent=TRUE)
  -get_cond_lik(fit_rok$fit,xdf$Y) 
  
}

tau_finder<-function(tau){
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,
                          plot.coef.path=FALSE,fit.frth=TRUE,tol=1e-5,maxIter = 15),silent=TRUE)
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


tauf<-opt_tau

li<-ee$values+opt_tau*(lm-ee$values)
psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
nu0=2*q-1

fit_rok_f<-try(rok.iter.2(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE,tol=1e-5,maxIter = 50)
             ,silent=TRUE)


##show the plot that illustrates how we set tau
taus<-seq(from=0,to=1,by=0.01)
liks<-rep(NA,length(taus))
zz=0
for (jjj in taus){
  zz=zz+1
  liks[zz]<-tau_finder_maxi(jjj)
}

plot(taus,-liks,type="l")
abline(h=-liks[1]+lim)
abline(v=tau2)

df_taus_f<-data.frame(taus=taus,liks=liks,lim=lim,opttau=tauf)

##one-step


fit_rok0<-rok.iter.2.onestep(xdf,nu1=NULL,psi1=NULL,nu2=NULL,psi2=NULL,weight.c = mv_multiplier(xdf)/2,plot.coef.path=TRUE,fit.frth=TRUE)


if (class(fit_rok0)=="try-error") D_est<-diag(1,2,2) else D_est<-VarCorr(fit_rok0$fit)$cond$grouping2[1:2,1:2]
q<-nrow(D_est)
tau=0
ee<-eigen(D_est)

lm<-1 #identity as the shrink target
li<-ee$values+tau*(lm-ee$values)
psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
nu0=2*q-1

fit_rok1<-rok.iter.2.onestep(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE)



get_cond_lik<-function(model,Y){
  
  p_hat<-predict(model,type="response")[1:length(Y)]
  
  sum(  Y*log(p_hat)+(1-Y)*log(1-p_hat)  )
  
}




get_lim<-function(model,Y,X){
  
  
  p_hat<-predict(model,type="response")[1:length(Y)]
  
  dl_deta<-Y-p_hat
  
  
  grad_beta <- t(X) %*% dl_deta
  
  
  vcov_beta <- vcov(model)$cond
  if (any(is.na(vcov_beta))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) 
  if (any(diag(vcov_beta<0))) vcov_beta<-4*diag(1,ncol(X),ncol(X)) #res from Kos
  sqrt(t(grad_beta)%*%vcov_beta%*%grad_beta)
  
  
}  





lik0<-get_cond_lik(fit_rok1$fit,xdf$Y)





lim<-get_lim(fit_rok1$fit,xdf$Y,xdf$X) #1SE rule 




tau_finder_maxi<-function(tau){
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter.2.onestep(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,
                          plot.coef.path=FALSE,fit.frth=TRUE)
               ,silent=TRUE)
  -get_cond_lik(fit_rok$fit,xdf$Y) 
  
}

tau_finder<-function(tau){
  li<-ee$values+tau*(lm-ee$values)
  psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
  nu0=2*q-1
  
  fit_rok<-try(rok.iter.2.onestep(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,
                          plot.coef.path=FALSE,fit.frth=TRUE),silent=TRUE)
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


tauf1<-opt_tau

li<-ee$values+opt_tau*(lm-ee$values)
psi0<-ee$vectors%*%diag(li)%*%t(ee$vectors)*3*q
nu0=2*q-1

fit_rok_f1<-try(rok.iter.2.onestep(xdf,nu1=NULL,psi1=NULL,nu2=nu0,psi2=psi0,weight.c = mv_multiplier(xdf)/2,plot.coef.path=FALSE,fit.frth=TRUE)
               ,silent=TRUE)


##show the plot that illustrates how we set tau
taus<-seq(from=0,to=1,by=0.01)
liks<-rep(NA,length(taus))
zz=0
for (jjj in taus){
  zz=zz+1
  liks[zz]<-tau_finder_maxi(jjj)
}

plot(taus,-liks,type="l")
abline(h=-liks[1]+lim)
abline(v=tau2)

df_taus_1<-data.frame(taus=taus,liks=liks,lim=lim,opttau=tauf1)





df_taus_1$algorithm<-"MPL-A2(1)"
df_taus$algorithm<-"MPL-A1"
df_taus_f$algorithm<-"MPL-A2"

df_taus_a<-rbind(df_taus_1,df_taus,df_taus_f)

library(dplyr)

df_ref <- df_taus_a %>%
  group_by(algorithm) %>%
  arrange(taus) %>%
  summarise(
    opttau = first(opttau),
    hline = -first(liks) + first(lim),
    .groups = "drop"
  )

library(ggplot2)

g_lim<-ggplot(df_taus_a, aes(x = taus, y = -liks)) +
  geom_line() +
  geom_vline(
    data = df_ref,
    aes(xintercept = opttau),
    linetype = "dashed",
    color = "blue"
  ) +
  geom_hline(
    data = df_ref,
    aes(yintercept = hline),
    linetype = "dashed",
    color = "red"
  ) +
  facet_wrap(~ algorithm, ncol = 3, scales = "free_y") +
  labs(
    x = expression(tau),
    y = "Conditional log-likelihood"
  ) +
  theme_minimal()

pdf("lims1.pdf",height=4,width=9)
g_lim
dev.off()



#summarize the data

mod<-fit_rok$fit
#sd for 1st RE
vr<-mod$sdr$cov.fixed[6,6]
est<-mod$sdr$par.fixed[6]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 2nd RE

vr<-mod$sdr$cov.fixed[7:8,7:8]
est<-mod$sdr$par.fixed[7:8]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[9,9]
est<-mod$sdr$par.fixed[9]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:5,1:5]
est<-mod$sdr$par.fixed[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_rok<-c(fxf,
            paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
            paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
            paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
            
            paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
            
)

mod<-fit_rok_f1$fit

#sd for 1st RE
vr<-mod$sdr$cov.fixed[6,6]
est<-mod$sdr$par.fixed[6]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 2nd RE

vr<-mod$sdr$cov.fixed[7:8,7:8]
est<-mod$sdr$par.fixed[7:8]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[9,9]
est<-mod$sdr$par.fixed[9]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:5,1:5]
est<-mod$sdr$par.fixed[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_rok2<-c(fxf,
             paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
             paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
             paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
             
             paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
             
)


mod<-fit_rok_f$fit
#sd for 1st RE
vr<-mod$sdr$cov.fixed[6,6]
est<-mod$sdr$par.fixed[6]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 2nd RE

vr<-mod$sdr$cov.fixed[7:8,7:8]
est<-mod$sdr$par.fixed[7:8]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[9,9]
est<-mod$sdr$par.fixed[9]
names(est)<-names(vr)<-paste0("theta",0)

s3.12<-deltaMethod(est,"theta0/sqrt(1+theta0^2)",vr)

#fixef
vr<-mod$sdr$cov.fixed[1:5,1:5]
est<-mod$sdr$par.fixed[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_rok3<-c(fxf,
            paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
            paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
            paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
            
            paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
            
)


###

#bglmer



mod<-fit_blm.t

#sd for 1st RE

h <- mod@optinfo$derivs$Hessian
h <- solve(h)
v <- forceSymmetric(h + t(h))

mod@optinfo$val

vr<-v[1,1]
est<-mod@optinfo$val[1]
names(est)<-"theta"

s1<-deltaMethod(est,"theta",vr)

#sd for 2nd RE

vr<-v[2:4,2:4]
est<-mod@optinfo$val[2:4]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"theta1",vr)
s3.2<-deltaMethod(est,"sqrt(theta2^2+theta3^2)",vr)

s3.12<-deltaMethod(est,"theta2/sqrt(theta2^2+theta3^2)",vr)


#fixef
vr<-vcov(mod)[1:5,1:5]
est<-fixef(mod)[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_bglmert<-c(fxf,
                paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
                paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
                paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
                
                paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
                
)


mod<-fit_blm.n

#sd for 1st RE

h <- mod@optinfo$derivs$Hessian
h <- solve(h)
v <- forceSymmetric(h + t(h))

mod@optinfo$val

vr<-v[1,1]
est<-mod@optinfo$val[1]
names(est)<-"theta"

s1<-deltaMethod(est,"theta",vr)

#sd for 2nd RE

vr<-v[2:4,2:4]
est<-mod@optinfo$val[2:4]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:3)

s3.1<-deltaMethod(est,"theta1",vr)
s3.2<-deltaMethod(est,"sqrt(theta2^2+theta3^2)",vr)

s3.12<-deltaMethod(est,"theta2/sqrt(theta2^2+theta3^2)",vr)


#fixef
vr<-vcov(mod)[1:5,1:5]
est<-fixef(mod)[1:5]

fxf<-paste0(round(est,3)," (",round(sqrt(diag(vr)),3),")")

resi_bglmern<-c(fxf,
                paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
                paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
                paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),
                
                paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")
                
)
########

res<-cbind(resi_rok,resi_rok3,resi_rok2,resi_bglmern,resi_bglmert)
colnames(res)<-c("MPL-A1","MPL-A2","MPL-A2(1)","BMN","BMT")
rownames(res)<-c("Intercept","migration","food[granivorous]","food[all]","food[insectivorous]","species (Intercept)"
                 ,"plyo (Intercept)","plyo (migration)","rho")

resx<-xtable(res)
resx



m1<-matrix(unlist(fit_rok_f$cfs),ncol=length(fit_rok_f$cfs[[1]]),byrow = TRUE)

m2<-matrix(unlist(fit_rok$cfs),ncol=length(fit_rok$cfs[[1]]),byrow = TRUE)

m3<-matrix( fit_rok_f1$fit$sdr$par.fixed,ncol=length(fit_rok_f1$fit$sdr$par.fixed),nrow=max(c(nrow(m1),nrow(m2))),byrow = TRUE)



library(ggplot2)
library(tidyr)
library(dplyr)

# Convert matrices to data frames
df1 <- as.data.frame(m1)
df2 <- as.data.frame(m2)
df3 <- as.data.frame(m3)

# Add row index
df1$row <- 1:nrow(df1)
df2$row <- 1:nrow(df2)
df3$row <- 1:nrow(df3)

# Add source label
df1$matrix <- "MPL-A2"
df2$matrix <- "MPL-A1"
df3$matrix <- "MPL-A2(1)"

# Combine
df <- bind_rows(df1, df2,df3 )

# Pivot to long format
df_long <- pivot_longer(
  df,
  cols = -c(row, matrix),
  names_to = "column",
  values_to = "value"
)


df_long$column[df_long$column=="V1"]<-"b0"
df_long$column[df_long$column=="V2"]<-"b1"
df_long$column[df_long$column=="V3"]<-"b2"
df_long$column[df_long$column=="V4"]<-"b3"
df_long$column[df_long$column=="V5"]<-"b4"
df_long$column[df_long$column=="V6"]<-"theta1"

df_long$column[df_long$column=="V7"]<-"theta21"
df_long$column[df_long$column=="V8"]<-"theta22"
df_long$column[df_long$column=="V9"]<-"theta23"



limits <- df_long %>%
  group_by(column) %>%
  summarise(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE)
  ) %>%
  mutate(
    center = (ymin + ymax) / 2,
    half_range = (ymax - ymin) / 2,
    ymin = center - 10 * half_range,
    ymax = center + 10* half_range
  )


#rather set to 0.1SE!

ses<-sqrt(diag(fit_rok$fit$sdr$cov.fixed))

limits$ymin<-limits$center-0.01*ses/2
limits$ymax<-limits$center+0.01*ses/2

dfpallete <- c(   "#CC6677", "#AA4499","#88CCEE")


g1<-ggplot(df_long, aes(x = row, y = value, color = matrix)) +
  geom_line(alpha=0.75) +
  geom_blank(
    data = limits,
    aes(y = ymin),
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values=dfpallete)+
  geom_blank(
    data = limits,
    aes(y = ymax),
    inherit.aes = FALSE
  ) +
  facet_wrap(~ column, scales = "free_y",nrow=3,ncol=3) +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)
  ) +
  labs(
    x = "Iteration",
    y = "",
    color = "Algorithm"
  ) +
  theme_minimal() 


pdf("figure_real1.pdf", width = 9, height = 6.5)
g1
dev.off()


