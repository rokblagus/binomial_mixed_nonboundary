
#######libraries

library(MASS)
library(logistf) #this is an issue on tools!
#library(brglm2)
library(optimx)
library(lme4)
library(numDeriv)
library(glmmTMB)
library(blme)
library(mvtnorm)

#######rok's functions and aux

expit<-function(x) 1/(1+exp(-x))
"%^%" <- function(x, n)   with(eigen(x), vectors %*% (values^n * t(vectors)))



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


  get_theta_tmb<-function(x){

    vr1=x[1]
    vr2=x[2]
    rho=x[3]

    theta1=log(vr1)/2
    theta2=log(vr2)/2
    theta0=sign(rho)*rho/sqrt(1-rho**2)
    c(theta1,theta2,theta0)
  }

beta_s<-c(b0,betas)
theta_s<-get_theta_tmb(c(diag(corr_mat),cov2cor(corr_mat)[1,2]))

  fit_ml1<-try(glmmTMB(Y~-1+X+(-1+Z|grouping),data=xdf,family="binomial",
                      start=list(beta=beta_s,theta=theta_s),
                      control = glmmTMBControl(optCtrl = list(iter.max=0, eval.max=0),rank_check ="skip",conv_check="skip")
  ),silent=TRUE)

  res_se<-sqrt(diag(vcov(fit_ml1)$cond))


 
 write(res_se,paste0("results/ses",nms,".txt"),ncolumns=length(res_se),append = TRUE,sep="\t")
 

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




