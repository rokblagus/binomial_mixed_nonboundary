
###should be ok

#setwd("D:/Mixed_separation_fixed_oct_25")


library(ggplot2)
library(gridExtra)


################# 1st - compare different implementations




zz=0
for (N in c(20,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.1,0.5)){
      for (mult in sqrt(c(0.01,0.1,0.25,0.5,1))){
        for (rh in c(-0.5,-0.2)){
          
          zz=zz+1
          
          sd_int<-sqrt(vr_int)
          multiplier_int<-round(mult,7)
          
          sd_slope<-sd_int*multiplier_int
          
          b0 <- -3.5
          
          # regression coefficients of fixed effects
          betas <- c(1, 1, 1, 0.35, -0.1)
          
          
          rho=rh
          
          # list of covariance matrices for random effects
          corr_mat <-  matrix(c(sd_int**2, rho*sd_int*sd_slope, rho*sd_int*sd_slope, sd_slope**2), ncol = 2, nrow = 2)
          true_sig<-corr_mat
          
          true_pars<-c(NA,min(eigen(true_sig)$values),max(eigen(true_sig)$values),mean(eigen(true_sig)$values),
                       determinant(true_sig)$modulus)
          
          nms<-paste(N,n,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")
          
          meth<-"pml"
          
          nmsfile<-paste0("results/",meth,nms,".txt") 
          
          
          dd.f<-try(read.table(nmsfile,header=FALSE,sep="\t"),silent=TRUE)
          if (class(dd.f)[1]!="try-error"){
          
          
          
          ddmf<-as.matrix(dd.f)
          if (nrow(ddmf)>1000) ddmf[1:1000,] 
          
          meth<-"pml1"
          
          nmsfile<-paste0("results/",meth,nms,".txt") #using tau=0.05, truncation for POR  
          
          
          dd.1<-read.table(nmsfile,header=FALSE,sep="\t")
          ddm1<-as.matrix(dd.1)
          if (nrow(ddm1)>1000) ddm1[1:1000,] 
          
          
          meth<-"pmla2"
          
          nmsfile<-paste0("results/",meth,nms,".txt") 
          
          
          dd.fa2<-read.table(nmsfile,header=FALSE,sep="\t")
          
          
          
          ddmfa2<-as.matrix(dd.fa2)
          if (nrow(ddmfa2)>1000) ddmfa2[1:1000,]
          
           
          par(mfrow=c(3,3))
          
          mad95<-mad<-rmse<-rep(NA,9)
          for (mm in 1:9){
             plot(ddmf[,mm],ddm1[,mm])
              abline(a=0,b=1)  
            mad[mm]<-mean(abs(ddmf[,mm]-ddm1[,mm])/sd(ddmf[,mm]))
            rmse[mm]<-sqrt(mean((ddmf[,mm]-ddm1[,mm])**2))/sd(ddmf[,mm])
            mad95[mm] <- quantile(abs(ddmf[,mm]-ddm1[,mm])/sd(ddmf[,mm]),probs=0.95)
            
            }
          
          sq_eucl_norm<-mean((apply(abs(ddmf[,1:9]-ddm1[,1:9])/matrix(apply(ddmf[,1:9],2,sd,na.rm=TRUE),ncol=9,nrow=nrow(ddmf),byrow = TRUE),1,mean,na.rm=TRUE)))
          sq_eucl_norm_95<-quantile((apply(abs(ddmf[,1:9]-ddm1[,1:9])/matrix(apply(ddmf[,1:9],2,sd,na.rm=TRUE),ncol=9,nrow=nrow(ddmf),byrow = TRUE),1,mean,na.rm=TRUE)),probs=0.95)
          
          
          par(mfrow=c(3,3))
          
          mad952<-mad2<-rmse2<-rep(NA,9)
          for (mm in 1:9){
            plot(ddmf[,mm],ddmfa2[,mm])
            abline(a=0,b=1)  
            mad2[mm]<-mean(abs(ddmf[,mm]-ddmfa2[,mm])/sd(ddmf[,mm]))
            rmse2[mm]<-sqrt(mean((ddmf[,mm]-ddmfa2[,mm])**2))/sd(ddmf[,mm])
            mad952[mm] <- quantile(abs(ddmf[,mm]-ddmfa2[,mm])/sd(ddmf[,mm]),probs=0.95)
            
          }
          
          sq_eucl_norm2<-mean((apply(abs(ddmf[,1:9]-ddmfa2[,1:9])/matrix(apply(ddmf[,1:9],2,sd,na.rm=TRUE),ncol=9,nrow=nrow(ddmf),byrow = TRUE),1,mean,na.rm=TRUE)))
          sq_eucl_norm2_95<-quantile((apply(abs(ddmf[,1:9]-ddmfa2[,1:9])/matrix(apply(ddmf[,1:9],2,sd,na.rm=TRUE),ncol=9,nrow=nrow(ddmf),byrow = TRUE),1,mean,na.rm=TRUE)),probs=0.95)
          
          meth<-rep(c("PML-A1 vs PML-A2(1)","PML-A1 vs PML-A2"),each=length(mad952))
          cfs<-c(paste0("beta",0:5),paste0("theta",1:3))
          dfi1<-data.frame(
            N=N,n=n,vr_int=vr_int,mult=mult,rho=rh,mad95=c(mad95,mad952),coef=rep(cfs,2),mad=c(mad,mad2),
            rmse=c(rmse,rmse2),meth=meth,useful=TRUE)
          
          dfi1_2<-data.frame(
            N=N,n=n,vr_int=vr_int,mult=mult,rho=rh,sqenorm=c(sq_eucl_norm,sq_eucl_norm2),sqenorm95=c(sq_eucl_norm_95,sq_eucl_norm2_95),meth=c("SAD(PML-A1,PML-A2(1))","SAD(PML-A1,PML-A2)"),useful=TRUE)
          
          ###oracle
          
          meth<-"o"
          
          nmsfile<-paste0("results/",meth,nms,".txt") 
          
          
          dd.f<-read.table(nmsfile,header=FALSE,sep="\t")
          
          
          
          ddmf<-as.matrix(dd.f)
          if (nrow(ddmf)>1000) ddmf[1:1000,] 
          
          meth<-"o1"
          
          nmsfile<-paste0("results/",meth,nms,".txt") #using tau=0.05, truncation for POR  
          
          
          dd.1<-read.table(nmsfile,header=FALSE,sep="\t")
          ddm1<-as.matrix(dd.1)
          if (nrow(ddm1)>1000) ddm1[1:1000,] 
          
          
          meth<-"oa2"
          
          nmsfile<-paste0("results/",meth,nms,".txt") 
          
          
          dd.fa2<-read.table(nmsfile,header=FALSE,sep="\t")
          
          
          
          ddmfa2<-as.matrix(dd.fa2)
          if (nrow(ddmfa2)>1000) ddmfa2[1:1000,]
          
          
          par(mfrow=c(3,3))
          
          mad95<-mad<-rmse<-rep(NA,9)
          for (mm in 1:9){
            plot(ddmf[,mm],ddm1[,mm])
            abline(a=0,b=1)  
            mad[mm]<-mean(abs(ddmf[,mm]-ddm1[,mm])/sd(ddmf[,mm]))
            rmse[mm]<-sqrt(mean((ddmf[,mm]-ddm1[,mm])**2))/sd(ddmf[,mm])
            mad95[mm] <- quantile(abs(ddmf[,mm]-ddm1[,mm])/sd(ddmf[,mm]),probs=0.95)
            
          }
          
          sq_eucl_norm<-mean((apply(abs(ddmf[,1:9]-ddm1[,1:9])/matrix(apply(ddmf[,1:9],2,sd,na.rm=TRUE),ncol=9,nrow=nrow(ddmf),byrow = TRUE),1,mean,na.rm=TRUE)))
          sq_eucl_norm_95<-quantile((apply(abs(ddmf[,1:9]-ddm1[,1:9])/matrix(apply(ddmf[,1:9],2,sd,na.rm=TRUE),ncol=9,nrow=nrow(ddmf),byrow = TRUE),1,mean,na.rm=TRUE)),probs=0.95)
          
          par(mfrow=c(3,3))
          
          mad952<-mad2<-rmse2<-rep(NA,9)
          for (mm in 1:9){
            plot(ddmf[,mm],ddmfa2[,mm])
            abline(a=0,b=1)  
            mad2[mm]<-mean(abs(ddmf[,mm]-ddmfa2[,mm])/sd(ddmf[,mm]))
            rmse2[mm]<-sqrt(mean((ddmf[,mm]-ddmfa2[,mm])**2))/sd(ddmf[,mm])
            mad952[mm] <- quantile(abs(ddmf[,mm]-ddmfa2[,mm])/sd(ddmf[,mm]),probs=0.95)
          }
          sq_eucl_norm2<-mean((apply(abs(ddmf[,1:9]-ddmfa2[,1:9])/matrix(apply(ddmf[,1:9],2,sd,na.rm=TRUE),ncol=9,nrow=nrow(ddmf),byrow = TRUE),1,mean,na.rm=TRUE)))
          sq_eucl_norm2_95<-quantile((apply(abs(ddmf[,1:9]-ddmfa2[,1:9])/matrix(apply(ddmf[,1:9],2,sd,na.rm=TRUE),ncol=9,nrow=nrow(ddmf),byrow = TRUE),1,mean,na.rm=TRUE)),probs=0.95)
          
          
          meth<-rep(c("ORACLE-A1 vs ORACLE-A2(1)","ORACLE-A1 vs ORACLE-A2"),each=length(mad952))
          cfs<-c(paste0("beta",0:5),paste0("theta",1:3))
          dfi2<-data.frame(
            N=N,n=n,vr_int=vr_int,mult=mult,rho=rh,mad95=c(mad95,mad952),coef=rep(cfs,2),mad=c(mad,mad2),
            rmse=c(rmse,rmse2),meth=meth,useful=FALSE)
          
          dfi2_2<-data.frame(
            N=N,n=n,vr_int=vr_int,mult=mult,rho=rh,sqenorm=c(sq_eucl_norm,sq_eucl_norm2),sqenorm95=c(sq_eucl_norm_95,sq_eucl_norm2_95),meth=c("SAD(ORACLE-A1,ORACLE-A2(1))","SAD(ORACLE-A1,ORACLE-A2)"),useful=FALSE)
          
          
          dfi<-rbind(dfi1,dfi2)
          dfi_2<-rbind(dfi1_2,dfi2_2)
          if(zz==1) df<-dfi else df<-rbind(df,dfi)
          if(zz==1) df_2<-dfi_2 else df_2<-rbind(df_2,dfi_2)
          
          print(zz) 
        }
        }
        
      }
    }
  }}

txt1<-"sigma[0]^2"
txt2<-"sigma[1]^2"

df$Nn<-paste0("N=",df$N,",n=",df$n)
df$vars<-paste0(txt1,"==",df$vr_int,"*','*~",txt2,"==",df$vr_int*df$mult**2)


df$Nn2<-as.numeric(as.factor(df$Nn))


df$r2<-paste0("rho==",df$rho)
df$group<-ifelse(df$N==20,"1","2")




df$vars<-factor(df$vars)

df_2$Nn<-paste0("N=",df_2$N,",n=",df_2$n)
df_2$vars<-paste0(txt1,"==",df_2$vr_int,"*','*~",txt2,"==",df_2$vr_int*df_2$mult**2)


df_2$Nn2<-as.numeric(as.factor(df_2$Nn))


df_2$r2<-paste0("rho==",df_2$rho)
df_2$group<-ifelse(df_2$N==20,"1","2")




df_2$vars<-factor(df_2$vars)


alphap=0.75



nms <- c(
  "atop(N==20, n[i]=='S')",
  "atop(N==20, n[i]=='M')",
  "atop(N==50, n[i]=='S')",
  "atop(N==50, n[i]=='M')"
)

df$r2<-factor(df$r2,levels=c("rho==-0.5","rho==-0.2"))
df_2$r2<-factor(df_2$r2,levels=c("rho==-0.5","rho==-0.2"))

p1b1<-ggplot(df[df$rho==-0.5,],aes(group=interaction(meth, group)))+
  geom_point(aes(x=Nn2,y=mad95,colour =  meth),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=mad95,colour =  meth),linewidth=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(
    coef+r2 ~ vars,
    labeller = labeller(vars = label_parsed,r2=label_parsed)#,
#    scales = "free_y"
  ) +
  theme_light()+
  scale_x_continuous(
    name = "",
    breaks = c(1,2,3,4),
    labels = function(x) parse(text = nms[x]),
    limits = c(0.75,4.25)
  ) +
  ylab("Standardized absolute difference (SAD), 95th percentile")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=c(0.01))+
  # scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p1b2<-ggplot(df[df$rho==-0.2,],aes(group=interaction(meth, group)))+
  geom_point(aes(x=Nn2,y=mad95,colour =  meth),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=mad95,colour =  meth),linewidth=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(
    coef+r2 ~ vars,
    labeller = labeller(vars = label_parsed,r2=label_parsed)#,
    #    scales = "free_y"
  ) +
  theme_light()+
  scale_x_continuous(
    name = "",
    breaks = c(1,2,3,4),
    labels = function(x) parse(text = nms[x]),
    limits = c(0.75,4.25)
  ) +
  ylab("Standardized absolute difference (SAD), 95th percentile")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=c(0.01))+
  # scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b1<-ggplot(df[df$rho==-0.5,],aes(group=interaction(meth, group)))+
  geom_point(aes(x=Nn2,y=mad,colour =  meth),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=mad,colour =  meth),linewidth=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(
    r2+coef  ~ vars,
    labeller = labeller(vars = label_parsed,r2=label_parsed)#,
    #    scales = "free_y"
  ) +
  theme_light()+
  scale_x_continuous(
    name = "",
    breaks = c(1,2,3,4),
    labels = function(x) parse(text = nms[x]),
    limits = c(0.75,4.25)
  ) +
  ylab("Standardized absolute difference (SAD), mean")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=c(0.01))+
  # scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b2<-ggplot(df[df$rho==-0.2,],aes(group=interaction(meth, group)))+
  geom_point(aes(x=Nn2,y=mad,colour =  meth),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=mad,colour =  meth),linewidth=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(
    r2+coef  ~ vars,
    labeller = labeller(vars = label_parsed,r2=label_parsed)#,
    #    scales = "free_y"
  ) +
  theme_light()+
  scale_x_continuous(
    name = "",
    breaks = c(1,2,3,4),
    labels = function(x) parse(text = nms[x]),
    limits = c(0.75,4.25)
  ) +
  ylab("Standardized absolute difference (SAD), mean")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=c(0.01))+
  # scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())


p4b<-ggplot(df_2,aes(group=interaction(meth, group)))+
  geom_point(aes(x=Nn2,y=sqenorm,colour =  meth),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=sqenorm,colour =  meth),linewidth=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(
    r2  ~ vars,
    labeller = labeller(vars = label_parsed,r2=label_parsed)#,
    #    scales = "free_y"
  ) +
  theme_light()+
  scale_x_continuous(
    name = "",
    breaks = c(1,2,3,4),
    labels = function(x) parse(text = nms[x]),
    limits = c(0.75,4.25)
  ) +
  ylab("Standardized absolute difference (SAD), mean")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=c(0.01))+
  # scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p4b2<-ggplot(df_2,aes(group=interaction(meth, group)))+
  geom_point(aes(x=Nn2,y=sqenorm95,colour =  meth),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=sqenorm95,colour =  meth),linewidth=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(
    r2  ~ vars,
    labeller = labeller(vars = label_parsed,r2=label_parsed)#,
    #    scales = "free_y"
  ) +
  theme_light()+
  scale_x_continuous(
    name = "",
    breaks = c(1,2,3,4),
    labels = function(x) parse(text = nms[x]),
    limits = c(0.75,4.25)
  ) +
  ylab("Standardized absolute difference (SAD), 95th percentile")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=c(0.01))+
  # scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_lin_p1.pdf",height=25,width=22)
p1b1
dev.off()

pdf("figs/supp_lin_p2.pdf",height=25,width=22)
p1b2
dev.off()


pdf("figs/supp_sadd_p1.pdf",height=25,width=22)
p2b1
dev.off()

pdf("figs/supp_sadd_p2.pdf",height=25,width=22)
p2b2
dev.off()


pdf("figs/supp_sad.pdf",height=20,width=22)
grid.arrange(p4b,p4b2,nrow=2)
dev.off()


################# 2nd - summarize the performance of all methods

######cnvergence params

##binom
lim_sing<-1e-5 #maybe we will need to change this! this might over-estimate the prevalence, but since its zero for our methods it should be okish??
#lim_sing<-sqrt(.Machine$double.eps) #this likely does not make sense
lim_sing_rho<-1e-5 #for this we could reference Chung et al. 2015, but this would miss cases where ests of diags are zero.
#lim_sing_s<--log(sqrt(.Machine$double.eps))
lim_sing_s<--log(sqrt(1e-5))


########################


###func to sum res

ent_loss<-function(true,est){
  q<-ncol(true)
  SS<-solve(true)%*%est
  if (det(SS)<1e-90) sum(diag(SS))-log(1e-90)-q else  sum(diag(SS))-log(det(SS))-q
}



#the same as above, maybe less prone to rounding error (but it will give Inf when the estimated matrix is not PD, could be an issue for PRIAL when ML has issues!)
stein_loss_reverse <- function(Sigma_true, Sigma_est) {
  p <- nrow(Sigma_true)
  inv_true <- solve(Sigma_true)
  trace_term <- sum(diag(Sigma_est %*% inv_true))
  logdet_est <- determinant(Sigma_est, logarithm = TRUE)$modulus
  logdet_true <- determinant(Sigma_true, logarithm = TRUE)$modulus
  loss <- trace_term - logdet_est + logdet_true - p
  return(as.numeric(loss))
}

#sum of squared differences of all elements
frobenius_norm<-function(A, B ) {
  diff <- A - B
  val <- sum(diff^2)
   
    return(sqrt(val))    # Frobenius norm
   
}

#as defined in https://www.sciencedirect.com/science/article/pii/S0047259X21000749?via%3Dihub
frobenius_norm<-function(A, B) {
  diff <- A - B
  sqrt(sum(diag(diff%*%t(diff))))
}

# (seveda sta obe isti!)

get_loss<-function(Sigma,ests){
  if (sum(is.na(ests))==3) NA else {
  Sigma_est<-matrix(c(ests[1],ests[2],ests[2],ests[3]),2,2)
  #ent_loss(Sigma,Sigma_est)
  #frobenius_norm(Sigma,Sigma_est)
  stein_loss_reverse(Sigma,Sigma_est)
  }
}

get_loss2<-function(Sigma,ests){
  if (sum(is.na(ests))==3) NA else {
    Sigma_est<-matrix(c(ests[1],ests[2],ests[2],ests[3]),2,2)
    #ent_loss(Sigma,Sigma_est)
    frobenius_norm(Sigma,Sigma_est)
    #stein_loss_reverse(Sigma,Sigma_est)
  }
}

min_eigen<-function(vr1,cvr,vr2){
  if (is.na(vr1)) NA else {
    Sigma<-matrix(c(vr1,cvr,cvr,vr2),2,2)
    min(eigen(Sigma)$values)
  }
}

max_eigen<-function(vr1,cvr,vr2){
  if (is.na(vr1)) NA else {
    Sigma<-matrix(c(vr1,cvr,cvr,vr2),2,2)
    max(eigen(Sigma)$values)
  }
}


mean_eigen<-function(vr1,cvr,vr2){
  if (is.na(vr1)) NA else {
    Sigma<-matrix(c(vr1,cvr,cvr,vr2),2,2)
    mean(eigen(Sigma)$values,na.rm=TRUE)
  }
}

get_bias<-function(est,true){
  mean(est,na.rm=TRUE)-true
}

get_rmse<-function(est,true){
  sqrt( mean( (est-true)**2,na.rm=TRUE )  )
}
dfmsi<-NULL

get_determinant<-function(ests){
  if (sum(is.na(ests))==3) NA else {
    Sigma_est<-matrix(c(ests[1],ests[2],ests[2],ests[3]),2,2)
    determinant(Sigma_est,log=TRUE)$mod
    }
}

###############
#get ses for each setting so that we can judge convergence of betas as in Geroldinger 
 
zz=0
for (N in c(20,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.1,0.5)){
      for (mult in sqrt(c(0.01,0.1,0.25,0.5,1))){
        for (rh in c(-0.5,-0.2)){
          
          sd_int<-sqrt(vr_int)
          multiplier_int<-round(mult,7)

          sd_slope<-sd_int*multiplier_int

          b0 <- -3.5

          # regression coefficients of fixed effects
          betas <- c(1, 1, 1, 0.35, -0.1)


          rho=rh

          # list of covariance matrices for random effects
          corr_mat <-  matrix(c(sd_int**2, rho*sd_int*sd_slope, rho*sd_int*sd_slope, sd_slope**2), ncol = 2, nrow = 2)
          true_sig<-corr_mat



          nms<-paste(N,n,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")

          nmsfile<-paste0("results/ses",nms,".txt")

          dd<-try(read.table(nmsfile,header=FALSE,sep="\t"),silent = TRUE)
          if (class(dd)!="try-error"){
            zz=zz+1
            if (nrow(dd)>1000) dd<-dd[1:1000,]
            ses<-apply(dd,2,mean,na.rm=TRUE)
            dfi<-data.frame(
              N=N,n=n,vr_int=vr_int,mult=mult,rho=rh,se0=ses[1],
              se1=ses[2],se2=ses[3],se3=ses[4],se4=ses[5],se5=ses[6],nsim=nrow(dd))
            if (zz==1) dfse<-dfi else dfse<-rbind(dfse,dfi)
          }

        }}}}}


###############
#par(mfcol=c(4,16),mar=c(2,2,3,1))
zz=0
for (N in c(20,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.1,0.5)){
      for (mult in sqrt(c(0.01,0.1,0.25,0.5,1))){
        for (rh in c(-0.5,-0.2)){
          


          for (meth in c("kos","kos2","ml","bmln","bmlt", "pml","o", "pml1","o1","pmla2","oa2")){
          #  for (meth in c("kos", "pml","ml")){
           #meth="ml"   

          sd_int<-sqrt(vr_int)
          multiplier_int<-round(mult,7)
          
          sd_slope<-sd_int*multiplier_int

          b0 <- -3.5

          # regression coefficients of fixed effects
          betas <- c(1, 1, 1, 0.35, -0.1)


          rho=rh

          # list of covariance matrices for random effects
          corr_mat <-  matrix(c(sd_int**2, rho*sd_int*sd_slope, rho*sd_int*sd_slope, sd_slope**2), ncol = 2, nrow = 2)
          true_sig<-corr_mat

true_pars<-c(NA,min(eigen(true_sig)$values),max(eigen(true_sig)$values),mean(eigen(true_sig)$values),
             determinant(true_sig)$modulus)

          nms<-paste(N,n,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")

            
          nmsfile<-paste0("results/",meth,nms,".txt")  
           
          
          dd<-try(read.table(nmsfile,header=FALSE,sep="\t"),silent = TRUE)

          if (class(dd)!="try-error"){
            zz=zz+1

            nsim<-nrow(dd)

            ddm<-as.matrix(dd)

            if (meth=="ml"){
              dfmse<-apply(ddm[,23:28,drop=FALSE],2,mean)
              if(is.null(dfmsi)) {dfmsi<-dfmse} else {dfmsi<-rbind(dfmsi,dfmse)}
            }


            is_sing_s<-function(v1,v2,tol){
              if (sum(is.na(c(v1,v2)))==2) NA else {
              if (abs(log(sqrt(v1)))>tol|abs(log(sqrt(v2)))>tol) 1 else 0
              }
            }

            det_s<-sing_s<-loss<-loss2<-min_e<-max_e<-mean_e<-rep(NA,nsim)
            for (i in 1:nsim){
              loss[i]<-get_loss(true_sig,ddm[i,7:9])
              loss2[i]<-get_loss2(true_sig,ddm[i,7:9])
              min_e[i]<- min_eigen(ddm[i,7],ddm[i,8],ddm[i,9])
              max_e[i]<- max_eigen(ddm[i,7],ddm[i,8],ddm[i,9])
              mean_e[i]<- mean_eigen(ddm[i,7],ddm[i,8],ddm[i,9])
              sing_s[i]<-is_sing_s(ddm[i,7],ddm[i,9],lim_sing_s)
              det_s[i]<-get_determinant(ddm[i,7:9])
            }
            get_bound_beta<-function(est,true,se){
             if (any(is.na(est))) 1 else {
              tr<- ifelse(abs(est-true)>10*se,1,0)
             if (sum(tr)>0) 1 else 0
            }
            }
           # if (meth%in%c("kos","pml","po","ml")){
            #  boxplot(cbind(ddm[,c(7:9)],loss),main=meth)
            #  boxplot(cbind(ddm[,c(7:9)]),main=meth)
            #abline(h=c(diag(true_sig),true_sig[1,2]))
            #}
            
            widthci<-apply(ddm[,17:22],2,mean,na.rm=TRUE)
            bound_beta<-rep(NA,nsim)
            for (i in 1:nsim){
              bound_beta[i]<-get_bound_beta(ddm[i,1:6],c(b0,betas),
                  dfse[dfse$N==N&dfse$n==n&dfse$vr_int==vr_int&dfse$mult==mult&dfse$rho==rh,6:11])
            }

            dfi<-data.frame(
              N=N,n=n,vr_int=vr_int,mult=mult,rho=rh,
              method=meth,
                mloss=mean(loss,na.rm=TRUE),
              mloss2=mean(loss2,na.rm=TRUE),
              #mloss=median(loss,na.rm=TRUE),
              na_mloss=mean(is.na(loss)),
            singular=mean(min_e<lim_sing,na.rm=TRUE),
            na_singular=mean(is.na(min_e),na.rm=TRUE),
            mean_eigen=mean(mean_e,na.rm=TRUE),
            na_mean_eigen=mean(is.na(mean_e),na.rm=TRUE),
            covr5=mean(ddm[,16],na.rm=TRUE),
            na_covr5=mean(is.na(ddm[,16]),na.rm=TRUE),

            covr4=mean(ddm[,15],na.rm=TRUE),
            na_covr4=mean(is.na(ddm[,15]),na.rm=TRUE),

            covr3=mean(ddm[,14],na.rm=TRUE),
            na_covr3=mean(is.na(ddm[,14]),na.rm=TRUE),

            covr2=mean(ddm[,13],na.rm=TRUE),
            na_covr2=mean(is.na(ddm[,13]),na.rm=TRUE),

            covr1=mean(ddm[,12],na.rm=TRUE),
            na_covr1=mean(is.na(ddm[,12]),na.rm=TRUE),

            covr0=mean(ddm[,11],na.rm=TRUE),
            na_covr0=mean(is.na(ddm[,11]),na.rm=TRUE),


            nsim=nsim,
            bound_rho=mean(1-abs(ddm[,10])<lim_sing_rho,na.rm=TRUE),
            na_bound_rho=mean(is.na(ddm[,10]),na.rm=TRUE),
            bound_s=mean(sing_s,na.rm=TRUE),
            na_bound_s=mean(is.na(sing_s),na.rm=TRUE),
            bound_t=mean(ifelse(1-abs(ddm[,10])<lim_sing_rho|sing_s==1,1,0),na.rm=TRUE),
            na_bound_t=mean(is.na(ddm[,10])|is.na(sing_s),na.rm=TRUE) ,na_fit=mean(is.na(ddm[,1]),na.rm=TRUE),
            bias_b5=get_bias(ddm[,6],-0.1),rmse_b5=get_rmse(ddm[,6],-0.1),
            bias_b1=get_bias(ddm[,2],1),rmse_b1=get_rmse(ddm[,2],1),
            bias_b2=get_bias(ddm[,3],1),rmse_b2=get_rmse(ddm[,3],1),
            bias_b3=get_bias(ddm[,4],1),rmse_b3=get_rmse(ddm[,4],1),
            bias_b4=get_bias(ddm[,5],0.35),rmse_b4=get_rmse(ddm[,5],0.35),
            bias_b0=get_bias(ddm[,1],-3.5),rmse_b0=get_rmse(ddm[,1],-3.5),
            bound_beta=mean(bound_beta,na.rm=TRUE),
            na_bound_beta=mean(is.na(bound_beta),na.rm=TRUE),
            bound_tot=mean(ifelse(bound_beta==1|1-abs(ddm[,10])<lim_sing_rho|sing_s==1,1,0),na.rm=TRUE),
            na_bound_tot=mean(is.na(bound_beta)|is.na(ddm[,10])|is.na(sing_s),na.rm=TRUE) ,na_fit=mean(is.na(ddm[,1]),na.rm=TRUE),
            bias_vr1=get_bias(ddm[,7],true_sig[1,1]),rmse_vr1=get_rmse(ddm[,7],true_sig[1,1]),
            bias_vr2=get_bias(ddm[,9],true_sig[2,2]),rmse_vr2=get_rmse(ddm[,9],true_sig[2,2]),
            bias_cvr=get_bias(ddm[,8],true_sig[1,2]),rmse_cvr=get_rmse(ddm[,8],true_sig[1,2]),
            bias_rho=get_bias(ddm[,10],rho),rmse_rho=get_rmse(ddm[,10],rho),
            det_est=mean(det_s,na.rm=TRUE),
            det_true=true_pars[5],
            det_est_o=mean(exp(det_s),na.rm=TRUE),
            det_true_o=exp(true_pars[5]),
            min_e=mean(min_e,na.rm=TRUE),
            max_e=mean(max_e,na.rm=TRUE),
            mean_e=mean(mean_e,na.rm=TRUE),
            true_min_e=true_pars[2],
            true_max_e=true_pars[3],
            true_mean_e=true_pars[4],widths_0=widthci[1]
            ,widths_1=widthci[2],widths_2=widthci[3]
            ,widths_3=widthci[4],widths_4=widthci[5]
            ,widths_5=widthci[6]
            )

            if (zz==1) df<-dfi else df<-rbind(df,dfi)

          }


          }
        }
    }}}
}

 
display_taus<-FALSE

if (display_taus){
##see taus - too many for all sims!
#par(mfrow=c(12,10),mar=c(3,2,3,1))
par(mfrow=c(8,10),mar=c(3,2,3,1))
zz=0
for (N in c(20,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.1,0.5)){
      for (mult in sqrt(c(0.01,0.1,0.25,0.5,1))){
        for (rh in c(-0.5,-0.2)){
          
          
          
             
            
            sd_int<-sqrt(vr_int)
            multiplier_int<-round(mult,7)
            
            sd_slope<-sd_int*multiplier_int
            
            b0 <- -3.5
            
            # regression coefficients of fixed effects
            betas <- c(1, 1, 1, 0.35, -0.1)
            
            
            rho=rh
            
            # list of covariance matrices for random effects
            corr_mat <-  matrix(c(sd_int**2, rho*sd_int*sd_slope, rho*sd_int*sd_slope, sd_slope**2), ncol = 2, nrow = 2)
            true_sig<-corr_mat
            
            
            
            nms<-paste(N,n,diag(corr_mat)[1],diag(corr_mat)[2],cov2cor(corr_mat)[1,2],sep="_")
            
                nmsfile<-paste0("results/","taus",nms,".txt") #using alpha=.5, no truncation for POR
            
             
            
            
            dd<-try(read.table(nmsfile,header=FALSE,sep="\t"),silent = TRUE)
            
            if (class(dd)!="try-error"){
              zz=zz+1
              
            
               
              ddm<-as.matrix(dd)
              colnames(ddm)<-c("PML-A1","PML-A2","PML-A2(1)")
            boxplot(ddm,main=nms,ylim=c(0,1))
          
             #   dfi<-data.frame(
            #    N=N,n=n,vr_int=vr_int,mult=mult,rho=rh
                
             #         )
              
            #  if (zz==1) df<-dfi else df<-rbind(df,dfi)
              
            }
            
            
          
        }
      }}}}
}



split(df$na_fit,df$method) #nas for -1 meths are since bglm2 is not working (it is not giving error, but it gives weird pis which cause error when inverting to get the hat matrix)!
#using logistf seems to remove this issue

#txt1<-expression(sigma[1]^2)
#txt2<-expression(sigma[2]^2)
txt1<-"sigma[0]^2"
txt2<-"sigma[1]^2"

df$Nn<-paste0("N=",df$N,",n=",df$n)
df$vars<-paste0(txt1,"==",df$vr_int,"*','*~",txt2,"==",df$vr_int*df$mult**2)


df$Nn2<-as.numeric(as.factor(df$Nn))


df$r2<-paste0("rho==",df$rho)




##num sims that are done
ggplot(df[df$method=="ml",])+
  geom_point(aes(x=Nn2,y=nsim),size=2)+
  geom_line(aes(x=Nn2,y=nsim),size=.7)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  #scale_x_continuous(name="", breaks=c(1,2,3,4), labels=levels(factor(dft$Nn)))+
  ylab("n sim")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))#+
  #geom_hline(yintercept=1000)





###reduced df (only for PRIAL)

zz=0
for (N in c(20,50)){
  for (n in c(10,20)){
    for (vr_int in c(0.1,0.5)){
      for (mult in sqrt(c(0.01,0.1,0.25,0.5,1))){
        for (rh in c(-0.5,-0.2)){
          

          for (meth in c("kos","kos2","bmln","bmlt", "pml","o", "pml1","o1","pmla2","oa2")){
            zz=zz+1
             
            dfi<-df[df$N==N&df$n==n&df$vr_int==vr_int&df$mult==mult&df$rho==rh&df$method%in%c(meth,"ml"),c(1:8,which(names(df)%in%c("Nn","vars","Nn2","r2")))]
            if (nrow(dfi)==2){
              prial<-(dfi$mloss[dfi$method=="ml"]-dfi$mloss[dfi$method==meth])/dfi$mloss[dfi$method=="ml"]*100
              prial2<-(dfi$mloss2[dfi$method=="ml"]-dfi$mloss2[dfi$method==meth])/dfi$mloss2[dfi$method=="ml"]*100
              dfii<-dfi[dfi$method==meth,]
              dfii$prial<-prial
              dfii$prial2<-prial2
              if (zz==1) df2<-dfii else df2<-rbind(df2,dfii)
            }
          }
        }
      }}}}


 
################final figs as for the paper: only show useful methods

df<-df[df$method%in%c("kos","ml","bmln","bmlt","pml","o","pml1","o1","pmla2","oa2"),]

#rename methods as in paper

##rename methods as in paper!
df$method<-as.character(df$method)
df$method[df$method=="kos"]<-"SPM"

df$method[df$method=="ml"]<-"ML"

df$method[df$method=="bmln"]<-"BMN"
df$method[df$method=="bmlt"]<-"BMT"



 df$method[df$method=="pml"]<-"MPL-A1"

df$method[df$method=="o"]<-"ORACLE-A1"

df$method[df$method=="pml1"]<-"MPL-A2(1)"

df$method[df$method=="o1"]<-"ORACLE-A2(1)"

df$method[df$method=="pmla2"]<-"MPL-A2"

df$method[df$method=="oa2"]<-"ORACLE-A2"

 


df$method<-factor(df$method,levels=c("ML","BMN","BMT","MPL-A1","MPL-A2","MPL-A2(1)","ORACLE-A1","ORACLE-A2","ORACLE-A2(1)","SPM"))
 
df2<-df2[df2$method%in%c("kos","bmln","bmlt","pml","o","pml1","o1","pmla2","oa2"),]
df2$method<-as.character(df2$method)
df2$method[df2$method=="kos"]<-"SPM"


df2$method[df2$method=="bmln"]<-"BMN"
df2$method[df2$method=="bmlt"]<-"BMT"


 df2$method[df2$method=="pml"]<-"MPL-A1"
 df2$method[df2$method=="o"]<-"ORACLE-A1"
 df2$method[df2$method=="o1"]<-"ORACLE-A2(1)"
 
 df2$method[df2$method=="pml1"]<-"MPL-A2(1)"
 
 df2$method[df2$method=="oa2"]<-"ORACLE-A2"
 
 df2$method[df2$method=="pmla2"]<-"MPL-A2"
 

df2$method<-factor(df2$method,levels=c("BMN","BMT","MPL-A1","MPL-A2","MPL-A2(1)","ORACLE-A1","ORACLE-A2","ORACLE-A2(1)","SPM"))
 
dfpallete <- c("black","#332288", "#117733", "#44AA99","#66BBCC", "#88CCEE",
                "#DDCC77","#D59977", "#CC6677", "#AA4499")

dfpallete2 <- c("#332288", "#117733", "#44AA99", "#66BBCC", "#88CCEE",
               "#DDCC77","#D59977", "#CC6677", "#AA4499")


 
df$group<-ifelse(df$N==20,"1","2")
df2$group<-ifelse(df2$N==20,"1","2")

 
nms <- c(
  "atop(N==20, n[i]==S)",
  "atop(N==20, n[i]==M)",
  "atop(N==50, n[i]==S)",
  "atop(N==50, n[i]==M)"
)


df$vars<-factor(df$vars)

df2$vars<-factor(df2$vars)

#supplement: boundary est all criteria, PRIAL, all methods
alphap=0.75

#######main document
which.rho=-0.2 #determines which rho we will show in the main doc
show.L2=FALSE #display figs with L2? (Set to FALSE to address Georg's comment)
offset.x=FALSE #offset the methods in direction of x? avoids overlap but seems weird!
ylab.prial="PRIAL (%)"
#scenarios to show in main file
#theta=0.25,2.25
mult2<-c(0.01,0.1,0.25)
vrint<-c(0.1,0.5)

mult2<-c(0.01,0.1,0.25,0.5,1)
vrint<-c(0.1)


plot_ind<-FALSE #this plots individual plots used in supp in rows
#methods to show (omit PML, ORACLE and say in the text they are the same (show them in supplement))
#methods_omit<-c("MPL-A1","ORACLE-A1","MPL-A2","ORACLE-A2")
methods_omit<-c("MPL-A2","ORACLE-A2")

df_main<-df[round(df$mult**2,3)%in%mult2&df$vr_int%in%vrint&!df$method%in%methods_omit&df$rho==which.rho,]
df2_main<-df2[round(df2$mult**2,3)%in%mult2&df2$vr_int%in%vrint&!df2$method%in%methods_omit&df2$rho==which.rho,]

#df_main<-df[df$mult%in%mult&!df$method%in%methods_omit,]
#df2_main<-df2[df2$mult%in%mult&!df2$method%in%methods_omit,]


dfpallete_main <- c("black","#332288", "#117733","#44AA99",  "#88CCEE",
                    "#DDCC77", "#CC6677", "#AA4499")

dfpallete2_main <- c("#332288", "#117733", "#44AA99", "#88CCEE",
                     "#DDCC77",  "#CC6677", "#AA4499")


if(offset.x){
offset<-seq(from=-0.2,to=0.2,length.out=6)
df_main$x_off<-df_main$Nn2+rep(offset,nrow(df_main)/6)
offset2<-seq(from=-0.2,to=0.2,length.out=5)
df2_main$x_off<-df2_main$Nn2+rep(offset2,nrow(df_main)/6)
} else {
    df_main$x_off<-df_main$Nn2 
   df2_main$x_off<-df2_main$Nn2 
  
  
}

###some desc for bound(s) to be mentioned in the text:

#number of scenarios with any boundary estimate (by type)

sapply(split(df$bound_tot,df$method),function(x) sum(x>0)) 
sapply(split(df$bound_beta,df$method),function(x) sum(x>0)) 
sapply(split(df$bound_t,df$method),function(x) sum(x>0)) 
sapply(split(df$bound_s,df$method),function(x) sum(x>0)) 
sapply(split(df$bound_r,df$method),function(x) sum(x>0)) 




#the largest proportion of boundary estimates


sapply(split(df$bound_tot,df$method),function(x) max(x)) 
sapply(split(df$bound_beta,df$method),function(x) max(x)) 
sapply(split(df$bound_t,df$method),function(x) max(x)) 
sapply(split(df$bound_s,df$method),function(x) max(x)) 
sapply(split(df$bound_r,df$method),function(x) max(x)) 


#total proportion of datasets with boundary estimates

round(sapply(split(df$bound_tot,df$method),function(x) sum(x*1000)/(80*1000)*100) ,3)
round(sapply(split(df$bound_beta,df$method),function(x) sum(x*1000)/(80*1000)*100) ,3)
      round(sapply(split(df$bound_t,df$method),function(x) sum(x*1000)/(80*1000)*100),3) 
#round(sapply(split(df$bound_s,df$method),function(x) sum(x*1000)/(80*1000)) ,2)
#round(sapply(split(df$bound_r,df$method),function(x) sum(x*1000)/(80*1000)),2) 


#total number of datasets with boundary estimates

sapply(split(df$bound_tot,df$method),function(x) sum(x*1000)) #all
sapply(split(df$bound_beta,df$method),function(x) sum(x*1000)) #betas
sapply(split(df$bound_t,df$method),function(x) sum(x*1000)) #sigma
#sapply(split(df$bound_s,df$method),function(x) sum(x*1000)) #sds
#sapply(split(df$bound_r,df$method),function(x) sum(x*1000)) #rho


###bound (overall only)

p_bound<-ggplot(df_main,aes(group=interaction(method, group)))+
  geom_point(aes(x=x_off,y=bound_tot*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=x_off,y=bound_tot*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+#, linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  geom_vline(xintercept = c(1.5,2.5, 3.5),  color = "grey90") +
  ylab("boundary estimate (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete_main)+
  theme(legend.title=element_blank())+
  theme(
    panel.grid.major.x = element_blank(),  # remove vertical major gridlines
    panel.grid.minor.x = element_blank()   # remove vertical minor gridlines
  )






###PRIAL

pp_prial1<-ggplot(df2_main,aes(group=interaction(method, group)))+
  geom_point(aes(x=x_off,y=prial,colour =  method),size=2,alpha=alphap )+
  geom_line(aes(x=x_off,y=prial,colour =  method),size=0.7,alpha=alphap,linetype = "dashed")+#, linetype = "dashed")+
  facet_grid(r2~vars,scales="free_y",labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  geom_vline(xintercept = c(1.5,2.5, 3.5),  color = "grey90") +
  ylab(ylab.prial)+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2_main)+
  theme(legend.title=element_blank())+
  theme(
    panel.grid.major.x = element_blank(),  # remove vertical major gridlines
    panel.grid.minor.x = element_blank()   # remove vertical minor gridlines
  )+
  ylim(c(-300,100))

pp_prial2<-ggplot(df2_main,aes(group=interaction(method, group)))+
  geom_point(aes(x=x_off,y=prial2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=x_off,y=prial2,colour =  method),size=0.7,alpha=alphap,linetype = "dashed")+#, linetype = "dashed")+
  facet_grid(r2~vars,scales="free_y",labeller = label_parsed)+
  geom_vline(xintercept = c(1.5,2.5, 3.5),  color = "grey90") +
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("PRIAL-L2 (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2_main)+
  theme(legend.title=element_blank())+
  theme(
    panel.grid.major.x = element_blank(),  # remove vertical major gridlines
    panel.grid.minor.x = element_blank()   # remove vertical minor gridlines
  )+
  ylim(c(-300,100))


#Fig1
if (show.L2){
pdf("figs/Fig1.pdf",height=15/5*3,width=15)
grid.arrange(p_bound,pp_prial1,pp_prial2)
dev.off()
} else {
  pdf("figs/Fig1.pdf",height=15/5*2,width=15/6*5)
  grid.arrange(p_bound,pp_prial1)
  dev.off()
}

#Fig2, beta5 properties

##bias,rmse, cover for beta5

p1b<-ggplot(df_main,aes(group=interaction(method, group)))+
  geom_point(aes(x=x_off,y=bias_b5,colour =  method),size=2,alpha=alphap )+
  geom_line(aes(x=x_off,y=bias_b5,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+#, linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  geom_vline(xintercept = c(1.5,2.5, 3.5),  color = "grey90") +
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete_main)+
  theme(legend.title=element_blank())+
  theme(
    panel.grid.major.x = element_blank(),  # remove vertical major gridlines
    panel.grid.minor.x = element_blank()   # remove vertical minor gridlines
  )

p1b2<-ggplot(df_main,aes(group=interaction(method, group)))+
  geom_point(aes(x=x_off,y=log(abs(bias_b5)),colour =  method),size=2,alpha=alphap )+
  geom_line(aes(x=x_off,y=log(abs(bias_b5)),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+#, linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  geom_vline(xintercept = c(1.5,2.5, 3.5),  color = "grey90") +
  ylab("log(|bias|)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete_main)+
  theme(legend.title=element_blank())+
  theme(
    panel.grid.major.x = element_blank(),  # remove vertical major gridlines
    panel.grid.minor.x = element_blank()   # remove vertical minor gridlines
  )


p2b<-ggplot(df_main,aes(group=interaction(method, group)))+
  geom_point(aes(x=x_off,y=rmse_b5,colour =  method),size=2,alpha=alphap )+
  geom_line(aes(x=x_off,y=rmse_b5,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+#, linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  geom_vline(xintercept = c(1.5,2.5, 3.5),  color = "grey90") +
  ylab("RMSE")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  scale_colour_manual(values=dfpallete_main)+
  theme(legend.title=element_blank())+
  theme(
    panel.grid.major.x = element_blank(),  # remove vertical major gridlines
    panel.grid.minor.x = element_blank()   # remove vertical minor gridlines
  )



p3b<-ggplot(df_main,aes(group=interaction(method, group)))+
  geom_point(aes(x=x_off,y=covr5*100,colour =  method),size=2,alpha=alphap )+
  geom_line(aes(x=x_off,y=covr5*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+#, linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  geom_vline(xintercept = c(1.5,2.5, 3.5),  color = "grey90") +
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete_main)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 0.75, xmax = 4.25, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)+
  theme(
    panel.grid.major.x = element_blank(),  # remove vertical major gridlines
    panel.grid.minor.x = element_blank()   # remove vertical minor gridlines
  )

 
pdf("figs/Fig2.pdf",height=15/5*3,width=15/6*5)
grid.arrange(p1b,p2b,p3b,nrow=3)
dev.off()








###supp, all scenarios, further results, no jittering

df$r2<-factor(df$r2,levels=c("rho==-0.5","rho==-0.2"))

###convergence



p1<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_rho*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_rho*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+ #this should be active when res are done!
  ylab("boundary correlation (%)")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_s*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_s*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("boundary standard deviations (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p0<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_t*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_t*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("boundary sigma (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p01<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_beta*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_beta*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("boundary beta (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p02<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_tot*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_tot*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("boundary estimate (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bound.pdf",height=25,width=22)
grid.arrange(p02,p01,p0,p1,p2,nrow=5)
dev.off()

if (plot_ind){
pdf("figs/supp_bound_all.pdf",height=20,width=24)
p02
dev.off()


pdf("figs/supp_bound_beta.pdf",height=20,width=24)
p01
dev.off()

pdf("figs/supp_bound_sigma.pdf",height=20,width=24)
p0
dev.off()

pdf("figs/supp_bound_rho.pdf",height=20,width=24)
p1
dev.off()

pdf("figs/supp_bound_sd.pdf",height=20,width=22)
p2
dev.off()
}

##estimation of D

df2$r2<-factor(df2$r2,levels=c("rho==-0.5","rho==-0.2"))


pp<-ggplot(df2,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=prial,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=prial,colour =  method),size=0.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,scales="free_y",labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(ylab.prial)+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2)+
  theme(legend.title=element_blank())+
  ylim(c(-300,100))

pp2<-ggplot(df2,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=prial2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=prial2,colour =  method),size=0.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,scales="free_y",labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("PRIAL-L2 (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2)+
  theme(legend.title=element_blank())+
  ylim(c(-300,100))







ploss<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y= log(mloss,10) ,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y= log(mloss,10) ,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10](loss)))+
  #ylab("loss")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+

  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

ploss2<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y= log(mloss2,10) ,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y= log(mloss2,10) ,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10](loss2)))+
  #ylab("loss")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

if (show.L2){
pdf("figs/supp_loss.pdf",height=25,width=22)
grid.arrange(ploss,pp,ploss2,pp2,nrow=4)
dev.off()
} else{
  pdf("figs/supp_loss.pdf",height=10,width=22)
  grid.arrange(ploss,pp,nrow=2)
  dev.off()
}

if (plot_ind){
pdf("figs/supp_loss_loss.pdf",height=20,width=24)
ploss
dev.off()

pdf("figs/supp_loss_prial.pdf",height=20,width=24)
pp
dev.off()
}

##bias,rmse, cover, width for beta5

p1b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b5,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b5,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p1b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log10(abs(bias_b5)),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log10(abs(bias_b5)),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10]*abs(bias)))+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  #geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b5,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b5,base=10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+

  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr5*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr5*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)

p4b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_5,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_5,base=10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
 # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
  #annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
   #        alpha = .2)


pdf("figs/supp_beta5.pdf",height=20,width=22)
grid.arrange(p1b,p2b,p3b,p4b,nrow=4)
dev.off()

if (plot_ind){
pdf("figs/supp_beta5_1.pdf",height=20,width=24)
p1b
dev.off()

pdf("figs/supp_beta5_2.pdf",height=20,width=24)
p2b
dev.off()

pdf("figs/supp_beta5_3.pdf",height=20,width=24)
p3b
dev.off()

pdf("figs/supp_beta5_4.pdf",height=20,width=24)
p4b
dev.off()
}

##bias, mse,...for other betas

p1b1<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b1,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b1,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p1b1<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log10(abs(bias_b1)),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log10(abs(bias_b1)),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10]*abs(bias)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  #geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())


p2b1<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b1,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b1,10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+

  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b1<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr1*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr1*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)

p4b1<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_1,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_1,base=10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)


pdf("figs/supp_beta1.pdf",height=20,width=22)
grid.arrange(p1b1,p2b1,p3b1,p4b1,nrow=4)
dev.off()

if (plot_ind){
pdf("figs/supp_beta1_1.pdf",height=20,width=24)
p1b1
dev.off()

pdf("figs/supp_beta1_2.pdf",height=20,width=24)
p2b1
dev.off()

pdf("figs/supp_beta1_3.pdf",height=20,width=24)
p3b1
dev.off()

pdf("figs/supp_beta1_4.pdf",height=20,width=24)
p4b1
dev.off()
}




p1b2<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b2,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p1b2<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log10(abs(bias_b2)),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log10(abs(bias_b2)),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10]*abs(bias)))+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  #geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())


p2b2<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b2,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b2,10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+

  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b2<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr2*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr2*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)

p4b2<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_2,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_2,base=10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)


pdf("figs/supp_beta2.pdf",height=20,width=22)
grid.arrange(p1b2,p2b2,p3b2,p4b2,nrow=4)
dev.off()

if (plot_ind){
pdf("figs/supp_beta2_1.pdf",height=20,width=24)
p1b2
dev.off()

pdf("figs/supp_beta2_2.pdf",height=20,width=24)
p2b2
dev.off()

pdf("figs/supp_beta2_3.pdf",height=20,width=24)
p3b2
dev.off()

pdf("figs/supp_beta2_4.pdf",height=20,width=24)
p4b2
dev.off()
}





p1b3<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b3,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b3,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())




p1b3<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log10(abs(bias_b3)),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log10(abs(bias_b3)),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10]*abs(bias)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  #geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())



p2b3<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b3,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b3,10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+

  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b3<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr3*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr3*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)

p4b3<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_3,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_3,base=10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)


pdf("figs/supp_beta3.pdf",height=20,width=22)
grid.arrange(p1b3,p2b3,p3b3,p4b3,nrow=4)
dev.off()

if (plot_ind){

pdf("figs/supp_beta3_1.pdf",height=20,width=24)
p1b3
dev.off()

pdf("figs/supp_beta3_2.pdf",height=20,width=24)
p2b3
dev.off()

pdf("figs/supp_beta3_3.pdf",height=20,width=24)
p3b3
dev.off()

pdf("figs/supp_beta3_4.pdf",height=20,width=24)
p4b3
dev.off()


}


p1b4<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b4,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b4,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p1b4<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log10(abs(bias_b4)),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log10(abs(bias_b4)),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10]*abs(bias)))+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  #geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())





p2b4<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b4,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b4,10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+

  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b4<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr4*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr4*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)


p4b4<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_4,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_4,base=10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)



pdf("figs/supp_beta4.pdf",height=20,width=22)
grid.arrange(p1b4,p2b4,p3b4,p4b4,nrow=4)
dev.off()


if (plot_ind){
pdf("figs/supp_beta4_1.pdf",height=20,width=24)
p1b4
dev.off()

pdf("figs/supp_beta4_2.pdf",height=20,width=24)
p2b4
dev.off()

pdf("figs/supp_beta4_3.pdf",height=20,width=24)
p3b4
dev.off()

pdf("figs/supp_beta4_4.pdf",height=20,width=24)
p4b4
dev.off()
}




p1b0<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b0,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b0,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p1b0<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log10(abs(bias_b0)),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log10(abs(bias_b0)),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10]*abs(bias)))+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  #geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())





p2b0<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b0,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b0,10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b0<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr0*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr0*100,colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)


p4b0<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_0,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_0,base=10),colour =  method),size=.7,alpha=alphap,linetype = "dashed")+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=parse(text = nms),limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)



pdf("figs/supp_beta0.pdf",height=20,width=22)
grid.arrange(p1b0,p2b0,p3b0,p4b0,nrow=4)
dev.off()


if (plot_ind){
  pdf("figs/supp_beta0_1.pdf",height=20,width=24)
  p1b0
  dev.off()
  
  pdf("figs/supp_beta0_2.pdf",height=20,width=24)
  p2b0
  dev.off()
  
  pdf("figs/supp_beta0_3.pdf",height=20,width=24)
  p3b0
  dev.off()
  
  pdf("figs/supp_beta0_4.pdf",height=20,width=24)
  p4b0
  dev.off()
}


####other plots, not needed

plot_others<-FALSE

if (plot_others){




#####the same, but only for PML, PML-1 and SPM


###convergence

dfs<-df[df$method%in%c("PML","PML-1","SPM"),]
df2s<-df2[df2$method%in%c("PML","PML-1","SPM"),]

dfpallete<-dfpallete[c(4,5,8)]
dfpallete2<-dfpallete2[c(4,5,8)-1]

p1<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_rho*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_rho*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+ #this should be active when res are done!
  ylab("boundary correlation (%)")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_s*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_s*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("boundary standard deviations (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p0<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_t*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_t*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("boundary sigma (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p01<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_beta*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_beta*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("boundary beta (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p02<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_tot*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_tot*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("boundary estimate (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())


##estimation of D

pp<-ggplot(df2s,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=prial,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=prial,colour =  method),size=0.7,alpha=alphap)+
  facet_grid(r2~vars,scales="free_y",labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("PRIAL-L1 (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2)+
  theme(legend.title=element_blank())#+
  #ylim(c(-300,100))

pp2<-ggplot(df2s,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=prial2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=prial2,colour =  method),size=0.7,alpha=alphap)+
  facet_grid(r2~vars,scales="free_y",labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("PRIAL-L2 (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2)+
  theme(legend.title=element_blank())#+
  #ylim(c(-300,100))


 

pdf("figs/supp_bin_boundary_2s.pdf",height=20,width=20)
grid.arrange(p02,p01,p0,p1,p2,nrow=5)
dev.off()


ploss<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y= log(mloss,10) ,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y= log(mloss,10) ,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](loss)))+
  #ylab("loss")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

ploss2<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y= log(mloss2,10) ,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y= log(mloss2,10) ,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](loss2)))+
  #ylab("loss")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())
 

pdf("figs/supp_bin_loss_2s.pdf",height=20/5*2,width=20)
grid.arrange(ploss,ploss2,nrow=2)
dev.off()

##bias,rmse  for vr1

p1b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_vr1,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_vr1,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_vr1,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_vr1,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_vr1s.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()

dfpalletem<-c("black", "#44AA99","#88CCEE",
               "#AA4499")

p1b<-ggplot(df[df$method%in%c("ML","PML-1","PML","SPM")&!(df$vr_int==0.5&df$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_vr1,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_vr1,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpalletem)+
  theme(legend.title=element_blank())

p2b<-ggplot(df[df$method%in%c("ML","PML-1","PML","SPM")&!(df$vr_int==0.5&df$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_vr1,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_vr1,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpalletem)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_vr1s_m.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()
##bias,rmse  for vr2

p1b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_vr2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_vr2,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_vr2,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_vr2,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_vr2s.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()


p1b<-ggplot(dfs[!(dfs$vr_int==0.5&dfs$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_vr2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_vr2,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(dfs[!(dfs$vr_int==0.5&dfs$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_vr2,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_vr2,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_vr2s_2.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()


p1b<-ggplot(df[df$method%in%c("ML","PML-1","PML","SPM")&!(df$vr_int==0.5&df$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_vr2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_vr2,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpalletem)+
  theme(legend.title=element_blank())

p2b<-ggplot(df[df$method%in%c("ML","PML-1","PML","SPM")&!(df$vr_int==0.5&df$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_vr2,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_vr2,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpalletem)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_vr2s_m.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()
##bias,rmse  for cvr

p1b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_cvr,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_cvr,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_cvr,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_cvr,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_cvrs.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()




p1b<-ggplot(dfs[!(dfs$vr_int==0.5&dfs$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_cvr,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_cvr,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(dfs[!(dfs$vr_int==0.5&dfs$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_cvr,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_cvr,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_cvrs_2.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()

p1b<-ggplot(df[df$method%in%c("ML","PML-1","PML","SPM")&!(df$vr_int==0.5&df$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_cvr,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_cvr,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpalletem)+
  theme(legend.title=element_blank())

p2b<-ggplot(df[df$method%in%c("ML","PML-1","PML","SPM")&!(df$vr_int==0.5&df$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_cvr,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_cvr,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpalletem)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_cvrs_m.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()

##bias,rmse  for rho

p1b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_rho,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_rho,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_rho,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_rho,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_rhos.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()



p1b<-ggplot(df[df$method%in%c("ML","PML-1","PML","SPM")&!(df$vr_int==0.5&df$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_rho,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_rho,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpalletem)+
  theme(legend.title=element_blank())

p2b<-ggplot(df[df$method%in%c("ML","PML-1","PML","SPM")&!(df$vr_int==0.5&df$mult==10),],aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_rho,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_rho,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpalletem)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_rhos_m.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()
##bias,rmse, cover for beta5

p1b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b5,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b5,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b5,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b5,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr5*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr5*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)

p4b<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_5,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_5,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)

pdf("figs/supp_bin_beta5_2s.pdf",height=20/5*4,width=20)
grid.arrange(p1b,p2b,p3b,p4b,nrow=4)
dev.off()


#main document, boundary

p0<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bound_tot*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bound_tot*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("boundary estimate (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/main_bin_boundary_2s.pdf",height=20/5,width=20)
p0
dev.off()

#tiff("figs/Fig3.tif",height=6,width=12,units="in",res=1200,compression = "lzw")
#p0
#dev.off()

#main document, prial

pdf("figs/main_bin_prial_2s.pdf",height=20/5*2,width=20)
grid.arrange(pp,pp2,nrow=2)
dev.off()


#tiff("figs/Fig4.tif",height=6,width=12,units="in",res=1200,compression = "lzw")
#pp
#dev.off()

##bias, mse,...for other betas

p1b1<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b1,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b1,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b1<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b1,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b1,10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b1<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr1*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr1*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)

p4b1<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_1,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_1,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)


pdf("figs/supp_bin_beta1_2s.pdf",height=20/5*4,width=20)
grid.arrange(p1b1,p2b1,p3b1,p4b1,nrow=4)
dev.off()

p1b2<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b2,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b2<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b2,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b2,10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b2<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr2*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr2*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)

p4b2<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_2,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_2,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)


pdf("figs/supp_bin_beta2_2s.pdf",height=20/5*4,width=20)
grid.arrange(p1b2,p2b2,p3b2,p4b2,nrow=4)
dev.off()



p1b3<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b3,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b3,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b3<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b3,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b3,10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b3<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr3*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr3*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)

p4b3<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_3,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_3,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)


pdf("figs/supp_bin_beta3_2s.pdf",height=20/5*4,width=20)
grid.arrange(p1b3,p2b3,p3b3,p4b3,nrow=4)
dev.off()



p1b4<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_b4,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_b4,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b4<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_b4,10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_b4,10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p3b4<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=covr4*100,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=covr4*100,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("coverage (%)")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+
  annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
           alpha = .2)


p4b4<-ggplot(dfs,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(widths_4,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(widths_4,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("log mean CI width")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  # geom_hline(yintercept=95)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())#+
#annotate("rect", xmin = 1, xmax = 4, ymin =95-sqrt(0.05*0.95/1000)*2*100 , ymax =95+sqrt(0.05*0.95/1000)*2*100 ,
#        alpha = .2)



pdf("figs/supp_bin_beta4_2s.pdf",height=20/5*4,width=20)
grid.arrange(p1b4,p2b4,p3b4,p4b4,nrow=4)
dev.off()




####used to understand what we have, not shown in paper


#estimation of the determinant, eigenvalues

d1log<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=det_est-det_true,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=det_est-det_true,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+ #this should be active when res are done!
  ylab("mean log determinant-true log determinant")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())


d1<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=det_est_o-det_true_o,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=det_est_o-det_true_o,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+ #this should be active when res are done!
  ylab("mean determinant-true determinant")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+ylim(c(-5,5))



emin<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=min_e-true_min_e,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=min_e-true_min_e,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+ #this should be active when res are done!
  ylab("mean min eigen-true min eigen")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+ylim(c(-0.5,0.5))

emax<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=max_e-true_max_e,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=max_e-true_max_e,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+ #this should be active when res are done!
  ylab("mean max eigen-true max eigen")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+ylim(c(-0.5,1))

emean<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=mean_e-true_mean_e,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=mean_e-true_mean_e,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+ #this should be active when res are done!
  ylab("mean mean eigen-true mean eigen")+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())+ylim(c(-0.5,0.5))


grid.arrange(d1log,emin,emax,emean,nrow=4)


##bias,rmse  for vr1

p1b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_vr1,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_vr1,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_vr1,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_vr1,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_vr1.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()

##bias,rmse  for vr2

p1b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_vr2,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_vr2,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_vr2,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_vr2,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_vr2.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()

##bias,rmse  for cvr

p1b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_cvr,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_cvr,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_cvr,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_cvr,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_cvr.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()


##bias,rmse  for rho

p1b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=bias_rho,colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=bias_rho,colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab("bias")+
  theme(axis.text.x = element_text(  vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

p2b<-ggplot(df,aes(group=interaction(method, group)))+
  geom_point(aes(x=Nn2,y=log(rmse_rho,base=10),colour =  method),size=2,alpha=alphap)+
  geom_line(aes(x=Nn2,y=log(rmse_rho,base=10),colour =  method),size=.7,alpha=alphap)+
  facet_grid(r2~vars,labeller = label_parsed,scales="free_y")+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=nms,limits = c(0.75,4.25))+
  ylab(expression(log[10](RMSE)))+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))+
  
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())

pdf("figs/supp_bin_rho.pdf",height=20/5*2,width=20)
grid.arrange(p1b,p2b,nrow=2)
dev.off()
}