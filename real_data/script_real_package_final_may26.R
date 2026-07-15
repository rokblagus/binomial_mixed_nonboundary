library(AUGglmmTMB)
library(glmmTMB)
library(MASS)
library(optimx)
library(lme4)
library(numDeriv)
library(blme)
library(mvtnorm)
library(xtable)
library(car)
library(ggplot2)
library(tidyr)
library(dplyr)
library(brglm2)
library(reformulas)

#dd<-read.table("dataBirds.txt",header=T,dec=",",sep="\t",fill=TRUE)
data(birds)
dd<-birds

#2 REs


fit_ml<-glmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=dd,family=binomial(link = "logit"))

#summarize the data

mod<-fit_ml

#sd for 2nd RE
vr<-mod$sdr$cov.fixed[9,9]
est<-mod$sdr$par.fixed[9]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 1st RE

vr<-mod$sdr$cov.fixed[6:7,6:7]
est<-mod$sdr$par.fixed[6:7]
names(est)<-colnames(vr)<-rownames(vr)<-paste0("theta",1:2)

s3.1<-deltaMethod(est,"exp(theta1)",vr)
s3.2<-deltaMethod(est,"exp(theta2)",vr)

#off diags

vr<-mod$sdr$cov.fixed[8,8]
est<-mod$sdr$par.fixed[8]
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


##lmer
fit_mle<-glmer(parasites~migration+food+(migration|phylogenetic)+(1|species),data=dd,family=binomial(link = "logit"))


mod<-fit_mle

#sd for 2nd RE (lme internaly shifts the order!)

h <- mod@optinfo$derivs$Hessian
h <- solve(h)
v <- forceSymmetric(h + t(h))

mod@optinfo$val

vr<-v[1,1]
est<-mod@optinfo$val[1]
names(est)<-"theta"

s1<-deltaMethod(est,"theta",vr)

#sd for 1st RE

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




res<-cbind(resi_glmer,resi_tmb)
colnames(res)<-c("lme4","glmmTMB")
rownames(res)<-c("Intercept","migration","food[granivourous]","food[all]","food[insectivourous]","species (Intercept)"
                 ,"plyo (Intercept)","plyo (migration)","rho")

resx<-xtable(res)
resx


####penalized methods

fit_blm.t<-bglmer(parasites~migration+food+(migration|phylogenetic)+(1|species),data=dd,family=binomial(link = "logit"),fixef.prior = "t")
fit_blm.n<-bglmer(parasites~migration+food+(migration|phylogenetic)+(1|species),data=dd,family=binomial(link = "logit"),fixef.prior = "normal")



#MPL-A1
fit_rok<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=dd,
                    link="logit",
                    penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=FALSE),
                    control=AUGglmmTMBControl(fit_pGLM = FALSE,maxiter = 50,tol=1e-5,save_coef=TRUE))




#summarize the data

mod<-fit_rok$fit
#sd for 2nd RE
vr<-mod$sdr$cov.fixed[6+3,6+3]
est<-mod$sdr$par.fixed[6+3]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 2nd RE

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

res<-cbind(resi_rok,resi_bglmern,resi_bglmert)
colnames(res)<-c("MPL","BMN","BMT")
rownames(res)<-c("Intercept","migration","food[granivorous]","food[all]","food[insectivorous]","species (Intercept)"
                 ,"plyo (Intercept)","plyo (migration)","rho")

resx<-xtable(res)
resx


####plots!


####plot coef paths for all methods


m2<-fit_rok$coefs$coefs_glmm

nm<-max(nrow(m2))







# Convert matrices to data frames

df2 <- as.data.frame(m2)

names(df2)<-c(paste0("beta",0:4),paste0("theta2",1:3),"theta1")
# Add row index
df2$row <- 1:nrow(df2)

# Add source label
df2$matrix <- "MPL"


# rescale thetas to sds(rho)

df2$theta21<-exp(df2$theta21)
df2$theta22<-exp(df2$theta22)
df2$theta23<-df2$theta23/sqrt(1+df2$theta23**2) 


df2$theta1<-exp(df2$theta1)

#rename cols to match table

colnames(df2)<-c("Intercept","migration","food[granivorous]","food[omnivorous]","food[insectivorous]",
                 "plyo(Intercept)","plyo(migration)","rho",
                 "species(Intercept)","row","matrix")



# Combine
df <- bind_rows( df2 )

# Pivot to long format
df_long <- pivot_longer(
  df,
  cols = -c(row, matrix),
  names_to = "column",
  values_to = "value"
)


df_long$column<-factor(df_long$column,levels=c("Intercept","migration","food[granivorous]","food[omnivorous]","food[insectivorous]",
                                               "species(Intercept)","plyo(Intercept)","plyo(migration)","rho"
                                               ))



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
  facet_wrap(~ column, scales = "free_y",nrow=3,ncol=3) +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1)
  ) +
  labs(
    x = "Iteration",
    y = "",
    color = "Algorithm"
  ) +
  theme_minimal()+
  theme(legend.position = "none")


pdf("figure_real1.pdf", width = 9, height = 6.5)
g1
dev.off()

####plot cloglik by tau

X <- model.matrix(parasites~migration+food, dd)
Z<-model.matrix(parasites~migration, dd)
Z1<-model.matrix(parasites~1, dd)

xdf22<-list(Y=dd$parasites,X=X,Z=list(Z,Z1),grouping=list(as.numeric(as.factor(dd$phylogenetic)),as.numeric(as.factor(dd$species))),M=rep(1,nrow(X)),W=rep(1,nrow(X)))

 

#MPL-A1
dm_a1<-get_data_plot_cloglik(xdf22,cfe=mv_multiplier(xdf22),nu=list(3),n_taus=50,fit_pGLM=FALSE,maxiter=15,tol=1e-5,link_fun="logit",save_coef=FALSE,inter_iter=1e3)
df_taus<-data.frame(taus=dm_a1$taus,liks=dm_a1$cloglik,lim=dm_a1$se,opttau=fit_rok$optre$tau)

 

 df_taus$algorithm<-"MPL"
 
df_taus_a<-rbind(df_taus)


df_ref <- df_taus_a %>%
  group_by(algorithm) %>%
  arrange(taus) %>%
  summarise(
    opttau = first(opttau),
    hline = first(liks) + first(lim),
    .groups = "drop"
  )


g_lim<-ggplot(df_taus_a, aes(x = taus, y = liks)) +
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
  theme_minimal()+
  theme(strip.text = element_blank())

pdf("lims1.pdf",height=4,width=13/3)
g_lim
dev.off()

#######single RE

 

#MPL-A1
fit_rok<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic),data=dd,
                    link="logit",
                    penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=FALSE),
                    control=AUGglmmTMBControl(fit_pGLM = FALSE,maxiter = 50,tol=1e-5,save_coef=TRUE))
 

source("mv_MSPAL_illustration_2.R")

fit_ml<-try(glmmTMB(parasites~migration+food+(migration|phylogenetic),data=dd,family="binomial"),silent=TRUE)
fit_blm.t<-try(bglmer(parasites~migration+food+(migration|phylogenetic),data=dd,
                      family="binomial",fixef.prior = "t"),silent=TRUE)
fit_blm.n<-try(bglmer(parasites~migration+food+(migration|phylogenetic),data=dd,
                      family="binomial",fixef.prior = "normal"),silent=TRUE)


np<-8

xdf<-list(Y=dd$parasites,X=X,Z=Z,grouping=as.numeric(as.factor(dd$phylogenetic)),M=rep(1,nrow(X)))

fit_kos <- try(mv_get_MSPAL(start = rep(0, np),
                            data = xdf,
                            mult=c(mv_multiplier(xdf)*2, mv_multiplier(xdf)*2),
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
res<- cbind(resi_ml,resi_rok, resi_bglmert,resi_bglmern,resi_kos)
colnames(res)<-c("ML","MPL","BMT","BMN","SPM")


rownames(res)<-c("Intercept","migration","food[granivourous]","food[omnivourous]","food[insectivourous]"
                 ,"plyo (Intercept)","plyo (migration)","rho")

resx<-xtable(res)
resx

#plots!


####plot coef paths for all methods

 
m2<-fit_rok$coefs$coefs_glmm

nm<-max( nrow(m2))


 



# Convert matrices to data frames
 
df2 <- as.data.frame(m2)
 
 names(df2)<-c(paste0("beta",0:4),paste0("theta",1:3))
# Add row index
 df2$row <- 1:nrow(df2)
 
# Add source label
 df2$matrix <- "MPL"

df2$theta1<-exp(df2$theta1)
df2$theta2<-exp(df2$theta2) 

df2$theta3<-df2$theta3/sqrt(1+df2$theta3**2)

names(df2)<-c("Intercept","migration","food[granivorous]","food[omnivourous]","food[insectivorous]","plyo(Intercept)","plyo(migration)","rho","row","matrix")
  
# Combine
df <- bind_rows(  df2  )

# Pivot to long format
df_long <- pivot_longer(
  df,
  cols = -c(row, matrix),
  names_to = "column",
  values_to = "value"
)

df_long$column<-factor(df_long$column,levels=c("Intercept","migration","food[granivorous]","food[omnivourous]","food[insectivorous]","plyo(Intercept)","plyo(migration)","rho"))




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
  theme_minimal()+ theme(legend.position = "none")



pdf("figure_real.pdf", width = 9/3*4, height = 6.5/3*2)
g1
dev.off()

###plot cloglik by tau
xdf22<-list(Y=dd$parasites,X=X,Z=list(Z),grouping=list(as.numeric(as.factor(dd$phylogenetic))),M=rep(1,nrow(X)),W=rep(1,nrow(X)))

 

#MPL-A1
dm_a1<-get_data_plot_cloglik(xdf22,cfe=mv_multiplier(xdf22),nu=list(3),n_taus=50,fit_pGLM=FALSE,maxiter=15,tol=1e-5,link_fun="logit",save_coef=FALSE,inter_iter=1e3)
df_taus<-data.frame(taus=dm_a1$taus,liks=dm_a1$cloglik,lim=dm_a1$se,opttau=fit_rok$optre$tau)

 

 df_taus$algorithm<-"MPL"
 
df_taus_a<-rbind( df_taus )


df_ref <- df_taus_a %>%
  group_by(algorithm) %>%
  arrange(taus) %>%
  summarise(
    opttau = first(opttau),
    hline = first(liks) + first(lim),
    .groups = "drop"
  )


g_lim<-ggplot(df_taus_a, aes(x = taus, y = liks)) +
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
  theme_minimal()+
  theme(strip.text = element_blank())

pdf("lims.pdf",height=4,width=13/3)
g_lim
dev.off()


