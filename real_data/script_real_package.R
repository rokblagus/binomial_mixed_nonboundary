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

#MPL-A2(1)

fit_rok_1<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=dd,
                      link="logit",
                      penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=FALSE),
                      control=AUGglmmTMBControl(fit_pGLM = TRUE,maxiter = 1,tol=1e-5,save_coef=TRUE))

#MPL-A1
fit_rok<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=dd,
                      link="logit",
                      penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=FALSE),
                      control=AUGglmmTMBControl(fit_pGLM = FALSE,maxiter = 50,tol=1e-5,save_coef=TRUE))

#MPL-A2

fit_rok_frth<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic)+(1|species),data=dd,
                      link="logit",
                      penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=FALSE),
                      control=AUGglmmTMBControl(fit_pGLM = TRUE,maxiter = 50,tol=1e-5,save_coef=TRUE))


#summarize the data

mod<-fit_rok$fit$fit
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

mod<-fit_rok_1$fit$fit

#sd for 2nd RE
vr<-mod$sdr$cov.fixed[6+3,6+3]
est<-mod$sdr$par.fixed[6+3]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 1st RE

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

resi_rok2<-c(fxf,
             paste0(round(s1$Estimate,3)," (",round(s1$SE,3),")"),
             paste0(round(s3.1$Estimate,3)," (",round(s3.1$SE,3),")"),
             paste0(round(s3.2$Estimate,3)," (",round(s3.2$SE,3),")"),

             paste0(round(s3.12$Estimate,3)," (",round(s3.12$SE,3),")")

)


mod<-fit_rok_frth$fit$fit
#sd for 2nd RE
vr<-mod$sdr$cov.fixed[6+3,6+3]
est<-mod$sdr$par.fixed[6+3]
names(est)<-"theta"

s1<-deltaMethod(est,"exp(theta)",vr)

#sd for 1st RE

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


####plots!


####plot coef paths for all methods

m1<-fit_rok_frth$fit$coefs

m2<-fit_rok$fit$coefs

nm<-max(nrow(m1),nrow(m2))


m3<- matrix(fit_rok_1$fit$coefs,nrow=nm,ncol=ncol(m1),byrow=TRUE)





# Convert matrices to data frames
df1 <- as.data.frame(m1)

df2 <- as.data.frame(m2)
df3 <- as.data.frame(m3)

names(df1)<-names(df2)<-names(df3)<-c(paste0("beta",0:4),paste0("theta2",1:3),"theta1")
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

ses<-sqrt(diag(fit_rok$fit$fit$sdr$cov.fixed))

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
  theme_minimal()


pdf("figure_real1.pdf", width = 9, height = 6.5)
g1
dev.off()

####plot cloglik by tau

X <- model.matrix(parasites~migration+food, dd)
Z<-model.matrix(parasites~migration, dd)
Z1<-model.matrix(parasites~1, dd)

xdf22<-list(Y=dd$parasites,X=X,Z=list(Z,Z1),grouping=list(as.numeric(as.factor(dd$phylogenetic)),as.numeric(as.factor(dd$species))),M=rep(1,nrow(X)),W=rep(1,nrow(X)))

#MPL-A2(1)
dm_1<-get_data_plot_cloglik(xdf22,cfe=mv_multiplier(xdf22),nu=list(3),n_taus=50,fit_pGLM=TRUE,maxiter=1,tol=1e-6,link_fun="logit",save_coef=FALSE,inter_iter=1e3)
df_taus_1<-data.frame(taus=dm_1$taus,liks=dm_1$cloglik,lim=dm_1$se,opttau=fit_rok_1$optre$tau)


#MPL-A1
dm_a1<-get_data_plot_cloglik(xdf22,cfe=mv_multiplier(xdf22),nu=list(3),n_taus=50,fit_pGLM=FALSE,maxiter=15,tol=1e-5,link_fun="logit",save_coef=FALSE,inter_iter=1e3)
df_taus<-data.frame(taus=dm_a1$taus,liks=dm_a1$cloglik,lim=dm_a1$se,opttau=fit_rok$optre$tau)

#MPL-A2
dm_a2<-get_data_plot_cloglik(xdf22,cfe=mv_multiplier(xdf22),nu=list(3),n_taus=50,fit_pGLM=TRUE,maxiter=15,tol=1e-5,link_fun="logit",save_coef = TRUE,use_previous = FALSE)
df_taus_f<-data.frame(taus=dm_a2$taus,liks=dm_a2$cloglik,lim=dm_a2$se,opttau=fit_rok_frth$optre$tau)


df_taus_1$algorithm<-"MPL-A2(1)"
df_taus$algorithm<-"MPL-A1"
df_taus_f$algorithm<-"MPL-A2"

df_taus_a<-rbind(df_taus_1,df_taus,df_taus_f)


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
  theme_minimal()

pdf("lims1.pdf",height=4,width=13)
g_lim
dev.off()

#######single RE



#MPL-A2(1)
fit_rok_1<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic),data=dd,
                      link="logit",
                      penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=FALSE),
                      control=AUGglmmTMBControl(fit_pGLM = TRUE,maxiter = 1,tol=1e-5,save_coef=TRUE))


#MPL-A1
fit_rok<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic),data=dd,
                      link="logit",
                      penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=FALSE),
                      control=AUGglmmTMBControl(fit_pGLM = FALSE,maxiter = 50,tol=1e-5,save_coef=TRUE))

#MPL-A2
fit_rok_frth<-AUGglmmTMB(parasites~migration+food+(migration|phylogenetic),data=dd,
                    link="logit",
                    penOpt = AUGglmmTMBPenalty(autrepen =TRUE,plot=FALSE),
                    control=AUGglmmTMBControl(fit_pGLM = TRUE,maxiter = 50,tol=1e-5,save_coef=TRUE))


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


mod<-fit_rok_1$fit$fit

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



mod<-fit_rok_frth$fit$fit

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


mod<-fit_rok$fit$fit

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
resx

#plots!


####plot coef paths for all methods

m1<-fit_rok_frth$fit$coefs

m2<-fit_rok$fit$coefs

nm<-max(nrow(m1),nrow(m2))


m3<- matrix(fit_rok_1$fit$coefs,nrow=nm,ncol=ncol(m1),byrow=TRUE)




# Convert matrices to data frames
df1 <- as.data.frame(m1)

df2 <- as.data.frame(m2)
df3 <- as.data.frame(m3)

names(df1)<-names(df2)<-names(df3)<-c(paste0("beta",0:4),paste0("theta",1:3))
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

ses<-sqrt(diag(fit_rok$fit$fit$sdr$cov.fixed))

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

###plot cloglik by tau
xdf22<-list(Y=dd$parasites,X=X,Z=list(Z),grouping=list(as.numeric(as.factor(dd$phylogenetic))),M=rep(1,nrow(X)),W=rep(1,nrow(X)))

#MPL-A2(1)
dm_1<-get_data_plot_cloglik(xdf22,cfe=mv_multiplier(xdf22),nu=list(3),n_taus=50,fit_pGLM=TRUE,maxiter=1,tol=1e-6,link_fun="logit",save_coef=FALSE,inter_iter=1e3)
df_taus_1<-data.frame(taus=dm_1$taus,liks=dm_1$cloglik,lim=dm_1$se,opttau=fit_rok_1$optre$tau)


#MPL-A1
dm_a1<-get_data_plot_cloglik(xdf22,cfe=mv_multiplier(xdf22),nu=list(3),n_taus=50,fit_pGLM=FALSE,maxiter=15,tol=1e-5,link_fun="logit",save_coef=FALSE,inter_iter=1e3)
df_taus<-data.frame(taus=dm_a1$taus,liks=dm_a1$cloglik,lim=dm_a1$se,opttau=fit_rok$optre$tau)

#MPL-A2
dm_a2<-get_data_plot_cloglik(xdf22,cfe=mv_multiplier(xdf22),nu=list(3),n_taus=50,fit_pGLM=TRUE,maxiter=15,tol=1e-5,link_fun="logit",save_coef = TRUE,use_previous = FALSE)
df_taus_f<-data.frame(taus=dm_a2$taus,liks=dm_a2$cloglik,lim=dm_a2$se,opttau=fit_rok_frth$optre$tau)


df_taus_1$algorithm<-"MPL-A2(1)"
df_taus$algorithm<-"MPL-A1"
df_taus_f$algorithm<-"MPL-A2"

df_taus_a<-rbind(df_taus_1,df_taus,df_taus_f)


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
  theme_minimal()

pdf("lims.pdf",height=4,width=13)
g_lim
dev.off()


