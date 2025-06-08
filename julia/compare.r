# load libraries
library(pacman)
p_load(dplyr,ggplot2,survival,survminer,ggsurvfit,RISCA)
# load df
df=read.csv("C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\df.csv")
# survival by groups
fit=survfit2(Surv(time,y=='true')~factor(trt),
             df)
plot(fit)
# cox model
model=coxph(Surv(time, y=='true') ~ trt + age,df,weights = weight )
summary(model)
# hazard
bh=basehaz(model)
bh$coeff=coef(model)
# create survival probabilities
bh=bh %>% 
  mutate(s0=exp(-hazard),
         prob=s0^(exp(coeff)))
# plot
plot(bh$time,bh$prob)
lines(bh$time,bh$prob)
# calculate predictions
newdata=data.frame(trt=as.factor('true'))
preds=survfit(model,df)
plot(preds)
# compare
bh$survift=preds$surv
bh=bh %>% 
  mutate(diff=survift-prob)
sum(bh$diff)
# create df for predictions
dfpred=data.frame(time=preds$time,
                  surv=preds$surv)
# create csv file
write.csv(dfpred,"C:\\Users\\P095206\\OneDrive - Amsterdam UMC\\Shared material with Elia\\PhD\\Project X - Matching methods in small sample sizes\\Data\\julia\\compare.csv")
### ipw
df$trt=ifelse(df$trt=='true',1,0)
weights=ipw::ipwpoint(exposure = trt,
                      family = "binomial",
                      link = "logit",
                      numerator= ~1,
                      denominator = ~ age,
                      data = df)
# model
model=coxph(Surv(time, y=='true') ~ trt,df,weights = weights$ipw.weights)
# yep
### gcomp
df=df %>% mutate(across(c(trt,y),~ifelse(.=='true',1,0)))
# model
model=coxph(Surv(time, y) ~ trt + age,df)
# risca
fit.risca=gc.survival(object=model,group='trt',times='time',failures='y',
                      max.time = max(df$time),data=df)
# results
ate.risca=as.numeric(fit.risca$logHR[1])
# by hand
# dataset
trt=df %>% mutate(trt=1)
ctrl=df %>% mutate(trt=0)
# predict
predtrt=survfit(model,trt); meantrt=rowMeans(as.matrix(predtrt$surv))
predctrl=survfit(model,ctrl); meanctrl=rowMeans(as.matrix(predctrl$surv))
# 
ate.hand=log(-log(meantrt[1000])/-log(meanctrl[1000]))

