##########imputation######
install.packages("mice")
install.packages("mitolls")
install.packages("jomo")
install.packages("survey")
install.packages("BaBooN")
install.packages("rms")
install.packages("latticeExtra")
install.packages("Hmisc", dependencies = T)
install.packages("ROCR")
install.packages("xfun")
install.packages("psfmi")
install.packages("Rcpp")
install.packages("psych")
library(haven)
library(mice)
library(mitools)
library(jomo)
library(survey)
library(BaBooN)
library(ggplot2)
library(rms)
library(Hmisc)
library(ROCR)
library(xfun)
library(psfmi)
library(Rcpp)
library(psych)
#########
#remove already installed package
#  remove.packages("") 

#################################
###use log xd and log xt#########
###Import Nha Trang's data######
setwd("C:/Users/Noriko/Documents/LSHTM/Data2019/2017_2019serology")
df <- read_dta("dtp2017_2019.dta")
df <- df[,c("id","xd1","xt1","lnxd1","sex","catdip1","cattet01","age1","dtp3card","dtp4card")]
df$lnxt <- log(df$xt1)
df$cattet01 <- factor(df$cattet01)
df$cattet2 <- ifelse(df$xt1>0.16, df$cattet2 <- 1, df$cattet2 <- 0)
df <- subset(df, df$age1< 10)
df$dtp3card <- factor(df$dtp3card)
df$dtp4card <- factor(df$dtp4card)
df$sex <- factor(df$sex)
df$cattet2 <- factor(df$cattet2)
attach(df)

#plot tetanus and diphtheria IgG over age# showing the correlation of age and tetanus IgG, and age and diphtheria IgG
scatter <- ggplot(subset(df, age1<20), aes(x=age1, y=lnxt)) + geom_point(colour ="red") + 
  geom_point(aes(x=age1, y=lnxd1),colour ="blue") + 
  labs(x= "age", y = "IgG (log_scale)", colours = c("tetanus", "diphtheria")) +
  theme_bw()+ theme(legend.position = "bottom") 
scatter  

dim(df)  #198 <15, 139<10yr
sum(is.na(df$dtp3card)) #71 < 15yr, 49 <10yr
sum(is.na(df$dtp4card)) #49 <10yr

cr.mod <- glm(dtp3card ~ lnxt + age1 + sex, family="binomial")
summary(cr.mod)

df.imp1 <-mice (data=df[,c("dtp3card","lnxt", "age1","sex")], m=35,
            method=c("logreg","norm","norm","logreg"),seed=8382,maxit=35)
plot(df.imp1)
xl1=(DTP3 ~ vaccination ~ status)
densityplot(df.imp1, ~dtp3card, xlab =xl1)
stripplot(df.imp1, dtp3card~.imp, col=c("blue",mdc(2)))
imptable = complete(df.imp1, "long", inc= F) #include = T then the raw data is included.
write.csv(imptable, file = "imptable.csv")
imptable$age <- floor(imptable$age1)
table(imptable$age,imptable$dtp3card)

###individual imputed DTP3 completed probability###
imptable$dtp3card <- as.numeric(imptable$dtp3card)
imptable_wide <- reshape(imptable, varying = NULL, timevar = ".imp", idvar =".id", times = seq_along(varying[[1]]), direction = "wide")
imptable_wide$mode <- imptable_wide$dtp3card.1 + imptable_wide$dtp3card.2 + imptable_wide$dtp3card.3 +imptable_wide$dtp3card.4 + imptable_wide$dtp3card.5 + imptable_wide$dtp3card.6 + imptable_wide$dtp3card.7 + 
      imptable_wide$dtp3card.8 + imptable_wide$dtp3card.9 + imptable_wide$dtp3card.10 + imptable_wide$dtp3card.11+ imptable_wide$dtp3card.12 + imptable_wide$dtp3card.13 + imptable_wide$dtp3card.14 +
      imptable_wide$dtp3card.15 + imptable_wide$dtp3card.16 + imptable_wide$dtp3card.17+imptable_wide$dtp3card.18 + imptable_wide$dtp3card.19 + imptable_wide$dtp3card.20+imptable_wide$dtp3card.21+
      imptable_wide$dtp3card.22 + imptable_wide$dtp3card.23 + imptable_wide$dtp3card.24 + imptable_wide$dtp3card.25 + imptable_wide$dtp3card.26 + imptable_wide$dtp3card.27 + imptable_wide$dtp3card.28 +
      imptable_wide$dtp3card.29 + imptable_wide$dtp3card.30 + imptable_wide$dtp3card.31 + imptable_wide$dtp3card.32 + imptable_wide$dtp3card.33 + imptable_wide$dtp3card.34 + imptable_wide$dtp3card.35 
#+ imptable_wide$dtp3card.36 + imptable_wide$dtp3card.37 + imptable_wide$dtp3card.38 + imptable_wide$dtp3card.39 + imptable_wide$dtp3card.40
##if frequency of dtp3 completed (=1) is more than half of the imputed number, that person is counted as completed DTP3##
imptable_wide$mode2 <- ifelse(imptable_wide$mode >= 53, 1, 0)
#mode:when imputation was 40 times >60, when imputation was 35 times >53...
table(imptable_wide$age.1, imptable_wide$mode2)
summary (imptable_wide$mode)

fit1<-with(data=df.imp1,exp=glm(dtp3card ~ lnxt + age1 + sex, family="binomial"))
summary(pool(fit1))
####################################################
####dtp4################
cr.mod <- glm(dtp4card ~ lnxt + age1 + sex, family="binomial")
summary(cr.mod)
df.imp4 <-mice (data=df[,c("dtp4card","lnxt", "age1","sex")], m=35,
                method=c("logreg","norm","norm","logreg"),seed=8382,maxit=35)
plot(df.imp4)
xl2 <- expression(DTP4~ vaccination ~ status)
densityplot(df.imp4, ~dtp4card, xlab=xl2)
stripplot(df.imp4, dtp4card~.imp, col=c("blue",mdc(2)))
imptable4 = complete(df.imp4, "long", inc= F) #include = T then the raw data is included.
write.csv(imptable4, file = "imptable4.csv")
imptable4$age <- floor(imptable4$age1)
table(imptable4$age,imptable4$dtp4card)

fit2<-with(data=df.imp4,exp=glm(dtp4card ~ lnxt + age1 + sex, family="binomial"))
summary(pool(fit2))
###individual imputed DTP4 completed probability###
imptable4$dtp4card <- as.numeric(imptable4$dtp4card)
imptable_wide4 <- reshape(imptable4, varying = NULL, timevar = ".imp", idvar =".id", times = seq_along(varying[[1]]), direction = "wide")
imptable_wide4$mode <- imptable_wide4$dtp4card.1 + imptable_wide4$dtp4card.2 + imptable_wide4$dtp4card.3 +imptable_wide4$dtp4card.4 + imptable_wide4$dtp4card.5 + imptable_wide4$dtp4card.6 + imptable_wide4$dtp4card.7 + 
  imptable_wide4$dtp4card.8 + imptable_wide4$dtp4card.9 + imptable_wide4$dtp4card.10 + imptable_wide4$dtp4card.11+ imptable_wide4$dtp4card.12 + imptable_wide4$dtp4card.13 + imptable_wide4$dtp4card.14 +
  imptable_wide4$dtp4card.15 + imptable_wide4$dtp4card.16 + imptable_wide4$dtp4card.17+imptable_wide4$dtp4card.18 + imptable_wide4$dtp4card.19 + imptable_wide4$dtp4card.20+imptable_wide4$dtp4card.21+
  imptable_wide4$dtp4card.22 + imptable_wide4$dtp4card.23 + imptable_wide4$dtp4card.24 + imptable_wide4$dtp4card.25 + imptable_wide4$dtp4card.26 + imptable_wide4$dtp4card.27 + imptable_wide4$dtp4card.28 +
  imptable_wide4$dtp4card.29 + imptable_wide4$dtp4card.30 + imptable_wide4$dtp4card.31 + imptable_wide4$dtp4card.32 + imptable_wide4$dtp4card.33 + imptable_wide4$dtp4card.34 + imptable_wide4$dtp4card.35 

##if frequency of dtp4 completed (=1) is more than half of the imputed number, that person is counted as completed DTP4##
imptable_wide4$mode2 <- ifelse(imptable_wide4$mode >= 53, 1, 0)
#mode:when imputation was 40 times >60, when imputation was 35 times >53...
table(imptable_wide4$age.1, imptable_wide4$mode2)
summary (imptable_wide$mode)
write.csv(imptable_wide4, file = "imp4")

####estimate probability of vaccinated 3 doses from each tetanus value
####using the pooled analysis from imputation
coverageest <- subset(imptable, imptable$.imp==0)
coverageest$sex <- as.numeric(as.character(coverageest$sex))
coverageest$est <- 4.0612973 + 0.9305096*coverageest$lnxt - 0.1142889*coverageest$age1 - 0.1134628*coverageest$sex
View(coverageest)
coverageest$est_probablity <- coverageest$est/(coverageest$est+1)
describeBy(coverageest$est_probablity, group = coverageest$age)

###how tetanus IgG wanes over time
summary(waningtet <- lm(lnxt ~ age1 + sex)) 
coverageest$est_lnxt <- waningtet[["coefficients"]][["(Intercept)"]] + waningtet[["coefficients"]][["age1"]]*coverageest$age1 + waningtet[["coefficients"]][["sex2"]]*coverageest$sex
coverageest$est2 <- 4.0102724 + 0.9199392*coverageest$est_lnxt -0.1115565*coverageest$age1 - 0.1062251*coverageest$sex
coverageest$est_probablity2 <- coverageest$est2/(coverageest$est2+1)
describeBy(coverageest$est_probablity2, group = coverageest$age)

###use df.imp1 for the results##
##dif. imp2 and df.imp3 was examined whether the model improved without sex 
## or with diphtheria IgG  it does not really improved ###
df.imp2 <-mice (data=df[,c("dtp3card","lnxt","age1")], m=40,
                method=c("logreg","norm","norm"),seed=8382,maxit=40)
plot(df.imp2)
densityplot(df.imp2, ~dtp3card)
stripplot(df.imp2, dtp3card~.imp, col=c("blue",mdc(2)))
imp2table = complete(df.imp2, "long", inc= TRUE)
imp2table$age <- floor(imp2table$age1)
table(imp2table$age, imp2table$dtp3card)
write.csv(imp2table, file = "imp2table.csv")

fit2<-with(data=df.imp2,exp=glm(dtp3card ~ lnxt + age1, family="binomial"))
summary(pool(fit2))

df.imp3<-mice(data=df[,c("dtp3card","lnxt","lnxd1","sex", "age1")],m=40,
              method=c("logreg","norm","norm","logreg", "norm"),seed=8382, maxit=40)
plot(df.imp3)
densityplot(df.imp3, ~dtp3card)
stripplot(df.imp3, dtp3card~.imp, col=c("blue",mdc(2)))

imp3table = complete(df.imp3, "long", inc= TRUE)
write.csv(imp3table, file = "imp3table.csv")

imp3table$age <- floor(imp3table$age1)
fit3<-with(data=df.imp3,exp=glm(dtp3card ~ lnxt + lnxd1 + sex + age1, family="binomial"))
summary(pool(fit3))


###create AUC curve####
dd <- datadist(df)
options(datadist = "dd")

model=lrm(formula=dtp3card ~ lnxt + age1 + sex,
           data=df,
           x=TRUE,
           y=TRUE)
val=validate(model,B=10000) ###validating model
c_opt_corr <- 0.5 * (val[1, 5] + 1) ##c statistic/AUC
cal <- calibrate(model,B=10000)   ##calibration plots(accuracy)
plot(cal,xlab="Predicted probability",ylab="Observed proportion")


#####
library(psfmi)
###calculate c statistic(AUC) and Brier score of the pooled model
#plot observed probabilities vs predicted probability in individual data##

perf <- pool_performance(data=imptable, nimp=35, impvar=".imp", 
                         Outcome = "dtp3card", predictors = c("lnxt","sex","age1"), 
                         cal.plot=FALSE, plot.indiv=FALSE)
perf

install.packages("CalibrationCurves")
install.packages("CalibrationCurves", dependencies=TRUE, libPaths=c("/Users/local/lib/R/site-library"))
install.packages("rms")
library(CalibrationCurves)

#CalibrationCurves are not available in my R version???
install.packages("devtools")
require(devtools)
#https://git-scm.com/download/win
install_git("https://github.com/BavoDC/CalibrationCurves")

install.packages("val.prob.ci.2")


#Generate the imputed dataset using MICE
df_minvar <- df[,c("lnxt","age1","sex","dtp3card")] #restrict the variable
imp <- mice(df_minvar, maxit=0)  
predM <- imp$predictorMatrix  
meth <- imp$method            

##set the different methods
meth[c("lnxt")]="norm"
meth[c("age1")]="norm"
meth[c("sex")]="logreg"
meth[c("dtp3card")]="logreg"

#predM[c("id", "xd1","xt1","catdip1","cattet0","dtp4","cattet2")]=0   ### keep variables which does not need to change to 0

##fit the imputation model
imp_auc3 <- mice(df_minvar, maxit =100, 
             predictorMatrix = predM, 
             method = meth, print =  FALSE, m=35)

#Fit a (logistic) model using MICE output
model_imp <- Hmisc::fit.mult.impute(dtp3card ~ age1 + sex + lnxt,
                                rms::lrm, #logistic regression model using rms package
                                xtrans = imp_auc3,
                                data=df,
                                x=T, 
                                y=T,
                                linear.predictors = T,
                                pr=F)

###HMISC
#Generate the imputed dataset using HMISC
Imputed_dataset <-  Hmisc::aregImpute(formula = ~dtp3card + lnxt + age1 + sex  ,
                                      data = df,
                                      burnin = 10,
                                      n.impute= 35,
                                      pr = F )

#Fit a (logistic) model using HMISC output
model_imp2 <- Hmisc::fit.mult.impute(dtp3card ~ age1 + sex + lnxt,
                                 rms::lrm, #logistic regression model using rms package
                                 xtrans = Imputed_dataset,
                                 data=df,
                                 x=T, 
                                 y=T,
                                 linear.predictors = T,
                                 pr=F)

###plot1 imputation
#use MICE package
source("mycalfunction2.R")
#source('~/LSHTM/Publication/Trivariate serology/mycalfunction2.R', echo=TRUE)
cal4 <- calibrate(model_imp,B=10000)

#get slope and intercept using RSM package
#DXY could be converted to c-statistics but no 95%CI, B =Brier score, 
cal5 <- validate(model_imp,B=10000) 
cal5

#create graph predicted probability vs observed proportion
myfun(cal4,xlab="Predicted probability",ylab="Observed proportion",main="Multiple imputation model",cex.axis = 1,cex.lab=1,cex.main=1,cex=0.7, group = 10) 
text(x=0.2,y=0.99,"C-statistic=0.835(0.687-0.920)",cex=1)
text(x=0.2,y=0.90,"Intercept=0.096",cex=1)
text(x=0.2,y=0.80,"Slope=0.894",cex=1)
text(x=0.2,y=0.70,"Brier=0.150",cex=1)

###plot2 imputation
source("mycalfunction.R")
cal4 <- calibrate(model_imp2,B=10000)
myfun2(cal4,xlab="Predicted probability",ylab="Observed proportion",main="Multiple imputation model",cex.axis = 1,cex.lab=1,cex.main=1,cex=0.7,group=10) 
text(x=0.53,y=0.99,"C-statistic=0.82(0.67-0.92)",cex=0.7)
text(x=0.51,y=0.95,"Intercept=0.093",cex=0.7)
text(x=0.5,y=0.90,"Slope=0.897",cex=0.7)
text(x=0.5,y=0.85,"Brier=0.230",cex=0.7)


########################################
###########use xt, (xd) age and sex#####
#### AUC : complete data analysis ######
model3=glm(formula=dtp3card ~ lnxt + age1+sex,
           data=df,
           family = binomial(link="logit"))
summary(model3)

##Internal validation 

##method 1.using bootstrapping method
##using the complete case only

#Step 1 calculate insample AUC
prob = predict(model3, type='response', newdata = df)
pred = prediction(prob, df$dtp3card)

# AUC
InsampleAuc = performance(pred,"auc")@y.values[[1]][1]
#not enough distinct predictions to compute area under the ROC curve(error)
# plot the ROC curve
perf <- performance(pred,"tpr","fpr")
plot(perf, col="navyblue", cex.main=1,
     main= paste("Logistic Regression ROC Curve: AUC =", round(auc,3)))
abline(a=0, b = 1, col='darkorange1')

#Step2  generate bootstraped aucs
R=10000
n=nrow(datanew)
B=matrix(nrow=R,ncol=2,dimnames = 
           list(paste("Sample",1:R),c("auc_orig","auc_boot")))

for(i in 1:R){
  obs.boot=sample(x=1:n,size=n,replace=T) ##random(bootstrapped) sample
  data.boot=datanew[obs.boot,]
}

  ##fit model on bootstrapped sample
  glm.boot=glm(model3$formula,
               data=data.boot,
               family=binomial)
               
##run model
  model3=glm(formula=dtp3card~lnxt + lnxd1 + sex + age1,
                          data=datanew,
                          family = binomial(link="logit"))
  summary(model3)



  
  #################################
  ### import Nha Trang's data #####
  setwd("C:/Users/Noriko/Documents/LSHTM/Data2019/2017_2019serology")
  df <- read_dta("dtp2017_2019.dta")
  df <- df[,c("id","xd1","xt1","lnxd1","sex","catdip1","cattet01","age1","dtp3card","dtp4")]
  df$lnxt <- log(df$xt1)
  df$cattet01 <- factor(df$cattet01)
  df$cattet2 <- ifelse(df$xt1>0.16, df$cattet2 <- 1, df$cattet2 <- 0)
  df <- subset(df, df$age1< 20) 
  df$dtp3card <- factor(df$dtp3card)
  df$sex <- factor(df$sex)
  df$cattet2 <- factor(df$cattet2)
  attach(df)
  cr.mod <- glm(dtp3card ~ xt1, family="binomial")
  cr.mod <- glm(dtp3card ~ xt1 + age1, family="binomial")
  cr.mod <- glm(dtp3card ~ xt1 + age1 + sex, family="binomial")
  cr.mod <- glm(dtp3card ~ cattet2, family="binomial")
  cr.mod <- glm(dtp3card ~ cattet2+ age1 + sex, family="binomial")
  cr.mod <- glm(dtp3card ~ lnxt  + age1 + sex, family="binomial")
  summary(cr.mod)
  
  dim(df)
  sum(is.na(df$dtp3card))
  #####check missing pattern
  df.complete <- rep(1,259) #198 <15yr
  df.complete[is.na(dtp3card)]<-0
  table(df.complete)
  df$df.complete <- factor(df$df.complete)
  
  
  summary(glm(df.complete~xt1 ))
  summary(glm(df.complete~age1 ))
  summary(glm(df.complete~sex))
  summary(glm(df.complete~xt1 + age1 + sex))
  summary(glm(df.complete~lnxt))
  summary(glm(df.complete~lnxt + age1))
  summary(glm(df.complete~cattet2))
  summary(glm(df.complete~cattet2 + age1))
  summary(glm(df.complete~cattet2 + age1 + sex))
  
  box1 <- ggplot(df, aes(x=df.complete, y=age1)) + geom_boxplot() + labs(x ="0 missing, 1 not-missing", y ="age")
  box1
  ggsave("nhatrang_age_complete.jpeg", width = 4, height = 3, dpi =500)
  
  box2 <- ggplot(df, aes(x=df.complete, y=xt1)) + geom_boxplot() + labs(x ="0 missing, 1 not-missing", y ="tetanus IgG")
  box2
  ggsave("nhatrang_xt_complete.jpeg", width = 4, height = 3, dpi =500)
  
  # Impute:
  df.imp2<-mice(data=df[,c("dtp3card","lnxt","sex", "age1")],m=40,
                method=c("logreg","norm","norm", "logreg", "norm"),seed=8382, maxit=40)
  plot(df.imp2)
  savePlot(filename = "nhatrang_imp2", type=c("jpeg"), device = dev.cur(), restoreConsole = TRUE)
  densityplot(df.imp2, ~dtp3card)
  stripplot(df.imp2, dtp3card~.imp, col=c("blue",mdc(2)))
  ###imputed data table
  imp2table = complete(df.imp2, "long", inc= TRUE)
  write.csv(imp2table, file = "imp2table.csv")
  
  ##pooled regression
  fit1<-with(data=df.imp2,exp=glm(dtp3card ~ lnxt + sex +age1 ,family="binomial"))
  summary(pool(fit1))
  
  #### categorical
  df.imp2<-mice(data=df[,c("dtp3card","cattet2","sex", "age1")],m=40,
                method=c("logreg","norm","logreg", "norm"),seed=8382, maxit=40)
  plot(df.imp2)
  densityplot(df.imp2, ~dtp3card)
  stripplot(df.imp2, dtp3card~.imp, col=c("blue",mdc(2)))
  
  fit1<-with(data=df.imp2,exp=glm(dtp3card ~ cattet2 + sex + age1, family="binomial"))
  summary(pool(fit1))
  
  #### export results
  imp1table = complete(df.imp1, "long", inc= TRUE)
  write.csv(imp1table, file = "imp1table.csv")
                         
