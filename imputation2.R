################################
####apply prediction model for the Quang Ngai data#######
####Quang Ngai Data########
setwd("C:/Users/Noriko/Documents/LSHTM/Data2019/QNgai2019")
df <- read_dta("qn2019analysis.dta")
df <- df[,c("id","xt","logxt","logxd","xp","sex","catxd","catxt","catxt2","agei","age2","bcg", "dtp3","dtp3mis","dtp4","dtp4mis","dtpdose","commune2","village","dur")]
df <- subset(df, df$agei< 10)
df$commune2 <- factor(df$commune2)
df$dtp3 <- factor(df$dtp3)
df$dtpdose <- factor(df$dtpdose)
df$dtp3mis <- factor(df$dtp3mis)
df$dtp4mis <- factor(df$dtp4mis)
df$catxt2 <- factor(df$catxt2)
df$catxd <- factor(df$catxd)
df$sex <- factor(df$sex)
#df = df %>% filter(!is.na(df$catxt2))
df <- subset(df, !is.na(df$catxt2))
#df <- subset(df, df$dtp4!=1)

scatter <- ggplot(subset(df, df$agei<10), aes(x=age2, y=logxt, colour = dtpdose)) + geom_point() + 
  labs(x= "age", y = "IgG (log_scale)", colours = c("tetanus", "diphtheria")) +
  theme_bw()+ theme(legend.position = "bottom") 
scatter  

attach(df)
box <- ggplot(df,aes(x=dtp3mis, y=agei)) + geom_boxplot()
box
ggsave("box.jpeg")

box2 <- ggplot(df,aes(x=dtp3mis, y=agei)) + geom_boxplot()
box2
##complete dataset model
cr.mod <- glm(dtp3 ~ logxt +age2 + sex , family="binomial")
summary(cr.mod)
cr.mod <- glm(dtp3 ~ logxt + age2 + sex + commune2, family="binomial")
summary(cr.mod)
cr.mod <- glm(dtp4 ~ logxt +age2 + sex , family="binomial")
summary(cr.mod)
cr.mod <- glm(dtp4 ~ logxt + age2 + sex + commune2, family="binomial")
summary(cr.mod)
cr.mod <- glm(dtp3mis ~ logxt +age2 + sex , family="binomial")
summary(cr.mod)
cr.mod <- glm(dtp4mis ~ logxt + age2+ sex , family="binomial")
summary(cr.mod)

summary(lm(xt ~ agei + dtp3mis + sex + commune2))
summary(lm(logxt ~ agei + dtp3mis + sex + commune2))

dim(df)
sum(is.na(df$dtp3))
#####check missing pattern
df.complete <- rep(1,556)
df.complete[is.na(dtp3)]<-0
table(df.complete)
summary(glm(df.complete~dtp3))
summary(glm(df.complete~logxt))
summary(glm(df.complete~agei))
summary(glm(df.complete~sex))
summary(glm(df.complete~commune2))
summary(glm(df.complete~logxt + agei))
summary(glm(df.complete~logxt + agei + sex + commune2))

# Impute:
df.imp1<-mice(data=df[,c("dtp3","logxt","agei","sex")],m=20,
              method=c("logreg","norm","norm","logreg"),seed=8382, maxit=20)
plot(df.imp1)
densityplot(df.imp1, ~dtp3)
stripplot(df.imp1, dtp3~.imp, col=c("blue",mdc(2)))

fit1<-with(data=df.imp1,exp=glm(dtp3 ~ logxt + agei +sex, family="binomial"))
summary(pool(fit1))
##logxt not significant###

#### export results
imp1table = complete(df.imp1, "long", inc= TRUE)
write.csv(imp1table, file = "imp1table.csv")

#external validation
####estimate probability of vaccinated 3 doses from each tetanus value
####using the pooled analysis from imputation
#coverageest <- subset(imp1table, imp1table$.imp==0)
coverageest <- df
coverageest$sex <- as.numeric(as.character(coverageest$sex))
coverageest$est <- 3.932 + 0.898*coverageest$logxt -0.094*coverageest$age2 -0.172*coverageest$sex
View(coverageest)
coverageest$est_probablity <- coverageest$est/(coverageest$est+1)
#mode
coverageest$est_dtp3[coverageest$est_probablity>0.5]<- 1
coverageest$est_dtp3[coverageest$est_probablity<=0.5]<- 0
describeBy(coverageest$est_dtp3, group = coverageest$agei)
table(coverageest$est_dtp3, coverageest$agei)
##estimate DTP4, using the pooled analysis from imputation###
coverageest$est4 <- 2.439 + 0.559*coverageest$logxt -0.270*coverageest$age2 + 0.335 *coverageest$sex
coverageest$est_probablity4 <- coverageest$est4/(coverageest$est4+1)
#mode
coverageest$est_dtp4[coverageest$est_probablity4>0.5]<- 1
coverageest$est_dtp4[coverageest$est_probablity4<=0.5]<- 0
table(coverageest$est_dtp4, coverageest$agei)
