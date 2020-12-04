#z.1

rm(list=ls())
data=read.csv2('Acemoglu.csv', header=TRUE, sep=',', dec='.')
model=lm(logpgp95~avexpr, data=data)
summary(model)

library(stargazer)
ols=lm(logpgp95~avexpr+logem4, data=data)
ehat<-ols$fitted.values
sst<-lm(logpgp95~avexpr+ehat,data=data)

library(AER)
IV=ivreg(logpgp95~avexpr+avexpr|logem4, data=data)
summary(IV)
summary(IV,diagnostics=T)
linearHypothesis(ols,"logem4=0")

model2=lm(logpgp95~avexpr+logem4+euro1900, data=data)
linearHypothesis(model2,c("logem4=0",'euro1900=0'))

IV2=ivreg(logpgp95~avexpr+avexpr|(logem4+euro1900),data=data)
summary(IV2,diagnostics=T)

model3=lm(logpgp95~avexpr+lat_abst, data=data)
summary(model3)

IV3=ivreg(logpgp95~avexpr+(avexpr+lat_abst)|(logem4+euro1900)+lat_abst,data=data)
summary(IV3,diagnostics=T)

library(wbstats)
dataGDP<-wb(indicator='NY.GDP.PCAP.PP.CD') #GDPdla2018roku
dataGDP<-subset(dataGDP, date==2018)
dataGDP<-dataGDP[,-c(4:7)]
dataGDP<-dataGDP[,-2]
data<-cbind(data,dataGDP[match(data[,'shortnam'],dataGDP[,'iso3c']),'value'])#wskaŸnikefektywnoœciinstytucjirz¹dowychdla2018roku
dataIns<-wb(indicator='GE.PER.RNK')
dataIns<-subset(dataIns, date==2018)
dataIns<-dataIns[,-c(4:7)]
dataIns<-dataIns[,-2]
data<-cbind(data,dataIns[match(data[,'shortnam'],dataIns[,'iso3c']),'value'])
colnames(data)[8]<-'GDP2018'
colnames(data)[9]<-'GovEff'

write.csv(data, file='dataRS82640PD.csv')

newmodel1=lm(log(GDP2018)~GovEff,data=data)
summary(newmodel1)
newmodel2=lm(log(GDP2018)~GovEff+lat_abst,data=data)
summary(newmodel2)
newIV1=ivreg(log(GDP2018)~GovEff+GovEff|(logem4+euro1900),data=data)
summary(newIV1, diagnostics=T)
newIV2=ivreg(log(GDP2018)~GovEff+lat_abst+(GovEff+lat_abst)|(logem4+euro1900),data=data)
summary(newIV2, diagnostics=T)


###########################################
#z.2
rm(list=ls())
data2=read.csv2('InternationalTradePoland.csv', header=TRUE, sep=',', dec='.')
data2=na.omit(data2)
N<-nrow(data2)

year<-data2$year
im<-log(data2$IM)
ex<-log(data2$EX)
dd<-log(data2$DD)
ddf<-log(data2$DD_F)
reer<-log(data2$REER)
library(Hmisc)
Lim<-Lag(im,-1)
Lex<-Lag(ex, -1)

library(systemfit)
library(dynlm)
modelim=dynlm(im~dd+reer+Lim)
modelex=dynlm(ex~ddf+reer+Lex)
tb=ex-im
summary(modelim)
summary(modelex)


FSTIM=dynlm(im~dd+reer+ddf)
FSTEX=dynlm(ex~dd+reer+ddf)
IMhat=FSTIM$fitted.values
EXhat=FSTEX$fitted.values

SSTIM=lm(im~dd+reer+EXhat+Lim)
SSTEX=lm(ex~ddf+reer+IMhat+Lex)

eqIM<-im ~ dd + reer + ex + Lim
eqEX<-ex ~ ddf + reer+ im + Lex
eqtb<-ex-im

inst<- ~ dd+ddf+reer+ex+im+Lim+Lex

system<-list(import=eqIM)
system2SLS<-systemfit(system,'2SLS',inst=inst)
summary(system2SLS)

system<-list(export=eqEX)
system2SLS2<-systemfit(system,'2SLS',inst=inst)
summary(system2SLS2)

system<-list(import=eqIM, export=eqEX)
system3SLS=systemfit(system,"3SLS",inst=inst,maxiter=1)
summary(system3SLS)
dwtest(modelex)
dwtest(modelim)

system<-list(IM=eqIM, EX=eqEX)
system2SLS<-systemfit(system,'2SLS',inst=inst)
bgtest(modelex,order=4)

names.endo=c('im','ex')
names.exo=c('Intercept','dd','ddf','reer','Lex','Lim')
K=length(names.endo) # liczba zamiennych endogneicznych
M=length(names.exo) # Liczba zmiennych egzogenicznych

#Macierze parametrów - generujemy puste macierze
Gamma=diag(K)
colnames(Gamma)=names.endo
rownames(Gamma)=names.endo
B=matrix(0,ncol=K,nrow=M)
colnames(B)=names.endo
rownames(B)=names.exo
coef<-as.data.frame(t(coef(system2SLS)))

Gamma["im","ex"]=-coef$IM_ex
Gamma["ex","im"]=-coef$EX_im
B['Intercept','ex']=-coef$`EX_(Intercept)`
B['Intercept','im']=-coef$`IM_(Intercept)`
B['dd','im']=-coef$IM_dd
B['ddf','ex']=-coef$EX_ddf
B['reer','ex']=-coef$EX_reer
B['reer','im']=-coef$IM_reer
B['Lex','ex']=-coef$EX_Lex
B['Lim','im']=-coef$IM_Lim

PI=-B%*%solve(Gamma)
PIX=PI[1:4,]

VARPHI=matrix(0,ncol=K,nrow=K)
colnames(VARPHI)=names.endo
rownames(VARPHI)=paste("L",names.endo,sep="")

VARPHI[rownames(PI)[5:6],]=PI[5:6,]
VARPHI

library(expm)

s <- 10
auuu<-PIX%*%(VARPHI%^%s)
auuu


# Mnoznik dlugookresowy/rownowagowy
LR=PIX%*%solve(diag(K)-VARPHI)
LR


#Skumulowany mnoznik dynamiczny
s <- 10
skumulowany<-LR%*%(diag(K)-VARPHI%^%s)
skumulowany


#Mnoznik dynamiczny - funkcje reakcji na impuls
H=20
IRF=array(dim=c(4,2,H+1))
IRF[,,1]=PIX
for(i in 2:(H+1)){
  IRF[,,i]=IRF[,,(i-1)]%*%VARPHI
}

test<-data.frame(IRF)
test1<-data.frame(t(test[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41)]))
test2<-data.frame(t(test[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42)]))
IM<-c(mean(test1$X1),mean(test1$X2),mean(test1$X3),mean(test1$X4))
EX<-c(mean(test2$X1),mean(test2$X2),mean(test2$X3),mean(test2$X4))
NAME<-c('Intercept','dd','ddf','reer')
au<-as.matrix(t(rbind(NAME,IM,EX)))
au

H=20
CIRF=array(dim=c(4,2,H+1))
for(i in 1:(H+1)){
  CIRF[,,i]=LR%*%(diag(K)-VARPHI%^%i)
}

test<-data.frame(CIRF)
test1<-data.frame(t(test[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41)]))
test2<-data.frame(t(test[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42)]))
IM<-c(mean(test1$X1),mean(test1$X2),mean(test1$X3),mean(test1$X4))
EX<-c(mean(test2$X1),mean(test2$X2),mean(test2$X3),mean(test2$X4))
NAME<-c('Intercept','dd','ddf','reer')
au2<-as.matrix(t(rbind(NAME,IM,EX)))
au2

CIRF202<-as.numeric(au2[,3])
CIRF202
CIRF201<-as.numeric(au2[,2])
IRF202<-as.numeric(au2[,3])
IRF201<-as.numeric(au[,2])
dynamiczny<-auuu[,2]-auuu[,1]
dlugookresowy<-LR[,2]-LR[,1]
skumulowany<-skumulowany[,2]-skumulowany[,1]
IRF20<-as.vector(t(IRF202-IRF201))
CIRF20<-as.vector(t(CIRF202-CIRF201))
saldo<-rbind(dynamiczny, dlugookresowy, skumulowany, IRF20, CIRF20)
saldo

###################

eqIM<-im ~ dd + reer + ex + Lim
eqEX<-ex ~ ddf + reer+ im + Lex

inst<- ~ reer+ex+im+Lim+Lex
names.endo=c('im','ex')
names.exo=c('Intercept','dd','ddf','reer','Lex','Lim')
K=length(names.endo) # liczba zamiennych endogneicznych
M=length(names.exo) # Liczba zmiennych egzogenicznych

s=10
H=20

data4<-list()
data5<-list()
data6<-list()
data7<-list()
data8<-list()
data9<-list()
for (i in 1:1000){
  data3<-data2[sample(nrow(data2), 95, replace=TRUE), ]

year<-log(data3$year)
im<-log(data3$IM)
ex<-log(data3$EX)
dd<-log(data3$DD)
ddf<-log(data3$DD_F)
reer<-log(data3$REER)

Lim<-Lag(im,-1)
Lex<-Lag(ex, -1)

system<-list(IM=eqIM, EX=eqEX)
system2SLS<-systemfit(system,'2SLS',inst=inst)

coef<-as.data.frame(t(coef(system2SLS)))

Gamma=diag(K)
colnames(Gamma)=names.endo
rownames(Gamma)=names.endo
B=matrix(0,ncol=K,nrow=M)
colnames(B)=names.endo
rownames(B)=names.exo

Gamma["im","ex"]=-coef$IM_ex
Gamma["ex","im"]=-coef$EX_im
B['Intercept','ex']=-coef$`EX_(Intercept)`
B['Intercept','im']=-coef$`IM_(Intercept)`
B['dd','im']=-coef$IM_dd
B['ddf','ex']=-coef$EX_ddf
B['reer','ex']=-coef$EX_reer
B['reer','im']=-coef$IM_reer
B['Lex','ex']=-coef$EX_Lex
B['Lim','im']=-coef$IM_Lim


PI=-B%*%solve(Gamma)
PIX=PI[1:4,]

VARPHI=matrix(0,ncol=K,nrow=K)
colnames(VARPHI)=names.endo
rownames(VARPHI)=paste("L",names.endo,sep="")

VARPHI[rownames(PI)[5:6],]=PI[5:6,]

a<-PIX%*%(VARPHI%^%s)

LR=PIX%*%solve(diag(K)-VARPHI)

skumulowany<-LR%*%(diag(K)-VARPHI%^%s)

data4[i]<-list(a)
data5[i]<-list(LR)
data6[i]<-list(skumulowany)
}

df<-data.frame(data4)
df<-data.frame(t(df))
df1 <- split(df, rep(1:2, length = length(df)))
df1 <- as.data.frame(df1)
df1 <- na.omit(df1)

df<-data.frame(data5)
df<-data.frame(t(df))
df2 <- split(df, rep(1:2, length = length(df)))
df2 <- as.data.frame(df2)
df2 <- na.omit(df2)

df<-data.frame(data6)
df<-data.frame(t(df))
df3 <- split(df, rep(1:2, length = length(df)))
df3 <- as.data.frame(df3)
df3 <- na.omit(df3)

mean1.1.1<-mean(df1$X1.Intercept)
mean1.1.2<-mean(df1$X1.dd)
mean1.1.3<-mean(df1$X1.ddf)
mean1.1.4<-mean(df1$X1.reer)
mean1.2.1<-mean(df1$X2.Intercept)
mean1.2.2<-mean(df1$X2.dd)
mean1.2.3<-mean(df1$X2.ddf)
mean1.2.4<-mean(df1$X2.reer)

mean2.1.1<-mean(df2$X1.Intercept)
mean2.1.2<-mean(df2$X1.dd)
mean2.1.3<-mean(df2$X1.ddf)
mean2.1.4<-mean(df2$X1.reer)
mean2.2.1<-mean(df2$X2.Intercept)
mean2.2.2<-mean(df2$X2.dd)
mean2.2.3<-mean(df2$X2.ddf)
mean2.2.4<-mean(df2$X2.reer)

mean3.1.1<-mean(df3$X1.Intercept)
mean3.1.2<-mean(df3$X1.dd)
mean3.1.3<-mean(df3$X1.ddf)
mean3.1.4<-mean(df3$X1.reer)
mean3.2.1<-mean(df3$X2.Intercept)
mean3.2.2<-mean(df3$X2.dd)
mean3.2.3<-mean(df3$X2.ddf)
mean3.2.4<-mean(df3$X2.reer)

q<-c(mean1.1.1,mean1.1.2,mean1.1.3,mean1.1.4)
qq<-c(mean1.2.1,mean1.2.2,mean1.2.3,mean1.2.4)
qqq<-c(mean2.1.1,mean2.1.2,mean2.1.3,mean2.1.4)
w<-c(mean2.2.1,mean2.2.2,mean2.2.3,mean2.2.4)
ww<-c(mean3.1.1,mean3.1.2,mean3.1.3,mean3.1.4)
www<-c(mean3.2.1,mean3.2.2,mean3.2.3,mean3.2.4)

oszacowania<-as.data.frame(rbind(q,qq,qqq,w,ww,www))
colnames(oszacowania)=c('Intercept', 'dd','ddf','reer')
rownames(oszacowania)=c('dynamiczny import','dynamiczny export','dlugookresowy import','dlugookresowy export','skumulowany import','skumulowany export')
oszacowania

library(Rmisc)
CI1.1.1<-CI(df1$X1.Intercept, ci=0.95)
CI1.1.2<-CI(df1$X1.dd, ci=0.95)
CI1.1.3<-CI(df1$X1.ddf, ci=0.95)
CI1.1.4<-CI(df1$X1.reer, ci=0.95)
CI1.2.1<-CI(df1$X2.Intercept, ci=0.95)
CI1.2.2<-CI(df1$X2.dd, ci=0.95)
CI1.2.3<-CI(df1$X2.ddf, ci=0.95)
CI1.2.4<-CI(df1$X2.reer, ci=0.95)

CI2.1.1<-CI(df2$X1.Intercept, ci=0.95)
CI2.1.2<-CI(df2$X1.dd, ci=0.95)
CI2.1.3<-CI(df2$X1.ddf, ci=0.95)
CI2.1.4<-CI(df2$X1.reer, ci=0.95)
CI2.2.1<-CI(df2$X2.Intercept, ci=0.95)
CI2.2.2<-CI(df2$X2.dd, ci=0.95)
CI2.2.3<-CI(df2$X2.ddf, ci=0.95)
CI2.2.4<-CI(df2$X2.reer, ci=0.95)

CI3.1.1<-CI(df3$X1.Intercept, ci=0.95)
CI3.1.2<-CI(df3$X1.dd, ci=0.95)
CI3.1.3<-CI(df3$X1.ddf, ci=0.95)
CI3.1.4<-CI(df3$X1.reer, ci=0.95)
CI3.2.1<-CI(df3$X2.Intercept, ci=0.95)
CI3.2.2<-CI(df3$X2.dd, ci=0.95)
CI3.2.3<-CI(df3$X2.ddf, ci=0.95)
CI3.2.4<-CI(df3$X2.reer, ci=0.95)

q<-rbind(CI1.1.1,CI1.1.2,CI1.1.3,CI1.1.4)
qq<-rbind(CI1.2.1,CI1.2.2,CI1.2.3,CI1.2.4)
qqq<-rbind(CI2.1.1,CI2.1.2,CI2.1.3,CI2.1.4)
w<-rbind(CI2.2.1,CI2.2.2,CI2.2.3,CI2.2.4)
ww<-rbind(CI3.1.1,CI3.1.2,CI3.1.3,CI3.1.4)
www<-rbind(CI3.2.1,CI3.2.2,CI3.2.3,CI3.2.4)

przedzialyufnosci<-cbind(q,qq,qqq,w,ww,www)
przedzialyufnosci<-as.data.frame(t(przedzialyufnosci))
colnames(przedzialyufnosci)=c('Intercept', 'dd','ddf','reer')
rownames(przedzialyufnosci)=c('gorny DI','sredni DI', 'dolny DI', 'gorny DE', 'sredni DE','dolny DE','gorny DOI','sredni DOI','dolny DOI','gorny DOE', 'sredni DOE','dolny DOE','gorny SI', 'sredni SI', 'dolny SI', 'gorny SE', 'sredni SE', 'dolny SE')
przedzialyufnosci

gornyDII<-przedzialyufnosci[1,1]
dolnyDII<-przedzialyufnosci[3,1]
gornyDIDD<-przedzialyufnosci[1,2]
dolnyDIDD<-przedzialyufnosci[3,2]
gornyDIDDF<-przedzialyufnosci[1,3]
dolnyDIDDF<-przedzialyufnosci[3,3]
gornyDIREER<-przedzialyufnosci[1,4]
dolnyDIREER<-przedzialyufnosci[3,4]
gornyDEI<-przedzialyufnosci[4,1]
dolnyDEI<-przedzialyufnosci[6,1]
gornyDEDD<-przedzialyufnosci[4,2]
dolnyDEDD<-przedzialyufnosci[6,2]
gornyDEDDF<-przedzialyufnosci[4,3]
dolnyDEDDF<-przedzialyufnosci[6,3]
gornyDEREER<-przedzialyufnosci[4,4]
dolnyDEREER<-przedzialyufnosci[6,4]
gornyDOI<-przedzialyufnosci[7,1]
dolnyDOI<-przedzialyufnosci[9,1]
gornyDOIDD<-przedzialyufnosci[7,2]
dolnyDOIDD<-przedzialyufnosci[9,2]
gornyDOIDDF<-przedzialyufnosci[7,3]
dolnyDOIDDF<-przedzialyufnosci[9,3]
gornyDOIREER<-przedzialyufnosci[7,4]
dolnyDOIREER<-przedzialyufnosci[9,4]
gornyDOEI<-przedzialyufnosci[10,1]
dolnyDOEI<-przedzialyufnosci[12,1]
gornyDOEDD<-przedzialyufnosci[10,2]
dolnyDOEDD<-przedzialyufnosci[12,2]
gornyDOEDDF<-przedzialyufnosci[10,3]
dolnyDOEDDF<-przedzialyufnosci[12,3]
gornyDOEREER<-przedzialyufnosci[10,4]
dolnyDOEREER<-przedzialyufnosci[12,4]
gornySII<-przedzialyufnosci[13,1]
dolnySII<-przedzialyufnosci[15,1]
gornySIDD<-przedzialyufnosci[13,2]
dolnySIDD<-przedzialyufnosci[15,2]
gornySIDDF<-przedzialyufnosci[13,3]
dolnySIDDF<-przedzialyufnosci[15,3]
gornySIREER<-przedzialyufnosci[13,4]
dolnySIREER<-przedzialyufnosci[15,4]
gornySEI<-przedzialyufnosci[16,1]
dolnySEI<-przedzialyufnosci[18,1]
gornySEDD<-przedzialyufnosci[16,2]
dolnySEDD<-przedzialyufnosci[18,2]
gornySEDDF<-przedzialyufnosci[16,3]
dolnySEDDF<-przedzialyufnosci[18,3]
gornySEREER<-przedzialyufnosci[16,4]
dolnySEREER<-przedzialyufnosci[18,4]

dolnyDII<df1$X1.Intercept<gornyDII
dolnyDII<df3$X1.Intercept & df3$X1.Intercept<gornyDII
#i tak dalej...
