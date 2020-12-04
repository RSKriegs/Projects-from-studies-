rm(list=ls())
# Pakiet do zapisu danych w formacie zoo
require(zoo)
require(ggplot2)
require(rugarch)
require(MASS)
require(moments)

raw.data <- read.csv2('BPH023.mst', head=TRUE, sep=",", dec=".")
head(raw.data)
dateForm <- "%Y%m%d"

class(raw.data$"X.DTYYYYMMDD.")
class(raw.data$"X.OPEN.")
prices <- raw.data$"X.OPEN."
dates  <- as.Date(as.character(raw.data$"X.DTYYYYMMDD."), format=dateForm)
data   <- zoo(prices , order.by = dates)
head(data)
plot(data)

P       <- data
r       <- diff(log(P))

setwd(*C:/Users/Rafal/Desktop/R/materialy/funds)
# source("LoadWigData.R")

# Panel wykresow: pelna proba 
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(P, col='red', main="Poziom", xlab="", ylab="")
plot(r, main="stopy zwrotu", xlab="", ylab="")

# Panel wykresow: pelna proba 
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(P, main="Poziom", xlab="", ylab="")
plot(r, main="stopy zwrotu", xlab="", ylab="")


##########################################################
# Ustawienia do analizy                                  #
##########################################################

T  <- length(r)                        # liczba obserwacji dla pelnej proby
N  <- 1250                             # ustalamy probe na podstawie ktorej liczymy statystyki
r  <- tail(r,N)                        # stopy zwrotu
R  <- coredata(r)

##########################################################
# Statystyki opisowe dla logarytmicznych stóp zwrotu WIG #
##########################################################

mu = sum(R)/N

R0    <- R - mu
M2    <- sum(R0^2)/N
M3    <- sum(R0^3)/N
M4    <- sum(R0^4)/N

sig <- sqrt(M2)       ## zmiennosc
S   <- M3/(sig^3)     ## skosnosc
K   <- M4/(sig^4)     ## kurtoza

# to samo z pakietem moments
require(moments)
M1   = moment(R, order = 1, central = FALSE, na.rm = TRUE)
M2   = moment(R, order = 2, central = TRUE, na.rm = TRUE)
M3   = moment(R, order = 3, central = TRUE, na.rm = TRUE)
M4   = moment(R, order = 4, central = TRUE, na.rm = TRUE)
mu0  = M1
sig0 = sqrt(M2)
S0   = M3/(sig0^3)
K0   = M4/(sig0^4)

# Annualizacja
Nyear <- 250
muA   <- mu*Nyear
sigA  <- sig*sqrt(Nyear)


##################################
# Test na normalnosc stop zwrotu #
##################################

#D'Agostino test for skewness
agostino.test(R)

#Anscombe-Glynn test of kurtosis for normal samples
anscombe.test(R)

#Jarque-Bera test of normality
jarque.test(R)

#######################################################
# Rozk³ad t-Studenta                                  #
#######################################################

x = seq(-5,5,0.01)
junk <- cbind(x,dnorm(x),ddist("std",x,shape=10),ddist("std",x,shape=5),ddist("std",x,shape=3))
colnames(junk) <- c("x","vInf","v10","v5","v3")

ggplot(data = data.frame(junk), aes(x = x)) + 
  geom_line(aes(y = vInf,colour="v=Inf")) + 
  geom_line(aes(y = v10,colour="v=10")) +
  geom_line(aes(y = v5,colour="v=5")) +
  geom_line(aes(y = v3,colour="v=3")) +
  labs(title="rozk³ad t-Studenta", y="", x="", caption="")+
  theme_bw()+
  scale_colour_manual("", 
                      breaks = c("v=Inf", "v=10","v=5","v=3"),
                      values = c("black", "red", "blue", "green"))


#######################################################
# Empiryczna funkcja gestosci (kernel density plot)   #
#######################################################

R0        <- (R-mu)/sig
bwdth     <- 0.1

ggplot(data.frame(R0), aes(x = R0)) +
  theme_bw() +
  geom_histogram(binwidth = bwdth, colour = "white", fill = "yellow4", size = 0.1) +
  stat_function(fun = function(x) dnorm(x)*N*bwdth, color = "red", size = 1)                        # rozklad normalny
# + stat_function(fun = function(x) ddist("std",x,shape=5)*N*bwdth, color = "black", size = 1)    # rozklad t-Studenta
# + xlim(-7,-2)


#######################################################
# Wykres kwantyl-kwantyl (QQplot)                     #
#######################################################

R0     <- (R-mu)/sig

# Tabela kwantyli
q              <- seq(0.01, 0.20, 0.01)
QQtable        <- data.frame(quantile(R0,q),qnorm(q))   # rozklad normalny
names(QQtable) <- c("kwantyl empiryczny", "kwantyl teoretyczny")

# QQ plot
q            <- seq(0.001, 0.999, 0.001)
Qemp         <- qnorm(q)    
#Qemp         <- quantile(rnorm(N),q)
#Qemp         <- quantile(R0,q)                    # kwantyl empiryczny
Qteo         <- qnorm(q)                          # kwantyl teoretyczny (rozklad normalny)
#v=5;   Qteo  <- qt(q,v)*sqrt((v-2)/v)             # rozklad t-Studenta (wariancja to v/v-2)
lim0    <- c(-5,5)                           # zakres na wykresie
par(mfrow=c(1,1), cex = 0.7, bty="l")
plot(Qemp,Qteo, main="QQplot", col="red", xlim = lim0, ylim = lim0,
     xlab="kwantyl empiryczny", ylab="kwantyl teoretyczny") # plots the results
abline(a=0,b=1, lwd=2)


#######################################################
# Estymacja parametrow  rozkladu t-Studenta           #
#######################################################

# Metoda momentow (K - kurtoza)
v0 <- 4 + 6/(K-3)

# Metoda najwiêkszej wiarygodnoœci
require(MASS)

d0 <- fitdistr(R0, "normal")
d0$loglik

d1 <- fitdistr(R0, "t", m = 0, start = list(s=sqrt((v0-2)/v0), df=v0), lower=c(0.001,3))
d1$loglik
v=d1$estimate[[2]]

# Comparison
(d1$loglik-d0$loglik)/N

# QQ plot
q        <- seq(0.001, 0.999, 0.001)
Qemp     <- quantile(R0,q)                     # kwantyl empiryczny
Qteo     <- qt(q,v)*sqrt((v-2)/v)              # rozklad t-Studenta (wariancja to v/v-2)
lim0    <- c(-5,5)                             # zakres na wykresie
par(mfrow=c(1,1), cex = 0.7, bty="l")
plot(Qemp,Qteo, main="QQplot", col="red", xlim = lim0, ylim = lim0,
     xlab="kwantyl empiryczny", ylab="kwantyl teoretyczny") # plots the results
abline(a=0,b=1, lwd=2)

#######################################################
# Autokorelacje, czyli zaleznosci w czasie            #
#######################################################

require(tseries)

# ACF dla stop zwrotu
z0 <- acf(R,20, plot=TRUE)
plot(z0$acf[2:20], type="h", main="ACF dla stop zwrotu", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))

# ACF dla kwadratow stop zwrotu
z1 <- acf(R^2,20, plot=FALSE)
plot(z1$acf[2:20], type="h", main="ACF dla kwadratow stop zwrotu", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.4), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))

# Kros-korelacje dla stop zwrotu i kwaratow stop zwrotu
z2 <- ccf(R^2, R, 20, type="correlation", plot=TRUE)
plot(z2$acf[22:40], type="h", main="korelacja kwadratow stop zwrotu i opoznionych st. zwrotu", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))

plot(z2$acf[20:1], type="h", main="korelacja st. zwrotu i opoznionych kwadratow st. zwrotu", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))

# panel wykresow
par(mfrow=c(2,2), cex = 0.7, bty="l")
plot(z0$acf[2:20], type="h", main="ACF dla stop zwrotu", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z1$acf[2:20], type="h", main="ACF dla kwadratow stop zwrotu", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.4), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z2$acf[22:40], type="h", main="korelacja kw. st. zwrotu i op. st. zw.", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))
plot(z2$acf[20:1], type="h", main="korelacja st. zw. i op. kw. st. zw.", xlab="Lag", ylab="ACF", ylim=c(-0.2,0.2), las=1)
abline(h=c(2/sqrt(length(R)),0,-2/sqrt(length(R))),lty=c(2,1,2))


# testy Ljunga-Boxa

Box.test(R, lag = 20, type = c("Ljung-Box"))
Box.test(R^2, lag = 20, type = c("Ljung-Box"))


##########################################################
# Ustawienia do modelu VaR                               #
##########################################################

T  <- length(r)                          # liczba obserwacji
p  <- 0.05                               # poziom tolerancji VaR
H  <- 1                                  # horyont

N  <- 1250                               # ustalamy probe na podstawie ktorej liczymy VaR
r  <- tail(r,N)                          # stopy zwrotu
R  <- coredata(r)

##########################################################
# Symulacja historyczna                                  #
##########################################################

# wartosc zagrozona VaR
R0       <- sort(R)              
N0       <- floor(N*p)                  # obserwacja dla p-tego kwintyla                   
VaR_HS  <- R0[N0]                      # porownaj z: quantile(R,p)
# Oczekiwana strata ES
ES_HS  <- mean(R0[1:N0])               # srednia wartosc stop zwrotu ponizej VaR

# wykres na histogramie
hist(R, main="Symulacja historyczna", xlab="stopa zwrotu", ylab="l. obserwacji", nclass = 500, col="black")
abline(v=c(VaR_HS,ES_HS) ,lwd=2, col=c("blue","red"),lty=c(2,2))
legend("left", c("VaR dla sym. hist.","ES dla sym. hist."), lty=c(2,2), col=c("blue","red"), bty="n")

# Przedzial ufnosci dla VaR 
d1   <- density(R)
x    <- d1$x                              # punkty obserwacji
temp <- tail(which(x<VaR_HS),1)
fhat <- d1$y[temp]

seVaRHS <- sqrt(p*(1-p)/(N*fhat^2))                        # odchylenie standardowe dla kwintyla
CI      <- c(VaR_HS-1.96*seVaRHS, VaR_HS+1.96*seVaRHS)   # 95% przedzial ufnosci
CI

##########################################################
# Model parametryczny                                    #
##########################################################

# rozklad normalny
m   <- mean(R) 
s   <- sd(R)

VaR_N     <- qnorm(p)*s + m            # VaR dla stopy zwrotu
ES_N      <- m - s*dnorm(qnorm(p))/p   # Expected shortfall 

# wykres na density plot
d2  <- density(rnorm(100000, mean=m,sd=s))
plot(d2, main="Rozklad normalny", xlab="stopa zwrotu", ylab="f. gestosci", col="black")
abline(v=c(VaR_N,ES_N) ,lwd=2, col=c("blue","red"),lty=c(2,2))
legend("left", c("VaR dla rozkl. norm.","ES dla rozkl. norm."), lty=c(2,2), col=c("blue","red"), bty="n")

##########################################################
# Symulacje Monete Carlo dla rozk³adu normalnego         #
##########################################################

M      <- 10000     # liczba losowan w symulacjach Monte Carlo
m      <- mean(R) 
s      <- sd(R)

draws0    <- rnorm(M, mean=m, sd=s)
draws1    <- sort(draws0)                 # uporzadkowane obserwacje
M0        <- floor(M*p)                   # obserwacja dla p-tego kwintyla                   
VaR_N0    <- draws1[M0]                   # porownaj z: quantile(draws1,p)

# Oczekiwana strata ES
ES_N0     <- mean(draws1[1:M0])           # srednia wartosc

# Porownanie uzyskanych wartosci
c(VaR_N,VaR_N0)
c(ES_N,ES_N0)

#########################
# Grube ogony           #
#########################

# Rozklad t-Studenta
# standaryzacja
m  <- mean(R)
s  <- sd(R)
R0 <- (R-m)/s

# estymacja liczby stopni swobody
K     <- kurtosis(R0)
v0    <- 4+6/(K-3)         # metoda momentow

dt    <- fitdistr(R0, "t", m = 0, start = list(s=sqrt((v0-2)/v0), df=v0), lower=c(0.001,3))
v1    <- dt$estimate[[2]]
# liczenie VaR
VaR_t <- m + s*qt(p,v1) * sqrt((v1-2)/v1)           
# Oczekiwana strata
qf       <- function(x) qt(x, v1)
ES_t     <- m + s*sqrt((v1-2)/v1)*(1/p * integrate(qf, 0, p)$value)   

# Uwaga. Alternatywa to skorzystanie z pakietu ugarch
# require(rugarch)
# temp  <- fitdist("std", R0)
# v1    <- temp$pars["shape"]
# VaR_t <- m + s*qdist("std",shape=v1)
# qf    <- function(x) qdist("std", p=x, shape=v1)
# ES_t  <- m + s*(1/p * integrate(qf, 0, p)$value) 

draws <- rt(100000,v1)*sqrt((v1-2)/v1)*s + m
d3    <- density(draws)
plot(d3, main="Rozklad t-Studenta", xlab="stopa zwrotu", ylab="f. gestosci", xlim=c(min(R),max(R)),col="black")
abline(v=c(VaR_t,ES_t) ,lwd=2, col=c("blue","red"),lty=c(2,2))
legend("left", c("VaR dla rozkl. t","ES dla rozkl. t"), lty=c(2,2), col=c("blue","red"), bty="n")

# Rozszerzenie Cornisha-Fischera    #
#####################################

require(moments)
m   = mean(R)
s   = sd(R)
S   = skewness(R)
K   = kurtosis(R)

PsiP    <- qnorm(p)
VaR_CF  <- m + s*(PsiP + (PsiP^2-1)/6*S + (PsiP^3-3*PsiP)/24*(K-3) - (2*PsiP^3-5*PsiP)/36*(S^2) ) 

# Porównanie 

# HS / NORM / t-Student na jednym wykresie #

plot(d1, main="Funkcje gestosci", xlab="stopa zwrotu", ylab="f. gestosci", xlim=c(min(R),max(R)),col="black", lwd=2)
lines(d2,col="red", lwd=2)
lines(d3,col="green", lwd=2)
abline(v=c(VaR_HS,VaR_N,VaR_t) ,lwd=2, col=c("black", "red","green"),lty=c(2,2,2))
abline(v=c(ES_HS,ES_N,ES_t) ,lwd=2, col=c("black", "red","green"),lty=c(3,3,3))
legend("left", c("sym. historyczna","rozklad normalny","rozklad t"), lty=c(1,1,1), col=c("black", "red","green"), bty="n")

ES_CF <- m+s*(1/p*integrate(qf,0,p)$value)

Tabela           <- matrix(c(VaR_HS,VaR_N,VaR_t,VaR_CF, ES_HS,ES_N,ES_t,NaN),4,2)
colnames(Tabela) <- c("VaR","ES")
rownames(Tabela) <- c("sym. hist.", "rozk. norm.", "rozklad t", "Cornish-Fischer")
round(-100*Tabela,2)

# Funkcje i wykres VaR dla ró¿nych p  #
#######################################
source("MRFzR_FunkcjeBlok1.R")
TabelaP = matrix(NA,100,7)
for(p in seq(0.001,0.1,0.001)){
  TabelaP[p*1000,] <- c(p,VaRhist(R,p)$VaR, VaRnorm(R,p)$VaR, VaRt(R,p,v1)$VaR, VaRCF(R,p)$VaR, VaREWMA(R,p,v,lambda)$VaR, VaRIGARCH(R,p,v,lambda)$VaR)
}

plot(TabelaP[,1],TabelaP[,2], main="Poziom tolerancji a VaR", 
     xlab="poziom tolerancji (p)", ylab="VaR", 
     ylim=c(min(TabelaP[,2:7]),max(TabelaP[,2:7])), xlim=c(0,0.1),col="black", lwd=2, type="l")
lines(TabelaP[,1],TabelaP[,3],col="red", lwd=2, type="l")
lines(TabelaP[,1],TabelaP[,4],col="green", lwd=2, type="l")
lines(TabelaP[,1],TabelaP[,5],col="orange", lwd=2, type="l")
lines(TabelaP[,1],TabelaP[,6],col="blue", lwd=2, type="l")
legend("right", c("sym. historyczna","rozklad normalny","rozklad t","Cornish-Fischer", 'EWMA','GARCH'), lty=c(1,1,1,1), col=c("black", "red","green","orange"), bty="n")
abline(v=c(0.05) ,lwd=2, lty=2)

##########################################################
# Ustawienia do modelu VaR                               #
##########################################################

T  <- length(r)                          # liczba obserwacji
p  <- 0.05                               # poziom tolerancji VaR
H  <- 1                                  # horyont

N     <- 1250                             # ustalamy probe na podstawie ktorej liczymy VaR (5 lat)
r     <- tail(r,N)                        # dane
R     <- coredata(r)

############################
# wspolzaleznosci w czasie #
############################
require(forecast)

# brak autokorelacja zwrotow
Acf(R)
Pacf(R)

# ale jest grupowanie zmiennosci
Acf(R*R)
Pacf(R*R)

# oraz efekt dzwigni
Ccf(R*R, R)

########################################################################
# Grupowanie wariancji a VaR / ES                                      #
# Dla uproszczenia zakladamy rozklad t studenta o 5 stopniach swobody  #
########################################################################

require(rugarch)

# Ustawienia, czy obliczenia dla rozk³adu normalnego czy t-Studenta z 5 stopniami swobody
q    <- qdist("std", p=p, shape=5)
qf   <- function(x) qdist("std", p=x, shape=5)
# UWAGA: w pakiecie rugarch automatyczne losowanie z rokladu o parametrach mu=0, sig=1. Porownaj z:
# q    <- qt(p,5)*sqrt((5-2)/5)
# q    <- qdist("norm",p)
# qf   <- function(x) qdist("norm", p=x)


# 1. Staly rozklad
# t-Student o 5 stopniach swobody / rozklad normalny
#######################################################

mu    <- mean(r)
sigma <- sd(r)

# VaR i ES (jako calka numeryczna)
VaR_const   <- mu + sigma*q
ES_const <- mu + sigma*(1/p * integrate(qf, 0, p)$value)   

# Wykresy
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(merge( r, mu + 2*sigma, mu - 2*sigma), plot.type="single", col=c(1,2,2),
     main="log-zwroty vs +-2sd", ylab="" )
plot(merge(r, VaR_const, ES_const), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR i ES z modelu o stalych param.",  sep=""), ylab="", xlab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# 2. Parametry zmienne w czasie: model sredniej ruchomej (MA) 
# rozk³ad t-Studenta z 5 stopniami swobody / rozk³ad normalny
###############################################################

# szerokosc okna
w_length = 20             

# rolowana srednia i odchylenie standardowe
MAmean  <- rollapply(r, width=w_length, mean, by=1, align="right")
MAstd   <- rollapply(r, width=w_length, sd,   by=1, align="right")
# opznione o jeden okres, aby przedzial ufnosci na kolejny okres by³ oparty o dane do momentu t-1
MAmean <- lag(MAmean, -1)
MAstd  <- lag(MAstd,  -1)
# VaR i ES (rozklad t-Studenta o 5 stopniach swobody)
MAvar <- MAmean + MAstd*q 
MAes  <- MAmean + MAstd*(1/p * integrate(qf, 0, p)$value)

# wykresy
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(merge( r, MAmean + 2*MAstd, MAmean - 2*MAstd), plot.type="single", col=c(1,2,2),
     main="log-zwroty vs +-2sd", ylab="" )
plot(merge( r, MAvar, MAes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR i ES z modelu MA",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# obliczenia VaR i ES dla ostatniego okresu
mT      <- mean(coredata(tail(r,w_length)))
sT      <- sd(coredata(tail(r,w_length)))
VaR_MA  <- mT + sT*q
ES_MA   <- mT + sT*(1/p * integrate(qf, 0, p)$value)

# 3. Kalibrowany EWMA (RiskMetrics) 
# rozk³ad t-Studenta z 5 stopniami swobody / rozk³ad normalny
###############################################################

lambda      <- 0.94                          # parametr wygladzajacy
EWMAsig2    <- rep(0,N)
EWMAsig2[1] <- var(R)                        # punkt startowy
for (t in 2:N){             
  EWMAsig2[t] = lambda * EWMAsig2[t-1]  + (1-lambda) * R[t-1]^2
}
EWMAstd      <- EWMAsig2^0.5
# VaR i ES 
EWMAvar <- EWMAstd*q 
EWMAes  <- EWMAstd*(1/p * integrate(qf, 0, p)$value)

# wykresy
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(merge( r, -2*EWMAstd, 2*EWMAstd), plot.type="single", col=c(1,2,2),
     main="log-zwroty vs EWMA +-2sd", ylab="" )
plot(merge( r, EWMAvar, EWMAes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR i ES z kalibrowanego modelu EWMA",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# obliczenia VaR i ES dla ostatniego okresu
mT  <- 0
s2T <- lambda * EWMAsig2[N]  + (1-lambda) * R[N]^2
sT  <- sqrt(s2T)
VaR_EWMA = mT + sT*q
ES_EWMA  = mT + sT*(1/p * integrate(qf, 0, p)$value)

# EWMA w pakiecie ugarch - wprowadzenie do modeli GARCH
###########################################################

lambda = 0.94

EWMAspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                       variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                       fixed.pars=list(alpha1=1-lambda, omega=0, shape=5), distribution.model = "std")
#                      fixed.pars=list(alpha1=1-lambda, omega=0),         distribution.model = "norm")

# przefiltrowanie danych
EWMAfit <- ugarchfit(data = r, spec = EWMAspec)
# nasze parametry
round(coef(EWMAfit),3)

# porownanie szeregow uzyskanych dwoma metodami
plot(merge(sigma(EWMAfit),EWMAstd), plot.type="single", col=c(1,2))
tail(merge(sigma(EWMAfit),EWMAstd))
plot(EWMAfit, which=9)

# 4. Estymacja parametru lambda w modelu EWMA: model IGARCH  
###########################################################

IGARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                         variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                         fixed.pars=list(omega=0,shape=5), start.pars = list(alpha1=0.06), distribution.model = "std")
#                      fixed.pars=list(omega=0), start.pars = list(alpha1=0.06), distribution.model = "norm")

IGARCHfit <- ugarchfit(data = r, spec = IGARCHspec, solver="nlminb")
round(coef(IGARCHfit),3)

# VaR i ES 
IGARCHvar <- sigma(IGARCHfit)*q 
IGARCHes  <- sigma(IGARCHfit)*(1/p * integrate(qf, 0, p)$value)

# wykresy
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(IGARCHfit, which=1)
plot(merge( r, IGARCHvar, IGARCHes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR i ES z estymowanego modelu EWMA (IGARCH)",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# Obliczenia VaR i ES dla ostatniego okresu
IGARCHfct   <- ugarchforecast(IGARCHfit,n.ahead = 1)
mT  <- as.numeric(fitted(IGARCHfct))
sT  <- as.numeric(sigma(IGARCHfct))

VaR_IGARCH  = mT + sT*q
ES_IGARCH   = mT + sT*(1/p * integrate(qf, 0, p)$value)

# 5. Estymacja wszystkich parametrow: GARCH(1,1)  
###########################################################

GARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                        variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                        fixed.pars=list(shape=5), distribution.model = "std")
#                        distribution.model = "norm")

GARCHfit <- ugarchfit(data = r, spec = GARCHspec, solver="hybrid")
round(coef(GARCHfit),5)

# VaR i ES (rozklad t-Studenta o 5 stopniach swobody)
GARCHvar   <- fitted(GARCHfit) + sigma(GARCHfit)*q 
GARCHes    <- fitted(GARCHfit) + sigma(GARCHfit)*(1/p * integrate(qf, 0, p)$value)

# wykresy
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(GARCHfit, which=1) 
plot(merge( r, GARCHvar, GARCHes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR i ES z modelu GARCH(1,1))",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# Obliczenia VaR i ES dla ostatniego okresu
GARCHfct   <- ugarchforecast(GARCHfit,n.ahead = 1)
mT  <- as.numeric(fitted(GARCHfct))
sT  <- as.numeric(sigma(GARCHfct))

VaR_GARCH  = mT + sT*q
ES_GARCH   = mT + sT*(1/p * integrate(qf, 0, p)$value)

# 6. Efekt d¿wigni: model eGARCH(1,1)  
###########################################################

EGARCHspec <- ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
                         variance.model=list(model="eGARCH", garchOrder=c(1,1)),
                         fixed.pars=list(shape=5), distribution.model = "std")
#                        distribution.model = "norm")

EGARCHfit <- ugarchfit(data = r, spec = EGARCHspec, solver="nlminb")
round(coef(EGARCHfit),5)

# VaR i ES (rozklad t-Studenta o 5 stopniach swobody)
EGARCHvar   <- fitted(EGARCHfit) + sigma(EGARCHfit)*q 
EGARCHes    <- fitted(EGARCHfit) + sigma(EGARCHfit)*(1/p * integrate(qf, 0, p)$value)

# wykresy
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(EGARCHfit, which=1)
plot(merge( r, EGARCHvar, EGARCHes), plot.type="single", col=c(1,2,3), main=paste(100*p,"% VaR i ES z modelu EGARCH(1,1))",  sep=""), ylab="" )
legend("bottomright", c("VaR", "ES"), lty=1, col=2:3)

# Obliczenia VaR i ES dla ostatniego okresu
EGARCHfct   <- ugarchforecast(EGARCHfit,n.ahead = 1)
mT  <- as.numeric(fitted(EGARCHfct))
sT  <- as.numeric(sigma(EGARCHfct))

VaR_EGARCH  = mT + sT*q
ES_EGARCH   = mT + sT*(1/p * integrate(qf, 0, p)$value)

# Wykresy dla zmiennosci sigma
##################################

windows()
labssd <- c("log-zwroty", "staly rozklad", "MA", "EWMA", "IGARCH", "GARCH", "eGARCH")
plot( merge(r,
            2*sigma,
            2*MAstd,
            2*EWMAstd,
            2*sigma(IGARCHfit),
            2*sigma(GARCHfit),
            2*sigma(EGARCHfit),
            -2*sigma,
            -2*MAstd,
            -2*EWMAstd,
            -2*sigma(IGARCHfit),
            -2*sigma(GARCHfit),
            -2*sigma(EGARCHfit)),
      facets=NULL,
      main="log-zwroty vs +-2sd",
      col= c(1:7, 2:7),
      plot.type="single",
      ylab="", xlab="")
legend("bottomright", labssd, lty=1, col=c(1:7, 2:7))


# Wykresy dla VaR 
###################
windows()
labsvar <- c("log-zwroty", "staly rozklad", "MA", "EWMA", "IGARCH", "GARCH", "eGARCH")
plot(merge(r,
           VaR_const,
           MAvar,
           EWMAvar,
           IGARCHvar,
           GARCHvar,
           EGARCHvar),
     main=paste("log-zwroty vs ", 100*p,"% VaR", sep=""),
     col=1:7,
     plot.type="single",
     ylab="", xlab="")
legend("bottomright", labsvar, lty=1, col=1:7)

# wszystko razem
-100*c(VaR_const,VaR_MA, VaR_EWMA, VaR_IGARCH, VaR_GARCH, VaR_EGARCH)
-100*c(ES_const,ES_MA, ES_EWMA, ES_IGARCH, ES_GARCH, ES_EGARCH)

# z wykorzystaniem funkcji w pliku MRFzR_FunkcjeBlok1.R
source("MRFzR_FunkcjeBlok1.R")
VaR_const1  <- VaRt(R,p,5)$VaR
VaR_MA1     <- VaRt(tail(R,w_length),p,5)$VaR
VaR_EWMA1   <- VaREWMA(r, p, v=5)$VaR  
VaR_IGARCH1 <- VaRIGARCH(r, p, v=5)$VaR
VaR_GARCH1  <- VaRGARCH(r, p, v=5)$VaR
VaR_EGARCH1 <- VaRGARCH(r, p, v=5, eGARCH=TRUE)$VaR


# Porownanie
-100*cbind(
  c(VaR_const,VaR_MA, VaR_EWMA, VaR_IGARCH, VaR_GARCH, VaR_EGARCH),
  c(VaR_const1,VaR_MA1, VaR_EWMA1, VaR_IGARCH1, VaR_GARCH1, VaR_EGARCH1))

c(M0,M1,M2,M3,M4)
