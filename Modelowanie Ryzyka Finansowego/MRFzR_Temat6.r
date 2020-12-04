rm(list=ls())
require(zoo)
require(ggplot2)
require(rugarch)
require(knitr)


source("LoadFundData.R")
#source("LoadWigData.R")


# Panel wykresow: pelna proba 
par(mfrow=c(2,1), cex = 0.7, bty="l")
plot(P, main="Poziom", xlab="", ylab="")
plot(r, main="stopy zwrotu", xlab="", ylab="")

##########################################################
# Ustawienia do testow warunkow skrajnych                #
##########################################################
 
# W przypadku testow warunkow skrajnych zazwyczaj wykorzystujemy jak najd³u¿sz¹ probe
T     <- length(r)
N     <- T                               
r     <- tail(r,N)                       
R     <- coredata(r)
T
p1<- P[2853]
p2<- P[2789]
p2
pc<-((p2-p1)/p2)*100
pc
##############################################################
# 1. Scenariusz historyczny                                  #
# Najgorsza stopa zwrotu w oknie w_length dni                #
##############################################################

# szerokosc okna
w_length = 250             

# rolowana minimum w ci¹gu nastêpnych w_length dni
MAmin     <- rollapply(P, width=w_length, min, by=1, align="left")
MAdecline <- MAmin/P[1:length(MAmin)]-1            # minimalna stopa zwrotu w ciagu kolejnych w_length dni

N0     <- which.min(MAdecline)
MAdecline[N0]

plot(MAdecline)
abline(v=index(MAdecline[N0]),col=2, lwd=3)

ScHist <- min(MAdecline)          # maksymalny spadek w ciagu roku

# data poczatkowa i koÅ„cowa spadku
N1        <- N0 + which.min(P[N0:(N0+w_length-1)]) - 1
P[c(N0,N1)]

##############################################################
# 2A. Standardowy scenariusz hipotetyczny                     #
##############################################################

# Sklad portfela
A1 = 0.01       # akcje krajowe
A2 = 0.03       # akcje zagraniczne
A3 = 0.05       # obligacje skarbowe, 5-letnie
A4 = 0.05       # obligacje korporacyjne, 3 -letnie
A5 = 0.86       # surowce
A = c(A1,A2,A3,A4,A5)
A
# Etap 1: analiza wrazliwosci wzgledem czynnikow ryzyka

RF1 = sum(c(0,0,0,0,1)*A)    #spadek cen surowców o 50%
# procentowa reakcja na 1% wzrost indeksow gieldowych (w kraju i zagranica) 
RF2 = sum(c(0,0,0,0,-0.5 )*A)   #wzroœcie cen z³ota o 20%
# procentowa reakcja na wzrost krzywej dochodowosci na calej dlugosci o 100 pb.
RF3 = sum(c(1.5,1.5,0,0,0 )*A)   #spadek cen akcji na rynkach rozwijajacych siê i rozwiniêtych o 15%
# procentowa reakcja na deprecjacje kursu o 1%
RF4 = sum(c(0,0,5,3,0)*A)   #przesuniêciu krzywej dochodowoœci o 50 pb w dó³
# procentowa reakcja na wzrost spreadu kredytowego o 100 pb.
RF5 = sum(c(0,0.5,0,0,0.5)*A)   # odp³ywie kapita³u z rynków wschodz¹cych i deprecjacji ich kursów o 10%
# procentowa reakcja na wzrost cen surowcow o 1%

# Etap 2: analizy scenariuszowe 
S1 = -abs(RF1*10); # increase and decrease in equity index values by 10 percent
S2 = -abs(RF2);    # parallel yield curve shifts of 100 basis points up and down
S3 = -abs(RF3*20);  # increase and decrease in the exchange value of foreign currencies by 20 percent 
S4 = -abs(RF4);     # increase and decrease in cedit spreads by 100 basis points.

ScHipA <- S1+S2+S3+S4;
ScHipA                       # hipotetyczny scenariusz skrajny       

##############################################################
# 2B. Unikatowy scenariusz hipotetyczny                      #
##############################################################

# Przyjmijmy, ze w scenariuszu nastapi zalamanie wzrostu w Chinach i wzrost globalnej awersji do ryzyka:
# 1. ucieczka kapitalu z rynkow wschodzacych, w tym PL: deprecjacja PLN o 10%
# 2. ucieczka kapitalu z papierow udzialowych: spadki indeksow gieldowych o 15%
# 3. naplyw kaptalu na rynek obligacji: spadek krzywej doch. o 100 pb w krajach rozwinietych
# 4. ale takze wzrost ryzyka kredytowego dla krajow rozwijajacych o 50 pb
# 5. spadek cen surowcow o 25%

# Wplyw na portfel
ScHipB <- RF1*(-50) + RF2*(20) + RF3*(-15) + RF4*-(0.5) + RF5*(-10)
ScHipB

#################################################################
# 3. Wartosc zagrozona w sytuacji skrajnej: stressed VaR, S-VaR #
#################################################################

p = 0.01   # poziom tolerancji
H = 250    # horyzont
w_length  <- 250            # okno dla rolowanej (stresowanej) statystyki
q         <- 0.05          # kwantyl dla (stresowanej) statystyki     

# przyjmijmy, ze r~N(mu,sig^2)

# 3.0. standardowy VaR/ES  
N    <- 1250               # liczba wykorzystanych obserwacji
R    <- tail(R,N)
# r    <- tail(r,N)

# por. Temat 5
m0   <- mean(R)*H
s0   <- sd(R)*sqrt(H)
S_VaR0    <- m0 + qnorm(p)*s0            
S_ES0     <- m0 - s0*dnorm(qnorm(p))/p   

# 3.1. zaburzamy wartosc sredniej stopy zwrotu
MAmean    <- rollapply(r, width=w_length, mean, by=1, align="right")
m1        <- quantile(MAmean,q, names=FALSE)*H
s1        <- s0
kable(c(m0,m1),digits=3, col.names = "stresses mean")
S_VaR1    <- m1 + qnorm(p)*s1            # stressed VaR 
S_ES1     <- m1 - s1*dnorm(qnorm(p))/p   # stressed Expected shortfall

# 3.2. zaburzamy wartosc dla odchylenia standardowego
MAsd      <- rollapply(r, width=w_length, sd, by=1, align="right")
m2        <- m0
s2        <- quantile(MAsd,1-q, names=FALSE)*sqrt(H)
kable(c(s0,s2),digits=3, col.names = "stresses std.")
S_VaR2    <- m2 + qnorm(p)*s2            # stressed VaR 
S_ES2     <- m2 - s2*dnorm(qnorm(p))/p   # stressed Expected shortfall

# 3.3. polaczenie
m3        <- m1
s3        <- s2
S_VaR3    <- m3 + qnorm(p)*s3            # stressed VaR 
S_ES3     <- m3 - s3*dnorm(qnorm(p))/p   # stressed Expected shortfall

kable(100*c(S_VaR0,S_VaR1,S_VaR2,S_VaR3), col.names = "porownanie SVaR")
kable(100*exp(c(S_VaR0,S_VaR1,S_VaR2,S_VaR3)), col.names = "Wartosc 100PLN")

