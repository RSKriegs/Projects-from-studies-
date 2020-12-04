#Usuniêcie starych zmiennych - dok³adnie jedno rozwi¹zanie
rm(list=ls())
require(lpSolveAPI)
require(zoo)

aktywo1<-"RANKPROGR.mst";n.aktywo1="Rank Progress"
aktywo2<-"ELKOP.mst";n.aktywo2="ELKOP"

benchmark <- "wig_nrchom_d.csv"
stopa.wolna <- "plbpln1m_m.csv" 

dane.aktywo1 <- read.csv2(file = aktywo1, head=TRUE, sep=",", dec=".")
data.format <- "%Y%m%d"
cena.otwarcia <- dane.aktywo1$"X.OPEN."
daty <- as.Date(as.character(dane.aktywo1$"X.DTYYYYMMDD."), format=data.format)
c.aktywo1 <- zoo(cena.otwarcia , order.by = daty)

dane.aktywo2 <- read.csv2(file = aktywo2, head=TRUE, sep=",", dec=".")
cena.otwarcia <- dane.aktywo2$"X.OPEN."
daty <- as.Date(as.character(dane.aktywo2$"X.DTYYYYMMDD."), format=data.format)
c.aktywo2 <- zoo(cena.otwarcia , order.by = daty)

dane.benchmark <- read.csv2(file = benchmark, head=TRUE, sep=",", dec=".")
data.format <- "%Y-%m-%d"
b.cena.otwarcia <- dane.benchmark$"Otwarcie"
daty <- as.Date(as.character(dane.benchmark$"Data"), format=data.format)
b.c <- zoo(b.cena.otwarcia , order.by = daty)

dane.stopa.wolna <- read.csv2(file = stopa.wolna, head=TRUE, sep=",", dec=".")
rf <- dane.stopa.wolna$"Otwarcie"
daty <- as.Date(as.character(dane.stopa.wolna$"Data"), format=data.format)
rf <- zoo(rf/1200 , order.by = daty) #dzielnie bo w procentach i skali rocznej

#Ustalenie daty pocz¹tku i koñca próby
data.start <- "2016-04-01"
data.koniec <-"2020-04-30" #"2020-02-28"

require(xts)
c1<-window(c.aktywo1, start=as.Date(data.start), end=as.Date(data.koniec)) 
c2<-window(c.aktywo2, start=as.Date(data.start), end=as.Date(data.koniec)) 
b.c1<-window(b.c, start=as.Date(data.start), end=as.Date(data.koniec)) 
rf<-window(rf, start=as.Date(data.start), end=as.Date(data.koniec))

#konwersja szeregu danych dziennych na dane miesiêczne
#aktywo
c1.m <- apply.monthly(c1, last) #first, mean, mojafunkcja
K<-49 # d³ugoœæ okresu próby, 6 lat 
c1.m <-tail(c1.m, n=K)
daty.miesiace<-time(c1.m)
c2.m <- apply.monthly(c2, last) #first, mean, mojafunkcja
c2.m <-tail(c2.m, n=K)
#benchmark
b.c1.m <- apply.monthly(b.c1, last) #first, mean, mojafunkcja
b.c1.m <-tail(b.c1.m, n=K)
#stopa wolna od ryzyka
rf <- apply.monthly(rf, last) #first, mean, mojafunkcja
rf <-tail(rf, n=K)
rf<-rf[-c(length(rf))]
rf <- zoo(rf, order.by = daty.miesiace[-1]) 

#logarytmiczna stopa zwrotu i jej charakterystyki
r1 <- log(as.vector(c1.m[-1]))-log(as.vector(c1.m[-K])) 
r1 <- zoo(r1 , order.by = daty.miesiace[-1])
r2 <- log(as.vector(c2.m[-1]))-log(as.vector(c2.m[-K])) 
r2 <- zoo(r2 , order.by = daty.miesiace[-1])

b.r <- log(as.vector(b.c1.m[-1]))-log(as.vector(b.c1.m[-K])) 
b.r <- zoo(b.r , order.by = daty.miesiace[-1])

c(mean(r1),mean(r2))

#portfel o minimalnym ryzyku z³o¿ony z akcji 1 i 2 (wagi: w1, w2)
portfel.markovitz<-function(r1,r2){

  T<-length(r1)
  r1.sr<-mean(r1)
  r2.sr<-mean(r2)
  sigma1<-sd(r1)
  sigma2<-sd(r2)
  ro12<-cor(r1,r2)

#portfel Markowitza (minimalne odchylenie standardowe); wagi: w1, w2
C<-matrix(c(2*sigma1^2,2*sigma1*sigma2*ro12,1,
            2*sigma2*sigma1*ro12,2*sigma2^2,1,
            1,1,0), byrow = TRUE, nrow = 3)
D<-solve(C)
w1<-D[1,3]
w2<-D[2,3]
return(c(w1,w2))
}

###################    portfel z³o¿ony z akcji 1 i akcji 2 (wagi: w5 i w6), 
###################    który minimalizuje  odchylenie przeciêtne MAD (Mean Average Deviation)
portfel.MAD<-function(r1,r2){
  T<-length(r1)
  r1.sr<-mean(r1)
  r2.sr<-mean(r2)
  sigma1<-sd(r1)
  sigma2<-sd(r2)
  ro12<-cor(r1,r2)
r0=-1
ZPL <- make.lp(nrow=0, ncol=(2*T+3))
lp.control(ZPL, sense="max",simplextype="primal")
#funkcja celu
set.objfn(ZPL, rep(-1, T), indices = (T+1):(2*T)) #y1...yT,d1...dT, \mu, w5,w6
#warunki ograniczaj¹ce
for (i in 1:T){add.constraint(ZPL, c(1,-1,-1), "<=", 0, indices=c(i,(T+i),(2*T+1)))} 
for (i in 1:T){add.constraint(ZPL, c(-1,-1,1), "<=", 0, indices=c(i,(T+i),(2*T+1)))} 
for (i in 1:T){add.constraint(ZPL, c(1,-exp(r1[i])+1,-exp(r2[i])+1), "=", 0, indices=c(i,(2*T+2),(2*T+3)))} 
add.constraint(ZPL, c(-1,exp(r1.sr)-1,exp(r2.sr)-1), "=", 0, indices=c((2*T+1),(2*T+2),(2*T+3)))
add.constraint(ZPL, c(1), ">=", r0, indices=c((2*T+1)))
#czy dla ka¿dej wartoœci r0 uda siê wyznaczyæ portfel przy za³o¿eniu braku krótkiej sprzeda¿y?
add.constraint(ZPL, c(1,1), "=", 1, indices=c((2*T+2),(2*T+3)))
set.bounds(ZPL, columns=c(2*T+2, 2*T+3), lower=c(0,0), upper=c(1,1))
set.bounds(ZPL, columns=1:T, lower=rep(-1, T))
set.bounds(ZPL, columns=2*T+1, lower=r0)
#gdy dopuœcimy krótk¹ sprzeda¿, to nale¿y do³aczyæ:
#set.bounds(ZPL, columns=2*T+1, lower=-1)
#set.bounds(ZPL, columns=c(2*T+2, 2*T+3), lower=c(-Inf, -Inf))
ZPL
write.lp(ZPL, "ZPL.txt")
#rozwi¹zujemy ZPL
solve(ZPL)
get.variables(ZPL)
#get.dual.solution(ZPL)
w5<-get.variables(ZPL)[2*T+2]
w6<-get.variables(ZPL)[2*T+3]
return(c(w5,w6))
}

###################   portfel z³o¿ony z akcji 1 i 2 (wagi: w7, w8) 
###################   maksymalizuj¹cy najwiêksz¹ stratê  WL (Worst Loss)
portfel.WL<-function(r1,r2){
  T<-length(r1)
  r1.sr<-mean(r1)
  r2.sr<-mean(r2)
  sigma1<-sd(r1)
  sigma2<-sd(r2)
  ro12<-cor(r1,r2)
  r0=-1
ZPL <- make.lp(nrow=0, ncol=4)
lp.control(ZPL, sense="max",simplextype="primal")
#funkcja celu
set.objfn(ZPL, c(0,0,1,0)) #w7, w8, y, \mu
for (i in 1:T) {add.constraint(ZPL, c(exp(r1[i])-1, exp(r2[i])-1,-1,0), ">=", 0)}
add.constraint(ZPL, c(1,1,0,0), "=", 1)
add.constraint(ZPL, c(0,0,0,1), ">=", r0)
add.constraint(ZPL, c(exp(r1.sr)-1, exp(r2.sr)-1,0,-1), "=", 0)
set.bounds(ZPL, columns=3, lower=-1)
set.bounds(ZPL, columns=c(1,2), lower=c(0,0), upper=c(1,1))
set.bounds(ZPL, columns=4, lower=-1)
# z krótk¹ sprzeda¿¹:
#set.bounds(ZPL, columns=c(1,2), lower=c(-Inf, -Inf))
ZPL
#rozwi¹zujemy ZPL
solve(ZPL)
get.variables(ZPL)
#get.dual.solution(ZPL)
select.solution(ZPL)
w7<-get.variables(ZPL)[1]
w8<-get.variables(ZPL)[2]
return(c(w7,w8))
}


#czy trzeba tutaj po prostu w r1[1:36] zamiast [1:36] daæ [1:12] i odpowiednio skalibrowaæ? Przysi¹dê.
r.portfela<-rep(0, times=36); #matrix(c(1,0,0,0), nrow=1)
for (i in 1:36) {
  if (i==1) {
            #portfel.wagi<-portfel.markovitz(as.vector(r1[1:12]), as.vector(r2[1:12])); 
            portfel.wagi<-portfel.MAD(as.vector(r1[1:12]), as.vector(r2[1:12])); 
            #portfel.wagi<-portfel.WL(as.vector(r1[1:12]), as.vector(r2[1:12])); 
            r.portfela[1]<-portfel.wagi[1]*(exp(as.numeric(r1[12]))-1)+portfel.wagi[2]*(exp(as.numeric(r2[12]))-1)
  }
  else {
            #portfel.wagi<-portfel.markovitz(as.vector(r1[(i:(12+i-1))]), as.vector(r2[(i:(12+i-1))])); 
            portfel.wagi<-portfel.MAD(as.vector(r1[(i:(12+i-1))]), as.vector(r2[(i:(12+i-1))])); 
            #portfel.wagi<-portfel.WL(as.vector(r1[(i:(12+i-1))]), as.vector(r2[(i:(12+i-1))])); 
            r.portfela[i]<-portfel.wagi[1]*(exp(as.numeric(r1[(12+i)]))-1)+portfel.wagi[2]*(exp(as.numeric(r2[(12+i)]))-1)}
}
r.portfela
#r.markovitz<-r.portfela
r.MAD<-r.portfela
#r.WL<-r.portfela

#r.markovitz<-zoo(r.markovitz, order.by=daty.miesiace[13:48])
r.MAD<-zoo(r.MAD, order.by=daty.miesiace[13:48])
#r.WL<-zoo(r.WL, order.by=daty.miesiace[13:48])
r.markovitz
r.MAD
r.WL

r.akcja1<-exp(r1[12:47])-1
r.akcja1<-zoo(r.akcja1, order.by=daty.miesiace[13:48])
r.akcja2<-exp(r2[12:47])-1
r.akcja2<-zoo(r.akcja2, order.by=daty.miesiace[13:48])

plot(r.akcja1, xlab=c("Miesi¹ce"), ylab=c("Prosta stopa zwrotu"))
lines(r.akcja2, col="grey")
lines(r.markovitz, col="red")
lines(r.MAD, col="orange")
lines(r.WL, col="blue")
legend("topleft",c("Akcja 1", "Akcja 2", "Markovitz", "MAD", "WL"),
      lty=c(1,1,1,1,1), col=c("black","grey","red","orange","blue"), cex=0.6)

#WskaŸniki
b.r<-exp(b.r[12:47])-1

a.r.markovitz<-as.vector(r.markovitz)-as.vector(b.r)
a.r.MAD<-as.vector(r.MAD)-as.vector(b.r)
a.r.WL<-as.vector(r.WL)-as.vector(b.r)

#TE
N.r<-length(b.r)

tracking.error.a.r.markovitz<-sd(a.r.markovitz)*sqrt((N.r-1)/N.r)
tracking.error.a.r.MAD<-sd(a.r.MAD)*sqrt((N.r-1)/N.r)
tracking.error.a.r.WL<-sd(a.r.WL)*sqrt((N.r-1)/N.r)

roczny.te.a.r.markovitz<-sqrt(12)*tracking.error.a.r.markovitz
roczny.te.a.r.MAD<-sqrt(12)*tracking.error.a.r.MAD
roczny.te.a.r.WL<-sqrt(12)*tracking.error.a.r.WL

#IR
roczny.ie.a.r.markovitz<-mean(a.r.markovitz)*(12)/roczny.te.a.r.markovitz
roczny.ie.a.r.MAD<-mean(a.r.MAD)*(12)/roczny.te.a.r.MAD
roczny.ie.a.r.WL<-mean(a.r.WL)*(12)/roczny.te.a.r.WL

#Sharpe
sigma.markovitz<-sd(r.markovitz)*sqrt((N.r-1)/N.r)*sqrt(12)
sigma.MAD<-sd(r.MAD)*sqrt((N.r-1)/N.r)*sqrt(12)
sigma.WL<-sd(r.WL)*sqrt((N.r-1)/N.r)*sqrt(12)
r.f<-mean(rf)*12

SR.markovitz<-(mean(r.markovitz)*(12)-r.f)/sigma.markovitz
SR.MAD<-(mean(r.MAD)*(12)-r.f)/sigma.MAD
SR.WL<-(mean(r.WL)*(12)-r.f)/sigma.WL

#M^2
sigma.b<-sd(b.r)*sqrt((N.r-1)/N.r)*sqrt(12)

m2.markovitz<-mean(r.markovitz)*(12)+SR.markovitz*(sigma.b-sigma.markovitz)
m2.MAD<-mean(r.MAD)*(12)+SR.MAD*(sigma.b-sigma.MAD)
m2.WL<-mean(r.WL)*(12)+SR.WL*(sigma.b-sigma.WL)

#Tabela
tabela<-data.frame(
  portfel=c('Markovitz','MAD','WL'),
  TE=c(roczny.te.a.r.markovitz, roczny.te.a.r.MAD, roczny.te.a.r.WL),
  IR=c(roczny.ie.a.r.markovitz, roczny.ie.a.r.MAD, roczny.ie.a.r.WL),
  SR=c(SR.markovitz, SR.MAD, SR.WL),            
  m2=c(m2.markovitz, m2.MAD, m2.WL)
)
tabela

tabela2<-data.frame(
  portfel=c('Markovitz','MAD','WL'),
  stopazwrotu=c(mean(r.markovitz),mean(r.MAD),mean(r.WL))
)
tabela2

