#Usuniêcie starych zmiennych
rm(list=ls())

#parametry potrzebne przy strategii œrednich ruchomych
p=1; q=5

#Wczytanie danych 
 aktywo <- "RANKPROGR.mst";n.aktywo="Rank Progress"
# aktywo <- "NETIA.mst";n.aktywo="Netia"
#aktywo <- "ORANGEPL.mst"; n.aktywo="Orange"

benchmark <- "WIG-NRCHOM.mst" #WIG WIG20
stopa.wolna <- "WIBID1M.csv" 

dane.aktywo <- read.csv2(file = aktywo, head=TRUE, sep=",", dec=".")
data.format <- "%Y%m%d"
require(zoo)
cena.otwarcia <- dane.aktywo$"X.OPEN."
daty <- as.Date(as.character(dane.aktywo$"X.DTYYYYMMDD."), format=data.format)
c <- zoo(cena.otwarcia , order.by = daty)

dane.benchmark <- read.csv2(file = benchmark, head=TRUE, sep=",", dec=".")
b.cena.otwarcia <- dane.benchmark$"X.OPEN."
daty <- as.Date(as.character(dane.benchmark$"X.DTYYYYMMDD."), format=data.format)
b.c <- zoo(b.cena.otwarcia , order.by = daty)

dane.stopa.wolna <- read.csv2(file = stopa.wolna, head=TRUE, sep=",", dec=".")
data.format <- "%Y-%m-%d"
rf <- dane.stopa.wolna$"Otwarcie"
daty <- as.Date(as.character(dane.stopa.wolna$"Data"), format=data.format)
rf <- zoo(rf/1200 , order.by = daty) #dzielnie bo w procentach i skali rocznej

#Ustalenie daty pocz¹tku i koñca próby
data.start <- "2017-01-01"
data.koniec <-"2020-03-31"

require(xts)
c1<-window(c, start=as.Date(data.start), end=as.Date(data.koniec)) 
b.c1<-window(b.c, start=as.Date(data.start), end=as.Date(data.koniec)) 
rf<-window(rf, start=as.Date(data.start), end=as.Date(data.koniec))

#konwersja szeregu danych dziennych na dane miesiêczne
#aktywo
c1.m <- apply.monthly(c1, last) #first, mean, mojafunkcja
K<-60+q # d³ugoœæ okresu próby, 5 lat + dane poprzedzaj¹ce ten okres potrzembne dla decyzji strategii w pierwszym miesi¹cu
c1.m <-c1.m[-c(1:(length(c1.m)-K))]
daty.miesiace<-time(c1.m)
N.c1.m<-length(c1.m) #liczba miesiêcy podlegaj¹cych analizie (dla cen)
#benchmark
b.c1.m <- apply.monthly(b.c1, last) #first, mean, mojafunkcja
K<-60 
b.c1.m <-b.c1.m[-c(1:(length(b.c1.m)-K-1))]
b.c1.m <-zoo(b.c1.m, order.by = daty.miesiace[-(1:q)])
N.b.c1.m<-length(b.c1.m)
#stopa wolna od ryzyka
rf <- apply.monthly(rf, last) #first, mean, mojafunkcja
K<-60 
rf <-rf[-c(1:(length(rf)-K-1))]
rf<-rf[-c(length(rf))]
rf <- zoo(rf, order.by = daty.miesiace[-(1:q)]) #time(r))

#logarytmiczna stopa zwrotu i jej charakterystyki
r <- log(as.vector(c1.m[-1]))-log(as.vector(c1.m[-(N.c1.m)])) 
r <- zoo(r , order.by = daty.miesiace[-1])
r<-tail(r, n=K)#r[-(1:q)]
N.r<-length(r)

b.r <- log(as.vector(b.c1.m[-1]))-log(as.vector(b.c1.m[-(N.b.c1.m)])) 
b.r <- zoo(b.r , order.by = time(r))
b.r<-tail(b.r, n=K)

dane.widok<-merge(c1.m,b.c1.m,rf,momentum)
head(dane.widok,n=12)
tail(dane.widok)

dane.widok<-merge(r,b.r,rf)
head(dane.widok,n=12)
tail(dane.widok)

#################################################################################


#strategia naiwna
s.naiwna<-rep(0, N.r)
s.naiwna <- zoo(s.naiwna , order.by = time(r))
delta=0.01
if(as.vector(c1.m[q])/as.vector(c1.m[q-1])>(1+delta))
  {s.naiwna[1]<-1} else {s.naiwna[1]<-0} 
for (i in index(as.vector(r))[-1]) {
  if(as.vector(c1.m[i-1+q])/as.vector(c1.m[i-2+q])>(1+delta)) 
    {s.naiwna[i]<-1} else 
    if(as.vector(c1.m[i-1+q])/as.vector(c1.m[i-2+q])<(1-delta)) 
      {s.naiwna[i]<-0} else 
      {s.naiwna[i]<-s.naiwna[i-1]}
}
r.s.naiwna<-rep(0, N.r)
r.s.naiwna[which(as.vector(s.naiwna)==1)]<-r[which(as.vector(s.naiwna)==1)]
r.s.naiwna[which(as.vector(s.naiwna)==0)]<-rf[which(as.vector(s.naiwna)==0)]
r.s.naiwna<-zoo(r.s.naiwna , order.by = time(r))


plot(r, type="l")
lines(r.s.naiwna, col="red")
abline(h=0)
mean(r)*12
mean(r.s.naiwna)*12
sd(r)*sqrt((N.r-1)/N.r)*sqrt(12)
sd(r.s.naiwna)*sqrt((N.r-1)/N.r)*sqrt(12)

# dane.widok<-merge(c1.m,r,b.c1.m,rf, s.naiwna, r.s.maiwna)
# head(dane.widok,n=12)
# tail(dane.widok)

#strategia moementum dwuokresowe (j=2)
j=2 #pamiêtajmy, ¿e j<=q
s.momentum<-rep(0, N.r)
s.momentum <- zoo(s.momentum , order.by = time(r))
momentum<-tail(as.vector(c1.m),n=(N.r+q-j)) - head(as.vector(c1.m),n=(N.r+q-j))
momentum<-zoo(momentum , order.by = time(c1.m)[-(1:(q-j))])
if(as.vector(momentum[1])<0 & as.vector(momentum[2])>0) {s.momentum[1]<-1} else {s.momentum[1]<-0} 
for (i in index(as.vector(r))[-1]) {
  if(as.vector(momentum[i])<0 & as.vector(momentum[i+1])>0) {s.momentum[i]<-1} else 
    if(as.vector(momentum[i])>0 & as.vector(momentum[i+1])<0) {s.momentum[i]<-0} else {s.momentum[i]<-s.momentum[i-1]}
}
r.s.momentum<-rep(0, N.r)
r.s.momentum[which(as.vector(s.momentum)==1)]<-r[which(as.vector(s.momentum)==1)]
r.s.momentum[which(as.vector(s.momentum)==0)]<-rf[which(as.vector(s.momentum)==0)]
r.s.momentum<-zoo(r.s.momentum , order.by = time(r))
# dane.widok<-merge(c1.m,r,b.c1.m,rf,momentum, s.momentum, r.s.momentum)
# head(dane.widok,n=12)
# tail(dane.widok)

#strategia œrednie ruchome p,q
r.s.ruchome<-rep(0, N.r)->p.ruchome->q.ruchome
r.s.ruchome <- zoo(r.s.ruchome , order.by = time(r))
#p.ruchome<-c()
#q.ruchome<-c()
for (i in index(as.vector(r))) {
  p.ruchome[i]<-mean(c1.m[(i+(q-p)):(i+q-1)])
  q.ruchome[i]<-mean(c1.m[(i):(i+q-1)])
}
p.ruchome <- zoo(p.ruchome , order.by = time(r))
q.ruchome <- zoo(q.ruchome , order.by = time(r))

if(as.vector(p.ruchome[1])>as.vector(q.ruchome[1])) {r.s.ruchome[1]<-r[1]} else 
  if(as.vector(p.ruchome[1])<as.vector(q.ruchome[1])) {r.s.ruchome[1]<-rf[1]}

for (i in index(as.vector(r))[-1]) {
  if(as.vector(p.ruchome[i])>as.vector(q.ruchome[i])) {r.s.ruchome[i]<-r[i]} else 
    if(as.vector(p.ruchome[i])<as.vector(q.ruchome[i])) {r.s.ruchome[i]<-rf[i]} else {r.s.ruchome[i]<-r.s.ruchome[i-1]}
}
# dane.widok<-merge(c1.m,r,b.c1.m,rf,p.ruchome, q.ruchome, r.s.ruchome)
# head(dane.widok,n=12)
# tail(dane.widok)

#Wykresy stóp zwrotu poszczególnych strategii
par(mfrow=c(3,1))
plot(r, type="l", main="strategia naiwna")
lines(r.s.naiwna, col="red")
abline(h=0)
plot(r, type="l", main="strategia momentum")
lines(r.s.momentum, col="blue")
abline(h=0)
plot(r, type="l", main="strategia œrednich ruchomych")
lines(r.s.ruchome, col="green")
abline(h=0)
par(mfrow=c(1,1))

#œrednia roczna stopa zwrotu i roczne odch. standardowe dla trzech strategii
r.aktywa<-mean(r)*12
r.naiwna<-mean(r.s.naiwna)*12
r.momentum<-mean(r.s.momentum)*12
r.ruchome<-mean(r.s.ruchome)*12
r.f<-mean(rf)*12
r.b<-mean(b.r)*12
sigma.aktywa<-sd(r)*sqrt((N.r-1)/N.r)*sqrt(12)
sigma.naiwna<-sd(r.s.naiwna)*sqrt((N.r-1)/N.r)*sqrt(12)
sigma.momentum<-sd(r.s.momentum)*sqrt((N.r-1)/N.r)*sqrt(12)
sigma.ruchome<-sd(r.s.ruchome)*sqrt((N.r-1)/N.r)*sqrt(12)
sigma.b<-sd(b.r)*sqrt((N.r-1)/N.r)*sqrt(12)

#mapa ryzyko-dochod dla strategii par(mfrow=c(1,1))
ryzyko<-c(sigma.aktywa,sigma.naiwna,sigma.momentum,sigma.ruchome)
dochod<-c(r.aktywa,r.naiwna,r.momentum,r.ruchome)

x.pocz<-min(ryzyko)-0.03
x.koniec<-max(ryzyko)+0.03
y.pocz<-min(dochod)-0.01
y.koniec<-max(dochod)+0.01

dane<-data.frame(ryzyko, dochod)
plot(dane, type = "p", col = "green", pch = 16, 
     xlim=c(x.pocz,x.koniec), ylim = c(y.pocz,y.koniec),
     xlab=c("annualizowane odchylenie standardowe"), ylab=c("annualizowana stopa zwrotu"),
     main=paste("Mapa ryzyko-dochód:", n.aktywo,", \n", "strategie:naiwna,", "momentum i œrednich ruchomych"))
     
text(dane, 
     labels = c(n.aktywo,"s. naiwna,", "s. momentum", paste("s.œrednich ruchomych \n p=", p, "q=", q)), 
                      pos = 4, cex = 0.5)
abline(h=0, col="red")

#zapisujemy ceny aktywa i stopy zwrotu strategii oraz benchmarku w tabeli
# dane.tabela<-merge(c1.m,r,b.c1.m,b.r,rf,r.s.naiwna, momentum, r.s.momentum,
#                    p.ruchome,q.ruchome, r.s.ruchome)
# dane.tabela[1:10]

#################################################################################################
#   Ocena strategii - relacja stopy zwrotu strategi benchmaku na rynku wzrostowym i spadkowym   #

b.obs.plus<-which(b.r>0)
b.obs.minus<-which(b.r<0)

UpCaptureIndicator.r<-mean(as.vector(r)[b.obs.plus])/mean(as.vector(b.r)[b.obs.plus])
DownCaptureIndicator.r<-mean(as.vector(r)[b.obs.minus])/mean(as.vector(b.r)[b.obs.minus])
UpNumberRatio.r<-length(which(as.vector(r)[b.obs.plus]>0))/length(b.obs.plus)
DownNumberRatio.r<-length(which(as.vector(r)[b.obs.minus]<0))/length(b.obs.minus)
PercentageGainRatio.r<-length(as.vector(r)[r>0])/length(b.obs.plus)

UpCaptureIndicator.naiwna<-mean(as.vector(r.s.naiwna)[b.obs.plus])/mean(as.vector(b.r)[b.obs.plus])
DownCaptureIndicator.naiwna<-mean(as.vector(r.s.naiwna)[b.obs.minus])/mean(as.vector(b.r)[b.obs.minus])
UpNumberRatio.naiwna<-length(which(as.vector(r.s.naiwna)[b.obs.plus]>0))/length(b.obs.plus)
DownNumberRatio.naiwna<-length(which(as.vector(r.s.naiwna)[b.obs.minus]<0))/length(b.obs.minus)
PercentageGainRatio.naiwna<-length(as.vector(r.s.naiwna)[r.s.naiwna>0])/length(b.obs.plus)

UpCaptureIndicator.momentum<-mean(as.vector(r.s.momentum)[b.obs.plus])/mean(as.vector(b.r)[b.obs.plus])
DownCaptureIndicator.momentum<-mean(as.vector(r.s.momentum)[b.obs.minus])/mean(as.vector(b.r)[b.obs.minus])
UpNumberRatio.momentum<-length(which(as.vector(r.s.momentum)[b.obs.plus]>0))/length(b.obs.plus)
DownNumberRatio.momentum<-length(which(as.vector(r.s.momentum)[b.obs.minus]<0))/length(b.obs.minus)
PercentageGainRatio.momentum<-length(as.vector(r.s.momentum)[r.s.momentum>0])/length(b.obs.plus)

UpCaptureIndicator.ruchome<-mean(as.vector(r.s.ruchome)[b.obs.plus])/mean(as.vector(b.r)[b.obs.plus])
DownCaptureIndicator.ruchome<-mean(as.vector(r.s.ruchome)[b.obs.minus])/mean(as.vector(b.r)[b.obs.minus])
UpNumberRatio.ruchome<-length(which(as.vector(r.s.ruchome)[b.obs.plus]>0))/length(b.obs.plus)
DownNumberRatio.ruchome<-length(which(as.vector(r.s.ruchome)[b.obs.minus]<0))/length(b.obs.minus)
PercentageGainRatio.ruchome<-length(as.vector(r.s.ruchome)[r.s.ruchome>0])/length(b.obs.plus)

tabela.podsumowanie<-data.frame(
  nazwa=rep(dane.aktywo$X.TICKER.[1], times=4),
  strategia=c("aktywo", "naiwna", "momentum", "œrednie ruchome"),
  s.zwrotu=c(r.aktywa*100, r.naiwna*100, r.momentum*100, r.ruchome*100),
  sigma.r=c(sigma.aktywa, sigma.naiwna, sigma.momentum, sigma.ruchome),
  UpCaptureIndicator=c(UpCaptureIndicator.r,UpCaptureIndicator.naiwna, UpCaptureIndicator.momentum, UpCaptureIndicator.ruchome),
  DownCaptureIndicator=c(DownCaptureIndicator.r,DownCaptureIndicator.naiwna, DownCaptureIndicator.momentum, DownCaptureIndicator.ruchome),
  UpNumberRatio=c(UpNumberRatio.r,UpNumberRatio.naiwna, UpNumberRatio.momentum, UpNumberRatio.ruchome),
  DownNumberRatio=c(DownNumberRatio.r,DownNumberRatio.naiwna, DownNumberRatio.momentum, DownNumberRatio.ruchome),
  PercentageGainRatio=c(PercentageGainRatio.r,PercentageGainRatio.naiwna, PercentageGainRatio.momentum, PercentageGainRatio.ruchome)
)




#################################################################################################
#   Ocena strategii - relacja stopy zwrotu strategi benchmaku na rynku wzrostowym i spadkowym   #
#                                z wykorzystaniem Excess Return                                 #

#definiujemy arytmetyczne nadwy¿kowe stopy zwrotu (Excess Return)
a.r<-as.vector(r)-as.vector(b.r)
a.r.s.naiwna<-as.vector(r.s.naiwna)-as.vector(b.r)
a.r.s.momentum<-as.vector(r.s.momentum)-as.vector(b.r)
a.r.s.ruchome<-as.vector(r.s.ruchome)-as.vector(b.r)


UpPercentageRatio.r<-length(which(a.r[b.obs.plus]>0))/length(b.obs.plus)
DownPercentageRatio.r<-length(which(a.r[b.obs.minus]>0))/length(b.obs.minus)

UpPercentageRatio.naiwna<-length(which(a.r.s.naiwna[b.obs.plus]>0))/length(b.obs.plus)
DownPercentageRatio.naiwna<-length(which(a.r.s.naiwna[b.obs.minus]>0))/length(b.obs.minus)

UpPercentageRatio.momentum<-length(which(a.r.s.momentum[b.obs.plus]>0))/length(b.obs.plus)
DownPercentageRatio.momentum<-length(which(a.r.s.momentum[b.obs.minus]>0))/length(b.obs.minus)

UpPercentageRatio.ruchome<-length(which(a.r.s.ruchome[b.obs.plus]>0))/length(b.obs.plus)
DownPercentageRatio.ruchome<-length(which(a.r.s.ruchome[b.obs.minus]>0))/length(b.obs.minus)


tabela.podsumowanie1<-data.frame(
  UpPercentageRatio=c(UpPercentageRatio.r,UpPercentageRatio.naiwna, UpPercentageRatio.momentum, UpPercentageRatio.ruchome),
  DownPercentageRatio=c(DownPercentageRatio.r,DownPercentageRatio.naiwna, DownPercentageRatio.momentum, DownPercentageRatio.ruchome
  )
)

#Doklejamy do poprzedniej tabeli podsumowuj¹cej nowe wyniki
tabela.podsumowanie<-cbind(tabela.podsumowanie,tabela.podsumowanie1)

#################################################################################################
#                           Ocena strategii - tracking error i information ratio                #


tracking.error.a.r<-sd(a.r)*sqrt((N.r-1)/N.r)
tracking.error.a.r.s.naiwna<-sd(a.r.s.naiwna)*sqrt((N.r-1)/N.r)
tracking.error.a.r.s.momentum<-sd(a.r.s.momentum)*sqrt((N.r-1)/N.r)
tracking.error.a.r.s.ruchome<-sd(a.r.s.ruchome)*sqrt((N.r-1)/N.r)

roczny.te.a.r<-sqrt(12)*tracking.error.a.r
roczny.te.a.r.s.naiwna<-sqrt(12)*tracking.error.a.r.s.naiwna
roczny.te.a.r.s.momentum<-sqrt(12)*tracking.error.a.r.s.momentum
roczny.te.a.r.s.ruchome<-sqrt(12)*tracking.error.a.r.s.ruchome

# i information ratio
roczny.ie.a.r<-mean(a.r)*(12)/roczny.te.a.r
roczny.ie.a.r.s.naiwna<-mean(a.r.s.naiwna)*(12)/roczny.te.a.r.s.naiwna
roczny.ie.a.r.s.momentum<-mean(a.r.s.momentum)*(12)/roczny.te.a.r.s.momentum
roczny.ie.a.r.s.ruchome<-mean(a.r.s.ruchome)*(12)/roczny.te.a.r.s.ruchome

tabela.podsumowanie2<-data.frame(
  te=c(roczny.te.a.r,roczny.te.a.r.s.naiwna, roczny.te.a.r.s.momentum, roczny.te.a.r.s.ruchome),
  ie=c(roczny.ie.a.r,roczny.ie.a.r.s.naiwna, roczny.ie.a.r.s.momentum, roczny.ie.a.r.s.ruchome)
)

# I znowu doklejamy do poprzedniej tabeli podsumowuj¹cej nowe wyniki
tabela.podsumowanie<-cbind(tabela.podsumowanie,tabela.podsumowanie2)

# 1. jeœli chcemy sklejaæ wyniki dla ró¿nych aktywów, to kolejn¹ linijkê kodu nale¿y wykonaæ dla pierwszego aktywa 
tabela.wynikow<-t(tabela.podsumowanie)
# 2. i nastêpne: sklejamy wyniki dla poprzedniego i obecnego aktywa
tabela.wynikow<-cbind(tabela.wynikow,t(tabela.podsumowanie))





