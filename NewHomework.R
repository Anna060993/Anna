rm(list=ls())
ls()
options(warn=-1)

#### GARCH(1,1) con media non nulla e innovazioni gaussiane ####

#######################
# Funzioni utilizzate nel corso del programma
#######################

# Funzione per simulare una serie storica di lunghezza arbitraria e con un generico vettore di parametri.
# In input: n (ampiezza della serie) e theta (vettore contente 4 parametri, nell'ordine: mi,omega,alpha,beta).
# In output: una matrice di dimensione nx2, la prima colonna contiene i valori della serie,
# la seconda quelli della varianza condizionata.

garch11<-function(n,n.start=200,theta=c(mi,omega,alpha,beta)) 
{
  mi<-theta[1]
  omega<-theta[2]
  alpha<-theta[3]
  beta<-theta[4]
  y<-rep(0,(n+n.start))        # Serie dei livelli
  h<-rep(0,(n+n.start))        # Serie delle arianze condizionate
  h[1]<-omega/(1-alpha-beta)   # Si inizializza con la varianza marginale
  y[1]<-rnorm(1,mi,sqrt(h[1])) # Coerentemente con l'informazione disponibile, si genera un primo valore per 
                               # inizializzare la serie dei livelli da una N(mi,h[1]) 
  for (i in 2:(n+n.start)) 
  {
    h[i]<-omega+alpha*(y[i-1]-mi)^2+beta*h[i-1]
    y[i]<-mi+rnorm(1)*sqrt(h[i]) 
  }
  
  serie<-y[(n.start+1):(n+n.start)]
  var.cond<-h[(n.start+1):(n+n.start)]
  return(cbind(serie,var.cond)) 
}

#### Funzione per il calcolo della log-verosimiglianza ####

lvgarch11<-function(lambda) 
{
  mi<-lambda[1]
  omega<-exp(lambda[2])      
  alpha<-exp(lambda[3])
  beta<-exp(lambda[4])
  e<-y-mi # Serie dei residui
  e2<-e^2
  h<-rep(0,n)
  h[1]<-omega/(1-alpha-beta)
  
  for (i in 2:n) 
  {
    h[i]<-omega+alpha*e2[i-1]+beta*h[i-1] 
  }
  
  h<-h[2:n]
  e2<-e2[2:n]
  # Stima della sola log-verosimiglianza condizionata:
  l<-sum(-0.5*log(h)-0.5*(e2/h))
  # Per la stazionarietà in senso debole: si assegnano valori estremamente bassi della log-verosimiglianza
  # a quelle coppie (alpha,beta) per cui tale condizione non è verificata.
  if (alpha+beta>=1) l<--10^(9)
  return(l) 
}

#### Funzione per il calcolo del gradiente: ####

grad11<-function(lambda)
{
  mi<-lambda[1]
  # La trasformazione esponenziale assicura il vincolo di positività per i parametri omega, alpha e beta:
  omega<-exp(lambda[2])
  alpha<-exp(lambda[3])
  beta<-exp(lambda[4])
  e<-y-mi
  e2<-e^2
  h<-rep(0,n)
  h[1]<-omega/(1-alpha-beta)
  for (i in 2:n) 
  {
    h[i]<-omega+alpha*e2[i-1]+beta*h[i-1] 
  }
 
  c1<-rep(0,n-1)
  c2<-cumsum(beta^seq(0,(n-2),1))
  c3<-rep(0,n-1)
  c4<-rep(0,n-1)
  
  for (i in 2:n) 
  {
    c1[i-1]<--2*alpha*sum((beta^seq(0,i-2,1))*e[(i-1):1])
    c3[i-1]<-sum((beta^seq(0,i-2,1))*e2[(i-1):1]) 
  }
  for (i in 3:n) 
  {
    c4[i-1]<-omega*sum(seq(1,i-2,1)*(beta^seq(0,i-3,1)))+alpha*sum(seq(1,i-2,1)*(beta^seq(0,i-3,1))*e2[(i-2):1]) 
  }  
  
  dlvmi<--0.5*sum((1/h[2:n])*c1)+sum(e[2:n]/h[2:n])+0.5*sum((e2[2:n]/(h[2:n]^2))*c1) 
  dlvomega<--0.5*sum((1/h[2:n])*c2)+0.5*sum((e2[2:n]/(h[2:n]^2))*c2) 
  dlvalpha<--0.5*sum(c3/h[2:n])+0.5*sum(c3*e2[2:n]/(h[2:n]^2))
  dlvbeta<--0.5*sum(c4/h[2:n])+0.5*sum(c4*e2[2:n]/(h[2:n]^2))
  
  dlvomega<-dlvomega*exp(lambda[2])
  dlvalpha<-dlvalpha*exp(lambda[3])
  dlvbeta<-dlvbeta*exp(lambda[4])
  return(c(dlvmi,dlvomega,dlvalpha,dlvbeta)) 
}

#### Funzione per il calcolo della matrice hessiana: ####

hessian<-function(theta) 
{
  h<-rep(0,n)
  mi<-theta[1]
  omega<-theta[2]
  alpha<-theta[3]
  beta<-theta[4]
  h[1]<-omega/(1-alpha-beta)
  e<-y-mi 
  e2<-e^2
  e1<-y1-mi
  e12<-e1^2
  
  for (i in 2:n) 
  { 
    h[i]<-omega+alpha*e2[i-1]+beta*h[i-1] 
  }  
  h1<-c(h[1],h[1:n-1])	        
  h1<-h1[2:n]
    
  c1<-rep(0,n-1) 
  c2<-cumsum(beta^seq(0,(n-2),1))
  c3<-rep(0,n-1)
  c4<-rep(0,n-1)
  c5<-2*alpha*cumsum(beta^seq(0,(n-2),1))
  c6<-rep(0,n-1)
  c7<-rep(0,n-1)
  c8<-rep(0,n-1)
  c9<-rep(0,n-1)
  c10<-rep(0,n-1)
    
  for (i in 2:n) 
  {
    c1[i-1]<--2*alpha*sum(beta^(seq(0,i-2,1))*e[(i-1):1])
    c3[i-1]<-sum((beta^seq(0,i-2,1))*e2[(i-1):1]) 
    c6[i-1]<--2*sum((beta^seq(0,i-2,1))*e[(i-1):1]) 
  }
  
  for (i in 3:n) 
  {
    c4[i-1]<-omega*sum(seq(1,i-2,1)*(beta^seq(0,i-3,1)))+alpha*sum(seq(1,i-2,1)*(beta^seq(0,i-3,1))*e2[(i-2):1])
    c7[i-1]<--2*alpha*sum(seq(1,i-2,1)*beta^seq(0,i-3,1)*e[(i-2):1])
    c8[i-1]<-sum(seq(1,i-2,1)*(beta^seq(0,i-3,1)))
    c9[i-1]<-sum(seq(1,i-2,1)*(beta^seq(0,i-3,1))*e2[(i-2):1]) 
  }
  
  for (i in 4:n) 
  {
    c10[i-1]<-omega*(sum(seq(1,i-3,1)*seq(2,i-2,1)*(beta^seq(0,i-4,1))))+
     alpha*(sum(seq(1,i-3,1)*seq(2,i-2,1)*(beta^seq(0,i-4,1)))*e2[(i-3):1]) 
  }
   
  h<-h[2:n]
  e<-e[2:n]
  e2<-e2[2:n]
  hess<-matrix(0,nrow=4,ncol=4)
  
  hess[1,1]<-0.5*sum((1/(h^2))*(c1^2))-0.5*sum(c5/h)-sum(1/h)-2*sum((c1*e)/(h^2))-sum((c1^2)*e2/(h^3))+
    sum(c5*e2/(h^2))
  hess[1,2]<-0.5*sum((c1*c2)/(h^2))-sum(c2*(e/(h^2)))-sum(c1*c2*(e2/(h^3)))
  hess[1,3]<-0.5*sum((1/(h^2))*c1*c3)-0.5*sum(c6/h)+sum(c6*e/h)-sum(c3*e/(h^2))-sum(c3*c1*(e2/(h^3)))+
    0.5*sum(c6*e2/(h^2))  
  hess[1,4]<-0.5*sum((1/(h^2))*c1*c4)-0.5*sum(c7/h)+sum((c7*e)/h)-sum((c4*e)/(h^2))-sum(c4*c1*(e2/(h^3)))+
    0.5*sum((c7*e2)/h^2)  
  hess[2,2]<-0.5*sum((c2^2)/(h^2))-sum((e2*(c2^2))/(h^3))  
  hess[2,3]<-0.5*sum((c2*c3)/(h^2))-sum((e2*(c2*c3))/(h^3))   
  hess[2,4]<-0.5*sum((c2*c4)/(h^2))-0.5*sum(c8/h)-sum((c2*c4*e2)/(h^3))+0.5*sum((c8*e2)/(h^2))
  hess[3,3]<-0.5*sum((c3^2)/(h^2))-sum((e2*(c3^2))/(h^3))
  hess[3,4]<-0.5*sum((c3*c4)/(h^2))-0.5*sum(c9/h)-sum((c3*c4*e2)/(h^3))+0.5*sum((c9*e2)/(h^2))
  hess[4,4]<-0.5*sum((c4^2)/(h^2))-0.5*sum(c10/h)-sum(((c4^2)*e2)/(h^3))+0.5*sum((c10*e2)/(h^2))
  hess[2,1]<-hess[1,2]   
  hess[3,1]<-hess[1,3]
  hess[4,1]<-hess[1,4] 
  hess[3,2]<-hess[2,3]
  hess[4,2]<-hess[2,4]
  hess[4,3]<-hess[3,4]
  return(hess) 
}

########################################################
#### Esempio numerico ####
########################################################
# Costruzione di un modello GARCH
# 0.Simulazione
# 1.Analisi preliminari 
# 2.Stima
# 3.Diagnostica


# Supponiamo di generare una serie di 500 osservazioni da un GARCH(1,1) con innovazioni gaussiane ed i 
# seguenti valori dei parametri:
mi<-0.10
omega<-0.05
alpha<-0.2
beta<-0.7

# Poichè la funzione "garch11" restituisce una matrice, per recuperare la serie dei livelli si seleziona 
# la prima colonna. La seconda colonna restituisce la varianza condizionata.
n<-500
set.seed(100)
dati<-garch11(n)
y<-dati[,1]
var.cond<-dati[,2] 

# E' utile avere un riscontro grafico del fatto che la serie simulata presenta le tipiche caratteristiche 
# empiriche della serie dei rendimenti.
par(mfrow=c(2,1),mar=c(4,4,2,1))
plot(y,type="l",main="Serie dei livelli",xlab="t")
plot(var.cond,type="l",main="Serie delle varianze condizionate",xlab="t",ylab="var. cond.")

# Operazioni preliminari:
y2<-y^2                  	# valori di y(t)^2
y1<-c(mean(y),y[1:n-1])   # valori di y(t)
y12<-y1^2 				        # valori di y(t-1)

# E' utile avere un riscontro grafico del fatto che la serie simulata presenti le tipiche caratteristiche 
# empiriche delle serie dei rendimenti. Si tratta di semplici analisi preliminari, indispensabili quando
# si vuole risalire al vero processo generatore dei dati non noto. In questo caso, in cui è noto dato che 
# la serie è stata simulata, tali analisi servono solo a confermare alcune note evidenze empiriche.
# I rendimenti risultano stazionari in media (pari a mi) ed eteroschedastici.
par(mfrow=c(2,1),mar=c(4,4,2,1))
plot(y,type="l",main="Serie dei livelli",xlab="t")
abline(h=.1,col=2)
plot(var.cond,type="l",main="Serie delle varianze condizionate",xlab="t",ylab="var. cond.")
# La serie dei livelli è incorrelata:
par(mfrow=c(2,1),mar=c(4,4,3,3))
acf(y)
pacf(y)
# Mentre non lo è quella dei quadrati dei livelli:
acf(y2)
pacf(y2)
# La distribuzione marginale del processo leptocurtica. 
asimm=function(x)
{
  n=length(x)
  m3=mean((x-mean(x))^3)
  m2=sd(x)
  return(m3/m2^3)
}
curtosi<-function(x) # Funzione per il calcolo della curtosi
{
  m4<-mean((x-mean(x))^4)
  m2<-var(x)
  return(m4/m2^2)
}
curtosi(y) # >3
par(mfrow=c(1,1))
# Distribuzione centrata sulla media, presenza di simmetria e leptocurtosi
hist(y)
asimm(y)
qqnorm(y)
qqline(y) # Code più pesanti di quelle che si avrebbero sotto Normalità.

#### Stima del modello con gradiente ed hessiano numerici ####

fit.y.num<-optim(par=c(0.15,log(0.05),log(0.1),log(0.7)),lvgarch11,gr=NULL,method="BFGS",
                 control=list(fnscale=-1,maxit=200),hessian=T)
# Per la proprietà di equivarianza della stima di massima verosimiglianza:
theta.hat.num<-c(fit.y.num$par[1],exp(fit.y.num$par[2:4]))
theta.hat.num          # Valori dei parametri stimati 
fit.y.num$value        # Valore massimo della log-verosimiglianza
fit.y.num$convergence  # Se =0, l'algoritmo è arrivato a convergenza
# Stima della matrice di varianza e covarianza:
q<-diag(theta.hat.num)
j<-(solve(q)%*%(-fit.y.num$hessian)%*%solve(t(q)))/n
varAs<-solve(j)/n

# Visualizzazione dei risultati:
parametri<-c("mi","omega","alpha","beta")
theta.hat.num<-as.matrix(theta.hat.num)
dimnames(theta.hat.num)<-list(parametri)
Std.Err<-sqrt(diag(varAs))
tstat2<-theta.hat.num/Std.Err
pval<-2*pnorm(tstat2,lower.tail=F)
tabella_num<-data.frame(theta.hat.num,Std.Err,tstat2,pval)
colnames(tabella)_num<-c("Coefficienti","Errori standard","Statistica t","p-value")
print(tabella_num)
cat("Verosimiglianza in theta.hat=",fit.y.num$value,"\n")
cat("Codice di convergenza=",fit.y.num$convergence,"\n")

# Controllo:
# library(fGarch)
# fit = garchFit(~ garch(1, 1)+1,data=y,trace=F)
# fit

#### Stima del modello con gradiente analitico ed hessiano numerico ####

fit.y.an1<-optim(par=c(0.15,log(0.05),log(0.1),log(0.7)),lvgarch11,gr=grad11,method="BFGS",
                 control=list(fnscale=-1,maxit=200),hessian=T)
theta.hat.an1<-c(fit.y.an1$par[1],exp(fit.y.an1$par[2:4]))
theta.hat.an1          # Valori dei parametri stimati 
fit.y.an1$value        # Valore massimo della log-verosimiglianza
fit.y.an1$convergence  # Se =0, l'algoritmo è arrivato a convergenza
q<-diag(c(1,theta.hat.an1[2:4]))
varAn1<-solve(solve(q)%*%(-fit.y.num$hessian)%*%solve(t(q)))
varAn1                 # Stima della matrice di varianza e covarianza
# Visualizzazione dei risultati:
parametri<-c("mi","omega","alpha","beta")
theta.hat.an1<-as.matrix(theta.hat.an1)
dimnames(theta.hat.an1)<-list(parametri)
Std.Err<-sqrt(diag(varAn1))
tstat2<-theta.hat.an1/Std.Err
pval<-2*pnorm(tstat2,lower.tail=F)
tabella_anum<-data.frame(theta.hat.an1,Std.Err,tstat2,pval)
colnames(tabella_anum)<-c("Coefficienti","Errori standard","Statistica t","p-value")
print(tabella_anum)
cat("Verosimiglianza in theta.hat=",fit.y.num$value,"\n")
cat("Codice di convergenza=",fit.y.num$convergence,"\n")

#### Calcolo della matrice hessiana per via analitica (tramite la funzione "Hessian") ####

j<--hessian(theta.hat.an1)/n
varas.ml.ana<-solve(j)/n
parametri<-c("mi","omega","alpha","beta")
theta.hat.an1<-as.matrix(theta.hat.an1)
dimnames(theta.hat.an1)<-list(parametri)
Std.Err<-sqrt(diag(varas.ml.ana))
tstat2<-theta.hat.an1/Std.Err
pval<-2*pnorm(tstat2,lower.tail=F)
tabella_ana<-data.frame(theta.hat.an1,Std.Err,tstat2,pval)
colnames(tabella_ana)<-c("Coefficienti","Errori standard","Statistica t","p-value")
print(tabella_ana)
cat("Verosimiglianza in theta.hat=",fit.y.num$value,"\n")
cat("Codice di convergenza=",fit.y.num$convergence,"\n")

#Confronto tra tabelle relative ai tre metodi utilizzati
print(tabella_num)
print(tabella_anum)
print(tabella_ana)


########################## 
# Diagnostica del modello: analisi dei residui.
# per procedura numerica
h<-rep(0,n)
omega.hat<-theta.hat.num[2]
alpha.hat<-theta.hat.num[3]
beta.hat<-theta.hat.num[4]
h[1]<-omega.hat/(1-alpha.hat-beta.hat)
for (i in 2:n) { h[i]=omega.hat+alpha.hat*y2[i-1]+beta.hat*h[i-1]  }  # valori di h(t)
# Residui standardizzati:
res.std.num<-(y-theta.hat.num[1])/(h[2:n])
head(res.std.num)
hist(res.std.num)
summary(res.std.num)
par(mfrow=c(1,1))
qqnorm(res.std.num)
qqline(res.std.num)
par(mfrow=c(2,1))
cfs(res.std.num,15,new=T)
cfs(res.std.num^2,15,new=T)
curtosi(res.std.num) 
asimm(res.std.num) 
jarque.bera.test(x)

# Diagnostica del modello: analisi dei residui.
# per procedura analitica
# Calcolo la serie delle varianze condizionate:
res.std.ana<-y[2:n]/sqrt(h[2:n]) # Residui semplici/sigma.t
# Correlazione dei residui semplici al quadrato:
source("sse")
cfs(res.std.ana,15,new=T)
cfs(res.std.ana^2,15,new=T)
qqnorm(res.std.ana)
qqline(res.std.ana)
curtosi(res.std.ana) 
cat("Curtosi dei residui standardizzati:", curtosi(res.std.ana), "\n")

###################################
###################################
# -------- analisi diagnostiche -----------------
cfs(res.std,20,new=F)
cfs(res.std^2,20,new=F)
hist(res.std,20)
qqnorm(res.std)
qqline(res.std)
curtosi(res.std)
garch.stats(res.std, 10)

###################################
#################################
#################################
#################################
# Genero M=100 serie storiche di lunghezza n=500 per confrontare la precisione delle stime 
# nei diversi approcci (approccio numerico e analitico)

n<-500
m<-100
theta.hat.num<-matrix(0,ncol=4,nrow=m)
theta.hat.an<-matrix(0,ncol=4,nrow=m)
st.error.num<-matrix(0,ncol=4,nrow=m)
st.error.an<-matrix(0,ncol=4,nrow=m)
for(i in 1:m)
{
  set.seed(i)
  print(i)
  y<-garch11(n)[,1]
  y2<-y^2                  		
  
  fit.y.num<-optim(par=c(0.15,log(0.05),log(0.1),log(0.7)),lvgarch11,gr=NULL,method="BFGS",
                   control=list(fnscale=-1,maxit=200),hessian=T)
  q<-diag(c(1,exp(fit.y.num$par[2:4])))
  varAsN<-solve(solve(q)%*%(-fit.y.num$hessian)%*%solve(t(q)))
  varAsN
  
  fit.y.an<-optim(par=c(0.15,log(0.05),log(0.1),log(0.7)),lvgarch11,gr=grad11, method="BFGS",
                 control=list(fnscale=-1,maxit=200),hessian=T)
  q<-diag(c(1,exp(fit.y.an$par[2:4])))
  varAsAn<-solve(solve(q)%*%(-fit.y.an$hessian)%*%solve(t(q)))
  varAsAn 
  # j<--hessian(theta.hat.an1)/n
  # varas.ml.ana<-solve(j)/n
  
  theta.hat.num[i,]<-c(fit.y.num$par[1],exp(fit.y.num$par[2:4]))
  theta.hat.an[i,]<-c(fit.y.an$par[1],exp(fit.y.an$par[2:4]))
  st.error.num[i,]<-c(sqrt(diag(varAsN)))
  st.error.an[i,]<-c(sqrt(diag(varAsAn)))
}

head(st.error.an)
head(st.error.num)
head(theta.hat.num)
head(theta.hat.an)
nulli<-rep(0,4)
for (i in 1:4)
{
  nulli[i]<-sum(is.na(st.error.an[,i]))
}
nulli

# per n=500---> nulli: 0 0 0 0

col1<-colMeans(theta.hat.num)
col2<-colMeans(theta.hat.an)
col3<-apply(theta.hat.num,2,sd)
col4<-apply(theta.hat.an,2,sd)
tab1<-matrix(c(col1,col2,col3,col4),ncol=4,byrow=F)
colnames(tab1)<-c("Media Num","Media An","Dev Std Num","Dev Std An")
row.names(tab1)<-c("Mi","Omega","Alpha","Beta")
print(tab1)


# Faccio lo stesso con n=1000 e n=2000: 
# come cambia la differenza tra varianza stimata e varianza teorica 
# all'aumentare della lunghezza della serie? 
# Si ottengono varianze sempre più precise
# La varianza campionaria tende a stimare sempre meglio la varianza 
# marginale del processo.



########################################
########################################
########################################
# QML
########################################
########################################
########################################
gradt<-function(theta)
{
n=length(y)
mi=(theta[1])
omega=(theta[2])
alpha=(theta[3])
beta=(theta[4])
 e<-y-mi # Serie dei residui
  e2<-e^2
h=matrix(0,nrow=n,ncol=1)
h[1]=omega/(1-alpha-beta)
for (i in 2:n)
{
 h[i]=omega+alpha*y2[i-1]+beta*h[i-1]
}
h1=c(h[1],h[1:n-1])
c1<-rep(0,n-1)
  c2<-cumsum(beta^seq(0,(n-2),1))
  c3<-rep(0,n-1)
  c4<-rep(0,n-1)
  
  for (i in 2:n) 
  {
    c1[i-1]<--2*alpha*sum((beta^seq(0,i-2,1))*e[(i-1):1])
    c3[i-1]<-sum((beta^seq(0,i-2,1))*e2[(i-1):1]) 
  }
  for (i in 3:n) 
  {
    c4[i-1]<-omega*sum(seq(1,i-2,1)*(beta^seq(0,i-3,1)))+alpha*sum(seq(1,i-2,1)*(beta^seq(0,i-3,1))*e2[(i-2):1]) 
  }  
  
  dlvmi<--0.5*sum((1/h[2:n])*c1)+sum(e[2:n]/h[2:n])+0.5*sum((e2[2:n]/(h[2:n]^2))*c1) 
  dlvomega<--0.5*sum((1/h[2:n])*c2)+0.5*sum((e2[2:n]/(h[2:n]^2))*c2) 
  dlvalpha<--0.5*sum(c3/h[2:n])+0.5*sum(c3*e2[2:n]/(h[2:n]^2))
  dlvbeta<--0.5*sum(c4/h[2:n])+0.5*sum(c4*e2[2:n]/(h[2:n]^2))
  
  return(c(dlvmi,dlvomega,dlvalpha,dlvbeta)) 
}

# -------- calcolo di I con gradiente analitico --------
#calcolo di I con gradiente analitico (slide 13 file QML)
#assegniamo a dlt l'output di gradt calcolata su theta.hat #dlt calcolato per t che parte da 2, per quello dopo nel ciclo mettiamo i-1!!
dlt=gradt(theta.hat.an1)    
k=4
n<-length(y)
i0=matrix(0,nrow=k,ncol=k)   #Ã¨ la matrice identitÃ¨ I di zeri
for (i in 2:n)   
{
  i0=i0+(dlt[i-1,])%*%t(dlt[i-1,])  #sto sommando n-1 termini
}
i0=i0/n


# -------- Calcolo della matrice di covarianza con QML -------
# Calcolo della matrice hessiana per via analitica (tramite la funzione "Hessian"):
j<--hessian(theta.hat.an1)/n
varas.ml.ana<-solve(j)/n

varAs.QML=(solve(j) %*% i0 %*% solve(j))/n


#STIMA QML CON DERIVATE ANALITICHE
# Visualizzazione dei risultati:  
parametri<-c("mi","omega","alpha","beta")
theta.hat.an<-as.matrix(theta.hat.an1)
dimnames(theta.hat.an)<-list(parametri)
Std.Err.QML<-sqrt(diag(varAs.QML))
tstat2<-theta.hat.an1/Std.Err.QML
pval<-2*pnorm(abs(tstat2),lower.tail=F)
tabella<-data.frame(theta.hat.an1,Std.Err.QML,tstat2,pval)
colnames(tabella)<-c("Coefficienti","Errori standard","Statistica t","p-value")
print(tabella)
cat("Verosimiglianza in theta.hat=",fit.y.num$value,"\n")
cat("Codice di convergenza=",fit.y.num$convergence,"\n")

print(tabella_ana)
# Ottengo stesse stime puntuali ma diversi errori standard
# con l'utilizzo della QML potrebbero quindi cambiare anche la conclusioni inferenziali

########################################
########################################
# Che valore di mi utilizzare?
########################################
m<-rep(c(0.05,0.10,0.15,0.20))
#theta.hat.num<-matrix(0,nrow=4,ncol=4)
#theta.hat.ana<-matrix(0,nrow=4,ncol=4)
for(i in 1:4)
{
y=garch11(n=500,theta=c(m[i],0.05,0.1,0.7))
fit.y.ana<-optim(par=c(0.15,log(0.05),log(0.1),log(0.7)),lvgarch11,gr=grad11,method="BFGS",
                 control=list(fnscale=-1,maxit=200),hessian=T)
theta.hat.ana<-c(fit.y.ana$par[1],exp(fit.y.ana$par[2:4]))
fit.y.num<-optim(par=c(0.15,log(0.05),log(0.1),log(0.7)),lvgarch11,gr=NULL,method="BFGS",
                 control=list(fnscale=-1,maxit=200),hessian=T)
theta.hat.num<-c(fit.y.num$par[1],exp(fit.y.num$par[2:4]))
# Calcolo standard error numerici
q<-(diag(theta.hat.num))
j<-(solve(q)%*%(-fit.y.num$hessian)%*%solve((q)))/n
varas.num<-solve(j)/n
se.num<-sqrt(varas.num)
# Calcolo standard error analitici
j<--hessian(theta.hat.ana)/n
varas.ana<-solve(j)/n
se.ana<-sqrt(varas.ana)
# Stampa dei risultati
print(i)
print(c(fit.y.num$par[1],exp(fit.y.num$par[2:4])))
print(diag(se.num))
print(c(fit.y.ana$par[1],exp(fit.y.ana$par[2:4])))
print(diag(se.ana))
}

#print(theta.hat.num)
#print(theta.hat.ana)

#Confronto tra tabelle relative ai tre metodi utilizzati
print(tabella_num)
print(tabella_anum)
print(tabella_ana)

