cov<-read.csv("C:\\Users\\David\\Downloads\\cov.csv",header=TRUE)
str(cov)
covar<-as.matrix(cov[,2:6])
colnames(covar)<-c("VIX","VWO","TLT","IEI","IEF")
rownames(covar)<-c("VIX","VWO","TLT","IEI","IEF")

meanas<-read.csv("C:\\Users\\David\\Downloads\\means.csv",header=TRUE)
meanas$AssetClass<-c("VIX","VWO","TLT","IEI","IEF")

Rvals=seq(min(meanas$m)+0.1^10,max(meanas$m)-0.1^10,length.out=14);

library(MASS)

tscalar <- 0.25
C <- covar

var1 <- ginv(tscalar * C)


P <- matrix( c(1,0,-1,-1,-1,
               -1,1,0,0,0,
               0,-1,1,1,1),3,5)

Omega <- matrix(c(0.049563538,-0.058822285,0.009258747,
                  -0.058822285,0.075653644,-0.016831359,
                  0.009258747,-0.016831359,0.007572612),3,3)

var2 <- t(P) %*% ginv(Omega) %*% P
Q <- c(0.01,0.05,0.03)

var12 <- ginv(var1+var2)
var3 <- (t(P) %*% ginv(Omega) %*% Q) + (var1 %*% meanas$m) 
mhat <- var12 %*% var3
mhat

#Define the QP
Dmat1 <- 2*covar
dvec1 <- rep(0,5)
Amat1 <- matrix(c(mhat,-mhat,rep(1,5),rep(-1,5),diag(length(mhat))),5,9)

# compute efficient frontier for eight stocks - BL
BvarP=vector()
BsigmaP=vector()
Bw1=vector()
Bw2=vector()
Bw3=vector()
Bw4=vector()
Bw5=vector()


#Expected Returns 14 values
BRvals=seq(min(mhat)+0.1^10,max(mhat)-0.1^10,length.out=14);


for (i in 1:length(Rvals)) {
  BR=BRvals[i]
  bvec1 <- c(BR,-BR,1,-1,0,0,0,0,0)
  BqpSol=solve.QP(Dmat1,dvec1,Amat1,bvec1)
  BvarP[i]=BqpSol$value
  BsigmaP[i]=sqrt(BvarP[i])
  Bw1[i]=BqpSol$solution[1];
  Bw2[i]=BqpSol$solution[2];
  Bw3[i]=BqpSol$solution[3];
  Bw4[i]=BqpSol$solution[4];
  Bw5[i]=BqpSol$solution[5];
}

#Portfolio weights
Bweightsoutput<-data.frame(Bw1,Bw2,Bw3,Bw4,Bw5)
#Checking whether summation is equal to one
weightsum<-apply(Bweightsoutput,1,sum)

#Efficient frontier Plot - BL
par(mfrow = c(1,1))
BLplot <- plot(BsigmaP,BRvals,type = 'l',lty = 1,lwd=3, xlab = 'Risk',ylab = 'Returns', main = 'Black litterman method',col = 'blue')

#Weights of portfolio assets vs expected returns in one plot - BL
par(mfrow = c(1,1))
plot(BRvals,Bw1,type = 'l', lty = 1,lwd=3, xlab ='Returns' ,ylab = 'Weights of Assets' , main = 'Visualization of Asset Weights vs Returns', col = 'red')
lines(BRvals,Bw2,"l",lty = 1,lwd=3,col = 'black')
lines(BRvals,Bw3,"l",lty = 1,lwd=3,col='green')
lines(BRvals,Bw4,"l",lty = 1,lwd=3,col = 'blue')
lines(BRvals,Bw5,"l",lty = 1,lwd=3,col='violet')

legend('topleft', c("VIX","VWO","TLT","IEI","IEF"), pch = 17, 
       col = c('red','black','green','blue','violet'), text.col = c('red','black','green','blue','violet','dark green','violet','violet'), cex = .6)

#Weights of portfolio assets vs expected returns in separate plots - BL
par(mfrow = c(2,3))
plot(BRvals,Bw1,type = 'l', lty = 1,lwd=3, xlab ='Returns' ,ylab = 'Weight' , main = 'VIX index', col = 'red')
plot(BRvals,Bw2,type = 'l', lty = 1,lwd=3, xlab ='Returns' ,ylab = 'Weight' , main = "Vanguard FTSE Emerging Markets ETF (VWO)", col = 'black')
plot(BRvals,Bw3,type = 'l', lty = 1,lwd=3, xlab ='Returns' ,ylab = 'Weight' , main = "iShares 20+ Year Treasury Bond ETF", col = 'green')
plot(BRvals,Bw4,type = 'l', lty = 1,lwd=3, xlab ='Returns' ,ylab = 'Weight' , main = "iShares 3-7 Year Treasury Bond ETF", col = 'blue')
plot(BRvals,Bw5,type = 'l', lty = 1,lwd=3, xlab ='Returns' ,ylab = 'Weight' , main = "iShares 7-10 Year Treasury Bond ETF", col = 'violet')
