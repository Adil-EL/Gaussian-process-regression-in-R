# Majeure Science des Données 2019-2020
# UP4 : exploitation de simulateurs numériques
# Gaussian process regression - Kriging
# Script du TP à compléter...

library(MASS) # fonction mvrnorm()

# Mouvement brownien

# Noyau brownien 
kBrown <- function(x,y,param=NULL){
  if(is.null(param)) param <-1
  param*outer(c(x),c(y),"pmin")
}
# Exemple d'utilisation de la fonction kBrown
X <- c(0.1,0.5,0.9)
kXX <- kBrown(X,X); print(kXX)
x <- 0.8
kxX <- kBrown(x,X); print(kxX)

# Simulation du mouvement brownien standard B

ngrid <- 400 # nombre de points de discrétisation
xgrid <- matrix(seq(from=0, to=1, length=ngrid),ncol=1) 
Kxgrid <- kBrown(xgrid,xgrid) 
muB <- 0*xgrid # moyenne du processus B
varB <- diag(Kxgrid) # variance de B

# Simulation avec mvrnorm = multivariate random normal
# attention Sigma = matrice de covariance

yB <- mvrnorm(n=5,mu=muB,Sigma=Kxgrid)
yB <- t(yB) # simulations en colonne

matplot(xgrid, yB, type='l', xlab="Temps t", ylab="B(t)", ylim=c(-3,3),
        main="Simulations du MB standard")
abline(h=0, lty=2)
lines(xgrid,1.96*sqrt(varB),lty=2,col="red")
lines(xgrid,-1.96*sqrt(varB),lty=2,col="red")

# Krigeage avec le noyau brownien

n <- 5
X <- matrix(seq(from=0.1, to=0.9, length=n), ncol=1)
yX <- matrix(c(0.5, 0, 2.5, 3, 2), ncol=1)
plot(X, yX, xlab="x", xlim=c(0,1), ylim=c(-0.5,3.5),ylab="y",
     main="f(x) connue en 6 points", cex.main=1)
points(0,0)
abline(h=0, lty=2)
abline(v=0, lty=2)

m <- 501
x <- matrix(seq(from=0, to=1, length=m),ncol=1)
kxx <- kBrown(x,x)
kxX <- kBrown(x,X)
kXX <- kBrown(X,X)
invkXX <- solve(kXX)
mu.krig <- kxX%*%invkXX%*%yX
K.krig <- kxx - kxX%*%invkXX%*%t(kxX)

y.krig <- mvrnorm(5,mu.krig,K.krig)

param <- 1
var.krig <- param*diag(K.krig)  
upp95 <- mu.krig + 1.96*sqrt(pmax(0,var.krig))
low95 <- mu.krig - 1.96*sqrt(pmax(0,var.krig))
par(mar=c(4.5,5.1,1.5,1.5))
plot(x, mu.krig, type="n", xlab="x",ylab="f(x) ?", ylim=range(low95-0.5,upp95+0.5),cex.axis=1,cex.lab=1)
lightblue <- rgb(114/255,159/255,207/255,.3) 
darkblue <- rgb(32/255,74/255,135/255,1)
darkbluetr <- rgb(32/255,74/255,135/255,.3)
polygon(c(x,rev(x)),c(upp95,rev(low95)),border=NA,col=lightblue)
lines(x,mu.krig,col=darkblue,lwd=3)
lines(x,low95,col=darkbluetr,lwd=2)  
lines(x,upp95,col=darkbluetr,lwd=2)
points(X, yX, pch=4, cex=1, lwd=3, col="red")
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(x, mu.krig, col="red",lwd=2)
for (simu in 1:5){
  lines(x, y.krig[simu,]) #col=darkbluetr)
}
title(main="Krigeage avec le noyau brownien",cex.main=1)



