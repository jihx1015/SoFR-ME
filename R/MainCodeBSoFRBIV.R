#--- This contains all the functions needed for estimation
#-------------------------------------------------------------------
require(LaplacesDemon)
require(GeneralizedHyperbolic)
require(mvtnorm)
require(splines)
require(nimble)
require(Matrix)
require(pracma)
require(truncnorm)
require(tmvtnorm)
require(emulator)
require(mclust)
require(Rcpp)
require(compiler)
require(TruncatedNormal)
require(inline)
require(data.table)
require(refund)
require(mixtools)
require(KScorrect)
require(magrittr)

#-------------------------------------------------------------------
#--- Convert C++ functions to help speed up sampling
#-------------------------------------------------------------------
{
  
  {Code <-
    '
    arma::mat CikW(arma::mat W, arma::mat X, arma::mat Mu, arma::mat Sig, arma::vec Pk){

    const double log2pi = std::log(2.0 * M_PI);

    int n = W.n_rows;
    int K = Mu.n_rows;
    int P = Mu.n_cols;
    double constants = -(static_cast<double>(P)/2.0) * log2pi;
    arma::mat sigma(P, P);
    arma::mat pik(n, K);
    arma::mat rooti(P, P);

    for(int k=0; k < K; k++){
    sigma = arma::reshape(Sig.row(k), P, P);
    rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    for(int i=0; i < n; i++){
    arma::rowvec xtp = W.row(i) - X.row(i);
    arma::vec z = rooti * arma::trans( xtp - Mu.row(k));
    pik(i,k) = Pk(k)*exp(constants - 0.5 * arma::sum(z%z) + rootisum);
    }
    }
    return(pik);
    }
    '
  }
  
  {Code1a <- '
    arma::mat CikM(arma::mat M, arma::mat X, arma::mat Mu, arma::mat Sig, arma::vec Pk, arma::mat bs2 ){

    const double log2pi = std::log(2.0 * M_PI);

    int n = M.n_rows;
    int K = Mu.n_rows;
    int P = M.n_cols;
    int pn = bs2.n_cols;
    double constants = -(static_cast<double>(P)/2.0) * log2pi;
    arma::mat sigma(pn, pn);
    arma::mat sigmaB(P, P);
    arma::mat pik(n, K);
    arma::mat rooti(P, P);

    for(int k=0; k < K; k++){
    sigma = arma::reshape(Sig.row(k), pn, pn);
    sigmaB = arma::mat (bs2 * sigma * bs2.t() + 0.0001 * arma::eye(P, P)); 
    rooti = arma::trans(arma::inv(trimatu(arma::chol(sigmaB))));
    double rootisum = arma::sum(log(rooti.diag()));
    for(int i=0; i < n; i++){
    arma::rowvec xtp = M.row(i) - X.row(i);
    arma::vec z = rooti * arma::trans( xtp - Mu.row(k));
    pik(i,k) = Pk(k)*exp(constants - 0.5 * arma::sum(z%z) + rootisum);
    }
    }
    return(pik);
    }
    '
  }
  
  
  {Code1 <- '
    arma::mat CikMvec(arma::mat M, arma::mat X, arma::mat Mu, arma::mat Sig, arma::vec Pk, arma::vec Delt){

    const double log2pi = std::log(2.0 * M_PI);

    int n = M.n_rows;
    int K = Mu.n_rows;
    int P = Mu.n_cols;
    double constants = -(static_cast<double>(P)/2.0) * log2pi;
    arma::mat sigma(P, P);
    arma::mat pik(n, K);
    arma::mat rooti(P, P);
    arma::mat diDelt = diagmat(Delt);

    for(int k=0; k < K; k++){
    sigma = arma::reshape(Sig.row(k), P, P);
    rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    for(int i=0; i < n; i++){
    arma::rowvec xtp = M.row(i) - X.row(i)*diDelt;
    arma::vec z = rooti * arma::trans( xtp - Mu.row(k));
    pik(i,k) = Pk(k)*exp(constants - 0.5 * arma::sum(z%z) + rootisum);
    }
    }
    return(pik);
    }
    '
  }
  
  {Code2 <- '
    arma::mat CikX(arma::mat X, arma::mat Mu, arma::mat Sig, arma::vec Pk){

    const double log2pi = std::log(2.0 * M_PI);

    int n = X.n_rows;
    int K = Mu.n_rows;
    int P = Mu.n_cols;
    double constants = -(static_cast<double>(P)/2.0) * log2pi;
    arma::mat sigma(P, P);
    arma::mat pik(n, K);
    arma::mat rooti(P, P);

    for(int k=0; k < K; k++){
    sigma = arma::reshape(Sig.row(k), P, P);
    rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    for(int i=0; i < n; i++){
    arma::rowvec xtp = X.row(i);
    arma::vec z = rooti * arma::trans( xtp - Mu.row(k));
    pik(i,k) = Pk(k)*exp(constants - 0.5 * arma::sum(z%z) + rootisum);
    }
    }
    return(pik);
    }
    '
  }
  #
  {Code3 <- '

    arma::colvec Y = Rcpp::as<arma::vec>(Ys);
    arma::mat X = Rcpp::as<arma::mat>(Xs);
    arma::colvec  Mu = Rcpp::as<arma::vec>(Mus);
    arma::colvec Sig = Rcpp::as<arma::vec>(Sigs);
    arma::colvec Pk  = Rcpp::as<arma::vec>(Pks);
    arma::colvec Gam = Rcpp::as<arma::vec>(Gams);


    int n = X.n_rows;
    int K = Mu.n_elem;

    arma::mat pik = arma::zeros<arma::mat>(n, K);
    arma::colvec Ytld = Y - X*Gam;

    for(int k=0; k < K; k++){
    for(int i=0; i < n; i++){
    pik(i,k) = Pk(k)*R::dnorm(Ytld(i), Mu(k), sqrt(Sig(k)), 0);
    }
    }
    return Rcpp::wrap(pik);
    '
  }
  #
  Cike <- cxxfunction(signature(Ys="numeric", Xs="numeric",Mus="numeric",Sigs="numeric",Pks="numeric",Gams="numeric"),
                      Code3, plugin="RcppArmadillo")
  
  cppFunction(Code,depends = "RcppArmadillo")
  cppFunction(Code1,depends = "RcppArmadillo")
  cppFunction(Code1a,depends = "RcppArmadillo")
  cppFunction(Code2,depends = "RcppArmadillo")
  
  {
    Code4 <- '

  arma::mat UpdateX(arma::vec Yd, arma::mat Wd, arma::mat Md, arma::mat Mukx, arma::vec Sige, arma::mat Sigw, arma::mat Sigm, arma::mat Sigx,
  arma::vec Cike, arma::vec Cikw, arma::vec Cikm, arma::vec Cikx,arma::vec Gamma, double Delta, arma::vec A, arma::vec B){
  int n = Yd.n_elem;
  int J = Wd.n_cols;
  arma::mat Sigt(J,J);
  arma::vec Mut;
  arma::mat Gamsq = Gamma*Gamma.t();
  double Deltsq = Delta*Delta;
  arma::mat Res(n,J);

  // Obtaining namespace of Matrix package
  Environment pkg1 = Environment::namespace_env("LaplacesDemon");
  Function F1 = pkg1["as.symmetric.matrix"];
  Function F2 = pkg1["as.inverse"];
  Environment pkg2 = Environment::namespace_env("TruncatedNormal");
  Function F3 = pkg2["rtmvnorm"];

  for(int i=0; i < n; i++){
  int wi = Cikw(i)-1;
  int mi = Cikm(i)-1;
  int xi = Cikx(i)-1;
  int ei = Cike(i)-1;
  arma::mat Sigw0 = as<arma::mat>(F2(F1(reshape(Sigw.row(wi), J, J))));
  arma::mat Sigm0 = as<arma::mat>(F2(F1(reshape(Sigm.row(mi), J, J))));
  arma::mat Sigx0 = as<arma::mat>(F2(F1(reshape(Sigx.row(xi), J, J))));

  Sigt = as<arma::mat>(F2(F1((Gamsq/Sige(ei) + Sigw0 + Sigx0 + Gamsq*Sigm0))));
  Mut = Sigt*((Gamma*Yd[i]/Sige(ei)) + Sigw0*Wd.row(i).t() + Sigm0*Md.row(i).t() + Sigx0*Mukx.row(xi).t());

  Res.row(i) = as<arma::rowvec>(F3(1, Mut, Sigt, A, B));

  }

  return Res;

  }
  '
  }
  cppFunction(Code4, depends = "RcppArmadillo")
  
}



#--- AR1 covariance function
#'@ n = size of the correlation matrix
#'@ rho = correlation value

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
#--- Distributon of errors for epsilon
#--  Options are: 
#'@ n = sample size
#'@ distF = "Rst", "Norm","MixNorm", "Gam"
#'@ sig.e = error standard deviation
DistEps <- function(n, distF, sig.e){
  pi.g <- .5
  ei = switch(distF,
              "Gam" = rgamma(n=n,shape=1,scale=1.5) - 1.5,
              "Rst"   = sn::rst(n = n, xi = 0, omega = .5, nu = 5, alpha = 2),
              "Norm" = rnorm(n=n, sd = sig.e),
  )
  ei
}

#--- Distributon of errors for U (the measurement error)
#'@ n = sample size 
#'@ distF = "Rmst" and "Mvnorm"
#'@ Sig is a positive define covariance matrix

DistUw <- function(n,distF, Sig){
  p0 <- nrow(Sig)
  p10 = round(0.5*p0); p11 <- p0 - p10
  
  ei = switch(distF,
              "Rmst" =  sn::rmst(n=n, xi = c(rep(-2,p10), rep(2,p11)), Omega = Sig, alpha = rep(2,p0), nu= 5),
              "Mvnorm" = mvtnorm::rmvnorm(n=n,mean = rep(0,nrow(Sig)), sigma = Sig)
  )
  scale(ei, center = T, scale = F)
}

#---- Cov matrix type
#'@ m denotes the size of the covariance matrix 
#'@ Met = "CS" (compound symmetry) or "Ar1" (Ar 1 type covariance)
#'@ sig = standard deviation 
#'@ rho =  correlation value
CovMat <- function(m, Met,sig,rho){
  
  ei = switch(Met,
              "CS" =   (ones(m,m) - diag(rep(1,m)))*((sig^2)*rho) + diag(rep(sig^2,m)),
              "Ar1" = ar1_cor(m, rho)*(sig^2) #CVTuningCov::AR1 
  )
  ei
}


#---  Beta function
#'@ t = scalar/vector of time point at which function is evaluated (0, 1)
#'@ met = c("f1","f2","f3","f4","f5")
Betafunc <- function(t, met){
  switch(met,
         f1 = 1.*sin(2*pi*t) + 1.5 ,
         f2 = 1/(1+exp(4*2*(t-.5))),
         f3 = 1.*sin(pi*(8*(t-.5))/2)/ (1 + (2*(8*(t-.5))^2)*(sign(8*(t-.5))+1)) + 1.5,
         f4 = sin(pi*(16*(t-.5))/2)/ (1 + (2*(16*(t-.5))^2)*(.5*sin(8*(t-.5))+1)),
         f5 = sin(pi*(16*(t-.5))/2)/ (1 + (2*(16*(t-.5))^2)*(sin(8*(t-.5))+1))
  )
}

#'@ t = scalar/vector of time point at which function is evaluated (0, 1)
#'@ met = c("f1","f2","f3","f4","f5")

BetafuncDel <- function(t, met){
  switch(met,
         f1 = sin(2*pi*t),
         f2 = 1/(1+exp(4*2*(t-.5))),
         f3 = sin(pi*(8*(t-.5))/2)/ (1 + (2*(8*(t-.5))^2)*(sign(8*(t-.5))+1)),
         f4 = sin(pi*(16*(t-.5))/2)/ (1 + (2*(16*(t-.5))^2)*(.5*sin(8*(t-.5))+1)),
         f5 = sin(pi*(16*(t-.5))/2)/ (1 + (2*(16*(t-.5))^2)*(sin(8*(t-.5))+1))
  )
}

#'@ function deal with values outside of the prescribed range during MCMC

smFun <- function(Mu, Lw, Uw){
  i0 = sum(c((Mu < Lw) + (Mu > Uw)))
  
  i0*(Lw*(Mu < Lw) + Uw*(Mu > Uw)) + (1 - i0)*Mu
}

#'@ a = denote the time points at whcih the function will be evaluated
#'@ Delta =  scalar scaling value 
#'@ met = c("f1","f2","f3","f4","f5")
Deltafun <- function(a, Delta, met){ Delta*(0.5 + 0.5*BetafuncDel(a, met)) + min(0.005, Delta) + 0.2 }


#-- Simulate data with constant variance
#'@ n = sample size
#'@ t = number of time points
#'@ j0 = number of error free covariates
#'@ rho.x = correlation in the process for X(t)
#'@ sig.x = standard deviation for X(t) at each time point
#'@ rho.u = correlation in the process for U(t) (measurement error)
#'@ sig.u = standard deviation for U(t) at each time point
#'@ rho.m = correlation in the process for omega(t) (error for the IV)
#'@ sig.w = standard deviation for Omega(t) at each time point (error for the IV)
#'@ sig.e =. error variance for the response Y
#'@ Delta : scalling value for the function delta(t)
#'@ distFep = c( "Norm","Rst","MixNorm") : distribution of the response Y 
#'@ distFUe = c("Mvnorm","Rmst") : distribution for the ME (U(t))
#'@ distFepWe = c("Mvnorm", "Rmst") : distribution for the IV error (Omega(t))  
#'@ idf : c(1:4) what function we want for the mean of X(t) (see the definition of `Betafunc` above)


#n = 500; sig.e=1; sig.x = 4; sig.u = 4; sig.w=1; rhox = 0.25; rho.u = 0.25; rho.m = 0.25;Estdetl=FALSE ; Delta=2; t = 50; pn = ceiling(n^{1/3.8})+1; idf = 1; distF <- "Norm"; distU = distU; distW=distW;Gx = 3}

DataSimFunc <- function(n = 300, t = 50, j0 = 2, rho.x = 0.25, sig.x = 4, rho.u = 0.25, sig.u = 4,
                        rho.m = 0.25, sig.w = 1, sig.e =1, Delta = 2,  distFep=c( "Norm","Rst","MixNorm")[1], distFUe = c("Mvnorm","Rmst")[1],
                        distFepWe=c("Mvnorm", "Rmst")[1], CovMet =c("CS", "Ar1")[1], idf=c(1:4)[1]){
  
  a <- seq(0, 1, length.out=t) ## time points at which obs. are takes
  
  met <- c("f1","f2","f3","f4","f5")
  
  P.mat <- function(K){
    # penalty matrix
    D <- diag(rep(1,K))
    D <- diff(diff(D))
    P <- t(D)%*%D 
    return(P)
  }
  
  #---------------------------------------------------------------------
  #---- Simulate X
  #---------------------------------------------------------------------
  #sig.x = .5
  SigX <- CovMat(t, CovMet,sig.x,rho.x) 
  MeanX <- Betafunc(a, met[idf]) 
  MeanX2 <- Betafunc(a, met[idf]) 
  X_t <- matrix(MeanX, ncol=t, nrow=n, byrow = T) + DistUw(n,distFUe, SigX)  
  plot(a,MeanX, type="l",lwd=3)
  
  
  #'@ Simulate the error free covariates
  Z <- matrix(rnorm(n=n*j0, sd = 1), ncol = j0); Z <- scale(Z); 
  
  #'@ Betaz denotes the effects for the error free covariates
  Betaz <- c( -1.05, 0.57)
  
  
  # We generate Yi the scalar response with error terms from Normal (0,.05) (Equation (1))  
  ei <- DistEps(n, distFep, sig.e)
  BetatF <- (Betafunc(a,met[idf]) - mean(Betafunc(a,met[idf])))
  fx0 = crossprod(t(X_t),BetatF)/length(a)
  
  Y <- Z%*%Betaz + fx0  + ei 
  
  
  #---------------------------------------------------------------------
  # We generate the observed surrogates values W 
  #---------------------------------------------------------------------
  SigU <- CovMat(t, CovMet,sig.u,rho.u) 
  MeanU <- numeric(t)
  U <- DistUw(n,distFUe, SigU) 
  W = X_t + U                 
  #---------------------------------------------------------------------
  # We generate the IV with errors from the Normal distribution (Equation (3))
  #---------------------------------------------------------------------
  SigW <- ones(t,t)*rho.m*sig.w*sig.w; diag(SigW) <- (sig.w^2) 
  MeanW <- numeric(t)
  omega <- DistUw(n,distFUe, SigW) #mvtnorm::rmvnorm(n=n,mean = MeanW, sigma = SigW)
  
  #--- Simulate delta(t)
  Delt = Deltafun(a, Delta, met[idf])
  
  
  M <- X_t%*%diag(Delt) + omega  # (n,t) M matrix
  
  Data <- list(a = a, Z = Z,  Y = Y, X = X_t, 
               Data=list(M= M, W = W, X = X_t), Delt = Delt, BetaZ = Betaz, Betat = BetatF)
  Data
}

#--- Input
#'@ Posterior sample analysis
#'@ SumRes  = output from  the main function BayesFIVSing
#'@ T0 = number of distinct time points
#'@  n = sample size from the data
#'@ dfm =. degree of smoothing to be passed to the `smooth.spline` function
#'@ Burn = How much burn in 
#'@ probs = probability for the posterior interval 

#'@ Output
#'
#'

PosteriorSum <- function(SumRes, T0 = 50, dfm = 8, Burn = NULL, probs = 0.95, Data, lwd=.5 ){

  if(is.null(Burn)){
    Burn = round(c(1:(0.3*nrow(SumRes$Gamma))))
  }
pn <- ncol(SumRes$Gamma)
a <- seq(0, 1, length.out = T0)
bs2 <- bs(a, df = pn, intercept = T)
bs2del <- qr.Q(qr(bs2)) 

GamPost <- SumRes$Gamma[-c(Burn),] %*% t(bs2)
GamPostSm <- matrix(NA, ncol = ncol(GamPost), nrow = nrow(GamPost))

DeltPost <- SumRes$Betdelt[-c(Burn),] %*% t(bs2del)

DeltPostSm <- matrix(NA, ncol = ncol(DeltPost), nrow = nrow(DeltPost))

for(i in 1: nrow(GamPost)){
  GamPostSm[i,] <- smooth.spline(a, GamPost[i,], df = dfm)$y
  DeltPostSm[i,] <- smooth.spline(a, DeltPost[i,], df = dfm)$y
}

BetaPost <- as.data.frame(rbind(colMeans(SumRes$BetaZ[-c(Burn),]),
                  round(apply(SumRes$BetaZ[-c(Burn),], 2, quantile, probs = 1 - (1 - probs)*0.5), 3),
                  round(apply(SumRes$BetaZ[-c(Burn),], 2, quantile, probs = (1 - probs)*0.5), 3)
))
colnames(BetaPost) <- paste0("Beta",1:ncol(SumRes$BetaZ))
rownames(BetaPost) <- paste0("Post_",c("Mean",round(1 - (1 - probs)*0.5, 5), round( (1 - probs)*0.5,3)))

Datgam <- data.frame(x = a, y = colMeans(GamPostSm), upy = apply(GamPostSm, 2, quantile, probs = 0.975) , 
                     lwy = apply(GamPostSm, 2, quantile, probs = 0.025),
                     ytrue = Data$Betat)

#Res <- list()
P1 = ggplot(Datgam, aes(x = x, y = y)) + geom_ribbon(aes(ymin = lwy, ymax = upy), alpha = 0.2) + geom_line(linewidth=lwd)
P11 <- P1 + geom_line(aes(x = x, y = ytrue, color = "#6BD7AF"),linewidth=lwd)
print(P11)

Datdel <- data.frame(x = a, y = colMeans(DeltPostSm), upy = apply(DeltPostSm, 2, quantile, probs = 0.975) , 
                     lwy = apply(DeltPostSm, 2, quantile, probs = 0.025),
                     ytrue = Data$Delt)

P2 = ggplot(Datdel, aes(x = x, y = y)) + geom_ribbon(aes(ymin = lwy, ymax = upy), alpha = 0.2) + geom_line(linewidth=lwd)
P22 <- P2 +  geom_line(aes(x = x, y = ytrue, color = "#6BD7AF"),linewidth=lwd)
print(P22)

# Res[[3]] = list(GamSm  = GamPostSm,GamRaw = GamPost ,
#      DeltatSm = DeltPostSm, DeltatRaw = DeltPost,
#      BetaZ = BetaPost
#      )
Res =  list(BetaZ = BetaPost)
Res

}




#@parms Dat: Is a list
#@parms pn: Number of basis used
SumeRiseX <- function(Tempdat, pn, bs2, burn0){
  Burn <- c(1:burn0)
  Res <- list()
  #-- Collect things related to X
  #-------------------------------------------------------------------
  #-- Mean 
  Res$Mu1.x <- (Tempdat$Muk.x[-Burn,1:pn])%*%t(bs2)*pn ; 
  Res$Mu2.x <- (Tempdat$Muk.x[-Burn,pn+c(1:pn)])%*%t(bs2)*pn; 
  Res$Mu3.x <- (Tempdat$Muk.x[-Burn,2*pn+c(1:pn)])%*%t(bs2)*pn; 
  
  #-- Covariance matrix
  CovFunc<- function(x, pn){matrix(c(x),nrow=pn, ncol=pn, byrow = T)}
  Res$Sig1.x <- Reduce("+", apply(Tempdat$Sigk.x[-Burn,1:(pn*pn)],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.x[-Burn,1:(pn*pn)])
  Res$Sig2.x <- Reduce("+", apply(Tempdat$Sigk.x[-Burn, pn*pn+c(1:(pn*pn))],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.x[-Burn,1:(pn*pn)])
  Res$Sig3.x <- Reduce("+", apply(Tempdat$Sigk.x[-Burn, (2*pn*pn)+c(1:(pn*pn))],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.x[-Burn,1:(pn*pn)])
  
  #-- Collect things related to W
  #-------------------------------------------------------------------
  #-- Mean 
  Res$Mu1.w <- (Tempdat$Muk.w[-Burn,1:pn])%*%t(bs2)*pn ; 
  Res$Mu2.w <- (Tempdat$Muk.w[-Burn,pn+c(1:pn)])%*%t(bs2)*pn; 
  Res$Mu3.w <- (Tempdat$Muk.w[-Burn,2*pn+c(1:pn)])%*%t(bs2)*pn; 
  #-- Covariance matrix
  CovFunc<- function(x, pn){matrix(c(x),nrow=pn, ncol=pn, byrow = T)}
  Res$Sig1.w <- Reduce("+", apply(Tempdat$Sigk.w[-Burn,1:(pn*pn)],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.w[-Burn,1:(pn*pn)])
  Res$Sig2.w <- Reduce("+", apply(Tempdat$Sigk.w[-Burn, pn*pn+c(1:(pn*pn))],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.w[-Burn,1:(pn*pn)])
  Res$Sig3.w <- Reduce("+", apply(Tempdat$Sigk.w[-Burn, (2*pn*pn)+c(1:(pn*pn))],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.w[-Burn,1:(pn*pn)])
  
  #-- Collect things related to M
  #-------------------------------------------------------------------
  #-- Mean 
  Res$Mu1.m <- (Tempdat$Muk.m[-Burn,1:pn])%*%t(bs2)*pn ; 
  Res$Mu2.m <- (Tempdat$Muk.m[-Burn,pn+c(1:pn)])%*%t(bs2)*pn; 
  Res$Mu3.m <- (Tempdat$Muk.m[-Burn,2*pn+c(1:pn)])%*%t(bs2)*pn; 
  #-- Covariance matrix
  CovFunc<- function(x, pn){matrix(c(x),nrow=pn, ncol=pn, byrow = T)}
  Res$Sig1.m <- Reduce("+", apply(Tempdat$Sigk.m[-Burn,1:(pn*pn)],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.w[-Burn,1:(pn*pn)])
  Res$Sig2.m <- Reduce("+", apply(Tempdat$Sigk.m[-Burn, pn*pn+c(1:(pn*pn))],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.w[-Burn,1:(pn*pn)])
  Res$Sig3.m <- Reduce("+", apply(Tempdat$Sigk.m[-Burn, (2*pn*pn)+c(1:(pn*pn))],1, CovFunc, pn=pn, simplify = F))/nrow(Tempdat$Sigk.w[-Burn,1:(pn*pn)])
  
  Res
  
}

#--- Obtain the counts/frequency of the each cluster visit for each individual
#-------------------------------------------------------------------------------
CountHr <- function(Vec,G){
  fx <- function(x,g) 1*(x==g)
  fx <- Vectorize(fx,"g")
  colMeans(fx(Vec,1:G))
}

CountHr(sample(c(1:4),size=10,replace=T),G=4)







#---- Main function to fit our approach
#'@ Nsim = number of MCMC samples
#'@ Y = a vector of length n (sample size)
#'@ W = matrix of n (rows) and t (columns) [surrogate of X]
#'@ M = matrix of n (rows) and t (columns) [ IV values of X]
#'@ Z : matrix of n (rows) and j0 (columns) for the error free 
#'@ sig02.del, = scaling factor for P-Spline (can be set of estimated)
#'@ alpx = 1  Dirichlet concentration parameter set at 1
#'@ Ge = chosen number of normal mixture needed to approximate the distribution of the response Y (integer 1 or above)
#'@ Gx = chosen number of normal mixture needed to approximate the distribution of the response X (integer 1 or above)
#'@ Gmw = chosen number of normal mixture needed to approximate the distribution of the response W and M (integer 1 or above)
#'@ BasFun = what basis function to choose from.Options are: c("Bspline", "Four")
#'@ PsiAccept = postive Scaling factor to improve MH performance 
#'@ sig2.propdelt  =  variance for the proposal distribution for delta
#'@ sig2.priordelt =  variance for the prior distribution for delta

#===============================================================================
#==== Exactly on component of the normal distribution
BayesFIVSing = function(Nsim = 100, Y, W, M, Z, sig02.del = NULL,  
                          alpx = 1, BasFun = "Bspline", PsiAccept = 0.05, 
                          sig2.propdelt = 2.75, sig2.priordelt = 5.0){
  #--- Functions
  func <- function(Pik, n=1, size=1){
    which.max(rmultinom(n=1,size=1,prob = Pik + abs(min(Pik))+ 0.001))
  }
  funcCmp <- cmpfun(func)
  
  F2 <- function(Pik, n=1, size=1){
    
    #print(Pik)
    apply(Pik,1,funcCmp)
  }
  F2Cmp <- cmpfun(F2)
  
  a <-  seq(0, 1, length.out = ncol(M))
  n = length(Y)
  cat("--- set some initial values ")
  #---- Parameters settings
  #-----------------------------------------------------------------------------
  pn <- round(n^{1/3}) + 2 ## min(15, ceiling(length(Y)^{1/3.8})+4)
  
  alpha.e = alpx
  alpha.w = alpha.W = alpx
  alpha.x = alpha.X = alpx
  alpha.m = alpx
  pn = pn 
  n = length(Y) 
  Zmod = Z
  
  #-----------------------------------------------------------------------------
  #--- Transform the data
  #-----------------------------------------------------------------------------
  Ge = 2
  bs2 <- bs(a, df = pn, intercept = T)
  bs2del <- qr.Q(qr(bs2))
  
  #--- Using different bases
  
  if(BasFun == "Four"){
    
    #Nbasis0 = c(3,4, 5, 7) + 2
    J = K = k  = pn = pn00 =  Nbasis0[which(c(100,200,500, 1000) == n)]
    #Nbasis0 =  c(3, 6, 15, 18)
    J = K = k  = pn = pn00 =  Nbasis0[which(c(100,200,500, 1000) == n)]
    
    
    Fourmeth <- fda::create.fourier.basis(c(0,1), nbasis = pn)
    bs2 = as.matrix(eval.basis(a, Fourmeth))
    colnames(bs2) <- NULL
    bs2del <- qr.Q(qr(bs2)) 
    
  }
  
  Mn <- (M%*%bs2)/length(a) #(M%*%bs2)/length(a) 
  Wn <- (W%*%bs2)/length(a)
  
  #--- Penalty matrix
  
  P.mat <- function(K){
    # penalty matrix
    D <- diag(rep(1,K))
    D <- diff(diff(D))
    P <- t(D)%*%D 
    return(P)
  }
  
  #--- Estimate initial values
  InitialValues <- function(Y, W, M, Z, df0=3, a, Ge, bs2){
    Gx = 3
    Res <- list()
    #bs.w <- splines::bs(a, df = pn+1, degree = min(c(df0,pn-1)), intercept=TRUE) 
    Zmod <- Z #cbind(1,Z)
    pn = ncol(bs2)
    
    
    bs2del <- qr.Q(qr(bs2))
    
    Wn <- crossprod(t(W), bs2)/length(a)
    Mn <- crossprod(t(M), bs2)/length(a)
    
    WsmB <- Wn
    WsM <- NULL
    for(i in 1:nrow(W)){ WsM <- rbind(WsM, smooth.spline(a, W[i,])$y)}
    W_i <- crossprod(t(WsM),bs2)/length(a)   #(n,k) matrix
    lmFit <- lm(Y ~ -1 + Zmod + W_i) 
    c_hatW <- lmFit$coefficients
    Res$Gamma <- as.numeric(c_hatW[-c(1:ncol(Zmod))])
    Res$BetaZ <- as.numeric(c_hatW[c(1:ncol(Zmod))])
    Res$Delta <- abs(mean(Mn)/mean(W_i))
    #Res$Betadel <- c(abs(solve(t(bs2del)%*%bs2del)%*%t(bs2del)%*%(colMeans(M)/colMeans(W))))
    
    #--- Estimate Delta
    Delhat <- function(bs, U0){
      
      Fe0 <- function(val){
        sum((bs%*%val - U0)^2)
      }
      
      k = ncol(bs)
      Res <- optim( rep(0, ncol(bs)), Fe0, method="L-BFGS-B", lower=rep(-60,k), upper=rep(60, k) )
      Res  
      
    }
    U0 = colMeans(M)/colMeans(W)
    R0d <- Delhat(bs2del, U0)
    
    Res$Betadel <- R0d$par 
    
    #---------------------------------------------------------------------------
    #--- Determine the number of clusters for epsilon_i
    #---------------------------------------------------------------------------
    
    reslmfit <- residuals.lm(lmFit)
    BIC <- mclustBIC(reslmfit)
    mod1 <- Mclust(reslmfit, x = BIC, G = Ge)
    Temp <- summary(mod1, parameters = TRUE)
    Res$Ke <- Temp$G #+ Gx
    Res$muk.e <- numeric(Res$Ke)
    Res$muk.e[1:Temp$G] = Temp$mean; 
    Res$sig2k.e <- rep(1, Res$Ke)*sum(lmFit$residuals^2)/lmFit$df.residual
    Res$sig2k.e[1:Temp$G] = Temp$variance
    Res$Cik.e <- as.numeric(Temp$classification)
    Res$Pik.e <- numeric(Res$Ke)
    Res$Pik.e[1:Temp$G] <- Temp$pro 
    #---------------------------------------------------------------------------
    #--- Determine the number of cluster (W)
    #---------------------------------------------------------------------------
    Ui <- Wn #- W_i 
    Mod1 <- Mclust(Ui, G = Gx)
    Tp <- summary(Mod1, parameters=T)
    Res$Kw <- Tp$G #+Gx
    Res$Muk.w <- matrix(0, ncol = pn, nrow=Res$Kw)
    Res$Muk.w <- as.matrix( rbind(t(Tp$mean))) #  ,  repmat(rep(0,pn), Gx, 1)))
    
    Res$Sigk.w <- matrix(0, nrow = Res$Kw, ncol = pn*pn)
    for(j in 1:Tp$G){Res$Sigk.w[j,] <- c(Tp$variance[,,j])}
    
    Res$Cik.w <- Tp$classification
    Res$Pik.w <- Tp$pro
    
    #---------------------------------------------------------------------------  
    #--- Determine the number of cluster (X)
    #---------------------------------------------------------------------------
    WsmB <- (W + M) /(1+Res$Delta)
    WsM <- NULL
    for(i in 1:nrow(W)){ WsM <- rbind(WsM, smooth.spline(a, WsmB[i,])$y)}
    W_i <- crossprod(t(WsM),bs2)/length(a)   #(n,k) matrix
    Ui <- W_i 
    Mod1 <- Mclust(Ui, G=3)
    Tp <- summary(Mod1, parameters=T)
    Res$Kx <- Tp$G #+ Gx
    Res$Muk.x <- matrix(0, ncol = pn, nrow=Res$Kw)
    Res$Muk.x <- as.matrix( rbind(t(Tp$mean),  repmat(colMeans(Ui), Gx, 1 )))
    
    Res$Sigk.x <- matrix(0, nrow = Res$Kx, ncol = pn*pn)
    for(j in 1:Tp$G){Res$Sigk.x[j,] <- c(Tp$variance[,,j])}
    #for(j in (Tp$G+1):Res$Kx){Res$Sigk.x[j,] <- c(.5*diag(pn))}
    Res$Cik.x <- Tp$classification
    Res$Pik.x <- Tp$pro  
    Res$X <- W_i
    #---------------------------------------------------------------------------
    #---- Determine the number of cluster (M)
    #---------------------------------------------------------------------------
    
    Ui <- Mn  #- Res$Delta*W_i 
    Mod1 <- Mclust(Ui, G=Gx)
    Tp <- summary(Mod1, parameters=T)
    Res$Km <- Tp$G #+ Gx
    #Res$Muk.m <- as.matrix(t(Tp$mean))
    Res$Muk.m <- as.matrix( rbind(t(Tp$mean))) #,  repmat(colMeans(Ui), Gx, 1 )))
    Res$Sigk.m <- matrix(0,nrow=Res$Km,ncol=pn*pn)
    for(j in 1:Tp$G){Res$Sigk.m[j,] <- c(Tp$variance[,,j])}
    #for(j in (Tp$G+1):Res$Km){Res$Sigk.m[j,] <- c(as.positive.definite(.5*diag(pn)))}
    Res$Cik.m <- Tp$classification
    Res$Pik.m <- Tp$pro
    Res$R0d <- R0d
    #------- return result 
    Res 
  }
  
  
  cat("\n #--- Obtain initial values...\n")
  
  Thet0 <- InitialValues(Y, W, M, Z, df0=3, a, Ge, bs2)
  
  #-----------------------------------------------------------------------------
  cat("\n#---- Initialize parameters ---> \n")
  #-----------------------------------------------------------------------------
  
  #Kx = Thet0$Kx               ## number of X  cluster
  Ke = Thet0$Ke               ## number of ei cluster
  #Km = Thet0$Km               ## number of ei cluster
  #Kw = Thet0$Kw               ## number of omega_i cluster 
  
  
  Mn0 <- (M%*%bs2)/length(a) #(M%*%bs2)/length(a) 
  Wn <- (W%*%bs2)/length(a)
  Mn <- ((M%*%diag(1/(c(bs2del%*%Thet0$Betadel))))%*%bs2)/length(a) 
  
  Thet <- list()
  
  Thet$Gamma <-  rnorm(n=length(Thet0$Gamma)) #as.numeric(Thet0$Gamma)   #- a vector of length J
  Thet$BetaZ <- Thet0$BetaZ
  Thet$siggam <-  ifelse(is.null(sig02.del), 1, sig02.del)
  Thet$bs2del <- bs2del
  Thet$Betadelt <- Thet0$Betadel
  
  
  Thet$Cik.x = rep(1, n) #Thet0$Cik.x #sample(x = 1:Kx,size = n, replace = T)    #- vector of length Kx
  Thet$Cik.e = rep(1, n) #Thet0$Cik.e #sample(x = 1:Ke,size = n, replace = T)     #- vector of length Ku
  Thet$Cik.w = rep(1, n)  #Thet0$Cik.w  #sample(x = 1:Kw,size = n, replace = T)     #- vector of length Kw
  Thet$Cik.m = rep(1, n)  #Thet0$Cik.m #sample(x = 1:Km,size = n, replace = T)     #- vector of length Km
  
  
  Thet$Pik.x = c(1,0) #rdirch(n=1,alpha = rep(alpha.X, Kx)/Kx)   #- vector of length Kx
  Thet$Pik.e = c(1,0) #rdirch(n=1,alpha = rep(alpha.e, Ke)/Ke)  #- vector of length Ku
  Thet$Pik.w = c(1,0) #rdirch(n=1,alpha = rep(alpha.W, Kw)/Kw)  #- vector of length Kw
  Thet$Pik.m = c(1,0) #rdirch(n=1,alpha = rep(alpha.m, Km)/Km)  #- vector of length Kw
  
  
  Thet$sig2k.e  =  Thet0$sig2k.e[1]  #rgamma(n=Ke, shape=1,scale=.1) # matrix dim Kx x J
  Thet$muk.e  =  0#Thet0$muk.e #rnorm(n=Ke)  #- Matrix of (Ku x J) and Ku x J - Note that $\Omega$ is a diagonal matrix
  
  Thet$Sigk.w <- .5*diag(pn)     #- vector of length Kw
  Thet$Muk.w <- rep(0,pn) #c(rnorm(n = pn))         #- vectors of length Kw and Kw respectively
  
  Thet$Sigk.x <- .5*diag(pn)     #-- vector of length Kx
  Thet$Muk.x  <- rep(0,pn) #rnorm(n=pn)   #- vectors of length Kx and Kx respectively
  
  Thet$Sigk.m <- .5*diag(pn)     #-- vector of length Kx
  Thet$Muk.m  <- rep(0,pn) #c(rnorm(n = pn))   #- vectors of length Kx and Kx respectively
  
  
  A <- min(c(c(Wn), c(Mn))) - 0.05*diff(range(c(Wn, Mn)))
  B <- max(c(c(Wn), c(Mn))) + 0.05*diff(range(c(Wn, Mn)))
  
  Thet$X <- TruncatedNormal::rtmvnorm(n=n,mu=rep(0,pn), sigma=diag(pn),lb = rep(A, pn), ub = rep(B, pn))

  
  cat("\n#---- Loading updating functions ---> \n")
  #-----------------------------------------------------------------------------
  # Update things related to Y
  #-----------------------------------------------------------------------------
  #-------- Update Gamma (vector of length pn + the Z together) 
  #--- Prior parms
  pn0 = ncol(Z)
  Pgam <- P.mat(pn)*1.0
  Sig0.gam0 = diag(pn0); Mu0.gam0 = numeric(pn0+pn);
  
  updateGam <- function(Thet, Y, Z, Sig0.gam0, Mu0.gam0){
    
    pn = ncol(Thet$X)
    Zmod <-  Z #cbind(1,Z)
    Y.td = Y - Thet$muk.e[Thet$Cik.e]
    InvSig0.gam0 <- 1.0*as.matrix(Matrix::bdiag(as.inverse(Sig0.gam0), P.mat(pn)/Thet$siggam))
    Xmod <- cbind(Zmod, Thet$X)
    X_sc <- sweep(Xmod, MARGIN = 1, Thet$sig2k.e[Thet$Cik.e], FUN="/") ## n*pn
    
    Sig.gam <- as.inverse(as.symmetric.matrix(crossprod(X_sc, Xmod) + InvSig0.gam0))
    Mu.gam <- c(Sig.gam%*%(crossprod(X_sc, Y.td) + InvSig0.gam0%*%Mu0.gam0))
    
    c(mvtnorm::rmvnorm(n=1, mean = Mu.gam, sigma = as.positive.definite(Sig.gam)))
    
  }
  updateGam <- cmpfun(updateGam)
  
  #Update siggma - siggam ~ IG(al0, bet0)
  al0 = 1
  bet0 = .005
  updateSigGam <- function(Bet, al0, bet0, order = 2){
    K = length(Bet)
    al = al0 + .5*(K - order)
    bet = .5*quad.form(P.mat(K), Bet) + bet0
    #1/rgamma(n = 1, shape = al, rate = bet)
    #nimble::rinvgamma(n = 1, shape = al, scale = bet)
    # tb = 10
    # 1/heavy::rtgamma(n=1,shape=al,scale = bet, t = tb)
    tb = 1
    1/cascsim::rtgamma(n=1,shape=al, scale = 1/bet, max = 10, min = 1/15)
    #tb = .01
    #ifelse( qgamma(.95, shape= al, rate = bet) > tb, 1/tb , 1/cascsim::rtgamma(n=1,shape=al, scale = 1/bet, max = tb))
  }
  
  #------- Update Pik.e 
  #alpha.e = 1.0
  updatePik.e <- function(Thet, alpha.e){
    
    nk <- numeric(length(Thet$muk.e))
    for(i in 1:length(Thet$muk.e)){ nk[i] = sum(Thet$Cik.e == i)}
    
    #update the Pis
    alpha.enew <- alpha.e/length(Thet$muk.e) + nk
    
    rdirichlet(n=1, alpha.enew)
    
  }
  
  # #--- Update muk.w and sig2k.w based on the mixture of normal prior Xtl
  # mu0.e = 0 ; sig20.e = 100
  # updatemuk.e <- function(Thet, Y, Z, mu0.e, sig20.e ){
  #   
  #   Ku = length(Thet$muk.e)
  #   muk <- numeric(Ku)
  #   sig2k <- numeric(Ku)
  #   Zmod <- Z #cbind(1,Z)
  #   #Ytl <- Y - Thet$X%*%Thet$Gamma
  #   Ytl = Y - Zmod%*%Thet$BetaZ -  Thet$X%*%Thet$Gamma
  #   for(k in 1:length(Thet$muk.e)){
  #     id <- which(Thet$Cik.e == k)
  #     nk = length(id)
  #     if(nk > 0){
  #       sig2k[k] <- 1/(nk/Thet$sig2k.e[k] + 1/sig20.e)
  #       muk[k] <- sig2k[k]*(mu0.e/sig20.e + sum(c(Ytl[id]))/Thet$sig2k.e[k]) 
  #     }else{
  #       sig2k[k] <- sig20.e
  #       muk[k] <- mu0.e 
  #     }
  #   }
  #   ## Force the zero mean constrains
  #   SigR0 <- Thet$Pik.e*sig2k
  #   SigRR <- sum((Thet$Pik.e^2)*sig2k)
  #   MuK0 = muk - SigR0*(1/SigRR)*sum(Thet$Pik.e*muk)
  #   SigK <- as.matrix(nearPD(diag(sig2k) - tcrossprod(SigR0)*(1/SigRR), doSym = T)$mat)
  #   
  #   if(Ku > 2){
  #     id0 = 1:(Ku-1)
  #     Muk_1 <- c(mvtnorm::rmvnorm(n=1,mean=MuK0[-Ku], sigma = SigK[id0,id0]))
  #     
  #     c(Muk_1, -sum(c(Thet$Pik.e[-Ku]*Muk_1)/ifelse(Thet$Pik.e[Ku] == 0.0, (.Machine$double.eps)^{.7}, Thet$Pik.e[Ku])))
  #   }else{
  #     
  #     id0 = 1:(Ku-1)
  #     Muk_1 <- c(rnorm(n=1,mean=MuK0[-Ku], sd = sqrt(SigK[id0,id0])))
  #     
  #     c(Muk_1, -(c(Thet$Pik.e[-Ku]*Muk_1)/ifelse(Thet$Pik.e[Ku] == 0.0, (.Machine$double.eps)^{.7}, Thet$Pik.e[Ku])))
  #   }
  #   
  # }
  # 
  # 
  
  gam0.e = 1 ; sig20.esig = 1
  #gam0.e = 1 ; sig20.e = .5
  updatesig2k.e <- function(Thet, Y, Z, gam0.e, sig20.e){
    
    Zmod <- Z #cbind(1,Z)
    Ku = length(Thet$muk.e)
    al <- numeric(Ku)
    bet <- numeric(Ku)
    Xtl <- Y - Thet$X%*%Thet$Gamma - Zmod%*%Thet$BetaZ  #- Thet$muk.e[Thet$Cik.e]
    for(k in 1:length(Thet$muk.e)){
      id <- which(Thet$Cik.e == k)
      nk = length(id)
      if(nk > 0){
        al[k] <- nk/2 + gam0.e
        bet[k] <- sig20.e + .5*sum(Xtl[id]^2) 
      }else{
        al[k] <- gam0.e
        bet[k] <- sig20.e 
      }
    }
    ## Force the zero mean constrains
    mapply(rinvgamma, n = 1, shape = al, scale = bet)
  }
  
  #-----------------------------------------------------------------------------
  # Update things related to W 
  #-----------------------------------------------------------------------------
  
  #--- Update Pik.x  
  #alpha.W = alpha.w = 1.0
  updatePik.w <- function(Thet, alpha.W){
    Ku = nrow(Thet$Muk.w)
    nk <- numeric(Ku)
    for(k in 1:Ku){ nk[k] = sum(Thet$Cik.w == k)}
    
    #update the Pis
    alpha.wnew <- alpha.W/length(Thet$muk.w) + nk
    
    rdirichlet(n=1, alpha.wnew)
    
  }
  
  #--- Update Ci.k mixture of normals 
  dmvnormVec <- function(x, Mean, SigmaVc){
    Ku <- nrow(Mean)
    Val <- numeric(Ku)
    for(k in 1:Ku){ Val[k] <- mvtnorm::dmvnorm(x = x, mean = c(Mean[k,]), sigma = as.symmetric.matrix(matrix(SigmaVc[k,],ncol=length(x),byrow=F))) }
    Val
  }
  
  updateCik.w <- function(i, Thet, W){ # mixture of normals 
    #Sig <- matrix(Thet$Sigk.w[id,], ncol = ncol(W), byrow = F)
    #kw = nrow(Thet$Muk.w)
    xi =  c(W[i,] - Thet$X[i,])
    #pi.k  = Thet$Pik.w*mapply(dmvnormVec, k = 1:Kw, MoreArgs = list(x = xi, Mean = Thet$Muk.w , SigmaVc = Thet$Sigk.w))
    pi.k  = Thet$Pik.w*dmvnormVec(x=xi, Mean = Thet$Muk.w , SigmaVc = Thet$Sigk.w) + .Machine$double.eps
    which(rmultinom(n=1,size=1,prob = pi.k) == 1)
  }
  
  #--- Update Cik.w 
  #--- Update muk.w and sig2k.w based on the mixture of normal prior
  Mu0.w = rep(0, pn);  Sig0.w = .5*diag(pn);inSig0.w = 2*diag(pn) # as.inverse(.5*cov(Wn)) #.5*diag(pn) #
  updateMuk.w <- function(Thet, W, Mu0.w, inSig0.w, Sig0.w){
    Ku = length(Thet$Pik.w)
    J = ncol(Thet$X)  
    Mu.k <- matrix(0, nrow = Ku, ncol = J) # K x J
    SigK <- list()  #-- store the covariances
    SigKWg <- matrix(0, ncol=J, nrow=J) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    Yi.tld = W - Thet$X 
    #Muk <- matrix(0, ncol=J,nrow=Ku)
    SigK <- list()
    for(k in 1:Ku){
      idk <- which(Thet$Cik.w == k)
      nk <- length(idk)
      if(length(idk) > 0){
        Yi0 = matrix(c(Yi.tld[idk,]), nrow= nk, byrow=F)
        TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.w[k,],ncol=J,byrow=F)))
        SigK[[k]] <- as.inverse(as.symmetric.matrix(nk*TpSi + inSig0.w)) ## sum of matrices inverse
        SigKWg <- SigKWg + (Thet$Pik.w[k]*Thet$Pik.w[k])*SigK[[k]] # J x J
        #Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi.tld[idk,])) ) + solve(Sig0.w)%*%Mu0.w)      ## some of mat
        Mu.k[k,] <- c(SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi0)) ) + inSig0.w%*%Mu0.w))      ## some of mat
        #Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(apply(t(Yi.tld[idk,]), 2, sum)) ) + solve(Sig0.w)%*%Mu0.w)      ## some of mat
      }else{
        SigK[[k]] <- Sig0.w
        SigKWg <- SigKWg + (Thet$Pik.w[k]*Thet$Pik.w[k])*SigK[[k]] # J x J
        Mu.k[k,] <- Mu0.w
      }
    }
    #--- we impose some constrainsts on the mean
    SIg0 <- bdiag(SigK) # create a block diagonal matrix
    SIGR0 <- NULL; for(k in 1:Ku){SIGR0 <- rbind(SIGR0,Thet$Pik.w[k]*SigK[[k]])} ## K*J x J
    MUR0 <-  colSums(diag(Thet$Pik.w) %*% Mu.k) # return a evector of length J
    
    InvSigKWg <- as.positive.definite(as.inverse(SigKWg))
    # print(dim(InvSigKWg))
    # print(dim(SIg0))
    # print(dim(quad.tform.inv(SigKWg,SIGR0)))
    # print(Ku)
    #MURFin <- c(t(Mu.k)) - SIGR0%*%solve(SigKWg)%*%MUR0 ## J*K
    MURFin <- c(t(Mu.k)) - SIGR0%*%InvSigKWg%*%MUR0 ## J*K
    #SIGR0fin <- 1.0*as.matrix(nearPD(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T)$mat) ## J*K x J*K
    #SIGR0fin <- 1.0*(as.positive.definite(as.symmetric.matrix(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T))) ## J*K x J*K
    SIGR0fin <- 1.0*(as.positive.definite(as.symmetric.matrix(as.matrix(SIg0 - SIGR0%*%InvSigKWg%*%t(SIGR0)))))
    #--- Simulate K-1 vectors from the degenerate dist.
    id <- 1:(J*(Ku-1))
    SimMU.k <- matrix(c(mvtnorm::rmvnorm(n=1, mean = MURFin[id], sigma = SIGR0fin[id,id])), ncol=J,byrow=T) ## (K-1) x J
    SimMU.k1 <- -colSums(diag(Thet$Pik.w[-Ku]) %*% SimMU.k)/Thet$Pik.w[Ku] ## get a K-1 x J matrix follow by colSums - a vector of J
    as.matrix(rbind(SimMU.k,SimMU.k1)) ## k * J
  }
  updateMuk.w <- cmpfun(updateMuk.w)
  
  #--- update Sigk.w
  nu0.w = pn + 2; Psi0.w = .5*diag(pn) 
  updateSigk.w <- function(Thet, W, nu0.w, Psi0.w){
    
    Ku = nrow(Thet$Muk.w)
    m = ncol(W)
    nk = nrow(W)
    
    Wtl <- W - Thet$X 
    nun <- nu0.w  + nk
    #Yi0 <- matrix(c(Wtl[id,]), nrow= nk, byrow=F)
    #Psin.w <- as.positive.definite(Psi0.w + crossprod(Yi0))
    Psin.w <- as.positive.definite(Psi0.w + crossprod(Wtl))
    rinvwishart(nun, Psin.w)
    
  }
  updateSigk.w <- cmpfun(updateSigk.w)
  
  #-----------------------------------------------------------------------------
  # Update things related to M 
  #-----------------------------------------------------------------------------
  
  #--- update Pik.m for epsilon_i
  #alpha.m = 1.0
  updatePik.m <- function(Thet, alpha.m){
    nk <- numeric(nrow(Thet$Muk.m))
    for(k in 1:nrow(Thet$Muk.m)){ nk[k] = sum(Thet$Cik.m == k)}
    
    #update the Pis
    alpha.mnew <- alpha.m/nrow(Thet$Muk.m) + nk
    
    c(rdirichlet(n=1, alpha.mnew))
  }
  
  #--- update Cik.m for epsilon_i based on mixture of normals
  updateCik.m <- function(i, Thet, M){ # mixture of normals 
    #id <- Thet$Cik.m[i]
    #Sig <- matrix(Thet$Sigk.m[id,], ncol = ncol(M), byrow = F)
    Km = nrow(Thet$Muk.m)
    xi = M[i,] - Thet$Delta*Thet$X[i,]
    #pi.k  = Thet$Pik.m*mapply(dmvnormVec, k = 1:Km, MoreArgs = list(x = xi , Mean = Thet$Muk.m , SigmaVc = Thet$Sigk.m))
    pi.k  = Thet$Pik.m*dmvnormVec(x = xi, Mean = Thet$Muk.m , SigmaVc = Thet$Sigk.m) + .Machine$double.eps
    which(rmultinom(n=1,size=1,prob = pi.k) == 1)
  }
  
  
  #--- Update Muk.w and sig2k.w based on the mixture of normal prior
  Mu0.m = rep(0,pn) ; Sig0.m = .5*diag(pn); inSig0.m = 2*diag(pn) 
  updateMuk.m <- function(Thet, M, Mu0.m, inSig0.m, Sig0.m){
    #cat(dim(inSig0.m),"(--1)\n")
    Ku = length(Thet$Pik.m)
    J = ncol(M)  
    Mu.k <- matrix(0, nrow = Ku, ncol = J) # K x J
    SigK <- list()  #-- store the covariances
    SigKWg <- matrix(0, ncol = J, nrow = J) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    Yi.tld = M - Thet$X%*%diag(Thet$Betadel) 
    Muk <- matrix(0, ncol=J, nrow=Ku)
    for(k in 1:Ku){
      idk <- which(Thet$Cik.m == k)
      nk <- length(idk)
      if(length(idk) > 0){
        Yi0 = matrix(c(Yi.tld[idk,]), nrow = nk, byrow=F)
        TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.m[k,],ncol=J,byrow=F)))
        SigK[[k]] <- as.inverse(as.symmetric.matrix(nk*TpSi + inSig0.m)) ## sum of matrices inverse
        #cat(dim(SigK[[k]]),"(--2 \n")
        SigKWg <- SigKWg + (Thet$Pik.m[k]*Thet$Pik.m[k])*SigK[[k]] # J x J
        Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi0)) ) + inSig0.m%*%Mu0.m)      ## some of mat
        #cat(dim(SigKWg),"(--3 \n")
      }else{
        SigK[[k]] <- Sig0.m
        SigKWg <- SigKWg + (Thet$Pik.m[k]*Thet$Pik.m[k])*SigK[[k]] # J x J
        Mu.k[k,] <- Mu0.m
      }
    }
    #--- we impose some constrainsts on the mean
    InvSigKWg <- as.positive.definite(as.inverse(SigKWg))
    SIg0 <- bdiag(SigK) # create a block diagonal matrix
    SIGR0 <- NULL;  for(k in 1:Ku){SIGR0 <- rbind(SIGR0,Thet$Pik.m[k]*SigK[[k]])} ## K*J x J
    MUR0 <-  colSums(diag(Thet$Pik.m) %*% Mu.k) # return a evector of length J
    
    # print(Ku)
    # print(dim(SIg0))
    # print(dim(SIGR0))
    #MURFin <- c(t(Mu.k)) - SIGR0%*%solve(SigKWg)%*%MUR0 ## J*K
    #SIGR0fin <- as.matrix(nearPD(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T)$mat) ## J*K x J*K
    
    MURFin <- c(t(Mu.k)) - SIGR0%*%InvSigKWg%*%MUR0 ## J*K
    SIGR0fin <- as.symmetric.matrix(as.positive.definite(as.matrix(SIg0 - SIGR0%*%InvSigKWg%*%t(SIGR0)))) ## J*K x J*K
    
    #--- Simulate K-1 vectors from the degenerate dist.
    id <- 1:(J*(Ku-1))
    SimMU.k <- matrix(c(mvtnorm::rmvnorm(n=1, mean = MURFin[id], SIGR0fin[id,id])), ncol=J,byrow=T) ## (K-1) x J
    SimMU.k1 <- -colSums(diag(Thet$Pik.m[-Ku]) %*% SimMU.k)/Thet$Pik.m[Ku] ## get a K-1 x J matrix follow by colSums - a vector of J
    as.matrix(rbind(SimMU.k,SimMU.k1)) ## k * J
  }
  updateMuk.m <- cmpfun(updateMuk.m)
  
  #--- Update Sig2k.w based on the mixture prior
  nu0.m = pn + 2; Psi0.m = .5*diag(pn)#as.symmetric.matrix(.25*cov(Wn))
  updateSigk.m <- function(Thet, M, nu0.m, Psi0.m, bs2=bs2){
    
    Ku = nrow(Thet$Muk.m)
    m = ncol(Thet$X)
    nk = nrow(M) 
    
    Mtl <- (M - (tcrossprod(Thet$X, bs2) * matrix(c(Thet$bs2del%*%Thet$Betadel), nrow=nk, ncol = ncol(M), byrow=T)))%*%bs2 
    
    nun <- nu0.m  + nk
    Psin.m <- as.positive.definite(Psi0.m + crossprod(Mtl))
    rinvwishart(nun, Psin.m)
  }
  updateSigk.m <- cmpfun(updateSigk.m)
  
  #---- Update things related to delta(t)
  #---- Function to Compute  hij matrix (N x T)
  LogSimDElVec <- function(M, Mn, X, Sigma, Val=1){
    #LogSimDElVec <- function(M, Mn, X, Sigma, Val){
    id <- Val
    #mvtnorm::dmvnorm(M[i,], mean = Mn[i,], sigma = Sigma, log = T)
    #LaplacesDemon::dmvnc(M[i,], mu = Mn[i,], U = Sigma, log = T)
    LaplacesDemon::dmvnc(M, mu = Mn, U = Sigma, log = T)
  }
  
  #-- Makes the hij (N x T) matrix
  Hij <- function(Betadel, Thet, bs2){
    
    (Thet$X%*% t(bs2)) * matrix(c(Thet$bs2del%*%matrix(Betadel,ncol=1)), nrow = nrow(Thet$X), ncol = nrow(bs2), byrow = T) ## Matrix NxK
  }
  Hij <- cmpfun(Hij)
  
  
  
  #'@ Muvec0 : prior mean vector
  #'@ Sig0 : prior variance
  #'@ likelihood + prior
  #--- Second version (without projecting delta(t)*X_{i}(t))
  LikHijMat <- function(Betadel, Thet, M, Muvec0, Sig0, bs2, Val=1){
    Sig = as.positive.definite(as.symmetric.matrix(quad.tform(Thet$Sigk.m, bs2)))
    Sigchol <- chol(Sig)
    MeanV = Hij(Betadel, Thet, bs2)
    #sum(LogSimDElVec(1:nrow(M), M, MeanV, Thet$X, Thet$Sigk.m, Val)) + mvtnorm::dmvnorm(Betadel, mean = Muvec0, sigma = Sig0, log = T)
    sum(LogSimDElVec(M, MeanV, Thet$X, Sigchol, Val=1)) + mvtnorm::dmvnorm(Betadel, mean = Muvec0, sigma = Sig0, log = T)
  }
  
  # 
  #--- Without projecting the function to lower space.
  LikHijMatOptim <- function(Thet, M, Muvec, Sig0, bs2, Val=1){
    #Thet$Betadel
    #MeanV = ((Thet$X%*% t(bs2) * matrix(c(bs2%*%matrix(Betadel,ncol=1)), nrow = nrow(Thet$X), ncol = nrow(bs2), byrow = T))%*%bs2)/nrow(bs2) ## Matrix NxK
    J = length(Muvec)
    
    #lb = Muvec - 0.5*abs(Muvec)
    #ub = Muvec + 0.5*abs(Muvec)
    
    parm0 <- c(rnorm(n = J))
    SigU = chol(as.positive.definite(quad.tform(Thet$Sigk.m, bs2)))
    #print(dim(SigU))
    
    funC <- function(parm){
      MeanV = Hij(parm, Thet, bs2)
      -sum(LogSimDElVec(M, MeanV, Thet$X, SigU, Val=1)) - mvtnorm::dmvnorm(parm, mean = Muvec,sigma = Sig0, log = T)
    }
    
    Res <- optim(parm0, funC, method="L-BFGS-B", lower=rep(-15,J), upper = rep(15,J), hessian = T)
    #Res <- optim(parm0, funC, method="L-BFGS-B", lower=c(rep(lb,J), 0.01), upper = c(rep(ub,J), 10.0), hessian = T)
    Res
  } 
  
  #--- Updated on 10/16/23 (update Deltat)
  #@' M : matrix of n (samples) x T(number of distinct time points)
  #@' Sigp0 = proposal covariance
  #@' Muvec : prior mean vector
  #@' Sig0 : prior variance
  #@' bs2 : are the basis function
  #mu0.del, Sig0.del
  #Muvec = c(t(bs2)%*%colMeans(M)/colMeans(W))
  #Muvec = mu0.del = rep(0, ncol(Mn)); Sig0 = Sig0.del = diag(rep(.3,ncol(Mn))); Sigp0 = diag(rep(.2,ncol(Mn)))
  #'@ Find the covariance matrix for the  proposal distribution
  
  Mu0dl = rep(0, ncol(Mn))
  Muvec = mu0.del = rep(0, ncol(Mn)); Sig0 = Sig0.del = diag(rep(sig2.priordelt,ncol(Mn))); Sigp0 = diag(rep(.2,ncol(Mn)))
  Muvec = mu0.del = Thet0$Betadel #c(t(bs2) %*% matrix(c(colMeans(M)/colMeans(W)), ncol=1))
  
  #--- expanded version 
  updateFunDeltaStNew <- function(Thet, M, Mu0dl, Muvec, Sig0, bs2, Sigp0){
    Mitl <- M 
    #-- Current value
    BetadelCur <- Thet$Betadelt  
    #-- Propose a value
    #BetadelNew <- c(mvtnorm::rmvnorm(n=1, mean = BetadelCur, sigma = Sigp0)) 
    BetadelNew <- c(mvtnorm::rmvnorm(n=1, mean = Muvec, sigma = Sigp0)) # simulate a proposal from the prior
    
    #-- Compute the likelihood (proposal)
    ernew <- LikHijMat(BetadelNew, Thet, Mitl, Mu0dl, Sig0, bs2,1) - mvtnorm::dmvnorm(BetadelNew, mean = Muvec, sigma = Sigp0, log=T)
    
    #-- Compute the likelihood (old)
    erold <- LikHijMat(BetadelCur, Thet, Mitl, Mu0dl, Sig0, bs2,1) - mvtnorm::dmvnorm(BetadelCur, mean = Muvec, sigma = Sigp0, log = T)
    
    # #-- Compute the likelihood (proposal)
    # ernew <- LikHijMat(BetadelNew, Thet, Mitl, Muvec, Sig0, bs2, Thet$Cik.m) - mvtnorm::dmvnorm(BetadelNew, mean = Muvec, sigma = Sigp0, log=T)
    # #-- Compute the likelihood (old)
    # erold <- LikHijMat(BetadelCur, Thet, Mitl, Muvec, Sig0, bs2, Thet$Cik.m) - mvtnorm::dmvnorm(BetadelCur, mean = Muvec, sigma = Sigp0, log = T)
    
    u0 <- runif(n=1)
    #--- 
    #Fval <- ifelse(log(u0) - (ernew-erold) < 0, ernew, erold)
    accept0 <- 1*(log(u0) < (ernew - erold))
    
    list(res=c(accept0, c(BetadelNew*(accept0) + (1 - accept0)*BetadelCur)), prop = BetadelNew, evid = ernew - erold)
  }
  updateFunDeltaStNew <- cmpfun(updateFunDeltaStNew)
  
  
  #-----------------------------------------------------------------------------
  # Update things related to X 
  #-----------------------------------------------------------------------------
  #--- update Pik.m for epsilon_i
  #alpha.x <- 1.0
  updatePik.x <- function(Thet, alpha.x){
    nk <- numeric(nrow(Thet$Muk.x))
    for(k in 1:nrow(Thet$Muk.x)){ nk[k] = sum(Thet$Cik.x == k)}
    
    #update the Pis
    alpha.xnew <- alpha.x/nrow(Thet$Muk.x) + nk
    
    c(rdirichlet(n=1, alpha.xnew))
  }
  
  #--- update Cik.m for epsilon_i based on mixture of normals
  updateCik.x <- function(i, Thet){ # mixture of normals 
    #Kx = ncol(Thet$X) # J x J
    #id <- Thet$Cik.x[i]
    #Sig <- matrix(Thet$Sigk.x[id,], ncol = length(Thet$Gamma), byrow = F)
    Xi = Thet$X[i,]
    #pi.k  = Thet$Pik.x*mapply(dmvnormVec, k = 1:Kx, x = Xi,  Mean = Thet$Muk.x , SigmaVc = Thet$Sigk.x)
    pi.k  = Thet$Pik.x*dmvnormVec(x=Xi,Mean = Thet$Muk.x , SigmaVc = Thet$Sigk.x) + .Machine$double.eps
    which(rmultinom(n=1,size=1,prob = pi.k) == 1)
  }
  #updateCik.x <- Vectorize(updateCik.x,"i")
  
  #--- Update muk.w based on the mixture of normal prior
  Mu0.x = numeric(pn); Sig0.x = .5*diag(pn); inSig0.x = 2*diag(pn)#as.inverse(.25*cov(Wn)) #.5*diag(pn) #
  updateMuk.x <- function(Thet, Mu0.x, inSig0.x, Sig0.x){
    K = length(Thet$Pik.x)
    J = ncol(Thet$X)
    Res0 <- matrix(0, nrow = K, ncol = J)
    Mu.k <- matrix(0, nrow = K, ncol = J) # K x J
    nk = nrow(Thet$X)
    #SigKWg <- matrix(0, ncol = ncol(Y),nrow=ncol(Y)) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    
    TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.x, ncol=J, byrow=F)))
    SigK <- as.inverse(as.symmetric.matrix((nk*TpSi + inSig0.x))) ## sum of matrices inverse
    Muk <- SigK%*%(crossprod(TpSi, c(colSums(Thet$X))) + inSig0.x%*%Mu0.x)      ## some of mat
    c(mvtnorm::rmvnorm(n=1, mean = Muk, sigma = SigK))
    
    #--- we impose some constrainsts on the mean
    
    
  }
  updateMuk.x <- cmpfun(updateMuk.x)
  
  #--- Update Sig2k.w based on the mixture of normal prior
  nu0.x <- pn + 2 ; Psi0.x <- .5*diag(pn) #as.symmetric.matrix(.25*cov(Wn))
  updateSigk.x <- function(Thet, nu0.x, Psi0.x){
    
    Ku = nrow(Thet$Muk.x)
    m = ncol(Thet$Muk.x)
    nk = nrow(Thet$X)
    
    Xtl <- Thet$X - matrix(c(Thet$Muk.x), ncol = ncol(Thet$X), nrow=nk, byrow = T)
    nux <- nu0.x  + nk
    #Yi0 <- matrix(c(Xtl[id,]), nrow=nk, byrow=F)
    #Psin.x <- Psi0.x + crossprod(Yi0)
    Psin.x <- Psi0.x + crossprod(Xtl)
    rinvwishart(nux, Psin.x)
    
  }
  updateSigk.x <- cmpfun(updateSigk.x)
  # 
  # A <- min(c(Wn, Mn)) - .2*diff(range(c(Wn, Mn))) 
  # B <- max(c(Wn, Mn)) + .2*diff(range(c(Wn, Mn))) 
  
  A <- min(c(Wn, Mn)) - .05*diff(range(c(Wn, Mn))) 
  B <- max(c(Wn, Mn)) + .05*diff(range(c(Wn, Mn))) 
  
  
  #--- Update X_i  
  
  udapateXi2 <- function(Thet, Y, W, M, A = A, B = B){
    
    J = ncol(Thet$X)
    n = length(Y)
    #Ytl = (Y - Thet$muk.e[Thet$Cik.e])/Thet$sig2k.e[Thet$Cik.e]
    Ytl = Y / Thet$sig2k.e[Thet$Cik.e]
    Wtl = W 
    Mtl = M #%*%diag(Thet$Betadel)
    Res = matrix(0, ncol=J, nrow = n)
    Sigw = matrix(0, J, J)
    Sigm = matrix(0, J, J)
    Sigx = matrix(0, J, J)
    GamMat= tcrossprod(c(Thet$Gamma))
    
    
    Sigw <- as.inverse(as.symmetric.matrix(Thet$Sigk.w))
    Sigm <- as.inverse(as.symmetric.matrix(Thet$Sigk.m))
    Sigx <- as.inverse(as.symmetric.matrix(Thet$Sigk.x))
    
    
    Mux <- matrix(c(Thet$Muk.x), ncol = ncol(W), nrow=length(Y), byrow=T)
    Sigxinv <- as.symmetric.matrix(Thet$Sigk.x)
    
    Sig <- as.positive.definite(as.inverse(as.symmetric.matrix(GamMat/Thet$sig2k.e + Sigw +  Sigm + Sigx))) + 0.0001*diag(J)
    
    Mu <- (matrix(c(Ytl), ncol=1) %*% matrix(Thet$Gamma,nrow=1) + tcrossprod(Wtl,Sigw) + tcrossprod(Mtl, Sigm) + 
             tcrossprod(Mux,Sigx)) %*% Sig
    
    for(i in 1:n){
      Res[i,] = tryCatch(expr={
        c(TruncatedNormal::rtmvnorm(n=1, mu = Mu[i,], sigma =  Sig,lb = rep(A,J), ub = rep(B,J)))
      },
      warning= function(W){
        if(grepl("Did not find a solution to the nonlinear system in `mvrandn`", W$message)){
          
          c(TruncatedNormal::rtmvnorm(n=1, mu = numeric(J), sigma =  Sig,lb = rep(A,J), ub = rep(B,J)))
          
        }
      }
      )  %>% as.numeric()
    }
    
    Res
    
  }
  udapateXi2 <- cmpfun(udapateXi2)
  
  
  #-- Update Xi new
  #'@ SigP0 = prior covariance matrix
  #'@ Mup0 = prior Mean 
  #'@ SigProp = proposal covariance matrix
  #'@ Muprop = proposal mean 
  
  ##----------------------------------------------------
  cat("\n #---- Create containers for the MCMC output ---> \n")
  #----------------------------------------
  #----- Divide these in sections - Construct the output
  #----- containers
  #--------------------------------------------------
  #Nsim <- 5000
  
  MCMCSample <- list()
  
  #-- for X
  Kx = 1
  MCMCSample$X <- matrix(0,nrow = Nsim,ncol=length(c(Thet$X))) ## c(t()) : how we transform the matrix (row combined) 
  MCMCSample$Muk.x <- matrix(0,nrow = Nsim, ncol = Kx*pn)
  MCMCSample$Sigk.x <- matrix(0,nrow = Nsim,ncol=Kx*pn*pn)
  #MCMCSample$Pik.x <- matrix(0,nrow = Nsim,ncol=Kx)
  #MCMCSample$Cik.x <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  
  #-- for W
  Kw = 1
  MCMCSample$Muk.w <- matrix(0,nrow = Nsim,ncol=Kw*pn)
  MCMCSample$Sigk.w <- matrix(0,nrow = Nsim,ncol=Kw*pn*pn)
  #MCMCSample$Pik.w <- matrix(0,nrow = Nsim,ncol=Kw)
  #MCMCSample$Cik.w <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  
  #-- for M
  Km=1
  MCMCSample$Muk.m <- matrix(0,nrow = Nsim,ncol=Km*pn)
  MCMCSample$Sigk.m <- matrix(0,nrow = Nsim,ncol=Km*pn*pn)
  #MCMCSample$Pik.m <- matrix(0,nrow = Nsim,ncol=Km)
  #MCMCSample$Cik.m <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  MCMCSample$Delta <- numeric(length(Thet$Delta))
  MCMCSample$Betdelt <- matrix(0, nrow=Nsim, ncol=length(Thet$Betadelt))
  MCMCSample$AcceptDelt <- numeric(Nsim)
  
  # For Y
  MCMCSample$muk.e <- matrix(0,nrow = Nsim, ncol=length(Thet0$muk.e)) ## c(t()) : how we transform the matrix (row combined)
  MCMCSample$sig2k.e <- matrix(0,nrow = Nsim, ncol=length(Thet0$muk.e))
  MCMCSample$Pik.e <- matrix(0,nrow = Nsim, ncol=Ke)
  MCMCSample$Cik.e <- matrix(0,nrow = Nsim, ncol=nrow(Wn))
  
  MCMCSample$Gamma <- matrix(0,nrow = Nsim, ncol=length(Thet0$Gamma))
  MCMCSample$BetaZ <- matrix(0,nrow = Nsim, ncol=length(Thet$BetaZ))
  MCMCSample$siggam <- numeric(Nsim)
  
  #------------------------------------------------------------------
  #Nsim = 1000 
  #library(profvis)
  #profvis(
  Sigp000 <- diag(rep(sig2.propdelt, length(Thet0$Betadel))) #solve(OptParm$hessian)
  Mu0dl <- Thet0$Betadel
  MCMCSample$Propval <- matrix(NA, nrow = Nsim, ncol = length(Thet0$Betadel))
  MCMCSample$evid <- NA*numeric(Nsim)
  Out0 <- NA
  
  cat("\n #--- MCMC is starting ----> \n")
  MCMCSample$PsiAcceptF <- numeric(6)
  for(iter in 1:Nsim){
    if(iter %% 500 == 0){cat("iter=",iter,"\n")}
    #----------------------------------------------
    # Update things related to X 
    #--------------------------------------------  
    Thet$Cik.x <-  rep(1, nrow(Wn)) 
    
    Thet$Muk.x <- updateMuk.x(Thet, Mu0.x, inSig0.x, Sig0.x)
    
    Thet$Sigk.x <- updateSigk.x(Thet, nu0.x, Psi0.x)
    
    Ytd = Y - Zmod%*%Thet$BetaZ
    
    #if( (iter %% lag == 0)*is.null(X)){
      Thet$X <- udapateXi2(Thet, Ytd, Wn, Mn, A=A, B=B)
    #}
    
    #---------------------------------------------
    #--- update things related to W
    #---------------------------------------------
    #--- Update Pik.w  
    # #updatePik.w <- function(Thet, alpha.W){
    # nk <- numeric(Kw)
    # for(k in 1:Kw){ nk[k] = sum(Thet$Cik.w == k)}
    # #update the Pis
    # alpha.wnew <- alpha.W/Kw + nk
    # #Thet$Pik.w <- c(1, rep(0, Kw-1)) #rdirch(n=1, alpha.wnew)
    # Thet$Pik.w <- rdirch(n=1, alpha.wnew)
    
    #--- Update Ci.k mixture of normals
    #Thet$Cik.w <- mapply(updateCik.w, i= 1:N, MoreArgs = list(Thet=Thet, W = Wn))
    #Thet$Cik.w <- rep(1, nrow(Wn)) #F2Cmp(CikW(Wn, Thet$X, Thet$Muk.w, Thet$Sigk.w, Thet$Pik.w)) # rep(1, nrow(Wn)) #
    # if(Gx0){Thet$Cik.w <- rep(1, length(Thet$Cik.w))}else{
    #   Thet$Cik.w <- F2Cmp(CikW(Wn, Thet$X, Thet$Muk.w, Thet$Sigk.w, Thet$Pik.w))} # rep(1, nrow(Wn)) #
    # 
    #--- Update muk.w and sig2k.w based on the mixture of normal prior
    #Thet$Muk.w <- updateMuk.w(Thet, Wn, Mu0.w, inSig0.w, Sig0.w)
    Thet$Sigk.w <- updateSigk.w(Thet, Wn, nu0.w, Psi0.w)
    
    #---------------------------------------------
    #--- update things related to M
    #---------------------------------------------
    Thet$Sigk.m <- updateSigk.m(Thet, M, nu0.m, Psi0.m, bs2=bs2)
    
    Out0 <- updateFunDeltaStNew(Thet, M, Mu0dl, mu0.del, Sig0, bs2, Sigp0 = PsiAccept * Sigp000) #Sigp000 PsiAccept
    
    #--- testing the code
    MCMCSample$Propval[iter, ] <- Out0$prop 
    Thet$Betadelt <- Out0$res[-1]
    MCMCSample$Betdelt[iter,] <- Thet$Betadelt
    MCMCSample$AcceptDelt[iter] <- Out0$res[1]
    MCMCSample$evid[iter] <- Out0$evid
    
    #---------------------------------------------------------------------------
    #--- update things related to Y or epsilon
    #--------------------------------------------------------------------------- 
    #---------------------------------------------
    #--- update Pik.e for epsilon_i
    #---------------------------------------------
    nk <- numeric(length(Thet$muk.e))
    for(k in 1:length(Thet$muk.e)){ nk[k] = sum(Thet$Cik.e == k)}
    #update the Pis
    alpha.enew <- alpha.e/length(Thet$muk.e) + nk
    Thet$Pik.e <- c(rdirch(n=1, alpha.enew))
    
    #update the Cik.e
    Thet$Cik.e <- F2Cmp(Cike(c(Y), Thet$X, Thet$muk.e, Thet$sig2k.e, Thet$Pik.e,Thet$Gamma))
    #update muk.e
    #Thet$muk.e <- updatemuk.e(Thet, Y, Z, mu0.e, sig20.e)
    
    #update Sigk.e
    Thet$sig2k.e <- updatesig2k.e(Thet,Y, Z, gam0.e, sig20.esig)
    #---------------------------------------------
    #--- update Gamma and BetaZ
    #---------------------------------------------
    Tp0 <- updateGam(Thet, Y, Z, Sig0.gam0, Mu0.gam0)
    Thet$BetaZ <- Tp0[c(1:ncol(Zmod))] #c(0, Tp0[c(2:ncol(Zmod))])  #Tp0[c(1:ncol(Zmod))]
    Thet$Gamma <- Tp0[-c(1:ncol(Zmod))]

    if(is.null(sig02.del)){
    Thet$siggam <- updateSigGam(Thet$Gamma, al0, bet0, order = 2)  # updateSigGam <- function(Bet,al0, bet0, order = 2) ##1/tauval #
    }
    #---------------------------------------------------------------------------
    #------------------------- Save output
    #---------------------------------------------------------------------------
    #-- For X
    MCMCSample$Muk.x[iter, ] <- c(t(Thet$Muk.x))  ## Stor rowwise
    MCMCSample$Sigk.x[iter, ] <- c(t(Thet$Sigk.x))
    #MCMCSample$Pik.x[iter, ] <- Thet$Pik.x
    #MCMCSample$Cik.x[iter, ] <- Thet$Cik.x
    MCMCSample$X[iter, ] <- c(t(Thet$X)) ## Store row-wise
    
    #-- for W
    MCMCSample$Muk.w[iter, ] <- c(t(Thet$Muk.w))
    MCMCSample$Sigk.w[iter, ] <- c(t(Thet$Sigk.w))
    #MCMCSample$Pik.w[iter, ] <- Thet$Pik.w
    #MCMCSample$Cik.w[iter, ] <- Thet$Cik.w
    
    #-- for M
    MCMCSample$Muk.m[iter, ] <- c(t(Thet$Muk.m))
    MCMCSample$Sigk.m[iter, ] <- c(t(Thet$Sigk.m))
    #MCMCSample$Pik.m[iter, ] <- Thet$Pik.m
    #MCMCSample$Cik.m[iter, ] <- Thet$Cik.m
    #MCMCSample$AcceptDelt[iter] <- Out0[1]
    #MCMCSample$Delta[iter] <- Thet$Delta
    #-- For Y
    MCMCSample$muk.e[iter, ] <- Thet$muk.e ## c(t()) : how we transform the matrix (row combined)
    MCMCSample$sig2k.e[iter, ] <- c(Thet$sig2k.e)
    MCMCSample$Pik.e[iter, ] <- Thet$Pik.e
    MCMCSample$Cik.e[iter, ] <- Thet$Cik.e
    
    MCMCSample$Gamma[iter, ] <- Thet$Gamma
    MCMCSample$BetaZ[iter, ] <- Thet$BetaZ
    MCMCSample$siggam[iter] <- Thet$siggam
  }  
  
  #---- Return the MCMC smaples
  cat("\n ---- MCMCs are done!! ---/n")
  #-- Same the raw (pointwise) estimate of delta(t)
  MCMCSample$Delthat <- colMeans(M)/colMeans(W) #Delthat
  
  #-- Same the smooth estimate of delta(t)
  MCMCSample$EstimDelthat <- (bs2del%*%Thet0$Betadel)
  
  #--- Save the initial values
  MCMCSample$Thet0 <- Thet0
  #MCMCSample$OptimOut <- OptParm
  MCMCSample
  
  
}
#===============================================================================

#===============================================================================
#==== Exactly on component of the normal distribution
BayesFIVSingV2 = function(Nsim = 100, Y, W, M, Z, sig02.del = NULL, Ge = 2,  
                        alpx = 1, BasFun = "Bspline", PsiAccept = 0.05, 
                        sig2.propdelt = 2.75, sig2.priordelt = 5.0){
  #--- Functions
  func <- function(Pik, n=1, size=1){
    which.max(rmultinom(n=1,size=1,prob = Pik + abs(min(Pik))+ 0.001))
  }
  funcCmp <- cmpfun(func)
  
  F2 <- function(Pik, n=1, size=1){
    
    #print(Pik)
    apply(Pik,1,funcCmp)
  }
  F2Cmp <- cmpfun(F2)
  
  a <-  seq(0, 1, length.out = ncol(M))
  n = length(Y)
  cat("--- set some initial values ")
  #---- Parameters settings
  #-----------------------------------------------------------------------------
  pn <- round(n^{1/3}) + 2 ## min(15, ceiling(length(Y)^{1/3.8})+4)
  
  alpha.e = alpx
  alpha.w = alpha.W = alpx
  alpha.x = alpha.X = alpx
  alpha.m = alpx
  pn = pn 
  n = length(Y) 
  Zmod = Z
  
  #-----------------------------------------------------------------------------
  #--- Transform the data
  #-----------------------------------------------------------------------------
  bs2 <- bs(a, df = pn, intercept = T)
  bs2del <- qr.Q(qr(bs2))
  
  #--- Using different bases
  
  if(BasFun == "Four"){
    
    #Nbasis0 = c(3,4, 5, 7) + 2
    J = K = k  = pn = pn00 =  Nbasis0[which(c(100,200,500, 1000) == n)]
    #Nbasis0 =  c(3, 6, 15, 18)
    J = K = k  = pn = pn00 =  Nbasis0[which(c(100,200,500, 1000) == n)]
    
    
    Fourmeth <- fda::create.fourier.basis(c(0,1), nbasis = pn)
    bs2 = as.matrix(eval.basis(a, Fourmeth))
    colnames(bs2) <- NULL
    bs2del <- qr.Q(qr(bs2)) 
    
  }
  
  Mn <- (M%*%bs2)/length(a) #(M%*%bs2)/length(a) 
  Wn <- (W%*%bs2)/length(a)
  
  #--- Penalty matrix
  
  P.mat <- function(K){
    # penalty matrix
    D <- diag(rep(1,K))
    D <- diff(diff(D))
    P <- t(D)%*%D 
    return(P)
  }
  
  #--- Estimate initial values
  InitialValues <- function(Y, W, M, Z, df0=3, a, Ge, bs2){
    Gx = 3
    Res <- list()
    #bs.w <- splines::bs(a, df = pn+1, degree = min(c(df0,pn-1)), intercept=TRUE) 
    Zmod <- Z #cbind(1,Z)
    pn = ncol(bs2)
    
    
    bs2del <- qr.Q(qr(bs2))
    
    Wn <- crossprod(t(W), bs2)/length(a)
    Mn <- crossprod(t(M), bs2)/length(a)
    
    WsmB <- Wn
    WsM <- NULL
    for(i in 1:nrow(W)){ WsM <- rbind(WsM, smooth.spline(a, W[i,])$y)}
    W_i <- crossprod(t(WsM),bs2)/length(a)   #(n,k) matrix
    lmFit <- lm(Y ~ -1 + Zmod + W_i) 
    c_hatW <- lmFit$coefficients
    Res$Gamma <- as.numeric(c_hatW[-c(1:ncol(Zmod))])
    Res$BetaZ <- as.numeric(c_hatW[c(1:ncol(Zmod))])
    Res$Delta <- abs(mean(Mn)/mean(W_i))
    #Res$Betadel <- c(abs(solve(t(bs2del)%*%bs2del)%*%t(bs2del)%*%(colMeans(M)/colMeans(W))))
    
    #--- Estimate Delta
    Delhat <- function(bs, U0){
      
      Fe0 <- function(val){
        sum((bs%*%val - U0)^2)
      }
      
      k = ncol(bs)
      Res <- optim( rep(0, ncol(bs)), Fe0, method="L-BFGS-B", lower=rep(-60,k), upper=rep(60, k) )
      Res  
      
    }
    U0 = colMeans(M)/colMeans(W)
    R0d <- Delhat(bs2del, U0)
    
    Res$Betadel <- R0d$par 
    
    #---------------------------------------------------------------------------
    #--- Determine the number of clusters for epsilon_i
    #---------------------------------------------------------------------------
    
    reslmfit <- residuals.lm(lmFit)
    BIC <- mclustBIC(reslmfit)
    mod1 <- Mclust(reslmfit, x = BIC, G = Ge)
    Temp <- summary(mod1, parameters = TRUE)
    Res$Ke <- Temp$G #+ Gx
    Res$muk.e <- numeric(Res$Ke)
    Res$muk.e[1:Temp$G] = Temp$mean; 
    Res$sig2k.e <- rep(1, Res$Ke)*sum(lmFit$residuals^2)/lmFit$df.residual
    Res$sig2k.e[1:Temp$G] = Temp$variance
    Res$Cik.e <- as.numeric(Temp$classification)
    Res$Pik.e <- numeric(Res$Ke)
    Res$Pik.e[1:Temp$G] <- Temp$pro 
    #---------------------------------------------------------------------------
    #--- Determine the number of cluster (W)
    #---------------------------------------------------------------------------
    Ui <- Wn #- W_i 
    Mod1 <- Mclust(Ui, G = Gx)
    Tp <- summary(Mod1, parameters=T)
    Res$Kw <- Tp$G #+Gx
    Res$Muk.w <- matrix(0, ncol = pn, nrow=Res$Kw)
    Res$Muk.w <- as.matrix( rbind(t(Tp$mean))) #  ,  repmat(rep(0,pn), Gx, 1)))
    
    Res$Sigk.w <- matrix(0, nrow = Res$Kw, ncol = pn*pn)
    for(j in 1:Tp$G){Res$Sigk.w[j,] <- c(Tp$variance[,,j])}
    
    Res$Cik.w <- Tp$classification
    Res$Pik.w <- Tp$pro
    
    #---------------------------------------------------------------------------  
    #--- Determine the number of cluster (X)
    #---------------------------------------------------------------------------
    WsmB <- (W + M) /(1+Res$Delta)
    WsM <- NULL
    for(i in 1:nrow(W)){ WsM <- rbind(WsM, smooth.spline(a, WsmB[i,])$y)}
    W_i <- crossprod(t(WsM),bs2)/length(a)   #(n,k) matrix
    Ui <- W_i 
    Mod1 <- Mclust(Ui, G=3)
    Tp <- summary(Mod1, parameters=T)
    Res$Kx <- Tp$G #+ Gx
    Res$Muk.x <- matrix(0, ncol = pn, nrow=Res$Kw)
    Res$Muk.x <- as.matrix( rbind(t(Tp$mean),  repmat(colMeans(Ui), Gx, 1 )))
    
    Res$Sigk.x <- matrix(0, nrow = Res$Kx, ncol = pn*pn)
    for(j in 1:Tp$G){Res$Sigk.x[j,] <- c(Tp$variance[,,j])}
    #for(j in (Tp$G+1):Res$Kx){Res$Sigk.x[j,] <- c(.5*diag(pn))}
    Res$Cik.x <- Tp$classification
    Res$Pik.x <- Tp$pro  
    Res$X <- W_i
    #---------------------------------------------------------------------------
    #---- Determine the number of cluster (M)
    #---------------------------------------------------------------------------
    
    Ui <- Mn  #- Res$Delta*W_i 
    Mod1 <- Mclust(Ui, G=Gx)
    Tp <- summary(Mod1, parameters=T)
    Res$Km <- Tp$G #+ Gx
    #Res$Muk.m <- as.matrix(t(Tp$mean))
    Res$Muk.m <- as.matrix( rbind(t(Tp$mean))) #,  repmat(colMeans(Ui), Gx, 1 )))
    Res$Sigk.m <- matrix(0,nrow=Res$Km,ncol=pn*pn)
    for(j in 1:Tp$G){Res$Sigk.m[j,] <- c(Tp$variance[,,j])}
    #for(j in (Tp$G+1):Res$Km){Res$Sigk.m[j,] <- c(as.positive.definite(.5*diag(pn)))}
    Res$Cik.m <- Tp$classification
    Res$Pik.m <- Tp$pro
    Res$R0d <- R0d
    #------- return result 
    Res 
  }
  
  
  cat("\n #--- Obtain initial values...\n")
  
  Thet0 <- InitialValues(Y, W, M, Z, df0=3, a, Ge, bs2)
  
  #-----------------------------------------------------------------------------
  cat("\n#---- Initialize parameters ---> \n")
  #-----------------------------------------------------------------------------
  
  #Kx = Thet0$Kx               ## number of X  cluster
  Ke = Thet0$Ke               ## number of ei cluster
  #Km = Thet0$Km               ## number of ei cluster
  #Kw = Thet0$Kw               ## number of omega_i cluster 
  
  
  Mn0 <- (M%*%bs2)/length(a) #(M%*%bs2)/length(a) 
  Wn <- (W%*%bs2)/length(a)
  Mn <- ((M%*%diag(1/(c(bs2del%*%Thet0$Betadel))))%*%bs2)/length(a) 
  
  Thet <- list()
  
  Thet$Gamma <-  rnorm(n=length(Thet0$Gamma)) #as.numeric(Thet0$Gamma)   #- a vector of length J
  Thet$BetaZ <- Thet0$BetaZ
  Thet$siggam <-  ifelse(is.null(sig02.del), 1, sig02.del)
  Thet$bs2del <- bs2del
  Thet$Betadelt <- Thet0$Betadel
  
  
  Thet$Cik.x = rep(1, n) #Thet0$Cik.x #sample(x = 1:Kx,size = n, replace = T)    #- vector of length Kx
  Thet$Cik.e = rep(1, n) #Thet0$Cik.e #sample(x = 1:Ke,size = n, replace = T)     #- vector of length Ku
  Thet$Cik.w = rep(1, n)  #Thet0$Cik.w  #sample(x = 1:Kw,size = n, replace = T)     #- vector of length Kw
  Thet$Cik.m = rep(1, n)  #Thet0$Cik.m #sample(x = 1:Km,size = n, replace = T)     #- vector of length Km
  
  
  Thet$Pik.x = Thet0$Pik.x #rdirch(n=1,alpha = rep(alpha.X, Kx)/Kx)   #- vector of length Kx
  Thet$Pik.e = Thet0$Pik.e #rdirch(n=1,alpha = rep(alpha.e, Ke)/Ke)  #- vector of length Ku
  Thet$Pik.w = Thet0$Pik.w #rdirch(n=1,alpha = rep(alpha.W, Kw)/Kw)  #- vector of length Kw
  Thet$Pik.m = Thet0$Pik.m #rdirch(n=1,alpha = rep(alpha.m, Km)/Km)  #- vector of length Kw
  
  
  Thet$sig2k.e  =  Thet0$sig2k.e  #rgamma(n=Ke, shape=1,scale=.1) # matrix dim Kx x J
  Thet$muk.e  =  Thet0$muk.e #rnorm(n=Ke)  #- Matrix of (Ku x J) and Ku x J - Note that $\Omega$ is a diagonal matrix
  
  Thet$Sigk.w <- .5*diag(pn)     #- vector of length Kw
  Thet$Muk.w <- rep(0,pn) #c(rnorm(n = pn))         #- vectors of length Kw and Kw respectively
  
  Thet$Sigk.x <- .5*diag(pn)     #-- vector of length Kx
  Thet$Muk.x  <- rep(0,pn) #rnorm(n=pn)   #- vectors of length Kx and Kx respectively
  
  Thet$Sigk.m <- .5*diag(pn)     #-- vector of length Kx
  Thet$Muk.m  <- rep(0,pn) #c(rnorm(n = pn))   #- vectors of length Kx and Kx respectively
  
  
  A <- min(c(c(Wn), c(Mn))) - 0.05*diff(range(c(Wn, Mn)))
  B <- max(c(c(Wn), c(Mn))) + 0.05*diff(range(c(Wn, Mn)))
  
  Thet$X <- TruncatedNormal::rtmvnorm(n=n,mu=rep(0,pn), sigma=diag(pn),lb = rep(A, pn), ub = rep(B, pn))
  
  
  cat("\n#---- Loading updating functions ---> \n")
  #-----------------------------------------------------------------------------
  # Update things related to Y
  #-----------------------------------------------------------------------------
  #-------- Update Gamma (vector of length pn + the Z together) 
  #--- Prior parms
  pn0 = ncol(Z)
  Pgam <- P.mat(pn)*1.0
  Sig0.gam0 = diag(pn0); Mu0.gam0 = numeric(pn0+pn);
  
  updateGam <- function(Thet, Y, Z, Sig0.gam0, Mu0.gam0){
    
    pn = ncol(Thet$X)
    Zmod <-  Z #cbind(1,Z)
    Y.td = Y - Thet$muk.e[Thet$Cik.e]
    InvSig0.gam0 <- 1.0*as.matrix(Matrix::bdiag(as.inverse(Sig0.gam0), P.mat(pn)/Thet$siggam))
    Xmod <- cbind(Zmod, Thet$X)
    X_sc <- sweep(Xmod, MARGIN = 1, Thet$sig2k.e[Thet$Cik.e], FUN="/") ## n*pn
    
    Sig.gam <- as.inverse(as.symmetric.matrix(crossprod(X_sc, Xmod) + InvSig0.gam0))
    Mu.gam <- c(Sig.gam%*%(crossprod(X_sc, Y.td) + InvSig0.gam0%*%Mu0.gam0))
    
    c(mvtnorm::rmvnorm(n=1, mean = Mu.gam, sigma = as.positive.definite(Sig.gam)))
    
  }
  updateGam <- cmpfun(updateGam)
  
  #Update siggma - siggam ~ IG(al0, bet0)
  al0 = 1
  bet0 = .005
  updateSigGam <- function(Bet, al0, bet0, order = 2){
    K = length(Bet)
    al = al0 + .5*(K - order)
    bet = .5*quad.form(P.mat(K), Bet) + bet0
    #1/rgamma(n = 1, shape = al, rate = bet)
    #nimble::rinvgamma(n = 1, shape = al, scale = bet)
    # tb = 10
    # 1/heavy::rtgamma(n=1,shape=al,scale = bet, t = tb)
    tb = 1
    1/cascsim::rtgamma(n=1,shape=al, scale = 1/bet, max = 10, min = 1/15)
    #tb = .01
    #ifelse( qgamma(.95, shape= al, rate = bet) > tb, 1/tb , 1/cascsim::rtgamma(n=1,shape=al, scale = 1/bet, max = tb))
  }
  
  #------- Update Pik.e 
  #alpha.e = 1.0
  updatePik.e <- function(Thet, alpha.e){
    
    nk <- numeric(length(Thet$muk.e))
    for(i in 1:length(Thet$muk.e)){ nk[i] = sum(Thet$Cik.e == i)}
    
    #update the Pis
    alpha.enew <- alpha.e/length(Thet$muk.e) + nk
    
    rdirichlet(n=1, alpha.enew)
    
  }
  
  #--- Update muk.w and sig2k.w based on the mixture of normal prior Xtl
  mu0.e = 0 ; sig20.e = 100
  updatemuk.e <- function(Thet, Y, Z, mu0.e, sig20.e ){

    Ku = length(Thet$muk.e)
    muk <- numeric(Ku)
    sig2k <- numeric(Ku)
    Zmod <- Z #cbind(1,Z)
    #Ytl <- Y - Thet$X%*%Thet$Gamma
    Ytl = Y - Zmod%*%Thet$BetaZ -  Thet$X%*%Thet$Gamma
    for(k in 1:length(Thet$muk.e)){
      id <- which(Thet$Cik.e == k)
      nk = length(id)
      if(nk > 0){
        sig2k[k] <- 1/(nk/Thet$sig2k.e[k] + 1/sig20.e)
        muk[k] <- sig2k[k]*(mu0.e/sig20.e + sum(c(Ytl[id]))/Thet$sig2k.e[k])
      }else{
        sig2k[k] <- sig20.e
        muk[k] <- mu0.e
      }
    }
    ## Force the zero mean constrains
    SigR0 <- Thet$Pik.e*sig2k
    SigRR <- sum((Thet$Pik.e^2)*sig2k)
    MuK0 = muk - SigR0*(1/SigRR)*sum(Thet$Pik.e*muk)
    SigK <- as.matrix(nearPD(diag(sig2k) - tcrossprod(SigR0)*(1/SigRR), doSym = T)$mat)

    if(Ku > 2){
      id0 = 1:(Ku-1)
      Muk_1 <- c(mvtnorm::rmvnorm(n=1,mean=MuK0[-Ku], sigma = SigK[id0,id0]))

      c(Muk_1, -sum(c(Thet$Pik.e[-Ku]*Muk_1)/ifelse(Thet$Pik.e[Ku] == 0.0, (.Machine$double.eps)^{.7}, Thet$Pik.e[Ku])))
    }else{

      id0 = 1:(Ku-1)
      Muk_1 <- c(rnorm(n=1,mean=MuK0[-Ku], sd = sqrt(SigK[id0,id0])))

      c(Muk_1, -(c(Thet$Pik.e[-Ku]*Muk_1)/ifelse(Thet$Pik.e[Ku] == 0.0, (.Machine$double.eps)^{.7}, Thet$Pik.e[Ku])))
    }

  }

  gam0.e = 1 ; sig20.esig = 1
  #gam0.e = 1 ; sig20.e = .5
  updatesig2k.e <- function(Thet, Y, Z, gam0.e, sig20.e){
    
    Zmod <- Z #cbind(1,Z)
    Ku = length(Thet$muk.e)
    al <- numeric(Ku)
    bet <- numeric(Ku)
    Xtl <- Y - Thet$X%*%Thet$Gamma - Zmod%*%Thet$BetaZ  #- Thet$muk.e[Thet$Cik.e]
    for(k in 1:length(Thet$muk.e)){
      id <- which(Thet$Cik.e == k)
      nk = length(id)
      if(nk > 0){
        al[k] <- nk/2 + gam0.e
        bet[k] <- sig20.e + .5*sum(Xtl[id]^2) 
      }else{
        al[k] <- gam0.e
        bet[k] <- sig20.e 
      }
    }
    ## Force the zero mean constrains
    mapply(rinvgamma, n = 1, shape = al, scale = bet)
  }
  
  #-----------------------------------------------------------------------------
  # Update things related to W 
  #-----------------------------------------------------------------------------
  
  #--- Update Pik.x  
  #alpha.W = alpha.w = 1.0
  updatePik.w <- function(Thet, alpha.W){
    Ku = nrow(Thet$Muk.w)
    nk <- numeric(Ku)
    for(k in 1:Ku){ nk[k] = sum(Thet$Cik.w == k)}
    
    #update the Pis
    alpha.wnew <- alpha.W/length(Thet$muk.w) + nk
    
    rdirichlet(n=1, alpha.wnew)
    
  }
  
  #--- Update Ci.k mixture of normals 
  dmvnormVec <- function(x, Mean, SigmaVc){
    Ku <- nrow(Mean)
    Val <- numeric(Ku)
    for(k in 1:Ku){ Val[k] <- mvtnorm::dmvnorm(x = x, mean = c(Mean[k,]), sigma = as.symmetric.matrix(matrix(SigmaVc[k,],ncol=length(x),byrow=F))) }
    Val
  }
  
  updateCik.w <- function(i, Thet, W){ # mixture of normals 
    #Sig <- matrix(Thet$Sigk.w[id,], ncol = ncol(W), byrow = F)
    #kw = nrow(Thet$Muk.w)
    xi =  c(W[i,] - Thet$X[i,])
    #pi.k  = Thet$Pik.w*mapply(dmvnormVec, k = 1:Kw, MoreArgs = list(x = xi, Mean = Thet$Muk.w , SigmaVc = Thet$Sigk.w))
    pi.k  = Thet$Pik.w*dmvnormVec(x=xi, Mean = Thet$Muk.w , SigmaVc = Thet$Sigk.w) + .Machine$double.eps
    which(rmultinom(n=1,size=1,prob = pi.k) == 1)
  }
  
  #--- Update Cik.w 
  #--- Update muk.w and sig2k.w based on the mixture of normal prior
  Mu0.w = rep(0, pn);  Sig0.w = .5*diag(pn);inSig0.w = 2*diag(pn) # as.inverse(.5*cov(Wn)) #.5*diag(pn) #
  updateMuk.w <- function(Thet, W, Mu0.w, inSig0.w, Sig0.w){
    Ku = length(Thet$Pik.w)
    J = ncol(Thet$X)  
    Mu.k <- matrix(0, nrow = Ku, ncol = J) # K x J
    SigK <- list()  #-- store the covariances
    SigKWg <- matrix(0, ncol=J, nrow=J) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    Yi.tld = W - Thet$X 
    #Muk <- matrix(0, ncol=J,nrow=Ku)
    SigK <- list()
    for(k in 1:Ku){
      idk <- which(Thet$Cik.w == k)
      nk <- length(idk)
      if(length(idk) > 0){
        Yi0 = matrix(c(Yi.tld[idk,]), nrow= nk, byrow=F)
        TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.w[k,],ncol=J,byrow=F)))
        SigK[[k]] <- as.inverse(as.symmetric.matrix(nk*TpSi + inSig0.w)) ## sum of matrices inverse
        SigKWg <- SigKWg + (Thet$Pik.w[k]*Thet$Pik.w[k])*SigK[[k]] # J x J
        #Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi.tld[idk,])) ) + solve(Sig0.w)%*%Mu0.w)      ## some of mat
        Mu.k[k,] <- c(SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi0)) ) + inSig0.w%*%Mu0.w))      ## some of mat
        #Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(apply(t(Yi.tld[idk,]), 2, sum)) ) + solve(Sig0.w)%*%Mu0.w)      ## some of mat
      }else{
        SigK[[k]] <- Sig0.w
        SigKWg <- SigKWg + (Thet$Pik.w[k]*Thet$Pik.w[k])*SigK[[k]] # J x J
        Mu.k[k,] <- Mu0.w
      }
    }
    #--- we impose some constrainsts on the mean
    SIg0 <- bdiag(SigK) # create a block diagonal matrix
    SIGR0 <- NULL; for(k in 1:Ku){SIGR0 <- rbind(SIGR0,Thet$Pik.w[k]*SigK[[k]])} ## K*J x J
    MUR0 <-  colSums(diag(Thet$Pik.w) %*% Mu.k) # return a evector of length J
    
    InvSigKWg <- as.positive.definite(as.inverse(SigKWg))
    # print(dim(InvSigKWg))
    # print(dim(SIg0))
    # print(dim(quad.tform.inv(SigKWg,SIGR0)))
    # print(Ku)
    #MURFin <- c(t(Mu.k)) - SIGR0%*%solve(SigKWg)%*%MUR0 ## J*K
    MURFin <- c(t(Mu.k)) - SIGR0%*%InvSigKWg%*%MUR0 ## J*K
    #SIGR0fin <- 1.0*as.matrix(nearPD(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T)$mat) ## J*K x J*K
    #SIGR0fin <- 1.0*(as.positive.definite(as.symmetric.matrix(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T))) ## J*K x J*K
    SIGR0fin <- 1.0*(as.positive.definite(as.symmetric.matrix(as.matrix(SIg0 - SIGR0%*%InvSigKWg%*%t(SIGR0)))))
    #--- Simulate K-1 vectors from the degenerate dist.
    id <- 1:(J*(Ku-1))
    SimMU.k <- matrix(c(mvtnorm::rmvnorm(n=1, mean = MURFin[id], sigma = SIGR0fin[id,id])), ncol=J,byrow=T) ## (K-1) x J
    SimMU.k1 <- -colSums(diag(Thet$Pik.w[-Ku]) %*% SimMU.k)/Thet$Pik.w[Ku] ## get a K-1 x J matrix follow by colSums - a vector of J
    as.matrix(rbind(SimMU.k,SimMU.k1)) ## k * J
  }
  updateMuk.w <- cmpfun(updateMuk.w)
  
  #--- update Sigk.w
  nu0.w = pn + 2; Psi0.w = .5*diag(pn) 
  updateSigk.w <- function(Thet, W, nu0.w, Psi0.w){
    
    Ku = nrow(Thet$Muk.w)
    m = ncol(W)
    nk = nrow(W)
    
    Wtl <- W - Thet$X 
    nun <- nu0.w  + nk
    #Yi0 <- matrix(c(Wtl[id,]), nrow= nk, byrow=F)
    #Psin.w <- as.positive.definite(Psi0.w + crossprod(Yi0))
    Psin.w <- as.positive.definite(Psi0.w + crossprod(Wtl))
    rinvwishart(nun, Psin.w)
    
  }
  updateSigk.w <- cmpfun(updateSigk.w)
  
  #-----------------------------------------------------------------------------
  # Update things related to M 
  #-----------------------------------------------------------------------------
  
  #--- update Pik.m for epsilon_i
  #alpha.m = 1.0
  updatePik.m <- function(Thet, alpha.m){
    nk <- numeric(nrow(Thet$Muk.m))
    for(k in 1:nrow(Thet$Muk.m)){ nk[k] = sum(Thet$Cik.m == k)}
    
    #update the Pis
    alpha.mnew <- alpha.m/nrow(Thet$Muk.m) + nk
    
    c(rdirichlet(n=1, alpha.mnew))
  }
  
  #--- update Cik.m for epsilon_i based on mixture of normals
  updateCik.m <- function(i, Thet, M){ # mixture of normals 
    #id <- Thet$Cik.m[i]
    #Sig <- matrix(Thet$Sigk.m[id,], ncol = ncol(M), byrow = F)
    Km = nrow(Thet$Muk.m)
    xi = M[i,] - Thet$Delta*Thet$X[i,]
    #pi.k  = Thet$Pik.m*mapply(dmvnormVec, k = 1:Km, MoreArgs = list(x = xi , Mean = Thet$Muk.m , SigmaVc = Thet$Sigk.m))
    pi.k  = Thet$Pik.m*dmvnormVec(x = xi, Mean = Thet$Muk.m , SigmaVc = Thet$Sigk.m) + .Machine$double.eps
    which(rmultinom(n=1,size=1,prob = pi.k) == 1)
  }
  
  
  #--- Update Muk.w and sig2k.w based on the mixture of normal prior
  Mu0.m = rep(0,pn) ; Sig0.m = .5*diag(pn); inSig0.m = 2*diag(pn) 
  updateMuk.m <- function(Thet, M, Mu0.m, inSig0.m, Sig0.m){
    #cat(dim(inSig0.m),"(--1)\n")
    Ku = length(Thet$Pik.m)
    J = ncol(M)  
    Mu.k <- matrix(0, nrow = Ku, ncol = J) # K x J
    SigK <- list()  #-- store the covariances
    SigKWg <- matrix(0, ncol = J, nrow = J) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    Yi.tld = M - Thet$X%*%diag(Thet$Betadel) 
    Muk <- matrix(0, ncol=J, nrow=Ku)
    for(k in 1:Ku){
      idk <- which(Thet$Cik.m == k)
      nk <- length(idk)
      if(length(idk) > 0){
        Yi0 = matrix(c(Yi.tld[idk,]), nrow = nk, byrow=F)
        TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.m[k,],ncol=J,byrow=F)))
        SigK[[k]] <- as.inverse(as.symmetric.matrix(nk*TpSi + inSig0.m)) ## sum of matrices inverse
        #cat(dim(SigK[[k]]),"(--2 \n")
        SigKWg <- SigKWg + (Thet$Pik.m[k]*Thet$Pik.m[k])*SigK[[k]] # J x J
        Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi0)) ) + inSig0.m%*%Mu0.m)      ## some of mat
        #cat(dim(SigKWg),"(--3 \n")
      }else{
        SigK[[k]] <- Sig0.m
        SigKWg <- SigKWg + (Thet$Pik.m[k]*Thet$Pik.m[k])*SigK[[k]] # J x J
        Mu.k[k,] <- Mu0.m
      }
    }
    #--- we impose some constrainsts on the mean
    InvSigKWg <- as.positive.definite(as.inverse(SigKWg))
    SIg0 <- bdiag(SigK) # create a block diagonal matrix
    SIGR0 <- NULL;  for(k in 1:Ku){SIGR0 <- rbind(SIGR0,Thet$Pik.m[k]*SigK[[k]])} ## K*J x J
    MUR0 <-  colSums(diag(Thet$Pik.m) %*% Mu.k) # return a evector of length J
    
    # print(Ku)
    # print(dim(SIg0))
    # print(dim(SIGR0))
    #MURFin <- c(t(Mu.k)) - SIGR0%*%solve(SigKWg)%*%MUR0 ## J*K
    #SIGR0fin <- as.matrix(nearPD(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T)$mat) ## J*K x J*K
    
    MURFin <- c(t(Mu.k)) - SIGR0%*%InvSigKWg%*%MUR0 ## J*K
    SIGR0fin <- as.symmetric.matrix(as.positive.definite(as.matrix(SIg0 - SIGR0%*%InvSigKWg%*%t(SIGR0)))) ## J*K x J*K
    
    #--- Simulate K-1 vectors from the degenerate dist.
    id <- 1:(J*(Ku-1))
    SimMU.k <- matrix(c(mvtnorm::rmvnorm(n=1, mean = MURFin[id], SIGR0fin[id,id])), ncol=J,byrow=T) ## (K-1) x J
    SimMU.k1 <- -colSums(diag(Thet$Pik.m[-Ku]) %*% SimMU.k)/Thet$Pik.m[Ku] ## get a K-1 x J matrix follow by colSums - a vector of J
    as.matrix(rbind(SimMU.k,SimMU.k1)) ## k * J
  }
  updateMuk.m <- cmpfun(updateMuk.m)
  
  #--- Update Sig2k.w based on the mixture prior
  nu0.m = pn + 2; Psi0.m = .5*diag(pn)#as.symmetric.matrix(.25*cov(Wn))
  updateSigk.m <- function(Thet, M, nu0.m, Psi0.m, bs2=bs2){
    
    Ku = nrow(Thet$Muk.m)
    m = ncol(Thet$X)
    nk = nrow(M) 
    
    Mtl <- (M - (tcrossprod(Thet$X, bs2) * matrix(c(Thet$bs2del%*%Thet$Betadel), nrow=nk, ncol = ncol(M), byrow=T)))%*%bs2 
    
    nun <- nu0.m  + nk
    Psin.m <- as.positive.definite(Psi0.m + crossprod(Mtl))
    rinvwishart(nun, Psin.m)
  }
  updateSigk.m <- cmpfun(updateSigk.m)
  
  #---- Update things related to delta(t)
  #---- Function to Compute  hij matrix (N x T)
  LogSimDElVec <- function(M, Mn, X, Sigma, Val=1){
    #LogSimDElVec <- function(M, Mn, X, Sigma, Val){
    id <- Val
    #mvtnorm::dmvnorm(M[i,], mean = Mn[i,], sigma = Sigma, log = T)
    #LaplacesDemon::dmvnc(M[i,], mu = Mn[i,], U = Sigma, log = T)
    LaplacesDemon::dmvnc(M, mu = Mn, U = Sigma, log = T)
  }
  
  #-- Makes the hij (N x T) matrix
  Hij <- function(Betadel, Thet, bs2){
    
    (Thet$X%*% t(bs2)) * matrix(c(Thet$bs2del%*%matrix(Betadel,ncol=1)), nrow = nrow(Thet$X), ncol = nrow(bs2), byrow = T) ## Matrix NxK
  }
  Hij <- cmpfun(Hij)
  
  
  
  #'@ Muvec0 : prior mean vector
  #'@ Sig0 : prior variance
  #'@ likelihood + prior
  #--- Second version (without projecting delta(t)*X_{i}(t))
  LikHijMat <- function(Betadel, Thet, M, Muvec0, Sig0, bs2, Val=1){
    Sig = as.positive.definite(as.symmetric.matrix(quad.tform(Thet$Sigk.m, bs2)))
    Sigchol <- chol(Sig)
    MeanV = Hij(Betadel, Thet, bs2)
    #sum(LogSimDElVec(1:nrow(M), M, MeanV, Thet$X, Thet$Sigk.m, Val)) + mvtnorm::dmvnorm(Betadel, mean = Muvec0, sigma = Sig0, log = T)
    sum(LogSimDElVec(M, MeanV, Thet$X, Sigchol, Val=1)) + mvtnorm::dmvnorm(Betadel, mean = Muvec0, sigma = Sig0, log = T)
  }
  
  # 
  #--- Without projecting the function to lower space.
  LikHijMatOptim <- function(Thet, M, Muvec, Sig0, bs2, Val=1){
    #Thet$Betadel
    #MeanV = ((Thet$X%*% t(bs2) * matrix(c(bs2%*%matrix(Betadel,ncol=1)), nrow = nrow(Thet$X), ncol = nrow(bs2), byrow = T))%*%bs2)/nrow(bs2) ## Matrix NxK
    J = length(Muvec)
    
    #lb = Muvec - 0.5*abs(Muvec)
    #ub = Muvec + 0.5*abs(Muvec)
    
    parm0 <- c(rnorm(n = J))
    SigU = chol(as.positive.definite(quad.tform(Thet$Sigk.m, bs2)))
    #print(dim(SigU))
    
    funC <- function(parm){
      MeanV = Hij(parm, Thet, bs2)
      -sum(LogSimDElVec(M, MeanV, Thet$X, SigU, Val=1)) - mvtnorm::dmvnorm(parm, mean = Muvec,sigma = Sig0, log = T)
    }
    
    Res <- optim(parm0, funC, method="L-BFGS-B", lower=rep(-15,J), upper = rep(15,J), hessian = T)
    #Res <- optim(parm0, funC, method="L-BFGS-B", lower=c(rep(lb,J), 0.01), upper = c(rep(ub,J), 10.0), hessian = T)
    Res
  } 
  
  #--- Updated on 10/16/23 (update Deltat)
  #@' M : matrix of n (samples) x T(number of distinct time points)
  #@' Sigp0 = proposal covariance
  #@' Muvec : prior mean vector
  #@' Sig0 : prior variance
  #@' bs2 : are the basis function
  #mu0.del, Sig0.del
  #Muvec = c(t(bs2)%*%colMeans(M)/colMeans(W))
  #Muvec = mu0.del = rep(0, ncol(Mn)); Sig0 = Sig0.del = diag(rep(.3,ncol(Mn))); Sigp0 = diag(rep(.2,ncol(Mn)))
  #'@ Find the covariance matrix for the  proposal distribution
  
  Mu0dl = rep(0, ncol(Mn))
  Muvec = mu0.del = rep(0, ncol(Mn)); Sig0 = Sig0.del = diag(rep(sig2.priordelt,ncol(Mn))); Sigp0 = diag(rep(.2,ncol(Mn)))
  Muvec = mu0.del = Thet0$Betadel #c(t(bs2) %*% matrix(c(colMeans(M)/colMeans(W)), ncol=1))
  
  #--- expanded version 
  updateFunDeltaStNew <- function(Thet, M, Mu0dl, Muvec, Sig0, bs2, Sigp0){
    Mitl <- M 
    #-- Current value
    BetadelCur <- Thet$Betadelt  
    #-- Propose a value
    #BetadelNew <- c(mvtnorm::rmvnorm(n=1, mean = BetadelCur, sigma = Sigp0)) 
    BetadelNew <- c(mvtnorm::rmvnorm(n=1, mean = Muvec, sigma = Sigp0)) # simulate a proposal from the prior
    
    #-- Compute the likelihood (proposal)
    ernew <- LikHijMat(BetadelNew, Thet, Mitl, Mu0dl, Sig0, bs2,1) - mvtnorm::dmvnorm(BetadelNew, mean = Muvec, sigma = Sigp0, log=T)
    
    #-- Compute the likelihood (old)
    erold <- LikHijMat(BetadelCur, Thet, Mitl, Mu0dl, Sig0, bs2,1) - mvtnorm::dmvnorm(BetadelCur, mean = Muvec, sigma = Sigp0, log = T)
    
    # #-- Compute the likelihood (proposal)
    # ernew <- LikHijMat(BetadelNew, Thet, Mitl, Muvec, Sig0, bs2, Thet$Cik.m) - mvtnorm::dmvnorm(BetadelNew, mean = Muvec, sigma = Sigp0, log=T)
    # #-- Compute the likelihood (old)
    # erold <- LikHijMat(BetadelCur, Thet, Mitl, Muvec, Sig0, bs2, Thet$Cik.m) - mvtnorm::dmvnorm(BetadelCur, mean = Muvec, sigma = Sigp0, log = T)
    
    u0 <- runif(n=1)
    #--- 
    #Fval <- ifelse(log(u0) - (ernew-erold) < 0, ernew, erold)
    accept0 <- 1*(log(u0) < (ernew - erold))
    
    list(res=c(accept0, c(BetadelNew*(accept0) + (1 - accept0)*BetadelCur)), prop = BetadelNew, evid = ernew - erold)
  }
  updateFunDeltaStNew <- cmpfun(updateFunDeltaStNew)
  
  
  #-----------------------------------------------------------------------------
  # Update things related to X 
  #-----------------------------------------------------------------------------
  #--- update Pik.m for epsilon_i
  #alpha.x <- 1.0
  updatePik.x <- function(Thet, alpha.x){
    nk <- numeric(nrow(Thet$Muk.x))
    for(k in 1:nrow(Thet$Muk.x)){ nk[k] = sum(Thet$Cik.x == k)}
    
    #update the Pis
    alpha.xnew <- alpha.x/nrow(Thet$Muk.x) + nk
    
    c(rdirichlet(n=1, alpha.xnew))
  }
  
  #--- update Cik.m for epsilon_i based on mixture of normals
  updateCik.x <- function(i, Thet){ # mixture of normals 
    #Kx = ncol(Thet$X) # J x J
    #id <- Thet$Cik.x[i]
    #Sig <- matrix(Thet$Sigk.x[id,], ncol = length(Thet$Gamma), byrow = F)
    Xi = Thet$X[i,]
    #pi.k  = Thet$Pik.x*mapply(dmvnormVec, k = 1:Kx, x = Xi,  Mean = Thet$Muk.x , SigmaVc = Thet$Sigk.x)
    pi.k  = Thet$Pik.x*dmvnormVec(x=Xi,Mean = Thet$Muk.x , SigmaVc = Thet$Sigk.x) + .Machine$double.eps
    which(rmultinom(n=1,size=1,prob = pi.k) == 1)
  }
  #updateCik.x <- Vectorize(updateCik.x,"i")
  
  #--- Update muk.w based on the mixture of normal prior
  Mu0.x = numeric(pn); Sig0.x = .5*diag(pn); inSig0.x = 2*diag(pn)#as.inverse(.25*cov(Wn)) #.5*diag(pn) #
  updateMuk.x <- function(Thet, Mu0.x, inSig0.x, Sig0.x){
    K = length(Thet$Pik.x)
    J = ncol(Thet$X)
    Res0 <- matrix(0, nrow = K, ncol = J)
    Mu.k <- matrix(0, nrow = K, ncol = J) # K x J
    nk = nrow(Thet$X)
    #SigKWg <- matrix(0, ncol = ncol(Y),nrow=ncol(Y)) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    
    TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.x, ncol=J, byrow=F)))
    SigK <- as.inverse(as.symmetric.matrix((nk*TpSi + inSig0.x))) ## sum of matrices inverse
    Muk <- SigK%*%(crossprod(TpSi, c(colSums(Thet$X))) + inSig0.x%*%Mu0.x)      ## some of mat
    c(mvtnorm::rmvnorm(n=1, mean = Muk, sigma = SigK))
    
    #--- we impose some constrainsts on the mean
    
    
  }
  updateMuk.x <- cmpfun(updateMuk.x)
  
  #--- Update Sig2k.w based on the mixture of normal prior
  nu0.x <- pn + 2 ; Psi0.x <- .5*diag(pn) #as.symmetric.matrix(.25*cov(Wn))
  updateSigk.x <- function(Thet, nu0.x, Psi0.x){
    
    Ku = nrow(Thet$Muk.x)
    m = ncol(Thet$Muk.x)
    nk = nrow(Thet$X)
    
    Xtl <- Thet$X - matrix(c(Thet$Muk.x), ncol = ncol(Thet$X), nrow=nk, byrow = T)
    nux <- nu0.x  + nk
    #Yi0 <- matrix(c(Xtl[id,]), nrow=nk, byrow=F)
    #Psin.x <- Psi0.x + crossprod(Yi0)
    Psin.x <- Psi0.x + crossprod(Xtl)
    rinvwishart(nux, Psin.x)
    
  }
  updateSigk.x <- cmpfun(updateSigk.x)
  # 
  # A <- min(c(Wn, Mn)) - .2*diff(range(c(Wn, Mn))) 
  # B <- max(c(Wn, Mn)) + .2*diff(range(c(Wn, Mn))) 
  
  A <- min(c(Wn, Mn)) - .05*diff(range(c(Wn, Mn))) 
  B <- max(c(Wn, Mn)) + .05*diff(range(c(Wn, Mn))) 
  
  
  #--- Update X_i  
  
  udapateXi2 <- function(Thet, Y, W, M, A = A, B = B){
    
    J = ncol(Thet$X)
    n = length(Y)
    Ytl = (Y - Thet$muk.e[Thet$Cik.e])/Thet$sig2k.e[Thet$Cik.e]
    #Ytl = Y / Thet$sig2k.e[Thet$Cik.e]
    Wtl = W 
    Mtl = M #%*%diag(Thet$Betadel)
    Res = matrix(0, ncol=J, nrow = n)
    Sigw = matrix(0, J, J)
    Sigm = matrix(0, J, J)
    Sigx = matrix(0, J, J)
    GamMat= tcrossprod(c(Thet$Gamma))
    
    
    Sigw <- as.inverse(as.symmetric.matrix(Thet$Sigk.w))
    Sigm <- as.inverse(as.symmetric.matrix(Thet$Sigk.m))
    Sigx <- as.inverse(as.symmetric.matrix(Thet$Sigk.x))
    
    
    Mux <- matrix(c(Thet$Muk.x), ncol = ncol(W), nrow=length(Y), byrow=T)
    Sigxinv <- as.symmetric.matrix(Thet$Sigk.x)
    
    Mu <- matrix(NA, ncol = ncol(Thet$X), nrow= nrow(Thet$X))
    
    nk <- unique(Thet$Cik.e)
    for(k in 1:length(nk)){
      id <- which(Thet$Cik.e == nk[k])
      if(length(id) > 0){
    Sig <- as.positive.definite(as.inverse(as.symmetric.matrix(GamMat/Thet$sig2k.e[k] + Sigw +  Sigm + Sigx))) + 0.0001*diag(J)
    
    Mu[id,] <- (matrix(c(Ytl[id]), ncol=1) %*% matrix(Thet$Gamma,nrow=1) + tcrossprod(Wtl[id,],Sigw) + tcrossprod(Mtl[id,], Sigm) + 
             tcrossprod(Mux[id,],Sigx)) %*% Sig
      }
    
    
    for(i in id){
      Res[i,] = tryCatch(expr={
        c(TruncatedNormal::rtmvnorm(n=1, mu = Mu[i,], sigma =  Sig,lb = rep(A,J), ub = rep(B,J)))
      },
      warning= function(W){
        if(grepl("Did not find a solution to the nonlinear system in `mvrandn`", W$message)){
          
          c(TruncatedNormal::rtmvnorm(n=1, mu = numeric(J), sigma =  Sig,lb = rep(A,J), ub = rep(B,J)))
          
        }
      }
      )  %>% as.numeric()
    }
    }
    
    Res
    
  }
  udapateXi2 <- cmpfun(udapateXi2)
  
  
  #-- Update Xi new
  #'@ SigP0 = prior covariance matrix
  #'@ Mup0 = prior Mean 
  #'@ SigProp = proposal covariance matrix
  #'@ Muprop = proposal mean 
  
  ##----------------------------------------------------------------------------
  cat("\n #---- Create containers for the MCMC output ---> \n")
  #-----------------------------------------------------------------------------
  #----- Divide these in sections - Construct the output
  #----- containers
  #--------------------------------------------------
  #Nsim <- 5000
  
  MCMCSample <- list()
  
  #-- for X
  Kx = 1
  MCMCSample$X <- matrix(0,nrow = Nsim,ncol=length(c(Thet$X))) ## c(t()) : how we transform the matrix (row combined) 
  MCMCSample$Muk.x <- matrix(0,nrow = Nsim, ncol = Kx*pn)
  MCMCSample$Sigk.x <- matrix(0,nrow = Nsim,ncol=Kx*pn*pn)
  #MCMCSample$Pik.x <- matrix(0,nrow = Nsim,ncol=Kx)
  #MCMCSample$Cik.x <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  
  #-- for W
  Kw = 1
  MCMCSample$Muk.w <- matrix(0,nrow = Nsim,ncol=Kw*pn)
  MCMCSample$Sigk.w <- matrix(0,nrow = Nsim,ncol=Kw*pn*pn)
  #MCMCSample$Pik.w <- matrix(0,nrow = Nsim,ncol=Kw)
  #MCMCSample$Cik.w <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  
  #-- for M
  Km=1
  MCMCSample$Muk.m <- matrix(0,nrow = Nsim,ncol=Km*pn)
  MCMCSample$Sigk.m <- matrix(0,nrow = Nsim,ncol=Km*pn*pn)
  #MCMCSample$Pik.m <- matrix(0,nrow = Nsim,ncol=Km)
  #MCMCSample$Cik.m <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  MCMCSample$Delta <- numeric(length(Thet$Delta))
  MCMCSample$Betdelt <- matrix(0, nrow=Nsim, ncol=length(Thet$Betadelt))
  MCMCSample$AcceptDelt <- numeric(Nsim)
  
  # For Y
  MCMCSample$muk.e <- matrix(0,nrow = Nsim, ncol=length(Thet0$muk.e)) ## c(t()) : how we transform the matrix (row combined)
  MCMCSample$sig2k.e <- matrix(0,nrow = Nsim, ncol=length(Thet0$muk.e))
  MCMCSample$Pik.e <- matrix(0,nrow = Nsim, ncol=Ke)
  MCMCSample$Cik.e <- matrix(0,nrow = Nsim, ncol=nrow(Wn))
  
  MCMCSample$Gamma <- matrix(0,nrow = Nsim, ncol=length(Thet0$Gamma))
  MCMCSample$BetaZ <- matrix(0,nrow = Nsim, ncol=length(Thet$BetaZ))
  MCMCSample$siggam <- numeric(Nsim)
  
  #------------------------------------------------------------------
  #Nsim = 1000 
  #library(profvis)
  #profvis(
  Sigp000 <- diag(rep(sig2.propdelt, length(Thet0$Betadel))) #solve(OptParm$hessian)
  Mu0dl <- Thet0$Betadel
  MCMCSample$Propval <- matrix(NA, nrow = Nsim, ncol = length(Thet0$Betadel))
  MCMCSample$evid <- NA*numeric(Nsim)
  Out0 <- NA
  
  cat("\n #--- MCMC is starting ----> \n")
  MCMCSample$PsiAcceptF <- numeric(6)
  for(iter in 1:Nsim){
    if(iter %% 500 == 0){cat("iter=",iter,"\n")}
    #----------------------------------------------
    # Update things related to X 
    #--------------------------------------------  
    Thet$Cik.x <-  rep(1, nrow(Wn)) 
    
    Thet$Muk.x <- updateMuk.x(Thet, Mu0.x, inSig0.x, Sig0.x)
    
    Thet$Sigk.x <- updateSigk.x(Thet, nu0.x, Psi0.x)
    
    Ytd = Y - Zmod%*%Thet$BetaZ
    
    #if( (iter %% lag == 0)*is.null(X)){
    Thet$X <- udapateXi2(Thet, Ytd, Wn, Mn, A=A, B=B)
    #}
    
    #---------------------------------------------
    #--- update things related to W
    #---------------------------------------------
    #--- Update Pik.w  
    # #updatePik.w <- function(Thet, alpha.W){
    # nk <- numeric(Kw)
    # for(k in 1:Kw){ nk[k] = sum(Thet$Cik.w == k)}
    # #update the Pis
    # alpha.wnew <- alpha.W/Kw + nk
    # #Thet$Pik.w <- c(1, rep(0, Kw-1)) #rdirch(n=1, alpha.wnew)
    # Thet$Pik.w <- rdirch(n=1, alpha.wnew)
    
    #--- Update Ci.k mixture of normals
    #Thet$Cik.w <- mapply(updateCik.w, i= 1:N, MoreArgs = list(Thet=Thet, W = Wn))
    #Thet$Cik.w <- rep(1, nrow(Wn)) #F2Cmp(CikW(Wn, Thet$X, Thet$Muk.w, Thet$Sigk.w, Thet$Pik.w)) # rep(1, nrow(Wn)) #
    # if(Gx0){Thet$Cik.w <- rep(1, length(Thet$Cik.w))}else{
    #   Thet$Cik.w <- F2Cmp(CikW(Wn, Thet$X, Thet$Muk.w, Thet$Sigk.w, Thet$Pik.w))} # rep(1, nrow(Wn)) #
    # 
    #--- Update muk.w and sig2k.w based on the mixture of normal prior
    #Thet$Muk.w <- updateMuk.w(Thet, Wn, Mu0.w, inSig0.w, Sig0.w)
    Thet$Sigk.w <- updateSigk.w(Thet, Wn, nu0.w, Psi0.w)
    
    #---------------------------------------------
    #--- update things related to M
    #---------------------------------------------
    Thet$Sigk.m <- updateSigk.m(Thet, M, nu0.m, Psi0.m, bs2=bs2)
    
    Out0 <- updateFunDeltaStNew(Thet, M, Mu0dl, mu0.del, Sig0, bs2, Sigp0 = PsiAccept * Sigp000) #Sigp000 PsiAccept
    
    #--- testing the code
    MCMCSample$Propval[iter, ] <- Out0$prop 
    Thet$Betadelt <- Out0$res[-1]
    MCMCSample$Betdelt[iter,] <- Thet$Betadelt
    MCMCSample$AcceptDelt[iter] <- Out0$res[1]
    MCMCSample$evid[iter] <- Out0$evid
    
    #---------------------------------------------------------------------------
    #--- update things related to Y or epsilon
    #--------------------------------------------------------------------------- 
    #---------------------------------------------
    #--- update Pik.e for epsilon_i
    #---------------------------------------------
    nk <- numeric(length(Thet$muk.e))
    for(k in 1:length(Thet$muk.e)){ nk[k] = sum(Thet$Cik.e == k)}
    #update the Pis
    alpha.enew <- alpha.e/length(Thet$muk.e) + nk
    Thet$Pik.e <- c(rdirch(n=1, alpha.enew))
    
    #update the Cik.e
    Thet$Cik.e <- F2Cmp(Cike(c(Y), Thet$X, Thet$muk.e, Thet$sig2k.e, Thet$Pik.e,Thet$Gamma))
    #update muk.e
    #Thet$muk.e <- updatemuk.e(Thet, Y, Z, mu0.e, sig20.e)
    
    #update Sigk.e
    Thet$sig2k.e <- updatesig2k.e(Thet,Y, Z, gam0.e, sig20.esig)
    #---------------------------------------------
    #--- update Gamma and BetaZ
    #---------------------------------------------
    Tp0 <- updateGam(Thet, Y, Z, Sig0.gam0, Mu0.gam0)
    Thet$BetaZ <- Tp0[c(1:ncol(Zmod))] #c(0, Tp0[c(2:ncol(Zmod))])  #Tp0[c(1:ncol(Zmod))]
    Thet$Gamma <- Tp0[-c(1:ncol(Zmod))]
    
    if(is.null(sig02.del)){
      Thet$siggam <- updateSigGam(Thet$Gamma, al0, bet0, order = 2)  # updateSigGam <- function(Bet,al0, bet0, order = 2) ##1/tauval #
    }
    #---------------------------------------------------------------------------
    #------------------------- Save output
    #---------------------------------------------------------------------------
    #-- For X
    MCMCSample$Muk.x[iter, ] <- c(t(Thet$Muk.x))  ## Stor rowwise
    MCMCSample$Sigk.x[iter, ] <- c(t(Thet$Sigk.x))
    #MCMCSample$Pik.x[iter, ] <- Thet$Pik.x
    #MCMCSample$Cik.x[iter, ] <- Thet$Cik.x
    MCMCSample$X[iter, ] <- c(t(Thet$X)) ## Store row-wise
    
    #-- for W
    MCMCSample$Muk.w[iter, ] <- c(t(Thet$Muk.w))
    MCMCSample$Sigk.w[iter, ] <- c(t(Thet$Sigk.w))
    #MCMCSample$Pik.w[iter, ] <- Thet$Pik.w
    #MCMCSample$Cik.w[iter, ] <- Thet$Cik.w
    
    #-- for M
    MCMCSample$Muk.m[iter, ] <- c(t(Thet$Muk.m))
    MCMCSample$Sigk.m[iter, ] <- c(t(Thet$Sigk.m))
    #MCMCSample$Pik.m[iter, ] <- Thet$Pik.m
    #MCMCSample$Cik.m[iter, ] <- Thet$Cik.m
    #MCMCSample$AcceptDelt[iter] <- Out0[1]
    #MCMCSample$Delta[iter] <- Thet$Delta
    #-- For Y
    MCMCSample$muk.e[iter, ] <- Thet$muk.e ## c(t()) : how we transform the matrix (row combined)
    MCMCSample$sig2k.e[iter, ] <- c(Thet$sig2k.e)
    MCMCSample$Pik.e[iter, ] <- Thet$Pik.e
    MCMCSample$Cik.e[iter, ] <- Thet$Cik.e
    
    MCMCSample$Gamma[iter, ] <- Thet$Gamma
    MCMCSample$BetaZ[iter, ] <- Thet$BetaZ
    MCMCSample$siggam[iter] <- Thet$siggam
  }  
  
  #---- Return the MCMC smaples
  cat("\n ---- MCMCs are done!! ---/n")
  #-- Same the raw (pointwise) estimate of delta(t)
  MCMCSample$Delthat <- colMeans(M)/colMeans(W) #Delthat
  
  #-- Same the smooth estimate of delta(t)
  MCMCSample$EstimDelthat <- (bs2del%*%Thet0$Betadel)
  
  #--- Save the initial values
  MCMCSample$Thet0 <- Thet0
  #MCMCSample$OptimOut <- OptParm
  MCMCSample
  
  
}
#===============================================================================
#==== Multiple on components of the truncated DP
BayesFIVMult = function(Nsim = 100, Y, W, M, Z, sig02.del = 1.0, tauval = 1, 
                        alpx = 1, Gx = 2, Ge = 2, Gmw = 2, BasFun = "Bspline", PsiAccept = 0.05, sig2.propdelt = 0.1, sig2.priordelt = 5.0 ,lag = 1, X = NULL){
  #--- Functions
  func <- function(Pik, n=1, size=1){
    which.max(rmultinom(n=1,size=1,prob = Pik + abs(min(Pik))+ 0.001))
  }
  funcCmp <- cmpfun(func)
  
  F2 <- function(Pik, n=1, size=1){
    
    #print(Pik)
    apply(Pik,1,funcCmp)
  }
  F2Cmp <- cmpfun(F2)
  
  a <-  seq(0, 1, length.out = ncol(M))
  n = length(Y)
  cat("--- set some initial values ")
  #---- Parameters settings
  #-------------------------------------
  pn <- round(n^{1/3}) + 2 ## min(15, ceiling(length(Y)^{1/3.8})+4)
  
  alpha.e = alpx
  alpha.w = alpha.W = alpx
  alpha.x = alpha.X = alpx
  alpha.m = alpx
  pn = pn 
  n = length(Y) 
  Zmod = Z
  #--- Transform the data
  #-----------------------------------------------------------------------------
  
  bs2 <- bs(a, df = pn, intercept = T)
  #bs2 <- qr.Q(qr(bs20))
  bs2del <- qr.Q(qr(bs2))
  
  #--- Using different bases
  
  if(BasFun == "Four"){
    
    #Nbasis0 = c(3,4, 5, 7) + 2
    J = K = k  = pn = pn00 =  Nbasis0[which(c(100,200,500, 1000) == n)]
    #Nbasis0 =  c(3, 6, 15, 18)
    J = K = k  = pn = pn00 =  Nbasis0[which(c(100,200,500, 1000) == n)]
    
    
    Fourmeth <- fda::create.fourier.basis(c(0,1), nbasis = pn)
    bs2 = as.matrix(eval.basis(a, Fourmeth))
    colnames(bs2) <- NULL
    bs2del <- qr.Q(qr(bs2)) 
    
  }
  
  Mn <- (M%*%bs2)/length(a) #(M%*%bs2)/length(a) 
  Wn <- (W%*%bs2)/length(a)
  if(!is.null(X)){
    Xn = (X%*%bs2)/length(a)
  }
  
  P.mat <- function(K){
    # penalty matrix
    D <- diag(rep(1,K))
    D <- diff(diff(D))
    P <- t(D)%*%D 
    return(P)
  }
  
  
  InitialValues <- function(Y, W, M, Z, df0=3, Gmw=Gmw, Gx, a, Ge, bs2){
    Res <- list()
    #bs.w <- splines::bs(a, df = pn+1, degree = min(c(df0,pn-1)), intercept=TRUE) 
    Zmod <- Z #cbind(1,Z)
    pn = ncol(bs2)
    
    
    bs2del <- qr.Q(qr(bs2))
    
    Wn <- crossprod(t(W), bs2)/length(a)
    Mn <- crossprod(t(M), bs2)/length(a)
    
    WsmB <- Wn
    WsM <- NULL
    for(i in 1:nrow(W)){ WsM <- rbind(WsM, smooth.spline(a, W[i,])$y)}
    W_i <- crossprod(t(WsM),bs2)/length(a)   #(n,k) matrix
    lmFit <- lm(Y ~ -1 + Zmod + W_i) 
    c_hatW <- lmFit$coefficients
    Res$Gamma <- as.numeric(c_hatW[-c(1:ncol(Zmod))])
    Res$BetaZ <- as.numeric(c_hatW[c(1:ncol(Zmod))])
    Res$Delta <- abs(mean(Mn)/mean(W_i))
    #Res$Betadel <- c(abs(solve(t(bs2del)%*%bs2del)%*%t(bs2del)%*%(colMeans(M)/colMeans(W))))
    
    #--- Estimate Delta
    Delhat <- function(bs, U0){
      
      Fe0 <- function(val){
        sum((bs%*%val - U0)^2)
      }
      
      k = ncol(bs)
      Res <- optim( rep(0, ncol(bs)), Fe0, method="L-BFGS-B", lower=rep(-60,k), upper=rep(60, k) )
      Res  
      
    }
    U0 = colMeans(M)/colMeans(W)
    R0d <- Delhat(bs2del, U0)
    
    Res$Betadel <- R0d$par 
    
    #--- Determine the number of clusters for epsilon_i
    reslmfit <- residuals.lm(lmFit)
    BIC <- mclustBIC(reslmfit)
    mod1 <- Mclust(reslmfit, x = BIC, G = Ge)
    Temp <- summary(mod1, parameters = TRUE)
    Res$Ke <- Temp$G #+ Gx
    Res$muk.e <- numeric(Res$Ke)
    Res$muk.e[1:Temp$G] = Temp$mean; 
    Res$sig2k.e <- rep(1, Res$Ke)*sum(lmFit$residuals^2)/lmFit$df.residual
    Res$sig2k.e[1:Temp$G] = Temp$variance
    Res$Cik.e <- as.numeric(Temp$classification)
    Res$Pik.e <- numeric(Res$Ke)
    Res$Pik.e[1:Temp$G] <- Temp$pro 
    #-----------------------------------------------------
    #--- Determine the number of cluster (W)
    #------ 
    Gmw0 <- ifelse(Gmw > 1, Gmw, 2)
    Ui <- Wn #- W_i 
    Mod1 <- Mclust(Ui, G = Gmw)
    Tp <- summary(Mod1, parameters=T)
    Res$Kw <- Tp$G #+Gx
    Res$Muk.w <- matrix(0, ncol = pn, nrow=Res$Kw)
    Res$Muk.w <- as.matrix( rbind(t(Tp$mean))) #  ,  repmat(rep(0,pn), Gx, 1)))
    
    Res$Sigk.w <- matrix(0, nrow = Res$Kw, ncol = pn*pn)
    for(j in 1:Tp$G){Res$Sigk.w[j,] <- c(Tp$variance[,,j])}
    
    Res$Cik.w <- Tp$classification
    Res$Pik.w <- Tp$pro
    
    #-------------------------------------------------   
    #--- Determine the number of cluster (X)
    #-----
    WsmB <- (W + M) /(1+Res$Delta)
    WsM <- NULL
    for(i in 1:nrow(W)){ WsM <- rbind(WsM, smooth.spline(a, WsmB[i,])$y)}
    W_i <- crossprod(t(WsM),bs2)/length(a)   #(n,k) matrix
    Ui <- W_i 
    Mod1 <- Mclust(Ui, G=Gx)
    Tp <- summary(Mod1, parameters=T)
    Res$Kx <- Tp$G #+ Gx
    Res$Muk.x <- matrix(0, ncol = pn, nrow=Res$Kw)
    Res$Muk.x <- as.matrix( rbind(t(Tp$mean),  repmat(colMeans(Ui), Gx, 1 )))
    
    Res$Sigk.x <- matrix(0, nrow = Res$Kx, ncol = pn*pn)
    for(j in 1:Tp$G){Res$Sigk.x[j,] <- c(Tp$variance[,,j])}
    #for(j in (Tp$G+1):Res$Kx){Res$Sigk.x[j,] <- c(.5*diag(pn))}
    Res$Cik.x <- Tp$classification
    Res$Pik.x <- Tp$pro  
    Res$X <- W_i
    #---------------------------------------------------
    #---- Determine the number of cluster (M)
    #------
    Gmw0 <- ifelse(Gmw > 1, Gmw, 2)
    Ui <- Mn  #- Res$Delta*W_i 
    Mod1 <- Mclust(Ui, G=Gmw0)
    Tp <- summary(Mod1, parameters=T)
    Res$Km <- Tp$G #+ Gx
    #Res$Muk.m <- as.matrix(t(Tp$mean))
    Res$Muk.m <- as.matrix( rbind(t(Tp$mean))) #,  repmat(colMeans(Ui), Gx, 1 )))
    Res$Sigk.m <- matrix(0,nrow=Res$Km,ncol=pn*pn)
    for(j in 1:Tp$G){Res$Sigk.m[j,] <- c(Tp$variance[,,j])}
    #for(j in (Tp$G+1):Res$Km){Res$Sigk.m[j,] <- c(as.positive.definite(.5*diag(pn)))}
    Res$Cik.m <- Tp$classification
    Res$Pik.m <- Tp$pro
    Res$R0d <- R0d
    #------- return result 
    Res 
  }
  
  
  cat("\n #--- Obtain initial values...\n")
  
  Thet0 <- InitialValues(Y, W, M, Z, df0=3, Gmw = Gmw, Gx=Gx, a, Ge, bs2)
  
  #---------------------------------------------
  cat("\n#---- Initialized parameters ---> \n")
  #---------------------------------------------
  
  Kx = Thet0$Kx               ## number of X  cluster
  Ke = Thet0$Ke               ## number of ei cluster
  Km = Thet0$Km               ## number of ei cluster
  Kw = Thet0$Kw               ## number of omega_i cluster 
  
  
  Mn0 <- (M%*%bs2)/length(a) #(M%*%bs2)/length(a) 
  Wn <- (W%*%bs2)/length(a)
  Mn <- ((M%*%diag(1/(c(bs2del%*%Thet0$Betadel))))%*%bs2)/length(a) 
  
  Thet <- list()
  Thet$Gamma <-  rnorm(n=length(Thet0$Gamma)) #as.numeric(Thet0$Gamma)   #- a vector of length J
  Thet$BetaZ <- Thet0$BetaZ
  Thet$siggam <- ifelse(is.null(sig02.del), 1, sig02.del)
  Thet$bs2del <- bs2del
  Thet$Betadelt <-  Thet0$Betadel
  
  Thet$Cik.x =  Thet0$Cik.x #sample(x = 1:Kx,size = n, replace = T)    #- vector of length Kx
  Thet$Cik.e =  Thet0$Cik.e #sample(x = 1:Ke,size = n, replace = T)     #- vector of length Ku
  Thet$Cik.w = Thet0$Cik.w  #sample(x = 1:Kw,size = n, replace = T)     #- vector of length Kw
  Thet$Cik.m = Thet0$Cik.m #sample(x = 1:Km,size = n, replace = T)     #- vector of length Km
  
  
  Thet$Pik.x = Thet0$Pik.x #rdirch(n=1,alpha = rep(alpha.X, Kx)/Kx)   #- vector of length Kx
  Thet$Pik.e = Thet0$Pik.e #rdirch(n=1,alpha = rep(alpha.e, Ke)/Ke)  #- vector of length Ku
  Thet$Pik.w = Thet0$Pik.w #rdirch(n=1,alpha = rep(alpha.W, Kw)/Kw)  #- vector of length Kw
  Thet$Pik.m = Thet0$Pik.m #rdirch(n=1,alpha = rep(alpha.m, Km)/Km)  #- vector of length Kw
  
  
  Thet$sig2k.e  =  Thet0$sig2k.e  #rgamma(n=Ke, shape=1,scale=.1) # matrix dim Kx x J
  Thet$muk.e  =  Thet0$muk.e #rnorm(n=Ke)  #- Matrix of (Ku x J) and Ku x J - Note that $\Omega$ is a diagonal matrix
  
  Thet$Sigk.w <- repmat(c(.5*diag(pn)), n = Kw, m = 1)     #- vector of length Kw
  Thet$Muk.w <- matrix(rnorm(n = Kw*pn), nrow = Kw)          #- vectors of length Kw and Kw respectively
  
  Thet$Sigk.x <- repmat(c(.5*diag(pn)), Kx, 1)     #-- vector of length Kx
  Thet$Muk.x  <- matrix(rnorm(n=Kx*pn), nrow=Kx)   #- vectors of length Kx and Kx respectively
  
  Thet$Sigk.m <- repmat(c(.5*diag(pn)), Km, 1)     #-- vector of length Kx
  Thet$Muk.m  <- matrix(rnorm(n=Km*pn), nrow = Km)   #- vectors of length Kx and Kx respectively
  
  
  A <- min(c(c(Wn), c(Mn))) - 0.05*diff(range(c(Wn, Mn)))
  B <- max(c(Wn, Mn)) + 0.05*diff(range(c(Wn, Mn)))
  
  Thet$X <- TruncatedNormal::rtmvnorm(n=n,mu=rep(0,pn), sigma=diag(pn),lb = rep(A, pn), ub = rep(B, pn))
  if(!is.null(X)){
    Thet$X <- Xn
  }
  
  cat("\n#---- Loading updating functions ---> \n")
  #-----------------------------------------------------------------------------
  # Update things related to Y
  #-----------------------------------------------------------------------------
  pn0 = ncol(Z)
  Pgam <- P.mat(pn)*tauval
  Sig0.gam0 = diag(pn0); Mu0.gam0 = numeric(pn0+pn);
  
  
  updateGam <- function(Thet, Y, Z, Sig0.gam0, Mu0.gam0){
    
    Zmod <-  Z #cbind(1,Z)
    Y.td = Y - Thet$muk.e[Thet$Cik.e]
    InvSig0.gam0 <- 1.0*as.matrix(bdiag(as.inverse(Sig0.gam0), P.mat(length(Thet$Gamma))/Thet$siggam))
    Xmod <- cbind(Zmod, Thet$X)
    X_sc <- sweep(Xmod, MARGIN = 1, Thet$sig2k.e[Thet$Cik.e], FUN="/") ## n*pn
    
    Sig.gam <- as.inverse(as.symmetric.matrix(crossprod(X_sc, Xmod) + InvSig0.gam0))
    Mu.gam <- c(Sig.gam%*%(crossprod(X_sc, Y.td) + InvSig0.gam0%*%Mu0.gam0))
    
    c(mvtnorm::rmvnorm(n=1, mean = Mu.gam, sigma = as.positive.definite(Sig.gam)))
    
  }
  updateGam <- cmpfun(updateGam)
  
  #siggam ~ IG(al0, bet0)
  al0 = 1
  bet0 = .005
  updateSigGam <- function(Bet, al0, bet0, order = 2){
    K = length(Bet)
    al = al0 + .5*(K - order)
    bet = .5*quad.form(P.mat(K), Bet) + bet0
    #1/rgamma(n = 1, shape = al, rate = bet)
    #nimble::rinvgamma(n = 1, shape = al, scale = bet)
    # tb = 10
    # 1/heavy::rtgamma(n=1,shape=al,scale = bet, t = tb)
    1/cascsim::rtgamma(n=1,shape=al, scale = 1/bet, max = 10, min = 1/15)
  }
  
  #------- Update Pik.e 
  #alpha.e = 1.0
  updatePik.e <- function(Thet, alpha.e){
    
    nk <- numeric(length(Thet$muk.e))
    for(i in 1:length(Thet$muk.e)){ nk[i] = sum(Thet$Cik.e == i)}
    
    #update the Pis
    alpha.enew <- alpha.e/length(Thet$muk.e) + nk
    
    rdirichlet(n=1, alpha.enew)
    
  }
  
  #--- Update muk.w and sig2k.w based on the mixture of normal prior Xtl
  mu0.e = 0 ; sig20.e = 100
  updatemuk.e <- function(Thet, Y, Z, mu0.e, sig20.e ){
    
    Ku = length(Thet$muk.e)
    muk <- numeric(Ku)
    sig2k <- numeric(Ku)
    Zmod <- Z #cbind(1,Z)
    #Ytl <- Y - Thet$X%*%Thet$Gamma
    Ytl = Y - Zmod%*%Thet$BetaZ -  Thet$X%*%Thet$Gamma
    for(k in 1:length(Thet$muk.e)){
      id <- which(Thet$Cik.e == k)
      nk = length(id)
      if(nk > 0){
        sig2k[k] <- 1/(nk/Thet$sig2k.e[k] + 1/sig20.e)
        muk[k] <- sig2k[k]*(mu0.e/sig20.e + sum(c(Ytl[id]))/Thet$sig2k.e[k]) 
      }else{
        sig2k[k] <- sig20.e
        muk[k] <- mu0.e 
      }
    }
    ## Force the zero mean constrains
    SigR0 <- Thet$Pik.e*sig2k
    SigRR <- sum((Thet$Pik.e^2)*sig2k)
    MuK0 = muk - SigR0*(1/SigRR)*sum(Thet$Pik.e*muk)
    SigK <- as.matrix(nearPD(diag(sig2k) - tcrossprod(SigR0)*(1/SigRR), doSym = T)$mat)
    
    if(Ku > 2){
      id0 = 1:(Ku-1)
      Muk_1 <- c(mvtnorm::rmvnorm(n=1,mean=MuK0[-Ku], sigma = SigK[id0,id0]))
      
      c(Muk_1, -sum(c(Thet$Pik.e[-Ku]*Muk_1)/ifelse(Thet$Pik.e[Ku] == 0.0, (.Machine$double.eps)^{.7}, Thet$Pik.e[Ku])))
    }else{
      
      id0 = 1:(Ku-1)
      Muk_1 <- c(rnorm(n=1,mean=MuK0[-Ku], sd = sqrt(SigK[id0,id0])))
      
      c(Muk_1, -(c(Thet$Pik.e[-Ku]*Muk_1)/ifelse(Thet$Pik.e[Ku] == 0.0, (.Machine$double.eps)^{.7}, Thet$Pik.e[Ku])))
    }
    
  }
  
  gam0.e = 1 ; sig20.esig = 1
  #gam0.e = 1 ; sig20.e = .5
  updatesig2k.e <- function(Thet, Y, Z, gam0.e, sig20.e){
    
    Zmod <- Z #cbind(1,Z)
    Ku = length(Thet$muk.e)
    al <- numeric(Ku)
    bet <- numeric(Ku)
    Xtl <- Y - Thet$X%*%Thet$Gamma - Zmod%*%Thet$BetaZ - Thet$muk.e[Thet$Cik.e]
    for(k in 1:length(Thet$muk.e)){
      id <- which(Thet$Cik.e == k)
      nk = length(id)
      if(nk > 0){
        al[k] <- nk/2 + gam0.e
        bet[k] <- sig20.e + .5*sum(Xtl[id]^2) 
      }else{
        al[k] <- gam0.e
        bet[k] <- sig20.e 
      }
    }
    ## Force the zero mean constrains
    mapply(rinvgamma, n = 1, shape = al, scale = bet)
  }
  
  #-----------------------------------------------------------------------------
  # Update things related to W 
  #-----------------------------------------------------------------------------
  
  #--- Update Pik.x  
  #alpha.W = alpha.w = 1.0
  updatePik.w <- function(Thet, alpha.W){
    Ku = nrow(Thet$Muk.w)
    nk <- numeric(Ku)
    for(k in 1:Ku){ nk[k] = sum(Thet$Cik.w == k)}
    
    #update the Pis
    alpha.wnew <- alpha.W/length(Thet$muk.w) + nk
    
    rdirichlet(n=1, alpha.wnew)
    
  }
  
  Mu0.w = rep(0, pn);  Sig0.w = .5*diag(pn);inSig0.w = 2*diag(pn) 
  updateMuk.w <- function(Thet, W, Mu0.w, inSig0.w, Sig0.w){
    Ku = length(Thet$Pik.w)
    J = ncol(Thet$X)  
    Mu.k <- matrix(0, nrow = Ku, ncol = J) # K x J
    SigK <- list()  #-- store the covariances
    SigKWg <- matrix(0, ncol=J, nrow=J) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    Yi.tld = W - Thet$X 
    #Muk <- matrix(0, ncol=J,nrow=Ku)
    SigK <- list()
    for(k in 1:Ku){
      idk <- which(Thet$Cik.w == k)
      nk <- length(idk)
      if(length(idk) > 0){
        Yi0 = matrix(c(Yi.tld[idk,]), nrow= nk, byrow=F)
        TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.w[k,],ncol=J,byrow=F)))
        SigK[[k]] <- as.inverse(as.symmetric.matrix(nk*TpSi + inSig0.w)) ## sum of matrices inverse
        SigKWg <- SigKWg + (Thet$Pik.w[k]*Thet$Pik.w[k])*SigK[[k]] # J x J
        #Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi.tld[idk,])) ) + solve(Sig0.w)%*%Mu0.w)      ## some of mat
        Mu.k[k,] <- c(SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi0)) ) + inSig0.w%*%Mu0.w))      ## some of mat
        #Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(apply(t(Yi.tld[idk,]), 2, sum)) ) + solve(Sig0.w)%*%Mu0.w)      ## some of mat
      }else{
        SigK[[k]] <- Sig0.w
        SigKWg <- SigKWg + (Thet$Pik.w[k]*Thet$Pik.w[k])*SigK[[k]] # J x J
        Mu.k[k,] <- Mu0.w
      }
    }
    #--- we impose some constrainsts on the mean
    SIg0 <- bdiag(SigK) # create a block diagonal matrix
    SIGR0 <- NULL; for(k in 1:Ku){SIGR0 <- rbind(SIGR0,Thet$Pik.w[k]*SigK[[k]])} ## K*J x J
    MUR0 <-  colSums(diag(Thet$Pik.w) %*% Mu.k) # return a evector of length J
    
    InvSigKWg <- as.positive.definite(as.inverse(SigKWg))
    # print(dim(InvSigKWg))
    # print(dim(SIg0))
    # print(dim(quad.tform.inv(SigKWg,SIGR0)))
    # print(Ku)
    #MURFin <- c(t(Mu.k)) - SIGR0%*%solve(SigKWg)%*%MUR0 ## J*K
    MURFin <- c(t(Mu.k)) - SIGR0%*%InvSigKWg%*%MUR0 ## J*K
    #SIGR0fin <- 1.0*as.matrix(nearPD(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T)$mat) ## J*K x J*K
    #SIGR0fin <- 1.0*(as.positive.definite(as.symmetric.matrix(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T))) ## J*K x J*K
    SIGR0fin <- 1.0*(as.positive.definite(as.symmetric.matrix(as.matrix(SIg0 - SIGR0%*%InvSigKWg%*%t(SIGR0)))))
    
    #--- Simulate K-1 vectors from the degenerate dist.
    if(Ku > 2){
      id <- 1:(J*(Ku-1))
      SimMU.k <- matrix(c(mvtnorm::rmvnorm(n=1, mean = MURFin[id], sigma = SIGR0fin[id,id])), ncol=J,byrow=T) ## (K-1) x J
      SimMU.k1 <- -colSums(diag(Thet$Pik.w[-Ku]) %*% SimMU.k)/Thet$Pik.w[Ku] ## get a K-1 x J matrix follow by colSums - a vector of J
      as.matrix(rbind(SimMU.k,SimMU.k1)) ## k * J
    }else{
      id <- 1:(J*(Ku-1))
      SimMU.k <- matrix(c(mvtnorm::rmvnorm(n=1, mean = MURFin[id], sigma = SIGR0fin[id,id])), ncol=J,byrow=T) ## (K-1) x J
      SimMU.k1 <- -(Thet$Pik.w[-Ku] * SimMU.k)/Thet$Pik.w[Ku] ## get a K-1 x J matrix follow by colSums - a vector of J
      as.matrix(rbind(SimMU.k,SimMU.k1)) ## k * J
      
    }
  }
  updateMuk.w <- cmpfun(updateMuk.w)
  
  #--- update Sigk.w
  nu0.w = pn + 2; Psi0.w = .5*diag(pn) 
  updateSigk.w <- function(Thet, W, nu0.w, Psi0.w){
    
    Ku = nrow(Thet$Muk.w)
    m = ncol(W)
    al <- numeric(Ku)
    bet <- numeric(Ku)
    Wtl <- W - Thet$X - Thet$Muk.w[Thet$Cik.w,]
    Res <- matrix(0,ncol = m*m ,nrow = Ku)
    for(i in 1:Ku){
      id <- which(Thet$Cik.w == i)
      nk = length(id)
      if(nk > 0){
        nun <- nu0.w  + nk
        Yi0 = matrix(c(Wtl[id,]), nrow= nk, byrow=F)
        Psin.w <- as.positive.definite(Psi0.w + crossprod(Yi0))
        Res[i, ] <- c(rinvwishart(nun, Psin.w))
      }else{
        Res[i, ] <- c(rinvwishart(nu0.w, Psi0.w))
      }
    }
    Res
  }
  updateSigk.w <- cmpfun(updateSigk.w)
  
  #-----------------------------------------------------------------------------
  # Update things related to M 
  #-----------------------------------------------------------------------------
  
  updatePik.m <- function(Thet, alpha.m){
    nk <- numeric(nrow(Thet$Muk.m))
    for(k in 1:nrow(Thet$Muk.m)){ nk[k] = sum(Thet$Cik.m == k)}
    
    #update the Pis
    alpha.mnew <- alpha.m/nrow(Thet$Muk.m) + nk
    
    c(rdirichlet(n=1, alpha.mnew))
  }
  
  Mu0.m = rep(0,pn) ; Sig0.m = .5*diag(pn); inSig0.m = 2*diag(pn) 
  updateMuk.m <- function(Thet, M, Mu0.m, inSig0.m, Sig0.m, bs2){
    Ku = length(Thet$Pik.m)
    J = ncol(bs2) #ncol(M)  
    Mu.k <- matrix(0, nrow = Ku, ncol = J) # K x J
    SigK <- list()  #-- store the covariances
    SigKWg <- matrix(0, ncol = J, nrow = J) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    
    Yi.tld <- ((M - (tcrossprod(Thet$X, bs2) * matrix(c(Thet$bs2del%*%Thet$Betadel), nrow=nrow(M), ncol = ncol(M), byrow=T))) %*% bs2) 
    
    Muk <- matrix(0, ncol=J, nrow=Ku)
    for(k in 1:Ku){
      idk <- which(Thet$Cik.m == k)
      nk <- length(idk)
      if(length(idk) > 0){
        Yi0 = matrix(c(Yi.tld[idk,]), nrow = nk, byrow=F)
        TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.m[k,],ncol=J, nrow=J,byrow=F)))
        SigK[[k]] <- as.inverse(as.symmetric.matrix(nk*TpSi + inSig0.m)) ## sum of matrices inverse
        #cat(dim(SigK[[k]]),"(--2 \n")
        SigKWg <- SigKWg + (Thet$Pik.m[k]*Thet$Pik.m[k])*SigK[[k]] # J x J
        Mu.k[k,] <- SigK[[k]]%*%(crossprod(TpSi, c(colSums(Yi0)) ) + inSig0.m%*%Mu0.m)      ## some of mat
        #cat(dim(SigKWg),"(--3 \n")
      }else{
        SigK[[k]] <- Sig0.m
        SigKWg <- SigKWg + (Thet$Pik.m[k]*Thet$Pik.m[k])*SigK[[k]] # J x J
        Mu.k[k,] <- Mu0.m
      }
    }
    #--- we impose some constrainsts on the mean
    InvSigKWg <- as.positive.definite(as.inverse(SigKWg))
    SIg0 <- bdiag(SigK) # create a block diagonal matrix
    SIGR0 <- NULL;  for(k in 1:Ku){SIGR0 <- rbind(SIGR0,Thet$Pik.m[k]*SigK[[k]])} ## K*J x J
    MUR0 <-  colSums(diag(Thet$Pik.m) %*% Mu.k) # return a evector of length J
    
    # print(Ku)
    # print(dim(SIg0))
    # print(dim(SIGR0))
    #MURFin <- c(t(Mu.k)) - SIGR0%*%solve(SigKWg)%*%MUR0 ## J*K
    #SIGR0fin <- as.matrix(nearPD(SIg0 - SIGR0%*%solve(SigKWg)%*%t(SIGR0),doSym = T)$mat) ## J*K x J*K
    
    MURFin <- c(t(Mu.k)) - SIGR0%*%InvSigKWg%*%MUR0 ## J*K
    SIGR0fin <- as.symmetric.matrix(as.positive.definite(as.matrix(SIg0 - SIGR0%*%InvSigKWg%*%t(SIGR0)))) ## J*K x J*K
    
    #--- Simulate K-1 vectors from the degenerate dist.
    if(Ku > 2){
      id <- 1:(J*(Ku-1))
      SimMU.k <- matrix(c(mvtnorm::rmvnorm(n=1, mean = MURFin[id], SIGR0fin[id,id])), ncol=J,byrow=T) ## (K-1) x J
      SimMU.k1 <- -colSums(diag(Thet$Pik.m[-Ku]) %*% SimMU.k)/Thet$Pik.m[Ku] ## get a K-1 x J matrix follow by colSums - a vector of J
      as.matrix(rbind(SimMU.k,SimMU.k1)) ## k * J
    }else{
      id <- 1:(J*(Ku-1))
      SimMU.k <- matrix(c(mvtnorm::rmvnorm(n=1, mean = MURFin[id], SIGR0fin[id,id])), ncol=J,byrow=T) ## (K-1) x J
      SimMU.k1 <- -(Thet$Pik.m[-Ku] * SimMU.k)/Thet$Pik.m[Ku] ## get a K-1 x J matrix follow by colSums - a vector of J
      as.matrix(rbind(SimMU.k,SimMU.k1)) ## k * J
    }
  }
  updateMuk.m <- cmpfun(updateMuk.m)
  
  #--- Update Sig2k.w based on the mixture prior
  nu0.m = pn + 2; Psi0.m = .5*diag(pn)#as.symmetric.matrix(.25*cov(Wn))  #diag(pn)
  updateSigk.m <- function(Thet, M, nu0.m, Psi0.m, bs2){
    
    Ku = nrow(Thet$Muk.m)
    m = ncol(Thet$X)
    al <- numeric(Ku)
    bet <- numeric(Ku)
    Mtl <- ((M - (tcrossprod(Thet$X, bs2) * matrix(c(Thet$bs2del%*%Thet$Betadel), nrow=nrow(M), ncol = ncol(M), byrow=T))) %*% bs2) 
    
    Res <- matrix(0,ncol = m*m ,nrow = Ku)
    for(i in 1:Ku){
      id <- which(Thet$Cik.m == i)
      nk = length(id)
      if(nk > 0){
        nun <- nu0.m  + nk
        Yi0 <- Mtl[id,]#matrix(c(Mtl[id,]), nrow=nk, byrow=F)
        Psin.m <- as.positive.definite(Psi0.m + crossprod(Yi0))
        Res[i, ] <- c(rinvwishart(nun, Psin.m))
      }else{
        Res[i, ] <- c(rinvwishart(nu0.m, Psi0.m))
      }
    }
    Res
  }
  updateSigk.m <- cmpfun(updateSigk.m)
  
  #-- Function to Compute  hij matrix (N x T)
  LogSimDElVec <- function(M, Mn,Thet,bs2){
    Ku <- length(Thet$Pik.m)
    K = ncol(Thet$X)
    Resf <- 0.0
    for(l in 1:Ku){
      id <- which(Thet$Cik.m == l)
      if(length(id) > 0){
        Sigma <- chol( as.positive.definite(quad.tform(as.symmetric.matrix(matrix(Thet$Sigk.w[l,], ncol = K, nrow = K, byrow=F)), bs2) + 0.0001*diag(nrow(bs2)) ))  #Thet$Sigk.w
        Resf = Resf + sum(LaplacesDemon::dmvnc(M[id,], mu = Mn[id,], U = Sigma, log = T))
      }
      
    }
    Resf
  }
  
  Hij <- function(Betadel, Thet, bs2){
    
    (Thet$X%*% t(bs2)) * matrix(c(Thet$bs2del%*%matrix(Betadel,ncol=1)), nrow = nrow(Thet$X), ncol = nrow(bs2), byrow = T) ## Matrix NxT
  }
  Hij <- cmpfun(Hij)
  
  #'@ Muvec0 : prior mean vector
  #'@ Sig0 : prior variance
  #'@ likelihood + prior
  #--- Second version (without projecting delta(t)*X_{i}(t))
  LikHijMat <- function(Betadel, Thet, M, Muvec0, Sig0, bs2){
    
    MeanV = Hij(Betadel, Thet, bs2)
    
    LogSimDElVec( M, MeanV, Thet, bs2) + mvtnorm::dmvnorm(Betadel, mean = Muvec0, sigma = Sig0, log = T)
  }
  
  
  #--- Updated on 10/16/23 (update Deltat)
  #@' Sigp0 = proposal covariance
  #@' Muvec : prior mean vector
  #@' Sig0 : prior variance
  #@' bs2 : are the basis function
  #--- Without projecting the function to lower space.
  LikHijMatOptim <- function(Thet, M, Muvec, Sig0, bs2){
    #Thet$Betadel
    #MeanV = ((Thet$X%*% t(bs2) * matrix(c(bs2%*%matrix(Betadel,ncol=1)), nrow = nrow(Thet$X), ncol = nrow(bs2), byrow = T))%*%bs2)/nrow(bs2) ## Matrix NxK
    J = length(Muvec)
    parm0 <- c(rnorm(n = J))
    #SigU = chol(as.positive.definite(quad.tform(matrix(Thet$Sigk.m[1,], ncol = J,byrow=F), bs2)))
    
    funC <- function(parm){
      MeanV = Hij(parm, Thet, bs2)
      #print(dim(MeanV))
      #-sum(LogSimDElVec(1, M, MeanV, Thet$X, SigU, Val)) - mvtnorm::dmvnorm(parm, mean = Muvec,sigma = Sig0, log = T)
      -LogSimDElVec(M, MeanV, Thet,bs2) - mvtnorm::dmvnorm(parm, mean = Muvec,sigma = Sig0, log = T)
    }
    #Res <- optim(parm0, funC, method="L-BFGS-B", lower=c(rep(-15,J), 0.01), upper = c(rep(15,J), 10.0), hessian = T)
    Res <- optim(parm0, funC, method="L-BFGS-B", lower=rep(-15,J), upper = rep(15,J), hessian = T)
    #Res <- optim(parm0, funC, method="L-BFGS-B", lower=c(rep(lb,J), 0.01), upper = c(rep(ub,J), 10.0), hessian = T)
    Res
  } 
  
  #--- Updated on 10/16/23 (update Deltat)
  #@' M : matrix of n (samples) x T(number of distinct time points)
  #@' Sigp0 = proposal covariance
  #@' Muvec : prior mean vector
  #@' Sig0 : prior variance
  #@' bs2 : are the basis function
  #mu0.del, Sig0.del
  #Muvec = c(t(bs2)%*%colMeans(M)/colMeans(W))
  #Muvec = mu0.del = rep(0, ncol(Mn)); Sig0 = Sig0.del = diag(rep(.3,ncol(Mn))); Sigp0 = diag(rep(.2,ncol(Mn)))
  #' Find the covariance matrix for the  proposal distribution
  
  Mu0dl = rep(0, ncol(Mn))
  Muvec = mu0.del = rep(0, ncol(Mn)); Sig0 = Sig0.del = diag(rep(sig2.priordelt,ncol(Mn))); Sigp0 = diag(rep(.2,ncol(Mn)))
  Muvec = mu0.del = Thet0$Betadel 
  
  #--- expanded version 
  updateFunDeltaStNew <- function(Thet, M, Mu0dl, Muvec, Sig0, bs2, Sigp0){
    Mitl <- M - Thet$Muk.m[Thet$Cik.m,] %*%t(bs2)
    #-- Current value
    BetadelCur <- Thet$Betadelt  
    #-- Propose a value
    #BetadelNew <- c(mvtnorm::rmvnorm(n=1, mean = BetadelCur, sigma = Sigp0)) 
    BetadelNew <- c(mvtnorm::rmvnorm(n=1, mean = Muvec, sigma = Sigp0)) # simulate a proposal from the prior
    
    #-- Compute the likelihood (proposal)
    ernew <- LikHijMat(BetadelNew, Thet, Mitl, Mu0dl, Sig0, bs2) - mvtnorm::dmvnorm(BetadelNew, mean = Muvec, sigma = Sigp0, log=T)
    
    #-- Compute the likelihood (old)
    erold <- LikHijMat(BetadelCur, Thet, Mitl, Mu0dl, Sig0, bs2) - mvtnorm::dmvnorm(BetadelCur, mean = Muvec, sigma = Sigp0, log = T)
    
    # #-- Compute the likelihood (proposal)
    # ernew <- LikHijMat(BetadelNew, Thet, Mitl, Muvec, Sig0, bs2, Thet$Cik.m) - mvtnorm::dmvnorm(BetadelNew, mean = Muvec, sigma = Sigp0, log=T)
    # #-- Compute the likelihood (old)
    # erold <- LikHijMat(BetadelCur, Thet, Mitl, Muvec, Sig0, bs2, Thet$Cik.m) - mvtnorm::dmvnorm(BetadelCur, mean = Muvec, sigma = Sigp0, log = T)
    
    u0 <- runif(n=1)
    #--- 
    #Fval <- ifelse(log(u0) - (ernew-erold) < 0, ernew, erold)
    accept0 <- 1*(log(u0) < (ernew - erold))
    
    list(res=c(accept0, c(BetadelNew*(accept0) + (1 - accept0)*BetadelCur)), prop = BetadelNew, evid = ernew - erold)
  }
  updateFunDeltaStNew <- cmpfun(updateFunDeltaStNew)
  
  
  #-----------------------------------------------------------------------------
  # Update things related to X 
  #-----------------------------------------------------------------------------
  
  updatePik.x <- function(Thet, alpha.x){
    nk <- numeric(nrow(Thet$Muk.x))
    for(k in 1:nrow(Thet$Muk.x)){ nk[k] = sum(Thet$Cik.x == k)}
    
    #update the Pis
    alpha.xnew <- alpha.x/nrow(Thet$Muk.x) + nk
    
    c(rdirichlet(n=1, alpha.xnew))
  }
  
  
  Mu0.x = numeric(pn); Sig0.x = .5*diag(pn); inSig0.x = 2*diag(pn)#as.inverse(.25*cov(Wn)) #.5*diag(pn) #
  updateMuk.x <- function(Thet, Mu0.x, inSig0.x, Sig0.x){
    K = length(Thet$Pik.x)
    J = ncol(Thet$X)
    Res0 <- matrix(0, nrow = K, ncol = J)
    Mu.k <- matrix(0, nrow = K, ncol = J) # K x J
    #SigKWg <- matrix(0, ncol = ncol(Y),nrow=ncol(Y)) #-- store the weights some of the covarianace - cov of \sum_{k=1}\pi_kMu_k
    Yi.tld = Thet$X 
    Muk <- matrix(0, ncol = J, nrow = K)
    for(k in 1:K){
      idk <- which(Thet$Cik.x == k)
      nk <- length(idk)
      if(length(idk) > 0){
        Yi0 <- matrix(c(Yi.tld[idk,]), nrow=nk, byrow=F)
        TpSi <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.x[k,], ncol=J, byrow=F)))
        SigK <- as.inverse(as.symmetric.matrix((nk*TpSi + inSig0.x))) ## sum of matrices inverse
        Muk <- SigK%*%(crossprod(TpSi, c(colSums(Yi0))) + inSig0.x%*%Mu0.x)      ## some of mat
        Mu.k[k,] <- c(mvtnorm::rmvnorm(n=1, mean = Muk, sigma = SigK))
      }else{
        Mu.k[k,] <- c(mvtnorm::rmvnorm(n=1, mean = Mu0.x, sigma = Sig0.x ))
      }
    }
    #--- we impose some constrainsts on the mean
    Mu.k
    
  }
  updateMuk.x <- cmpfun(updateMuk.x)
  
  #--- Update Sig2k.w based on the mixture of normal prior
  nu0.x <- pn + 2 ; Psi0.x <- .5*diag(pn) #as.symmetric.matrix(.25*cov(Wn))
  updateSigk.x <- function(Thet, nu0.x, Psi0.x){
    
    Ku = nrow(Thet$Muk.x)
    m = ncol(Thet$Muk.x)
    Xtl <- Thet$X - Thet$Muk.x[Thet$Cik.x,]
    Res <- matrix(0,ncol = m*m ,nrow = Ku)
    for(i in 1:Ku){
      id <- which(Thet$Cik.x == i)
      nk = length(id)
      if(nk > 0){
        nux <- nu0.x  + nk
        Yi0 <- matrix(c(Xtl[id,]), nrow=nk, byrow=F)
        Psin.x <- Psi0.x + crossprod(Yi0)
        Res[i, ] <- c(rinvwishart(nux, Psin.x))
      }else{
        Res[i, ] <- c(rinvwishart(nu0.x, Psi0.x))
      }
    }
    Res
  }
  updateSigk.x <- cmpfun(updateSigk.x)
  
  A <- min(c(Wn, Mn)) - .05*diff(range(c(Wn, Mn))) 
  B <- max(c(Wn, Mn)) + .05*diff(range(c(Wn, Mn))) 
  
  
  udapateXi2 <- function(Thet, Y, W, M, A = A, B = B){
    
    J = ncol(Thet$X)
    n = length(Y)
    Ytl = (Y - Thet$muk.e[Thet$Cik.e])/Thet$sig2k.e[Thet$Cik.e]
    Wtl = W - Thet$Muk.w[Thet$Cik.w,]
    Mtl = (M - Thet$Muk.m[Thet$Cik.m,]) #%*%diag(Thet$Betadel)
    Res = matrix(0, ncol=J, nrow = n)
    Sigw = matrix(0, J, J)
    Sigm = matrix(0, J, J)
    Sigx = matrix(0, J, J)
    GamMat= tcrossprod(c(Thet$Gamma))
    
    #for(i in 1:n){
    i = 1
    ie = Thet$Cik.e[i]
    iw = Thet$Cik.w[i]
    im = Thet$Cik.m[i]
    ix = Thet$Cik.x[i]
    
    Sigw <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.w[iw,],ncol=J, byrow=F)))
    Sigm <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.m[im,],ncol=J, byrow=F)))
    Sigx <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.x[ix,],ncol=J, byrow=F)))
    
    Sigxinv <- as.symmetric.matrix(matrix(Thet$Sigk.x[ix,],ncol=J, byrow=F))
    
    Sig <- as.positive.definite(as.inverse(as.symmetric.matrix(GamMat/Thet$sig2k.e[ie] + Sigw +  Sigm + Sigx))) + 0.0001*diag(J)
    
    Mu <- (matrix(c(Ytl), ncol=1) %*% matrix(Thet$Gamma,nrow=1) + tcrossprod(Wtl,Sigw) + tcrossprod(Mtl, Sigm) + 
             tcrossprod(Thet$Muk.x[Thet$Cik.x,],Sigx)) %*% Sig
    
    for(i in 1:n){
      Res[i,] = tryCatch(expr={
        c(TruncatedNormal::rtmvnorm(n=1, mu = Mu[i,], sigma =  Sig,lb = rep(A,J), ub = rep(B,J)))
      },
      warning= function(W){
        if(grepl("Did not find a solution to the nonlinear system in `mvrandn`", W$message)){
          
          c(TruncatedNormal::rtmvnorm(n=1, mu = numeric(J), sigma =  Sig,lb = rep(A,J), ub = rep(B,J)))
          
        }
      }
      )  %>% as.numeric()
    }
    
    Res
    
  }
  udapateXi2 <- cmpfun(udapateXi2)
  
  
  
  udapateXi2 <- function(Thet, Y, W, M, A = A, B = B){
    
    J = ncol(Thet$X)
    n = length(Y)
    #Ytl = (Y - Thet$muk.e[Thet$Cik.e])/Thet$sig2k.e[Thet$Cik.e]
    #Ytl = Y / Thet$sig2k.e[Thet$Cik.e]
    Ytl = (Y - Thet$muk.e[Thet$Cik.e])/Thet$sig2k.e[Thet$Cik.e]
    Wtl = W - Thet$Muk.w[Thet$Cik.w,]
    Mtl = (M - Thet$Muk.m[Thet$Cik.m,]) #%*%diag(Thet$Betadel)
    Res = matrix(0, ncol=J, nrow = n)
    Sigw = matrix(0, J, J)
    Sigm = matrix(0, J, J)
    Sigx = matrix(0, J, J)
    GamMat= tcrossprod(c(Thet$Gamma))
    Mu <- matrix(NA, ncol = ncol(Thet$X), nrow= nrow(Thet$X))
    Mux <- matrix(c(Thet$Muk.x), ncol = ncol(W), nrow=length(Y), byrow=T)
    
    
    nke <- unique(Thet$Cik.e)
    nkx <- unique(Thet$Cik.x)
    nkw <- unique(Thet$Cik.w)
    nkm <- unique(Thet$Cik.m)
    for(kx in 1:length(nkx)){
      Sigx <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.x[kx, ], ncol = J, nrow=J)))
      for(kw in 1:length(nkw)){
        Sigw <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.w[kw, ], ncol = J, nrow=J)))
        for(km in 1:length(nkm)){
          Sigm <- as.inverse(as.symmetric.matrix(matrix(Thet$Sigk.m[km, ], ncol = J, nrow=J)))
          for(ke in 1:length(nke)){
            
            id <- which(as.logical((Thet$Cik.e == nke[ke])*(Thet$Cik.x == nkx[kx])*(Thet$Cik.w == nkw[kw])*(Thet$Cik.m == nkm[km])))
            #id <- which(as.logical((Thet$Cik.e == nke[ke])&(Thet$Cik.x == nkx[kx])&(Thet$Cik.w == nkw[kw])&(Thet$Cik.m == nkm[km])))
            if(length(id) > 0){
              Sig <- as.positive.definite(as.inverse(as.symmetric.matrix(GamMat/Thet$sig2k.e[ke] + Sigw +  Sigm + Sigx))) + 0.0001*diag(J)
              
              Mu[id,] <- (matrix(c(Ytl[id]), ncol=1) %*% matrix(Thet$Gamma,nrow=1) + tcrossprod(Wtl[id,],Sigw) + tcrossprod(Mtl[id,], Sigm) + 
                            tcrossprod(Mux[id,],Sigx)) %*% Sig
              
              for(i in id){
                Res[i,] = tryCatch(expr={
                  c(TruncatedNormal::rtmvnorm(n=1, mu = Mu[i,], sigma =  Sig,lb = rep(A,J), ub = rep(B,J)))
                },
                warning= function(W){
                  if(grepl("Did not find a solution to the nonlinear system in `mvrandn`", W$message)){
                    
                    c(TruncatedNormal::rtmvnorm(n=1, mu = numeric(J), sigma =  Sig,lb = rep(A,J), ub = rep(B,J)))
                    
                  }
                }
                )  %>% as.numeric()
              }
            }
          }
        }
      }
    }
    Res
    
  }
  udapateXi2 <- cmpfun(udapateXi2)
  
  
  
  
  ##----------------------------------------------------
  cat("\n #---- Create containers for the MCMC output ---> \n")
  #----------------------------------------
  #----- Divide these in sections - Construct the output
  #----- containers
  #--------------------------------------------------
  #Nsim <- 5000
  
  MCMCSample <- list()
  
  #-- for X
  MCMCSample$X <- matrix(0,nrow = Nsim,ncol=length(c(Thet$X))) ## c(t()) : how we transform the matrix (row combined) 
  MCMCSample$Muk.x <- matrix(0,nrow = Nsim, ncol = Kx*pn)
  MCMCSample$Sigk.x <- matrix(0,nrow = Nsim,ncol=Kx*pn*pn)
  MCMCSample$Pik.x <- matrix(0,nrow = Nsim,ncol=Kx)
  MCMCSample$Cik.x <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  
  #-- for W
  MCMCSample$Muk.w <- matrix(0,nrow = Nsim,ncol=Kw*pn)
  MCMCSample$Sigk.w <- matrix(0,nrow = Nsim,ncol=Kw*pn*pn)
  MCMCSample$Pik.w <- matrix(0,nrow = Nsim,ncol=Kw)
  MCMCSample$Cik.w <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  
  #-- for M
  MCMCSample$Muk.m <- matrix(0,nrow = Nsim,ncol=Km*pn)
  MCMCSample$Sigk.m <- matrix(0,nrow = Nsim,ncol=Km*pn*pn)
  MCMCSample$Pik.m <- matrix(0,nrow = Nsim,ncol=Km)
  MCMCSample$Cik.m <- matrix(0,nrow = Nsim,ncol=nrow(Wn))
  MCMCSample$Delta <- numeric(length(Thet$Delta))
  MCMCSample$Betdelt <- matrix(0, nrow=Nsim, ncol=length(Thet$Betadelt))
  MCMCSample$AcceptDelt <- numeric(Nsim)
  
  # For Y
  MCMCSample$muk.e <- matrix(0,nrow = Nsim, ncol=length(Thet0$muk.e)) ## c(t()) : how we transform the matrix (row combined)
  MCMCSample$sig2k.e <- matrix(0,nrow = Nsim, ncol=length(Thet0$muk.e))
  MCMCSample$Pik.e <- matrix(0,nrow = Nsim, ncol=Ke)
  MCMCSample$Cik.e <- matrix(0,nrow = Nsim, ncol=nrow(Wn))
  
  MCMCSample$Gamma <- matrix(0,nrow = Nsim, ncol=length(Thet0$Gamma))
  MCMCSample$BetaZ <- matrix(0,nrow = Nsim, ncol=length(Thet$BetaZ))
  MCMCSample$siggam <- numeric(Nsim)
  
  # MCMCSample$AcceptBetaS <- matrix(0,nrow = Nsim, ncol=length(Thet$BetaS))
  # 
  # MCMCSample$BetaX <- matrix(0,nrow = Nsim, ncol=length(Thet$BetaX))
  # 
  # MCMCSample$BetaZ <- matrix(0,nrow = Nsim, ncol=length(c(Thet$BetaZ))) ## c(t()) : how we transform the matrix (row combined)
  # 
  # MCMCSample$BetaV <- matrix(0,nrow = Nsim, ncol=length(Thet$BetaV))
  # MCMCSample$AcceptBetaV <- numeric(Nsim)
  # 
  # MCMCSample$AcceptX <- matrix(0,nrow = Nsim, ncol=nrow(Y))
  # MCMCSample$X <- matrix(0,nrow = Nsim, ncol=nrow(Y))
  #stop()
  #-----------------------------------------------------------------------------
  #------   MCMC sampling 
  #-----------------------------------------------------------------------------
  #Nsim = 1000 
  #library(profvis)
  #profvis(
  #OptParm <- LikHijMatOptim(Thet, Mn, Muvec, Sig0, bs2, Thet$Cik.w)
  OptParm <- LikHijMatOptim(Thet, M, Muvec, Sig0, bs2)
  #Sigp000 <- diag(rep(0.1, length(Thet0$Betadel))) #solve(OptParm$hessian)
  Sigp000 <- diag(rep(sig2.propdelt, length(Thet0$Betadel))) #solve(OptParm$hessian)
  Mu0dl <- Thet0$Betadel
  
  MCMCSample$Propval <- matrix(NA, nrow = Nsim, ncol = length(Thet0$Betadel))
  MCMCSample$evid <- NA*numeric(Nsim)
  
  cat("\n #--- MCMC is starting ----> \n")
  MCMCSample$PsiAcceptF <- numeric(6)
  for(iter in 1:Nsim){
    if(iter %% 500 == 0){cat("iter=",iter,"\n")}
    #----------------------------------------------
    # Update things related to X 
    #--------------------------------------------  
    
    Thet$Pik.x <- updatePik.x(Thet, alpha.x)
    
    Thet$Cik.x <- F2Cmp(CikX(Thet$X, Thet$Muk.x, Thet$Sigk.x, Thet$Pik.x)) 
    
    Thet$Muk.x <- updateMuk.x(Thet, Mu0.x, inSig0.x, Sig0.x)
    
    Thet$Sigk.x <- updateSigk.x(Thet, nu0.x, Psi0.x)
    
    Ytd = Y - Zmod%*%Thet$BetaZ
    
    if( (iter %% lag == 0)*is.null(X)){
      Thet$X <- udapateXi2(Thet, Ytd, Wn, Mn, A=A, B=B)
    }
    
    #---------------------------------------------
    #--- update things related to W
    #---------------------------------------------
    nk <- numeric(Kw)
    for(k in 1:Kw){ nk[k] = sum(Thet$Cik.w == k)}
    alpha.wnew <- alpha.W/Kw + nk
    
    Thet$Pik.w <- rdirch(n=1, alpha.wnew)
    Thet$Cik.w <- F2Cmp(CikW(Wn, Thet$X, Thet$Muk.w, Thet$Sigk.w, Thet$Pik.w))
    
    Thet$Muk.w <- updateMuk.w(Thet, Wn, Mu0.w, inSig0.w, Sig0.w)
    Thet$Sigk.w <- updateSigk.w(Thet, Wn, nu0.w, Psi0.w)
    
    #---------------------------------------------
    #--- update things related to M
    #---------------------------------------------
    nk <- numeric(Km)
    for(k in 1:nrow(Thet$Muk.m)){ nk[k] = sum(Thet$Cik.m == k)}
    alpha.mnew <- alpha.m/Kw + nk
    Thet$Pik.m <- rdirch(n=1, alpha.mnew) 
    
    Xfn = (Thet$X%*% t(bs2) * matrix(c(Thet$bs2del%*%matrix(Thet$Betadelt,ncol=1)), nrow = nrow(Thet$X), ncol = nrow(bs2), byrow = T))
    Thet$Cik.m <- F2Cmp(CikM(M, Xfn, Thet$Muk.m%*%t(bs2), Thet$Sigk.m, Thet$Pik.m, bs2)) 
    
    Thet$Muk.m <- updateMuk.m(Thet, M, Mu0.m, inSig0.m, Sig0.m, bs2)
    Thet$Sigk.m <- updateSigk.m(Thet, M, nu0.m, Psi0.m, bs2)
    
    Out0 <- NA
    Out0 <- updateFunDeltaStNew(Thet, M, Mu0dl, mu0.del, Sig0, bs2, Sigp0 = PsiAccept * Sigp000)
    
    #--- testing the code
    MCMCSample$Propval[iter, ] <- Out0$prop 
    Thet$Betadelt <- Out0$res[-1]
    MCMCSample$Betdelt[iter,] <- Thet$Betadelt
    MCMCSample$AcceptDelt[iter] <- Out0$res[1]
    MCMCSample$evid[iter] <- Out0$evid
    
    Mn <- ((M%*%diag(1/(c(Thet$bs2del%*%Thet$Betadelt))))%*%bs2)/length(a) 
    
    #---------------------------------------------
    #--- update things related to Y or epsilon
    #--------------------------------------------- 
    nk <- numeric(length(Thet$muk.e))
    for(k in 1:length(Thet$muk.e)){ nk[k] = sum(Thet$Cik.e == k)}
    alpha.enew <- alpha.e/length(Thet$muk.e) + nk
    Thet$Pik.e <- c(rdirch(n=1, alpha.enew))
    
    Thet$Cik.e <- F2Cmp(Cike(c(Y), Thet$X, Thet$muk.e, Thet$sig2k.e, Thet$Pik.e,Thet$Gamma))
    Thet$muk.e <- updatemuk.e(Thet, Y, Z, mu0.e, sig20.e)
    Thet$sig2k.e <- updatesig2k.e(Thet,Y, Z, gam0.e, sig20.esig)
    
    Tp0 <- updateGam(Thet, Y, Z, Sig0.gam0, Mu0.gam0)
    Thet$BetaZ <- Tp0[c(1:ncol(Zmod))] 
    Thet$Gamma <- Tp0[-c(1:ncol(Zmod))]
    
    if(is.null(sig02.del)){
      Thet$siggam <- updateSigGam(Thet$Gamma, al0, bet0, order = 2) 
    }
    
    #---------------------------------------------------------------------------
    #--- Save output
    #---------------------------------------------------------------------------
    #-- For X
    MCMCSample$Muk.x[iter, ] <- c(t(Thet$Muk.x))  ## Stor rowwise
    MCMCSample$Sigk.x[iter, ] <- c(t(Thet$Sigk.x))
    MCMCSample$Pik.x[iter, ] <- Thet$Pik.x
    MCMCSample$Cik.x[iter, ] <- Thet$Cik.x
    MCMCSample$X[iter, ] <- c(t(Thet$X)) ## Store row-wise
    
    #-- for W
    MCMCSample$Muk.w[iter, ] <- c(t(Thet$Muk.w))
    MCMCSample$Sigk.w[iter, ] <- c(t(Thet$Sigk.w))
    MCMCSample$Pik.w[iter, ] <- Thet$Pik.w
    MCMCSample$Cik.w[iter, ] <- Thet$Cik.w
    
    #-- for M
    MCMCSample$Muk.m[iter, ] <- c(t(Thet$Muk.m))
    MCMCSample$Sigk.m[iter, ] <- c(t(Thet$Sigk.m))
    MCMCSample$Pik.m[iter, ] <- Thet$Pik.m
    MCMCSample$Cik.m[iter, ] <- Thet$Cik.m
    #-- For Y
    MCMCSample$muk.e[iter, ] <- Thet$muk.e ## c(t()) : how we transform the matrix (row combined)
    MCMCSample$sig2k.e[iter, ] <- c(Thet$sig2k.e)
    MCMCSample$Pik.e[iter, ] <- Thet$Pik.e
    MCMCSample$Cik.e[iter, ] <- Thet$Cik.e
    
    MCMCSample$Gamma[iter, ] <- Thet$Gamma
    MCMCSample$BetaZ[iter, ] <- Thet$BetaZ
    MCMCSample$siggam[iter] <- Thet$siggam
  }  
  
  #---- Return the MCMC smaples
  cat("\n ---- MCMCs are done!! ---/n")
  MCMCSample$Delthat <- colMeans(M)/colMeans(W) #Delthat
  MCMCSample$EstimDelthat <- (bs2del%*%Thet0$Betadel)
  MCMCSample$Thet0 <- Thet0
  MCMCSample$OptimOut <- OptParm
  #MCMCSample$PsiAcceptF <- PsiAccept
  MCMCSample
  #list(Out = MCMCSample, Thet0 = Thet0)
  
}



