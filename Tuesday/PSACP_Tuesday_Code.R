### At the very start, let's install (if necessary) and load the packages needed

if (system.file(package = "FKSUM") == "") install.packages("FKSUM")
library(FKSUM)

if (system.file(package = "fastICA") == "") install.packages("fastICA")
library(fastICA)

if (system.file(package = "ProDenICA") == "") install.packages("ProDenICA")
library(ProDenICA)

if (system.file(package = "rARPACK") == "") install.packages("rARPACK")
library(rARPACK)




whiten <- function(X, ncomp){
  p <- ncol(X)
  n <- nrow(X)
  E <- eigs_sym(cov(X), ncomp)
  (X-matrix(colMeans(X), n, p, byrow = TRUE))%*%(E$vectors/matrix(sqrt(E$values), nrow = p, ncol = ncomp, byrow = TRUE))
}

ica <- function(X, ncomp, itermax = 200){
  ### we begin by capturing the size of the data matrix
  p <- ncol(X)
  n <- nrow(X)
  
  ### next we whiten the data. Most (all I have looked into)
  ### implementations of ICA operate only on the first ncomp
  ### principal components, and not the complete data set.
  ### This cannot find the optimal solution unless all
  ### PCs are included. To be consistent with existing methods
  ### we'll do that here too.
  Z <- whiten(X, ncomp)
  V <- diag(ncomp)
  
  ### The bandwidth is a simple heuristic from kernel density
  ### estimation.
  h <- .5/n^.2
  
  ### We now define the objective function to be optimised,
  ### and then its gradient.
  f_ica <- function(v, X, h){
    ### the start is as it was with PCA
    p <- X%*%v/sqrt(sum(v^2))
    n <- length(p)
    
    ### the function fk_sum computes the kernel sums of
    ### the form at the top of slide 12. The first argument
    ### is the sample of points whose differences determine
    ### the kernel weights, while the second is the value
    ### for the coefficients ''omega''. The KDE has as
    ### coefficients just the constant value 1/nh. Or
    ### just 1, and then the entire output can be divided
    ### by nh
    -mean(log(fk_sum(p, rep(1, n), h)/n/h))
  }
  df_ica <- function(v, X, h){
    ### again the start is as it was for PCA
    nv <- sqrt(sum(v^2))
    p <- X%*%v/nv
    n <- length(p)
    
    ### fhats stores the KDE values for the projected points
    fhats <- fk_sum(p, rep(1, n), h)/n/h
    
    ### the partial derivatives of the objective w.r.t. the
    ### individual projections, as in the final slide, is
    ### a sum of kernel derivatives with "omegas" equal to
    ### one divided by the density estimates and a sum of
    ### kernel derivatives with "omegas" equal to 1, all divided
    ### by the corresponding density estimate(s). The fk_sum function
    ### can also compute these sums of kernel derivatives by simply
    ### specifying type = 'dksum'.
    dp <- (fk_sum(p, 1/fhats, h, type = 'dksum') + fk_sum(p, rep(1, n), h, type = 'dksum')/fhats)/n^2/h^2
    1/nv*t(X-p%*%t(v)/nv)%*%dp
  }
  
  ### We now, as we did in the PP formulation of PCA, estimate the ICs
  ### sequentially. We cannot use the same deflation scheme as we did
  ### before (why do you think that is?). Instead, we can use a simple
  ### gradient descent approach and project the initialisations and all
  ### gradients into the null space of the IC vectors found so far.
  ### Since the solution will be the sum of the initialisation and a
  ### weighted sum of these gradients, it will too be in this null space.
  for(j in 1:ncomp){
    ### The first IC vector doesn't need to be projected into the
    ### null space of anything, but subsequent ones will
    if(j > 1) V[,j] <- V[,j] - V[,1:(j-1)]%*%t(V[,1:(j-1)])%*%V[,j]
    
    ### We will now use our own implementation of a first order
    ### gradient descent. This is done using backtracking line-search
    ### to determine the step-size for each iteration, with a
    ### minimum step-size of 10^-9.
    for(iter in 1:itermax){
      fval <- f_ica(V[,j], Z, h)
      df <- df_ica(V[,j], Z, h)
      if(j > 1) dir <- df - V[,1:(j-1)]%*%t(V[,1:(j-1)])%*%df
      else dir <- df
      m <- sum(dir*df)
      stp <- 1
      repeat{
        fnew <- f_ica(V[,j] - stp*dir, Z, h)
        if((fval - fnew) > (stp*m/2)){
          V[,j] <- V[,j] - stp*dir
          break
        }
        else stp <- stp/3
        if(stp < 1e-9) break
      }
      if(stp < 1e-9) break
    }
    V[,j] <- V[,j]/sqrt(sum(V[,j]^2))
  }
  Z%*%V
}


### We will look at the ECG example from the slides. We first download the data from a JSS paper
### in which they are used. The original source for the data is cited in that paper.

dataset <- matrix(scan("https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v076i02/foetal_ecg.dat"), 2500, 9, byrow = TRUE)
X <- dataset[ , 2:9]

### We use our method from above, as well as two other existing implementations of ICA.

system.time(ica_us <- ica(X, 8))

system.time(ica_proden <- ProDenICA(X, 8, whiten = TRUE, W0 = diag(8)))

system.time(ica_fast <- fastICA(X, 8))

### We now plot the original data (left) column, and the our method's outputs
### (second column), ProDenICA (third column) and fastICA (final column)

par(mfrow = c(8, 4))
par(mar = c(0, 3, 3, 0))
plot(X[,j], type = 'l', main = 'Original')
plot(ica_us[,j], type = 'l', main = 'Ours')
plot(ica_proden$s[,j], type = 'l', main = 'ProDen')
plot(ica_fast$S[,j], type = 'l', main = 'fast')
for(j in 2:8){
  plot(X[,j], type = 'l')
  plot(ica_us[,j], type = 'l')
  plot(ica_proden$s[,j], type = 'l')
  plot(ica_fast$S[,j], type = 'l')
}

### The plots are difficult to compare by eye, and we have no ground truth to determine
### accuracy. We can try to estimate the Shannon entropy of the components discovered
### by each method. We will use kernel estimates for the density, which arguably biases
### the comparison to favour our own implementation. However, we will use a much more
### sophisticated bandwidth estimation method, based on leave-one-out cross-validation,
### which should make the results more fair.
### Note that the ProDenICA and fastICA fix the maximum likelihood estimates for the
### standard deviations to be one, where we have used the unbiased sample covariance.
### To adjust this, we simple need to divide each of the estimates from the other two
### methods by sqrt(2500/2499) = 1.0002.

sum(sapply(1:8, function(i) sum(log(fk_density(X[,i], h = 'mlcv', x_eval = X[,i])$y))))
sum(sapply(1:8, function(i) sum(log(fk_density(ica_us[,i], h = 'mlcv', x_eval = ica_us[,i])$y))))
sum(sapply(1:8, function(i) sum(log(fk_density(ica_proden$s[,i]/1.0002, h = 'mlcv', x_eval = ica_proden$s[,i]/1.0002)$y))))
sum(sapply(1:8, function(i) sum(log(fk_density(ica_fast$S[,i]/1.0002, h = 'mlcv', x_eval = ica_fast$S[,i]/1.0002)$y))))

### If we are not convinced by using the LOOCV estimates for the bandwidth, we can plot
### the estimated Shannon information for a range of bandwidth values:

bandwidths <- .5/2500^.2*seq(.05, 2, length = 20)

dev.new()

plot(bandwidths, sapply(bandwidths, function(h) sum(sapply(1:8, function(i) sum(log(fk_density(ica_us[,i], h = h, x_eval = ica_us[,i])$y))))), ylim = c(-24472.13, -20000), type = 'l', main = "Ours (solid line), ProDen (dashed line), fast (dotted line)")
lines(bandwidths, sapply(bandwidths, function(h) sum(sapply(1:8, function(i) sum(log(fk_density(ica_proden$s[,i]/1.0002, h = h, x_eval = ica_proden$s[,i]/1.0002)$y))))), lty = 2)
lines(bandwidths, sapply(bandwidths, function(h) sum(sapply(1:8, function(i) sum(log(fk_density(ica_fast$S[,i]/1.0002, h = h, x_eval = ica_fast$S[,i]/1.0002)$y))))), lty = 3)