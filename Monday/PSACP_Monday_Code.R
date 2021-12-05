set.seed(1)

### At the very start, let's install (if necessary) and load a package
### which contains the data set we will be looking at a bit later.

if (system.file(package = "PPCI") == "") install.packages("PPCI")
library(PPCI)


### We begin by defining two functions which can be used to
### estimate various principal components models based on the
### formulations in the slides.

### PCA_reconstruction takes four arguments:
# X = data matrix with observations row-wise
# ncom = integer number of components desired
# loss = a real valued function on the reals (i.e. loss: R -> R)
#         which determines the loss used to capture
#         reconstruction quality/ penalise the size of
#         residuals.
# dloss = a real valued function on the reals which
#         returns the derivative of the loss function
# V0 = optional initialisation for the optimisation

PCA_reconstruction <- function(X, ncomp, loss, dloss, V0 = NULL){
  ### start with some capturing of parameters
  n <- nrow(X)
  p <- ncol(X)
  
  ### Xc is the (column) centered matrix of 
  ### observations
  Xc <- X - matrix(colMeans(X), n, p, byrow = TRUE)
  
  ### V is the randomly determined initialisation
  ### for the projection matrix
  if(is.null(V0)) V <- matrix(rnorm(p*ncomp), ncol = ncomp)
  else V <- V0
  
  ### obj and dobj are the objective function
  ### and its gradient. R's base function
  ### optim will vectorise the matrix
  ### argument V, and hence it needs to
  ### be turned back into a matrix (Vuse)
  ### within the evaluation of both functions.
  obj <- function(V){
    Vuse <- matrix(V, nrow = ncol(Xc))
    sum(loss(Xc-Xc%*%Vuse%*%t(Vuse)))
  }
  dobj <- function(V){
    Vuse <- matrix(V, nrow = ncol(Xc))
    R <- Xc-Xc%*%Vuse%*%t(Vuse)
    dR <- dloss(R)
    -t(Xc)%*%dR%*%Vuse - t(dR)%*%X%*%Vuse
  }
  
  optim(V, obj, dobj, method = 'BFGS')$par
}

### PCA_proj_pursuit takes the same four arguments

PCA_proj_pursuit <- function(X, ncomp, loss, dloss, V0 = NULL){
  ### The initial stages are the same as in the
  ### previous.
  n <- nrow(X)
  p <- ncol(X)
  Xc <- X - matrix(colMeans(X), n, p, byrow = TRUE)
  if(is.null(V0)) V <- matrix(rnorm(p*ncomp), ncol = ncomp)
  else V <- V0
  
  ### The objective and gradient differ from before.
  ### By default R's optim minimises, hence our
  ### function and gradient are given a negative
  ### sign. Here, because we are using a deflation
  ### scheme, the data matrix will change after each
  ### component is added, hence we give the objective
  ### and gradient the argument Z to represent this
  ### matrix.
  obj <- function(v, Z) -sum(loss(Z%*%v/sqrt(sum(v^2))))
  dobj <- function(v, Z){
    nv <- sqrt(sum(v^2))
    p <- Z%*%v/nv
    dp <- dloss(p)
    -1/nv*t(Z-p%*%t(v)/nv)%*%dp
  }
  ### With the deflation scheme we don't optimise
  ### over all columns in V simultaneously, but
  ### rather sequentially
  for(j in 1:ncomp){
    ### One of the limitations of the way I handle normalisation of the
    ### projection v is that the objective is constant along radii.
    ### This sometimes makes the BFGS algorithm, which is a quasi-Newtom
    ### method converge prematurely. I prefer to run a second instance of
    ### starting from the solution obtained in the first run.
    v <- optim(V[,j]/sqrt(sum(V[,j]^2)), obj, dobj, Xc, method = 'BFGS')$par
    v <- optim(v/sqrt(sum(v^2)), obj, dobj, Xc, method = 'BFGS')$par
    V[,j] <- v/sqrt(sum(v^2))
    ### after determining the j-th column in V, we project the
    ### centered data into its nullspace
    Xc <- Xc-Xc%*%V[,j]%*%t(V[,j])
  }
  V
}

##################################################################################

### Let's now continue with a brief experiment
### We will consider three loss functions: the
### squared loss, the absolute loss, and a non-convex loss function
### with logarithmic growth given by loss(e) = log(1+e^2)
### We begin by plotting these loss functions. In the
### plot they have been scaled to hit the same value
### at x = +- 10 for better visualisation.

xs <- seq(-10, 10, length = 1000)
plot(xs, xs^2, type = 'l', main = 'Loss functions: solid line = squared, dashed line = absolute and dotted line = logarithmic growth')
lines(xs, abs(xs)*100/abs(xs[1]), lty = 2)
lines(xs, log(1+xs^2)*100/log(1+xs[1]^2), lty = 3)

### Next we will estimate principal components using
### these three loss functions, and the two different
### approaches (reconstruction and projection pursuit).
### We will be adding outliers to the data to obscure
### the structure and see how the robustness approaches
### fare.

### First we load the dataset, which is a compressed version
### of the Yale Faces Database B.

data('yale')

### We can begin by splitting the data into a training and
### test set. Although we are not using a supervised
### method, we can assess the performance in the presence of
### outliers by projecting the unseen (uncorrupted) test set
### onto the projection matrices estimated from the corrupted
### training set.

train_ids <- sample(1:nrow(yale$x), ceiling(nrow(yale$x)/2))
Xtr <- yale$x[train_ids,]
Xte <- yale$x[-train_ids,]

### It will also be useful to create function objects which
### represent the loss functions and their derivatives

loss_sq <- function(x) x^2
dloss_sq <- function(x) 2*x

loss_abs <- function(x) abs(x)
dloss_abs <- function(x) sign(x)

loss_log <- function(x) log(1+x^2)
dloss_log <- function(x) 2*x/(1+x^2)

### First we use the data as they are (without added outliers)
### For consistency we will initialise in all cases with the
### same matrix

V0 <- matrix(rnorm(2*ncol(Xtr)), ncol = 2)/sqrt(ncol(Xtr))

system.time(PC_recon_sq <- PCA_reconstruction(Xtr, 2, loss_sq, dloss_sq, V0))

system.time(PC_pp_sq <- PCA_proj_pursuit(Xtr, 2, loss_sq, dloss_sq, V0))

system.time(PC_recon_abs <- PCA_reconstruction(Xtr, 2, loss_abs, dloss_abs, V0))

system.time(PC_pp_abs <- PCA_proj_pursuit(Xtr, 2, loss_abs, dloss_abs, V0))

system.time(PC_recon_log <- PCA_reconstruction(Xtr, 2, loss_log, dloss_log, V0))

system.time(PC_pp_log <- PCA_proj_pursuit(Xtr, 2, loss_log, dloss_log, V0))



### We now plot all 2-dimensional projections obtained
dev.new()
par(mar = c(0, 0, 4, 0))
par(mfrow = c(2, 3))

plot(Xte%*%PC_recon_sq, xaxt = 'n', yaxt = 'n', main = 'Reconstruction. Squared loss')
plot(Xte%*%PC_recon_abs, xaxt = 'n', yaxt = 'n', main = 'Reconstruction. Absolute loss')
plot(Xte%*%PC_recon_log, xaxt = 'n', yaxt = 'n', main = 'Reconstruction. Log loss')
plot(Xte%*%PC_pp_sq, xaxt = 'n', yaxt = 'n', main = 'Projection pursuit. Squared loss')
plot(Xte%*%PC_pp_abs, xaxt = 'n', yaxt = 'n', main = 'Projection pursuit. Absolute loss')
plot(Xte%*%PC_pp_log, xaxt = 'n', yaxt = 'n', main = 'Projection pursuit. Log loss')


dev.new()
par(mar = c(0, 0, 4, 0))
par(mfrow = c(2, 4))

image(matrix(Xte[100,600:1], ncol = 30, byrow = T), main = 'Original', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_recon_sq%*%t(PC_recon_sq))[600:1], ncol = 30, byrow = T), main = 'Reconstruction. Squared loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_recon_abs%*%t(PC_recon_abs))[600:1], ncol = 30, byrow = T), main = 'Reconstruction. Abs loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_recon_log%*%t(PC_recon_log))[600:1], ncol = 30, byrow = T), main = 'Reconstruction. Log loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix(Xte[100,600:1], ncol = 30, byrow = T), main = 'Original', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_pp_sq%*%t(PC_pp_sq))[600:1], ncol = 30, byrow = T), main = 'Projection pursuit. Squared loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_pp_abs%*%t(PC_pp_abs))[600:1], ncol = 30, byrow = T), main = 'Projection pursuit. Abs loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_pp_log%*%t(PC_pp_log))[600:1], ncol = 30, byrow = T), main = 'Projection pursuit. Log loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')


### Now we add outliers by adding Cauchy distributed
### corruptions to 10% of the entries in the training
### data

corr_id <- sample(1:length(Xtr), ceiling(.1*length(Xtr)))

Xtr[corr_id] <- Xtr[corr_id] + rt(length(corr_id), 1)

system.time(PC_recon_sq <- PCA_reconstruction(Xtr, 2, loss_sq, dloss_sq, V0))

system.time(PC_pp_sq <- PCA_proj_pursuit(Xtr, 2, loss_sq, dloss_sq, V0))

system.time(PC_recon_abs <- PCA_reconstruction(Xtr, 2, loss_abs, dloss_abs, V0))

system.time(PC_pp_abs <- PCA_proj_pursuit(Xtr, 2, loss_abs, dloss_abs, V0))

system.time(PC_recon_log <- PCA_reconstruction(Xtr, 2, loss_log, dloss_log, V0))

system.time(PC_pp_log <- PCA_proj_pursuit(Xtr, 2, loss_log, dloss_log, V0))


### We now plot all 2-dimensional projections obtained
dev.new()
par(mar = c(0, 0, 4, 0))
par(mfrow = c(2, 3))

plot(Xte%*%PC_recon_sq, xaxt = 'n', yaxt = 'n', main = 'Reconstruction. Squared loss')
plot(Xte%*%PC_recon_abs, xaxt = 'n', yaxt = 'n', main = 'Reconstruction. Absolute loss')
plot(Xte%*%PC_recon_log, xaxt = 'n', yaxt = 'n', main = 'Reconstruction. Log loss')
plot(Xte%*%PC_pp_sq, xaxt = 'n', yaxt = 'n', main = 'Projection pursuit. Squared loss')
plot(Xte%*%PC_pp_abs, xaxt = 'n', yaxt = 'n', main = 'Projection pursuit. Absolute loss')
plot(Xte%*%PC_pp_log, xaxt = 'n', yaxt = 'n', main = 'Projection pursuit. Log loss')


dev.new()
par(mar = c(0, 0, 4, 0))
par(mfrow = c(2, 4))

image(matrix(Xte[100,600:1], ncol = 30, byrow = T), main = 'Original', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_recon_sq%*%t(PC_recon_sq))[600:1], ncol = 30, byrow = T), main = 'Reconstruction. Squared loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_recon_abs%*%t(PC_recon_abs))[600:1], ncol = 30, byrow = T), main = 'Reconstruction. Abs loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_recon_log%*%t(PC_recon_log))[600:1], ncol = 30, byrow = T), main = 'Reconstruction. Log loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix(Xte[100,600:1], ncol = 30, byrow = T), main = 'Original', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_pp_sq%*%t(PC_pp_sq))[600:1], ncol = 30, byrow = T), main = 'Projection pursuit. Squared loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_pp_abs%*%t(PC_pp_abs))[600:1], ncol = 30, byrow = T), main = 'Projection pursuit. Abs loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')
image(matrix((Xte[100,]%*%PC_pp_log%*%t(PC_pp_log))[600:1], ncol = 30, byrow = T), main = 'Projection pursuit. Log loss', col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')


