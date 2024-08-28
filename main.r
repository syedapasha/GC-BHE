#------------------------------------------------------------
# Granger causality via directed likelihood ratio test (DLRT)
# and directed mutual information (DMI)
#------------------------------------------------------------

rm(list = ls(all = TRUE))
# set working directory, include packages and source files
setwd(getwd())

library(doParallel)
library(foreach)

source("data_gen.r")
source("em.r")
source("plots.r")
source("residual.r")

 
#debug(em_sHL)
ncpu <- c(16)


# simulation setting 
d <- c(2)									# point process dimension
p <- c(3)									# number of Laguerre basis

cc <- rep(.1,d)								# background rates
be <- c(1)									# 1/time constant
rho <- c(.9)								# spectral radius of HIR matrix

T <- c(1000)								# simulation time
del <- c(.1)								# discretization step
s <- seq(0,T,del)							# time steps
len <- length(s)

nR <- c(500)								# number of repeats to construct histogram
Rs <- c(250)								# number of self-simulations
eps <- c(1e-9)								# stopping threshold 


#----------------------
# simulation under null 
#----------------------

# model specification
Wl <- array(0,c(d,d,p))						# array of p non-negative matrices
Wl[1,1,] <- c(1.0, 0.5, 3)
Wl[1,2,] <- c(0.9, 2.0, 3)
Wl[2,2,] <- c(1.0, 0.8, 2)

W <- apply(Wl,1:2,sum)						# non-negative matrix with arbitrary spectral radius
deg <- rowSums(W)
D <- diag(deg)

Al <- array(0,c(d,d,p))						# array of p non-negative matrices satisfying Hawkes stability 
for (i in 1:p) {
	Al[,,i] <- rho*solve(D,Wl[,,i])
}


# #*** uncomment for plots
# load("res.RData")
# source("plots.r")
# source("gen_gc_plots.r")
# 
# browser()
# #*** 




# DLRT statistic
source("test.r")

dim(TS) <- c(d, nR)

DLRT_n <- sort(TS[1,], decreasing = FALSE)


# FPR values
m <- 50										# number of FPR grid points
alf <- c(1:m) / m							# FPR values


r <- nR * (1 - alf) + 1						# index
t <- DLRT_n[r]								# threshold


DMI_n <- sort(TS[2,], decreasing = FALSE)
t_mi <- DMI_n[r]



#-----------------------------
# simulation under alternative 
#-----------------------------

# model specification
Wl[2,1,] <- seq(.05,by=.05,length.out=p)

W <- apply(Wl,1:2,sum)						# non-negative matrix with arbitrary spectral radius
deg <- rowSums(W)
D <- diag(deg)

Al <- array(0,c(d,d,p))						# array of p non-negative matrices satisfying Hawkes stability 
for (i in 1:p) {
	Al[,,i] <- rho*solve(D,Wl[,,i])
}


# DLRT statistic
source("test.r")

dim(TS) <- c(d, nR)

DLRT_a <- sort(TS[1,], decreasing = FALSE)

S <- 1 - sapply(t, function(t_j) which(DLRT_a > t_j, arr.ind = TRUE)[1]) / nR


DMI_a <- sort(TS[2,], decreasing = FALSE)
S_mi <- 1 - sapply(t_mi, function(t_j) which(DMI_a > t_j, arr.ind = TRUE)[1]) / nR


save(DLRT_n, DLRT_a, DMI_n, DMI_a, file="res.RData")


#*** Null distributions ***#
plot_null(cbind(DLRT_n, DMI_n), 1)
plot_null(cbind(DLRT_n, DMI_n), 2)

#*** QQ-plot of Null distributions ***#
plot_QQ_null(DLRT_n, DMI_n) 
  

#*** ROC ***#
plot_roc(alf, S, S_mi)
