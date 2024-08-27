#---------------
# Generate plots 
#---------------
source("plots.r")

# simulate a bivariate Hawkes-Laguerre process
pp <- sim_vHL(T, cc, Al, be)		

# fit a bivariate Hawkes-Laguerre model via EM algorithm
parJ <- em_vHL(p, pp, be, T, eps)

# fit a scalar Hawkes-Laguerre model via EM algorithm
parM <- em_sHL(p, pp[[2]], be, T, eps)	


#*** bivariate Hawkes-Laguerre QQ plot ***# 
tau <- residual_vHL(pp, parJ$c, parJ$alf, be)
plot_vQQ(tau)


#*** scalar Hawkes-Laguerre QQ plot ***#
tau <- residual_sHL(pp[[2]], parM$c, parM$alf, be)
plot_sQQ(tau)



#*** HIR for the bivariate Hawkes-Laguerre model ***# 
u <- seq(0, 10, by=.1)
phi <- sapply(1:p, function(el) exp(-be*u) * ((be*u)^(el-1)) * be / gamma(el))
HIR <- sapply(1:d, function(k) sapply(1:d, function(j) phi %*% parJ$alf[k,j,]), simplify = "array")
dim(HIR) <- c(length(u), d^2)

plot_hir(u, HIR)


#*** HIR for the scalar Hawkes-Laguerre model ***# 
u <- seq(0, 10, by=.1)
phi <- sapply(1:p, function(el) exp(-be*u) * ((be*u)^(el-1)) * be / gamma(el))
HIR <- phi %*% parM$alf

plot_hir(u, HIR, "hir1")


#*** Null distributions ***#
plot_null(cbind(DLRT_n, DMI_n), 1)
plot_null(cbind(DLRT_n, DMI_n), 2)

#*** QQ-plot of Null distributions ***#
plot_QQ_null(DLRT_n, DMI_n) 
  

#*** ROC ***#
plot_roc(alf, S, S_mi)
