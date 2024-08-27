cl <- makeCluster(ncpu) 
registerDoParallel(cl)

TS <- foreach (R = 1:nR, .combine='c') %dopar% {

	# simulate a bivariate Hawkes-Laguerre process
	pp <- sim_vHL(T, cc, Al, be)		

	# fit a bivariate Hawkes-Laguerre model via EM algorithm
	parJ <- em_vHL(p, pp, be, T, eps)

	# fit a scalar Hawkes-Laguerre model via EM algorithm
	parM <- em_sHL(p, pp[[2]], be, T, eps)	
	
	# compute LR test statistic	
	LRT <- tail(parJ$L[[2]], 1) - tail(parM$L, 1)
	
	
	
	I <- rep(0, Rs)
	
	for (r in 1:Rs) {
	
		# self-simulation
		pp <- sim_vHL(T, parJ$c, parJ$alf, be)		


		x <- sapply(1:d, function(j) sapply(1:p, function(el) sapply(s, function(tm) sum(exp(-be*(tm - pp[[j]][pp[[j]]<tm])) * (be*(tm - pp[[j]][pp[[j]]<tm]))^(el-1)) * be / gamma(el))), simplify="array")

		xi <- as.vector(rep(1, len))
		for (j in 1:d) {
			xi <- rbind(xi, t(x[,,j]))
		}

		lam_b <- c(parJ$c[2], t(parJ$alf[2,,])) %*% xi
		
		
		x <- sapply(1:p, function(el) sapply(s, function(tm) sum(exp(-be*(tm - pp[[2]][pp[[2]]<tm])) * (be*(tm - pp[[2]][pp[[2]]<tm]))^(el-1)) * be / gamma(el)))
		
		xi <- as.vector(rep(1, len))
		xi <- rbind(xi, t(x))

		mu_b <- c(parM$c, parM$alf) %*% xi

		I[r] <- sum(lam_b * log(lam_b / mu_b) - lam_b + mu_b) * del

	}

	DMI <- sum(I) / Rs

	return(c(LRT, DMI))
}

stopCluster(cl)
