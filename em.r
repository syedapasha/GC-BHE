#----------------------------------------------
# EM algorithm for vector Hawkes-Laguerre model 
#----------------------------------------------

em_vHL = function(p, pp, be, T, eps) {

	M <- sapply(pp, length)					# vector of counts
	d <- length(pp)							# point process dimension
	
	cc <- rep(0,d)							# background rates
	alf <- array(0,c(d,d,p))				# HIR coefficients
	
	L_arr <- vector("list",d)				# log-likelihood (ratio) iterates


	for (k in 1:d) {	

		# compute Riemann-Stieltjes integral of phi at event times
		x <- sapply(1:d, function(j) sapply(1:p, function(el) sapply(pp[[k]], function(tm) sum(exp(-be*(tm - pp[[j]][pp[[j]]<tm])) * (be*(tm - pp[[j]][pp[[j]]<tm]))^(el-1)) * be / gamma(el))), simplify="array")

		xi <- as.vector(rep(1, M[k]))
		for (j in 1:d) {
			xi <- rbind(xi, t(x[,,j]))
		}


		# recursive formula to compute integral of x on interval [0,T]
		# cannot be implemented using sapply
		# B <- matrix(0,p,d)
		# B[1,] <- sapply(1:d, function(j) M[j] - sum(exp(-be*(T - pp[[j]]))))
		# B[2:p,] <- sapply(1:d, function(j) sapply(2:p, function(el) B[el-1,j] - sum(exp(-be*(T - pp[[j]])) * (be*(T - pp[[j]]))^(el-1)) / gamma(el)))

		# faster way to compute integral of x on interval [0,T]
		B <- sapply(1:d, function(j) sapply(1:p, function(el) (-1)^(2*(el-1))*(M[j] - sum(sapply(0:(el-1), function(q) exp(-be*(T - pp[[j]])) * ((be*(T - pp[[j]]))^q) / factorial(q))))))
		
		grad_B <- c(T, c(B))				# second piece of score function


		#*** EM algorithm ***#

		# starting values of parameter vector (c, alf)
		tht <- as.matrix(c(runif(1,0,.1), rep(1/(d*p+1), d*p)))		
		
		mu <- t(tht) %*% xi					# k-th intensity function evaluated at event times of k-th process
		L <- sum(log(mu)) + t(tht) %*% grad_B
		
		err <- 1							# relative error


		while (err > eps) {

			mu <- t(tht) %*% xi				
			grad_A <- sapply(1:(1+d*p), function(i) sum(xi[i,] / mu))		# first piece of score function

			tht_old <- tht
			# multiplicative update
			tht <- tht * grad_A / grad_B
			
			# compute log-likelihood (ratio) iterate
			L_old <- L
			L <- sum(log(mu)) + t(tht) %*% grad_B	
			err <- abs((L - L_old) / L)
			
			L_arr[[k]] <- c(L_arr[[k]], L)
		}

		# parameter estimates
		cc[k] <- tht[1]
		
		temp <- tht[2:(1+d*p)]
		dim(temp) <- c(p,d)
		
		for (j in 1:d) {
			alf[k,j,] <- temp[,j]
		}

	}	
	
	pars <- NULL
	pars$c <- cc
	pars$alf <- alf
	pars$L <- L_arr							# log-likelihood estimates
	return(pars)

}




#----------------------------------------------
# EM algorithm for scalar Hawkes-Laguerre model 
#----------------------------------------------

em_sHL = function(p, pp, be, T, eps) {

	M <- length(pp)

	L_arr <- NULL


	# compute Riemann-Stieltjes integral of phi at event times
	x <- sapply(1:p, function(el) sapply(pp, function(tm) sum(exp(-be*(tm - pp[pp<tm])) * (be*(tm-pp[pp<tm]))^(el-1)) * be / gamma(el)))
	
	xi <- as.vector(rep(1,M))
	xi <- rbind(xi, t(x))


	# recursive formula to compute integral of x on interval [0,T]
	# cannot be implemented using sapply
	# B <- rep(0,p)
	# B[1] <- M - sum(exp(-be*(T - pp)))
	# B[2:p] <- sapply(2:p, function(el) B[el-1] - sum(exp(-be*(T - pp)) * (be*(T - pp))^(el-1)) / gamma(el))

	# faster way to compute integral of x on interval [0,T]
	B <- sapply(1:p, function(el) (-1)^(2*(el-1)) * (M - sum(sapply(0:(el-1), function(q) exp(-be*(T-pp)) * ((be*(T-pp))^q) / factorial(q)))))

	grad_B <- c(T, B)						# second piece of score function


	#*** EM algorithm ***#

	# starting values of parameter vector (c, alf)
	tht <- as.matrix(c(runif(1,0,.1), rep(1/(p+1), p)))		

	mu <- t(tht) %*% xi					# intensity function evaluated at event times 
	L <- sum(log(mu)) + t(tht) %*% grad_B
	
	err <- 1							# relative error


	while (err > eps) {

		mu <- t(tht) %*% xi				
		grad_A <- sapply(1:(1+p), function(i) sum(xi[i,] / mu))		# first piece of score function

		tht_old <- tht
		# multiplicative update
		tht <- tht * grad_A / grad_B
		
		# compute log-likelihood (ratio) iterate
		L_old <- L
		L <- sum(log(mu)) + t(tht) %*% grad_B	
		err <- abs((L - L_old) / L)
		
		L_arr <- c(L_arr, L)
	}
	
	# parameter estimates
	cc <- tht[1]
	alf <- tht[2:(1+p)]

	pars <- NULL
	pars$c <- cc
	pars$alf <- alf
	pars$L <- L_arr							# log-likelihood estimates
	return(pars)

}





#------------------------------------------
# EM algorithm for a vector Poisson process 
#------------------------------------------

em_vPoi = function(pp, T) {

	M <- sapply(pp, length)
	d <- length(pp)

	cc <- sapply(1:d, function(k) M[k] / T)
	L  <- sapply(1:d, function(k) M[k]*log(cc[k]) - cc[k]*T)
	
	pars <- NULL
	pars$c <- cc
	pars$L <- L
	return(pars)
}
