#------------------------------------------
# simulate a vector Hawkes-Laguerre process
#------------------------------------------
 
sim_vHL = function(T, cc, A, be) {

	dimA <- dim(A)
	d <- dimA[1]
	p <- dimA[3]
	
	pp <- vector("list",d)					# event times of scalar processes

	#*** initialization ***
	lam <- cc								# intensity at time 0
	Is <- sum(lam)

	#*** first event ***
	s <- -log(runif(1)) / Is
	
	#*** attribution test ***
	k <- sum(runif(1) > (cumsum(lam)/Is)) + 1
	pp[[k]] <- c(pp[[k]], s)


	bProceed <- 1
	while (s<T & bProceed==1) {
	
		bProceed <- 0	
		
		#*** general routine ***
		#*** update maximum intensity ***
		
		# compute Riemann-Stieltjes (RS) integral of phi at event time s
		sum_phi <- apply(expand.grid(1:d,1:p), 1, function(x) sum(exp(-be*(s - pp[[x[1]]][pp[[x[1]]]<s])) * ((be*(s - pp[[x[1]]][pp[[x[1]]]<s]))^(x[2]-1)) * be / gamma(x[2])))
		dim(sum_phi) <- c(d,p)
		
		lam <- cc + sapply(1:d, function(x) sum(A[x,,]*sum_phi))

		I <- sum(lam)
		Is <- I + sum(A[,k,1]) * be
		
		while (s<T & bProceed==0) {

			#*** new event ***
			s <- s - log(runif(1)) / Is

			#*** attribution-rejection test ***

			# compute the RS integral of phi at the new event time s
			sum_phi <- apply(expand.grid(1:d,1:p), 1, function(x) sum(exp(-be*(s - pp[[x[1]]][pp[[x[1]]]<s])) * ((be*(s - pp[[x[1]]][pp[[x[1]]]<s]))^(x[2]-1)) * be / gamma(x[2])))
			dim(sum_phi) <- c(d,p)

			lam <- cc + sapply(1:d, function(x) sum(A[x,,]*sum_phi))
			I <- sum(lam)

			D <- runif(1)
			if (s<T & D <= I/Is) {
			
				k <- sum(D > (cumsum(lam)/Is)) + 1
				pp[[k]] <- c(pp[[k]], s)
				
				bProceed <- 1
			}
			else Is = I
		}
	}
	return(pp)
}



#---------------------------
# simulate a Poisson process
#---------------------------

sim_pois = function(T, lam) {
 
	pp <- 0
	len <- 1

	while (pp[len] < T) {
		U <- runif(1)
		pp[len+1] <- pp[len] - log(U)/lam
		len <- len + 1
	}
	
	pp <- pp[2:(length(pp)-1)]
	
	return(pp)	
}
