#-------------------------------------------------------
# time transformation for vector Hawkes-Laguerre process
#-------------------------------------------------------

residual_vHL <- function(pp, cc, Al, be) {

	dimA <- dim(Al)
	d <- dimA[1]
	p <- dimA[3]
	M <- sapply(pp, length)
	
	ncpu <- c(d)						# number of cpus for parallel processing



	part_sum <- function(pp, Tm, be, el) {

		res <- 0
		
		T_list <- pp[pp < Tm]
		len <- length(T_list)
		
		if (len) {
		
			dT <- Tm - T_list
			res <- 1 - exp(-be*dT) * .rowSums(sapply(1:el, function(q) ((be * dT)^(q-1)) / gamma(q), simplify="array"), len, el)
		}
		
		return(sum(res))
	}



	cl<-makeCluster(ncpu) 
	registerDoParallel(cl)

	tau <- foreach (k = 1:d) %dopar% {
	
		sapply(pp[[k]], function(Tm) cc[k]*Tm + sum(t(Al[k,,]) * sapply(1:d, function(j) sapply(1:p, function(el) part_sum(pp[[j]], Tm, be, el) ))))
	}

	stopCluster(cl)
	
	return(tau)
}



residual_sHL <- function(pp, cc, Al, be) {

	p <- length(Al)
	M <- length(pp)
	

	tau <- c(cc*pp[1], sapply(2:M, function(m) cc*pp[m] + sum(Al * sapply(1:p, function(el) sum(1 - exp(-be*(pp[m] - pp[1:(m-1)])) * .rowSums(sapply(1:el, function(q) ((be*(pp[m] - pp[1:(m-1)]))^(q-1)) / gamma(q), simplify="array"), m-1, el) )))))

	return(tau)
}
