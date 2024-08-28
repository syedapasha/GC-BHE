#------
# Plots
#------
plot_eps = function(fname) {
	postscript(file=paste(fname,".eps",sep=""), onefile=FALSE, horizontal=FALSE, width=4, height=4, paper="special", family="Times")
	# Trim off excess margin space (bottom, left, top, right)
	par(mar=c(3, 2.9, 0.2, 0.7))
	# Trim off excess outer margin space (bottom, left, top, right)
	par(oma=c(0,0,0,0))
	# Trim off excess space for label and ticks (label, ticks, line)
	par(mgp=c(1.8,0.6,0))
	# lty: line styles (1=solid, 2=dash, 3=dot, 4=dash-dot)
	# lab: (# of x-ticks, # of y-ticks, len of ticks), approximately
	# lwd: line-width
	# cex.lab: fontsize scaling-factor for labels
}
# plot(xxyyzz, xlab="x-label", ylab="y-label", xlim=c(0, 120), ylim=c(0, 50), lty=1:4, lab=c(10, 7, 5), lwd=2, cex.lab=1.3)
# # cex: fontsize scaling-factor for legends
# legend("topright", c("Legend 1", "Legend 2", "Legend 3", "Legend 4"), lty=1:4, lwd=2, cex=1.05)
# dev.off()



#---------------
# Vector QQ Plot
#---------------
plot_vQQ = function(pp) {

	plot_eps("qqJ")
	par(mar = c(3, 3, 1.5, 0.2))
	s <- seq(0.01, 0.99, 0.01)
	len <- length(s)
  
	d <- length(pp)

	qp  <- matrix(0,len,d)
	qpp <- matrix(0,len,d)
	poiss <- list("vector",d) 

	for (k in 1:d) {
		
		# simulate unit rate Poisson process
		poiss[[k]] <- sim_pois(pp[[k]][length(pp[[k]])], 1)
		qp[,k]  <- quantile(poiss[[k]], probs = s)
		qpp[,k] <- quantile(pp[[k]], probs = s) 
	}
	
	plot(qp[,1], qpp[,1], type="l", lwd=1, cex.axis=1.2, cex.lab=1.5, cex.main=1.8, xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="Bivariate Q-Q Plot")

	for (k in 2:d) {
		lines(qp[,k], qpp[,k], type="l", lwd=1, cex.axis=1.2)
	}

	temp <- max(sapply(1:d, function(k) tail(poiss[[k]], 1)))
	lines(0:temp, 0:temp, type="l", lty=2, lwd=1)
	
	dev.off()
}


#---------------
# Scalar QQ Plot
#---------------
plot_sQQ = function(pp) {

	plot_eps("qqM")
	par(mar = c(3, 3, 1.5, 0.8))

	s <- seq(0.01, 0.99, 0.01)
	len <- length(s)

	# simulate unit rate Poisson process
	poiss <- sim_pois(pp[length(pp)], 1)		
	qp <- quantile(poiss, probs = s)
	qpp <- quantile(pp, probs = s) 
	
	plot(qp, qpp, type="l", lwd=1, cex.axis=1.2, cex.lab=1.5, cex.main=1.8, xlab="Theoretical Quantiles", ylab="Sample Quantiles", main="Scalar Q-Q Plot")
	lines(0:tail(poiss,1), 0:tail(poiss,1), type="l", lty=2, lwd=1)  
	
	dev.off()
}


#---------------------------
# QQ Plot of Test Statistics 
#---------------------------
plot_QQ_null = function(DLRT, DMI) {

	plot_eps("qq_null")
	par(mar=c(3, 3, 1.5, 0.8))
	
	s <- seq(0.01, 0.99, 0.01)
	len <- length(s)
	
	qDLRT <- quantile(DLRT, probs = s) 
	qDMI <- quantile(DMI, probs = s)

	M <- max(DLRT[length(DLRT)], DMI[length(DLRT)])

	plot(qDLRT, qDMI, type="l", lwd=1, cex.axis=1.2, cex.lab=1.5, cex.main=1.8, xlab="DLRT Quantiles", ylab="DMI Quantiles", main="Null Q-Q Plot")
	lines(0:M, 0:M, type="l", lty=2, lwd=1)  
	
	dev.off()
}


#---------
# HIR Plot
#---------
plot_hir = function(u, HIR, fname="hir") {
  
	plot_eps(fname)
	par(mar=c(2.8, 1.6, 1.2, 0.2))
  
	if (dim(HIR)[2] > 1) {

		plot(u, HIR[,1], type="l", lty=1, ylim=c(0,max(HIR)), cex.axis=1.2, cex.lab=1.5, cex.main=1.8, xlab="Time (s)", ylab="", main="Bivariate HL HIRs")

		for (k in 2:dim(HIR)[2]) {
			lines(u, HIR[,k], lty=k)	
		}
		
		legend("topright", c(expression(h[aa](t)), expression(h[ab](t)), expression(h[ba](t)), expression(h[bb](t))), lty=1:(dim(HIR)[2]), lwd=1, cex=1.5)
	}
	
	else { 	
	
		plot(u, HIR[,1], type="l", lty=1, ylim=c(0,max(HIR)), cex.axis=1.2, cex.lab=1.5, cex.main=1.8, xlab="Time (s)", ylab="", main="Scalar HL HIR")
	}

	dev.off()
}


#----------
# ROC Curve
#----------
plot_roc = function(alf, S1, S2) {
	
	plot_eps("roc")
	par(mar=c(2.8, 3, 1.5, 0.2))
	
	plot(c(0,1), c(0,1), type="l", lty=3, lwd=1, cex.axis=1.2, cex.lab=1.5, cex.main=1.8, xlab="False Positive Rate", ylab="True Positive Rate", main="ROC")
	lines(alf, S1, type="l", lty=2)
	lines(alf, S2, type="l", lty=1)
	
	legend("bottomright", c("DLRT", "DMI"), lty=c(2,1))
  
	dev.off()
}


#-------------------
# Null Distributions
#-------------------
plot_null = function(dist, ix) {

	h1 <- hist(dist[,1], breaks=50, plot=FALSE)
	h2 <- hist(dist[,2], breaks=50, plot=FALSE)

	if (ix==1) test= "DLRT"		else test= "DMI"
	
	plot_eps(paste0(test,"null"))
	par(mar=c(1.6, 1.6, 1.5, 0.1))

	hist(dist[,ix], breaks=31, xlim=c(min(dist),max(dist)), ylim=c(0,max(h1$counts, h2$counts)), cex.axis=1.2, cex.lab=1.5, cex.main=1.8, xlab="", ylab="", main=test)
	box()
	
	dev.off()
}
