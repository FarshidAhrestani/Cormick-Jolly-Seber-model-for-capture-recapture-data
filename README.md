# Cormick-Jolly-Seber-for-capture-recapture-data
This JAGS code, initially developed by James Saracco and Andy Royle, can be used to analyze capture-recatpur data from the Monitoring for Avian Population Survivorship (MAPS) program coordinated by the Institute of Bird Populations
 model{

    #### year-specific intercepts for survival (phi), residency (pi),
    #### prob.of predetermining a resident (rho), and recapture prob. (p)
    for(t in 1:nyears-1){
    	p0[t] ~ dunif(0,1) # Priors for mean recapture probability
    	lp0[t] <- log(p0[t]/(1-p0[t])) # Logit transformation of recapture probability 
    	phi0[t] ~ dunif(0,1) # Priors for mean survival probability
    	lphi0[t] <- log(phi0[t]/(1-phi0[t])) # Logit transformation of survival probability
    	rho0[t] ~ dunif(0,1) # Priors for residency probability                       
    	lrho0[t] <- log(rho0[t]/(1-rho0[t]))
    	pi0[t] ~ dunif(0,1) 
    	lpi0[t] <- log(pi0[t]/(1-pi0[t]))
    }
    
    #### random station effects for p model ####
    sigma.p ~ dunif(0, 10)
    tau.p  <- pow(sigma.p, -2)
    for(j in 1:nsta){
    	alpha[j] ~ dnorm(0,tau.p)
    }
    

    for(i in 1:nind){
    	for(t in 1:first[i]){
    		## Define latent state at first capture
    		z[i,t] ~ dbern(1)
    	}
    	logit(rho[i]) <- lrho0[first[i]]
    	R[i] ~ dbern(pi[i,first[i]])
    	mu[i] <- R[i]*rho[i]
    	r[i] ~ dbern(mu[i])
    	for(t in first[i]:nyears-1){
    		logit(p[i,t]) <- lp0[t]  + alpha[sta[i]] 
    		logit(pi[i,t]) <- lpi0[t] 
    		logit(phi[i,t]) <- lphi0[t] 
    	}
    	## State process
    	for (t in first[i]+1:nyears){
    		mu2[i,t] <- z[i,t-1]*phi[i,t-1]*R[i]
    		z[i,t] ~ dbern(mu2[i,t])
    		# Observation process
    		mu1[i,t] <- p[i,t-1]*z[i,t]
    		y[i,t] ~ dbern(mu1[i,t])
    	}
    }
    
    }
