# CBJ April 2018 (Coryn Bailer-Jones, calj@mpia.de)
# Functions for distance inference

library(PolynomF) # required by mode.post3()
library(MASS)     # required by quantile.post3()
source("metropolis.R") # provides metrop(), required by quantile.post3()

# Return 1 if rlo <= r <= rhi, otherwise zero. Vectorized in any one parameter.
lim <- function(r, rlo, rhi) ifelse(r<rlo | r>rhi, 0, 1)

# Normalized Gaussian likelihood in w. Vectorized in any one parameter.
d.like <- function(w, r, wsd) dnorm(x=w, mean=1/r, sd=wsd)

# Unnormalized posterior in distance using 
# exponentially decreasing space density prior (with length scale rlen)
# Vectorized in any one parameter.
ud.distpost3 <- function(r, w, wsd, rlen) {
  d.like(w, r, wsd)*lim(r, 0, Inf)*exp(-r/rlen)*r^2
}

# Return integral of function. Is set to NA if it cannot be calculated (includes case Z<=0)
func.int <- function(func, lower, upper, ...) {
  integ <- integrate(f=func, lower=lower, upper=upper, ..., subdivisions=1e3, rel.tol=1e-5, 
                     abs.tol=0, stop.on.error=FALSE)
  return( ifelse(integ$message=="OK" && integ$value>0 && integ$abs.error>0, integ$value, NA) )
}

# Return posterior mode - a single real number (or NA if something is wrong).
# If retall=TRUE, return all three roots of derivative of PDF (default is FALSE)
# sorted in increasing order of real part
# Inputs:
# w    - parallax,             unrestricted
# wsd  - parallax uncertainty, must be >0
# rlen - prior length scale,   must be >0 and <Inf
# Two types of solutions of the cubic equation: 
# 1. only one root is real => one maxima. Take this as the mode.
# 2. all three roots are real => two maxima
#    if w>=0, take the smallest.
#    if w<0   take the positive one (there should be only one).
# The other two possibilities, 0 or 2 real roots, should not occur.
# If they do, return NA. 
mode.post3 <- function(w, wsd, rlen, retall=FALSE) {
  # special cases:
  # w<=0 works normally, except for w=-Inf
  if(any(is.na(c(w, wsd, rlen)))) return(NA)
  if(w==-Inf)  return(Inf)
  if(w==Inf)   return(0)
  if(wsd==0)   return(1/w)
  if(wsd==Inf) return(2*rlen)
  if(wsd<0)    return(NA) # makes no sense
  if(rlen<=0)  return(NA) # makes no sense
  r <- PolynomF::polynom()
  p <- r^3/rlen - 2*r^2 + (w/wsd^2)*r - 1/wsd^2
  roots <- solve(p) # complex vector of length 3, sorted by increasing size of real part
  rMode <- switch(EXPR = toString(length(which(Im(roots)==0))),
                  "0" = NA,
                  "1" = Re(roots[which(Im(roots)==0)]),
                  "2" = NA,
                  "3" = ifelse(w>0, min(roots), roots[roots>0]) # should be real and unique
  )
  if(retall) {
    return(roots)
  } else {
    return(rMode)
  }
}

# Return the cubic equation discriminant for distpost3
# -ve => one real root (a maximum), two complex roots
#   0 => three identical real roots (a maximum)
# +ve => three distinct real roots
#        select smallest positive root
disc.distpost3 <- function(w, wsd, rlen) {
  return( 36*w/(rlen*wsd^4) - 32/wsd^2 + 4*w^2/wsd^4 - 
            4*w^3/(rlen*wsd^6) - 27/(rlen^2*wsd^4) )
}

# Derivative of ud.distpost3 w.r.t r. Vectorized in any one parameter.
# It calls ud.distpost3 to ensure same scaling constant is used.
grad.ud.distpost3 <- function(r, w, wsd, rlen) {
  ifelse(r<=0, 0, ud.distpost3(r, w, wsd, rlen)*(2/r - 1/rlen - (w-1/r)/(r*wsd)^2))
}

# Second derivative of ud.distpost3 w.r.t r. Vectorized in any one parameter.
# It calls ud.distpost3 to ensure same scaling constant is used.
gradgrad.ud.distpost3 <- function(r, w, wsd, rlen) {
  x <- w-1/r
  ud.distpost3(r, w, wsd, rlen)*
    (1/rlen^2 + 2/r^2 - 4/(r*rlen) - 1/(r^4*wsd^2) + 
       2*x*(1/rlen-1/r)/(r*wsd)^2 + (x/r^2*wsd^2)^2)
}

# Find HDI for distpost3 for given data about specified rMode
# (provider has to ensure this is a mode!)
# Return three element list:
# 1. code:     0=failed due to incorrect input or results which should not occur
#             -1=HDI not found (valid failure)
#             +1=HDI found
# 2. val:     if code=+1 then it is a
#               4-element vector of rMode, lower and upper bounds, 
#               computed probability between them, no. iterations
#             else NA
# 3. message: if code=0 or -1 it is a string error message 
#             else NA
hdi.distpost3 <- function(w, wsd, rlen, rMode, HDIprob, verbose=FALSE) {

  # Trap bogus inputs and return code 0
  if(any(is.na(c(w, wsd, rlen, rMode, HDIprob)))) 
    return(list(code=0, val=NA, message="some inputs NA"))
  if(!is.finite(w)) 
    return(list(code=0, val=NA, message="parallax not finite"))
  if(!(wsd>0 && is.finite(wsd))) 
    return(list(code=0, val=NA, message="parallax uncertainty not (finite and positive)"))
  if(!(rlen>0 && is.finite(rlen))) 
    return(list(code=0, val=NA, message="rlen not (finite and positive)"))
  if(!(rMode>0 && is.finite(rMode))) 
    return(list(code=0, val=NA, message="rMode not (finite and positive)"))
  if(HDIprob<=0 || HDIprob>=1) 
    return(list(code=0, val=NA, message="HDIprob not 0-1"))
  
  # Compute normalization constant - try different upper and lower limits for integration
  rIntHi <- sort(c(1/w + 20*wsd/w^2, 20*rlen)) # candidate upper limits (smaller could be negative)
  rIntLo <- 0 # initial lower limit; must be zero
  sel <- 2 # try largest upper limit
  Z <- func.int(func=ud.distpost3, lower=rIntLo, upper=rIntHi[sel], w, wsd, rlen)
  if(is.na(Z) && rIntHi[1]>0) { # try smaller upper limit, if it's positive
    sel <- 1
    Z <- func.int(func=ud.distpost3, lower=rIntLo, upper=rIntHi[sel], w, wsd, rlen)
      if(is.na(Z)) { # try larger lower limit too, if it's positive and less than upper.
                     # This would be used if parallax is very precise
      rIntLo <- 1/w - 20*wsd/w^2
      if(rIntLo>0 && rIntLo<rIntHi) {
        Z <- func.int(func=ud.distpost3, lower=rIntLo, upper=rIntHi[sel], w, wsd, rlen)
      }
    }
  }
  if(is.na(Z))       return(list(code=-1, val=NA, message="All attempts to find Z via Gaussian quadrature return func.int=NA")) # could happen
  if(is.infinite(Z)) return(list(code=0, val=NA, message="Z computed to be infinite")) # should never happen
  if(Z<=0)           return(list(code=0, val=NA, message="Z computed to be finite but non-positive")) # should never happen
  if(verbose) cat("Z integration limits: ", rIntLo, rIntHi[sel], "Z =", Z, "\n")

  # Iteratively search for HDI starting at rMode
  dP <- 0.01*ud.distpost3(r=rMode, w=w, wsd=wsd, rlen=rlen) # target step size in ud.distpost3
                                                            # (relative to maximum) for search (>0)
  uphill <- FALSE # will be set if search goes uphill (which stops search)
  iMax <- 150     # maximum number of steps for search (>=0)
  rb      <- matrix(nrow=iMax, ncol=2) # lower and upper bounds on r
  udPostb <- matrix(nrow=iMax, ncol=2) # values of unnormalized PDF at rb
  fprob   <- matrix(nrow=iMax, ncol=2) # fraction of probability below and above rMode
  # First step computed using second derivative, as first derivate is zero at mode
  # I don't bother to trap possibility that first step takes us uphill,
  # as this will be handled at next step inside the loop (excluding the
  # very unlikely possibility that next step goes down again).
  dr <- sqrt(-2*dP/gradgrad.ud.distpost3(r=rMode, w=w, wsd=wsd, rlen=rlen))
  if(is.na(dr)) return(list(code=0, val=NA, 
                            message="dr^2 was negative, probably because 2nd derivative was"))
  i <- 1
  rb[i,] <- rMode + c(-1,+1)*dr
  if(rb[i,1]<=0) { # very unlikely, provided dP is sensibly small
    rb[i,1]=0
  }
  udPostb[i,] <- ud.distpost3(r=rb[i,], w=w, wsd=wsd, rlen=rlen)
  fprob[i,]   <- (ud.distpost3(r=rMode, w=w, wsd=wsd, rlen=rlen) + udPostb[i,])*dr/(2*Z)
  if(verbose) cat(i, rb[i,], fprob[i,], sum(fprob[i,]), "\n")
  for(i in 2:iMax) {
    rb[i,] <- rb[i-1,] - dP/grad.ud.distpost3(r=rb[i-1,], w=w, wsd=wsd, rlen=rlen)
    if(rb[i,1]<=0) {
      rb[i,1]=0
    }
    udPostb[i,] <- ud.distpost3(r=rb[i,], w=w, wsd=wsd, rlen=rlen)
    fprob[i,]   <- (udPostb[i,]+udPostb[i-1,])*abs(rb[i,]-rb[i-1,])/(2*Z) + fprob[i-1,]
    if(verbose) cat(i, rb[i,], udPostb[i,], fprob[i,], sum(fprob[i,]), "\n")
    # Trap uphill moves. Either function increases, or gradient changed sign
    # so step was in wrong direction.
    # Also trap a crazy step having led to a non-finite value
    if(!all(is.finite(c(rb[i,], udPostb[i,])))
       || any(udPostb[i,] > udPostb[i-1,]) 
       || (rb[i,1]>rb[i-1,1]) || (rb[i,2]<rb[i-1,2])) { 
      uphill <- TRUE
      break()
    }
    if(sum(fprob[i,])>=HDIprob) break()
  }
  
  # Returns
  if(uphill) {
    return(list(code=-1, val=NA, message="HDI search went uphill"))
  }
  if(i==iMax && sum(fprob[i,])<HDIprob) { # could happen
    return(list(code=-1, val=NA, 
                message="HDIprob not reached within iteration limit"))
  } else {
    return(list(code=+1, val=c(rb[i,], sum(fprob[i,]), i), message=NA))
  }
}

# Define func() required by metrop() for posterior sampling
func.post3 <- function(r, w, wsd, rlen) {
  return( c(log10(ud.distpost3(r, w, wsd, rlen)), 0) )
}

# Find quantiles at specified probs for distpost3 for given data
# using MCMC with specified initialization and step size.
# Return three element list:
# 1. code:     0=failed due to incorrect input or results which should not occur
#             -1=<not used>
#             +1=quantiles found
# 2. val:     if code=+1 then it is a
#               length(prob)-element vector the quantiles
#             else NA
# 3. message: if code=0 or -1 it is a string error message 
#             else NA
quantile.distpost3 <- function(w, wsd, rlen, rInit, rStep, probs, verbose=FALSE) {
  
  # Trap bogus inputs and return code 0
  if(any(is.na(c(w, wsd, rlen, rInit, rStep, probs)))) 
    return(list(code=0, val=NA, message="some inputs NA"))
  if(!is.finite(w)) 
    return(list(code=0, val=NA, message="parallax not finite"))
  if(!(wsd>0 && is.finite(wsd))) 
    return(list(code=0, val=NA, message="parallax uncertainty not (finite and positive)"))
  if(!(rlen>0 && is.finite(rlen))) 
    return(list(code=0, val=NA, message="rlen not (finite and positive)"))
  if(!(rInit>0 && is.finite(rInit)))
    return(list(code=0, val=NA, message="rInit not (finite and positive)"))
  if(!(rStep>0 && is.finite(rStep)))
    return(list(code=0, val=NA, message="rStep not (finite and positive)"))
  if(any(probs<0) || any(probs>1) )
    return(list(code=0, val=NA, message="probs not in range 0-1"))
  
  # metrop must be intialized with positive log density
  if(ud.distpost3(r=rInit, w=w, wsd=wsd, rlen=rlen)<=0) {
    return(list(code=0, val=NA, 
                message="metrop fails as posterior=zero at rInit"))
  } 
  samp <- metrop(func=func.post3, thetaInit=rInit, Nburnin=1e3, Nsamp=2e4,
                 verbose=ifelse(verbose, 1e3, Inf), sampleCov=rStep^2, 
                 w=w, wsd=wsd, rlen=rlen)[,3]

  if(verbose) {
    cat(quantile(samp, probs=probs), "\n")
    MASS::truehist(samp, col="white", xlab="r")
    abline(v=quantile(samp, probs=probs), col="blue")
    Z <- func.int(func=ud.distpost3, lower=0, upper=Inf, w, wsd, rlen)
    if(!is.na(Z)) {
      r <- seq(from=min(samp), to=max(samp), length.out=1e4)
      lines(r, ud.distpost3(r=r, w=w, wsd=wsd, rlen=rlen)/Z)
    }
  }
  
  return(list(code=+1, val=quantile(samp, probs=probs), message=NA))
}

# Plot one PDF specified by w, wsd, rlen (all scalar) in open device.
# If provided, rRes[1] is point estimate and is plotted in blue,
# and rRes[2:3] are CIs and plotted as as dashed lines.
# If modality!=1, plot PDF in red.
plot.ud.distpost3 <- function(w, wsd, rlen, rRes=rep(NA, 3), modality=1, xlim=NA) {
  if(is.na(xlim)) {
    rplotlo <- 0
    rplothi <- ifelse(any(is.na(rRes)), 5*rlen, 2*rRes[3])
  } else {
    rplotlo <- xlim[1]
    rplothi <- xlim[2]
  }
  Nplot   <- 1e3
  s <- seq(from=1/(2*Nplot), by=1/Nplot, length.out=Nplot+1)
  rplot <- s*(rplothi-rplotlo) + rplotlo
  plot(rplot, ud.distpost3(r=rplot, w=w, wsd=wsd, rlen=rlen), 
       col=ifelse(modality==1, "black", "red"),
       xaxs="i", yaxs="i", xlim=c(rplotlo,rplothi), xlab="", ylab="", type="l")
  abline(v=rRes[2:3], lty=2)
  abline(v=rRes[1], col="blue")
}

