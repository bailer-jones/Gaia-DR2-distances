# CBJ April 2018 (Coryn Bailer-Jones, calj@mpia.de)
# Compute distance estimators given parallax information

# Input data must be specified before calling this (done in the notebook)

source("functions.R")
source("plot_posterior_prior.R")
source("./length_scale_model/sphharm_basis.R") # provides sphharm.basis

# Setups
set.seed(12345) # only used by metropolis.R
HDIprob <- 0.6827
HDIprobMax <- 0.9 # if computed prob exceeds this, revert to quantiles
wzp <- -0.029e-3  # parallax ZP in arcseconds: True = GDR2 - ZP

# Check inputs
if(!is.finite(w)) {stop("parallax must be finite")}
if(!is.finite(wsd) || wsd<=0) {stop("parallax sd must be finite and positive")}
if(!is.na(rlen)) {
  if(!is.finite(rlen) || rlen<=0) stop("rlen must be finite and positive")
} else {
  if(is.na(glon) || is.na(glat)) {
    stop("either glon,glat or rlen must be specified")
    } else {
      if(glon < 0   || glon > 360) {stop("glon must be in range 0,360")}
      if(glat < -90 || glat > +90) {stop("glat must be in range -90,+90")}
    }
}

# Convert inputs
w   <- 1e-3*w
wsd <- 1e-3*wsd
if(is.na(rlen)) {
  cat("Using Galactic coordinates to compute prior length scale\n\n")
  # To be consistent with model "mod", and for plotting purposes,
  # fold glon from (0,+360) to (-180,+180). glat remains at (-90,+90)
  if(glon>180) {
    glon <- glon-360
  }
  # Determine rlen using fitted Galactic model
  if(!exists("modelisloaded")) { # only load on first run in a session
    load("./length_scale_model/K6_100p_rlenfit_logmedian_sphharm_30_mod_small.Robj") # provides mod
    modelisloaded <- TRUE
  }
  sphharm <- sphharm.basis(glon, glat, lmax=sqrt(mod$rank)-1)
  rlen <- 10^(predict(mod, data.frame(sphharm))) # call *must* have "data.frame(sphharm)" in it
}

# Determine modality and estimate
# If bimodal, set rMode to the one with largest posterior density
if(w>=0 && disc.distpost3(w=w, wsd=wsd, rlen=rlen)>0) {
  modality <- 2
  modes   <- mode.post3(w=w, wsd=wsd, rlen=rlen, retall=TRUE)[c(1,3)]
  rMode   <- modes[which.max(ud.distpost3(r=modes, w, wsd, rlen))]
} else {
  modality <- 1
  rMode <- mode.post3(w=w, wsd=wsd, rlen=rlen, retall=FALSE)
}
  
# Attempt to find HDI about specified mode
failMessage <- NA # set to failure message only if there is a failure 
rHDI <- hdi.distpost3(w=w, wsd=wsd, rlen=rlen, rMode=rMode, HDIprob=HDIprob)
if(rHDI$code==+1) { # HDI found
  if(rHDI$val[3]>HDIprobMax) { # nullify success if HDI search overstep bound too much
    rHDI$code <- -1
  } else {
    rRes <- c(rMode, rHDI$val[1:2], 1, modality)
    HDIinfo <- rHDI$val[3:4]
  }
}
if(rHDI$code==0) { # HDI failed due to pathology: don't find quantiles
  failMessage <- rHDI$message
}
if(rHDI$code==-1) { # valid HDI failure, find quantiles
  rQuant <- quantile.distpost3(w=w, wsd=wsd, rlen=rlen, 
                     rInit=rMode, rStep=rMode/2, 
                     probs=c( 0.5, (1-HDIprob)/2, (1+HDIprob)/2 ), 
                     verbose=TRUE)
  if(rQuant$code==1) {
    rRes <- c(rQuant$val, 2, modality)
  } else {
    failMessage <- rQuant$message
  }
}
if(!is.na(failMessage)) {
  rRes <- c(rep(NA, 3), 0, modality)
}

# Print summary statistics
cat("w[mas] wsd[mas] wsd/w  glon[deg] glat[deg]\n")
cat(1e3*w, 1e3*wsd, wsd/w, " ", glon, glat, "\n\n")
cat("rest[pc] rlo[pc] rhi[pc] rlen[pc] result_flag modality_flag\n")
cat(rRes[1:3], rlen, as.integer(rRes[4]), as.integer(rRes[5]), "\n\n")
if(rRes[4]==1) {
  cat("HDI: probability contained, #steps to find:", HDIinfo, "\n")
}
if(rRes[5]==2) {
  cat("Bimodal: modes at", modes, "\n")
}
if(!is.na(failMessage)) { 
  cat(failMessage, "\n")
}

if(rRes[4]==0) {
  cat("result_flag=0. Attemping to plot...\n")
}

# Plot
if(is.na(rplotlo)) {
  rplotlo <- 0
}
if(is.na(rplothi)) {
  rplothi <- ifelse(any(is.na(rRes)), 5*rlen, 2*rRes[3])
}
plot.post.prior(w=w, wsd=wsd, rlen=rlen, rRes=rRes, rlo, rhi, rplotlo, rplothi)

