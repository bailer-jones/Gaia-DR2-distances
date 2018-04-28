# CBJ April 2018 (Coryn Bailer-Jones, calj@mpia.de)
# Function to plot posterior and prior

# rlo,rhi are range for normalizing posterior. Prior is analytically normalized.
# rplotlo, rplothi is plotting range.
# Assumes astrometry is in arcseconds and distances in pc
plot.post.prior <- function(w, wsd, rlen, rRes, rlo, rhi, rplotlo, rplothi) {
  Z <- func.int(func=ud.distpost3, lower=rlo, upper=rhi, w, wsd, rlen)
  Nplot   <- 1e3
  s <- seq(from=1/(2*Nplot), by=1/Nplot, length.out=Nplot+1)
  rplot <- s*(rplothi-rplotlo) + rplotlo
  dprior <- (1/(2*rlen^3))*exp(-rplot/rlen)*rplot^2
  dpost  <- ud.distpost3(r=rplot, w=w, wsd=wsd, rlen=rlen)/Z
  plot(1e-3*rplot, dprior, lwd=3, col="green3",
       xaxs="i", yaxs="i", yaxt="n",
       xlim=1e-3*c(rplotlo,rplothi), ylim=c(0,1.05*max(c(dprior,dpost))),
       xlab="distance [kpc]", ylab="", type="l")
  if(rRes[5]==2) {mycol="red"} else {mycol="black"}
  lines(1e-3*rplot, dpost, lwd=3, col=mycol)
  abline(v=1e-3*rRes[1], lwd=2)
  abline(v=1e-3*rRes[2:3], lwd=2, lty=2)
}
