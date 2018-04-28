# CBJ April 2018 (Coryn Bailer-Jones, calj@mpia.de)
# Function to produce spherical harmonic basis set

library(pracma)  # provides legendre

# Return spherical harmonic basis with specified lmax for each coordinate (glon,glat)
# Returns a matrix with #rows=length(glon)=length(glat) and #cols=(lmax+1)^2 - 1, 
# i.e. glon,glat projected onto the basis. 
# Note that the l=0 term in the basis (which is just 1) is exclude (hence -1 in the #cols)
# because it can be accommodated by an intercept in a fit.
# Assumes glon=(-180,+180), glat=(-90,+90). Note in degrees, not radian
# In conventional definition of spherical harmonics theta is (0,pi), but here
# I use glat with (-pi/2,+pi/2), I have cos(phi) = sin(glat)
sphharm.basis <- function(glon, glat, lmax) {
  if(!is.finite(lmax) || lmax<1) return(NA)
  if(length(glon)!=length(glat)) return(NA)
  conv <- pi/180
  sphharm <- NULL
  for(l in 1:lmax) {
    cos_mphi <- cos(outer(0:l, conv*glon)) # (l+1) by N matrix of evaluations of cos(m*glon) for m=0..l
    sin_mphi <- sin(outer(1:l, conv*glon)) #     l by N matrix of evaluations of sin(m*glon) for m=1..l
    P_l <- legendre(l, sin(conv*glat))     # (l+1) by N matrix of evaluations of P^m_l(sin glat) for m=0..l
    Ylm_sin <- t(P_l[2:(l+1),] * sin_mphi) # N by (l+1) matrix with elements P^m_l(sin glat)*cos(m*glon) for m=0..l
    Ylm_cos <- t(P_l * cos_mphi)           # N by l     matrix with elements P^m_l(sin glat)*sin(m*glon) for m=1..l
    # cbind(Ylm_sin, Ylm_cos)              # N by 2l+1  matrix of odd and even spherical harmonics of order l
    sphharm <- cbind(sphharm, Ylm_sin, Ylm_cos)
  }
  return(sphharm)
}
