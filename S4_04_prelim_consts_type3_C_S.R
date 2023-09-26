########################################################
###
### Parameters for model, i.e. knots, consts
###
########################################################
########################################################
###
### Function for determining basis function
### for Constrained GAMs
### Convex shape neither decreasing or increasing
###
########################################################

# Decreasing and Convex
decconvex=function(x,t)
{
  n=length(x)
  k=length(t)-2
  m=k+3
  sigma=matrix(1:(m*n)*0,nrow=m,ncol=n)
  for(j in 1:k){
    i1=x<=t[j]
    sigma[j,i1] = x[i1]-t[1]
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = t[j]-t[1]+((t[j+1]-t[j])^3-(t[j+1]-x[i2])^3)/3/(t[j+1]-t[j])/(t[j+2]-t[j]) +(x[i2]-t[j])*(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = t[j]-t[1] + (t[j+1]-t[j])^2/3/(t[j+2]-t[j]) + (t[j+2]-t[j+1])*(t[j+1]-t[j])/(t[j+2]-t[j]) +((t[j+2]-t[j+1])^3-(t[j+2]-x[i3])^3)/3/(t[j+2]-t[j+1])/(t[j+2]-t[j])
    i4=x>=t[j+2]
    sigma[j,i4] = t[j]-t[1] + (t[j+1]-t[j])^2/3/(t[j+2]-t[j]) + (t[j+2]-t[j+1])*(t[j+1]-t[j])/(t[j+2]-t[j]) +(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[2]
  sigma[k+1,i1]=-(t[2]-x[i1])^3/3/(t[2]-t[1])^2
  i2=x>t[2]
  sigma[k+1,i2]=0
  i1=x<=t[k+1]
  sigma[k+2,i1]=x[i1]-t[1]
  i2=x>t[k+1]&x<=t[k+2]
  sigma[k+2,i2]=t[k+1]-t[1]+((t[k+2]-t[k+1])^2*(x[i2]-t[k+1])-(x[i2]-t[k+1])^3/3)/(t[k+2]-t[k+1])^2
  i3=x>t[k+2]
  sigma[k+2,i3]=t[k+1]-t[1]+((t[k+2]-t[k+1])^2*(t[k+2]-t[k+1])-(t[k+2]-t[k+1])^3/3)/(t[k+2]-t[k+1])^2
  sigma[k+3,]=x
  
  center.vector=apply(sigma,1,mean)
  
  list(sigma=-sigma, center.vector=-center.vector)
}



##############################################################
###
### Basis calculated from the BCGAM Meyer (2008) and bcgam package
###
##############################################################

quant_age <- .2
knots_age <- c(1,round(quantile(right_age-1,c(seq(quant_age,.99,by=quant_age),.99))))
knots_age <- unique(knots_age)
delta.i <- decconvex(1:nT_age,knots_age)
delta <- t(delta.i$sigma-delta.i$center.vector)
Z_age <- delta/max(delta)
nknots_age <- dim(Z_age)[2]

##############################################################
###
### plot of the basis functions
###
##############################################################
pdf("basis_function_age.pdf")
plot(1:nT_age,
     Z_age[, 1],
     type = "l",
     ylim = c(-1, 1),
     main = "Basis Function Age Effect")
for(i in 2:nknots_age){
  lines(1:nT_age,delta[,i])
}
dev.off()

########################################
###
### Spline basis matrix for Period
###
##########################################

intvl_period <- 7
knots_period <- unique(seq(1, nT_period, by = intvl_period))

##Basis for time-varying cause
splinebasis <- bs(1:nT_period, knots = knots_period)

##A constraint matrix so period_effects = 0
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis

##Get a basis for null space of constraint
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc, 
          complete = TRUE)[, (nrow(constr_sumzero) + 1):ncol(constr_sumzero)]
Z_period <- splinebasis %*% Z
nknots_period <- dim(Z_period)[2]

########################################################
###
### Number of MCMC iterations, Chains, Burn-in, Thining
###
########################################################

reps <- 50000
bin <- reps * .5
n_chains <- 3
n_thin <- 1