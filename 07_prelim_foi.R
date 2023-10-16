####################################################################################
###
### calculating age_lookup for foi model, 
### which is the ageclass given the number of weeks of age
###
####################################################################################

age_lookup_f <- c(rep(1:4, each = 52),
                       rep(5, 2 * 52),
                       rep(6, 3 * 52))
age_lookup_f <- c(age_lookup_f,
                  rep(7, nT_age - length(age_lookup_f)))

age_lookup_m <- c(rep(1:4, each = 52),
                       rep(5, 2 * 52),
                       rep(6, 3 * 52))
age_lookup_m <- c(age_lookup_m,
                  rep(6, nT_age - length(age_lookup_m)))

age_lookup <- age_lookup_m

# ################################################################
# ###
# ### Function for Basis Expansion - 
# ### Convex shape neither decreasing or increasing
# ### From Meyer et al (2008), and BCGAM R package
# ###
# ################################################################

# #' @keywords internal
# incconvex=function(x,t)
# {
# 	n=length(x)
# 	k=length(t)-2
# 	m=k+3
# 	sigma=matrix(1:m*n,nrow=m,ncol=n)
# 	for(j in 1:(k-1)){
# 		i1=x<=t[j]
# 		sigma[j,i1] = 0
# 	 	i2=x>t[j]&x<=t[j+1]
# 	 	sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
# 	    i3=x>t[j+1]&x<=t[j+2]
# 	    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
# 	    i4=x>t[j+2]
# 	    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
# 	}
# 	i1=x<=t[k]
# 	sigma[k,i1] = 0
# 	i2=x>t[k]&x<=t[k+1]
# 	sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
# 	i3=x>t[k+1]
# 	sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
# 	i1=x<=t[2]
# 	sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3
# 	i2=x>t[2]
# 	sigma[k+1,i2]=x[i2]-t[1]
# 	i1=x<=t[k+1]
# 	sigma[k+2,i1]=0
# 	i2=x>t[k+1]
# 	sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
# 	sigma[k+3,]=x

# 	center.vector=apply(sigma,1,mean)
	
# 	list(sigma=sigma, center.vector=center.vector)
# }


# ##############################################################
# ###
# ### Basis calculated from Meyer (2008) and
# ### bcgam R-package
# ###
# ##############################################################

# knots_foi_cgam <- round(seq(2, n_year, length = 6))
# delta_i <- incconvex(1:n_year, knots_foi_cgam)
# delta <- t(delta_i$sigma - delta_i$center.vector)
# Z_foi_cgam <- delta / max(delta)
# nknots_foi_cgam <- dim(Z_foi_cgam)[2]


# #############################################
# ###
# ### Spline basis matrix for Period Effects
# ###
# #############################################

# knots_foi_spline <- round(seq(2, round(n_year * .5), length = 6))
# splinebasis <- bs(1:(n_year * .5), knots = knots_foi_spline)
# constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis
# qrc <- qr(t(constr_sumzero))
# Z <- qr.Q(qrc,complete=TRUE)[,(nrow(constr_sumzero)+1):ncol(constr_sumzero)]
# Z_foi_spline <- splinebasis%*%Z
# nknots_foi_spline <- dim(Z_foi_spline)[2]

