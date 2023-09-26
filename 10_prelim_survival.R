

######################################################################################
###
### Preliminary constants for running in the model
###
######################################################################################

### Number of age effects for survival
nT_age_surv <- max(d_surv$right_age_s, na.rm = TRUE) - 1
nT_age_surv_aah_f <- intvl_step_yr_weekly * n_agef + 2
nT_age_short_f <- intvl_step_yr_weekly * (n_agef - 1) + 2
nT_age_surv_aah_m <- intvl_step_yr_weekly * n_agem + 1 
nT_age_short_m <- intvl_step_yr_weekly * (n_agem - 1) + 1
### adding 1 and 2 to account for the fact that deer that live 
### to the 9 or 7 years or older actually live more than
### 52 weeks in a year, due to leap years and timing of weeks beginning into end

########################################################################
### matrix with indexes for averaging period effects from collar data
### to account for within year annual variation when 
### there is no collar data
########################################################################

# interval("2016-05-15","2017-01-09") %/% weeks(1)
indx_mat_pe_surv <- matrix(NA, nrow = 6, ncol = 52)
indx_mat_pe_surv[1, 35:52] <- 1:18
indx_mat_pe_surv[2,] <- 19:(19 + 52 - 1)
indx_mat_pe_surv[3,] <- (19 + 52):(19 + 2 * 52 - 1)
indx_mat_pe_surv[4,] <- (19 + 2 * 52):(19 + 3 * 52 - 1)
indx_mat_pe_surv[5,] <- (19 + 3 * 52):(19 + 4 * 52 - 1)
indx_mat_pe_surv[6,] <- (19 + 4 * 52):(19 + 5 * 52 - 1)

########################################################################
### calibrating age of deer with the study time 
### for indexing in the likelihood for loops
### age2date = left_period - left_age
########################################################################

sus_age2date <- d_fit_sus$left_period_e - d_fit_sus$left_age_e
icap_cens_age2date <- d_fit_icap_cens$left_period_e - d_fit_icap_cens$left_age_e
icap_mort_age2date <- d_fit_icap_mort$left_period_e - d_fit_icap_mort$left_age_e
idead_age2date <- d_fit_idead$left_period_e - d_fit_idead$left_age_e
rec_neg_cens_posttest_age2date <- d_fit_rec_neg_cens_posttest$left_period_e - d_fit_rec_neg_cens_posttest$left_age_e
rec_neg_cens_postno_age2date <- d_fit_rec_neg_cens_postno$left_period_e - d_fit_rec_neg_cens_postno$left_age_e
rec_neg_mort_age2date <- d_fit_rec_neg_mort$left_period_e - d_fit_rec_neg_mort$left_age_e
rec_pos_cens_age2date <- d_fit_rec_pos_cens$left_period_e - d_fit_rec_pos_cens$left_age_e
rec_pos_mort_age2date <- d_fit_rec_pos_mort$left_period_e - d_fit_rec_pos_mort$left_age_e
sus_cens_postno_age2date <- d_fit_sus_cens_postno$left_period_e - d_fit_sus_cens_postno$left_age_e
sus_cens_posttest_age2date <- d_fit_sus_cens_posttest$left_period_e - d_fit_sus_cens_posttest$left_age_e
sus_mort_posttest_age2date <- d_fit_sus_mort_posttest$left_period_e - d_fit_sus_mort_posttest$left_age_e
sus_mort_postno_age2date <- d_fit_sus_mort_postno$left_period_e - d_fit_sus_mort_postno$left_age_e
endlive_age2date <- d_fit_endlive$left_period_e - d_fit_endlive$left_age_e

################################################################
###
### Calculating basis expansions for survival hazards
###
################################################################

####################################
###
### Basis function for Age Effects
### using natural splines
###
####################################

# quant_age <- .05
# knots_age <- unique(c(1,round(quantile(d_fit_sus$right_age_r,c(seq(quant_age,.99, by=quant_age),.99)))))
# nknots_age <- length(knots_age)
# splinebasis <- ns(1:nT_age_surv, knots = knots_age)#,intercept=TRUE,
# constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis
# qrc <- qr(t(constr_sumzero))
# Z <- qr.Q(qrc,complete=TRUE)[,(nrow(constr_sumzero)+1):ncol(constr_sumzero)]
# Z_age <- splinebasis%*%Z
# nknots_age <- dim(Z_age)[2]

#############################################
###
### Spline basis matrix for Period Effects
###
#############################################

intvl_period <- 13
knots_period <- seq(1,nT_period_collar, by = intvl_period)
splinebasis <- bs(1:nT_period_collar, knots = knots_period)
constr_sumzero <- matrix(1, 1, nrow(splinebasis)) %*% splinebasis
qrc <- qr(t(constr_sumzero))
Z <- qr.Q(qrc,complete=TRUE)[,(nrow(constr_sumzero)+1):ncol(constr_sumzero)]
Z_period <- splinebasis%*%Z
nknots_period <- dim(Z_period)[2]


# pdf("figures/basis_function_time_bs.pdf")
# plot(1:nT_period_collar,Z_period[,1],type="l",ylim=c(-1,1),main="Basis Function Time Effect")
# for(i in 2:nknots_period){
#   lines(1:nT_period_collar,Z_period[,i])
# }
# dev.off()


################################################################
###
### Functions for Basis Expansion - 
### From Meyer et al (2008)
###
################################################################

################################################################
### Convex shape neither decreasing or increasing
################################################################

convex = function(x, t, pred.new=TRUE){
  n=length(x)
  k=length(t)-2
  m=k+2
  sigma=matrix(1:m*n,nrow=m,ncol=n)
  for(j in 1:(k-1)){
    i1=x<=t[j]
    sigma[j,i1] = 0
    i2=x>t[j]&x<=t[j+1]
    sigma[j,i2] = (x[i2]-t[j])^3 / (t[j+2]-t[j]) / (t[j+1]-t[j])/3
    i3=x>t[j+1]&x<=t[j+2]
    sigma[j,i3] = x[i3]-t[j+1]-(x[i3]-t[j+2])^3/(t[j+2]-t[j])/(t[j+2]-t[j+1])/3+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
    i4=x>t[j+2]
    sigma[j,i4]=(x[i4]-t[j+1])+(t[j+1]-t[j])^2/3/(t[j+2]-t[j])-(t[j+2]-t[j+1])^2/3/(t[j+2]-t[j])
  }
  i1=x<=t[k]
  sigma[k,i1] = 0
  i2=x>t[k]&x<=t[k+1]
  sigma[k,i2] = (x[i2]-t[k])^3 / (t[k+2]-t[k]) / (t[k+1]-t[k])/3
  i3=x>t[k+1]
  sigma[k,i3] = x[i3]-t[k+1]-(x[i3]-t[k+2])^3/(t[k+2]-t[k])/(t[k+2]-t[k+1])/3+(t[k+1]-t[k])^2/3/(t[k+2]-t[k])-(t[k+2]-t[k+1])^2/3/(t[k+2]-t[k])
  i1=x<=t[2]
  sigma[k+1,i1]=x[i1]-t[1]+(t[2]-x[i1])^3/(t[2]-t[1])^2/3
  i2=x>t[2]
  sigma[k+1,i2]=x[i2]-t[1]
  i1=x<=t[k+1]
  sigma[k+2,i1]=0
  i2=x>t[k+1]
  sigma[k+2,i2]=(x[i2]-t[k+1])^3/(t[k+2]-t[k+1])^2/3
  v1=1:n*0+1
  v2=x
  x.mat=cbind(v1,v2)
  if(pred.new==TRUE){
    list(sigma=sigma,x.mat=x.mat)}
  else{
    if(pred.new==FALSE){
      coef=solve(t(x.mat)%*%x.mat)%*%t(x.mat)%*%t(sigma)
      list(sigma=sigma, x.mat=x.mat, center.vector=coef)}
  }
}


################################################################
### Decreasing Convex shape 
################################################################

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
### Caculating the basis expansion CONVEX
##############################################################

# quant_age <- .2
# knots_age <- c(1, round(quantile(d_surv$right_age_r,
#                        c(seq(quant_age, .99, by = quant_age),
#                        .99))))
# knots_age <- unique(knots_age)
# delta_i <- convex(1:nT_age_surv, knots_age, pred.new = FALSE)
# delta <- t(rbind(delta_i$sigma - 
#                 t(delta_i$x.mat %*%
#                 delta_i$center.vector)))
# delta <- delta / max(delta)
# Z_age <- delta
# nknots_age <- dim(Z_age)[2]

##############################################################
### Caculating the basis expansion CONVEX DECREASING
##############################################################

quant_age <- .2
knots_age <- c(1, round(quantile(d_surv$right_age_r,
                       c(seq(quant_age, .99, by = quant_age),
                       .99))))
knots_age <- unique(knots_age)
delta_i <- decconvex(1:nT_age_surv, knots_age)
delta <- t(delta_i$sigma - delta_i$center.vector)
Z_age <- delta / max(delta)
nknots_age <- dim(Z_age)[2]


# #############################################################
# ###
# ### plot of the basis functions
# ###
# ##############################################################
# pdf("basis_function_age.pdf")
# plot(1:nT_age_surv_aah_f,
#      Z_age[, 1],
#      ylim = c(-1, 1),
#      type = "l",
#      main = "Basis Function Age Effect")
# for (i in 2:nknots_age) {
#   lines(1:nT_age_surv_aah_f, Z_age[, i])
# }
# dev.off()



################################################################
###
### Function for calculating kernel convolution basis expansion
### For period effects
###
#################################################################

# kernel_conv <- nimbleFunction(
#   run = function(nT = double(0),
#                  Z = double(2),
#                  stauk = double(0),
#                  nconst = double(0),
#                  tauk = double(0),
#                  nknots = double(0),
#                  alphau = double(1)
#   ){
#     temp <- nimMatrix(value = 0, nrow = nT, ncol = nknots)
#     temp1 <- nimMatrix(value = 0, nrow = nT, ncol = nknots)
#     temp2 <- nimNumeric(nknots)
#     KA <- nimNumeric(nT)

#     for (i in 1:nT) {
#       for (j in 1:nknots) {
#         temp1[i, j] <- stauk * nconst * exp(-0.5 * Z[i, j]^2 * tauk)
#       }
#     }

#     for (j in 1:nknots) {
#       temp2[j] <- sum(temp1[1:nT, j])
#     }

#     for (i in 1:nT) {
#       for (j in 1:nknots) {
#         temp[i, j] <- (temp1[i, j] / temp2[j]) * alphau[j]
#       }
#       KA[i] <- sum(temp[i, 1:nknots])
#     }
#     muKA <- mean(KA[1:nT])
#     KA[1:nT] <- KA[1:nT] - muKA

#     returnType(double(1))
#     return(KA[1:nT])
#   })

# Ckernel_conv <- compileNimble(kernel_conv)


#########################################
###
### Setting up the distance matrix
### for the kernel convolution on Time
###
#########################################

# intvl_period <- 1
# knots_period <- c(seq(1,
#                       nT_period_collar,
#                       by = intvl_period),
#                 nT_period_collar)
# knots_period <- unique(knots_period)
# nknots_period <- length(knots_period)

# Z_period <- matrix(0, nT_period_collar, nknots_period)
# for (i in 1:nrow(Z_period)) {
#   for (j in 1:nknots_period) {
#     Z_period[i, j] <- abs(i - knots_period[j])
#   }
# }


#############################################################
###
### plot of the basis functions
###
##############################################################
# pdf("figures/basis_function_age.pdf")
# plot(1:nT_age_surv,delta[,1],ylim=c(-1,1),type="l",main="Basis Function Age Effect")
# for(i in 2:nknots_age){
#   lines(1:nT_age_surv,delta[,i])
# }
# dev.off()
