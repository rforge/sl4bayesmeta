
####---- epsilon-grid for scaled A*|X| HN, EXP, HC, LMX heterogeneity priors ----####


pri_par_epsilon_grid <- function(AA0_HN, AA0_EXP, AA0_HC, AA0_LMX, grid_epsilon=0.00354){
  # Function for automatic epsilon-grid computation for HN, EXP, HC, LMX applying the A*|X| methodology
  # input:
  # AA0_HN (EXP, HC, LMX): A*|X| scaling of the base heterogeneity prior HN (EXP, HC, LMX) 
  # grid_epsilon: value of the grid_epsilon for epsilon-grid (only two values lower and upper) computation
  # see Roos et al. (2015) for details
  # output:
  # lower and upper parameters for HN, EXP, HC, LMX
  
 
  # supporting functions
  
  
  # HN
  HN_A0_2_Al_Au<-function(AA0, eps=grid_epsilon){
    ## AA0: scaling of the base half normal distribution
    ## eps: epsilon for the epsilon local grid
    ## output AAl, AAu: epsilon local grid for the scaled half normal distribution
    AAl<-AA0*(1/(1-eps^2)^2-sqrt(1/(1-eps^2)^4-1))
    AAu<-AA0*(1/(1-eps^2)^2+sqrt(1/(1-eps^2)^4-1))
    return(c(AAl,AAu))
  }
  
  p_HN_l<-p_HN_l<-HN_A0_2_Al_Au(AA0=AA0_HN)[1] 
  p_HN_u<-p_HN_u<-HN_A0_2_Al_Au(AA0=AA0_HN)[2]
  
  
  # EXP
  EXP_A0_2_Al_Au<-function(AA0, eps=grid_epsilon){
    ## AA0: scaling of the base exponential distribution
    ## eps: epsilon for the epsilon local grid
    ## output AAl, AAu: epsilon local grid for the scaled exponential distribution
    AAl<-2*AA0*(1-(1-eps^2)^2/2-sqrt(1-(1-eps^2)^2))/(1-eps^2)^2
    AAu<-2*AA0*(1-(1-eps^2)^2/2+sqrt(1-(1-eps^2)^2))/(1-eps^2)^2
    return(c(AAl,AAu))
  }
  
  p_EXP_l<-EXP_A0_2_Al_Au(AA0=AA0_EXP)[1] 
  p_EXP_u<-EXP_A0_2_Al_Au(AA0=AA0_EXP)[2]
  
  
  # HC
  obj_HC <- function(x, AA0, eps=grid_epsilon){
    # Function for numerical search for the epsilon local grid for scaled half cauchy distributions
    return(integrate(function(t) {sqrt(dhalfcauchy(t, scale=x)*dhalfcauchy(t, scale=AA0))}, lower = 0, upper = Inf)$value-(1-eps^2))
  }
  
  p_HC_l <- uniroot(obj_HC, lower=0.0001, upper=AA0_HC, tol= 1e-9, AA0=AA0_HC, eps=grid_epsilon)$root
  p_HC_u <- uniroot(obj_HC, lower=AA0_HC, upper=100, tol= 1e-9, AA0=AA0_HC, eps=grid_epsilon)$root
  
  
  # Lomax (LMX)
  obj_LMX <- function(x, AA0, eps=grid_epsilon){
    # Function for numerical search for the epsilon local grid for scaled Lomax distributions
    return(integrate(function(t) {sqrt(dlomax(t, scale=x)*dlomax(t, scale=AA0))}, lower = 0, upper = Inf)$value-(1-eps^2))
  }
  
  p_LMX_l <- uniroot(obj_LMX, lower=0.0001, upper=AA0_LMX, tol= 1e-9, AA0=AA0_LMX, eps=grid_epsilon)$root
  p_LMX_u <- uniroot(obj_LMX, lower=AA0_LMX, upper=100, tol= 1e-9, AA0=AA0_LMX, eps=grid_epsilon)$root
  
  
  return(list(p_HN_l=p_HN_l, p_HN_u=p_HN_u, 
              p_EXP_l=p_EXP_l, p_EXP_u=p_EXP_u, 
              p_HC_l=p_HC_l, p_HC_u=p_HC_u, 
              p_LMX_l=p_LMX_l, p_LMX_u=p_LMX_u))
  
}