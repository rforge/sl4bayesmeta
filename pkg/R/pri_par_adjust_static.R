
####---- function for a static adjustment of prior parameters for HN, EXP, HC, LMX ----####

pri_par_adjust_static <- function(UU=1, tail_prob=0.05){
  # Function for a static adjustment of parameters for HN, EXP, HC, LMX 
  # if tail_prob=0.05 it corresponds to a 5%-tail adjustment suggested by Neuenschwander et al. (2010)
  # input:
  # UU: threshold for P[tau > UU] = tail_prob
  # tail_prob: defined the probability of the tail abbove the threshold UU
  # output:
  # parameters for HN, EXP, HC, LMX
  
  
  # supporting functions
  
  # Analytical formulae to find the scaling factor for a prior distribution fulfilling tail-adjustment 
  # given a threshold UU and an alpha tail-probability
  # They answer the question: Given U and alpha fixed -> how much is AA?
  
  AA_from_Ualpha_HN <- function(UU, alpha){
    return(UU/qnorm(1-alpha/2, mean = 0, sd = 1, lower.tail = TRUE))
  }
  
  AA_from_Ualpha_Exp <- function(UU, alpha){
    return(-UU/log(alpha))
  }
  
  AA_from_Ualpha_HC <- function(UU, alpha){
    return(UU/tan(pi*(1-alpha)/2))
  }
  
  AA_from_Ualpha_Lomax <- function(UU, alpha){
    return(UU*alpha/(1-alpha))
  }
  
  # compuation of scaling parameters
  p_HN <- AA_from_Ualpha_HN(UU=UU, alpha=tail_prob)
  p_EXP <- AA_from_Ualpha_Exp(UU=UU, alpha=tail_prob)
  p_HC <- AA_from_Ualpha_HC(UU=UU, alpha=tail_prob)
  p_LMX <- AA_from_Ualpha_Lomax(UU=UU, alpha=tail_prob)
  
  return(list(p_HN=p_HN, p_EXP=p_EXP, p_HC=p_HC, p_LMX=p_LMX))
  
}
