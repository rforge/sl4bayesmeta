
pri_par_adjust_dynamic <- function(df, rlmc=0.5, tail_prob=0.5){
  # function for a dynamic median RLMC-based scaling parameters adjustment for HN, EXP, HC, LMX
  # input:
  # df: data frame
  # rlmc: target relative latent model complexity 
  # output:
  # parameters for HN, EXP, HC, LMX
  
  # ### Analytical formulae to find the scaling factor for a prior distribution fulfilling tail-adjustment given a threshold UU and an alpha tail-probability
  # ### U and alpha fixed -> how much is AA?
  
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
  
  
  # computation of the reference threshold U_ref
  # P[tau>U_ref]=tail_prob
  U_ref<-sqrt(rlmc/(1-rlmc))*sigma_ref(df)
  
  # tail_prob adjusting
  p_HN <- AA_from_Ualpha_HN(U_ref, alpha=tail_prob)
  p_EXP <- AA_from_Ualpha_Exp(U_ref, alpha=tail_prob)
  p_HC <- AA_from_Ualpha_HC(U_ref, alpha=tail_prob)
  p_LMX <- AA_from_Ualpha_Lomax(U_ref, alpha=tail_prob)
  
  return(list(p_HN=p_HN, p_EXP=p_EXP, p_HC=p_HC, p_LMX=p_LMX))
  
}
