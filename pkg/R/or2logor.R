
or2logor <- function(tab){
  # Function transforming the results on the OR (exp(mu)) scale to results on the logOR (mu) scale
  # Remark: The invariance of the Hellinger distance to any monotone transformation guarantees that S and L values do not change.
  
  # Transform from the OR-scale (exp(mu)) to the log-OR-scale (mu)
  tab_logOR_scale <- tab
  tab_logOR_scale[, "median_post_OR"] <- log(tab_logOR_scale[, "median_post_OR"])
  tab_logOR_scale[, "95CrI_post_OR_low"] <- log(tab_logOR_scale[, "95CrI_post_OR_low"])
  tab_logOR_scale[, "95CrI_post_OR_up"] <- log(tab_logOR_scale[, "95CrI_post_OR_up"])
  tab_logOR_scale[, "length_95CrI_post_OR"] <- tab_logOR_scale[, "95CrI_post_OR_up"]-tab_logOR_scale[, "95CrI_post_OR_low"]
  
  # rename the transformed entries
  colnames(tab_logOR_scale)[5:8] <- c("median_post_mu", "95CrI_post_mu_low", "95CrI_post_mu_up",
                                      "length_95CrI_post_mu") 
  return(tab_logOR_scale)
}



