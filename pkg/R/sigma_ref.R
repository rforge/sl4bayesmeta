
# computation of the reference standard deviation
sigma_ref<-function(df, type_sigma_ref="geometric"){
  # input:
  # df: data frame with one column "sigma" containing the standard deviations sigmai in each study
  # type_sigma_ref: type of computation of the sigma_ref:
  ## "geometric" geometric mean of within-study sds as suggested in equation (7) by Sorbye and Rue (scaling)
  ## "harmonic" weighted harmonic mean as suggested on p. 490 by Hoaglin (2016) "Misunderstandings about Q..."
  # output:
  # reference standard deviation 
  
  sigma<-df$sigma
  
  if (type_sigma_ref=="geometric"){
     val <- exp(mean(log(sigma)))
  }
  
  if (type_sigma_ref=="harmonic")
  {
    # number of studies
    kk<-length(sigma)
    # weights
    wi<-1/(sigma^2)
    # variance
    ss2<-((kk-1)*sum(wi))/((sum(wi))^2-sum(wi^2))
    # standard deviation
    val <- sqrt(ss2)
  }
  
  return(val)
}