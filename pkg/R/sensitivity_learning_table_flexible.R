
sensitivity_learning_table_flexible <- function(df, tail_alpha_static=0.04550026, U1=1, U2=2, 
                                               tail_alpha_dynamic=0.5, rlmc1=0.25, rlmc2=0.5, 
                                               type_sigma_ref="geometric", mu_mean=0, mu_sd=4, grid_epsilon=0.00354){
  # This flexible function enables a flexible change of all parameters and provides redundant programming.
  # If necessary, this function can be split into parts without much recoding.
  # It takes ca. 10 min to run this function.
  # learning_sensitivity_table function (multiple calls of bayesmeta)
  # computes a table containing the paremeters of heterogeneity priors HN, EXP, HC, LMX, 
  # descriptive statistics of posteriors,
  # learning and sensitivity estimates for static 5%-tail HN(0.5)(U=1)- and HN(1)(U=2)-like scale parameters for for HN, EXP, HC, LMX heterogeneity priors 
  # learning and sensitivity for dynamically assigned scaling parameters for HN, EXP, HC, LMX for rlmc1 and rlmc2 adjusted at the median
  # input:
  # df: data frame in form compatible with data objects in the bayesmeta environment
  ## tail_alpha_static: P[HN(0.5)>1]=0.04550026, P[HN(1)>2]=0.04550026
  # tail_alpha_static: =0.04550026 (recommended value)
  ## U1: threshold for  P[HN(0.5)>1]=0.04550026
  # U1: =1  (recommended value)
  ## U1: threshold for P[HN(1)>2]=0.04550026
  # U2: =2  (recommended value)
  ## tail_alpha_dynamic=0.5: centering of Uref in the median P[tau>=Uref]=tail_alpha_tynamic=0.5 for a heterogeneity prior for tau
  # tail_alpha_dynamic: =0.5  (recommended value)
  ## rlmc1=0.25: The Relative Latent Model Complexity (RLMC) set to 0.25 (empirically corresponding to HN(0.5))
  #rlmc1: =0.25  (recommended value)
  # rlmc2=0.5: The Relative Latent Model Complexity (RLMC) set to 0.5 (empirically corresponding to HN(1))
  # rlmc2: =0.5  (recommended value)
  ##  Method for the computation of the reference within-study standard deviation: either "geometric" or "harmonic"
  # type_sigma_ref: ="geometric"  (recommended value)
  ## grid_epsilon for sensitivity computation
  ## grid_epsilon=0.00354: grid_epsilon as suggested by Roos et al. (2015) corresponds to H(N(0,1),N(0.01,1))=0.00354
  # grid_epsilon: =0.00354  (recommended value)
  # output:
  # table with the following columns:
  # "U"
  # "tail_prob"
  # "par_val" 
  # "MRLMC" 
  # "median_post_OR"
  # "95CrI_post_OR_low"
  # "95CrI_post_OR_up"
  # "length_95CrI_post_OR" 
  # "median_post_tau"
  # "95CrI_post_tau_low"
  # "95CrI_post_tau_up"
  # "length_95CrI_post_tau"
  # "L_mu"
  # "L_tau"
  # "S_mu"
  # "S_tau"
  # for HN, EXP, HC and LMX heterogeneity priors with the following scaling parameters (indicated in rows)
  # "HN(0.5)_U1tail", "EXP_U1tail", "HC_U1tail", "LMX_U1tail": tail_alpha_static %-tail adjusted with threshold U=U1 (static)
  # "HN_rlmc1", "EXP_rlmc1", "HC_rlmc1", "LMX_rlmc1": dynamically adjusted with target RLMC=rlmc1 and tail probability tail_alpha_dynamic
  # "HN(1)_U2tail", "EXP_U2tail", "HC_U2tail", "LMX_U2tail":: tail_alpha_static %-tail adjusted with threshold U=U2 (static)
  # "HN_rlmc2", "EXP_rlmc2", "HC_rlmc2", "LMX_rlmc2": dynamically adjusted with target RLMC=rlmc2 and tail probability tail_alpha_dynamic
  
  
  # with bayesmeta
  # library(bayesmeta)
  # Version downloaded on 20180928
  
  ####---- parameters for leraning and sensitivity computation (are now provided in the call of the flexible function) ----####
  
  # # mu_mean: mean of the Normal prior for mu
  # mu_mean <- 0
  # # mu_sd: sd of the Normal prior for mu (according to Roever(2018) Bayesmeta-Paper (unit-information prior for log-odds-ratios))
  # mu_sd <- 4
  # ## static choice of the scale parameter for the heterogeneity prior for tau
  # # 1-phalfnormal(q=1, scale=0.5) # P[HN(0.5)>1]=0.04550026
  # # 1-phalfnormal(q=2, scale=1) # P[HN(1)>2]=0.04550026
  # # tail_alpha_static: P[HN(0.5)>1]=0.04550026, P[HN(1)>2]=0.04550026
  # tail_alpha_static <- 0.04550026
  # # U1: thershold for  P[HN(0.5)>1]=0.04550026
  # U1 <- 1
  # # U1: thershold for P[HN(1)>2]=0.04550026
  # U2 <- 2
  # ## dynamic choice of the scale parameter for the heterogeneity prior for tau
  # # tail_alpha_dynamic=0.5: centering of Uref in the median P[tau>=Uref]=tail_alpha_tynamic=0.5 for a heterogeneity prior for tau
  # tail_alpha_dynamic <- 0.5
  # # rlmc1=0.25: The Relative Latent Model Complexity (RLMC) set to 0.25 (empirically corresponding to HN(0.5))
  # rlmc1 <- 0.25
  # # rlmc2=0.5: The Relative Latent Model Complexity (RLMC) set to 0.5 (empirically corresponding to HN(1))
  # rlmc2 <- 0.5
  # ## grid_epsilon for sensitivity computation
  # # grid_epsilon=0.00354: grid_epsilon as suggested by Roos et al. (2015) corresponds to H(N(0,1),N(0.01,1))=0.00354
  # grid_epsilon <- 0.00354
  # 
  
  ####---- epsilon-local grid for mu ----####
  
  
  # normal (for mu)
  # upper and lower epsilon grid values for a normal distribution with mean=0 according to equation (20) in Roos et al. (2015) expressed in terms of mu_sd
  AAu_normal_mu <- mu_sd*(1+sqrt(1-(1-grid_epsilon^2)^4))/(1-grid_epsilon^2)^2
  AAl_normal_mu <- mu_sd*(1-sqrt(1-(1-grid_epsilon^2)^4))/(1-grid_epsilon^2)^2
  
  
  ####---- create an object to collect the results ----####
  
  
  tres<-matrix(NA, nrow=16, ncol=16)
  colnames(tres)<-c("U", "tail_prob", "par_val", "MRLMC", 
                    "median_post_OR", "95CrI_post_OR_low", "95CrI_post_OR_up", "length_95CrI_post_OR", 
                    "median_post_tau", "95CrI_post_tau_low", "95CrI_post_tau_up", "length_95CrI_post_tau",
                    "L_mu", "L_tau", "S_mu", "S_tau")
  rownames(tres)<-c("HN(0.5)_U1tail", "EXP_U1tail", "HC_U1tail", "LMX_U1tail",
                    "HN_rlmc1", "EXP_rlmc1", "HC_rlmc1", "LMX_rlmc1",
                    "HN(1)_U2tail", "EXP_U2tail", "HC_U2tail", "LMX_U2tail",
                    "HN_rlmc2", "EXP_rlmc2", "HC_rlmc2", "LMX_rlmc2")
  
  
  
  ####----  (HN(0.5)-like) static 0.05-tail adjustment with threshold U=1 (HN(0.5)) ----####
  
  # prior parameters 
  
  tres[c(1:4),1]<-U1
  tres[c(1:4),2]<-tail_alpha_static
  tres[1,3]<-pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HN
  tres[2,3]<-pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_EXP
  tres[3,3]<-pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HC
  tres[4,3]<-pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_LMX
  
  # effective mrlmc
  
  tres[1,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HN), 
                                     MM=1000000, output="sample", step=0.03))
  tres[2,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rexp(n=MM, rate=1/pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_EXP), 
                                     MM=1000000, output="sample", step=0.03))
  tres[3,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rhalfcauchy(n=MM, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HC), 
                                     MM=1000000, output="sample", step=0.03))
  tres[4,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rlomax(n=MM, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_LMX, shape=1), 
                                     MM=1000000, output="sample", step=0.03))
  
  
  # posteriors for base HN, EXP, HC and LMX heterogeneity priors
  
  fit.bayesmeta.HN05 <- bayesmeta(y=df[,"y"], 
                                  sigma=df[,"sigma"],
                                  labels=df[,"labels"],
                                  mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                  tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HN)})
  tres[1,c(5:7)]<-exp(fit.bayesmeta.HN05$summary[c("median","95% lower","95% upper"),"mu"])
  tres[1,c(9:11)]<-fit.bayesmeta.HN05$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.EXP_HN05adj <- bayesmeta(y=df[,"y"], 
                                         sigma=df[,"sigma"],
                                         labels=df[,"labels"],
                                         mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                         tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_EXP)})
  tres[2,c(5:7)]<-exp(fit.bayesmeta.EXP_HN05adj$summary[c("median","95% lower","95% upper"),"mu"])
  tres[2,c(9:11)]<-fit.bayesmeta.EXP_HN05adj$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.HC_HN05adj <- bayesmeta(y=df[,"y"], 
                                        sigma=df[,"sigma"],
                                        labels=df[,"labels"],
                                        mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                        tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HC)})
  tres[3,c(5:7)]<-exp(fit.bayesmeta.HC_HN05adj$summary[c("median","95% lower","95% upper"),"mu"])
  tres[3,c(9:11)]<-fit.bayesmeta.HC_HN05adj$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.LMX_HN05adj <- bayesmeta(y=df[,"y"], 
                                         sigma=df[,"sigma"],
                                         labels=df[,"labels"],
                                         mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                         tau.prior=function(t){dlomax(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_LMX, shape=1)})
  tres[4,c(5:7)]<-exp(fit.bayesmeta.LMX_HN05adj$summary[c("median","95% lower","95% upper"),"mu"])
  tres[4,c(9:11)]<-fit.bayesmeta.LMX_HN05adj$summary[c("median","95% lower","95% upper"),"tau"]
  
  
  # learning quantification
  
  # mu
  
  integrand_HN05_mu <- function(x) {sqrt(fit.bayesmeta.HN05$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_EXP_HN05adj_mu <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_HC_HN05adj_mu <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_LMX_HN05adj_mu <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  
  tres[1,13]<-sqrt(1-integrate(integrand_HN05_mu, lower = -Inf, upper = Inf)$value)
  tres[2,13]<-sqrt(1-integrate(integrand_EXP_HN05adj_mu, lower = -Inf, upper = Inf)$value)
  tres[3,13]<-sqrt(1-integrate(integrand_HC_HN05adj_mu, lower = -Inf, upper = Inf)$value) 
  tres[4,13]<-sqrt(1-integrate(integrand_LMX_HN05adj_mu, lower = -Inf, upper = Inf)$value) 
  
  # tau
  
  integrand_HN05_tau <- function(x) {sqrt(fit.bayesmeta.HN05$dposterior(tau=x)*dhalfnormal(x, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HN))}
  integrand_EXP_HN05adj_tau <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj$dposterior(tau=x)*dexp(x, rate=1/pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_EXP))}
  integrand_HC_HN05adj_tau <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj$dposterior(tau=x)*dhalfcauchy(x, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HC))}
  integrand_LMX_HN05adj_tau <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj$dposterior(tau=x)*dlomax(x, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_LMX, shape=1))}
  
  tres[1,14]<-sqrt(1-integrate(integrand_HN05_tau, lower = 0, upper = Inf)$value)
  tres[2,14]<-sqrt(1-integrate(integrand_EXP_HN05adj_tau, lower = 0, upper = Inf)$value)
  tres[3,14]<-sqrt(1-integrate(integrand_HC_HN05adj_tau, lower = 0, upper = Inf)$value) 
  tres[4,14]<-sqrt(1-integrate(integrand_LMX_HN05adj_tau, lower = 0, upper = Inf)$value) 
  
  
  # epsilon-local grid for A*|X| scaled distributions for tau
  
  
  gr2p_HN05adj<-pri_par_epsilon_grid(AA0_HN=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HN, 
                                     AA0_EXP=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_EXP,
                                     AA0_HC=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HC,
                                     AA0_LMX=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_LMX,
                                     grid_epsilon=grid_epsilon)
  
  
  # sensitivity quantification
  
  # HN
  
  fit.bayesmeta.HN05 <- bayesmeta(y=df[,"y"], 
                                  sigma=df[,"sigma"],
                                  labels=df[,"labels"],
                                  mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                  tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HN)})
  
  fit.bayesmeta.HN05_mu_l <- bayesmeta(y=df[,"y"], 
                                       sigma=df[,"sigma"],
                                       labels=df[,"labels"],
                                       mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                       tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HN)})
  
  fit.bayesmeta.HN05_mu_u <- bayesmeta(y=df[,"y"], 
                                       sigma=df[,"sigma"],
                                       labels=df[,"labels"],
                                       mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                       tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HN)})
  
  fit.bayesmeta.HN05_tau_l <- bayesmeta(y=df[,"y"], 
                                        sigma=df[,"sigma"],
                                        labels=df[,"labels"],
                                        mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                        tau.prior=function(t){dhalfnormal(t, scale=gr2p_HN05adj$p_HN_l)})
  
  fit.bayesmeta.HN05_tau_u <- bayesmeta(y=df[,"y"], 
                                        sigma=df[,"sigma"],
                                        labels=df[,"labels"],
                                        mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                        tau.prior=function(t){dhalfnormal(t, scale=gr2p_HN05adj$p_HN_u)})
  
  
  integrand_HN05_l_base_mu <- function(x) {sqrt(fit.bayesmeta.HN05$dposterior(mu=x)*fit.bayesmeta.HN05_mu_l$dposterior(mu=x))}
  sens_HN05_l_base_mu<-sqrt(1-integrate(integrand_HN05_l_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_HN05_u_base_mu <- function(x) {sqrt(fit.bayesmeta.HN05$dposterior(mu=x)*fit.bayesmeta.HN05_mu_u$dposterior(mu=x))}
  sens_HN05_u_base_mu<-sqrt(1-integrate(integrand_HN05_u_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon  
  worst_sens_HN05_mu<-max(sens_HN05_l_base_mu, sens_HN05_u_base_mu)
  
  integrand_HN05_l_base_tau <- function(x) {sqrt(fit.bayesmeta.HN05$dposterior(tau=x)*fit.bayesmeta.HN05_tau_l$dposterior(tau=x))}
  sens_HN05_l_base_tau<-sqrt(1-integrate(integrand_HN05_l_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon  
  integrand_HN05_u_base_tau <- function(x) {sqrt(fit.bayesmeta.HN05$dposterior(tau=x)*fit.bayesmeta.HN05_tau_u$dposterior(tau=x))}
  sens_HN05_u_base_tau<-sqrt(1-integrate(integrand_HN05_u_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon  
  worst_sens_HN05_tau<-max(sens_HN05_l_base_tau, sens_HN05_u_base_tau)
  
  
  # EXP
  
  fit.bayesmeta.EXP_HN05adj <- bayesmeta(y=df[,"y"], 
                                         sigma=df[,"sigma"],
                                         labels=df[,"labels"],
                                         mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                         tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_EXP)})
  
  fit.bayesmeta.EXP_HN05adj_mu_l <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                              tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_EXP)})
  
  fit.bayesmeta.EXP_HN05adj_mu_u <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                              tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_EXP)})
  
  fit.bayesmeta.EXP_HN05adj_tau_l <- bayesmeta(y=df[,"y"], 
                                               sigma=df[,"sigma"],
                                               labels=df[,"labels"],
                                               mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                               tau.prior=function(t){dexp(t, rate=1/gr2p_HN05adj$p_EXP_l)})
  
  fit.bayesmeta.EXP_HN05adj_tau_u <- bayesmeta(y=df[,"y"], 
                                               sigma=df[,"sigma"],
                                               labels=df[,"labels"],
                                               mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                               tau.prior=function(t){dexp(t, rate=1/gr2p_HN05adj$p_EXP_u)})
  
  
  integrand_EXP_HN05adj_l_base_mu <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj$dposterior(mu=x)*fit.bayesmeta.EXP_HN05adj_mu_l$dposterior(mu=x))}
  sens_EXP_HN05adj_l_base_mu<-sqrt(1-integrate(integrand_EXP_HN05adj_l_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_EXP_HN05adj_u_base_mu <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj$dposterior(mu=x)*fit.bayesmeta.EXP_HN05adj_mu_u$dposterior(mu=x))}
  sens_EXP_HN05adj_u_base_mu<-sqrt(1-integrate(integrand_EXP_HN05adj_u_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  worst_sens_EXP_HN05adj_mu<-max(sens_EXP_HN05adj_l_base_mu, sens_EXP_HN05adj_u_base_mu)
  
  integrand_EXP_HN05adj_l_base_tau <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj$dposterior(tau=x)*fit.bayesmeta.EXP_HN05adj_tau_l$dposterior(tau=x))}
  sens_EXP_HN05adj_l_base_tau<-sqrt(1-integrate(integrand_EXP_HN05adj_l_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_EXP_HN05adj_u_base_tau <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj$dposterior(tau=x)*fit.bayesmeta.EXP_HN05adj_tau_u$dposterior(tau=x))}
  sens_EXP_HN05adj_u_base_tau<-sqrt(1-integrate(integrand_EXP_HN05adj_u_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon
  worst_sens_EXP_HN05adj_tau<-max(sens_EXP_HN05adj_l_base_tau, sens_EXP_HN05adj_u_base_tau)
  
  # HC
  
  fit.bayesmeta.HC_HN05adj <- bayesmeta(y=df[,"y"], 
                                        sigma=df[,"sigma"],
                                        labels=df[,"labels"],
                                        mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                        tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HC)})
  
  fit.bayesmeta.HC_HN05adj_mu_l <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                             tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HC)})
  
  fit.bayesmeta.HC_HN05adj_mu_u <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                             tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_HC)})
  
  fit.bayesmeta.HC_HN05adj_tau_l <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                              tau.prior=function(t){dhalfcauchy(t, scale=gr2p_HN05adj$p_HC_l)})
  
  fit.bayesmeta.HC_HN05adj_tau_u <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                              tau.prior=function(t){dhalfcauchy(t, scale=gr2p_HN05adj$p_HC_u)})
  
  
  integrand_HC_HN05adj_l_base_mu <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj$dposterior(mu=x)*fit.bayesmeta.HC_HN05adj_mu_l$dposterior(mu=x))}
  sens_HC_HN05adj_l_base_mu<-sqrt(1-integrate(integrand_HC_HN05adj_l_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_HC_HN05adj_u_base_mu <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj$dposterior(mu=x)*fit.bayesmeta.HC_HN05adj_mu_u$dposterior(mu=x))}
  sens_HC_HN05adj_u_base_mu<-sqrt(1-integrate(integrand_HC_HN05adj_u_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  worst_sens_HC_HN05adj_mu<-max(sens_HC_HN05adj_l_base_mu, sens_HC_HN05adj_u_base_mu)
  
  integrand_HC_HN05adj_l_base_tau <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj$dposterior(tau=x)*fit.bayesmeta.HC_HN05adj_tau_l$dposterior(tau=x))}
  sens_HC_HN05adj_l_base_tau<-sqrt(1-integrate(integrand_HC_HN05adj_l_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_HC_HN05adj_u_base_tau <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj$dposterior(tau=x)*fit.bayesmeta.HC_HN05adj_tau_u$dposterior(tau=x))}
  sens_HC_HN05adj_u_base_tau<-sqrt(1-integrate(integrand_HC_HN05adj_u_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  worst_sens_HC_HN05adj_tau<-max(sens_HC_HN05adj_l_base_tau, sens_HC_HN05adj_u_base_tau)
  
  # LMX
  
  fit.bayesmeta.LMX_HN05adj <- bayesmeta(y=df[,"y"], 
                                         sigma=df[,"sigma"],
                                         labels=df[,"labels"],
                                         mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                         tau.prior=function(t){dlomax(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN05adj_mu_l <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                              tau.prior=function(t){dlomax(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN05adj_mu_u <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                              tau.prior=function(t){dlomax(t, scale=pri_par_adjust_static(UU=U1, tail_prob=tail_alpha_static)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN05adj_tau_l <- bayesmeta(y=df[,"y"], 
                                               sigma=df[,"sigma"],
                                               labels=df[,"labels"],
                                               mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                               tau.prior=function(t){dlomax(t, scale=gr2p_HN05adj$p_LMX_l, shape=1)})
  
  fit.bayesmeta.LMX_HN05adj_tau_u <- bayesmeta(y=df[,"y"], 
                                               sigma=df[,"sigma"],
                                               labels=df[,"labels"],
                                               mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                               tau.prior=function(t){dlomax(t, scale=gr2p_HN05adj$p_LMX_u, shape=1)})
  
  
  integrand_LMX_HN05adj_l_base_mu <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj$dposterior(mu=x)*fit.bayesmeta.LMX_HN05adj_mu_l$dposterior(mu=x))}
  sens_LMX_HN05adj_l_base_mu<-sqrt(1-integrate(integrand_LMX_HN05adj_l_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_LMX_HN05adj_u_base_mu <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj$dposterior(mu=x)*fit.bayesmeta.LMX_HN05adj_mu_u$dposterior(mu=x))}
  sens_LMX_HN05adj_u_base_mu<-sqrt(1-integrate(integrand_LMX_HN05adj_u_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon
  worst_sens_LMX_HN05adj_mu<-max(sens_LMX_HN05adj_l_base_mu, sens_LMX_HN05adj_u_base_mu)
  
  integrand_LMX_HN05adj_l_base_tau <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj$dposterior(tau=x)*fit.bayesmeta.LMX_HN05adj_tau_l$dposterior(tau=x))}
  sens_LMX_HN05adj_l_base_tau<-sqrt(1-integrate(integrand_LMX_HN05adj_l_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_LMX_HN05adj_u_base_tau <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj$dposterior(tau=x)*fit.bayesmeta.LMX_HN05adj_tau_u$dposterior(tau=x))}
  sens_LMX_HN05adj_u_base_tau<-sqrt(1-integrate(integrand_LMX_HN05adj_u_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  worst_sens_LMX_HN05adj_tau<-max(sens_LMX_HN05adj_l_base_tau, sens_LMX_HN05adj_u_base_tau)
  
  
  
  tres[1,15]<-worst_sens_HN05_mu
  tres[2,15]<-worst_sens_EXP_HN05adj_mu
  tres[3,15]<-worst_sens_HC_HN05adj_mu 
  tres[4,15]<-worst_sens_LMX_HN05adj_mu 
  
  tres[1,16]<-worst_sens_HN05_tau
  tres[2,16]<-worst_sens_EXP_HN05adj_tau
  tres[3,16]<-worst_sens_HC_HN05adj_tau  
  tres[4,16]<-worst_sens_LMX_HN05adj_tau
  
  
  
  
  
  
  ####----  (HN(1)-like) static 0.05-tail adjustment with threshold U=2 (HN(1)) ----####
  
  # prior parameters 
  
  tres[c(9:12),1]<-U2
  tres[c(9:12),2]<-tail_alpha_static
  tres[9,3]<-pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HN
  tres[10,3]<-pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_EXP
  tres[11,3]<-pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HC
  tres[12,3]<-pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_LMX
  
  
  # effective mrlmc
  
  
  tres[9,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HN), 
                                     MM=1000000, output="sample", step=0.03))
  tres[10,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rexp(n=MM, rate=1/pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_EXP), 
                                      MM=1000000, output="sample", step=0.03))
  tres[11,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rhalfcauchy(n=MM, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HC), 
                                      MM=1000000, output="sample", step=0.03))
  tres[12,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rlomax(n=MM, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_LMX, shape=1), 
                                      MM=1000000, output="sample", step=0.03))
  
  # posteriors for base HN, EXP, HC and LMX heterogeneity priors
  
  
  fit.bayesmeta.HN1 <- bayesmeta(y=df[,"y"], 
                                 sigma=df[,"sigma"],
                                 labels=df[,"labels"],
                                 mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                 tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HN)})
  tres[9,c(5:7)]<-exp(fit.bayesmeta.HN1$summary[c("median","95% lower","95% upper"),"mu"])
  tres[9,c(9:11)]<-fit.bayesmeta.HN1$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.EXP_HN1adj <- bayesmeta(y=df[,"y"], 
                                        sigma=df[,"sigma"],
                                        labels=df[,"labels"],
                                        mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                        tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_EXP)})
  tres[10,c(5:7)]<-exp(fit.bayesmeta.EXP_HN1adj$summary[c("median","95% lower","95% upper"),"mu"])
  tres[10,c(9:11)]<-fit.bayesmeta.EXP_HN1adj$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.HC_HN1adj <- bayesmeta(y=df[,"y"], 
                                       sigma=df[,"sigma"],
                                       labels=df[,"labels"],
                                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                       tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HC)})
  tres[11,c(5:7)]<-exp(fit.bayesmeta.HC_HN1adj$summary[c("median","95% lower","95% upper"),"mu"])
  tres[11,c(9:11)]<-fit.bayesmeta.HC_HN1adj$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.LMX_HN1adj <- bayesmeta(y=df[,"y"], 
                                        sigma=df[,"sigma"],
                                        labels=df[,"labels"],
                                        mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                        tau.prior=function(t){dlomax(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_LMX, shape=1)})
  tres[12,c(5:7)]<-exp(fit.bayesmeta.LMX_HN1adj$summary[c("median","95% lower","95% upper"),"mu"])
  tres[12,c(9:11)]<-fit.bayesmeta.LMX_HN1adj$summary[c("median","95% lower","95% upper"),"tau"]
  
  
  # learning quantification
  
  # mu
  
  integrand_HN1_mu <- function(x) {sqrt(fit.bayesmeta.HN1$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_EXP_HN1adj_mu <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_HC_HN1adj_mu <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_LMX_HN1adj_mu <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  
  tres[9,13]<-sqrt(1-integrate(integrand_HN1_mu, lower = -Inf, upper = Inf)$value)
  tres[10,13]<-sqrt(1-integrate(integrand_EXP_HN1adj_mu, lower = -Inf, upper = Inf)$value)
  tres[11,13]<-sqrt(1-integrate(integrand_HC_HN1adj_mu, lower = -Inf, upper = Inf)$value) 
  tres[12,13]<-sqrt(1-integrate(integrand_LMX_HN1adj_mu, lower = -Inf, upper = Inf)$value) 
  
  # tau
  
  integrand_HN1_tau <- function(x) {sqrt(fit.bayesmeta.HN1$dposterior(tau=x)*dhalfnormal(x, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HN))}
  integrand_EXP_HN1adj_tau <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj$dposterior(tau=x)*dexp(x, rate=1/pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_EXP))}
  integrand_HC_HN1adj_tau <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj$dposterior(tau=x)*dhalfcauchy(x, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HC))}
  integrand_LMX_HN1adj_tau <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj$dposterior(tau=x)*dlomax(x, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_LMX, shape=1))}
  
  tres[9,14]<-sqrt(1-integrate(integrand_HN1_tau, lower = 0, upper = Inf)$value)
  tres[10,14]<-sqrt(1-integrate(integrand_EXP_HN1adj_tau, lower = 0, upper = Inf)$value)
  tres[11,14]<-sqrt(1-integrate(integrand_HC_HN1adj_tau, lower = 0, upper = Inf)$value) 
  tres[12,14]<-sqrt(1-integrate(integrand_LMX_HN1adj_tau, lower = 0, upper = Inf)$value) 
  
  
  # epsilon-local grid for A*|X| scaled distributions for tau
  
  
  gr2p_HN1adj<-pri_par_epsilon_grid(AA0_HN=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HN, 
                                    AA0_EXP=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_EXP,
                                    AA0_HC=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HC,
                                    AA0_LMX=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_LMX,
                                    grid_epsilon=grid_epsilon)
  
  # sensitivity quantification
  
  # HN
  
  fit.bayesmeta.HN1 <- bayesmeta(y=df[,"y"], 
                                 sigma=df[,"sigma"],
                                 labels=df[,"labels"],
                                 mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                 tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HN)})
  
  fit.bayesmeta.HN1_mu_l <- bayesmeta(y=df[,"y"], 
                                      sigma=df[,"sigma"],
                                      labels=df[,"labels"],
                                      mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                      tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HN)})
  
  fit.bayesmeta.HN1_mu_u <- bayesmeta(y=df[,"y"], 
                                      sigma=df[,"sigma"],
                                      labels=df[,"labels"],
                                      mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                      tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HN)})
  
  fit.bayesmeta.HN1_tau_l <- bayesmeta(y=df[,"y"], 
                                       sigma=df[,"sigma"],
                                       labels=df[,"labels"],
                                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                       tau.prior=function(t){dhalfnormal(t, scale=gr2p_HN1adj$p_HN_l)})
  
  fit.bayesmeta.HN1_tau_u <- bayesmeta(y=df[,"y"], 
                                       sigma=df[,"sigma"],
                                       labels=df[,"labels"],
                                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                       tau.prior=function(t){dhalfnormal(t, scale=gr2p_HN1adj$p_HN_u)})
  
  
  integrand_HN1_l_base_mu <- function(x) {sqrt(fit.bayesmeta.HN1$dposterior(mu=x)*fit.bayesmeta.HN1_mu_l$dposterior(mu=x))}
  sens_HN1_l_base_mu<-sqrt(1-integrate(integrand_HN1_l_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_HN1_u_base_mu <- function(x) {sqrt(fit.bayesmeta.HN1$dposterior(mu=x)*fit.bayesmeta.HN1_mu_u$dposterior(mu=x))}
  sens_HN1_u_base_mu<-sqrt(1-integrate(integrand_HN1_u_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon  
  worst_sens_HN1_mu<-max(sens_HN1_l_base_mu, sens_HN1_u_base_mu)
  
  integrand_HN1_l_base_tau <- function(x) {sqrt(fit.bayesmeta.HN1$dposterior(tau=x)*fit.bayesmeta.HN1_tau_l$dposterior(tau=x))}
  sens_HN1_l_base_tau<-sqrt(1-integrate(integrand_HN1_l_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon  
  integrand_HN1_u_base_tau <- function(x) {sqrt(fit.bayesmeta.HN1$dposterior(tau=x)*fit.bayesmeta.HN1_tau_u$dposterior(tau=x))}
  sens_HN1_u_base_tau<-sqrt(1-integrate(integrand_HN1_u_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon  
  worst_sens_HN1_tau<-max(sens_HN1_l_base_tau, sens_HN1_u_base_tau)
  
  
  # EXP
  
  fit.bayesmeta.EXP_HN1adj <- bayesmeta(y=df[,"y"], 
                                        sigma=df[,"sigma"],
                                        labels=df[,"labels"],
                                        mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                        tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_EXP)})
  
  fit.bayesmeta.EXP_HN1adj_mu_l <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                             tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_EXP)})
  
  fit.bayesmeta.EXP_HN1adj_mu_u <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                             tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_EXP)})
  
  fit.bayesmeta.EXP_HN1adj_tau_l <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                              tau.prior=function(t){dexp(t, rate=1/gr2p_HN1adj$p_EXP_l)})
  
  fit.bayesmeta.EXP_HN1adj_tau_u <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                              tau.prior=function(t){dexp(t, rate=1/gr2p_HN1adj$p_EXP_u)})
  
  
  integrand_EXP_HN1adj_l_base_mu <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj$dposterior(mu=x)*fit.bayesmeta.EXP_HN1adj_mu_l$dposterior(mu=x))}
  sens_EXP_HN1adj_l_base_mu<-sqrt(1-integrate(integrand_EXP_HN1adj_l_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_EXP_HN1adj_u_base_mu <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj$dposterior(mu=x)*fit.bayesmeta.EXP_HN1adj_mu_u$dposterior(mu=x))}
  sens_EXP_HN1adj_u_base_mu<-sqrt(1-integrate(integrand_EXP_HN1adj_u_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  worst_sens_EXP_HN1adj_mu<-max(sens_EXP_HN1adj_l_base_mu, sens_EXP_HN1adj_u_base_mu)
  
  integrand_EXP_HN1adj_l_base_tau <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj$dposterior(tau=x)*fit.bayesmeta.EXP_HN1adj_tau_l$dposterior(tau=x))}
  sens_EXP_HN1adj_l_base_tau<-sqrt(1-integrate(integrand_EXP_HN1adj_l_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_EXP_HN1adj_u_base_tau <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj$dposterior(tau=x)*fit.bayesmeta.EXP_HN1adj_tau_u$dposterior(tau=x))}
  sens_EXP_HN1adj_u_base_tau<-sqrt(1-integrate(integrand_EXP_HN1adj_u_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon
  worst_sens_EXP_HN1adj_tau<-max(sens_EXP_HN1adj_l_base_tau, sens_EXP_HN1adj_u_base_tau)
  
  
  # HC
  
  fit.bayesmeta.HC_HN1adj <- bayesmeta(y=df[,"y"], 
                                       sigma=df[,"sigma"],
                                       labels=df[,"labels"],
                                       mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                       tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HC)})
  
  fit.bayesmeta.HC_HN1adj_mu_l <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                            tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HC)})
  
  fit.bayesmeta.HC_HN1adj_mu_u <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                            tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_HC)})
  
  fit.bayesmeta.HC_HN1adj_tau_l <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                             tau.prior=function(t){dhalfcauchy(t, scale=gr2p_HN1adj$p_HC_l)})
  
  fit.bayesmeta.HC_HN1adj_tau_u <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                             tau.prior=function(t){dhalfcauchy(t, scale=gr2p_HN1adj$p_HC_u)})
  
  
  integrand_HC_HN1adj_l_base_mu <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj$dposterior(mu=x)*fit.bayesmeta.HC_HN1adj_mu_l$dposterior(mu=x))}
  sens_HC_HN1adj_l_base_mu<-sqrt(1-integrate(integrand_HC_HN1adj_l_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_HC_HN1adj_u_base_mu <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj$dposterior(mu=x)*fit.bayesmeta.HC_HN1adj_mu_u$dposterior(mu=x))}
  sens_HC_HN1adj_u_base_mu<-sqrt(1-integrate(integrand_HC_HN1adj_u_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  worst_sens_HC_HN1adj_mu<-max(sens_HC_HN1adj_l_base_mu, sens_HC_HN1adj_u_base_mu)
  
  integrand_HC_HN1adj_l_base_tau <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj$dposterior(tau=x)*fit.bayesmeta.HC_HN1adj_tau_l$dposterior(tau=x))}
  sens_HC_HN1adj_l_base_tau<-sqrt(1-integrate(integrand_HC_HN1adj_l_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_HC_HN1adj_u_base_tau <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj$dposterior(tau=x)*fit.bayesmeta.HC_HN1adj_tau_u$dposterior(tau=x))}
  sens_HC_HN1adj_u_base_tau<-sqrt(1-integrate(integrand_HC_HN1adj_u_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  worst_sens_HC_HN1adj_tau<-max(sens_HC_HN1adj_l_base_tau, sens_HC_HN1adj_u_base_tau)
  
  
  # LMX
  
  fit.bayesmeta.LMX_HN1adj <- bayesmeta(y=df[,"y"], 
                                        sigma=df[,"sigma"],
                                        labels=df[,"labels"],
                                        mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                        tau.prior=function(t){dlomax(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN1adj_mu_l <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                             tau.prior=function(t){dlomax(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN1adj_mu_u <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                             tau.prior=function(t){dlomax(t, scale=pri_par_adjust_static(UU=U2, tail_prob=tail_alpha_static)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN1adj_tau_l <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                              tau.prior=function(t){dlomax(t, scale=gr2p_HN1adj$p_LMX_l, shape=1)})
  
  fit.bayesmeta.LMX_HN1adj_tau_u <- bayesmeta(y=df[,"y"], 
                                              sigma=df[,"sigma"],
                                              labels=df[,"labels"],
                                              mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                              tau.prior=function(t){dlomax(t, scale=gr2p_HN1adj$p_LMX_u, shape=1)})
  
  
  integrand_LMX_HN1adj_l_base_mu <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj$dposterior(mu=x)*fit.bayesmeta.LMX_HN1adj_mu_l$dposterior(mu=x))}
  sens_LMX_HN1adj_l_base_mu<-sqrt(1-integrate(integrand_LMX_HN1adj_l_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_LMX_HN1adj_u_base_mu <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj$dposterior(mu=x)*fit.bayesmeta.LMX_HN1adj_mu_u$dposterior(mu=x))}
  sens_LMX_HN1adj_u_base_mu<-sqrt(1-integrate(integrand_LMX_HN1adj_u_base_mu, lower = -Inf, upper = Inf)$value)/grid_epsilon
  worst_sens_LMX_HN1adj_mu<-max(sens_LMX_HN1adj_l_base_mu, sens_LMX_HN1adj_u_base_mu)
  
  integrand_LMX_HN1adj_l_base_tau <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj$dposterior(tau=x)*fit.bayesmeta.LMX_HN1adj_tau_l$dposterior(tau=x))}
  sens_LMX_HN1adj_l_base_tau<-sqrt(1-integrate(integrand_LMX_HN1adj_l_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_LMX_HN1adj_u_base_tau <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj$dposterior(tau=x)*fit.bayesmeta.LMX_HN1adj_tau_u$dposterior(tau=x))}
  sens_LMX_HN1adj_u_base_tau<-sqrt(1-integrate(integrand_LMX_HN1adj_u_base_tau, lower = 0, upper = Inf)$value)/grid_epsilon 
  worst_sens_LMX_HN1adj_tau<-max(sens_LMX_HN1adj_l_base_tau, sens_LMX_HN1adj_u_base_tau)
  
  
  tres[9,15]<-worst_sens_HN1_mu
  tres[10,15]<-worst_sens_EXP_HN1adj_mu
  tres[11,15]<-worst_sens_HC_HN1adj_mu 
  tres[12,15]<-worst_sens_LMX_HN1adj_mu
  
  
  tres[9,16]<-worst_sens_HN1_tau
  tres[10,16]<-worst_sens_EXP_HN1adj_tau
  tres[11,16]<-worst_sens_HC_HN1adj_tau 
  tres[12,16]<-worst_sens_LMX_HN1adj_tau 
  
  
  
  ####----  dynamic 0.5-tail adjustment with rlmc=0.25 ----####
  
  # prior parameters 
  
  U_ref025<-sqrt(rlmc1/(1-rlmc1))*sigma_ref(df, type_sigma_ref=type_sigma_ref)
  
  tres[c(5:8),1] <- U_ref025
  tres[c(5:8),2] <- tail_alpha_dynamic
  tres[5,3]<-pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN
  tres[6,3]<-pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP
  tres[7,3]<-pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC
  tres[8,3]<-pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX
  
  # effective mrlmc
  
  tres[5,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN), 
                                     MM=1000000, output="sample", step=0.03))
  tres[6,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rexp(n=MM, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP), 
                                     MM=1000000, output="sample", step=0.03))
  tres[7,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rhalfcauchy(n=MM, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC), 
                                     MM=1000000, output="sample", step=0.03))
  tres[8,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rlomax(n=MM, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1), 
                                     MM=1000000, output="sample", step=0.03))
  
  
  # posteriors for base HN, EXP, HC and LMX heterogeneity priors
  
  
  fit.bayesmeta.HN05_dyn <- bayesmeta(y=df[,"y"], 
                                      sigma=df[,"sigma"],
                                      labels=df[,"labels"],
                                      mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                      tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN)})
  tres[5,c(5:7)]<-exp(fit.bayesmeta.HN05_dyn$summary[c("median","95% lower","95% upper"),"mu"])
  tres[5,c(9:11)]<-fit.bayesmeta.HN05_dyn$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.EXP_HN05adj_dyn <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                             tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP)})
  tres[6,c(5:7)]<-exp(fit.bayesmeta.EXP_HN05adj_dyn$summary[c("median","95% lower","95% upper"),"mu"])
  tres[6,c(9:11)]<-fit.bayesmeta.EXP_HN05adj_dyn$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.HC_HN05adj_dyn <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                            tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC)})
  tres[7,c(5:7)]<-exp(fit.bayesmeta.HC_HN05adj_dyn$summary[c("median","95% lower","95% upper"),"mu"])
  tres[7,c(9:11)]<-fit.bayesmeta.HC_HN05adj_dyn$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.LMX_HN05adj_dyn <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                             tau.prior=function(t){dlomax(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1)})
  tres[8,c(5:7)]<-exp(fit.bayesmeta.LMX_HN05adj_dyn$summary[c("median","95% lower","95% upper"),"mu"])
  tres[8,c(9:11)]<-fit.bayesmeta.LMX_HN05adj_dyn$summary[c("median","95% lower","95% upper"),"tau"]
  
  
  # learning quantification
  
  # mu
  
  integrand_HN05_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HN05_dyn$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_EXP_HN05adj_mu_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj_dyn$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_HC_HN05adj_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj_dyn$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_LMX_HN05adj_mu_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj_dyn$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  
  
  tres[5,13]<-sqrt(1-integrate(integrand_HN05_mu_dyn, lower = -Inf, upper = Inf)$value)
  tres[6,13]<-sqrt(1-integrate(integrand_EXP_HN05adj_mu_dyn, lower = -Inf, upper = Inf)$value)
  tres[7,13]<-sqrt(1-integrate(integrand_HC_HN05adj_mu_dyn, lower = -Inf, upper = Inf)$value) 
  tres[8,13]<-sqrt(1-integrate(integrand_LMX_HN05adj_mu_dyn, lower = -Inf, upper = Inf)$value) 
  
  # tau
  
  integrand_HN05_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HN05_dyn$dposterior(tau=x)*dhalfnormal(x, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN))}
  integrand_EXP_HN05adj_tau_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj_dyn$dposterior(tau=x)*dexp(x, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP))}
  integrand_HC_HN05adj_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj_dyn$dposterior(tau=x)*dhalfcauchy(x, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC))}
  integrand_LMX_HN05adj_tau_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj_dyn$dposterior(tau=x)*dlomax(x, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1))}
  
  
  tres[5,14]<-sqrt(1-integrate(integrand_HN05_tau_dyn, lower = 0, upper = Inf)$value)
  tres[6,14]<-sqrt(1-integrate(integrand_EXP_HN05adj_tau_dyn, lower = 0, upper = Inf)$value)
  tres[7,14]<-sqrt(1-integrate(integrand_HC_HN05adj_tau_dyn, lower = 0, upper = Inf)$value) 
  tres[8,14]<-sqrt(1-integrate(integrand_LMX_HN05adj_tau_dyn, lower = 0, upper = Inf)$value) 
  
  # epsilon-local grid for A*|X| scaled distributions for tau
  
  gr2p_HN05adj_dyn<-pri_par_epsilon_grid(AA0_HN=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN, 
                                         AA0_EXP=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP,
                                         AA0_HC=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC,
                                         AA0_LMX=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX,
                                         grid_epsilon=grid_epsilon)
  
  
  # sensitivity quantification
  
  # HN
  
  fit.bayesmeta.HN05_dyn <- bayesmeta(y=df[,"y"], 
                                      sigma=df[,"sigma"],
                                      labels=df[,"labels"],
                                      mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                      tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN)})
  
  fit.bayesmeta.HN05_mu_l_dyn <- bayesmeta(y=df[,"y"], 
                                           sigma=df[,"sigma"],
                                           labels=df[,"labels"],
                                           mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                           tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN)})
  
  fit.bayesmeta.HN05_mu_u_dyn <- bayesmeta(y=df[,"y"], 
                                           sigma=df[,"sigma"],
                                           labels=df[,"labels"],
                                           mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                           tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN)})
  
  fit.bayesmeta.HN05_tau_l_dyn <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                            tau.prior=function(t){dhalfnormal(t, scale=gr2p_HN05adj_dyn$p_HN_l)})
  
  fit.bayesmeta.HN05_tau_u_dyn <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                            tau.prior=function(t){dhalfnormal(t, scale=gr2p_HN05adj_dyn$p_HN_u)})
  
  
  integrand_HN05_l_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HN05_dyn$dposterior(mu=x)*fit.bayesmeta.HN05_mu_l_dyn$dposterior(mu=x))}
  sens_HN05_l_base_mu_dyn<-sqrt(1-integrate(integrand_HN05_l_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_HN05_u_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HN05_dyn$dposterior(mu=x)*fit.bayesmeta.HN05_mu_u_dyn$dposterior(mu=x))}
  sens_HN05_u_base_mu_dyn<-sqrt(1-integrate(integrand_HN05_u_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon  
  worst_sens_HN05_mu_dyn<-max(sens_HN05_l_base_mu_dyn, sens_HN05_u_base_mu_dyn)
  
  integrand_HN05_l_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HN05_dyn$dposterior(tau=x)*fit.bayesmeta.HN05_tau_l_dyn$dposterior(tau=x))}
  sens_HN05_l_base_tau_dyn<-sqrt(1-integrate(integrand_HN05_l_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon  
  integrand_HN05_u_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HN05_dyn$dposterior(tau=x)*fit.bayesmeta.HN05_tau_u_dyn$dposterior(tau=x))}
  sens_HN05_u_base_tau_dyn<-sqrt(1-integrate(integrand_HN05_u_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon  
  worst_sens_HN05_tau_dyn<-max(sens_HN05_l_base_tau_dyn, sens_HN05_u_base_tau_dyn)
  
  # EXP
  
  fit.bayesmeta.EXP_HN05adj_dyn <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                             tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP)})
  
  fit.bayesmeta.EXP_HN05adj_mu_l_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                                  tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP)})
  
  fit.bayesmeta.EXP_HN05adj_mu_u_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                                  tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP)})
  
  fit.bayesmeta.EXP_HN05adj_tau_l_dyn <- bayesmeta(y=df[,"y"], 
                                                   sigma=df[,"sigma"],
                                                   labels=df[,"labels"],
                                                   mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                   tau.prior=function(t){dexp(t, rate=1/gr2p_HN05adj_dyn$p_EXP_l)})
  
  fit.bayesmeta.EXP_HN05adj_tau_u_dyn <- bayesmeta(y=df[,"y"], 
                                                   sigma=df[,"sigma"],
                                                   labels=df[,"labels"],
                                                   mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                   tau.prior=function(t){dexp(t, rate=1/gr2p_HN05adj_dyn$p_EXP_u)})
  
  
  integrand_EXP_HN05adj_l_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj_dyn$dposterior(mu=x)*fit.bayesmeta.EXP_HN05adj_mu_l_dyn$dposterior(mu=x))}
  sens_EXP_HN05adj_l_base_mu_dyn<-sqrt(1-integrate(integrand_EXP_HN05adj_l_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_EXP_HN05adj_u_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj_dyn$dposterior(mu=x)*fit.bayesmeta.EXP_HN05adj_mu_u_dyn$dposterior(mu=x))}
  sens_EXP_HN05adj_u_base_mu_dyn<-sqrt(1-integrate(integrand_EXP_HN05adj_u_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  worst_sens_EXP_HN05adj_mu_dyn<-max(sens_EXP_HN05adj_l_base_mu_dyn, sens_EXP_HN05adj_u_base_mu_dyn)
  
  integrand_EXP_HN05adj_l_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj_dyn$dposterior(tau=x)*fit.bayesmeta.EXP_HN05adj_tau_l_dyn$dposterior(tau=x))}
  sens_EXP_HN05adj_l_base_tau_dyn<-sqrt(1-integrate(integrand_EXP_HN05adj_l_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_EXP_HN05adj_u_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN05adj_dyn$dposterior(tau=x)*fit.bayesmeta.EXP_HN05adj_tau_u_dyn$dposterior(tau=x))}
  sens_EXP_HN05adj_u_base_tau_dyn<-sqrt(1-integrate(integrand_EXP_HN05adj_u_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon
  worst_sens_EXP_HN05adj_tau_dyn<-max(sens_EXP_HN05adj_l_base_tau_dyn, sens_EXP_HN05adj_u_base_tau_dyn)
  
  # HC
  
  fit.bayesmeta.HC_HN05adj_dyn <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                            tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC)})
  
  fit.bayesmeta.HC_HN05adj_mu_l_dyn <- bayesmeta(y=df[,"y"], 
                                                 sigma=df[,"sigma"],
                                                 labels=df[,"labels"],
                                                 mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                                 tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC)})
  
  fit.bayesmeta.HC_HN05adj_mu_u_dyn <- bayesmeta(y=df[,"y"], 
                                                 sigma=df[,"sigma"],
                                                 labels=df[,"labels"],
                                                 mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                                 tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC)})
  
  fit.bayesmeta.HC_HN05adj_tau_l_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                  tau.prior=function(t){dhalfcauchy(t, scale=gr2p_HN05adj_dyn$p_HC_l)})
  
  fit.bayesmeta.HC_HN05adj_tau_u_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                  tau.prior=function(t){dhalfcauchy(t, scale=gr2p_HN05adj_dyn$p_HC_u)})
  
  
  integrand_HC_HN05adj_l_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj_dyn$dposterior(mu=x)*fit.bayesmeta.HC_HN05adj_mu_l_dyn$dposterior(mu=x))}
  sens_HC_HN05adj_l_base_mu_dyn<-sqrt(1-integrate(integrand_HC_HN05adj_l_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_HC_HN05adj_u_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj_dyn$dposterior(mu=x)*fit.bayesmeta.HC_HN05adj_mu_u_dyn$dposterior(mu=x))}
  sens_HC_HN05adj_u_base_mu_dyn<-sqrt(1-integrate(integrand_HC_HN05adj_u_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  worst_sens_HC_HN05adj_mu_dyn<-max(sens_HC_HN05adj_l_base_mu_dyn, sens_HC_HN05adj_u_base_mu_dyn)
  
  integrand_HC_HN05adj_l_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj_dyn$dposterior(tau=x)*fit.bayesmeta.HC_HN05adj_tau_l_dyn$dposterior(tau=x))}
  sens_HC_HN05adj_l_base_tau_dyn<-sqrt(1-integrate(integrand_HC_HN05adj_l_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_HC_HN05adj_u_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN05adj_dyn$dposterior(tau=x)*fit.bayesmeta.HC_HN05adj_tau_u_dyn$dposterior(tau=x))}
  sens_HC_HN05adj_u_base_tau_dyn<-sqrt(1-integrate(integrand_HC_HN05adj_u_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  worst_sens_HC_HN05adj_tau_dyn<-max(sens_HC_HN05adj_l_base_tau_dyn, sens_HC_HN05adj_u_base_tau_dyn)
  
  # LMX
  
  fit.bayesmeta.LMX_HN05adj_dyn <- bayesmeta(y=df[,"y"], 
                                             sigma=df[,"sigma"],
                                             labels=df[,"labels"],
                                             mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                             tau.prior=function(t){dlomax(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN05adj_mu_l_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                                  tau.prior=function(t){dlomax(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN05adj_mu_u_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                                  tau.prior=function(t){dlomax(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc1, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN05adj_tau_l_dyn <- bayesmeta(y=df[,"y"], 
                                                   sigma=df[,"sigma"],
                                                   labels=df[,"labels"],
                                                   mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                   tau.prior=function(t){dlomax(t, scale=gr2p_HN05adj_dyn$p_LMX_l, shape=1)})
  
  fit.bayesmeta.LMX_HN05adj_tau_u_dyn <- bayesmeta(y=df[,"y"], 
                                                   sigma=df[,"sigma"],
                                                   labels=df[,"labels"],
                                                   mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                   tau.prior=function(t){dlomax(t, scale=gr2p_HN05adj_dyn$p_LMX_u, shape=1)})
  
  
  integrand_LMX_HN05adj_l_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj_dyn$dposterior(mu=x)*fit.bayesmeta.LMX_HN05adj_mu_l_dyn$dposterior(mu=x))}
  sens_LMX_HN05adj_l_base_mu_dyn<-sqrt(1-integrate(integrand_LMX_HN05adj_l_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_LMX_HN05adj_u_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj_dyn$dposterior(mu=x)*fit.bayesmeta.LMX_HN05adj_mu_u_dyn$dposterior(mu=x))}
  sens_LMX_HN05adj_u_base_mu_dyn<-sqrt(1-integrate(integrand_LMX_HN05adj_u_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon
  worst_sens_LMX_HN05adj_mu_dyn<-max(sens_LMX_HN05adj_l_base_mu_dyn, sens_LMX_HN05adj_u_base_mu_dyn)
  
  integrand_LMX_HN05adj_l_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj_dyn$dposterior(tau=x)*fit.bayesmeta.LMX_HN05adj_tau_l_dyn$dposterior(tau=x))}
  sens_LMX_HN05adj_l_base_tau_dyn<-sqrt(1-integrate(integrand_LMX_HN05adj_l_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_LMX_HN05adj_u_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN05adj_dyn$dposterior(tau=x)*fit.bayesmeta.LMX_HN05adj_tau_u_dyn$dposterior(tau=x))}
  sens_LMX_HN05adj_u_base_tau_dyn<-sqrt(1-integrate(integrand_LMX_HN05adj_u_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  worst_sens_LMX_HN05adj_tau_dyn<-max(sens_LMX_HN05adj_l_base_tau_dyn, sens_LMX_HN05adj_u_base_tau_dyn)
  
  
  tres[5,15]<-worst_sens_HN05_mu_dyn
  tres[6,15]<-worst_sens_EXP_HN05adj_mu_dyn
  tres[7,15]<-worst_sens_HC_HN05adj_mu_dyn 
  tres[8,15]<-worst_sens_LMX_HN05adj_mu_dyn
  
  tres[5,16]<-worst_sens_HN05_tau_dyn
  tres[6,16]<-worst_sens_EXP_HN05adj_tau_dyn
  tres[7,16]<-worst_sens_HC_HN05adj_tau_dyn  
  tres[8,16]<-worst_sens_LMX_HN05adj_tau_dyn
  
  
  
  
  
  ####----  dynamic 0.5-tail adjustment with rlmc=0.5 ----####
  
  # prior parameters 
  
  U_ref05<-sqrt(rlmc2/(1-rlmc2))*sigma_ref(df, type_sigma_ref=type_sigma_ref)
  
  tres[c(13:16),1]<-U_ref05
  tres[c(13:16),2]<-tail_alpha_dynamic
  tres[13,3]<-pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN
  tres[14,3]<-pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP
  tres[15,3]<-pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC
  tres[16,3]<-pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX
  
  # effective mrlmc
  
  tres[13,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rhalfnormal(n=MM, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN), 
                                      MM=1000000, output="sample", step=0.03))
  tres[14,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rexp(n=MM, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP), 
                                      MM=1000000, output="sample", step=0.03))
  tres[15,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rhalfcauchy(n=MM, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC), 
                                      MM=1000000, output="sample", step=0.03))
  tres[16,4]<- median(effective_rlmc(df=df, r.tau.prior=function(MM)rlomax(n=MM, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1), 
                                      MM=1000000, output="sample", step=0.03))
  
  
  # posteriors for base HN, EXP, HC and LMX heterogeneity priors
  
  
  fit.bayesmeta.HN1_dyn <- bayesmeta(y=df[,"y"], 
                                     sigma=df[,"sigma"],
                                     labels=df[,"labels"],
                                     mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                     tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN)})
  tres[13,c(5:7)]<-exp(fit.bayesmeta.HN1_dyn$summary[c("median","95% lower","95% upper"),"mu"])
  tres[13,c(9:11)]<-fit.bayesmeta.HN1_dyn$summary[c("median","95% lower","95% upper"),"tau"]
  
  
  fit.bayesmeta.EXP_HN1adj_dyn <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                            tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP)})
  tres[14,c(5:7)]<-exp(fit.bayesmeta.EXP_HN1adj_dyn$summary[c("median","95% lower","95% upper"),"mu"])
  tres[14,c(9:11)]<-fit.bayesmeta.EXP_HN1adj_dyn$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.HC_HN1adj_dyn <- bayesmeta(y=df[,"y"], 
                                           sigma=df[,"sigma"],
                                           labels=df[,"labels"],
                                           mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                           tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC)})
  tres[15,c(5:7)]<-exp(fit.bayesmeta.HC_HN1adj_dyn$summary[c("median","95% lower","95% upper"),"mu"])
  tres[15,c(9:11)]<-fit.bayesmeta.HC_HN1adj_dyn$summary[c("median","95% lower","95% upper"),"tau"]
  
  fit.bayesmeta.LMX_HN1adj_dyn <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                            tau.prior=function(t){dlomax(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1)})
  tres[16,c(5:7)]<-exp(fit.bayesmeta.LMX_HN1adj_dyn$summary[c("median","95% lower","95% upper"),"mu"])
  tres[16,c(9:11)]<-fit.bayesmeta.LMX_HN1adj_dyn$summary[c("median","95% lower","95% upper"),"tau"]
  
  
  # learning quantification
  
  # mu
  
  integrand_HN1_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HN1_dyn$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_EXP_HN1adj_mu_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj_dyn$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_HC_HN1adj_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj_dyn$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  integrand_LMX_HN1adj_mu_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj_dyn$dposterior(mu=x)*dnorm(x, mean=mu_mean, sd=mu_sd))}
  
  
  tres[13,13]<-sqrt(1-integrate(integrand_HN1_mu_dyn, lower = -Inf, upper = Inf)$value)
  tres[14,13]<-sqrt(1-integrate(integrand_EXP_HN1adj_mu_dyn, lower = -Inf, upper = Inf)$value)
  tres[15,13]<-sqrt(1-integrate(integrand_HC_HN1adj_mu_dyn, lower = -Inf, upper = Inf)$value) 
  tres[16,13]<-sqrt(1-integrate(integrand_LMX_HN1adj_mu_dyn, lower = -Inf, upper = Inf)$value) 
  
  # tau
  
  integrand_HN1_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HN1_dyn$dposterior(tau=x)*dhalfnormal(x, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN))}
  integrand_EXP_HN1adj_tau_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj_dyn$dposterior(tau=x)*dexp(x, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP))}
  integrand_HC_HN1adj_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj_dyn$dposterior(tau=x)*dhalfcauchy(x, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC))}
  integrand_LMX_HN1adj_tau_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj_dyn$dposterior(tau=x)*dlomax(x, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1))}
  
  
  tres[13,14]<-sqrt(1-integrate(integrand_HN1_tau_dyn, lower = 0, upper = Inf)$value)
  tres[14,14]<-sqrt(1-integrate(integrand_EXP_HN1adj_tau_dyn, lower = 0, upper = Inf)$value)
  tres[15,14]<-sqrt(1-integrate(integrand_HC_HN1adj_tau_dyn, lower = 0, upper = Inf)$value) 
  tres[16,14]<-sqrt(1-integrate(integrand_LMX_HN1adj_tau_dyn, lower = 0, upper = Inf)$value) 
  
  
  # epsilon-local grid for A*|X| scaled distributions for tau
  
  
  gr2p_HN1adj_dyn<-pri_par_epsilon_grid(AA0_HN=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN, 
                                        AA0_EXP=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP,
                                        AA0_HC=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC,
                                        AA0_LMX=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX,
                                        grid_epsilon=grid_epsilon)
  
  
  # sensitivity quantification
  
  # HN
  
  fit.bayesmeta.HN1_dyn <- bayesmeta(y=df[,"y"], 
                                     sigma=df[,"sigma"],
                                     labels=df[,"labels"],
                                     mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                     tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN)})
  
  fit.bayesmeta.HN1_mu_l_dyn <- bayesmeta(y=df[,"y"], 
                                          sigma=df[,"sigma"],
                                          labels=df[,"labels"],
                                          mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                          tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN)})
  
  fit.bayesmeta.HN1_mu_u_dyn <- bayesmeta(y=df[,"y"], 
                                          sigma=df[,"sigma"],
                                          labels=df[,"labels"],
                                          mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                          tau.prior=function(t){dhalfnormal(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HN)})
  
  fit.bayesmeta.HN1_tau_l_dyn <- bayesmeta(y=df[,"y"], 
                                           sigma=df[,"sigma"],
                                           labels=df[,"labels"],
                                           mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                           tau.prior=function(t){dhalfnormal(t, scale=gr2p_HN1adj_dyn$p_HN_l)})
  
  fit.bayesmeta.HN1_tau_u_dyn <- bayesmeta(y=df[,"y"], 
                                           sigma=df[,"sigma"],
                                           labels=df[,"labels"],
                                           mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                           tau.prior=function(t){dhalfnormal(t, scale=gr2p_HN1adj_dyn$p_HN_u)})
  
  
  integrand_HN1_l_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HN1_dyn$dposterior(mu=x)*fit.bayesmeta.HN1_mu_l_dyn$dposterior(mu=x))}
  sens_HN1_l_base_mu_dyn<-sqrt(1-integrate(integrand_HN1_l_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_HN1_u_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HN1_dyn$dposterior(mu=x)*fit.bayesmeta.HN1_mu_u_dyn$dposterior(mu=x))}
  sens_HN1_u_base_mu_dyn<-sqrt(1-integrate(integrand_HN1_u_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon  
  worst_sens_HN1_mu_dyn<-max(sens_HN1_l_base_mu_dyn, sens_HN1_u_base_mu_dyn)
  
  integrand_HN1_l_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HN1_dyn$dposterior(tau=x)*fit.bayesmeta.HN1_tau_l_dyn$dposterior(tau=x))}
  sens_HN1_l_base_tau_dyn<-sqrt(1-integrate(integrand_HN1_l_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon  
  integrand_HN1_u_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HN1_dyn$dposterior(tau=x)*fit.bayesmeta.HN1_tau_u_dyn$dposterior(tau=x))}
  sens_HN1_u_base_tau_dyn<-sqrt(1-integrate(integrand_HN1_u_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon  
  worst_sens_HN1_tau_dyn<-max(sens_HN1_l_base_tau_dyn, sens_HN1_u_base_tau_dyn)
  
  # EXP
  
  fit.bayesmeta.EXP_HN1adj_dyn <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                            tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP)})
  
  fit.bayesmeta.EXP_HN1adj_mu_l_dyn <- bayesmeta(y=df[,"y"], 
                                                 sigma=df[,"sigma"],
                                                 labels=df[,"labels"],
                                                 mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                                 tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP)})
  
  fit.bayesmeta.EXP_HN1adj_mu_u_dyn <- bayesmeta(y=df[,"y"], 
                                                 sigma=df[,"sigma"],
                                                 labels=df[,"labels"],
                                                 mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                                 tau.prior=function(t){dexp(t, rate=1/pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_EXP)})
  
  fit.bayesmeta.EXP_HN1adj_tau_l_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                  tau.prior=function(t){dexp(t, rate=1/gr2p_HN1adj_dyn$p_EXP_l)})
  
  fit.bayesmeta.EXP_HN1adj_tau_u_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                  tau.prior=function(t){dexp(t, rate=1/gr2p_HN1adj_dyn$p_EXP_u)})
  
  
  integrand_EXP_HN1adj_l_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj_dyn$dposterior(mu=x)*fit.bayesmeta.EXP_HN1adj_mu_l_dyn$dposterior(mu=x))}
  sens_EXP_HN1adj_l_base_mu_dyn<-sqrt(1-integrate(integrand_EXP_HN1adj_l_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_EXP_HN1adj_u_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj_dyn$dposterior(mu=x)*fit.bayesmeta.EXP_HN1adj_mu_u_dyn$dposterior(mu=x))}
  sens_EXP_HN1adj_u_base_mu_dyn<-sqrt(1-integrate(integrand_EXP_HN1adj_u_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  worst_sens_EXP_HN1adj_mu_dyn<-max(sens_EXP_HN1adj_l_base_mu_dyn, sens_EXP_HN1adj_u_base_mu_dyn)
  
  integrand_EXP_HN1adj_l_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj_dyn$dposterior(tau=x)*fit.bayesmeta.EXP_HN1adj_tau_l_dyn$dposterior(tau=x))}
  sens_EXP_HN1adj_l_base_tau_dyn<-sqrt(1-integrate(integrand_EXP_HN1adj_l_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_EXP_HN1adj_u_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.EXP_HN1adj_dyn$dposterior(tau=x)*fit.bayesmeta.EXP_HN1adj_tau_u_dyn$dposterior(tau=x))}
  sens_EXP_HN1adj_u_base_tau_dyn<-sqrt(1-integrate(integrand_EXP_HN1adj_u_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon
  worst_sens_EXP_HN1adj_tau_dyn<-max(sens_EXP_HN1adj_l_base_tau_dyn, sens_EXP_HN1adj_u_base_tau_dyn)
  
  
  # HC
  
  fit.bayesmeta.HC_HN1adj_dyn <- bayesmeta(y=df[,"y"], 
                                           sigma=df[,"sigma"],
                                           labels=df[,"labels"],
                                           mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                           tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC)})
  
  fit.bayesmeta.HC_HN1adj_mu_l_dyn <- bayesmeta(y=df[,"y"], 
                                                sigma=df[,"sigma"],
                                                labels=df[,"labels"],
                                                mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                                tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC)})
  
  fit.bayesmeta.HC_HN1adj_mu_u_dyn <- bayesmeta(y=df[,"y"], 
                                                sigma=df[,"sigma"],
                                                labels=df[,"labels"],
                                                mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                                tau.prior=function(t){dhalfcauchy(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_HC)})
  
  fit.bayesmeta.HC_HN1adj_tau_l_dyn <- bayesmeta(y=df[,"y"], 
                                                 sigma=df[,"sigma"],
                                                 labels=df[,"labels"],
                                                 mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                 tau.prior=function(t){dhalfcauchy(t, scale=gr2p_HN1adj_dyn$p_HC_l)})
  
  fit.bayesmeta.HC_HN1adj_tau_u_dyn <- bayesmeta(y=df[,"y"], 
                                                 sigma=df[,"sigma"],
                                                 labels=df[,"labels"],
                                                 mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                 tau.prior=function(t){dhalfcauchy(t, scale=gr2p_HN1adj_dyn$p_HC_u)})
  
  
  integrand_HC_HN1adj_l_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj_dyn$dposterior(mu=x)*fit.bayesmeta.HC_HN1adj_mu_l_dyn$dposterior(mu=x))}
  sens_HC_HN1adj_l_base_mu_dyn<-sqrt(1-integrate(integrand_HC_HN1adj_l_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_HC_HN1adj_u_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj_dyn$dposterior(mu=x)*fit.bayesmeta.HC_HN1adj_mu_u_dyn$dposterior(mu=x))}
  sens_HC_HN1adj_u_base_mu_dyn<-sqrt(1-integrate(integrand_HC_HN1adj_u_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  worst_sens_HC_HN1adj_mu_dyn<-max(sens_HC_HN1adj_l_base_mu_dyn, sens_HC_HN1adj_u_base_mu_dyn)
  
  integrand_HC_HN1adj_l_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj_dyn$dposterior(tau=x)*fit.bayesmeta.HC_HN1adj_tau_l_dyn$dposterior(tau=x))}
  sens_HC_HN1adj_l_base_tau_dyn<-sqrt(1-integrate(integrand_HC_HN1adj_l_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_HC_HN1adj_u_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.HC_HN1adj_dyn$dposterior(tau=x)*fit.bayesmeta.HC_HN1adj_tau_u_dyn$dposterior(tau=x))}
  sens_HC_HN1adj_u_base_tau_dyn<-sqrt(1-integrate(integrand_HC_HN1adj_u_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  worst_sens_HC_HN1adj_tau_dyn<-max(sens_HC_HN1adj_l_base_tau_dyn, sens_HC_HN1adj_u_base_tau_dyn)
  
  
  # LMX
  
  fit.bayesmeta.LMX_HN1adj_dyn <- bayesmeta(y=df[,"y"], 
                                            sigma=df[,"sigma"],
                                            labels=df[,"labels"],
                                            mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                            tau.prior=function(t){dlomax(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN1adj_mu_l_dyn <- bayesmeta(y=df[,"y"], 
                                                 sigma=df[,"sigma"],
                                                 labels=df[,"labels"],
                                                 mu.prior.mean=mu_mean, mu.prior.sd=AAl_normal_mu,
                                                 tau.prior=function(t){dlomax(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN1adj_mu_u_dyn <- bayesmeta(y=df[,"y"], 
                                                 sigma=df[,"sigma"],
                                                 labels=df[,"labels"],
                                                 mu.prior.mean=mu_mean, mu.prior.sd=AAu_normal_mu,
                                                 tau.prior=function(t){dlomax(t, scale=pri_par_adjust_dynamic(df=df, rlmc=rlmc2, tail_prob=tail_alpha_dynamic, type_sigma_ref=type_sigma_ref)$p_LMX, shape=1)})
  
  fit.bayesmeta.LMX_HN1adj_tau_l_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                  tau.prior=function(t){dlomax(t, scale=gr2p_HN1adj_dyn$p_LMX_l, shape=1)})
  
  fit.bayesmeta.LMX_HN1adj_tau_u_dyn <- bayesmeta(y=df[,"y"], 
                                                  sigma=df[,"sigma"],
                                                  labels=df[,"labels"],
                                                  mu.prior.mean=mu_mean, mu.prior.sd=mu_sd,
                                                  tau.prior=function(t){dlomax(t, scale=gr2p_HN1adj_dyn$p_LMX_u, shape=1)})
  
  
  integrand_LMX_HN1adj_l_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj_dyn$dposterior(mu=x)*fit.bayesmeta.LMX_HN1adj_mu_l_dyn$dposterior(mu=x))}
  sens_LMX_HN1adj_l_base_mu_dyn<-sqrt(1-integrate(integrand_LMX_HN1adj_l_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon 
  integrand_LMX_HN1adj_u_base_mu_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj_dyn$dposterior(mu=x)*fit.bayesmeta.LMX_HN1adj_mu_u_dyn$dposterior(mu=x))}
  sens_LMX_HN1adj_u_base_mu_dyn<-sqrt(1-integrate(integrand_LMX_HN1adj_u_base_mu_dyn, lower = -Inf, upper = Inf)$value)/grid_epsilon
  worst_sens_LMX_HN1adj_mu_dyn<-max(sens_LMX_HN1adj_l_base_mu_dyn, sens_LMX_HN1adj_u_base_mu_dyn)
  
  integrand_LMX_HN1adj_l_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj_dyn$dposterior(tau=x)*fit.bayesmeta.LMX_HN1adj_tau_l_dyn$dposterior(tau=x))}
  sens_LMX_HN1adj_l_base_tau_dyn<-sqrt(1-integrate(integrand_LMX_HN1adj_l_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  integrand_LMX_HN1adj_u_base_tau_dyn <- function(x) {sqrt(fit.bayesmeta.LMX_HN1adj_dyn$dposterior(tau=x)*fit.bayesmeta.LMX_HN1adj_tau_u_dyn$dposterior(tau=x))}
  sens_LMX_HN1adj_u_base_tau_dyn<-sqrt(1-integrate(integrand_LMX_HN1adj_u_base_tau_dyn, lower = 0, upper = Inf)$value)/grid_epsilon 
  worst_sens_LMX_HN1adj_tau_dyn<-max(sens_LMX_HN1adj_l_base_tau_dyn, sens_LMX_HN1adj_u_base_tau_dyn)
  
  
  tres[13,15]<-worst_sens_HN1_mu_dyn
  tres[14,15]<-worst_sens_EXP_HN1adj_mu_dyn
  tres[15,15]<-worst_sens_HC_HN1adj_mu_dyn 
  tres[16,15]<-worst_sens_LMX_HN1adj_mu_dyn
  
  tres[13,16]<-worst_sens_HN1_tau_dyn
  tres[14,16]<-worst_sens_EXP_HN1adj_tau_dyn
  tres[15,16]<-worst_sens_HC_HN1adj_tau_dyn 
  tres[16,16]<-worst_sens_LMX_HN1adj_tau_dyn 
  
  
  #### Computation of the length of 95%CrI
  
  tres[,8]<-tres[,7]-tres[,6] #length_95CrI_post_OR
  tres[,12]<-tres[,11]-tres[,10] #length_95CrI_post_tau
  
  return(tres)
}
