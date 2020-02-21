\name{sensitivity_learning_table_flexible}
\alias{sensitivity_learning_table_flexible}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior, sensitivity and learning estimates for Bayesian meta-analysis: Flexible version of the basic \code{sensitivity_learning_table} function
%%  ~~function to do ... ~~
}
\description{
Generates a table containing posterior, sensitivity and learning estimates and
effective RLMC 
values
for a meta-analysis data set with odds (ratios) as effects.
The effect estimates returned are either given on the odds (ratio) scale (the default)
or on the log-odds (ratio) scale.
This function supports computation of the reference within-study standard deviation
based on both the geometric mean (the default) and a weighted harmonic mean.
Assumes a Bayesian normal-normal hierarchical model with different 
priors (half-normal, exponential, half-Cauchy and Lomax) for the between-study standard deviation.
These priors are either \code{tail_alpha_static}-tail adjusted or adjusted to a target RLMC value
(with respect to a tail probability \code{tail_alpha_dynamic}).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sensitivity_learning_table_flexible(df, tail_alpha_static=0.04550026, 
                                    U1=1, U2=2, 
                                    tail_alpha_dynamic=0.5, 
                                    rlmc1=0.25, rlmc2=0.5,
                                    scale_effect="OR",
                                    type_sigma_ref="geometric", 
                                    mu_mean=0, mu_sd=4, 
                                    grid_epsilon=0.00354)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  a data frame containing a column "y" with the estimates of the log-odds (ratios) for the individual studies, 
  a column "sigma" with the corresponding standard errors and a column "labels" with labels for the studies
%%     ~~Describe \code{x} here~~
}
  \item{tail_alpha_static}{
  tail probability for the static heterogeneity prior tail-adjustment
}
\item{U1}{
  first reference threshold for the static tail-adjustment
}
\item{U2}{
  second reference threshold for the static tail-adjustment
}
\item{tail_alpha_dynamic}{
  tail probability for the dynamic heterogeneity prior adjustment based on RLMC. The default 0.5 corresponds to a median-based adjustment.
}
\item{rlmc1}{
  first target value for RLMC
}
\item{rlmc2}{
  second target value for RLMC
}
\item{scale_effect}{
  either \code{"OR"} or \code{"logOR"}. Defaults to \code{"OR"}.
  Specifies the scale for effect estimates returned: either the odds (ratio)
  or log-odds (ratio) scale.
 % Note that OR=exp(mu).
}
\item{type_sigma_ref}{ 
  either \code{"geometric"} or \code{"harmonic"}. Defaults to \code{"geometric"}.
   Specifies if the geometric mean or a weighted harmonic mean is used to compute
   the reference standard deviation. See \code{\link{sigma_ref}} for details.
}
\item{mu_mean}{
  mean of the normal prior for the effect mu on the log-odds (ratio) scale. Defaults to 0.
}
\item{mu_sd}{
  standard deviation of the normal prior for the effect mu on the log-odds (ratio) scale. Defaults to 4.
}
\item{grid_epsilon}{
%  value for epsilon, the Hellinger distance between the base prior and the
% "lower"/"upper" prior on the grid
  step size of the grid for the epsilon-local sensitivity computation. \cr
  See \code{\link{pri_par_epsilon_grid}} for details.
}
}
\details{
The sensitivity measure used is epsilon-local sensitivity (Roos et al. 2015), which is based on the
Hellinger distance.
Here, learning refers to the ability of the data to modify the prior.
It is quantified by the Hellinger distance between the prior and its marginal posterior.
See Ott et al. (2020) for a description of the implemented methodology and some examples.
Since the sensitivity and learning measures are invariant
with respect to monotone transformations of the parameter,
the sensitvity and learning estimates do not change if
we switch from the odds (ratio) scale to the log-odds (ratio) scale or
vice versa,
i.e. S(OR)=S(mu) and L(OR)=L(mu).

The shape parameter of the Lomax prior is fixed at 1.
%The reference thresholds used for the 5 \% tail-adjustment are U=1 and 
%U=2 and the target RLMC values considered are RLMC=0.25 and RLMC=0.5.
The posterior estimates are obtained from the \code{bayesmeta()} function
in the package \code{bayesmeta}.

This function is a generalized version of the basic \code{sensitivity_learning_table()}
function, which enables a flexible choice of multiple parameters
that are fixed in the latter function.
The default values of the parameters correspond
to the values used in the basic \cr
\code{sensitivity_learning_table}()
function.
These values have been chosen as follows:
for \code{tail_alpha_static}, the default value is based 
on the half-normal prior with scale parameter 0.5 and the reference threshold
U1=1: we have P[HN(0.5)>1]=0.04550026.
For the mean and the standard deviation 
of the normal prior for the effect mu, we use the default values
 recommended by Roever (2018) for outcomes on the log-odds ratio scale,
which are also suitable for the log-odds scale.
For \code{grid_epsilon}, the default value 0.00354 corresponds to the
Hellinger distance between two unit-variance normal distributions
with means 0 and 0.01 (Roos et al., 2015).

%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with the different estimates in the columns and
one row per heterogeneity prior (i.e. the prior for the between-study standard deviation).
If \code{scale_effect="OR"},
the quantities given in the columns are as follows, where 
OR=exp(mu) denotes the effect on the odds ratio or odds scale and 
tau the heterogeneity standard deviation:
\item{U }{reference threshold for prior adjustment}
\item{tail_prob }{tail probability for prior adjustment}
\item{par_val }{scale parameter value of the heterogeneity prior}
\item{MRLMC }{median relative latent model complexity estimated from MC sampling}
\item{median_post_OR }{posterior median for the effect OR-exp(mu)}
\item{95CrI_post_OR_low }{lower end point of the 95 \% shortest credible interval (CrI) for the effect OR}
\item{95CrI_post_OR_up }{upper end point of the 95 \% shortest CrI for the effect OR}
\item{length_95CrI_post_OR }{length of the 95 \% shortest CrI for the effect OR}
\item{median_post_tau }{posterior median for tau}
\item{95CrI_post_tau_low}{lower end point of the 95 \% shortest CrI for tau}
\item{95CrI_post_tau_up }{upper end point of the 95 \% shortest CrI for tau}
\item{length_95CrI_post_tau }{length of the 95 \% shortest CrI for tau}
\item{L_mu }{learning estimate for the effect mu (\code{L_mu=L_OR})}
\item{L_tau }{learning estimate for tau}
\item{S_mu }{sensitivity estimate for the effect mu (\code{S_mu=S_OR})}
\item{S_tau }{sensitivity estimate for tau}
If \code{scale_effect="logOR"},
then the estimates in columns 5-8 and the corresponding column names are adapted as follows
(the remaining columns remain unchanged), where 
mu denotes the effect on the
log-odds ratio or log-odds scale:
\item{median_post_mu }{posterior median for the effect mu}
\item{95CrI_post_mu_low }{lower end point of the 95 \% shortest credible interval (CrI) for the effect mu}
\item{95CrI_post_mu_up }{upper end point of the 95 \% shortest CrI for the effect mu}
\item{length_95CrI_post_mu }{length of the 95 \% shortest CrI for the effect mu}


The following heterogeneity priors are given in the rows (from top to bottom):
  \item{HN_U1tail, EXP_U1tail, HC_U1tail, LMX_U1tail }{half-normal, exponential, half-Caucy and Lomax prior, all \cr \code{tail_alpha_static}-tail adjusted with threshold U=\code{U1} (static)}
 \item{HN_rlmc1, EXP_rlmc1, HC_rlmc1, LMX_rlmc1 }{dynamically adjusted with target RLMC=\code{rlmc1} and tail probability \cr \code{tail_alpha_dynamic}}
   \item{HN_U2tail, EXP_U2tail5, HC_U2tail5, LMX_U2tail }{\code{tail_alpha_static}-tail adjusted with threshold U=\code{U2} (static)}
   \item{HN_rlmc2, EXP_rlmc2, HC_rlmc2, LMX_rlmc2 }{dynamically adjusted with target RLMC=\code{rlmc2} and tail probability \cr \code{tail_alpha_dynamic}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. Sensitivity quantification in Bayesian meta-analysis. Manuscript revised for \emph{Research Synthesis Methods}. 2020.

Roever C. Bayesian random-effects meta-analysis using the bayesmeta R package. \emph{Journal of Statistical Software} (accepted), 2018. 

Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. 
\emph{Bayesian Analysis} \bold{10(2)}, 321--349.
\url{https://projecteuclid.org/euclid.ba/1422884977}
}
\note{
For effects which are not on the log-odds (ratio) scale, the prior for the effect
needs to be adjusted, but otherwise the same code can be used
if a Bayesian normal-normal hierarchical model is appropriate.
%%  ~~further notes~~
}
\section{Warning }{
This function takes ca. 5-10 minutes to run
on the acute graft rejection data set given in the example below.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[bayesmeta]{bayesmeta}} in package bayesmeta,
\code{\link{pri_par_adjust_static}}, \cr
\code{\link{pri_par_adjust_dynamic}},
\code{\link{effective_rlmc}},
\code{\link{pri_par_epsilon_grid}},
\code{\link{or2logor}}
}
\examples{
# Acute Graft rejection (AGR) data analyzed in Friede et al. (2017),  
# Sect. 3.2, URL: https://doi.org/10.1002/bimj.201500236
# First study: experimental group: 14 cases out of 61; 
# control group: 15 cases out of 20 
# Second study: experimental group: 4 cases out of 36; 
# control group: 11 cases out of 36 
rT <- c(14,4)
nT <- c(61,36)
rC <- c(15,11)
nC <- c(20,36)
df <- data.frame(y = log((rT*(nC-rC))/(rC*(nT-rT))), # log-OR
                 sigma = sqrt(1/rT+1/(nT-rT)+1/rC+1/(nC-rC)), # SE(log-OR)
                 labels = c(1:2))
  
# compute the table for the AGR data with median (alpha=50\%) tail-
# adjustment with thresholds U1=0.337 and U2=0.674
# and a weighted harmonic mean to compute the reference sd
# warning: it takes ca. 5-10 minutes to run this function
# on the above data set!
sensitivity_learning_table_flexible(df, tail_alpha_static=0.5,
                                    U1=0.337, U2=0.674,
                                    type_sigma_ref="harmonic")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")

