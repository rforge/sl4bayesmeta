\name{sensitivity_learning_table}
\alias{sensitivity_learning_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior, sensitivity and learning estimates for Bayesian meta-analysis
%%  ~~function to do ... ~~
}
\description{
Generates a table containing posterior, sensitivity and learning estimates and
effective RLMC 
values
for a meta-analysis data set with effects expressed as log-odds ratios or log-odds.
%The effect estimates are also given on the log-odds (ratio) scale.
Assumes a Bayesian normal-normal hierarchical model with different 
priors (half-normal, exponential, half-Cauchy and Lomax) for the between-study standard deviation.
These priors are either 5 \% tail-adjusted or adjusted to a target RLMC value.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sensitivity_learning_table(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  a data frame containing a column "y" with the estimates of the log-odds (ratios) for the individual studies, 
  a column "sigma" with the corresponding standard errors and a column "labels" with labels for the studies
%%     ~~Describe \code{x} here~~
}
}
\details{
The sensitivity measure used is epsilon-local sensitivity (Roos et al. 2015), which is based on the
Hellinger distance.
Here, learning refers to the ability of the data to modify the prior.
It is quantified by the Hellinger distance between the prior and its marginal posterior.
See Ott et al. (2020) for a description of the implemented methodology and some examples.

The shape parameter of the Lomax prior is fixed at 1.
The reference thresholds used for the 5 \% tail-adjustment are U=1 and 
U=2 and the target RLMC values considered are RLMC=0.25 and RLMC=0.5.
The posterior estimates are obtained from the \code{bayesmeta()} function
in the package \code{bayesmeta}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with the different estimates in the columns and
one row per heterogeneity prior (i.e. the prior for the between-study standard deviation).
The quantities given in the columns are as follows, where 
mu denotes the effect on the log-odds ratio or log-odds scale and 
tau the heterogeneity standard deviation:
\item{U }{reference threshold for prior adjustment}
\item{tail_prob }{tail probability for prior adjustment}
\item{par_val }{scale parameter value of the heterogeneity prior}
\item{MRLMC }{median relative latent model complexity estimated from MC sampling}
\item{median_post_mu }{posterior median for the effect mu}
\item{95CrI_post_mu_low }{lower end point of the 95 \% shortest credible interval (CrI) for the effect mu}
\item{95CrI_post_mu_up }{upper end point of the 95 \% shortest CrI for the effect mu}
\item{length_95CrI_post_mu }{length of the 95 \% shortest CrI for the effect mu}
\item{median_post_tau }{posterior median for tau}
\item{95CrI_post_tau_low}{lower end point of the 95 \% shortest CrI for tau}
\item{95CrI_post_tau_up }{upper end point of the 95 \% shortest CrI for tau}
\item{length_95CrI_post_tau }{length of the 95 \% shortest CrI for tau}
\item{L_mu }{learning estimate for the effect mu}
\item{L_tau }{learning estimate for tau}
\item{S_mu }{sensitivity estimate for the effect mu}
\item{S_tau }{sensitivity estimate for tau}

The following heterogeneity priors are given in the rows (from top to bottom):
  \item{HN(0.5)_U1tail5perz, EXP_U1tail5perz, HC_U1tail5perz, LMX_U1tail5perz }{half-normal, exponential, half-Caucy and Lomax prior, all 5 \%-tail adjusted with threshold U=1 (static)}
 \item{HN_mrlmc025, EXP_mrlmc025, HC_mrlmc025, LMX_mrlmc025 }{median-adjusted with target RLMC=0.25 (dynamic)}
   \item{HN(1)_U2tail5perz, EXP_U2tail5perz, HC_U2tail5perz, LMX_U2tail5perz }{5 \%-tail adjusted with threshold U=2 (static)}
   \item{HN_mrlmc050, EXP_mrlmc050, HC_mrlmc050, LMX_mrlmc050 }{median-adjusted with target RLMC=0.50 (dynamic)}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. Sensitivity quantification in Bayesian meta-analysis. Manuscript revised for \emph{Research Synthesis Methods}. 2020.

Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. 
\emph{Bayesian Analysis} \bold{10(2)}, 321--349.
\url{https://projecteuclid.org/euclid.ba/1422884977}
}
\note{
This function covers the 5\%-tail and the RLMC-based adjustment
of heterogeneity priors, but not
 the 50\%-tail adjustment (which is not based on RLMC) considered in
Ott et al. (2020).
That adjustment can be studied by using the function
\code{sensitivity_learning_table_flexible}.
%It focusses on the 5\%-tail and the RLMC-based adjustment.

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
\code{\link{effective_rlmc}}, \cr
\code{\link{sensitivity_learning_table_flexible}}
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
  
# compute the table for the AGR data
# warning: it takes ca. 5-10 minutes to run this function
# on the above data set!
sensitivity_learning_table(df)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")

