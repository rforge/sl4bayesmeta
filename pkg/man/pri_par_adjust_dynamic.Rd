\name{pri_par_adjust_dynamic}
\alias{pri_par_adjust_dynamic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic heterogeneity prior adjustment based on the relative latent model complexity (RLMC)
%%  ~~function to do ... ~~
}
\description{
Adjusts the scale parameter of 
the half-normal, exponential, half-Cauchy and Lomax prior
for the between-studies standard deviation
such that the probability mass above the reference threshold
equals the given tail probability.
The reference threshold used depends on the target RLMC and the with-study
standard errors in the data set.
The shape parameter of the Lomax distribution is fixed at 1.
}
\usage{
pri_par_adjust_dynamic(df, rlmc=0.5, tail_prob=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "sigma" containing the standard errors of the estimates for the individual studies
}
\item{rlmc}{
  target relative latent model complexity. Real number in [0,1].
}
\item{tail_prob}{
  probability mass of the prior above the reference threshold (which depends on \code{rlmc} and \code{df$sigma})
}
}
\details{
The reference threshold U is given by
\deqn{U = \sigma_{ref} \sqrt{rlmc/(1-rlmc)},}
where \eqn{\sigma_{ref}} is the reference standard deviation of the data set,
i.e. the geometric mean of \code{df$sigma}.
Then, the static prior tail-adjustment is applied for this reference threshold
and the specified tail probability.
This prior adjustment applies to Bayesian meta-analysis
expressed by a normal-normal hierarchical model.

Ott et al. (2019) suggest to use \code{tail_prob=0.5} as default,
so that the medians of the priors will be aligned with the reference threshold.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of four scale parameter values (one for each prior considered):
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{p_HN }{parameter of half-normal prior}
  \item{p_EXP }{parameter of exponential prior}
  \item{p_HC }{parameter of half-Cauchy prior}
  \item{p_LMX }{scale parameter for Lomax prior with shape parameter=1}
%  \item{p_HN }{scale parameter of half-normal prior}
%  \item{p_EXP }{scale parameter of exponential prior}
%  \item{p_HC }{scale parameter of half-Cauchy prior}
%  \item{p_LMX }{scale parameter of Lomax prior with shape parameter=1}
%% ...
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. The relative latent model complexity adjustment for heterogeneity prior specification in Bayesian meta-analysis. Manuscript submitted to \emph{Research Synthesis Methods}. 2019.
}

\seealso{
\code{\link{pri_par_adjust_static}} 
%for the alternative static tail adjustment of priors
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Acute Graft rejection data analyzed in Friede et al. (2017), Sect. 3.2, 
# URL: https://doi.org/10.1002/bimj.201500236
df <- data.frame(y = c(-2.310, -1.258), # log-odds-ratio
                 sigma = c(0.599, 0.642), # SE(log-odds-ratio)
                 labels = c(1:2))
                  
pri_par_adjust_dynamic(df=df, rlmc=0.25)
pri_par_adjust_dynamic(df=df, rlmc=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
