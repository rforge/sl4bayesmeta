\name{pri_par_adjust_static}
\alias{pri_par_adjust_static}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Static prior tail-adjustment 
%%  ~~function to do ... ~~
}
\description{
Adjusts the scale parameter of the half-normal, exponential, half-Cauchy
and Lomax distributions
such that the probability mass above the reference threshold
equals the given tail probability.
The shape parameter of the Lomax distribution is fixed at 1.
}
\usage{
pri_par_adjust_static(UU=1, tail_prob=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UU}{
  reference threshold
}
\item{tail_prob}{
  probability mass of the prior above the reference threshold
}
}
\details{
The condition imposed on the distribution is
\eqn{Pr[x > UU] = \alpha}{Pr[x > UU] = alpha},
where \eqn{\alpha}{alpha} denotes the tail probability \code{tail_prob}.

% In the literature, it is common to choose \code{tail_prob=0.05}.
To find a prior for the heterogenity standard deviation,
the standard approach for outcomes on the log-odds ratio scale
is to take the tail probability \eqn{\alpha = 5 \%}{alpha = 5 \%}
(\code{tail_prob=0.05}) and
thresholds
\eqn{UU=1} or \eqn{UU=2}
(Spiegelhalter et al., 2004; Neuenschwander et al., 2010).
For the half-normal distribution,
these conditions lead to the HN(0.5) and HN(1) prior,
respectively.
%If the outcome is on the log-odds (ratio) scale,
%Neuenschwander et al. (2010) recommend the thresholds
%\eqn{UU=1} or \eqn{UU=2} for the 5 \% tail adjustment 
%of the half-normal prior
%for the heterogeneity standard deviation.
}
\value{
A list of four scale parameter values (one for each prior considered):
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{p_HN }{parameter of half-normal prior}
  \item{p_EXP }{parameter of exponential prior}
  \item{p_HC }{parameter of half-Cauchy prior}
  \item{p_LMX }{scale parameter for Lomax prior with shape parameter=1}
%% ...
}
\references{
Neuenschwander, B., Capkun-Niggli, G., Branson, M., Spiegelhalter, D. (2010). Summarizing historical information on controls in clinical trials. \emph{Clinical Trials} \bold{7(1)}, 5--18.
\url{https://doi.org/10.1177/1740774509356002}

Spiegelhalter, D., Abrams, K., Myles, J. (2004). Bayesian Approaches to Clinical Trials and Health-Care Evaluation. John Wiley & Sons, Ltd.
%% ~put references to the literature/web site here ~
}


\seealso{
\code{\link{pri_par_adjust_dynamic}} 
%for the alternative dynamic prior adjustment
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# 5 % tail-adjustment
pri_par_adjust_static(UU=1)
pri_par_adjust_static(UU=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
