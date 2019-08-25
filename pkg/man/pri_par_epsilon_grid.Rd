\name{pri_par_epsilon_grid}
\alias{pri_par_epsilon_grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grid computation for epsilon-local sensitivity
%%  ~~function to do ... ~~
}
\description{
Starting with a base prior from a one-parameter family, this function determines the parameters of two
priors from the same parametric family which have the given Hellinger distance to the base prior.
The available one-parameter distributions are half-normal, exponential, half-Cauchy and Lomax
with shape parameter fixed at 1.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pri_par_epsilon_grid(AA0_HN, AA0_EXP, AA0_HC, AA0_LMX, 
                     grid_epsilon=0.00354)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{AA0_HN}{
 scale parameter of the base half-normal prior
}
\item{AA0_EXP}{
  scale parameter of the base exponential prior
}
\item{AA0_HC}{
  scale parameter of the base half-Cauchy prior
}
\item{AA0_LMX}{
  scale parameter of the base Lomax prior with shape parameter=1
}
\item{grid_epsilon}{
  value for epsilon, the Hellinger distance between the base prior
  and the \cr
  "lower"/"upper" prior on the grid
}
}
\details{
This 1-dimensional epsilon grid is needed for computation of the
epsilon-local sensitivity in the function \code{sensitivity_learning_table}.
See Ott et al. (2019) for a brief description of this methodology and Roos et al. (2015) for more details.
The default value for \code{grid_epsilon} corresponds to the 
Hellinger distance between two normal distribution with unit variance and
a difference in means of 0.01,
see Roos et al. (2015, Sect. 2.2) for calibration and interpretation of Hellinger distance values.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of the following eight scale parameter values:
\item{p_HN_l }{parameter of the "lower" half-normal prior on the grid}
\item{p_HN_u }{parameter of the "upper" half-normal prior on the grid}
\item{p_EXP_l }{parameter of the "lower" exponential prior on the grid}
\item{p_EXP_u }{parameter of the "upper" exponential prior on the grid}
\item{p_HC_l }{parameter of the "lower" half-Cauchy prior on the grid}
\item{p_HC_u }{parameter of the "upper" half-Cauchy prior on the grid}
\item{p_LMX_l }{scale parameter of the "lower" Lomax prior (with shape parameter=1) on the grid}
\item{p_LMX_u }{scale parameter of the "upper" Lomax prior (with shape parameter=1) on the grid}
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. The relative latent model complexity adjustment for heterogeneity prior specification in Bayesian meta-analysis. Manuscript submitted to \emph{Research Synthesis Methods}. 2019.

Roos, M., Martins, T., Held, L., Rue, H. (2015). Sensitivity analysis for Bayesian hierarchical models. 
\emph{Bayesian Analysis} \bold{10(2)}, 321--349.
\url{https://projecteuclid.org/euclid.ba/1422884977}
%% ~put references to the literature/web site here ~
}

\seealso{
\code{\link{sensitivity_learning_table}},
\code{\link{cal_h_dist}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# epsilon grid for 5 \% tail-adjusted priors with 
# reference threshold UU=1 and epsilon = 0.00354 (default)
pri_par_epsilon_grid(AA0_HN=pri_par_adjust_static(UU=1)$p_HN, 
                     AA0_EXP=pri_par_adjust_static(UU=1)$p_EXP,
                     AA0_HC=pri_par_adjust_static(UU=1)$p_HC,
                     AA0_LMX=pri_par_adjust_static(UU=1)$p_LMX)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

