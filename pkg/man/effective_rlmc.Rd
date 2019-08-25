\name{effective_rlmc}
\alias{effective_rlmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effective relative latent model complexity computation
%%  ~~function to do ... ~~
}
\description{
Generates a Monte Carlo sample of relative latent model complexity (RLMC) values
for the given meta-analysis data set and the prior distribution for the
between-studies standard deviation.
Can also return estimates for plotting a smoothed histogram of effective RLMC values.
Assumes a Bayesian normal-normal hierarchical model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
effective_rlmc(df, r.tau.prior, MM=10^6, output="sample", 
               step=0.03)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
  data frame with one column "sigma" containing the standard errors of the estimates for the individual studies
}
\item{r.tau.prior}{
  randomisation function for the prior, \cr
  e.g. \code{function(n)rhalfnormal(n=n, scale=1)}
  for the HN prior with scale parameter 1
}
\item{MM}{
  number of MC samples. Defaults to 10^6.
}
\item{output}{
  either "sample" or "summary" or "prob". Specifies the type of output.
  "sample" returns the MC samples of RLMC values, "summary" gives summary statistics
  of the MC samples and "prob" provides x- and y-values for plotting a smoothed histogram
}
\item{step}{
  bin width for \code{output="prob"}
}
}
\details{
For the definition of RLMC and further information, see Ott et al. (2019).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If \code{output=="sample"}: a vector of size \code{MM} of RLMC values.

If \code{output=="summary"}: summary statistics for the MC samples (minimum, first quantile, median, mean, third quantile and maximum).

If \code{output=="prob"}: a data frame with columns "x" and "y", giving the values in the middle of the bins ("x") and the corresponding "density estimates" ("y").
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. The relative latent model complexity adjustment for heterogeneity prior specification in Bayesian meta-analysis. Manuscript submitted to \emph{Research Synthesis Methods}. 2019.
%% ~put references to the literature/web site here ~
}

\examples{
# Acute Graft rejection data analyzed in Friede et al. (2017), Sect. 3.2, 
# URL: https://doi.org/10.1002/bimj.201500236
df <- data.frame(y = c(-2.310, -1.258), # log-odds-ratio
                  sigma = c(0.599, 0.642), # SE(log-odds-ratio)
                  labels = c(1:2))

# MC sample of RLMC values
effective_rlmc(df=df, r.tau.prior=function(n)rhalfnormal(n=n, scale=1))  
# summary statistics for MC sample of RLMC values
effective_rlmc(df=df, r.tau.prior=function(n)rhalfnormal(n=n, scale=1), 
               output="summary")
# values for smoothed histogram of RLMC values
effective_rlmc(df=df, r.tau.prior=function(n)rhalfnormal(n=n, scale=1), 
               output="prob")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")

