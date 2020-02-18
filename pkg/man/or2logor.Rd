\name{or2logor}
\alias{or2logor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform the effect estimates from the odds (ratio) to the
log-odds (ratio) scale
%%  ~~function to do ... ~~
}
\description{
Transforms the posterior, sensitivity and learning estimates 
for the odds (ratio) in the given table to the log-odds (ratio) scale.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
or2logor(tab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tab}{
  % TODO: make more general?
  A matrix with 16 rows and 16 columns,
  which has the same structure as the output of the
  \code{sensitivity_learning_table_flexibe} function (including the column names)
%%     ~~Describe \code{x} here~~
}
}
\details{
  Since the sensitivity and learning measures are invariant
with respect to monotone transformations of the parameter,
the sensitvity and learning estimates for the effect do not change
under the log transformation.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
% TODO: adapt
A matrix with 16 rows and 16 columns,
which is identical to the input matrix \code{tab},
except that columns 5-8 have been recalculated and renamed as follows,
where mu=log(OR) denotes the effect on the log-odds (ratio) scale:
\item{median_post_mu }{posterior median for the effect mu}
\item{95CrI_post_mu_low }{lower end point of the 95 \% shortest credible interval (CrI) for the effect mu}
\item{95CrI_post_mu_up }{upper end point of the 95 \% shortest CrI for the effect mu}
\item{length_95CrI_post_mu }{length of the 95 \% shortest CrI for the effect mu}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
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
tab.OR <- sensitivity_learning_table_flexible(df)
tab.logOR <- or2logor(tab.OR)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ htest }% use one of  RShowDoc("KEYWORDS")

