\name{sl4bayesmeta-package}
\alias{sl4bayesmeta-package}
\docType{package}
\title{
\packageTitle{sl4bayesmeta}
}
\description{
The main function \code{\link{sensitivity_learning_table}}() provides posterior, sensitivity and learning estimates
for the Bayesian normal-normal hierarchical model used for Bayesian meta-analysis
under four different heterogeneity priors (half-normal, half-Cauchy, exponential, Lomax).
The more advanced function \code{\link{sensitivity_learning_table_flexible}}()
enables a flexible choice of several parameters
and supports computation of the reference within-study standard deviation
based on both the geometric mean and a weighted harmonic mean.

In order to unify notation, the heterogeneity priors are defined
as scaled distributions \eqn{\tau \sim A_0 |X|}{tau ~ A_0 |X|},
where \eqn{A_0} is a scale parameter and
\eqn{X} is the standard form of the distribution.
The methodology implemented is proposed in Ott et al. (2020).
The function \code{\link{pri_par_adjust_dynamic}}() implements the novel heterogeneity prior adjustment
with respect to the relative latent model complexity (RLMC).
The function \code{\link{pri_par_adjust_static}}() implements the standard approach to
heterogeneity prior tail adjustment (Spiegelhalter et al. 2004).
}
\details{
Package: sl4bayesmeta

Type: Package

Title: Sensitivity and learning for Bayesian meta-analysis

Version: 0.3-1

Date: 2020-02-18

Author: Manuela Ott [aut, cre], Malgorzata Roos [aut] 

Maintainer: Manuela Ott <manuela.ott@uzh.ch>

Depends: bayesmeta

License: GPL (>=2)
%The DESCRIPTION file:
%\packageDESCRIPTION{pCalibrate}
%\packageIndices{pCalibrate}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Manuela Ott, Malgorzata Roos
%\packageAuthor{pCalibrate}
Maintainer: \packageMaintainer{sl4bayesmeta}
}
\references{
Ott, M., Hunanyan, S., Held, L., Roos, M. Sensitivity quantification in Bayesian meta-analysis. Manuscript revised for \emph{Research Synthesis Methods}. 2020.

Spiegelhalter, D., Abrams, K., Myles, J. (2004). Bayesian Approaches to Clinical Trials and Health-Care Evaluation. John Wiley & Sons, Ltd.
}
\keyword{ package
}
\keyword{htest}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
# Acute Graft rejection (AGR) data analyzed in Friede et al. (2017), 
# Sect. 3.2, URL: https://doi.org/10.1002/bimj.201500236
df <- data.frame(y = c(-2.310, -1.258), # log-odds-ratio
                  sigma = c(0.599, 0.642), # SE(log-odds-ratio)
                  labels = c(1:2))
  
# compute posterior, sensitivity and learning estimates for AGR data
# warning: it takes ca. 5-10 minutes to run this function 
# on the above data set!
sensitivity_learning_table(df)

# dynamic prior adjustement based on RLMC
pri_par_adjust_dynamic(df=df, rlmc=0.25)
# static 5 % prior tail adjustement with reference threshold UU=1
pri_par_adjust_static(UU=1)
}
