%- $Author: sinnwell $
%- $Date: 2011/01/25 21:31:09 $
%- $Header: /projects/genetics/cvs/cvsroot/trex/man/optimalDesign.Rd,v 1.1.1.1 2011/01/25 21:31:09 sinnwell Exp $
%- $Locker:  $
% 
% $Log: optimalDesign.Rd,v $
% Revision 1.1.1.1  2011/01/25 21:31:09  sinnwell
% initial for trex package
%
% Revision 1.7  2010/04/08 20:06:27  sinnwell
% rm keyword
%
% Revision 1.6  2010/04/08 18:51:30  sinnwell
% *** empty log message ***
%
% Revision 1.5  2010/04/06 15:30:32  sinnwell
% *** empty log message ***
%
% Revision 1.4  2010/04/05 14:30:09  sinnwell
% fix value
%
% Revision 1.1  2010/04/02 21:38:06  sinnwell
% Initial revision
%
%

\name{optimalDesign}
\alias{optimalDesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal design for a two-stage case-control study for evaluating rare genetic variants
}
\description{
  Find a series of potential designs, including and an optimal design
  for a two-stage case-control study for evaluating rare genetic variants
}
\usage{
optimalDesign(probCarrier, or, alpha, power, maxThreshold = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probCarrier}{
Probability that a person is a carrier of any of the rare genetic variants
}
  \item{or}{
Odds ratio of disease for carriers of any of the rare genetic variants
}
  \item{alpha}{
Type-I Error
}
  \item{power}{
Power of the two-stage design
}
  \item{maxThreshold}{
Maximum allowable threshold over all potential designs
}
}
\value{
A list, which contains the input parameters and a data.frame with the following columns:
 \item{Ntotal}{total number of cases and controls, assuming Ntotal/2 = Ncases=Ncontrols}
 \item{Power.cond}{conditional power, conditional on continuing to stage-2}
 \item{Ncase1}{number of cases for stage-1}
 \item{Threshold}{continuation threshold at stage-1 (Ncarriers in stage-1 >= threshold)}
 \item{ProbStopNull}{probability of stopping at stage-1 under the null (OR=1)}
 \item{ProbStopAlt}{probability of stopping at stage-1 under alternative (OR = input or)}
 \item{ESN.Null}{expected sample size under the null}
 \item{ESN.Alt}{expected sample size under the alternative specified OR}
 \item{Optimal}{ TRUE for the design with minimum ESN.Null}
}
\references{
Schaid DJ, Sinnwell JP, "Two-Stage Case-Control Designs for Rare Genetic 
Variants".  Hum Genetics. In Press.
}
\examples{

## Example with pr(carrier)=.05, OR=5, alpha=.05, power=.9
design.pr05 <- optimalDesign(probCarrier=.05, or=5, alpha=.05, power=.9)

design.pr05

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

