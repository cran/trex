%- $Author: sinnwell $
%- $Date: 2010/04/05 13:25:06 $
%- $Header: /people/biostat3/sinnwell/Projects/TwoStage/Make/RCS/trex-internal.Rd,v 1.2 2010/04/05 13:25:06 sinnwell Exp $
%- $Locker:  $
% 
% $Log: trex-internal.Rd,v $
% Revision 1.2  2010/04/05 13:25:06  sinnwell
% *** empty log message ***
%
% Revision 1.1  2010/04/02 21:37:39  sinnwell
% Initial revision
%
% 
\name{trex-internal}
\title{Internal trex functions}
\alias{dbin}
\alias{probStop }
\alias{powerCaseControl}
\alias{opt.min.esn}
\alias{design.stage1}
\alias{possibleDesign}
\alias{bisect2}
\alias{powBisect2.func}
\alias{findN}
\description{Internal trex functions}
\usage{
 findN(pCarrier, or, alpha, pow.conditional)
 bisect2(fn, pars, xL, xU)
 possibleDesign(pCarrier, or, nTotal, pow.min, alpha=.05, dat=NULL)
}
\author{Dan Schaid, Jason Sinnwell}
\details{These are not intended for use by users}
\keyword{internal}
