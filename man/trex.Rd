\name{trex}
\alias{trex}
\alias{print.trex}
\title{
Truncated Exact Test
}
\description{
Truncated exact test for two-stage case-control designs for evaluating 
rare genetic variants
}
\usage{
trex(tab32, threshold=2)
\method{print}{trex}(x, ...)

}
\arguments{
\item{tab32}{
A matrix of integer counts of cases and controls with 3 rows and 2 columns. 
Rows are: stage-1 cases, stage-2 cases, and stage-2 controls.
Columns are carrier status of any rare genetic variants (yes in first
column, no in second column)
}
\item{threshold}{
  Total number of cases at stage-1 that carry any rare variants required
  to continue to stage-2.
}
\item{x}{
An object of trex class, returned from the trex() function.
}
\item{...}{
Miscellaneous parameters for the print method
}
}
\value{
A trex class object with components:
\item{tab32}{
The 3x2 matrix passed to trex
}
\item{threshold}{
The treshold value passed to trex
}
\item{probExcluded}{
Null probability of excluded 3x2 tables
}
\item{chistat}{
Chi-square test statistic
}
\item{chi2}{
Two-sided chi-square test exact p-value
}
\item{chi1}{
One-sided chi-square test exact p-value
}
\item{fisher2}{
Two-sided Fisher's Exact test p-value
}
\item{fisher1}{
One-sided Fisher's Exact test p-value
}
}
\details{
The two-stage design uses diseased cases to screen for rare variants at
stage-1, then compares the frequency of variants among all cases versus
controls at stage-2 by an exact test that corrects for the stage-1
ascertainment.
}
\source{
Schaid DJ, Sinnwell JP, "Two-Stage Case-Control Designs for Rare Genetic 
Variants".  Hum Genetics. Published online, 30 Mar 2010.
}
\examples{

## Example 1
tab11 <- cbind(c(2,1,8), c(46,46,98))

out11 <- trex(tab11, threshold=2)
print.trex(out11)

## Example 2
tab10 <- cbind(c(8, 0, 2), c(92, 0, 98))
out10 <- trex(tab10, threshold=2)

print(out10)

out10.t4 <- trex(tab10, threshold=4)

out10.t4

}

