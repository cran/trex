#$Author: sinnwell $
#$Date: 2010/03/23 13:12:50 $
#$Header: /people/biostat3/sinnwell/Projects/TwoStage/Make/RCS/print.trex.q,v 1.2 2010/03/23 13:12:50 sinnwell Exp $
#$Locker:  $
#$Log: print.trex.q,v $
#Revision 1.2  2010/03/23 13:12:50  sinnwell
#changes suggested by Dan, to match stand-alone output
#
#Revision 1.1  2010/03/17 20:51:06  sinnwell
#Initial revision
#

## Purpose: Print a trex (truncated exact test) object
## Author: Jason Sinnwell and Dan Schaid

print.trex <- function(x, ...) { 

  ## print a trex object, similar to format in stand-alone C code

  df32 <- as.data.frame(x$tab32, row.names=c("case1", "case2", "control2"))
  df22 <- as.data.frame(rbind(x$tab32[1,]+x$tab32[2,], x$tab32[3,]), row.names=c("case", "control"))
  names(df32) <- names(df22) <- c("positive", "negative")

  cat("\nThreshold = ", x$threshold, "\n")

  cat("\nObserved 3x2 table:\n")
  print(df32)
  
  cat("\nObserved 2x2 table:\n")
  print(df22)

  cat("\nNull probability of excluded tables = ", x$probExcluded, "\n\n")

  cat("Chi-square statistic:\n")
  cat("  chi-square = ", x$chistat, "\n") 
  cat("pval 2-sided = ", x$chi2, "\n")
  cat("pval 1-sided = ", x$chi1, "\n\n")

  cat("Fisher's p-values:\n")
  cat("pval 2-sided = ", x$fisher2, "\n")
  cat("pval 1-sided = ", x$fisher1, "\n\n")

  invisible(x)

}

