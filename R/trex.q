#$Author: sinnwell $
#$Date: 2011/01/25 21:31:09 $
#$Header: /projects/genetics/cvs/cvsroot/trex/R/trex.q,v 1.1.1.1 2011/01/25 21:31:09 sinnwell Exp $
#$Locker:  $
#$Log: trex.q,v $
#Revision 1.1.1.1  2011/01/25 21:31:09  sinnwell
#initial for trex package
#
#Revision 1.3  2010/03/24 15:02:24  sinnwell
#remove chiSign and fisherSign from result
#
#Revision 1.2  2010/03/17 20:40:17  sinnwell
#update with more output from C, and for print method
#
#Revision 1.1  2010/03/10 20:25:40  sinnwell
#Initial revision
#

##########################################################
## Purpose: Calculate truncated exact test p-values for two-stage
##         case-control designs of sampling only cases in the 1st stage,
##         and cases+controls in the 2nd stage if there are
##         greater than threshold variants in case1 sample
##
## Authors:  Dan Schaid and Jason Sinnwell
## Created On: 3/4/2010
##
###########################################################

## Input:
## threshold: number of rare variants present in stage 1
##            to continue to 2nd stage
## tab32: 3x2 table of counts with col-1 being variant+, 2nd column variant-
##      row1 is 1st-stage cases
##      row2 is 2nd-stage cases
##      row3 is 2nd-stage controls

trex <- function(tab32, threshold=2) {

  ## check the 3x2 table for only integers and dims
  if(nrow(tab32) != 3 || ncol(tab32) != 2) {
    stop("Error, 3x2 table expected\n")
  }   
  
  ## make a vector from table by rows, check for integers only
  tablevec <- c(tab32[1,], tab32[2,], tab32[3,])
    
  if(!all(tablevec == as.integer(tablevec)))
    stop("tab32 has non-integer values.")
  
  chistat <- chi1 <-  chi2 <- fish1 <- fish2 <- nullprob <- 0.0
  chiSign <- fisherSign <- 0
  
  ## check the threshold against stage 1 variants
  if(threshold > tab32[1,1])
    stop("Stage 1 does not meet threshold.")

  if(threshold != as.integer(threshold))
    stop("threshold is not an integer.")
  
  save <- .C("trexR",
             threshold=as.integer(threshold),
             table=as.integer(tablevec),
             chistat =as.double(chistat),
             chi2sided =as.double(chi2),
             chi1sided =as.double(chi1),
             chiSign = as.integer(chiSign),
             fisher2sided =as.double(fish2),
             fisher1sided =as.double(fish1),
             fisherSign = as.integer(fisherSign),
             probExcluded=as.double(nullprob),
             PACKAGE="trex")
  
  out <- list(tab32=tab32,
              threshold=threshold,
              probExcluded=save$probExcluded,
              chistat=save$chistat,
              chi2=save$chi2sided,
              chi1=save$chi1sided,
              fisher2=save$fisher2sided,
              fisher1=save$fisher1sided)

  class(out) <- "trex"
  return(out)
  
}

