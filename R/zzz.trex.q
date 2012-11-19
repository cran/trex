#$Author: sinnwell $
#$Date: 2011/01/25 21:31:09 $
#$Header: /projects/genetics/cvs/cvsroot/trex/R/zzz.trex.q,v 1.1.1.1 2011/01/25 21:31:09 sinnwell Exp $
#$Locker:  $
#$Log: zzz.trex.q,v $
#Revision 1.1.1.1  2011/01/25 21:31:09  sinnwell
#initial for trex package
#
#Revision 1.1  2010/03/10 20:25:40  sinnwell
#Initial revision
#
#
#
#

.onLoad <- function(lib, pkg) {
   library.dynam("trex", pkg, lib)
}

.onAttach <- function(lib, pkg) {
   library.dynam("trex", pkg, lib)
}

##.First.lib <- function(lib, pkg) {
##
##library.dynam("trex", pkg, lib)
##}


.Last.lib <- function(libpath) {
  library.dynam.unload("trex", libpath)
}

.onUnload <- function(libpath) {
  library.dynam.unload("trex", libpath)
}
