#$Author: sinnwell $
#$Date: 2010/03/10 20:25:40 $
#$Header: /people/biostat3/sinnwell/Projects/TwoStage/Make/RCS/zzz.trex.q,v 1.1 2010/03/10 20:25:40 sinnwell Exp $
#$Locker:  $
#$Log: zzz.trex.q,v $
#Revision 1.1  2010/03/10 20:25:40  sinnwell
#Initial revision
#
#
#
#

.First.lib <- function(lib, pkg) {

   if(is.R()) library.dynam("trex", pkg, lib)

   
 }
