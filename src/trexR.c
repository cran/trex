/* $Author: sinnwell $ */
/* $Date: 2010/03/23 13:37:32 $ */
/* $Header: /people/biostat3/sinnwell/Projects/TwoStage/Make/RCS/trexR.c,v 1.2 2010/03/23 13:37:32 sinnwell Exp $ */
/* $Locker:  $ */
/*
 * $Log: trexR.c,v $
 * Revision 1.2  2010/03/23 13:37:32  sinnwell
 * update with errmsg
 *
 * Revision 1.1  2010/03/17 20:44:18  sinnwell
 * Initial revision
 *
 * Revision 1.2  2010/03/11 19:37:53  sinnwell
 * change back to malloc and free(), instead of Calloc and Free
 *
 * Revision 1.1  2010/03/10 20:25:40  sinnwell
 * Initial revision
 * * 
 */
#include "trexdriver.h"

void trexR(int *threshold,  // input, threshold count of rare variants
	   int *tablevec,   // input, vector to fill 3x2 obsTable
	   double *chistatObs,
	   double *chi2sided,
	   double *chi1sided,
	   int *chistatSign,
	   double *fisher2sided,
	   double *fisher1sided,
	   int *fisherSign,
	   double *probExcluded)
{

  // DECLARE OBJECTS FOR INPUT AND OUTPUT FROM PROGRAM
  static int verbose=0;

  // PREPARE DATA FOR TREX DRIVER
  int nrow=3;
  int ncol=2;
  int **obsTable = imatrix(1,nrow, 1,ncol);
  if(!obsTable)
    errmsg("Memory allocation failure for obsTable\n");

  fillTable(tablevec, obsTable, nrow, ncol);


  // CALL TREX DRIVER 
  trexDriver(*threshold, obsTable, 
	     chistatObs, chi2sided, chi1sided, chistatSign, 
             fisher2sided, fisher1sided, fisherSign,
	     probExcluded, verbose);

  // CLEAN MEMORY AND RETURN
  free_imatrix(obsTable, 1, nrow, 1, ncol);
  
  return;

}

/***********************************************************************************/
