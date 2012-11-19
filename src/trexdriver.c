/* $Author: sinnwell $ */
/* $Date: 2011/01/25 21:31:09 $ */
/* $Header: /people/biostat3/sinnwell/Projects/TwoStage/Make/RCS/trexdriver.c,v 
1.1 2010/03/17 20:44:18 sinnwell Exp sinnwell $ */
/* $Locker:  $ */
/*
 * $Log: trexdriver.c,v $
 * Revision 1.1.1.1  2011/01/25 21:31:09  sinnwell
 * initial for trex package
 *
 * Revision 1.4  2010/03/31 18:22:12  sinnwell
 * declare flogm and factlm before any goto m600
 *
 * Revision 1.3  2010/03/29 21:02:47  sinnwell
 * move some functions to trexC.cpp
 *
 * Revision 1.1  2010/03/17 20:44:18  sinnwell
 * Initial revision
 * * 
 */
/*  License: 
# 
# Copyright 2003 Mayo Foundation for Medical Education and Research. 
# 
# This program is free software; you can redistribute it and/or modify it under the terms of 
# the GNU General Public License as published by the Free Software Foundation; either 
# version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
# more details.
# 
# You should have received a copy of the GNU General Public License along with this 
# program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
# Boston, MA 02111-1307 USA
# 
# For other licensing arrangements, please contact Daniel J. Schaid.
# Daniel J. Schaid, Ph.D.
# Division of Biostatistics
# Harwick Building  Room 775
# Mayo Clinic
# 200 First St., SW
# Rochester, MN 55905
# 
# phone: 507-284-0639
# fax:      507-284-9542
# email: schaid@mayo.edu
# 

*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "trexdriver.h"


// 
// Authors: Dan Schaid and Jason Sinnwell
// Purpose: Driver function and supporting functsion for trex source code.
// Date Created: 3/2010
//


void trexDriver(int threshold, 
		int **obsTable,
		double *chistatObs,
		double *chi2sided,
		double *chi1sided,
		int *chiSign,
		double *fisher2sided,
		double *fisher1sided,
		int *fisherSign,
		double *probExcluded, 
		int verbose)
{

  int r1, r2, r3, r4;
  int c1, c2, c3, c4;
  int ifault = 0;
  int ntables, ntablesExcluded;
  int nenum;
  int ngroups;
  int i,j,k;
  int sign;
  int indexObs;

  int nrow = 3;
  int ncol= 2;

  double probIncluded;
  double chistat;
  
  /* Compute marginal row totals */

  int *rowTot = ivector(1,nrow);
  if(!rowTot)
    errmsg("Memory allocation failure for rowTot\n");

  //    errmsg("Memory allocation failure for rowTot\n");

  int *colTot = ivector(1,ncol);
  if(!colTot)
    errmsg("Memory allocation failure for colTot\n");

  for(i=1; i<=nrow; i++){
    rowTot[i] = 0  ;
    for(j=1; j<=ncol; j++)
      {
	rowTot[i] += obsTable[i][j];
      }
  }
 
  /* Compute marginal col totals */

  for(j=1; j<=ncol; j++)
    {
      colTot[j] = 0  ;
      for(i=1; i<=nrow; i++)
	{
	  colTot[j] += obsTable[i][j];
	}
    }

  for(j=1; j<=ncol; j++)
    {
      if(colTot[j] <= 0 )
	{
	  REprintf("Error: total col[%d] = %d\n", j, colTot[j]);
	  errmsg("Ending program\n");
	}
    }


  // check if rows 1,3  totals <= 0
  for(i=1; i<=nrow; i++)
    {      
      /* skip over row 2 to check later */
      if(i==2) continue;
     
      if(rowTot[i] <= 0 )
	{
	  REprintf("Error: total row[%d] = %d\n", i, rowTot[i]);	  
          errmsg("Row totals not sufficient\n");
	}
    }
      

  /* Account for row 2 of obsTable having all 0's, implying need to 
     consider 2x2 table instead of 3x2 table */

  if(rowTot[2] == 0)
    {
      nrow = 2;
      rowTot[2] = rowTot[3];
    }
  
  /* count of all potential enumerated tables */
  ntables = countTables(colTot[1], rowTot, nrow,0, 0);
  

  TABLE * tableVec;
  tableVec = (TABLE *)malloc((size_t) ((ntables)*sizeof(TABLE)));

  if(!tableVec)
    errmsg("Memory allocation failure for tableVec\n");


  /**** Enumerate all possible 3x2 tables ********/

  nenum = enumTableAll(nrow, ncol, rowTot, colTot, &ifault, tableVec, ntables);
  if(ifault != 0)
    {
      // REprintf("Error: ifault = %d, nenum = %d, ntables = %d\n", ifault, nenum, ntables);
      // exit(1);
      errmsg("Error: ifault problem in enumTableAll\n");
    }


  if(nenum != ntables){
    // REprintf("Error: no. enumerated tables (%d) != number expected (%d)\n", nenum, ntables);
    // exit(1);
    errmsg("Error, number of enumerated tables not same as expected.\n");
  }
  
 

  /* Create exclude indicator for tables that have first row/col cell count < threshold */

  ntablesExcluded = excludeTables(tableVec, ntables, threshold);

  /* Convert 3x2 table to 2x2, by summing rows 1 & 2 and copying row 2 onto row 2 */

  make22(tableVec, ntables);
 
 
  /* sort 2x2 tables to create groups (equivalance classes) */ 

  qsort(tableVec, ntables, sizeof(TABLE), compareTable22);
 
  /* code groups (indices for equivalence classes) */

  groupCode(tableVec, ntables);
  ngroups = tableVec[ntables-1].group;

  //if(verbose) {  /* Output */
  //  REprintf("Counts of tables\n");
  //  REprintf("  enumerated 3x2: %d\n", ntables);
  //  REprintf("    excluded 3x2: %d\n", ntablesExcluded);
  //  REprintf("    included 3x2: %d\n", ntables - ntablesExcluded);
  //  REprintf("      unique 2x2: %d\n", ngroups);
  
    // REprintf("\n\n=========== After code group ===================\n\n");
    // printTableVec(tableVec, ntables);  
    // REprintf("Number of tables = %d, number of groups for 2x2 tables = %d\n", ntables, ngroups);
  // }

  /* compute probabilities of equivalence classes */
 
  groupProb(tableVec, ntables);

  /* sort in order to get unique "representative" tables (one from each group) to top of vector */ 
  qsort(tableVec, ntables, sizeof(TABLE), compareGroupLast);


  probIncluded = 0.0;
  for(k=0; k<ngroups; k++){
    probIncluded += tableVec[k].groupCumProb;
  }

  *probExcluded = 1.0 - probIncluded;


  /* standardized representative tables to have prob sum = 1, by dividing by
     probIncluded tables */

  for(k=0; k<ngroups; k++){
    tableVec[k].groupCumProb =  tableVec[k].groupCumProb/probIncluded ;
  }
  
  // compute stats for each table 
  for(k=0; k<ngroups; k++){
    chi22Table(tableVec[k].t32, &chistat, &sign);
    tableVec[k].chistat = chistat;
    tableVec[k].chistatSign = sign;
  }

  /* convert observed 3x2 to 2x2 table */

  for(j=1; j<=ncol; j++)
    {
      obsTable[1][j] += obsTable[2][j];
      obsTable[2][j]  = obsTable[3][j];
    }

  indexObs = findObs22Table(obsTable, tableVec, ngroups);

  if(indexObs < 0 || indexObs >= ntables)
    {
      //     REprintf("Error: indexObs (%d) for observed table out of range (0, %d)\n", indexObs, ntables - 1);
      // exit(1);
      errmsg("Error: table index out of range\n");
    }

  if(verbose) {
    REprintf("\nObserved 2x2 table\n");
    for(i=1; i<=2; i++){
      REprintf(" ");
      for(j=1; j<=2; j++){
	REprintf("%d ", tableVec[indexObs].t32[i-1][j-1]);
      }
      REprintf("\n");
    }
  }

  computeChiSqPvals(tableVec, ngroups, indexObs, 
                    chistatObs, chi2sided, chi1sided, chiSign);

  computeFisherPvals(tableVec, ngroups, indexObs, 
		     fisher2sided, fisher1sided, fisherSign);

  // if(verbose) {
  //  REprintf("\n\n===================  Possible 2x2 Tables  ========================\n\n");
  //  print22Tables(tableVec,ngroups);
  // }

  // free allocated objects
  // REprintf("Freeing memory in trexC...\n");
  free_ivector(rowTot, 1, nrow);
  free_ivector(colTot, 1, ncol);

  free((FREE_ARG) tableVec);

  return;

}





/************************ OTHER FUNCTIONS FOR DRIVER *******************/


/***********************************************************************************/
void make22(TABLE * tableVec, int ntable){

  /* make a 2x2 table, assuming 0-offset, by overwriting 3x2 table row 1 with sum of rows 1 and 2; 
     copy row 3 onto row 2  */

  int j,k;
  for(k=0; k<ntable; k++){
    for(j=0; j<2; j++){
      tableVec[k].t32[0][j] += tableVec[k].t32[1][j];
      tableVec[k].t32[1][j]  = tableVec[k].t32[2][j];

    }
  }
  return;
}


/***********************************************************************************/

int enumTableAll(int nrow,int ncol,int *rowmarg, int *colmarg, int *ifault,
		  TABLE * tableVec, int lengthTableVec){


  /* 
     Saunders, IW. 1984. Enumeration of R x C tables with repeated row totals.
     Applied Statistics. 1984: 33(3):340-352.

	ALGORITHM AS 205 APPL. STATIST. (1984) VOL 33, NO. 3.
       ENUMERATES ALL R*C CONTINGENCY TABLES WITH GIVEN ROW TOTALS N(I)
       AND COLUMN TOTALS M(J) AND CALCULATES THE HYPERGEOMETRIC
       PROBABILITY OF EACH TABLE.

       FOR TABLES HAVING TWO OR MORE ROW SUMS REPEATED, EQUIVALENT
       TABLES DIFFERING ONLY BY A ROW PERMUTATION ARE NOT SEPARATELY
       ENUMERATED. A REPRESENTATIVE OF EACH EQUIVALENCE CLASS IS ENUMERATED
       AND THE MULTIPLICITY OF EACH CLASS CALCULATED.

       FOR EACH TABLE ENUMERATED, SUBROUTINE EVAL IS CALLED TO CARRY OUTC
       CALCULATIONS ON THE TABLE.


  DAN'S NOTES:

  Translated from Fortran 77 to C by Dan Schaid, June, 2009
  
  Used vectors and matrices with 1-offset, intead of traditional C
  0-offset, to make it easier to translate from Fortan to C, keeping
  track of limits of arrays.

  --> rowmarg and colmarg input vectors have index starting at 1

  But, when copy table to tableVec, the table in the strutc of tableVec
  uses 0-offsets.

  IMPORTANT: this version forces all multiplicities = 1, so that
  all tables are created, including all in the same equivalence class.

  Equivalence class tables occur when there are ties in the row margin totals.
  Swapping rows (permuting rows) doesn't change row totals for these tied rows,
  nor col totals. So, the original version creates only one representative table 
  from an equivalence class, its probability, and it multiplicity. The probabilty for the
  entire equivalence class  = multiplicity * prob.

  However, when creating r x 2 tables, and then collapsing to 2x2 table, different tables
  within an equivalence class can give different 2x2 tables, so need to be 
  careful about handling multiplicities when collapsing tables, eventhough
  all tables within an equivalence class have the same probabilities; their
  resulting 2x2 tables will not have equivalent probabilities, because they
  have different values in the 2x2 tables.

  So, this version (enumTableAll) enumerates all tables, so the prob of each table
  is simply the output prob (i.e., multiplicity = 1 for all tables).


  JASON'S NOTES
  
  Key steps to make enumTableAll work with trexC being called from R:

  1) Every exit point with "return" was replaced by setting indexTable to a valid value, 
     then jumping to c600 where memory is free()-ed and indexTable is returned

  2) Every object that is created as an ivector, dvector, imatrix, or dmatrix 
     must have a free_imatrix method at the end.  The constructors and free() methods for those 
     objects are also used in trexC, where Calloc() and Free() are used instead 
     of malloc() and free()

  3) As in trexC, all REprintf statements are within if(DEBUG) {} 
     and all stop-errors are handled by errmsg()
 
  */

  int r = nrow;
  int c = ncol;
  long indexTable = -1;

  int repsr, repsc, multr, multc;
  int ntotal, ntot, max, nmax, rm, cm, left, rowbnd, rowsum;
  int i, j, ii, jj, iip, iim, jp, keep;
  int k, ip, ij, jnext, jm, jjm, jkeep;
  int TRUE = 1;
  int FALSE = 0;
  int ieqim;

  int totTables = 0;

  double prob0;
  double probTable = 0.0;

  int ONE = 1;
  double ZERO = 0.0;
 
   /* max of r & c for allocation of memory, accounting for possible need
    to transpose table */

   int maxrc = (r > c)? r : c;

   double* flogm;
   double* factlm;

   int *z     = ivector(1, maxrc);
   if(!z)
     errmsg("Memory allocation failure for z\n");

   int *rept  = ivector(1, maxrc);
   if(!rept)
     errmsg("Memory allocation failure for rept\n");

   int *reps  = ivector(1, maxrc);
   if(!reps)
     errmsg("Memory allocation failure for reps\n");

   int *mult  = ivector(1, maxrc);
   if(!mult)
     errmsg("Memory allocation failure for mult\n");

   int *n     = ivector(1, maxrc);
   if(!n)
     errmsg("Memory allocation failure for n\n");

   int *m     = ivector(1, maxrc);
   if(!m)
     errmsg("Memory allocation failure for m\n");

   int    **table = imatrix(1,maxrc, 1,maxrc);
   if(!table)
     errmsg("Memory allocation failure for table\n");

   int    **bound = imatrix(1,maxrc, 1,maxrc);
   if(!bound)
     errmsg("Memory allocation failure for bound\n");

   double **prob  = dmatrix(1,maxrc, 1,maxrc);
   if(!prob)
     errmsg("Memory allocation failure for prob\n");

   /* assume input vecs are 1-offset, and copy to 1-offset vecs internally */

   for(i=1; i <= r; i++){
     n[i] = rowmarg[i];
   }
   for(i =1; i <= c; i++){
     m[i] = colmarg[i];
   }
 
   /* ORIGINAL FORTRAN COMMENTS :   LOGICAL REPT(10),REPTC(10),IEQIM
       LOCAL VARIABLES -
         TABLE(I,J)   - (I,J)-TH ENTRY OF CURRENT TABLE
         NTOTAL       - TOTAL OF TABLE ENTRIES
         BOUND(I,J)   - CURRENT UPPER BOUND ON TABLE(I,J) TO SATISFY
                        ROW AND COLUMN TOTALS
         REPT(I)      - LOGICAL = TRUE IF ROW TOTALS N(I), N(I-1) ARE EQUAL
                                = FALSE OTHERWISE
         REPS(I)      - NUMBER OF PREVIOUS ROWS EQUAL TO ROW I
         MULT(I)      - MAXIMUM NUMBER OF EQUIVALENT TABLES
                        GIVEN FIRST I ROWS
         Z(J)         - LOWER BOUND ON SUM OF ENTRIES USED BY ALGORITHM C
         PROB(I,J)    - PARTIAL SUM OF TERMS IN LOG(P)

         MAX  - MAXIMUM DIMENSION OF TABLE ( DS removed this)
         NMAX - MAXIMUM NUMBER OF OBSERVATIONS + 1
   */

   /*  check input */
 
   /* ifault values:
      0 = no fault
      1 = row or col dims out of range (should be bounded by 1)
      2 = row total != col total
      3 = row totals or col totals <= 0
   */

   *ifault = 1;
 
   if(r <=0 || c <=0)
     { 
       //return 0;
       indexTable=0;  // jps added
       goto c600;     // jps added
     }
 
   *ifault = 3;
   /* total = total of cells in table, so is the sum of  row totals */
   ntotal = 0;
   for(i=1; i <= r; i++){
     if(n[i] <= 0){
       indexTable=0;  // jps added
       goto c600;     // jps added
       //return 0;
     }

     ntotal += n[i];
   }


   nmax = ntotal + 1;

   flogm  = dvector(1, nmax);
   factlm = dvector(1, nmax);

   ntot = 0;
   for(j=1; j <= c; j++){
     if(m[j] <= 0){
       indexTable = 0; // jps added
       goto c600;  // jps added
       //return 0;
     }
     ntot += m[j];
   }
   *ifault = 2;

   if(ntot != ntotal){
     indexTable = 0; // jps added
     goto c600; // jps added
     //     return 0;
   }

   *ifault = 0;

   /*INITIALISE FLOGM(K)=LOG(K-1), FACTLM(K)=LOG(K-1 FACTORIAL) */

   flogm[1] = ZERO;
   factlm[1] = ZERO;
   for(k=1;  k<=ntotal; k++){
     flogm[k+1]= log((double)(k));
     factlm[k+1]= factlm[k]+flogm[k+1];
   }

   rm = r - 1;
   cm = c - 1;

   /* Changes from FORTRAN algorithm:
      - No longer SORT ROWS AND COLUMNS INTO ASCENDING ORDER 
      - Force MULTIPLICITIES OF ROWS AND COLUMNS = 1
      - No longer transpose table
   */
 

   multr = ONE;
   repsr = ONE;
   for(i=1; i <= r; i++){
     rept[i] = FALSE;
   }


   mult[1] = multr;
   reps[1] = ONE;

   /*      CONSTANT TERM IN PROBABILITY */

   prob0 = - factlm[ntotal + 1];

   for(i=1; i <= r; i++)
     {
       ii = n[i];
       prob0 = prob0 + factlm[ii + 1];
     }

   for(j = 1; j <= c; j++)
     {
       jj = m[j];
       prob0 = prob0 + factlm[jj + 1];
     }

   /*       CALCULATE BOUNDS  ON ROW 1 */

   for(j = 1; j <= c; j++)
     {
       bound[1][j] = m[j];
     }

   /*       FOR EACH I FIND GREATEST I-TH ROW SATISFYING BOUNDS */

   for(i=1; i <= r; i++)
     {
       if(i != 1) prob0 = prob[i-1][c];
       left = n[i];

       /* ELEMENTS OF ROW I */

       ieqim = rept[i];
       
       for(j=1; j <= cm; j ++)
	 {
	   
	   ij = (bound[i][j] < left)? bound[i][j] : left;

	   table[i][j] = ij;

	   if(j == 1) prob[i][j] = prob0 - factlm[ij + 1];
	   if(j != 1) prob[i][j] = prob[i][j-1] - factlm[ij + 1];
	   left = left - table[i][j];
	   if(i < r) bound[i+1][j] = bound[i][j] - table[i][j];

	   if(left == 0) goto c121;
	   
	   if(ieqim) ieqim = (table[i][j] == table[i-1][j]) ? TRUE : FALSE; 

	 }
       
       table[i][c] = left;
       prob[i][c] = prob[i][cm] - factlm[left + 1];
       if(i < r) bound[i+1][c] = bound[i][c] - left;

       goto c123;

       c121: 
       jp = j + 1;
       
       for(jj = jp; jj <= c; jj++)
	 {
	   table[i][jj] = 0;
	   prob[i][jj] = prob[i][jj - 1];
	   bound[i+1][jj] = bound[i][jj];
	 }
   
       c123: 
       if(i == 1) continue;
       
       mult[i] = mult[i-1];
       reps[i] = ONE;

       if(ieqim == FALSE){
	   continue;
	 }
     
       reps[i] = reps[i-1] + ONE;
       mult[i] = mult[i] / reps[i];

     }

   /*       CALL EVAL FOR TABLE 1
	    CALL EVAL(TABLE,R,C,N,M,PROB(R,C),MULT(R))
   */

   /* copy table to tableVec */

   indexTable ++;
   
   if(indexTable > lengthTableVec){
     errmsg("Error: index for tableVec > length tableVec\n");
     //     exit(1);
   }

   probTable = exp(prob[r][c]);
   totTables += mult[r];
   tableVec[indexTable].index = indexTable + 1;
   tableVec[indexTable].prob = probTable;

   if(nrow == 2)
     {
       for(j=1; j<=2; j++)
	 {
	   tableVec[indexTable].t32[0][j-1] = table[1][j];
	   tableVec[indexTable].t32[1][j-1] = 0 ;
	   tableVec[indexTable].t32[2][j-1] = table[2][j];

	 }
     } 
   else
     {
       for(i=1; i<=3; i++){
	 for(j=1; j<=2; j++){
	   tableVec[indexTable].t32[i-1][j-1] = table[i][j];
	 }
       }
     }

   /*  COMMENCE ENUMERATION OF REMAINING TABLES */

   /******************       START OF MAIN LOOP ****************/

   c200:  i = r;
   c210:  i = i-1;

   /*       IF I = 0 NO MORE TABLES ARE POSSIBLE */
   
   if(i == 0) {
     indexTable +=1;  //jps added
     goto c600;  // jps added
     //return indexTable + 1;
   }

   j = cm;
   left = table[i][c];
   rowbnd = bound[i][c];

   /*       TRY TO DECREASE ELEMENT (I,J) */


   c220: if(table[i][j] > 0 && left < rowbnd) goto c230;

   /*       ELEMENT (I,J) CANNOT BE DECREASED - TRY (I,J-1) */

   if(j == 1) goto c210;

   left = left + table[i][j];
   rowbnd = rowbnd + bound[i][j];
   j = j - 1; 

   goto c220;

   /*       DECREASE ELEMENT (I,J) */

   c230: 
   ij = table [i][j];
   prob[i][j] = prob[i][j] + flogm[ij + 1];
   table[i][j] = table[i][j] - 1;
   bound[i+1][j] = bound[i+1][j] + 1;

   /*     IF ROW I WAS THE SAME AS ROW I-1 IT IS NO LONGER */

   if(reps[i] == ONE) goto c270;

   reps[i] = ONE;
   mult[i] = mult[i-1];

   /*       COMPLETE ROW I WITH THE LARGEST POSSIBLE VALUES */

   c270:   
   ii = i;
   iip = ii + 1;
   iim = ii - 1;
   jnext = j + 1;
   left = left + 1;
   
   goto c380;

   /*      FILL UP REMAINING ROWS */

   c300: 
   ii = ii + 1;

   /*       THE LAST ROW IS TREATED SEPARATELY */

   if(ii == r) goto c400;
   iip = ii + 1;
   iim = ii - 1;

   /* if(rept[ii]) goto c310;*/
   /* never TRUE, so commented out here and  below at c310 */


   /*       ROW TOTAL N(II) IS NOT A REPEAT - MAKE ROW II AS LARGE AS POSSIBLE */

   left = n[ii];
   jnext = 1;
   goto c380;
  

   c330: 
   if(j > 1) goto c331;

   /*      IF J=1 THE BOUNDS ARE SATISFIED AUTOMATICALLY */

   ij = bound[ii][1];
   table[ii][1] = ij; 
   prob[ii][1] = prob[iim][c] - factlm[ij + 1];
   jnext = 2;
   left = n[ii] - table[ii][1];
   bound[iip][1] = 0;
   
   goto c380;


   c331: 
   z[j] = n[ii]; 
   jm = j - 1;
   if(j == c) goto c336;
   jp = j + 1; 

   for(jj = jp; jj <= c; jj++)
     {
       z[j] = z[j] - bound[ii][jj];
     }


   c336:
   for(jjm = 1; jjm <= jm; jjm++){
     jj = j - jjm;
     z[jj] = z[jj + 1] - bound[ii][jj + 1];
   }

   /*     (II) IF THE CUMULATIVE TOTALS OF ROW II-1 ALL EXCEED THE BOUNDS Z(J)
	  MAKE ELEMENT (II,J) EQUAL TO ITS BOUND */

   rowsum = 0;
   jkeep = 0;

   for(jj = 1; jj <= jm; jj++)
     {
       rowsum=rowsum+table[iim][jj];
       if(rowsum < z[jj]) goto c360;
       if(rowsum > z[jj] && table[iim][jj] > 0) jkeep = jj; 
     }

   table[ii][j] = bound[ii][j];

   bound[iip][j] = 0;
   ij = table[ii][j]; 
   prob[ii][j] = prob[ii][jm] - factlm[ij + 1];
   reps[ii] = ONE; 
   mult[ii] = mult[iim]; 

   /*       COMPLETE ROW II WITH THE LARGEST POSSIBLE ELEMENTS */
   jnext = jp;
   left = n[ii];

   for(jj = 1; jj <= j; jj++)
     {
       left = left - table[ii][jj];
     }

   goto c380;

   /*
       (III) THE CUMULATIVE SUMS VIOLATE THE BOUNDS
       IF NO ELEMENT OF ROW II-1 CAN BE CHANGED TO SATISFY THE BOUNDS
                 NO SUITABLE ROW II IS POSSIBLE
       IN THAT CASE GO BACK AND TRY DECREASING ROW II-1
   */

   c360: 
   if (jkeep != 0) goto c370; 

   i = ii;
   goto c210;
 
   /*      ELEMENT (II,JKEEP) CAN BE DECREASED */

   c370:  
   bound[iip][jkeep] = bound[iip][jkeep] + 1;
   
   ij = table[ii][jkeep];
   prob[ii][jkeep] = prob[ii][jkeep] + flogm[ij + 1];
   table[ii][jkeep] = table[ii][jkeep] - 1;

   /*       COMPLETE THE ROW */

   jnext = jkeep + 1;
   left = n[ii];
   for(jj=1; jj <= jkeep; jj++)
     {
       left = left - table[ii][jj];
     }


   /*      ROW II IS COMPLETE UP TO ELEMENT JNEXT-1
	   MAKE THE REMAINING ELEMENTS AS LARGE AS POSSIBLE
	   (THIS SECTION OF CODE IS USED FOR EVERY ROW, REPEATED OR NOT)
   */


   c380: 
   if(jnext == c) goto c390;

   for(j=jnext; j <= cm; j++)
     {
       table[ii][j] = (left < bound[ii][j]) ?  left : bound[ii][j];
       left = left - table[ii][j];
       bound[iip][j] = bound[ii][j] - table[ii][j];
       ij = table[ii][j];
       if(j == 1) prob[ii][j] = prob[iim][c]  - factlm[ij + 1];
       if(j != 1) prob[ii][j] = prob[ii][j-1] - factlm[ij + 1];
       if(left == 0) goto c391;
     }

   c390: 
   table[ii][c] = left;
   prob[ii][c] = prob[ii][cm] - factlm[left + 1];
   bound[iip][c] = bound[ii][c] - left;

   goto c393;

   c391: 
   jp = j + 1; 

   for(jj = jp; jj <= c; jj++)
     {
       table[ii][jj] = 0;
       prob[ii][jj] = prob[ii][jj - 1];
       bound[iip][jj] = bound[ii][jj];
     }

   c393: 
   reps[ii] = ONE; 

   if(ii > 1) mult[ii] = mult[iim];
   goto c300;

   /*       THE FINAL ROW */

   c400: 

   /* if(rept[r]) goto c420; */
   /* never true, so commented out here and below at c420 */

   /*       NOT A REPEAT - SET ROW R EQUAL TO ITS BOUNDS */


   ij = bound[r][1];
   table[r][1] = ij;
   prob[r][1] = prob[rm][c] - factlm[ij + 1]; 


   for(j=2; j <= c; j++)
     {
       ij = bound[r][j];
       table[r][j] = ij;
       prob[r][j] = prob[r][j - 1] - factlm[ij + 1];
     }

   mult[r] = mult[rm];
   goto c500;

   /*       ROW TOTAL R IS A REPEAT - ENSURE THAT IT IS LESS THAN ROW R-1 */

   /*
     c420:
  
     for(j=1; j <= c; j++)
     {
     if(bound[r][j] > table[rm][j]) goto c440;
     ij = bound[r][j];
     table[r][j] = ij;
     if(j == 1) prob[r][j] = prob[rm][c]  - factlm[ij + 1];
     if(j != 1) prob[r][j] = prob[r][j-1] - factlm[ij + 1];
     if(table[r][j] != table[rm][j]) goto c450;
     }
   */

   /*       ROW R IS A REPEAT OF ROW R-1 */
   /*
     reps[r] = reps[rm] + ONE;
     mult[r] = mult[rm] / reps[r];
     goto c500;

   */

   /*       IF ROW R WOULD BE BIGGER THAN ROW R-1 GO BACK AND TRY
	    DECREASING ROW R-2
   */

   c440:
   i = rm;
   goto c210;

   /*      ROW R IS ALREADY LESS THEN ROW R-1 SO NO MORE CHECKS ARE NEEDED */

   c450:
   jp = j + 1;

   for(jj = jp; jj <= c; jj++)
     {
       ij = bound[r][jj];
       table[r][jj] = ij;
       prob[r][jj] = prob[r][jj - 1] - factlm[ij + 1];
     }

   mult[r] = mult[rm];

   /*       THE TABLE IS COMPLETE - CALL SUBROUTINE EVAL */


   c500:     

    indexTable ++;

    if(indexTable > lengthTableVec){
      errmsg("Error: index for tableVec > length tableVec\n");
      //      exit(1);
    }

    probTable = exp(prob[r][c]);
    tableVec[indexTable].index = indexTable + 1;
    tableVec[indexTable].prob = probTable;

    if(nrow == 2)
      {
	for(j=1; j<=2; j++)
	  {
	    tableVec[indexTable].t32[0][j-1] = table[1][j];
	    tableVec[indexTable].t32[1][j-1] = 0 ;
	    tableVec[indexTable].t32[2][j-1] = table[2][j];
	  }
      } 
    else
      {
	for(i=1; i<=3; i++){
	  for(j=1; j<=2; j++){
	    tableVec[indexTable].t32[i-1][j-1] = table[i][j];
	  }
	}
      }


   /* CALL EVAL(TABLE,R,C,N,M,PROB(R,C),MULT(R))
       END OF MAIN LOOP
   */
   
   goto c200;


   // Every exit point goes to this spot: c600 and returns the value of 
   // indexTable at that point

   c600:
   // Free allocated containers, 
   // which may not get done with all the goto-s
   //REprintf("Freeing memory in enumTables...\n");
   free_ivector(z, 1, maxrc);
   free_ivector(rept, 1, maxrc);
   free_ivector(reps, 1, maxrc);
   free_ivector(mult, 1, maxrc);
   free_ivector(n, 1, maxrc);
   free_ivector(m, 1, maxrc);

   free_dvector(flogm, 1, nmax);
   free_dvector(factlm, 1, nmax);

   free_imatrix(table, 1, maxrc, 1, maxrc);
   free_imatrix(bound, 1, maxrc, 1, maxrc);

   free_dmatrix(prob, 1, maxrc, 1, maxrc);

   return indexTable;

}


/***********************************************************************************/

void print_imatrix(int ** mat, int nrow, int ncol, int transpose){

  int i,j;

  if(transpose)
    {
      for(i=1; i <= ncol; i++){
	for(j=1; j <= nrow; j++){
	  REprintf("%6d ", mat[j][i]);
	}
	REprintf("\n");
      }
    } 
  else
    {
      for(i=1; i <= nrow; i++){
	for(j=1; j <= ncol; j++){
	  REprintf("%6d ", mat[i][j]);
	}
	REprintf("\n");
      }
    }

}

/***********************************************************************************/

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	//v = (int *)Calloc((nh-nl+1+NR_END), int);
	if (!v) errmsg("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
  //Free(v);
}

/***********************************************************************************/

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	//v = (double *)Calloc((nh-nl+1+NR_END), double);
	if (!v) errmsg("allocation failure in dvector()"); 
	return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
  //Free(v);
}


/***********************************************************************************/

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	//m = (double **) Calloc((nrow+NR_END), double *);

	if (!m) errmsg("allocation failure 1 in matrix()");  
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	// m[nrl] = (double *) Calloc((nrow*ncol+NR_END), double);
	if (!m[nrl]) errmsg("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
  //Free(m[nrl]);
  //Free(m);
}

/***********************************************************************************/

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	//m=(int **) Calloc((nrow+NR_END), int*);

	if (!m)
	    errmsg("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	//m[nrl]=(int *) Calloc((nrow*ncol*NR_END), int);
	if (!m[nrl]) 
	     errmsg("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));

}


/*********************************************************************************/
  
void errmsg(char *string){

  /* Function to emulate "stop" of S+ - see page 134, S Programing, by
     Venables and Ripley */
  
  PROBLEM "%s", string RECOVER(NULL_ENTRY);
  
  // PROBLEM and RECOVER are defined as below in RS.h, 
  // define them in trexdriver.h if USING_R is not defined.
  // #define PROBLEM  {char R_problem_buf[R_PROBLEM_BUFSIZE];(sREprintf)(R_problem_buf,
  // #define RECOVER(x)  ),error(R_problem_buf);}
  // #define WARNING(x)	),warning(R_problem_buf);}

}


/***********************************************************************************/

int countTables(int t1, int *m, int nrow, int n, int approximation){
  /* 
     Count number of r x 2 contingency tables with
     fixed margins, by equation (2.1), or approximation (3.2).
     Gai l M,  Mantel, N. Counting the number of rxc contingency
     tables with fixed margins. JASA 1977 72:859-862

     t1 = marginal total for column one
     m = vector of rowmargin totals, array index beginning at 1
     nrow = nrows
     n set equal to 0 when called.

     Return = count of number of tables.

  */

  int j, jUpper;

  if(nrow < 1)
    return 0;

  if(approximation)
    {

      if(nrow == 1)
	{
	  n = (t1 > m[1]) ? 0 : 1;
	  return n;
	}
  
      double mu = 0.0;
      double v = 0.0;
      double q;
      double napprox;
      double twopi = 6.283185;

      for(j=1; j<=nrow; j++){
	mu += m[j];
	v  += m[j]*(m[j] + 2.0);
      }
      mu /= 2.0;
      v /= 12.0;

      napprox = 1.0;
      for(j=1; j<=nrow; j++){
	napprox = napprox*(m[j] + 1);
      }
      q = (t1 - mu)*(t1 - mu)/v;

      napprox = napprox * exp(-q / 2.0)/ sqrt(twopi * v);
      n = (int) (napprox + 1);
      return(n);
    }
  else
    {
  
      if(nrow == 1)
	{
	  n = (t1 > m[1]) ? 0 : 1;
	  return n;
	}
      
      jUpper = (t1 < m[nrow]) ? t1 : m[nrow];
      int tmpn = 0;
      for(j=0; j <= jUpper; j++){
	tmpn += countTables(t1-j, m, nrow-1, n, approximation);
      }
      return tmpn;
    }
}


/***********************************************************************************/

void printTableVec(TABLE * tableVec, int ntables){

  /* note 0-offset for tableVec and each t32 table in tableVec */

  int i, j, k;
  for(k=0; k < ntables; k++){

    REprintf("table[%d]: group = %d, last = %d, exlude = %d\n", tableVec[k].index, tableVec[k].group, tableVec[k].last, tableVec[k].exclude);
    REprintf("           prob = %f, groupCumProb = %f\n", tableVec[k].prob, tableVec[k].groupCumProb);

    for(i=0; i<3; i++){
      for(j=0; j<2; j++){
	REprintf("%6d ",tableVec[k].t32[i][j]);
      }
      REprintf("\n");
    }
  }

  return;
}

/***********************************************************************************/

void print22Tables(TABLE * tableVec, int ntables){

  /* note 0-offset for tableVec and each t32 table in tableVec */

  int i, j, k;
  for(k=0; k < ntables; k++){

    REprintf("table[%3d]: chistat = %f, sign = %d,  prob = %f\n", k+1, tableVec[k].chistat, tableVec[k].chistatSign,tableVec[k].groupCumProb);

    for(i=0; i<2; i++){
      for(j=0; j<2; j++){
	REprintf("%6d ",tableVec[k].t32[i][j]);
      }
      REprintf("\n");
    }
  }

  return;
}


/***********************************************************************************/

int  compareTable22(const void *to_one, const void *to_two){

  /* note 0-offset for each t32 table in tableVec, but using only 2x2 portion
   (i.e., skip 3rd row)
   This function is used to sort 2x2 tables, to subsequently group 2x2 tables with
   equivalent counts */

  int i,j;
  TABLE *one, *two;
  one = (TABLE *) to_one;
  two = (TABLE *) to_two;
 
  for (i=0; i<2; i++) { /* using rows 1 and 2, ignoring row 3 of 3x2 table */
    for(j=0; j<2; j++);
    if (one->t32[i][j] < two->t32[i][j])  return -1;
    if (one->t32[i][j] > two->t32[i][j])  return +1;
  }
  return 0;

}

/***********************************************************************************/

void groupCode(TABLE * tableVec, int ntables) {

  /* assume tableVec (vector of table structs)  is sorted */

  int k;
  int code;
 
  /* group code for 1st table */
  code = 1;
  tableVec[0].group = code;
 
  for(k=1; k<ntables; k++)
    {

      if(compareTable22(&tableVec[k-1], &tableVec[k]) != 0)
	{
	  /* if transition between tables, then last prior table
	     (last in its group) and increment group code */
	  code ++;
	  tableVec[k-1].last = 1;
	}
      else
	{
	  tableVec[k-1].last = 0;
	}

      /* group code for table */
      tableVec[k].group = code; 

    }
  
  /* always use last table */
  tableVec[ntables - 1].last = 1;
  
  return;
}

/***********************************************************************************/

void groupProb(TABLE * tableVec, int ntables) {

  /* assume tableVec is sorted according to group */

  int k;
  double probsum;

  if(tableVec[0].exclude == 1)
    {
      probsum =  0.0;
    } else
    {
      probsum = tableVec[0].prob;
    }

  tableVec[0].groupCumProb = probsum;

  for(k=1; k<ntables; k++)
    {

      if(tableVec[k-1].group != tableVec[k].group)
	{
	  probsum = (tableVec[k].exclude == 1) ? 0.0 : tableVec[k].prob;
	}
      else
	{
	  probsum += (tableVec[k].exclude == 1) ? 0.0 : tableVec[k].prob;
	}

      tableVec[k].groupCumProb = probsum;
    }

  return;
}

/***********************************************************************************/

int excludeTables(TABLE * tableVec, int ntables, int threshold_cell_11){

  /* exclude 3x2 tables that have cell count table[0][0] < thresheshold,
     so keep tables with cell count >= threshold */

  int k;
  int countExclude = 0;

  for(k=0; k < ntables; k++){
    tableVec[k].exclude = (tableVec[k].t32[0][0] < threshold_cell_11) ? 1 : 0;
    countExclude += tableVec[k].exclude;
  }

  return countExclude;
}

/***********************************************************************************/

int compareGroupLast(const void *to_one, const void *to_two){
  int i,j;
  TABLE *one, *two;
  one = (TABLE *) to_one;
  two = (TABLE *) to_two;
 
  /* sort such that last tables of each group are at beginning of vector,
     and sorted next according to group code */

  /* reverse sort on last */ 
  if(one->last > two->last) return -1;
  if(one->last < two->last) return +1;
  
  /* for ties on last, sort on group */
  if(one->group < two->group) return -1;
  if(one->group > two->group) return 1;
  return 0;
}

/***********************************************************************************/

void chi22Table(int table32[3][2], double *chistat, int *sign){
  double a,b,c,d;
  int i,j;

  /* compute chisquare stat for 2x2 table, first 2 rows and cols of 3x3 table 

     assumes 0-offset table32

     return -1.0 if any cells < 0 */

  for(i=0; i<2; i++){
    for(j=0; j<2; j++){
      if(table32[i][j] < 0)
	{
	  *chistat = -1.0;
	  *sign = 0;
	  return;
	}
    }
  }

  a = (double) table32[0][0];
  b = (double) table32[0][1];
  c = (double) table32[1][0];
  d = (double) table32[1][1];

  double row1, row2, col1, col2, tot;
  row1 = a + b;
  row2 = c + d;
  col1 = a + c;
  col2 = b + d;
  tot = row1 + row2;
  
  double temp  = (a*d - b*c);
  *chistat = temp*temp*tot / (row1 * row2 * col1 * col2);

  int ad = table32[0][0] * table32[1][1];
  int bc = table32[0][1] * table32[1][0];

  if(ad > bc)
    {
      *sign = 1;
    } else if(ad < bc)
    {
      *sign = -1;
    } else
    { 
      *sign = 0;
    }


  return;
}

/***********************************************************************************/

int findObs22Table(int **obsTable, TABLE * tableVec, int ntables){

  /* assume 1-offset for obsTable, but 0-offset for t32 in tableVec,
     and use only rows 1 & 2 from t32 (row indices 0 & 1) */
 
  int k;
  
  for(k=0; k < ntables; k++){
 
    if( (obsTable[1][1] == tableVec[k].t32[0][0]) &&
	(obsTable[1][2] == tableVec[k].t32[0][1]) &&
	(obsTable[2][1] == tableVec[k].t32[1][0]) &&
	(obsTable[2][2] == tableVec[k].t32[1][1]) ) 
      {
	return k;
      } 

  }

  /* return -1 if obsTable not found in tableVec */

  return -1;
}


/***********************************************************************************/  

void fillTable(int *rtablevec, int **table, int nrow, int ncol){

  int r, c;
  int k = 0;
  for(r=1; r<=nrow; r++) {
    for(c=1; c<=ncol; c++) {
      table[r][c] = rtablevec[k];
      k++;
    }
  }

  return;

}



/***********************************************************************************/

void computeChiSqPvals(TABLE *tableVec, int ngroups, int indexObs, 
    double *chistatObs, double *pvalTwoSided, double *pvalOneSided, int *signObs) {
  
  int k=0;
  
  //  double pvalChi = 0.0;
  
  *chistatObs = tableVec[indexObs].chistat;
  *signObs = tableVec[indexObs].chistatSign;
  
  //  double pvalOneSided = 0.0;
  double tol = 1e-10;
  double chiSigned, obsSigned;
  obsSigned = (*signObs) * (*chistatObs);
  *pvalTwoSided=0.0;
  *pvalOneSided=0.0;  
  for(k=0; k<ngroups; k++){

    if ( (tableVec[k].chistat > *chistatObs)  || (fabs(tableVec[k].chistat - *chistatObs) < tol ) )
      {
	*pvalTwoSided += tableVec[k].groupCumProb;
      }
    
    chiSigned = tableVec[k].chistatSign * tableVec[k].chistat;
    
    if ( (chiSigned > obsSigned) || (fabs(chiSigned - obsSigned) < tol) )
      {
	*pvalOneSided += tableVec[k].groupCumProb;
      }
  }
  return;  
}

/***********************************************************************************/

void computeFisherPvals(TABLE *tableVec, int ngroups, int indexObs, 
          double *pvalTwoSided, double *pvalOneSided, int *signObs){

  int k=0;
  //double pvalTwoSided = 0.0;

  double probObs = tableVec[indexObs].groupCumProb;
  *signObs = tableVec[indexObs].chistatSign;

  //double pvalOneSided = 0.0;
  int obsCell = tableVec[indexObs].t32[0][0];
  double tol = 1e-10;
  *pvalTwoSided = 0.0;
  *pvalOneSided = 0.0;
  for(k=0; k<ngroups; k++)
    {      
      /* two-sided p-val if prob of table same or smaller than obs table prob  */

      if ( (tableVec[k].groupCumProb < probObs)  || (fabs(tableVec[k].groupCumProb - probObs) < tol ) )
	{
	  *pvalTwoSided += tableVec[k].groupCumProb;
	}

      /* one-sided p-val in direction of cases having same or more variants than controls */

      if(tableVec[k].t32[0][0] >= obsCell)
	{
	  *pvalOneSided += tableVec[k].groupCumProb;
	}
    }

  return;
}

