#define MAXLINE 10000
#define MAXFILENAME 25
#define MAXLABELSIZE 100
#define SEPCHAR " ,\t\n"

#define NR_END 1
#define FREE_ARG char*

#ifdef USING_R
#include <R.h>
#else
#define RECOVER(x)   ); exit(1);}
#define PROBLEM     {printf(
#define NULL_ENTRY  /**/
#endif

#ifndef TREX_DRIVER_H
#define TREX_DRIVER_H

typedef struct TABLE_T {
  int index;
  int group;
  int exclude;
  int last;
  double groupCumProb;
  double prob;
  double chistat;
  int chistatSign;
  int t32[3][2];
  int valid;
} TABLE;


/* main computational functions for enumerating tables */

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
		int verbose);

int countTables(int t1, int * m, int nrow, int n, int approximation);

int enumTableAll(int nrow,int ncol,int *rowmarg, int *colmarg, int *ifault,
		 TABLE * tableVec, int ntables);

/* functions to populate 3x2 table */
void fillTable(int *rtablevec, int **table, int nrow, int ncol);


/* memory allocation functions */
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int *ivector(long nl, long nh);
double *dvector(long nl, long nh);


/*  free() methods for use in R */
void free_dvector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);

void print_imatrix(int ** mat, int nrow, int ncol, int transpose);


/* function used to sort 2x2 tables in tableVec array */

int compareTable22(const void *to_one, const void *to_two);
void groupCode(TABLE * tableVec, int ntables);
void groupProb(TABLE * tableVec, int ntables);
int excludeTables(TABLE * tableVec, int ntables, int threshold_cell_11);
int compareGroupLast(const void *to_one, const void *to_two);

void make22(TABLE * tableVec, int ntable);
void print22Tables(TABLE * tableVec, int ntables);
void chi22Table(int table32[3][2], double *chistat, int *sign);
int findObs22Table(int **obsTable, TABLE * tableVec, int ntables);


/*  error handler, like R's stop() function in R
    Set macros above for non-R to exit program with message */
void errmsg(char *string);

/* a few print functions */
void printTableVec(TABLE * tableVec, int ntables);


/* add pvalue args so they can be returned to R instead of printed out */
void computeChiSqPvals(TABLE *tableVec, int ngroups, int indexObs, 
     double *chistatObs, double *pvalTwoSided, double *pvalOneSided, int *signObs);

void computeFisherPvals(TABLE *tableVec, int ngroups, int indexObs, 
     double *pvalTwoSided, double *pvalOneSided, int *signObs);

#endif
