#ifndef _BINARYS_H_
#define _BINARYS_H_

#include "gene.h"


#define SIGNED      1
#define UNSIGNED    0
#define NO_CHANGE 1.0
#define HAM_LIMIT 1
#define SAME 1
#define APART 0



struct p_list {
     int lngth;
     int sign;
     int posn;
     double scale;
     double wrap_low;
     double wrap_high; };


#ifdef ANSI

void
allocate_bin_params (double **prb_params,
                     struct p_list **Work_params,
                     int Nb_params,
		     int Length,
		     char **gray_wrk);
double
ctod ( char   *string,
       int    length,
       double step_value,
       int    sign);

void
transform( GENE_DATA         *buff,
              double         *Vals,
              int            Nb_params,
              struct p_list  *Work_params);


int
ham_dist( GENE_DATA *bf1,
	  GENE_DATA *bf2,
          int Nb_params,
	  struct p_list *Work_params);

int
hamming_diff( GENE_DATA *bf1,
	      GENE_DATA *bf2,
	      int       length);

int
hamming_distance(GENE_DATA *bf1,
                 GENE_DATA *bf2,
                 int       length);

void 
gray (char *inbuf,
      char *outbuf,
      int  Nb_params,
      struct p_list *Work_params);

void 
degray (char *inbuf,
        char *outbuf,
        int  Nb_params,
	struct p_list *Work_params);

#else

void allocate_bin_params( );
double ctod( );
void transform( );
int ham_dist( );
int hamming_diff( );
int hamming_distance( );
void gray ( );
void degray ( );

#endif

#endif
