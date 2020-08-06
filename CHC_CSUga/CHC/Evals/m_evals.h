#ifndef _M_EVALS_H_
#define _M_EVALS_H_

#include "gene.h"
#include "binarys.h"


#ifdef ANSI

float 
get_tru_fit(void);

void 
set_tru_fit(float val);

void
prblm_init_dejong(char   *NodeFile,
                  int    *Nb_params,
                  double **Params,
                  struct p_list **Work_params,
		  int     Length,
		  char   **gray_wrk);

void
raw_init_dejong(char   *NodeFile,
              int    *Nb_params,
              double **Params,
              struct p_list **Work_params,
	      int Length, FILE **parmfil,
	      int *n_degray, char *val_mode);

void
sym_trnsfrm( GENE_DATA *buff,
             double *Vals,
             int Nb_params,
             struct p_list *Work_params);

float  
mstr_eval( GENE_DATA *buff,
           int       length);

float  
shrt_eval( GENE_DATA *buff,
           int       length,
	   int       opt);

float  
raw_eval( GENE_DATA *buff,
	  double *Params,
	  int    n_degray,
	  char   val_mode ); 

float
cplx_fx( double *Params,
	 int    Nb_Params,
         float  (*eval)(double *Params, int Nb_Params),
	 char   cplxty);

float
compos_fx( double *Params,
	   int    Nb_Params,
	   int    cf1,
	   int    cf2,
	   char   cplxty,
	   char   cform);

void f_fx ( int    cf, 
	    float  (**eval)( ));

void  
mstr_rpt( GENE_DATA *buff);

void  
raw_rpt(double *params);

float 
gaussian( );

float  
f1( double *Params,
    int    Nb_params);

float  
f2( double *Params,
    int    Nb_params);

float  
f3( double *Params,
    int    Nb_params);

float  
f4( double *Params,
    int    Nb_params);

float  
f5( double *Params,
    int    Nb_params);

float  
f6( double *Params,
    int    Nb_params);

float  
f7( double *Params,
    int    Nb_params);

float  
f8( double *Params,
    int    Nb_params);

float  
f9( double *Params,
    int    Nb_params);

float  
f10( double *Params,
    int    Nb_params);

float  
f20( double *Params,
    int    Nb_params);

float  
f21( double *Params,
    int    Nb_params);

float  
f22( double *Params,
    int    Nb_params);

float
f97( GENE_DATA  *buff,
     double     *Params,
     int        Nb_params,
     int        c1,
     int        num_dgry);

float
f98( GENE_DATA  *buff,
     double     *Params,
     int        Nb_params,
     int        c1,
     int        terms);

float
f99( GENE_DATA  *buff,
     double     *Params,
     int        Nb_params,
     int        c1);

float
f100( double *Params,
    int    Nb_params); 

#else

float get_tru_fit( );
void set_tru_fit( );
void prblm_init_dejong( );
void raw_init_dejong( );
void sym_trnsfrm( );
float mstr_eval( );
float shrt_eval( );
float raw_eval( );
float cplx_fx( );
float compos_fx( );
void f_fx ( );
void mstr_rpt( );
void raw_rpt( );
float gaussian( );
float f1( );
float f2( );
float f3( );
float f4( );
float f5( );
float f6( );
float f7( );
float f8( );
float f9( );
float f10( );
float f20( );
float f21( );
float f22( );
float f97( );
float f98( );
float f99( );
float f100( );
#endif

#endif
