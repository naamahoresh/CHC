/*************************************************************/
/*                                                           */
/*  Copyright (c) 1990                                       */
/*  Darrell L. Whitley                                       */
/*  Computer Science Department                              */
/*  Colorado State University                                */
/*                                                           */
/*  Permission is hereby granted to copy all or any part of  */
/*  this program for free distribution.   The author's name  */
/*  and this copyright notice must be included in any copy.  */
/*                                                           */
/*************************************************************/

#ifndef _GA_POOL_H_
#define _GA_POOL_H_

#ifdef ANSI

void    
free_pool (POOLPTR pool,
	   int start_pt);
			   
GENEPTR 
get_gene (int string_length);

GENEPTR
get_pool_data (int pool_size,
	       int string_length);
POOLPTR 
get_pool (int  pool_size,
	  int  string_length);
int
compare_genes (GENEPTR gene1,
	       GENEPTR gene2,
	       int     length,
	       int     sequence_flag);

int
init_pool (char   *seed_file,
	   POOLPTR pool,
           int     strt, 
	   int     stop,
	   float   (*eval_fun)(GENE_DATA *tour, int length) );

int
uniq_edge_init_pool (char   *seed_file,
	   	     POOLPTR pool,
           	     int     strt, 
	   	     int     stop,
	   	     float   (*eval_fun)(GENE_DATA *tour, int length) );

int 
isprime(int n);

int nextprime(int n);

void 
insert_gene (GENEPTR newgene, 
	     POOLPTR pool);

void
insert_unique_gene (GENEPTR newgene, 
		    POOLPTR pool, 
		    int     sequence_flag);

int
random_init_pool (POOLPTR pool, 
		  int     strt, 
		  int     stop,
		  float   (*eval_fun)(GENE_DATA *tour, int length) );

int
bias_init_pool (  float   bias,
		  POOLPTR pool, 
		  int     strt, 
		  int     stop,
		  float   (*eval_fun)(GENE_DATA *tour, int length) );

int
seed_pool (FILE    *fp,
	   POOLPTR pool, 
	   int     istart, 
	   int     istop, 
	   float   (*eval_fun)(GENE_DATA *tour, int length) );

void 
sort_pool (POOLPTR pool);


#else

int     compare_genes ();
void    free_pool ();
GENEPTR get_gene ();
GENEPTR get_pool_data ();
POOLPTR get_pool (); 
int     init_pool ();
void    insert_gene (); 
void    insert_unique_gene (); 
int     random_init_pool (); 
int     seed_pool (); 
void    sort_pool ();
int     uniq_edge_init_pool ();
int     isprime();
int     nextprime();

#endif


#endif
