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

#ifndef _GA_TOURNAMENT_H_
#define _GA_TOURNAMENT_H_


#ifdef ANSI

void
tourny_insert_gene (GENEPTR   newgene,
	            POOLPTR   pool,
	            int       replace);

void
get_tourny_parents ( GENEPTR  mom, 
	             GENEPTR  dad,
	             POOLPTR  pool,
	             float    tourn_sz,
	             int      *replace);

float
tourny_status (POOLPTR   pool,
	       float     *best,
	       float     *worst);

void
print_tourny_head (FILE     *fp);

void
tourny_progress (FILE     *fp,
	         POOLPTR   pool,
	         int       current_generation,
		 float     *bst);

#else

void tourny_insert_gene ( );
void get_tourny_parents ( );
float tourny_status ( );
void print_tourny_head ( );
void tourny_progress ( );

#endif

#endif
