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
#include <stdio.h>
#include "gene.h"
#include "ga_random.h"
#include "ga_copy.h"
#include "ga_tournament.h"

#define TRUE  1
#define FALSE 0

/****************************************************************************
 * FUNCTION: tourny_insert_gene
 *
 * DESCRIBE: inserts a new gene into pool, displacing string in population
 *           pointed to by replace index.
 *
 * INPUT PARAMETERS: gene structure;
 *                   pool of genes;
 *                   index to gene to replace
 *
 * RETURN VALUE: none
 ****************************************************************************/
#ifdef ANSI
void
tourny_insert_gene (GENEPTR   newgene,
	            POOLPTR   pool,
	            int       replace)
#else
void
tourny_insert_gene (newgene, pool, replace)
GENEPTR      newgene;
POOLPTR      pool;
int          replace;
#endif
{
 
 /* new gene is worse so don't use it */
 if (newgene->worth > pool->data[replace].worth)
    return;

  gene_copy (&pool->data[replace], newgene, pool->string_length);

}


/***************************************************************************
 * FUNCTION: get_tourny_parents
 *
 * DESCRIBE: according to tournament size, select strings from pool
 *           and then find the best two genes in the tournament
 *
 * INPUT PARAMETERS: GENEPTRs to parents, 
 *                   POOLPTR  to pool,
 *                   float    tournament_size
 *                   int      pointer to population index for replacement
 *
 * RETURN VALUE: none
 *
 * RESULT: input mom and dad pointers are given values of two genes in pool
 *         and replace is set to the index of the gene in the population
 *         that is the worst in the tournament and should be replaced.
 *
 ****************************************************************************/
#ifdef ANSI
void
get_tourny_parents ( GENEPTR  mom, 
	             GENEPTR  dad,
	             POOLPTR  pool,
	             float    tourn_sz,
	             int      *replace)
#else
void
get_tourny_parents (mom, dad, pool, tourn_sz, replace)
GENEPTR  mom, dad;
POOLPTR  pool;
float    tourn_sz;
int      *replace;
#endif
{
  int      k, m, indx, wst_idx, Finished; 
  int      dxs[50];
  float    best, secnd, last;

  indx = (int) (fracrand() * pool->size);    /* choose the first gene */
  best = pool->data[indx].worth;
  dxs[0] = indx;

  indx = (int) (fracrand() * pool->size);    /* choose a second gene */
  while (indx == dxs[0])                     /* make sure it is not the same */
     indx = (int) (fracrand() * pool->size);
 
  if (pool->data[indx].worth < best)
  {
     secnd = best;
     best = pool->data[indx].worth;
     dxs[1] = dxs[0];
     dxs[0] = indx;
  }
  else
  {
     secnd = pool->data[indx].worth;
     dxs[1] = indx;
  }

  last = secnd;   wst_idx = dxs[1];
  for(k=2; k<tourn_sz; k++)
  {
     Finished = FALSE;
     while (!Finished)
     {
        indx = (int) (fracrand() * pool->size);   /* get next gene */
	Finished = TRUE;
	for(m=0; m<k; m++)
	{
	   if (indx == dxs[m])             /* make sure it's not a duplicate */
	   {
	      Finished = FALSE;    break;
	   }
	}
     }
     if (pool->data[indx].worth < best)
     {
	secnd = best;
	best = pool->data[indx].worth;
	dxs[k] = dxs[1];
	dxs[1] = dxs[0];
	dxs[0] = indx;
     }
     else
     {
	if (pool->data[indx].worth < secnd)
	{
	   secnd = pool->data[indx].worth;
	   dxs[k] = dxs[1];
	   dxs[1] = indx;
	}
	else
	{
	   if (pool->data[indx].worth > last)
	   {
	      last = pool->data[indx].worth;
	      wst_idx = indx;
	   }
	}
     }
  }

  gene_copy (mom, &pool->data[dxs[0]], pool->string_length);
  gene_copy (dad, &pool->data[dxs[1]], pool->string_length);
  *replace = wst_idx;
}


/***************************************************************************
 * FUNCTION: tourny_status
 *
 * DESCRIBE: Computes average value of pool worths as well as best and worst.
 *
 * INPUT PARAMETERS: (POOLPTR) pool
 *                   pointer to best value in pool variable
 *                   pointer to worst value in pool variable
 *
 * RETURN VALUE: (float) average value of pool
 *               (float) the value of the best worth in the pool
 *               (float) the value of the worst worth in the pool
 *
 * CALLS: fatal_error
 ****************************************************************************/
#ifdef ANSI
float
tourny_status (POOLPTR   pool,
	       float     *best,
	       float     *worst)
#else
float
tourny_status (pool, best, worst)
     POOLPTR   pool;
     float     *best, *worst;
#endif
{
 int     i;
 float   bst, wst;
 double  cumulative = 0.0;
 
 if (pool->size==0)
	fatal_error ("avg_pool: pool_size of zero.\n");

 wst = bst = pool->data[0].worth;
 for (i=0; i<pool->size; i++)
 {
     cumulative = cumulative + pool->data[i].worth;
     if (pool->data[i].worth < bst)
	bst = pool->data[i].worth;
     if (pool->data[i].worth > wst)
	wst = pool->data[i].worth;
 }

 *best = bst;    *worst = wst;
 return ((float) cumulative/pool->size);
}



/***************************************************************************
 * FUNCTION: print_tourny_head
 *
 * DESCRIBE: Used to print a brief progress header for the tournament report
 *
 * INPUT PARAMETERS: FILE ptr for printing out
 *
 * RETURN VALUE: none
 *
 * CALLS:
 ****************************************************************************/
#ifdef ANSI
void
print_tourny_head (FILE     *fp)
#else
void
print_tourny_head(fp)
FILE         *fp;
#endif
{
  fprintf (fp, "\nTrial   Bst           Wst           Avg\n");
}



/***************************************************************************
 * FUNCTION: tourny_progress
 *
 * DESCRIBE: Used to print following information:
 *                Best 
 *                Worst
 *                Average
 *
 * INPUT PARAMETERS: FILE ptr for printing out
 *                   genetic pool
 *                   current generation number
 *
 * RETURN VALUE: none
 *
 * CALLS:
 ****************************************************************************/
#ifdef ANSI
void
tourny_progress (FILE     *fp,
	         POOLPTR   pool,
	         int       current_generation,
		 float     *bst)
#else
void
tourny_progress(fp, pool, current_generation, bst)
FILE         *fp;
POOLPTR       pool;
int           current_generation;
float         *bst;
#endif
{
  float avg, best, worst;
 
  avg = tourny_status(pool, &best, &worst);
  fprintf (fp, "%5d  %f   %f   %f\n",
	   current_generation, best,  worst, avg);

  *bst = best; 
}

