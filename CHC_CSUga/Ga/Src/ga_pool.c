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
#include <malloc.h>
#include "gene.h"
#include "ga_status.h"
#include "ga_random.h"
#include "ga_pool.h"
#include "ga_copy.h"

/****************************************************************************
 * FUNCTION: free_pool
 *
 * DESCRIBE: Deallocates memory for the data structures of a genetic pool
 *           which was created by get_pool().
 *           NOTE: if start_pt != 0 (i.e. we will not be freeing the
 *                 whole pool, just reducing the pool_size) then the
 *                 pool ptr will not be freed
 *                 
 *
 * INPUT PARAMETERS: pointer to a pool of GENE structures
 *                   the location at which to begin freeing genes                   
 *
 * RETURN VALUE: none
 ****************************************************************************/
#ifdef ANSI
void
free_pool ( POOLPTR    pool,
	    int        start_pt)
#else
void
free_pool (pool, start_pt)
POOLPTR    pool;
int        start_pt;
#endif
{
 int i;  /* loop counter */

 for (i=start_pt; i<pool->size; i++)
     if (pool->data[i].string)
         free (pool->data[i].string);

 /************************************
  The pool structure should only be
  freed when all genes have been freed
  ************************************/
 if (!start_pt)
    {
    free (pool->data);
    free (pool);
    }

 /******************
  Update pool's size 
  ******************/
 else
    {
    pool->size = start_pt;
    }
}

/***************************************************************************
 * FUNCTION: get_gene
 *
 * DESCRIBE: allocates a gene and string space
 *
 * INPUT PARAMETERS: length of string
 *
 * RETURN VALUE: geneptr
 ****************************************************************************/
#ifdef ANSI
GENEPTR
get_gene ( int   string_length)
#else
GENEPTR
get_gene (string_length)
int       string_length;
#endif
{
 GENEPTR gene;

 if (!(gene = (GENEPTR) malloc (sizeof(GENE)))
     ||
     !(gene->string = 
      (GENE_DATAPTR) malloc ((string_length+1)*sizeof(GENE_DATA)))
    )
    fatal_error ("get_gene(): malloc failure\n");
    
 return (gene);
}


/***************************************************************************
 * FUNCTION: get_pool_data
 *
 * DESCRIBE: Allocates memory for the gene data structures of a genetic pool
 *
 * INPUT PARAMETERS: pool_size, length of string
 *
 * RETURN VALUE: GENEPTR
 ****************************************************************************/
#ifdef ANSI
GENEPTR
get_pool_data ( int    pool_size,
		int    string_length)
#else
GENEPTR
get_pool_data (pool_size, string_length)
int            pool_size;
int            string_length;
#endif
{
 int i;
 GENEPTR pool_data;

 
 if ( !( pool_data = (GENEPTR) malloc (pool_size * sizeof(GENE))) )
    fatal_error("get_pool_data(): Malloc failure.\n");

 for (i=0; i<pool_size; i++)
 {
     if ( !(pool_data[i].string = 
            (GENE_DATAPTR) malloc ((string_length+1)*sizeof(GENE_DATA))) )
        fatal_error("get_pool_data(): Malloc failure.\n");
 }

 return (pool_data);
}


/****************************************************************************
 * FUNCTION: get_pool
 *
 * DESCRIBE: Allocates memory for the pool data structure and 
 *           the data structures of its genetic pool
 *
 * INPUT PARAMETERS: number of genes in pool;
 *                   length of a gene;
 *
 * RETURN VALUE: NULL if sufficient memory was not available,
 *               a pointer to a new pool if was available.
 ****************************************************************************/
#ifdef ANSI
POOLPTR
get_pool ( int     pool_size,
	   int     string_length)
#else
POOLPTR
get_pool (pool_size, string_length)
int        pool_size;
int        string_length;
#endif
{
 POOLPTR new_pool;

 if ( !(new_pool = (POOLPTR) malloc (sizeof(POOL))) )
    fatal_error ("get_pool(): Malloc failure.\n");

 new_pool->size = pool_size;
 new_pool->string_length = string_length;


 if ( !(new_pool->data = get_pool_data (pool_size, string_length)) )
    fatal_error ("get_pool(): Malloc failure.\n");

 return (new_pool);
}


/***************************************************************************
 * FUNCTION: init_pool
 *
 * DESCRIBE: gives initial values to genetic pool
 *
 * INPUT PARAMETERS: seed file name; 
 *                   ptr to pool of genes;
 *                   start and stop position in gene pool;
 *                   evaluation function
 *
 * RETURN VALUE: number of genes initialized
 *
 * RESULT: changes values of genes strings & worths in input pool
 *
 * CALLS: random_init_pool
 *        seed_pool
 ****************************************************************************/
#ifdef ANSI
int
init_pool ( char      *seed_file,
	    POOLPTR    pool,
	    int        strt, 
	    int        stp,
	    float      (*eval)(GENE_DATA *tour, int length) )
#else
int
init_pool (seed_file, pool, strt, stp, eval)
char      *seed_file;
POOLPTR    pool;
int        strt, stp;
float      (*eval)();
#endif
{
 int num_init;

 if (strlen (seed_file) == 0)
    num_init =  random_init_pool (pool, strt, stp, eval);

 else
    {
     FILE *fp;
     char mssg[80];

     if (!(fp = fopen (seed_file, "r")))
        { sprintf (mssg, 
                   "Cannot open pool initialization file '%s'\n", 
                   seed_file);
          fatal_error (mssg);
        }
        

     num_init = seed_pool (fp, pool, strt, stp, eval);
   
     if (num_init < 0)
        fatal_error ("init_pool: bad initialization indices");

     /*********************************************
      If file did not contain enough initialization
      values, finish off with random initialization
      *********************************************/
     if (num_init != stp-strt)
        {
        char mssg[80];
        int  rand_num_init;

        sprintf (mssg, 
    "\n%d genes read from file '%s';\n%d genes will be randomly generated.\n",
                num_init, seed_file, stp-num_init);
        warning(mssg);

        rand_num_init = random_init_pool (pool, num_init-1, stp, eval);
        num_init = num_init + rand_num_init;
        }

     fclose (fp);
    }
 return (num_init);
}


/****************************************************************************
 * FUNCTION: insert_gene
 *
 * DESCRIBE: inserts a new gene into pool, displacing worst gene in pool
 *
 *           NOTE: assumes best->worst = smallest->largest
 *                 (see sort_pool() for further details)
 *
 * INPUT PARAMETERS: gene structure;
 *                   pool of genes;
 *
 * RETURN VALUE: none
 ****************************************************************************/
#ifdef ANSI
void
insert_gene (GENEPTR   newgene,
	     POOLPTR   pool)
#else
void
insert_gene (newgene, pool)
GENEPTR      newgene;
POOLPTR      pool;
#endif
{
 int top, mid, bot;
 int rank;

 /* new gene is so bad we can't use it */
 if (newgene->worth > pool->data[pool->size-1].worth)
    return;


 /* do a binary search to find the rank of the new gene */

 top = 0;
 mid = pool->size/2;
 bot = pool->size-1;
 rank = -1;      

 while (rank == -1)
    {
     /* these 4 cases find a new location */

     if (newgene->worth <= pool->data[top].worth)
        rank = top;
     else
     if (newgene->worth == pool->data[mid].worth)
        rank = mid;
     else
     if (newgene->worth == pool->data[bot].worth)
        rank = bot;
     else
     if (bot-top <=1)
        rank = bot;


     /* These 2 cases move the search indices since
        a new location has not yet been found. */

     else if (newgene->worth < pool->data[mid].worth)
        {
         bot = mid;
         mid = top+((bot-top)/2);
        }
     else /* (newgene->worth > pool->data[mid].worth) */
        {
         top = mid;
         mid = top+((bot-top)/2);
        }
    }

 /*
  * Move every gene from rank on down
  * one position to make room for newgene.
  */

 {
  GENE travel_down, temp;
  int  i;

  /* copy new gene into pool storage;
   * always replace worst gene in pool 
   */
  gene_copy (&pool->data[pool->size-1], newgene, pool->string_length);

  travel_down.string = pool->data[pool->size-1].string;
  travel_down.worth  = pool->data[pool->size-1].worth;


  for (i=rank; i<pool->size; i++)
     {
      temp.string = pool->data[i].string;
      temp.worth  = pool->data[i].worth;

      pool->data[i].string = travel_down.string;
      pool->data[i].worth  = travel_down.worth;

      travel_down.string = temp.string;
      travel_down.worth  = temp.worth;
     }
 }
}

/****************************************************************************
 * FUNCTION: compare_genes
 *
 * RETURN VALUE: 1 if strings are different, 0 if same
 ****************************************************************************/
#ifdef ANSI
int 
compare_genes( GENEPTR  gene1,
	       GENEPTR  gene2,
	       int      length,
	       int      sequence_flag)    
#else
int compare_genes(gene1, gene2, length, sequence_flag) 
GENEPTR      gene1;
GENEPTR      gene2;
int          length;
int sequence_flag;    /* must be non-zero for a sequencing problem */
#endif
{
   int i;
   int return_value = 0;

   if (!sequence_flag) {
      for (i=0; i<length; i++) {
         if (gene1->string[i] != gene2->string[i]) {
	    return_value = 1;
	    break;
	 }
      }
   }
   else {          /* SequenceFlag is set-- its a sequencing problem */
      if (gene1->worth != gene2->worth)
	 return_value = 1;
   }
   return(return_value);
}

/****************************************************************************
 * FUNCTION: insert_unique_gene
 *
 * DESCRIBE: inserts a new gene into pool, displacing worst gene in pool,
 *           as long as it is unique.
 *
 *           NOTE: assumes best->worst = smallest->largest
 *                 (see sort_pool() for further details)
 *
 * INPUT PARAMETERS: gene structure;
 *                   pool of genes;
 *
 * RETURN VALUE: none
 ****************************************************************************/
#ifdef ANSI
void
insert_unique_gene ( GENEPTR      newgene,
		     POOLPTR      pool,
		     int          sequence_flag)    
#else
void
insert_unique_gene (newgene, pool, sequence_flag)
GENEPTR      newgene;
POOLPTR      pool;
int sequence_flag;    /* must be non-zero for a sequencing problem */
#endif
{
 int top, mid, bot;
 int rank;
 int temp_rank; 
 int different;

 /* new gene is so bad we can't use it */
 if (newgene->worth > pool->data[pool->size-1].worth)
    return;


 /* do a binary search to find the rank of the new gene */

 top = 0;
 mid = pool->size/2;
 bot = pool->size-1;
 rank = -1;      

 while (rank == -1)
    {
     /* these 4 cases find a new location */

     if (newgene->worth <= pool->data[top].worth)
        rank = top;
     else
     if (newgene->worth == pool->data[mid].worth)
        rank = mid;
     else
     if (newgene->worth == pool->data[bot].worth)
        rank = bot;
     else
     if (bot-top <=1)
        rank = bot;


     /* These 2 cases move the search indices since
        a new location has not yet been found. */

     else if (newgene->worth < pool->data[mid].worth)
        {
         bot = mid;
         mid = top+((bot-top)/2);
        }
     else /* (newgene->worth > pool->data[mid].worth) */
        {
         top = mid;
         mid = top+((bot-top)/2);
        }
    }

 /*
  * Move every gene from rank on down
  * one position to make room for newgene,
  * unless not a unique string.
  */

 different = 1;
 temp_rank = rank;
 while ( different && (newgene->worth == pool->data[temp_rank].worth) ) {
    different = compare_genes(newgene, &pool->data[temp_rank], 
			      pool->string_length, sequence_flag);
    temp_rank++;
    if (temp_rank == pool->size) 
       break;
 }

 if (different) {

 {
  GENE travel_down, temp;
  int  i;

  /* copy new gene into pool storage;
   * always replace worst gene in pool 
   */
  gene_copy (&pool->data[pool->size-1], newgene, pool->string_length);

  travel_down.string = pool->data[pool->size-1].string;
  travel_down.worth  = pool->data[pool->size-1].worth;


  for (i=rank; i<pool->size; i++)
     {
      temp.string = pool->data[i].string;
      temp.worth  = pool->data[i].worth;

      pool->data[i].string = travel_down.string;
      pool->data[i].worth  = travel_down.worth;

      travel_down.string = temp.string;
      travel_down.worth  = temp.worth;
     }
 }
 }
}




/****************************************************************************
 * FUNCTION: random_init_pool
 *
 * DESCRIBE: Randomly initializes the string data structures of a genetic pool;
 *           registers the worth of each randomly created string.
 *
 * INPUT PARAMETERS: pointer to genetic pool;
 *                   start & stop position in pool;
 *                   pointer to evaluation function;  
 *
 * RETURN VALUE: number of genes initialized
 ****************************************************************************/
#ifdef ANSI
int
random_init_pool ( POOLPTR   pool,
		   int       strt, 
		   int       stp,
		   float     (*eval_fun)(GENE_DATA *tour, int length))
#else
int
random_init_pool (pool, strt, stp, eval_fun)
POOLPTR    pool;
int        strt, stp;
float      (*eval_fun)();
#endif
{
 int i; /* loop counter */

 for (i=strt; i<stp; i++)
     {
      INIT(pool->data[i].string, pool->string_length);
      pool->data[i].worth = 
       (*eval_fun)(pool->data[i].string, pool->string_length);
     }
 return (stp-strt);
}




/****************************************************************************
 * FUNCTION: bias_init_pool
 *
 * DESCRIBE: Initializes the string data structures of a binary genetic pool
 *           according to some biased percentage of 1's;
 *           registers the worth of each randomly created string.
 *
 * INPUT PARAMETERS: bias fraction;
 *                   pointer to genetic pool;
 *                   start & stop position in pool;
 *                   pointer to evaluation function;  
 *
 * RETURN VALUE: number of genes initialized
 ****************************************************************************/
#ifdef ANSI
int
bias_init_pool ( float     bias,
		 POOLPTR   pool,
		 int       strt, 
		 int       stp,
		 float     (*eval_fun)(GENE_DATA *tour, int length))
#else
int
bias_init_pool (bias, pool, strt, stp, eval_fun)
float      bias;
POOLPTR    pool;
int        strt, stp;
float      (*eval_fun)();
#endif
{
 int i, k;                                                 /* loop counter */

 for (i=strt; i<stp; i++)
     {
	for (k = 0; k < pool->string_length; k++)
	{
	   if (fracrand() < bias)
	      pool->data[i].string[k] = '1';
	   else
	      pool->data[i].string[k] = '0';
	}

        pool->data[i].worth = 
                      (*eval_fun)(pool->data[i].string, pool->string_length);
     }
 return (stp-strt);
}

/***************************************************************************
 * FUNCTION: seed_pool
 *
 * DESCRIBE: Sets ga pool structures with intial values and worths.
 *           File format as follows:  string worth
 *                                       :     :
 *
 * INPUT PARAMETERS: (FILE *) ptr to read from, 
 *                   (GENEPTR) pool to write to,
 *                   (int) place in pool to begin writing;
 *                         place in pool to quit writing
 *                   ptr to eval function (if null, use seeded worth)
 *
 * RETURN VALUE: number of genes initialized
 *
 * CALLS: copy_string
 ****************************************************************************/
#ifdef ANSI
int
seed_pool ( FILE     *fp,
	    POOLPTR   pool,
	    int       istart, 
	    int       istop,
	    float     (*eval_fun)(GENE_DATA *tour, int length))
#else
int
seed_pool (fp, pool, istart, istop, eval_fun)
FILE     *fp;
POOLPTR   pool;
int       istart, istop;
float      (*eval_fun)();
#endif
{
 int  i,j;
 float val;

 for (i=istart; i<istop; i++)
     {
      for (j=0; j<pool->string_length; j++)
          {
          /*
          if (fscanf(fp, GENE_DATA_IN_FORMAT, &temp) != 1)
             return (i);
          pool->data[i].string[j] = GENE_DATA_IN_TRANS(temp);
          */
          if (fscanf(fp, GENE_DATA_IN_FORMAT, &pool->data[i].string[j]) != 1)
             return (i);
          }

      if (fscanf(fp, "%f", &val) != 1)
         return(i);

      if (eval_fun)
          pool->data[i].worth = 
           (*eval_fun)(pool->data[i].string, pool->string_length);
      else  
          pool->data[i].worth = val;
     }
 return (istop-istart);
}

/****************************************************************************
 * FUNCTION: sort_pool
 *
 * DESCRIBE: Sorts input pool according to worth, from smallest to largest
 *           
 *           NOTE: To affect the opposite sort (from largest to smallest)
 *                 the recommended method is to make the evaluation
 *                 function return a negative value.
 *
 *                 ( If you change sort_pool() itself (NOT recommended)
 *                   don't forget to change insert_gene() )
 *                 
 * INPUT PARAMETERS: number of genes in pool;
 *                   length of a gene;
 *
 * RETURN VALUE: none
 ****************************************************************************/
#ifdef ANSI
void
sort_pool( POOLPTR   pool)
#else
void
sort_pool(pool)
POOLPTR   pool;
#endif
{ 
 int  ndx, next, min;
 GENE temp;
        
 for (ndx=0; ndx<pool->size-1; ndx++)
    { 
     min = ndx;

     for (next=(ndx+1); next<pool->size; next++)
         if (pool->data[next].worth < pool->data[min].worth) 
             min = next;
           
     /* swap */
     if (min != ndx)
        {
	  temp.string = pool->data[ndx].string;
	  temp.worth = pool->data[ndx].worth;
	  pool->data[ndx].string = pool->data[min].string;
          pool->data[ndx].worth = pool->data[min].worth;
	  pool->data[min].string = temp.string;
          pool->data[min].worth = temp.worth;
        }
    }
}



/***************************************************************************
 * FUNCTION: nextprime
 *
 * DESCRIBE: Finds the next higher number that is a prime
 *
 * INPUT PARAMETERS: number to start search from
 *
 * RETURN VALUE: next higher integer that is a prime
 *
 * RESULT: 
 *
 * CALLS: 
 ****************************************************************************/
#ifdef ANSI
int nextprime(int n)
#else
int nextprime(n)
   int n;
#endif
{
     for (n=n+1; !isprime(n); n++) {};
     return n;
}



/***************************************************************************
 * FUNCTION: is_prime
 *
 * DESCRIBE: Determines if a number is a prime
 *
 * INPUT PARAMETERS: number to be checked
 *
 * RETURN VALUE: 1 if it is prime
 *               0 if not
 *
 * RESULT: 
 *
 * CALLS: 
 ****************************************************************************/
#ifdef ANSI
int isprime(int n)
#else
int isprime(n)
   int     n;
#endif
{
     int fact;

     if (n%2==0) return (n==2);
     if (n%3==0) return (n==3);
     if (n%5==0) return (n==5);

     fact=7;
     while (fact*fact <= n) {
          if (n%fact==0) return 0;
          fact+=4;
          if (n%fact==0) return 0;
          fact+=2;
          if (n%fact==0) return 0;
          fact+=4;
          if (n%fact==0) return 0;
          fact+=2;
          if (n%fact==0) return 0;
          fact+=4;
          if (n%fact==0) return 0;
          fact+=6;
          if (n%fact==0) return 0;
          fact+=2;
          if (n%fact==0) return 0;
          fact+=6;
     }
     return 1;
}


#ifdef INT
/***************************************************************************
 * FUNCTION: uniq_edge_init_pool
 *
 * DESCRIBE: gives initial values to genetic pool which are generated using
 *           a unique edge covering technique guaranteeing unique edges for
 *           each tour up to (Nb-1)/2 tours.
 *           
 *
 * INPUT PARAMETERS: seed file name; 
 *                   ptr to pool of genes;
 *                   start and stop position in gene pool;
 *                   evaluation function
 *
 * RETURN VALUE: number of genes initialized
 *
 * RESULT: changes values of genes strings & worths in input pool
 *
 * CALLS: seed_pool
 ****************************************************************************/
#ifdef ANSI
int
uniq_edge_init_pool ( char      *seed_file,
	    	      POOLPTR    pool,
	    	      int        strt, 
	    	      int        stp,
	    	      float      (*eval)(GENE_DATA *tour, int length) )
#else
int
uniq_edge_init_pool (seed_file, pool, strt, stp, eval)
char      *seed_file;
POOLPTR    pool;
int        strt, stp;
float      (*eval)();
#endif
{
 int num_init;
 int i, k, p, next, arc_lngth, cnt, idx, stg_idx, length;

 int *tour;

 if (strlen (seed_file) == 0)
 {
     length = pool->string_length;
     tour = (int *) malloc(length * sizeof(int));
     p = nextprime(length-1);
     cnt = (stp-strt) / ((p-1)/2);
     idx = strt;
     for (i=0; i<cnt; i++)
     {
         maketour(tour, length);
         for (arc_lngth=1; arc_lngth<=((p-1)/2); arc_lngth++)
         {
             next = 0;  stg_idx = 0;
             for (k=0; k<p; k++)
             {
                 next = (next + arc_lngth)%p;
                 if (next<length)                /* store new city "next" */
		 {
		    pool->data[idx].string[stg_idx] = tour[next];
		    stg_idx++;
		 }
             }
             idx++;
          } 
      }     

     cnt = (stp-strt) % ((p-1)/2);
     maketour(tour, length);
     for (arc_lngth=1; arc_lngth<=cnt; arc_lngth++)
       {
	 next = 0;  stg_idx = 0;
	 for (k=0; k<p; k++)
	   {
	     next = (next + arc_lngth)%p;
	     if (next<length)                 /* store new city "next" */
	       {
		 pool->data[idx].string[stg_idx] = tour[next];
		 stg_idx++;
	       }
	   }
	 idx++;
       } 

     
     free(tour); 
     for (i=strt; i<stp; i++)
         pool->data[i].worth = 
                        (*eval)(pool->data[i].string, pool->string_length);

      num_init = stp - strt; 
 
 }
 else
    {
     FILE *fp;
     char mssg[80];

     if (!(fp = fopen (seed_file, "r")))
        { sprintf (mssg, 
                   "Cannot open pool initialization file '%s'\n", 
                   seed_file);
          fatal_error (mssg);
        }
        

     num_init = seed_pool (fp, pool, strt, stp, eval);
   
     if (num_init < 0)
        fatal_error ("init_pool: bad initialization indices");

     /*********************************************
      If file did not contain enough initialization
      values, finish off with random initialization
      *********************************************/
     if (num_init != stp-strt)
        {
        char mssg[80];
        int  rand_num_init;

        sprintf (mssg, 
    "\n%d genes read from file '%s';\n%d genes will be randomly generated.\n",
                num_init, seed_file, stp-num_init);
        warning(mssg);

        rand_num_init = random_init_pool (pool, num_init-1, stp, eval);
        num_init = num_init + rand_num_init;
        }

     fclose (fp);
    }
 return (num_init);
}
#endif

