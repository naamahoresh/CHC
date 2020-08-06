
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
#include "ga_status.h"
#include "ga_params.h"


/***************************************************************************
 * FUNCTION: avg_pool
 *
 * DESCRIBE: Computes average value of pool worths.
 *
 * INPUT PARAMETERS: (POOLPTR) pool
 *
 * RETURN VALUE: (float) average value of pool
 *
 * CALLS: fatal_error
 ****************************************************************************/
#ifdef ANSI
float
avg_pool ( POOLPTR   pool)
#else
float
avg_pool (pool)
POOLPTR   pool;
#endif
{
 int i;
 double   cumulative = 0.0;
 
 if (pool->size==0)
	fatal_error ("avg_pool: pool_size of zero.\n");

 for (i=0; i<pool->size; i++)
	 cumulative = cumulative + pool->data[i].worth;
 
 return ((float) cumulative/pool->size);
}

/***************************************************************************
 * FUNCTION: print_pool
 *
 * DESCRIBE: prints contents of ga gene pool to input file from start_pt
 *           to start_pt+count
 *
 * INPUT PARAMETERS: file pointer,
 *                   POOLPTR,
 *                   start point, count
 *
 * RETURN VALUE: none
 *
 * CALLS:
 ****************************************************************************/
#ifdef ANSI
void
print_pool ( FILE       *fp,
	     POOLPTR     pool,
	     int         start_pt, 
	     int         count)
#else
void
print_pool (fp, pool, start_pt, count)
FILE       *fp;
POOLPTR     pool;
int         start_pt, count;
#endif
{
 int i, j;

 /* be extra careful that start and count are valid inputs */

 if (start_pt < 0)
	start_pt = 0;

 if (count > pool->size)
	count = pool->size;

 if (start_pt+count > pool->size)
	{ start_pt = 0; 
	  count = pool->size;
	}

 for (i=start_pt; i<count; i++)
	 {
	  for (j=0; j<pool->string_length; j++)
          fprintf (fp, GENE_DATA_OUT_FORMAT(pool->data[i].string[j]));
      fprintf (fp, " %f\n", pool->data[i].worth);
	 }
}

/***************************************************************************
 * FUNCTION: show_progress
 *
 * DESCRIBE: Used to print following information:
 *                Best 
 *                Worst
 *                Mean 
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
show_progress( FILE      *fp,
	       POOLPTR    pool,
	       int        current_generation)
#else
void
show_progress(fp, pool, current_generation)
FILE         *fp;
POOLPTR       pool;
int           current_generation;
#endif
{
 int lowest;
 
 /* Get index to lowest ranking gene in poplulation. */
 /* Use 2nd to last since last is buffer. */ 
 lowest = pool->size > 1 ? pool->size-2 : 0;

 fprintf (fp, 
		  "%5d | Bst: %f  Wst: %f  Mean: %f  Avg: %f\n",
		   current_generation,
		   pool->data[0].worth,
		   pool->data[lowest].worth, 
		   pool->data[pool->size/2].worth,
		   avg_pool(pool));
}



/***************************************************************************
 * FUNCTION: print_prog_head
 *
 * DESCRIBE: Used to print a brief progress header
 *
 * INPUT PARAMETERS: FILE ptr for printing out
 *
 * RETURN VALUE: none
 *
 * CALLS:
 ****************************************************************************/
#ifdef ANSI
void
print_prog_head( FILE     *fp)
#else
void
print_prog_head(fp)
FILE         *fp;
#endif
{
  fprintf (fp, "\nTrial   Bst           Wst           Median        Avg\n");
}



/***************************************************************************
 * FUNCTION: show_progress_brief
 *
 * DESCRIBE: Used to print following information:
 *                Best 
 *                Worst
 *                Median 
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
show_progress_brief( FILE     *fp,
		     POOLPTR   pool,
		     int       current_generation)
#else
void
show_progress_brief(fp, pool, current_generation)
FILE         *fp;
POOLPTR       pool;
int           current_generation;
#endif
{
 int lowest;
 
 /* Get index to lowest ranking gene in poplulation. */
 /* Use 2nd to last since last is buffer. */ 
 lowest = pool->size > 1 ? pool->size-2 : 0;

 fprintf (fp, 
		  "%5d  %f   %f   %f   %f\n",
		   current_generation,
		   pool->data[0].worth,
		   pool->data[lowest].worth, 
		   pool->data[pool->size/2].worth,
		   avg_pool(pool));
}



/***************************************************************************
 * FUNCTION: show_progress_exp
 *
 * DESCRIBE: Used to print following information:
 *                Best 
 *                Worst
 *                Median 
 *                Average
 *           in scientific format
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
show_progress_exp( FILE     *fp,
		   POOLPTR   pool,
		   int       current_generation)
#else
void
show_progress_exp(fp, pool, current_generation)
FILE         *fp;
POOLPTR       pool;
int           current_generation;
#endif
{
 int lowest;
 
 /* Get index to lowest ranking gene in poplulation. */
 /* Use 2nd to last since last is buffer. */ 
 lowest = pool->size > 1 ? pool->size-2 : 0;

 fprintf (fp, 
		  "%5d  %e   %e   %e   %e\n",
		   current_generation,
		   pool->data[0].worth,
		   pool->data[lowest].worth, 
		   pool->data[pool->size/2].worth,
		   avg_pool(pool));
}


/***************************************************************************
 * FUNCTION: final_pool
 *
 * DESCRIBE: after GA main loop finishes executing, this function is called
 *           to summarize results
 *
 * INPUT PARAMETERS: (char *)  filename for saving final pool (or NULL)
 *                   (GENEPTR) to pool structure
 *                   (int)     number of genes in pool
 *                   (int)     number of positions/gene
 *
 * RETURN VALUE: none
 *
 * CALLS: print_pool
 ****************************************************************************/
#ifdef ANSI
void
final_pool ( char        filename[],
	     POOLPTR     pool,
	     int         current_generation)
#else
void
final_pool (filename, pool, current_generation)
char        filename[];
POOLPTR     pool;
int         current_generation;
#endif
{
 FILE *fp;

 show_progress (stderr, pool, current_generation);
  
 if (filename
	 &&
	 (fp = fopen (filename, "w")))
	 {
	 print_pool (fp, pool, 0, pool->size);
	 fclose (fp);
	 }
}


/***************************************************************************
 * FUNCTION: dump_status
 *
 * DESCRIBE: Used to print a snapshot of the current process, which can
 *           later be used to initialize a new process.  Snapshot exists
 *           in two files: config file (.config) for ga_params,
 *                         pool file (.pool) for population
 *
 * INPUT PARAMETERS: genetic pool
 *                   dump file basename
 *
 * RETURN VALUE: none
 *
 * CALLS: print_params
 *        print_pool
 ****************************************************************************/
#ifdef ANSI
void
dump_status( POOLPTR   pool,
	     char      dump_base[])
#else
void
dump_status(pool, dump_base)
POOLPTR     pool;
char        dump_base[];
#endif
{
 char config_file[80];
 char pool_file[80];
 FILE *fp;
 char mssg[100];
 
 sprintf (config_file, "%s.config", dump_base);
 sprintf (pool_file,   "%s.pool",   dump_base);

 
 if (!(fp = fopen (pool_file, "w")))
	{ sprintf (mssg, 
			   "dump_status: unable to open pool file '%s' for write.\n",
			   pool_file
              );
      warning (mssg);

	  /* i stands for InitializationFile ; must be set so that the config
		 file contains correct pool file name for a restart */
	  set_parameter ('i', NULL);
    }
 else
	{ print_pool (fp, pool, 0, pool->size);

	  /* i stands for InitializationFile ; must be set so that the config
		 file contains correct pool file name for a restart */
	  set_parameter ('i', pool_file);
	  fclose (fp);
    }

 if (!(fp = fopen (config_file, "w")))
	{ sprintf (mssg, 
			   "dump_status: unable to open config file '%s' for write.\n",
			   config_file
              );
      warning (mssg);
    }

 else
	{ print_params(fp);
	  fclose (fp);
    }
}


/***************
 * FATAL ERROR *
 ***************/
#ifdef ANSI
void
fatal_error( char  *mssg)
#else
void
fatal_error(mssg)
char       *mssg;
#endif
{
 if (mssg)
	fprintf (stderr, "\n%s", mssg);

 fprintf (stderr, "\nFATAL ERROR: Exiting.\n\n");
 exit (-1);
}


/***********
 * WARNING *
 ***********/
#ifdef ANSI
void
warning( char  mssg[])
#else
void
warning(mssg)
char    mssg[];
#endif
{
 if (mssg)
	fprintf (stderr, "\n%s\n", mssg);
}

/**********************
 * for DEBUG or TRACE *
 **********************/
#ifdef ANSI
void
pause_it(void)
#else
void
pause_it()
#endif
{
 char inp[50];

 fprintf (stderr, "Paused...type 'c <CR>' to continue.\n");
 fscanf (stdin, "%s", inp);
}

