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

#ifndef _GA_STATUS_H_
#define _GA_STATUS_H_


#ifdef ANSI

float
avg_pool (POOLPTR pool);

void
print_pool (FILE *fp, 
	    POOLPTR pool,
	    int start_pt,
	    int count);
void
show_progress (FILE *fp, 
	       POOLPTR pool,
	       int current_generation);

void
print_prog_head (FILE *fp); 

void
show_progress_brief (FILE * fp, 
		     POOLPTR pool,
		     int current_generation);

void
show_progress_exp (FILE * fp, 
		   POOLPTR pool,
		   int current_generation);

void
final_pool (char filename[],
	    POOLPTR pool,
	    int current_generation);

void
dump_status (POOLPTR pool,
             char dump_base[]);
void
fatal_error (char *mssg);

void
warning (char mssg[]);

void
pause_it (void);

#else

float avg_pool ();
void  print_pool ();
void  show_progress ();
void  print_prog_head (); 
void  show_progress_brief ();
void  show_progress_exp ();
void  final_pool ();
void  dump_status ();
void  fatal_error ();
void  warning ();
void  pause_it ();

#endif



#endif
