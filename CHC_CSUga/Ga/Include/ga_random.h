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

#ifndef _GA_RANDOM_H_
#define _GA_RANDOM_H_


#include <math.h>        /* Needed for the floor function */


#define bitgen()         ((fracrand()>.5)? 1 : 0 )
#define bitgenc()        ((fracrand()>.5)?'1':'0')
#define randomain(hi,lo) ((int) floor(fracrand()*((hi-lo)+0.999999))+lo)


#ifdef ANSI

void 
rand_sequence (int seq[],
	       int num_points);
float fracrand (void);
void  srandom (long random_seed);

#else

void  rand_sequence ();
float fracrand ();
void  srandom ();
#endif


#endif
