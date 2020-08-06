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
#ifndef _GA_SELECTION_H_
#define _GA_SELECTION_H_


#ifdef ANSI

int 
linear (int max,
	float bias);

void 
get_parents (GENEPTR mom, 
	     GENEPTR dad, 
	     POOLPTR pool,
	     float bias,
	     int (*bias_fun)(int max, float bias));
#else

int  linear ();
void get_parents ();

#endif


#endif
