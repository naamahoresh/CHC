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

#ifndef _GA_COPY_H_
#define _GA_COPY_H_

#ifdef ANSI

void 
gene_copy (GENEPTR gene1, 
	   GENEPTR gene2,
	   int     string_length);
#else

void
gene_copy ();

#endif


#endif
