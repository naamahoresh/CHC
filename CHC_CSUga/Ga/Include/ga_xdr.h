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

#ifndef _GA_XDR_H_
#define _GA_XDR_H_


#ifdef ANSI

int
xdr_GENE (XDR *xdrs,
          GENEPTR  gene);

int
process_xdr (XDR     *xdrs,
             POOLPTR  pool,
             unsigned int pool_size,
	     unsigned int string_len);

int
print_xdr (FILE      *fp,
           POOLPTR      pool);

int
read_xdr (FILE        *fp,
          POOLPTR      pool);

#else

int xdr_GENE ();
int process_xdr ();
int print_xdr ();
int read_xdr ();


#endif

#endif
