
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
/****************************************************************************
                                 ga_xdr.c

 This file contains routines which allow genetic pool data to be written to
 and read from (binary) disk files using the "External Data Representation"
 (XDR) standard.  The XDR standard is most frequently used in conjunction
 with the RPC utility, but need not be and, indeed, is useful outside the
 domain of RPC procedures.  XDR was designed to work across different
 languages, operating systems, and machine architectures.  Thus, an XDR 
 format binary file may be written on one machine, say a VAX running VMS,
 and read on another machine, say a SUN Sparc running Unix. Even when machine
 portability is not an issue, XDR is still a good choice because its protocols
 are standardized and reliable.

 The implementation found in this file was constructed using :
 --> Sun Technical Notes, Revision A, of May 9, 1988.
 ****************************************************************************/

#include <stdio.h>
#include <rpc/rpc.h>
#include "gene.h"           /* this order is not negotiable.  gene.h must*/
#include "ga_xdr.h"         /* proceed ga_xdr.h for protecting AIX builds */
			    /* and ga_xdr.h must come after rpc.h         */


unsigned int xdrStringLen;  /* use of this global should be
			       restricted to this file. */

/***************************************************************************
 * FUNCTION: xdr_GENE
 *
 * DESCRIBE: XDR provides primitive routines for basic data types:
 *           char,int,long,short,float,double,string,array,union,pointer
 *           These primitives can be used to construct XDR routines for
 *           arbitrary C data structures, as is done here for the GENE
 *           structures which make up the genetic pool.
 *
 *           XDR routines are "direction independent"; that is, the same
 *           routines are called to write (serialize) and read (deserialize)
 *           data. This represents good software engineering; it guarantees
 *           that serialized data can be deserialized.  
 *
 * INPUT PARAMETERS: pointer to XDR structure,
 *                   pointer to GENE structure,
 *                  
 *                   length of GENE string
 *
 * RETURN VALUE: 1 for success, 0 for failure
 *
 * NOTES: XDR_DATA is dependent on the actual data type of the string
 *        and thus, it is defined in gene.h
 *
 *        The global xdrStringLen is required for xdr_array.  It is not
 *        passed to this routine, however, because the input arguments
 *        conform to the XDR standard.
 ****************************************************************************/
#ifdef ANSI
int xdr_GENE( XDR     *xdrs,
	  GENEPTR  gene)
#else
int xdr_GENE(xdrs, gene)
XDR     *xdrs;
GENEPTR  gene;
#endif
{
 return(xdr_array(xdrs, &gene->string, &xdrStringLen, 
		  xdrStringLen, sizeof(GENE_DATA), XDR_DATA) &&
		  xdr_float(xdrs, &gene->worth)
	   );
}


/***************************************************************************
 * FUNCTION: process_xdr
 *
 * DESCRIBE: Traverses the pool structure, calling xdr_GENE to read or 
 *           write each GENE.
 *
 * INPUT PARAMETERS: activated xdrs (direction can be XDR_ENCODE or XDR_DECODE)
 *                   pointer to pool of GENEs
 *                   number of GENEs in pool
 *                   length of GENE strings
 *
 * RETURN VALUE: 1 for success, 0 for failure
 ****************************************************************************/
#ifdef ANSI
int
process_xdr ( XDR         *xdrs,
	      POOLPTR      pool,
	      unsigned int pool_size, 
	      unsigned int string_len)
#else
int
process_xdr (xdrs, pool, pool_size, string_len)
XDR         *xdrs;
POOLPTR      pool;
unsigned int pool_size, string_len;
#endif
{
 int i;

 xdrStringLen = string_len;

 for (i=0; i<pool_size; i++)
      if (!xdr_GENE(xdrs, (GENEPTR)(&pool[i])))
		 { fprintf (stderr, "xdr_GENE failure\n");
		   return (0);
         }	

 return (1);
}


/***************************************************************************
 * FUNCTION: print_xdr
 *
 * DESCRIBE: prints contents of ga gene pool to input file 
 *           in XDR format
 *
 * INPUT PARAMETERS: file pointer,
 *                   POOLPTR,
 *
 * RETURN VALUE: 1 for success, 0 for failure
 ****************************************************************************/
#ifdef ANSI
int
print_xdr ( FILE        *fp, 
	    POOLPTR      pool)
#else
int
print_xdr (fp, pool)
FILE        *fp; 
POOLPTR      pool;
#endif
{
 XDR  xdrs;

 xdrstdio_create (&xdrs, fp, XDR_ENCODE);
 return (process_xdr (&xdrs, pool, pool->size, pool->string_length));
}

/***************************************************************************
 * FUNCTION: read_xdr
 *
 * DESCRIBE: reads contents of ga gene pool from input file 
 *           in XDR format
 *
 * INPUT PARAMETERS: file pointer,
 *                   POOLPTR,
 *
 * RETURN VALUE: 1 for success, 0 for failure
 ****************************************************************************/
#ifdef ANSI
int
read_xdr ( FILE        *fp,
	   POOLPTR      pool)
#else
int
read_xdr (fp, pool)
FILE        *fp;
POOLPTR      pool;
#endif
{
 XDR      xdrs;

 xdrstdio_create (&xdrs, fp, XDR_DECODE);
 return (process_xdr (&xdrs, pool, pool->size, pool->string_length));
}
