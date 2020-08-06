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
#include "ga_copy.h"

/***************************************************************************
 * FUNCTION: gene_copy
 *
 * DESCRIBE: This routine copies one gene to another
 *           (genes are of the same type)
 *
 * INPUT PARAMETERS: 2 gene pointers, string_length
 *
 * RETURN VALUE: none
 ****************************************************************************/
#ifdef ANSI
void
gene_copy ( GENEPTR    gene1, 
	    GENEPTR    gene2,
            int        string_length)
#else
void
gene_copy (gene1, gene2, string_length)
GENEPTR    gene1, gene2;
int        string_length;
#endif
{
 int i;

 for (i=0; i<string_length; i++)
	 gene1->string[i] = gene2->string[i];
 gene1->worth = gene2->worth;
}
