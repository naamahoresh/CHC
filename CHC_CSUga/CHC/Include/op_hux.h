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
#ifndef _OP_HUX_H_
#define _OP_HUX_H_


#ifdef ANSI

int * 
allocate_hux_params(int stringlength);

void 
hux( char    *mom,
     char    *dad,
     int     *dif_ary,
     int     length);

#else

extern int * allocate_hux_params( );
extern void hux( );
#endif

#endif
