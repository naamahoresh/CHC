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
#ifndef _OP_CATACLYSM_H_
#define _OP_CATACLYSM_H_


#ifdef ANSI

void 
cataclysm( POOLPTR    pool,
	   int        size,
	   int        length,
	   int        seed_idx,
	   float      pct_mut);

void 
parm_cataclysm( POOLPTR    pool,
	        int        size,
	        int        nb_parms,
                int        sparse,
                int        parm_lngth,
	        int        seed_idx,
	        float      pct_mut);
#else

extern void cataclysm( );
extern void parm_cataclysm( );

#endif
#endif
