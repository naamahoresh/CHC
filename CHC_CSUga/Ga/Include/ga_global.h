
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
/*********************************************************************
 * These REQUIRED and OPTIONAL parameters can be chosen by the user. *
 *********************************************************************/

 /* REQUIRED */
 int     PoolSize;
 int     StringLength;
 int     NumberTrials;

 char    NodeFile[80];       /* contains coordinates of tsp nodes */
 char    UserFile[80];       /* contains name of user specified file */

 /* OPTIONAL */
 long    RandomSeed;
 float   SelectionBias;
 float   MutateRate = 0.0;
 int     StatusInterval = 0; /* dump status every nth generation */
 int     DumpInterval = 0;   /* save state of every nth generation */
 char    SeedPool[80];       /* file containing init data */
 char    FinalPool[80];      /* file containing final results */
 char    DumpBase[80];       /* basename of file(s) into
                                which to dump population */

 int     NumPop;             /* Number of Subpopulations */
 int     SwapInterval;       /* Trials between Swapping of Subpopulations */
 int     SwapNumber;         /* Number of Strings Swapped between Subpops */
 int     Experiments;        /* Experiments must be used in main() */
 float   CutOff;             /* Cutoff value for a given experiment */

/***************************************************/

 int SequenceFlag = 0;       /* if set to 1 in main, insert_unique_gene() */
			     /* will use different criteria for sameness  */
 int     CurrentGeneration = 0;  
 POOLPTR Pool;
