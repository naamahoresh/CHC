
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
#include "ga_random.h"

/*************************************************************************/
/*  The following are static variables and functions needed for the      */
/*  random number generator                                              */
/*************************************************************************/

#define M 714025
#define IA 1366
#define IC 150889

static long Seed = 0;

#ifdef ANSI
static float random0_1 (long *idum);
#else
static float random0_1 ();
#endif



/***************************************************************************
 * FUNCTION: rand_sequence
 *
 * DESCRIBE: Randomly generates a random sequence of numbers from 
 *           0 to num_points-1
 *           (i.e. where each point is enumerated only once.)
 *
 *           Essentially, this routine fills an array with all possible
 *           numbers in the set and randomly chooses the 'next' number from
 *           this array.  When a number is chosen, the array is shortened
 *           and the procedure repeated.
 *
 * INPUT PARAMETERS: an empty tour, max number
 *
 * RETURN VALUE: none
 *
 * RESULT: fills input array
 ****************************************************************************/
#ifdef ANSI
void
rand_sequence( int      seq[],
               int      num_points)
#else
void
rand_sequence( seq, num_points)
int      seq[];
int      num_points;
#endif
{
static int          first_time = 1;
static int         *temp;
int                 remainder;
int                 next, i;

if (first_time)
   {
    if (!(temp = (int *) malloc (num_points*sizeof(int))))
 	   fatal_error ("rand_sequence: malloc failure\n");
    first_time = 0;
   }
   
for(i = 0; i < num_points; i++)
   temp[i] = (int) i;

remainder = num_points - 1;

for(i = 0; i < num_points; i++)
   {
    next = randomain(remainder, 0);
    seq[i] = temp[next];
    temp[next] = temp[remainder];
    remainder--;
   }
} 



/***************************************************************************
 * FUNCTION: random0_1
 *
 * DESCRIBE: Randomly generates numbers between 0.0 and 1.0.  This function
 *           is local to this file and it maintains the random number seed.
 *
 * INPUT PARAMETERS: A seed value that gets updated during each call.
 *
 * RETURN VALUE: float
 *
 * RESULT: random number between 0.0 and 1.0   
 *         an updated seed
 ****************************************************************************/
#ifdef ANSI
static float random0_1 (long *idum)
#else
static float random0_1 (idum)
long *idum;
#endif
{
      static long iy,ir[98];
      static int iff=0;
      int j;
      void nrerror();

      if (*idum < 0 || iff == 0) {
              iff=1;
              if ((*idum=(IC-(*idum)) % M) < 0) *idum = -(*idum);
              for (j=1;j<=97;j++) {
                      *idum=(IA*(*idum)+IC) % M;
                      ir[j]=(*idum);
              }
              *idum=(IA*(*idum)+IC) % M;
              iy=(*idum);
    }
      j=1 + 97.0*iy/M;
      if (j > 97 || j < 1) exit(0);
      iy=ir[j];
      *idum=(IA*(*idum)+IC) % M;
      ir[j]=(*idum);
      return (float) iy/M;
}


/***************************************************************************
 * FUNCTION: fracrand
 *
 * DESCRIBE: Randomly generates numbers between 0.0 and 1.0.  This is the
 *           function exported from this file.  It is does not have to
 *           worry about maintaining the random number seed.
 *
 * INPUT PARAMETERS: none
 *
 * RETURN VALUE: float
 *
 * RESULT: random number between 0.0 and 1.0   
 ****************************************************************************/
#ifdef ANSI
float fracrand (void)
#else
float fracrand ()
#endif
{
  return random0_1(&Seed);
}



/***************************************************************************
 * FUNCTION: srandom
 *
 * DESCRIBE: Set the random number seed
 *
 * INPUT PARAMETERS: A long integer to set the random number generator.
 *
 * RETURN VALUE: none
 *
 * RESULT: sets the random number seed
 ****************************************************************************/
#ifdef ANSI
void  srandom (long random_seed)
#else
void  srandom (random_seed)
long random_seed;
#endif
{
  Seed = random_seed;
}

