  /* **************************************************************** */
  /*          Parameters needed by Statistical Info routines          */
  /* **************************************************************** */

#ifndef _STATS_PARAMS_H_
#define _STATS_PARAMS_H_

  /*  The following are used to keep track of how things do in the long
      run such as how many successful solutions are found etc.  */

int    Bin_Tries[1000];
float  Bin_Solns[1000];
float  Extra_Solns[1000];
int    cur_trials     = 0;
long   running_trials = 0;
int    successes      = 0;
float  solution       = 0.0;
long   lowest_trials  = 10000000;
long   highest_trials = 0;
float  low_soln       = 10000000.0;


#endif
