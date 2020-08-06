#ifndef _STATS_H_
#define _STATS_H_

#ifdef ANSI

void
final_rpt (int    expt,
           int    successes,
           long   running_trials,
           long   highest_trials,
           long   lowest_trials,
           float  solution,
           float  low_soln);
void
stats_bin (int    *trials,
           float  error,
	   float  extra,
           int    *successes,
           long   *running_trials,
           long   *highest_trials,
           long   *lowest_trials,
           float  *solution,
           float  *low_soln,
           float  Cutoff,
	   int    exp);

#else

void final_rpt( );
void stats_bin( );

#endif

#endif
