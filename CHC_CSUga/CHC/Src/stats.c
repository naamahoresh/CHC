/**************    Binary Generic Functions ***************/


#include <stdio.h>
#include <math.h>
#include <malloc.h>

extern   int    Bin_Tries[];
extern   float  Bin_Solns[];
extern   float  Extra_Solns[];



/**************    Final Report and Averages  *****************/
#ifdef ANSI
void
final_rpt (int    expt,
           int    successes,
           long   running_trials,
           long   highest_trials,
           long   lowest_trials,
           float  solution,
           float  low_soln)
#else
void final_rpt (expt, successes, running_trials, highest_trials, lowest_trials,
                solution, low_soln)
   int    expt;
   int    successes;
   long   running_trials;
   long   highest_trials;
   long   lowest_trials;
   float  solution;
   float  low_soln;
#endif
{
    float avg_trial;
    float avg_soln;
    float avg_extra;
    float low_ex;
    int   k;
    float accum, std_dev, inter;


    /* ************************************************************ */
    /*   Deal with the Solutions First                              */
    /* ************************************************************ */
    avg_soln = solution / (float) expt;

    /* Added to calculate Standard Deviation for Solutions */
    if (successes == expt)
       fprintf(stdout, "\nMean Soln: %e   All Solved", avg_soln);
    else
    {

/*     fprintf(stdout, "\nSOLUTIONS");
       for (k = 0; k < expt; k++)
          fprintf(stdout, "\nSol. %d   %e", k, Bin_Solns[k]);
*/

       accum = 0.0;
       for (k = 0; k < expt; k++)
       {
          inter = Bin_Solns[k] - avg_soln;
          accum += (inter * inter);
       }
       std_dev = (float) sqrt((double)(accum/(expt - 1)));
       fprintf (stdout, "\nMean Soln: %e    Sigma: %e", avg_soln, std_dev);
       fprintf (stdout, "\nLowest:  %e\n", low_soln);
    }


    /* ************************************************************ */
    /*   Now Deal with the Other Data You Want to Track             */
    /* ************************************************************ */
    solution = 0;
    low_ex = Extra_Solns[0];
    for (k = 0; k < expt; k++)
    {
       solution += Extra_Solns[k];
       if (Extra_Solns[k] < low_ex)
	  low_ex = Extra_Solns[k];
    }

/*    fprintf(stdout, "\nEXTRA DATA");
    for (k = 0; k < expt; k++)
       fprintf(stdout, "\nData %d   %f", k, Extra_Solns[k]);
*/
    avg_extra = solution / (float) expt;

    /* Added to calculate Standard Deviation for Extra Data */
    accum = 0.0;
    for (k = 0; k < expt; k++)
    {
       inter = Extra_Solns[k] - avg_extra;
       accum += (inter * inter);
    }
    std_dev = (float) sqrt((double)(accum/(expt - 1)));
    fprintf (stdout, "\nMean Extra (EXPTS): %f    Sigma: %f", 
	     avg_extra, std_dev);
    fprintf (stdout, "\nLowest Extra:  %f\n", low_ex);


    /* ************************************************************ */
    /*   Now deal with the trials                                   */
    /* ************************************************************ */
    if (successes != 0)
       avg_trial = (float) ( (float)running_trials / (float)successes);
    else 
       avg_trial = 0.0;

    /* Added to calculate Standard Deviation for the Trials to Solve */
    if (avg_trial > 0.0)
    {
       accum = 0.0;
       for (k = 0; k < successes; k++)
       {
          inter = Bin_Tries[k] - avg_trial;
          accum += (inter * inter);
       }
       std_dev = (float) sqrt((double)(accum/(successes - 1)));
       fprintf (stdout, "\nMean Trials: %f    Sigma: %f", avg_trial, std_dev);
       fprintf (stdout, "\nHighest: %ld   Lowest: %ld", 
		highest_trials, lowest_trials);
    }
    else
       fprintf (stdout, "\nMean Trials: %f    None Solved", avg_trial);


    fprintf (stdout, "\n\n *** THERE WERE %d SUCCESSES.\n\n", successes);
}
/*****************************************************/
/* end of final_rpt()                                */
/*****************************************************/



/**************    Final Report Statistics Update   *****************/
#ifdef ANSI
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
           float  CutOff,
	   int    exp)
#else
void stats_bin (trials, error, extra, successes, running_trials, 
                highest_trials, lowest_trials, solution, low_soln, 
		CutOff, exp)
   int    *trials;
   float  error;
   float  extra;
   int    *successes;
   long   *running_trials;
   long   *highest_trials;
   long   *lowest_trials;
   float  *solution;
   float  *low_soln;
   float  CutOff;
   int    exp;
#endif
{


  *solution += error;
  Bin_Solns[exp] = error;
  Extra_Solns[exp] = extra;

  if (error < *low_soln)
     *low_soln = error;

  if (error <= CutOff)
  {
     Bin_Tries[(*successes)] = (*trials);
     (*successes)++;
     *running_trials += (long)(*trials);
     if ((long)(*trials) > *highest_trials)
        *highest_trials = (long)(*trials);
     if ((long)(*trials) < *lowest_trials)
        *lowest_trials = (long)(*trials);
  }

  *trials = 0;

}
/*****************************************************/
/* end of stats_bin()                                */
/*****************************************************/

