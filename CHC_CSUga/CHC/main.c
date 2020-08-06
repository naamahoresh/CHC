/*  **********************  CHC Model ********************* */

#include <stdio.h>
#include <ctype.h>

#include "ga_random.h"
#include "gene.h"
#include "ga_global.h"
#include "ga_params.h"
#include "ga_pool.h"
#include "ga_copy.h"
#include "ga_selection.h"
#include "ga_status.h"
#include "ga_signals.h"


/*******************************************
 Declare evaluation function which tells you
 how "good" a particular gene's solution is.
 NOTE: the input parameters for the eval
       function should always be a gene and
       the gene's length
 *******************************************/
#include "binarys.h"
#include "stats.h"
#include "m_evals.h"
#include "op_hux.h"
#include "op_cataclysm.h"
#include "binarys_params.h"
#include "stats_params.h"
#include "op_hux_params.h"



/**************************************/
/* main program for running CHC GA    */
/**************************************/
int main (argc, argv)
int   argc;
char *argv[];
{
 int        k, m;
 int        experiment_cnt;
 int        threshold, org_threshold, eval_count, produced, insert;
 
 GENEPTR    mom, dad;
 POOLPTR    *parent_array, *child_array;

 int        subpop = 0;

 setup_signal();                 /* Setup signal handlers */

/*************************************************************
 Set the global parameters according to command line arguments.
 *************************************************************/
 argc--;    argv++;
 parse_command_line (argc, argv);

 srandom(RandomSeed);          /*  Seed the Random Number Generator */
/******************************************************
 Allocate a parent and child genetic pool 
 ******************************************************/
 if (!(parent_array = (POOLPTR *)malloc(sizeof(POOLPTR) * NumPop)))
    fatal_error(NULL);
 if (!(child_array = (POOLPTR *)malloc(sizeof(POOLPTR) * NumPop)))
    fatal_error(NULL);
    
 for (k=0; k<NumPop; k++) 
    if ( !(parent_array[k] = get_pool(PoolSize, StringLength)) )
       fatal_error(NULL);

 for (k=0; k<NumPop; k++) 
    if ( !(child_array[k] = get_pool(PoolSize, StringLength)) )
       fatal_error(NULL);


/******************************************************
 PUT IN ANY INITIALIZATION CALLS HERE
  ******************************************************/
  dif_ary = allocate_hux_params(StringLength);
  prblm_init_dejong(NodeFile, &Nb_params, &Params, &Work_params,
		    StringLength, &Grey_gene);

 fprintf (stdout, "\n");  
 SelectionBias = 1.0;            /* Make sure of Random selection */
 /* StatusInterval = 1;             Report every generation */
 print_params(stdout);           /* Print Parameter Values */
 print_prog_head(stdout);

/*******************************************************
 Allocate temporary storage for parents of reproduction.
 *******************************************************/
 mom = get_gene (StringLength);
 dad = get_gene (StringLength);


/* MAIN LOOP FOR EXPERIMENTS */
for (experiment_cnt=0; experiment_cnt<Experiments; experiment_cnt++) 
{
    if (experiment_cnt)
       CurrentGeneration = 0;  
 
    for (k=0; k<NumPop; k++)             /* Initialize parental genetic pool */
         init_pool (SeedPool, parent_array[k], 0, 
		    parent_array[k]->size, mstr_eval); 

    for (k=0; k<NumPop; k++)             /* sort the pool - easier this way */
        sort_pool(parent_array[k]);

    org_threshold = threshold = (int)(StringLength/4);
    eval_count = 0;    

     /********** Optimize ! **********/
     for (/* CurrentGeneration either intialized to 0 in declaration OR
             initialized by restart of previous experiment */;
             CurrentGeneration < NumberTrials; 
             CurrentGeneration++, cur_trials++)
     {
 	 produced = 0;

         if (StatusInterval && !(CurrentGeneration % StatusInterval)) 
	  show_progress_exp (stdout, parent_array[subpop], eval_count);

	 for (m=0; m<PoolSize; m+=2)
	 {
	    get_parents(mom, dad, parent_array[subpop], SelectionBias, linear);
            
	    if ((int)(hamming_distance(mom->string,dad->string,StringLength)/2) 
		    > threshold)
	    {
	       hux(mom->string, dad->string, dif_ary, StringLength);
	       mom->worth = mstr_eval(mom->string, StringLength);
	       dad->worth = mstr_eval(dad->string, StringLength);
	       gene_copy(&(child_array[subpop]->data[produced]), mom, 
			 StringLength);
	       gene_copy(&(child_array[subpop]->data[produced+1]), dad, 
			 StringLength);
	       produced += 2;    eval_count += 2;
	    }
	 }
	
	 insert = 0;
	 for (m=0; m<produced; m++)
	 {
	    if (child_array[subpop]->data[m].worth < 
		parent_array[subpop]->data[PoolSize-1].worth)
            {
	     insert++;
	     insert_gene(&(child_array[subpop]->data[m]), parent_array[subpop]);
	    }
	 }

         if (parent_array[subpop]->data[0].worth <= CutOff)
	    break;

	 if(!insert)
	 {
	    threshold--;
	    if (threshold < 0)
	    {
	       cataclysm(parent_array[subpop], PoolSize-1, StringLength, 
			 0, MutateRate);
               threshold = org_threshold;

	       for (m=0; m<PoolSize; m++)
	       {
	          parent_array[subpop]->data[m].worth = 
		    mstr_eval(parent_array[subpop]->data[m].string, StringLength);
                  eval_count++;
	       }
    
	       for (k=0; k<NumPop; k++)               /* sort the pool again */
                   sort_pool(parent_array[k]);
	    }
	 }

     if (eval_count >= NumberTrials)
     {
	show_progress_exp (stdout, parent_array[subpop], eval_count);
	break;
     }

    /*****************************************************************
     If the DumpInterval parameter was set and this is the appropriate
     time, save the population and key parameters to disk for later
     reference (or to restart execution later.
     *****************************************************************/
     if (DumpInterval && !(CurrentGeneration % DumpInterval))
        dump_status(parent_array[subpop], DumpBase);

    }  /* End of WHILE NUMBER of Trials/Generations in an experiment */
	  
	  
    show_progress_exp (stdout, parent_array[subpop], eval_count);
	    
 
/*****************
 Summarize Results
 *****************/
  stats_bin (&eval_count, parent_array[subpop]->data[0].worth, 
	     parent_array[subpop]->data[0].worth, &successes, 
	     &running_trials, &highest_trials, &lowest_trials, 
	     &solution, &low_soln, CutOff, experiment_cnt);

  mstr_rpt(parent_array[0]->data[0].string);

 }  /* end of FOR NUMBER of Experiments loop */

  final_rpt(Experiments, successes, running_trials, highest_trials, 
	    lowest_trials, solution, low_soln);

} /* End of Main Program */
