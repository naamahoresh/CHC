
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
#include <ctype.h>
#include <time.h>
#include "gene.h"
#include "ga_global_extern.h"
#include "ga_status.h"
#include "ga_params.h"
    
/***************************************************************************
 * FUNCTION: parse_command_line
 *
 * DESCRIBE: accepts the main() command line arguments (minus the program
 *           name itself) and sets the appropriate options. 
 *           Command line format is as follows: -x value
 *           where 'x' is a one-character tag and a space must appear between
 *           the tag and associated value.
 *
 * INPUT PARAMETERS: integer;
 *                   array of pointers to character strings
 *
 * RETURN VALUE: none; in case of error, the program is exited
 *
 * CALLS: fatal_error
 *        set_parameter
 ****************************************************************************/
#ifdef ANSI
void
parse_command_line ( int     argc,
 		     char    *argv[])
#else
void
parse_command_line (argc, argv)
int                 argc;
char               *argv[];
#endif
{
 int i;
 char mssg[80];

 if (argc == 0)
    {
     usage();
     fatal_error ("No GA parameters supplied.");
    }

 if ((argc % 2) != 0)
    fatal_error ("Illegal command line syntax.");


 /* Initially, calling set_parameter with a 'z' value of "reset"
  * makes sure all required parameters have old values cleared.
  */
 set_parameter ('z', "reset");

 for (i=0; i< argc; i+=2)
     {
      /* required format of args is "-x val"
       * where x is a recognized tag char,
       * and a space between tag and value required
       */

      if (strlen(argv[i]) != 2 ||
          *argv[i] != '-' ||
          !isalpha(*(argv[i]+1)))
         {
          sprintf (mssg, "Illegal syntax: Command line argument %d.", i+1);
          fatal_error (mssg);
         }
      
       if (!set_parameter (*(argv[i]+1), argv[i+1]))
          {
           sprintf (mssg, 
                    "-%c %s is an unacceptable tag-value pair.", 
                    *(argv[i]+1), argv[i+1]);

           fatal_error (mssg);
          }
      }


 /* Finally, calling set_parameter with a 'z' value of "check"
  * makes sure all required parameters have values
  */
 if (!set_parameter ('z', "check"))
    fatal_error ("All required parameters were not set");
}

/***************************************************************************
 * FUNCTION: parse_config_file
 *
 * DESCRIBE: accepts as input a pointer to a file containing ga parameter
 *           settings.  The file format is as follows:
 *
 *           PoolSize       120
 *           Length_String  25
 *           NumberTrials   25000
 *           RandomSeed     12345678
 *           Bias           1.900000
 *           MutateRate     0.000000
 *           InitPool       dump.pool
 *           FinalPool      final.pool
 *           StatusInterval 1000
 *           DumpInterval   5000
 *           OutfileBase    my_dump
 *           NodeFile       nodes.tsp
 *           UCutoff         any appropriate value
 *           WNumPop        10
 *           XSwapInterval  1000
 *           YSwapNumber    5
 *           Experiments    5
 * 	     File           userdefined file name
 *	    
 *           the following requirements apply:
 *           1) only 1 'tag value' pair per line
 *           2) the tag and its value must appear on the same line,
 *              separated by at least one space
 *           3) the tag can only be one word
 *           4) only the first letter of the tag will be used by this
 *              parser...and it must agree with the command line option
 *              tag.
 *
 *           NOTE: the 'c' tag is disallowed in the config file since
 *                 it could cause infinite recursion
 *
 * INPUT PARAMETERS: FILE ptr
 *
 * RETURN VALUE: 1 on success; 0 on failure 
 *
 * CALLS: fatal_error
 *        warning
 *        set_parameter
 ****************************************************************************/
#ifdef ANSI
int
parse_config_file (FILE   *fp)
#else
int
parse_config_file (fp)
FILE              *fp;
#endif
{
 char  tag[80];
 char  value[80];

 while (fscanf (fp, "%s %s", tag, value) == 2)
       {
        /* avoid infinite recursion at all costs! */
        if (tag[0] == 'c' || tag[0] == 'C')
           {
            warning ("Configuration files may not specify other configuration files.");
            return (0);
           }
        
        else
        if (!set_parameter (tag[0], value))
            return (0);
       }
 fclose (fp);
 return (1);

}

/***************************************************************************
 * FUNCTION: print_params
 *
 * DESCRIBE: Prints the current values of the GA global params
 *           (see ga_params.h)
 *           Note that the output format is compatible with tag-value pairs
 *           expected in config files and used for dumps.
 *
 * INPUT PARAMETERS: none
 *
 * RETURN VALUE: none
 ****************************************************************************/
#ifdef ANSI
void
print_params( FILE   *fp)
#else
void
print_params(fp)
FILE        *fp;
#endif
{
 fprintf(fp, "\n");
 fprintf (fp, "%-19s %d\n",    "PoolSize",          PoolSize);
 fprintf (fp, "%-19s %d\n",    "LengthString",      StringLength);
 fprintf (fp, "%-19s %d\n",    "Trials",            NumberTrials);
 fprintf (fp, "%-19s %d\n",    "GenerationCurrent", CurrentGeneration);
 fprintf (fp, "%-19s %ld\n",   "RandomSeed",        RandomSeed);
 fprintf (fp, "%-19s %f\n",    "Bias",              SelectionBias);
 fprintf (fp, "%-19s %f\n",    "MutateRate",        MutateRate);
 fprintf (fp, "%-19s %s\n",    "OutfileDump",       DumpBase);
 if (DumpInterval > 0)
    fprintf (fp, "%-19s %d\n", "DumpInterval", DumpInterval);
 if (StatusInterval > 0)
    fprintf (fp, "%-19s %d\n", "StatusInterval", StatusInterval);
 if (strlen (SeedPool) > 0)
    fprintf (fp, "%-19s %s\n", "InitPoolFile",  SeedPool);
 if (strlen (FinalPool) > 0)
    fprintf (fp, "%-19s %s\n", "FinalPoolFile", FinalPool);
 if (strlen (NodeFile) > 0)
    fprintf (fp, "%-19s %s\n", "NodeFile", NodeFile);
 if (Experiments > 0)
    fprintf (fp, "%-19s %d\n", "Experiments", Experiments);

 fprintf (fp, "%-19s %f\n",    "CutOff", CutOff);

 if (NumPop > 0) {
    fprintf (fp, "%-19s %d\n", "NumPop", NumPop);
    fprintf (fp, "%-19s %d\n", "SwapInterval", SwapInterval);
    fprintf (fp, "%-19s %d\n", "SwapNumber", SwapNumber);
 if (strlen (UserFile) > 0)
    fprintf (fp, "%-19s %s\n", "UserFile", UserFile);
 }
}

/***************************************************************************
 * FUNCTION: set_parameter
 *
 * DESCRIBE: accepts as input a tag associated with a ga parameter and a value
 *           for that parameter. Sets the parameter and keeps a record which
 *           parameters have been set. 
 *           TO RESET: (such that it thinks NO parameters have been set)
 *                     tag = 'z', value = "reset"
 *           TO CHECK: (make sure all required parameters have been set)
 *                     tag = 'z'  value = "check"
 *                     Note that this function will set defaults when necessary;
 *                     a warning message is printed in this event.
 *
 * INPUT PARAMETERS: character; character string
 *
 * RETURN VALUE: 1 for success; 0 for failure
 *
 * RESULT: the following ga global variables may be set:
 *         SelectionBias, StringLength, MutateRate,
 *         PoolSize, RandomSeed, NumberTrials, SeedPool
 *
 * CALLS: parse_config_file
 ****************************************************************************/
#ifdef ANSI
int
set_parameter ( char    tag,
 		char    value[])
#else
int
set_parameter (tag, value)
char           tag;
char           value[];
#endif
{

 /* these static variables are used to signal that a value has been set */
 static int          pool_size = 0;
 static int      string_length = 0;
 static int      number_trials = 0;
 static int        random_seed = 0;
 static int        select_bias = 0;
 static int             mutate = 0;
 static int          seed_pool = 0; 
 static int    status_interval = 0;
 static int      dump_interval = 0;
 static int          dump_base = 0;
 static int         final_pool = 0;
 static int generation_current = 0;
 static int          node_file = 0;
 static int            num_pop = 0;
 static int      swap_interval = 0;
 static int        swap_number = 0;
 static int        experiments = 0;
 static int             cutoff = 0;
 static int          user_file = 0;

 switch (tag)
    {
     /* a for user specified name of file in UserFile */
     case 'a':
     case 'A':
        strcpy (UserFile, value);
        return (user_file = 1);
        break;

     /* b for bias as in SelectionBias */
     case 'b':
     case 'B':
        if (sscanf (value, "%f", &SelectionBias) == 1)
           return (select_bias = 1);
        break;

     /* c for configuration file */
     case 'c':
     case 'C':
         {
         FILE *fp;

         if (!(fp = fopen (value, "r")))
            {
             char mssg[80];
             sprintf (mssg, "Cannot open configuration file %s.", value);
             warning (mssg);
            }
         else
            {
             if (parse_config_file (fp))
                {
                fclose (fp);
                return (1);
                }
             fclose (fp);
            }
         }
        break;

     /* d for dump interval */
     case 'd':
     case 'D':
        if (sscanf (value, "%d", &DumpInterval) == 1)
           return (dump_interval = 1);
        break;
        
     /* e for number of experiments */
     case 'e':
     case 'E':
        if (sscanf (value, "%d", &Experiments) == 1)
           return (experiments = 1);
        break;
        
     /* f for final population file */
     case 'f':
     case 'F':
        strcpy (FinalPool, value);
        return (final_pool = 1);
        break;

     /* g for generation to begin with */
     case 'g':
     case 'G':
        if (sscanf (value, "%d", &CurrentGeneration) == 1)
          {
           return (generation_current = 1);
          }
        break;

     /* i for init population */
     case 'i':
     case 'I':
        strcpy (SeedPool, value);
        return (seed_pool = 1);
        break;


     /* l for length as in StringLength */
     case 'l':
     case 'L':
        if (sscanf (value, "%d", &StringLength) == 1)
           return (string_length = 1);
        break;

     /* m for mutate as in MutateRate */
     case 'm':
     case 'M':
        if (sscanf (value, "%f", &MutateRate) == 1)
           return (mutate = 1);
        break;

     /* n for nodes as in NodeFile */
     case 'n':
     case 'N':
        strcpy (NodeFile, value);
        return (node_file = 1);
        break;

     /* o for output dump file basename as in DumpBase */
     case 'o':
     case 'O':
        strcpy (DumpBase, value);
        return (dump_base = 1);
        break;

     /* p for pool as in PoolSize */
     case 'P':
     case 'p':
        if (sscanf (value, "%d", &PoolSize) == 1)
           return (pool_size = 1);
        break;

     /* r for random as in RandomSeed */
     case 'R':
     case 'r':
        if (sscanf (value, "%ld", &RandomSeed) == 1)
           return (random_seed = 1);
        break;
      

     /* s for status as in StatusInterval */
     case 'S':
     case 's':
        if (sscanf (value, "%d", &StatusInterval) == 1)
           return (status_interval = 1);
        break;


     /* t for trials as in NumberTrials */
     case 't':
     case 'T':
        if (sscanf (value, "%d", &NumberTrials) == 1)
           return (number_trials = 1);
        break;

     /* u for cutoff value for a run */
     case 'u':
     case 'U':
        if (sscanf (value, "%f", &CutOff) == 1)
           return (cutoff = 1);
        break;

     /* w for number of subpopulations as in NumPop */
     case 'w':
     case 'W':
        if (sscanf (value, "%d", &NumPop) == 1)
           return (num_pop = 1);
        break;

     /* x for swap interval between subpopulations */
     case 'x':
     case 'X':
        if (sscanf (value, "%d", &SwapInterval) == 1)
           return (swap_interval = 1);
        break;

     /* y for number of strings swapped between subpopulations */
     case 'y':
     case 'Y':
        if (sscanf (value, "%d", &SwapNumber) == 1)
           return (swap_number = 1);
        break;

     /* z should never be used as a command line option;
      * it signals administrative actions for this routine.
      */
     case 'z':
     case 'Z':
        if (!strcmp ("reset", value))
           {
                    pool_size = 0;
                string_length = 0;
                number_trials = 0;
                  random_seed = 0;
                  select_bias = 0;
                       mutate = 0;
                    seed_pool = 0;
                   final_pool = 0;
                    dump_base = 0;
              status_interval = 0;
                dump_interval = 0;
           generation_current = 0;
                    node_file = 0;
                      num_pop = 0;
                swap_interval = 0;
                  swap_number = 0;
		  experiments = 0;
		       cutoff = 0;
           }

        else
        if (!strcmp ("check", value) & 
            pool_size &                 /* these 3 values MUST be set */
            string_length &             /* by the user; no defaults allowed */
            number_trials                    
           )
           {
            /* these 3 values CAN be set by the user; 
             * but will default if not set 
             */
            if (!random_seed)
               {
                char rand_seed[80];
                RandomSeed = (long) 123456789;

                sprintf (rand_seed, 
                         "Using (Default) Random Number Seed of %ld\n",
                         RandomSeed);
                warning (rand_seed);
               }
            if (!select_bias)
               {
                warning ("Using (Default)  Selection Bias of 2.0.\n");
                SelectionBias = 2.0;
               }
            if (!mutate)
                warning ("(By Default) Adaptive Mutation Will Not Be Used.\n");
             
             /* this parameter is required for TSP; default is nil */
             if (!node_file)
                NodeFile[0] = '\0';

             /* this parameter is optional for users discretion */
             if (!user_file)
                UserFile[0] = '\0';

             /* this parameter is optional; default is nil */
             if (!seed_pool)
                SeedPool[0]='\0';

             /* this parameter is optional; default is nil */
             if (!final_pool)
                FinalPool[0] ='\0';

             /* this parameter is optional; default is "dump" */
             if (!dump_base)
                strcpy (DumpBase, "dump");

             /* this parameter is optional; default is 0 */
             if (!generation_current)
                CurrentGeneration = 0;

             /* this parameter is optional; default is 0 */
             if (!status_interval)
                StatusInterval = 0;

             /* this parameter is optional; default is 0 */
             if (!dump_interval)
                dump_interval = 0;

             /* this parameter is optional; default is 0 */
	     if (!num_pop)
		NumPop = 0;

             /* this parameter is optional; default is 1 */
	     if (!experiments)
                Experiments = 1;

             /* this parameter is optional; not sure about default-
		remove from main if not using!!! */
	     if (!cutoff)
		CutOff = -9999;

            /* print_params(stderr); */
            return (1);
           }
        break;

     default:
        break;
    }

 return (0);  /* error */
}


/***************************************************************************
 * FUNCTION: usage
 *
 * DESCRIBE: prints a usage message to the screen
 *
 * INPUT PARAMETERS: none
 *
 * RETURN VALUE: none
 ****************************************************************************/
#ifdef ANSI
void
usage ( void )
#else
void
usage ( )
#endif
 {
  fprintf (stderr, "\nUSAGE:\n");
  fprintf (stderr, "\nGenitor [-option value]\n");
  fprintf (stderr, "\n-a file < User file; string >");
  fprintf (stderr, "\n-b 1.25 < Bias selection pressure; float >");
  fprintf (stderr, "\n-c file < Config Filename; string >");
  fprintf (stderr, "\n-d int  < Dump Interval; integer >");
  fprintf (stderr, "\n-e 5    < Number of Experiments; integer>");
  fprintf (stderr, "\n        < variable 'Experiments' must be in main()>");
  fprintf (stderr, "\n-f file < Final pool output dump file; string >");
  fprintf (stderr, "\n-g int  < Current generation to start with; integer >");
  fprintf (stderr, "\n-i file < Initial Population Filename; string >");
  fprintf (stderr, "\n-l 25   < String Length; integer >");
  fprintf (stderr, "\n-m .15  < Mutation Rate; floating point (0-1.0) >");
  fprintf (stderr, "\n-n file < Node file (TSP only); string >");
  fprintf (stderr, "\n-o file < Output Dump File Basename; string >");
  fprintf (stderr, "\n-p 100  < Population Size; integer >");
  fprintf (stderr, "\n-r 8765 < Random Seed; (long) integer >");
  fprintf (stderr, "\n-s 500  < Status Interval; integer >");
  fprintf (stderr, "\n-t 1000 < Number Trials; integer >");
  fprintf (stderr, "\n-u .003 < Cutoff point for a run; float >");
  fprintf (stderr, "\n-w 10   < Number Subpopulations; integer >");
  fprintf (stderr, "\n-x 1000 < Trials Between Swapping; integer >");
  fprintf (stderr, "\n-y 5    < Number Of Strings Swapped; integer >");

  fprintf (stderr, "\n\n");
 }

