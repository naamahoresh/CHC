/**************     DeJong Evaluation Functions    ***************/
/*                    and other Eval Functions     ***************/
#include <stdio.h>
//#include <stdlib.h> 
#include <math.h>
#include <string.h>

#include "gene.h" 
#include "ga_random.h" 
#include "binarys.h"
#include "m_evals.h"

 /* ******************************************************************** */
 /* Global Parameters needed by the Test routines but defined elsewhere  */
 /* ******************************************************************** */
  extern int Nb_params;                  /* number of parameters in the run */
  extern struct p_list *Work_params;     /* Parameters for run              */
  extern double *Params;                 /* Work space for binary params    */
  extern char *Grey_gene;                /* Work space for de-graying a gene */
 /* **************************************************************** */
 /* **************************************************************** */
 
 /* ******************************************************************** */
 /*                Defines Unique to this file                           */
 /*        These defines describe the various parameter options          */
 /* available to drive experiments in the .params file                   */
 /* ******************************************************************** */
#define GREY          'g'
#define BCD           'b'
#define FUNCTION      'f'
#define COMPOSE       'c'
#define DUAL          'd'
#define DISTRIBUTED   'd'
#define MULTIPLICTIV  'm'
#define ADDITIV       'a'
#define SIMPLE        's'
#define FULL          'f'
#define LOWERLEFT     'l'
#define WRAP          'w'
#define STANDARD      's'

#define TERMINAL      't'
#define INNER         'i'
#define OUTTER        'o'
#define BOTH          'b'
#define ALL           'a'

#define XTERMINAL     0x01
#define XINNER        0x02
#define XOUTTER       0x04
#define XBOTH         0x06
#define XALL          0x07


 /* ******************************************************************** */
 /*                                Macros                                */
 /* ******************************************************************** */

#define LOG2(a)        log(a)/log(2.0)

 /* ******************************************************************** */
 /*                Global Parameters Unique to this file                 */
 /* ******************************************************************** */

  static char   Coding[20], Cplxty[20], Cform[20];
  static int    Func, Cfunc1, Cfunc2;
  static int    NoiseMask;
  static float  RealFit;

#include "dejong_params.h"

/*
using namespace std;
int myproc;
extern "C"
{
    void shgval_(int*,int*,int*,double *, double[]);
}
*/
int myproc;
void shgval_(int*,int*,int*,double *, double[]);

#ifdef ANSI
float get_tru_fit(void)
#else
float get_tru_fit()
#endif
{
   return (RealFit);
};

#ifdef ANSI
void set_tru_fit(float  val)
#else
void
set_tru_fit(val)
   float   val;
#endif
{
   RealFit = val;
};





/* ***************************************************************** */
/* Problem initialization routine                                    */
/* ***************************************************************** */
#ifdef ANSI
void
prblm_init_dejong(char   *NodeFile,
              int    *Nb_params,
              double **Params,
              struct p_list **Work_params,
	      int Length,
	      char **gray_wrk)
#else
void prblm_init_dejong(NodeFile, Nb_params, Params, Work_params, Length,
		       gray_wrk)
   char  *NodeFile;
   int  *Nb_params;
   struct p_list **Work_params;
   double **Params;
   int Length;
   char **gray_wrk;
#endif
{
  FILE  *fopen(), *parmfil;
  int i, postn, count;
  char ParmFil[50];
  char sign[2], par[10], note[20];


  strcpy(ParmFil, NodeFile);         /* open file with delta configuration */
  strcat(ParmFil, ".params");
  if ( (parmfil = fopen(ParmFil, "r")) == NULL){
     fprintf(stdout, "Could not open parameter configuration file %s\n",
                      NodeFile);
     exit(1); }

  
  count = fscanf (parmfil, "%1c\n", par);
  if ((par[0] == 'n') || (par[0] == 'N'))
  {
     count = fscanf (parmfil, "%1c\n", par);
     switch (par[0])
     {
        case TERMINAL:
	case 'T':
	   NoiseMask = XTERMINAL;
           fprintf(stdout, "Terminal noise has been selected. This is only valid for function composition of the form f8(f2) + f8(f2) + f8(f2)\n");
	   break;

        case INNER:
	case 'I':
	   NoiseMask = XINNER;
           fprintf(stdout, "The INNER noise option has been selected.\n");
	   break;

        case OUTTER:
	case 'O':
	   NoiseMask = XOUTTER;
           fprintf(stdout, "The OUTTER noise option has been selected.\n");
	   break;

        case BOTH:
	case 'B':
	   NoiseMask = XBOTH;
           fprintf(stdout, "Both the INNER and OUTTER noise options have been selected.\n");
	   break;

        case ALL:
	case 'A':
	   NoiseMask = XALL;
           fprintf(stdout, "All noise options have been selected (INNER, OUTTER, and TERMINAL). This is only valid for function composition of the form f8(f2) + f8(f2) + f8(f2)\n");
	   break;

        default:
           fprintf(stdout, "Invalid Noise Option!\n");
           exit(7);
           break;
     }
     count = fscanf (parmfil, "%1c\n", par);
  }
  else
     NoiseMask = 0; 

  switch (par[0])
  {
     case COMPOSE:
     case 'C':
        count = fscanf (parmfil, "%d %d %s\n", &Cfunc1, &Cfunc2, Cform);
	Func = 0;               /* function composition indicated by type 0 */
        
	switch (Cform[0])
	{
	   case ADDITIV:
	   case 'A':
	      Cform[0] = ADDITIV;
              fprintf(stdout, "Composition: f%d(f%d + f%d + f%d)\n", 
	              Cfunc1, Cfunc2, Cfunc2, Cfunc2);
              break;

	   case DISTRIBUTED:
	   case 'D':
	      Cform[0] = DISTRIBUTED;
              fprintf(stdout, "Composition: f%d(f%d, f%d, f%d)\n", 
	              Cfunc1, Cfunc2, Cfunc2, Cfunc2);
              break;

	   case MULTIPLICTIV:
	   case 'M':
	      Cform[0] = MULTIPLICTIV;
              fprintf(stdout, "Composition: f%d(f%d) + f%d(f%d) + f%d(f%d)\n", 
	              Cfunc1, Cfunc2, Cfunc1, Cfunc2, Cfunc1, Cfunc2);
              break;

	   case SIMPLE:
	   case 'S':
	      Cform[0] = SIMPLE;
              fprintf(stdout, "Composition: f%d(f%d)\n", Cfunc1, Cfunc2);
              break;

           default:
              fprintf(stdout, "Invalid composition!\n");
              exit(7);
	      break;
        }
	break;
     
     case DUAL:
     case 'D':
        count = fscanf (parmfil, "%d %d\n", &Cfunc1, &Cfunc2);
	Func = Cfunc1;                 /* DUAL function composition - Unique */
        fprintf(stdout, "DUAL Function Composition: f%d(f%d)\n", 
	        Cfunc1, Cfunc2);
	break;

     case FUNCTION:
     case 'F':
        count = fscanf (parmfil, "%d\n", &Func);
        fprintf(stdout, "Test function selected: f%d\n", Func);
	break;
     
     default:
        fprintf(stdout, "Invalid function type!\n");
        exit(7);
	break;
  }
  count = fscanf (parmfil, "%s\n", Coding);
  count = fscanf (parmfil, "%s\n", Cplxty);


  switch (Coding[0])
  {
     case BCD:
     case 'B':
        Coding[0] = BCD;
        fprintf(stdout, "Encoding Selected: BCD\n");
	break;
  
     case GREY:
     case 'G':
        Coding[0] = GREY;
        fprintf(stdout, "Encoding Selected: Grey\n");
	break;
  
     default:
        fprintf(stdout, "Invalid string encoding!\n");
        exit(7);
	break;
  }
  
  switch (Cplxty[0])
  {
     case FULL:
     case 'F':
        Cplxty[0] = FULL;
        fprintf(stdout, "Complexity Selected: Full\n");
	break;
     
     case LOWERLEFT:
     case 'L':
        Cplxty[0] = LOWERLEFT;
        fprintf(stdout, "Complexity Selected: LowerLeft\n");
	break;
  
     case STANDARD:
     case 'S':
     case ' ':
        Cplxty[0] = STANDARD;
        fprintf(stdout, "Standard function selected. No Complexity Option.\n");
	break;
  
     case WRAP:
     case 'W':
        Cplxty[0] = WRAP;
        fprintf(stdout, "Complexity Selected: Wrap\n");
	break;
  
     default:
        fprintf(stdout, "Invalid Function Complexity!\n");
        exit(7);
	break;
  }

  if (!Func)     /* then we are doing composition */
  {
      if (Cform[0] == SIMPLE)
      {
	 if (Cplxty[0] != STANDARD)
	 {
	    fprintf(stdout, "Invalid function composition/complexty combination.\n");
	    exit(7);
	 }
      }
  }

  count = fscanf (parmfil, "%d\n", Nb_params);
  fprintf(stdout, "Number of parameters: %d\n",(*Nb_params));

  allocate_bin_params(Params, Work_params, (*Nb_params), Length, gray_wrk);

  /* the syntax for addressing the items inside the array of structure
     pointers Work_params inside this routine is (*Work_params)[index].
     This is because of the **Work_params definition inside this routine. */
  postn = 0;
  for (i = 0; i < *Nb_params; i++)
  {
     count = fscanf(parmfil, "%s %d %c %lf %lf %lf %s\n",
                    par, &((*Work_params)[i].lngth), sign,
                    &((*Work_params)[i].scale), &((*Work_params)[i].wrap_low),
                    &((*Work_params)[i].wrap_high), note);


     (*Work_params)[i].scale = 1.0 / (*Work_params)[i].scale;
     (*Work_params)[i].posn = postn;
     postn += (*Work_params)[i].lngth;

     switch (sign[0])
     {
        case 's' :
        case 'S' :
                   (*Work_params)[i].sign = SIGNED;
                   break;
  
        case 'u' :
        case 'U' :
                   (*Work_params)[i].sign = UNSIGNED;
                   break;

     }  /* end of sign switch */
  
     fprintf(stdout, "P:%d Length: %d\n",i, (*Work_params)[i].lngth);

  }  /* end of parameter config reading */

  fclose(parmfil);
  fflush(stdout);

}
/***********************************/
/* end of prblm_init_dejong()      */
/***********************************/

/* ***************************************************************** */
/* Evaluation initialization routine for raw_eval                    */
/* ***************************************************************** */
#ifdef ANSI
void
raw_init_dejong(char   *NodeFile,
              int    *Nb_params,
              double **Params,
              struct p_list **Work_params,
	      int Length, FILE **parmfil, 
              int *n_degray, char *val_mode)
#else
     void raw_init_dejong(NodeFile, Nb_params, Params, Work_params, Length,
		     parmfil, n_degray, val_mode)
   char  *NodeFile;
   int  *Nb_params;
   struct p_list **Work_params;
   double **Params;
   int Length;
   FILE **parmfil;
   int *n_degray;
   char *val_mode;
#endif
{
  FILE  *fopen();
  int i, postn, count;
  char ParmFil[50];
  char sign[2], par[10], note[20];
  char *gray_wrk;
/**
  NOTE:  The input file format for raw_eval is identical to the
         mstr_eval params file except that there will be data 
	 following the params and the noise options are not implemented
	 in the raw_eval function.  For this reason, the initialization
	 routine is identical to prblm_init_dejong except that it 
	 does not close the input file (which should be input to the
	 raw_eval routine
**/

  strcpy(ParmFil, NodeFile);         /* open file with delta configuration */
  strcat(ParmFil, ".rparams");
  if ( (*parmfil = fopen(ParmFil, "r")) == NULL){
     fprintf(stdout, "Could not open parameter configuration file %s\n",
                      NodeFile);
     exit(1); }

  /* these were added to the .rparams file to read to generalize the
     function evaluation selection process for a single executable   */
  count = fscanf (*parmfil, "%1c\n", par);
  switch (par[0])
  {
     case COMPOSE:
     case 'C':
        count = fscanf (*parmfil, "%d %d %s\n", &Cfunc1, &Cfunc2, Cform);
	Func = 0;               /* function composition indicated by type 0 */
        
	switch (Cform[0])
	{
	   case ADDITIV:
	   case 'A':
	      Cform[0] = ADDITIV;
              fprintf(stdout, "Composition: f%d(f%d + f%d + f%d)\n", 
	              Cfunc1, Cfunc2, Cfunc2, Cfunc2);
              break;

	   case DISTRIBUTED:
	   case 'D':
	      Cform[0] = DISTRIBUTED;
              fprintf(stdout, "Composition: f%d(f%d, f%d, f%d)\n", 
	              Cfunc1, Cfunc2, Cfunc2, Cfunc2);
              break;

	   case MULTIPLICTIV:
	   case 'M':
	      Cform[0] = MULTIPLICTIV;
              fprintf(stdout, "Composition: f%d(f%d) + f%d(f%d) + f%d(f%d)\n", 
	              Cfunc1, Cfunc2, Cfunc1, Cfunc2, Cfunc1, Cfunc2);
              break;

	   case SIMPLE:
	   case 'S':
	      Cform[0] = SIMPLE;
              fprintf(stdout, "Composition: f%d(f%d)\n", Cfunc1, Cfunc2);
              break;

           default:
              fprintf(stdout, "Invalid composition!\n");
              exit(7);
	      break;
        }
	break;
     
     case DUAL:
     case 'D':
        count = fscanf (*parmfil, "%d %d\n", &Cfunc1, &Cfunc2);
	Func = Cfunc1;                 /* DUAL function composition - Unique */
        fprintf(stdout, "DUAL Function Composition: f%d(f%d)\n", 
	        Cfunc1, Cfunc2);
	break;

     case FUNCTION:
     case 'F':
        count = fscanf (*parmfil, "%d \n", &Func);
	fprintf(stdout, "Test function selected: f%d\n", Func);
	break;
     
     default:
        fprintf(stdout, "Invalid function type!\n");
        exit(7);
	break;
  }

  count = fscanf (*parmfil, "%s\n", Coding);
  count = fscanf (*parmfil, "%s\n", Cplxty);

  switch (Coding[0])
  {
     case BCD:
     case 'B':
        Coding[0] = BCD;
        fprintf(stdout, "Encoding Selected: BCD\n");
	break;
  
     case GREY:
     case 'G':
        Coding[0] = GREY;
        fprintf(stdout, "Encoding Selected: Grey\n");
	break;
  
     default:
        fprintf(stdout, "Invalid string encoding!\n");
        exit(7);
	break;
  }

  switch (Cplxty[0])
  {
     case FULL:
     case 'F':
        Cplxty[0] = FULL;
        fprintf(stdout, "Complexity Selected: Full\n");
	break;
     
     case LOWERLEFT:
     case 'L':
        Cplxty[0] = LOWERLEFT;
        fprintf(stdout, "Complexity Selected: LowerLeft\n");
	break;
  
     case STANDARD:
     case 'S':
     case ' ':
        Cplxty[0] = STANDARD;
        fprintf(stdout, "Standard function selected. No Complexity Option.\n");
	break;
  
     case WRAP:
     case 'W':
        Cplxty[0] = WRAP;
        fprintf(stdout, "Complexity Selected: Wrap\n");
	break;
  
     default:
        fprintf(stdout, "Invalid Function Complexity!\n");
        exit(7);
	break;
  }

  if (!Func)     /* then we are doing composition */
  {
      if (Cform[0] == SIMPLE)
      {
	 if (Cplxty[0] != STANDARD)
	 {
	    fprintf(stdout, "Invalid function composition/complexty combination.\n");
	    exit(7);
	 }
      }
  }

  count = fscanf (*parmfil, "%d\n", Nb_params);
  fprintf(stdout, "Number of parameters: %d\n",(*Nb_params));

  allocate_bin_params(Params, Work_params, (*Nb_params), Length, &gray_wrk);

  /* the syntax for addressing the items inside the array of structure
     pointers Work_params inside this routine is (*Work_params)[index].
     This is because of the **Work_params definition inside this routine. */
  postn = 0;
  for (i = 0; i < *Nb_params; i++)
  {
     count = fscanf(*parmfil, "%s %d %c %lf %lf %lf %s\n",
                    par, &((*Work_params)[i].lngth), sign,
                    &((*Work_params)[i].scale), &((*Work_params)[i].wrap_low),
                    &((*Work_params)[i].wrap_high), note);


     (*Work_params)[i].scale = 1.0 / (*Work_params)[i].scale;
     (*Work_params)[i].posn = postn;
     postn += (*Work_params)[i].lngth;

     switch (sign[0])
     {
        case 's' :
        case 'S' :
                   (*Work_params)[i].sign = SIGNED;
                   break;
  
        case 'u' :
        case 'U' :
                   (*Work_params)[i].sign = UNSIGNED;
                   break;

     }  /* end of sign switch */
  
     fprintf(stdout, "P:%d Length: %d\n",i, (*Work_params)[i].lngth);

  }  /* end of parameter config reading */

  fflush(stdout);

  fscanf(*parmfil,"%s\n", val_mode); 

  switch (val_mode[0])
  {
  case 'b':
  case 'B':
    fprintf(stdout,"Values in binary format\n");
    break;
  case 'f':
  case 'F':
    fprintf(stdout,"Values in floating point format\n");
    break;
  default:
    fprintf(stdout, "Invalid Value Format!\n");
    exit(7);
    break;
  }

  if (Func==97)
   {
      fscanf(*parmfil,"%d\n",n_degray);
      fprintf(stdout,"# Degrays %d\n",*n_degray);
   }

  /* The file pointer for parmfil should be located at the data section
     of the initialization file */

}
/***********************************/
/* end of raw_init_dejong()      */
/***********************************/

/* ********************************************* */
/* Symmetric Parameter Transformation Code       */
/* ********************************************* */
#ifdef ANSI
void  
sym_trnsfrm( GENE_DATA *buff,
	     double *Vals,
             int Nb_params,                  
             struct p_list *Work_params)    
#else
void 
sym_trnsfrm(buff, Vals, Nb_params, Work_params) 
  GENE_DATA *buff;
  double    *Vals;               
  int       Nb_params;
  struct p_list *Work_params;    
#endif
{

  int    i;
  float  max;

  /* decode the bit patterns for the parameters */
  for (i = 0; i < Nb_params; i++)
  {
       /* This change will only screw things up that decoded things
	  piecemeal on their own by depending on the addresses/offsets
          of the buff and Work_params[].posn to somehow sync up.    */

/*     This is the way of the old code which probably should not have
       been done this way.  Rather buff is a pointer that should be 
       adjusted each time through the loop rather than depending on 
       absolute offsets as this will make things relative each time.
       Resulting in much more flexibility.
       Vals[i] = ctod(&buff[Work_params[i].posn], Work_params[i].lngth,
                    Work_params[i].scale, Work_params[i].sign);
*/
     Vals[i] = ctod(buff, Work_params[i].lngth, Work_params[i].scale, 
		    Work_params[i].sign);
     buff += Work_params[i].lngth;
     max = (float) pow(2, (double) Work_params[i].lngth);
     Vals[i] = ((Vals[i] * Work_params[i].wrap_high) / (max / 2.0)) - 
	       Work_params[i].wrap_high;
  }

}
/*********************************/
/* end of sym_transform()        */
/*********************************/

/* ********************************************* */
/* Generic Parameter Transformation Code       */
/* ********************************************* */
#ifdef ANSI
void  
gen_trnsfrm( GENE_DATA *buff,
	     double *Vals,
             int Nb_params,                  
             struct p_list *Work_params)    
#else
void 
gen_trnsfrm(buff, Vals, Nb_params, Work_params) 
  GENE_DATA *buff;
  double    *Vals;               
  int       Nb_params;
  struct p_list *Work_params;    
#endif
{

  int    i;
  double  max,temp;

  static  char vals_set=0;
  static  double scale = 1.0;

  if (!vals_set) /* This assumes that all parameters will have the same range and same #bits  */
    {            /* should be made to compute the scaling factor independently if that changes*/

      max = (double) pow(2, (double) Work_params[0].lngth);
      temp = (Work_params[0].wrap_high - Work_params[0].wrap_low);
      temp = (temp + temp/max)/max;

      while (temp < 1.0)
	{
	  temp *= 10.0;
	  scale /= 10.0;
	}

      if (temp >= 5.0) /* temp should actually be 1 or 9, but 5 handles all possible cases */
	scale*= 10.0;  /* if we need to round up, multiply scale */
      
      vals_set=1;
    }

  /* decode the bit patterns for the parameters */
  for (i = 0; i < Nb_params; i++)
  {
     Vals[i] = ctod(buff, Work_params[i].lngth, Work_params[i].scale, 
		    Work_params[i].sign);
     buff += Work_params[i].lngth;
     Vals[i] = ((Vals[i]*scale)+ Work_params[i].wrap_low);
  }

}
/*********************************/
/* end of gen_transform()        */
/*********************************/


/* **************************************** */
/* Master Function Translation Code         */
/* **************************************** */
#ifdef ANSI
float  
mstr_eval( GENE_DATA *buff,
	   int       length)
#else
float 
mstr_eval(buff, length)
  GENE_DATA *buff;
  int       length;
#endif
{

  GENE_DATA *nw_buf;
  float     sum;
  float     (*eval_func)();
  static int  pow2_length=0;                  /* nearest power of 2 of length */

  if (Coding[0] == GREY)                      /* if this a grey coded string */
  {
      degray(buff, Grey_gene, Nb_params, Work_params);
      nw_buf = Grey_gene;
  }
  else
     nw_buf = buff;
 
  /* transform the binary to input parameters */
  sym_trnsfrm(nw_buf, Params, Nb_params, Work_params);


  if (Func)
     f_fx(Func, &eval_func);    /* look up the address of the function */
  
  switch (Func)
  {
     case 0:                              /* handle function de-composition */
	sum = compos_fx(Params, Nb_params, Cfunc1, Cfunc2, Cplxty[0], Cform[0]);
	break;

     case 1:
     case 100:
     case 3:
     case 4:
     case 5:
     case 6:
     case 7:
     case 20:
     case 21:
     case 22:
	sum = (*eval_func)(Params, Nb_params);
        if (NoiseMask & XOUTTER)
        {
           set_tru_fit(sum);
           sum += gaussian();
        }
        break;
     
     case 2:
     case 8:
     case 9:
     case 10:
	sum = cplx_fx(Params, Nb_params, eval_func, Cplxty[0]);
        break;

     case 98:                    /* handle general function de-composition */
        if (pow2_length == 0) 
	 {
	   pow2_length = (int)pow(2.0,(double)
				  ((int)(LOG2((double)Work_params[0].lngth))));
	   
	   if (pow2_length < (int)Work_params[0].lngth) pow2_length*=2;

	 }
        
	sum = (*eval_func)(nw_buf, Params, Nb_params, Cfunc2, pow2_length);
	break;
     
     case 99:                    /* handle speical function de-composition */
	sum = (*eval_func)(nw_buf, Params, Nb_params, Cfunc2);
	break;

     default: 
	printf("Should never reach default statement of mstr_eval\n");
	exit(99);
  }

  return(sum);
}
/**************************************************** */
/* end of mstr_eval()                                 */
/**************************************************** */


/* **************************************** */
/* Short Function Translation Code          */
/* **************************************** */
#ifdef ANSI
float  
shrt_eval( GENE_DATA *buff,
	   int       length,
	   int       opt)
#else
float 
shrt_eval(buff, length, opt)
  GENE_DATA *buff;
  int       length;
  int       opt;
#endif
{

  GENE_DATA *nw_buf;
  float     sum;
  float     (*eval_func)();


  if (Coding[0] == GREY)                      /* if this a grey coded string */
  {
      if (opt == -1)
         degray(buff, Grey_gene, Nb_params, Work_params);
      else
	 degray(&buff[Work_params[opt].posn], &Grey_gene[Work_params[opt].posn],
	        1, &Work_params[opt]);

      nw_buf = Grey_gene;
  }
  else
     nw_buf = buff;


  if (opt == -1)
     sym_trnsfrm(nw_buf, Params, Nb_params, Work_params);
  else
     sym_trnsfrm(&nw_buf[Work_params[opt].posn], &Params[opt], 
		 1, &Work_params[opt]);


  if (Func)
     f_fx(Func, &eval_func);

  
  switch (Func)
  {
     case 0:                              /* handle function de-composition */
	sum = compos_fx(Params, Nb_params, Cfunc1, Cfunc2, Cplxty[0], Cform[0]);
	break;

     case 1:
     case 100:
     case 3:
     case 4:
     case 5:
     case 6:
     case 7:
     case 20:
     case 21:
     case 22:
	sum = (*eval_func)(Params, Nb_params);
        if (NoiseMask & XOUTTER)
        {
           set_tru_fit(sum);
           sum += gaussian();
        }
        break;
     
     case 2:
     case 8:
     case 9:
     case 10:
	sum = cplx_fx(Params, Nb_params, eval_func, Cplxty[0]);
        break;

/* NOTE: this short evaluation strategy won't work with f98 and f99 */

     default: 
	printf("Should never reach default statement of shrt_eval\n");
	exit(99);
     
  }

  return(sum);
}
/*****************************************************/
/* end of shrt_eval()                                 */
/*****************************************************/



/**************************************************** */
/* Function to generate binary string from an int     */
/**************************************************** */
#ifdef ANSI
static void 
  FillBuff(GENE_DATA *str, int value, int length)
#else
static void 
FillBuff(str, value, length)
     GENE_DATA *str;
     int  value;
     int  length;
#endif
{
   int  tmp;
   int  i;

   for (i=length; i>0; i--)   /* Encode value as a bitstring */
    {
       tmp= (int)pow(2.0, (double)(i-1));
       if (value < tmp)
	 str[length-i] = '0';
       else
	{
	   value-=tmp;
	   str[length-i] = '1';
	}
    }
}	  
/**************************************************** */
/* end of FillBuff()                                 */
/**************************************************** */
   

/**************************************************** */
/* Construct a binary string from a double            */
/**************************************************** */
#ifdef ANSI
static void
  BuildBinString(GENE_DATA *strbuf, double value)
#else
static void
  BuildBinString(strbuf, value)
     GENE_DATA *strbuf;
     float value;
#endif
{
   int ivalue = 0;
   int ioffset;
   int bufflen=Work_params[0].lngth;
   float round = (value < 0.0? -0.9 : 0.9);
   float temp, max;

   static  char vals_set=0;
   static  double scale = 1.0;
   
   if (!vals_set) /* This assumes that all parameters have the same range 
		     and same #bits  */
    {            /* should be made to compute the scaling factor 
		    independently if that changes*/
       
       max = (double) pow(2, (double) Work_params[0].lngth);
       temp = (Work_params[0].wrap_high - Work_params[0].wrap_low);
       temp = (temp + temp/max)/max;
       
      while (temp < 1.0)
       {
	  temp *= 10.0;
	  scale /= 10.0;
       }
       
       if (temp >= 5.0) /* temp should actually be 1 or 9, but 5 
			   handles all possible cases */
	 scale*= 10.0;  /* if we need to round up, multiply scale */
       
       vals_set=1;
    }
   
   ivalue = (int)(value/scale + round);
   
   if (Work_params[0].sign == SIGNED)
    {
       /* signed representation */
       if (value < 0) 
	{
	   strbuf[0]='1';
	   ivalue*=(-1);
	}
       else    
	 strbuf[0]='0';
       
       FillBuff(&(strbuf[1]), ivalue, bufflen-1);
    }
   else
    {  /* unsigned represention with shift */
       round = (Work_params[0].wrap_low < 0.0? -0.9 : 0.9);
       ioffset = (int)-1*(Work_params[0].wrap_low/scale+round);/* add offset */
       FillBuff(strbuf, ivalue+ioffset, bufflen);
    }   
}
/**************************************************** */
/* end of BuildBinString()                            */
/**************************************************** */


/* *************************************************** */
/* Raw (floating pt) Function Translation Code         */
/* *************************************************** */
#ifdef ANSI
float  
raw_eval( GENE_DATA *buff,
	  double *Params,
	  int    n_degray,  /* n_degray is non-zero ONLY for f97 */
	  char   val_mode)
#else
float                       /* n_degray is non_zero ONLY for f97 */
raw_eval(buff, Params, n_degray, val_mode)
  GENE_DATA *buff;
  double    *Params;
  int       n_degray;
  char      val_mode;
#endif
{

  GENE_DATA tmp_buf[1000];
  GENE_DATA *nw_buf=NULL;
  int       i;
  float     sum;
  float     (*eval_func)();
  static int  pow2_length=0;                  /* nearest power of 2 of length */

  if (val_mode == 'b' || val_mode == 'B')     /* if we are looking at strings */
    {
      if (Coding[0] == GREY)
	{
	  degray(buff, tmp_buf, Nb_params, Work_params);
	  nw_buf=tmp_buf;
	}
      else
	nw_buf=buff;

      /* transform the binary to input parameters */
      sym_trnsfrm(nw_buf, Params, Nb_params, Work_params);
 
    }

  if (Func)
     f_fx(Func, &eval_func);    /* look up the address of the function */
  
  switch (Func)
  {
     case 0:                       /* handle function de-composition */
	sum = compos_fx(Params, Nb_params, Cfunc1, Cfunc2, Cplxty[0], Cform[0]);
	break;

     case 1:
     case 100:
     case 3:
     case 4:
     case 5:
     case 6:
     case 7:
     case 20:
     case 21:
     case 22:
	sum = (*eval_func)(Params, Nb_params);
        break;
     
     case 2:
     case 8:
     case 9:
     case 10:
	sum = cplx_fx(Params, Nb_params, eval_func, Cplxty[0]);
        break;

     case 97:                    /* handle function n-de-composition */
	if (val_mode=='f' || val_mode=='F')
	  {
	    for (i=0; i< Nb_params; i++)
	      BuildBinString(&(tmp_buf[i*Work_params[0].lngth]),Params[i]);
	    nw_buf=tmp_buf;
	  }

	sum = (*eval_func)(nw_buf, Params, Nb_params, Cfunc2, n_degray);
	break;

     case 98:                    /* handle general function de-composition */
	if (val_mode=='f' || val_mode=='F')
	  {
	    for (i=0; i< Nb_params; i++)
	      BuildBinString(&(tmp_buf[i*Work_params[0].lngth]),Params[i]);
	    nw_buf=tmp_buf;
	  }

        if (pow2_length == 0) 
	 {
	   pow2_length = (int)pow(2.0,(double)
				  ((int)(LOG2((double)Work_params[0].lngth))));
	   
	   if (pow2_length < (int)Work_params[0].lngth) pow2_length*=2;

	 }
        
	sum = (*eval_func)(nw_buf, Params, Nb_params, Cfunc2, pow2_length);
	break;
     
     default: 
	printf("Should never reach default statement of raw_eval\n");
	exit(99);
  }

  return(sum);
}
/*****************************************************/
/* end of raw_eval()                                 */
/*****************************************************/



/* **************************************** */
/* Complex Function Distinction Code        */
/* **************************************** */
#ifdef ANSI
float  
cplx_fx( double *Params,
	 int    Nb_Params,
	 float  (*eval)(double *Params, int Nb_Params),
	 char   cplxty)
#else
float 
cplx_fx(Params, Nb_Params, eval, cplxty)
    double *Params;
    int    Nb_Params;
    float  (*eval)( );
    char   cplxty;

#endif
{
  int    k, m;
  float  sum, rl_sum, interim;
  double t_params[2];

  sum = rl_sum = 0.0;
  switch (cplxty)
  {
     case FULL:
        for (k = 0; k < (Nb_Params-1); k++)
        {
           for (m = k+1; m < Nb_Params; m++)
           {
              t_params[0] = Params[k];
              t_params[1] = Params[m];
              interim = (*eval)(t_params, 2);
	      if (NoiseMask & XINNER)
	      {
		 rl_sum += interim;
		 sum += gaussian() + interim;
	      }
	      else
		 sum += interim;

              /* now do it in reverse to cover other triangle */
              t_params[0] = Params[m];
              t_params[1] = Params[k];
              interim = (*eval)(t_params, 2);
	      if (NoiseMask & XINNER)
	      {
		 rl_sum += interim;
		 sum += gaussian() + interim;
	      }
	      else
		 sum += interim;
           }
        }
	break;
     
     case LOWERLEFT:
        for (k = 0; k < (Nb_Params-1); k++)
        {
           for (m = k+1; m < Nb_Params; m++)
           {
              t_params[0] = Params[k];
              t_params[1] = Params[m];
              interim = (*eval)(t_params, 2);
	      if (NoiseMask & XINNER)
	      {
		 rl_sum += interim;
		 sum += gaussian() + interim;
	      }
	      else
		 sum += interim;
           }
        }
	break;
    
     case STANDARD:
	sum = (*eval)(Params, Nb_Params);
	break;

     case WRAP:
        for (k = 0; k < Nb_Params; k++)
        {
           t_params[0] = Params[k];
           m = (k + 1) % Nb_Params;
           t_params[1] = Params[m];
           interim = (*eval)(t_params, 2);
	   if (NoiseMask & XINNER)
	   {
	      rl_sum += interim;
	      sum += gaussian() + interim;
           }
	   else
	     sum += interim;
        }
	break;

     default:
        fprintf(stdout, "Invalid Function Complexity!\n");
        exit(7);
	break;
  }

  if (NoiseMask & XINNER)
     set_tru_fit(rl_sum);

  if (XOUTTER == (NoiseMask & XBOTH))
     set_tru_fit(sum);

  if (NoiseMask & XOUTTER)
     sum += gaussian();

  return(sum);
}
/*****************************************************/
/* end of cplx_fx()                                  */
/*****************************************************/



/* **************************************** */
/* Function Composition Code                */
/* **************************************** */
#ifdef ANSI
float  
compos_fx( double *Params,
	   int    Nb_Params,
	   int    cf1,
	   int    cf2,
	   char   cplxty,
	   char   cform)
#else
float 
compos_fx(Params, Nb_Params, cf1, cf2, cplxty, cform)
    double *Params;
    int    Nb_Params;
    int    cf1;
    int    cf2;
    char   cplxty;
    char   cform;

#endif
{
  int    k, m;
  int    cnt_tmp;
  float  sum;
  double t_params[2];
  double tmp_stor[100000];
  float  noise_stor[100000];
  float  (*eval1)();
  float  (*eval2)();

  f_fx(cf1, &eval1);
  f_fx(cf2, &eval2);

  cnt_tmp = 0;
  switch (cplxty)
  {
     case FULL:
        for (k = 0; k < (Nb_Params-1); k++)
        {
           for (m = k+1; m < Nb_Params; m++)
           {
              t_params[0] = Params[k];
              t_params[1] = Params[m];
	      tmp_stor[cnt_tmp] = (*eval2)(t_params, 2);
	      if (NoiseMask & XINNER)
		 noise_stor[cnt_tmp] = gaussian();
	      cnt_tmp++;
              /* now do it in reverse to cover other triangle */
              t_params[0] = Params[m];
              t_params[1] = Params[k];
	      tmp_stor[cnt_tmp] = (*eval2)(t_params, 2);
	      if (NoiseMask & XINNER)
		 noise_stor[cnt_tmp] = gaussian();
	      cnt_tmp++;
           }
        }
	break;
     
     case LOWERLEFT:
        for (k = 0; k < (Nb_Params-1); k++)
        {
           for (m = k+1; m < Nb_Params; m++)
           {
              t_params[0] = Params[k];
              t_params[1] = Params[m];
	      tmp_stor[cnt_tmp] = (*eval2)(t_params, 2);
	      if (NoiseMask & XINNER)
		 noise_stor[cnt_tmp] = gaussian();
	      cnt_tmp++;
           }
        }
	break;
    
     case STANDARD:
	tmp_stor[0] = (*eval2)(Params, Nb_Params);
	if (NoiseMask & XINNER)
	   noise_stor[cnt_tmp] = gaussian();
	cnt_tmp++;
	break;

     case WRAP:
        for (k = 0; k < Nb_Params; k++)
        {
           t_params[0] = Params[k];
           m = (k + 1) % Nb_Params;
           t_params[1] = Params[m];
	   tmp_stor[cnt_tmp] = (*eval2)(t_params, 2);
	   if (NoiseMask & XINNER)
	      noise_stor[cnt_tmp] = gaussian();
	   cnt_tmp++;
        }
	break;

     default:
        fprintf(stdout, "Invalid Function Complexity!\n");
        exit(7);
	break;
  }

  sum = 0.0;
  switch (cform)
  {
     case ADDITIV:
	for (k = 1; k < cnt_tmp; k++)
	   tmp_stor[0] += tmp_stor[k];
	   
	sum = (*eval1)(tmp_stor, 1);

	if (NoiseMask & XBOTH)
	{
	   set_tru_fit(sum);

           if (NoiseMask & XINNER)
	   {
	      for (k = 0; k < cnt_tmp; k++)
	         tmp_stor[0] += (double)noise_stor[k];
	      
	      sum = (*eval1)(tmp_stor, 1);
	   }
	
	   if (NoiseMask & XOUTTER)
	      sum += gaussian();
        }
	break;
     
     case DISTRIBUTED:
	sum = (*eval1)(tmp_stor, cnt_tmp);
	
	if (NoiseMask & XBOTH)
	{
	   set_tru_fit(sum);
           
	   if (NoiseMask & XINNER)
	   {
	      for (k = 0; k < cnt_tmp; k++)
	         tmp_stor[k] += (double)noise_stor[k];
	      
	      sum = (*eval1)(tmp_stor, cnt_tmp);
	   }
	   
	   if (NoiseMask & XOUTTER)
	      sum += gaussian();
	}
	break;
    
     case MULTIPLICTIV:
	for (k = 0; k < cnt_tmp; k++)
	   sum += (*eval1)(&tmp_stor[k], 1);
	
	if (NoiseMask & XBOTH)      /* really asks if either in or out is set */
	{
	   set_tru_fit(sum);        /* valid as long as terminal is not added */
	   
	   if (NoiseMask & XINNER)
	   {
	      sum = 0.0;
	      for (k = 0; k < cnt_tmp; k++)
	      {
	         tmp_stor[k] += (double)noise_stor[k];
	         sum += (*eval1)(&tmp_stor[k], 1);

	         if (NoiseMask & XOUTTER)
	            sum += gaussian();
	      }
	   }
	   
	   if (XOUTTER == (NoiseMask & XBOTH))      /* only outter noise term */
           {
	      sum = 0.0;
	      for (k = 0; k < cnt_tmp; k++)
	         sum += (*eval1)(&tmp_stor[k], 1) + gaussian();
           }
        }	

	if (NoiseMask & XTERMINAL)
        {
	   set_tru_fit(sum);
	   sum += gaussian();
        }

	break;
     
     case SIMPLE:
	sum = (*eval1)(tmp_stor, 1);
        break;

     default:
        fprintf(stdout, "Invalid Composition!\n");
        exit(7);
	break;
  }

  return(sum);
}
/*****************************************************/
/* end of compos_fx()                                */
/*****************************************************/

/* **************************************** */
/* Function Address Lookup Code             */
/* **************************************** */
#ifdef ANSI
void  
f_fx( int    cf,
      float  (**eval)( ))
#else
void 
f_fx(cf, eval)
    int    cf;
    float  (**eval)( );
#endif
{

  switch (cf)
  {
     case 1:
	*eval = f1;
	break;
     
     case 2:
	*eval = f2;
	break;
     
     case 3:
	*eval = f3;
	break;
     
     case 4:
	*eval = f4;
	break;
     
     case 5:
	*eval = f5;
	break;
     
     case 6:
	*eval = f6;
	break;
     
     case 7:
	*eval = f7;
	break;
     
     case 8:
	*eval = f8;
	break;
     
     case 9:
	*eval = f9;
	break;
     
     case 10:
	*eval = f10;
	break;
     
     case 20:
	*eval = f20;
	break;
     
     case 21:
	*eval = f21;
	break;
     
     case 22:
	*eval = f22;
	break;
     
     case 97:
	*eval = f97;
	break;

     case 98:
	*eval = f98;
	break;

     case 99:
	*eval = f99;
	break;
     case 100:
	*eval = f100;
     	break;

     default:
	fprintf(stdout, "Function chosen is not supported!\n");
	exit(7);
	break;
  }

}
/*****************************************************/
/* end of f_fx()                                     */
/*****************************************************/



/* **************************************** */
/* Master Function Parameter Report Code    */
/* **************************************** */
#ifdef ANSI
void  
mstr_rpt( GENE_DATA *buff)
#else
void 
mstr_rpt(buff)
  GENE_DATA *buff;
#endif
{

  GENE_DATA  *nw_buf;
  GENE_DATA  gbuf[5000], dbuf[5000], ddbuf[5000];
  int        buf_toggle=0;
  GENE_DATA  buf[2][5000];

  int        i, j, k, prms, psn;
  int        terms;

  if (Coding[0] == GREY)                    /* if this a grey coded string */
  {
      degray(buff, Grey_gene, Nb_params, Work_params);
      nw_buf = Grey_gene;
  }
  else
     nw_buf = buff;


  if (Func < 97)
  {
     sym_trnsfrm(nw_buf, Params, Nb_params, Work_params);
     for (k=0; k<Nb_params; k++)
     {
        fprintf(stdout, "P%d: %f    ", k, Params[k]);
        if ((k+1)%5 == 0)
           fprintf(stdout, "\n");
     }
  }
  else
  {
     switch (Func)
      {
       case 98:
	 k=0;
	 terms = (int)pow(2.0,(double)
			  ((int)(LOG2((double)Work_params[0].lngth))));
	 if (terms < (int)Work_params[0].lngth) terms*=2;
	 
	 prms = Nb_params / terms;
	 psn = 0;
	 
	 gen_trnsfrm(nw_buf, Params, prms, &Work_params[psn]);
	 for (j=0; j<prms; j++)
	   fprintf(stdout,"P%d, %f\t",k++, Params[j]);
	 
	 for (i=1; i<(terms); i++)            /* generate inner terms */
	  {
	     psn += prms;
	     /* Generate a degray term */
	     degray(&nw_buf[Work_params[psn].posn], buf[buf_toggle], prms, 
		    &Work_params[psn]);
	     for (j=1; j<i; j++)
	      {
		 degray(buf[buf_toggle], buf[!buf_toggle], prms, 
			&Work_params[psn]);
		 buf_toggle=!buf_toggle;
	      }
	     
	     gen_trnsfrm(buf[buf_toggle], Params, prms, &Work_params[psn]);
	     for (j=0; j<prms; j++)
	      {
		 fprintf(stdout,"P%d, %f\t",k++, Params[j]);
		 
		 if ((k+1)%5 == 0)
		   fprintf(stdout,"\n");
	      }
	  }
	 break;
       case 99:
	 k = 0;
	 prms = Nb_params / 4.0;
	 sym_trnsfrm(nw_buf, Params, prms, Work_params);
	 fprintf(stdout, "P%d: %f    P%d: %f", k, Params[0], k, Params[1]);
	 
	 k += 2;
	 psn = prms;
	 gray(&buff[Work_params[psn].posn], gbuf, prms, &Work_params[psn]);
	 sym_trnsfrm(gbuf, Params, prms, &Work_params[psn]);
	 fprintf(stdout, "P%d: %f    P%d: %f", k, Params[0], k, Params[1]);
	 psn += prms;
	 
	 k += 2;
	 degray(&buff[Work_params[psn].posn], dbuf, prms, &Work_params[psn]);
	 sym_trnsfrm(dbuf, Params, prms, &Work_params[psn]);
	 fprintf(stdout, "P%d: %f    P%d: %f", k, Params[0], k, Params[1]);
	 psn += prms;
	 
	 
	 k += 2;
	 degray(&buff[Work_params[psn].posn], dbuf, prms, &Work_params[psn]);
	 degray(dbuf, ddbuf, prms, &Work_params[psn]);
	 sym_trnsfrm(ddbuf, Params, prms, &Work_params[psn]);
	 fprintf(stdout, "P%d: %f    P%d: %f", k, Params[0], k, Params[1]);
         break;
      }
  }
  
  fprintf(stdout, "\n");
  fflush(stdout);
}
/*****************************************************/
/* end of mstr_rpt()                                 */
/*****************************************************/


/* **************************************** */
/* Raw Parameter Report Code                */
/* **************************************** */
#ifdef ANSI
void  
raw_rpt(double *params)
#else
void 
raw_rpt(params)
  double *params;
#endif
{

  int k;

  for (k=0; k<Nb_params; k++)
  {
     fprintf(stdout, "P%d: %f    ", k, Params[k]);
     if ((k+1)%5 == 0)
        fprintf(stdout, "\n");
  }
  
  fprintf(stdout, "\n\n");
  fflush(stdout);
}
/*****************************************************/
/* end of raw_rpt()                                  */
/*****************************************************/

/* *********************************** */
/* Gaussian Noise Function Code        */
/* *********************************** */
#ifdef ANSI
float  
gaussian( )
#else
float 
gaussian( )
#endif
{
  int k;
  float gaus;


  gaus = 0.0;
  for (k=0; k<12; k++)
      gaus += fracrand();

  gaus -= 6.0;

  return(gaus);
}
/*********************************/
/* end of gaussian()             */
/*********************************/



/* ******************************* */
/* DeJong Function F1 Code         */
/* ******************************* */
#ifdef ANSI
float  
f1( double *Params,
    int     Nb_params)
#else
float 
f1(Params, Nb_params)
  double *Params;
  int     Nb_params;
#endif
{
  int k;
  double sum;
  float rtn;

  sum = 0;
  for (k=0; k < Nb_params; k++)
     sum += (Params[k] * Params[k]);
 
  rtn = (float) sum;
  return(rtn);
}
/*********************************/
/* end of f1()                   */
/*********************************/


/* ********************************** */
/* DeJong Function F2 Code            */
/* ********************************** */
#ifdef ANSI
float  
f2( double *Params,
    int    Nb_params)
#else
float 
f2(Params, Nb_params)
  double *Params;
  int    Nb_params;
#endif
{
  float rtn;
  double sum;
  double x1, x2, sq_x1, diff_x1;

  x1 = Params[0];        x2 = Params[1];
  sq_x1 = x1 * x1;       
  diff_x1 = 1.0 - x1;
  sum = (100.0 * ((sq_x1 - x2) * (sq_x1 - x2)) ) + 
	(diff_x1 * diff_x1);

  rtn = (float)(sum);
  return(rtn);
}
/***********************************/
/* end of f2()                     */
/***********************************/



/* *********************************** */
/* DeJong Function F3 Code             */
/* *********************************** */
#ifdef ANSI
float  
f3( double *Params,
    int    Nb_params)
#else
float 
f3(Params, Nb_params)
  double *Params;
  int    Nb_params;
#endif
{
  int k;
  int sum, part;
  float rtn;

  sum = 0;
  for (k=0; k<Nb_params; k++)
  {
      part = (int) (Params[k]);
      if (part > Params[k])
	 part -= 1;
      sum += part;
  }

  rtn = (float)(sum);
  return(rtn);
}
/**************************************/
/* end of f3()                        */
/**************************************/


/* *********************************** */
/* DeJong Function F4 Code             */
/* *********************************** */
#ifdef ANSI
float  
f4( double *Params,
    int    Nb_params)
#else
float 
f4(Params, Nb_params)
  double *Params;
  int    Nb_params;
#endif
{
  int k;
  float rtn;
  double sum, quad;

  sum = 0.0;
  for (k=0; k<Nb_params; k++)
  {
     quad = Params[k] * Params[k] * Params[k] * Params[k];
     quad *= (k+1);
     sum += quad;
  }

  rtn = (float)(sum);
  set_tru_fit(rtn);
  rtn += gaussian();

  return(rtn);
}
/*********************************/
/* end of f4()                   */
/*********************************/

/* *********************************** */
/* DeJong Function F5 Code             */
/* *********************************** */
#ifdef ANSI
float  
f5( double *Params,
    int     Nb_params)
#else
float 
f5(Params, Nb_params)
  double *Params;
  int     Nb_params;
#endif
{
  int k;
  double sum, subsum, diff1, diff2, invrs;
  float rtn;

  sum = 0.0;
  for (k=0; k < 25; k++)
  {
     diff1 = Params[0] - f5_arr[0][k];
     diff2 = Params[1] - f5_arr[1][k];
     subsum = (diff1 * diff1 * diff1 * diff1 * diff1 * diff1) +
              (diff2 * diff2 * diff2 * diff2 * diff2 * diff2);
     subsum += (double) (k + 1.0);
     subsum = (double) (1.0) / subsum;
     sum += subsum;
  } 

  sum += (double) (0.002);
  invrs = (double) (1.0) / sum;
  rtn = (float) invrs;
  return(rtn);
}
/*********************************/
/* end of f5()                   */
/*********************************/


/* *********************************** */
/* Rastrigin Function F6 Code          */
/* *********************************** */
#ifdef ANSI
float  
f6( double *Params,
    int     Nb_params)
#else
float 
f6(Params, Nb_params)
  double *Params;
  int     Nb_params;
#endif
{
  int k;
  double sum;
  float rtn;
  double two;
  double A = 10.0;

  two = (double)(2.0);
  sum = 0.0;
  for (k=0; k < Nb_params; k++)
     sum += ((Params[k] * Params[k]) - (A * (cos(two * M_PI * Params[k]))));

  sum += (double) (Nb_params * A);
  rtn = (float) sum;
  return(rtn);
}
/*********************************/
/* end of f6()                   */
/*********************************/




/* *********************************** */
/* Schewfel Function F7 Code           */
/* *********************************** */
#ifdef ANSI
float  
f7( double *Params,
    int     Nb_params)
#else
float 
f7(Params, Nb_params)
  double *Params;
  int     Nb_params;
#endif
{
  int k;
  double sum;
  float rtn;

  sum = 0.0;
  for (k=0; k < Nb_params; k++)
     sum += (-Params[k] * (sin( sqrt( fabs(Params[k])))));

  rtn =(float) sum;
  return(rtn);
}
/*********************************/
/* end of f7()                   */
/*********************************/



/* *********************************** */
/* Griewangk Function F8 Code          */
/* *********************************** */
#ifdef ANSI
float  
f8( double *Params,
    int     Nb_params)
#else
float 
f8(Params, Nb_params)
  double *Params;
  int     Nb_params;
#endif
{
  int k;
  double sum, prod;
  float rtn;

  sum = 0.0;   prod = 1.0;
  for (k=0; k < Nb_params; k++)
  {
     sum += ((Params[k] * Params[k]) / 4000.0);
     prod *= (cos(Params[k]/sqrt(k+1)));
  }

  sum -= prod;
  sum += (double) 1.0;
  rtn = (float) sum;
  return(rtn);
}
/*********************************/
/* end of f8()                   */
/*********************************/





/* *********************************** */
/* Eshelmans Function F9 Code          */
/* This is the Sine envelope sine wave */
/* *********************************** */
#ifdef ANSI
float  
f9( double *Params,
    int    Nb_params)
#else
float 
f9(Params, Nb_params)
  double *Params;
  int    Nb_params;
#endif
{
  float rtn;
  double sq1, sq2;
  double sum, sq_sum, temp, temp2;
		
    sq1 = Params[0] * Params[0];
    sq2 = Params[1] * Params[1];
    sq_sum = sq1 + sq2;

    temp = 1.0 + (0.001 * sq_sum);
    temp2 = sin(sqrt(sq_sum));
    sum = 0.5 + (temp2*temp2 - 0.5) /(temp*temp);  

  rtn = (float)(sum);
  return(rtn);
}
/*********************************/
/* end of f9()                   */
/*********************************/



/* ********************************************* */
/* Eshelmans Function F10 Evaluation Code        */
/* This is the Streched V sine wave              */
/* ********************************************* */
#ifdef ANSI
float  
f10( double *Params,
     int    Nb_params)
#else
float 
f10(Params, Nb_params)
  double *Params;
  int    Nb_params;
#endif
{
  float rtn;
  double sq1, sq2;
  double sum, sq_sum, temp1, temp2;


    sq1 = Params[0] * Params[0];
    sq2 = Params[1] * Params[1];
    sq_sum = sq1 + sq2;
	
    temp1 = sqrt(sq_sum);
    temp2 = sin(pow(temp1,0.2) * 50.0);
    sum = sqrt(temp1) * (temp2*temp2 + 1.0); 


  rtn = (float)(sum);
  return(rtn);
}
/*********************************/
/* end of f10()                  */
/*********************************/


/* *********************************** */
/* Corana Function F20 Code            */
/* *********************************** */
#ifdef ANSI
float  
f20( double *Params,
    int     Nb_params)
#else
float 
f20(Params, Nb_params)
  double *Params;
  int     Nb_params;
#endif
{
  int k;
  double sum, dev, z, x, run;
  float rtn;

  sum = 0.0;
  for (k=0; k < Nb_params; k++)
  {
     x = Params[k];
     z = floor((fabs(x / 0.20)) + 0.49999999);
     if (x < 0.0)
	z = 0.0 - z;
     z *= 0.20;

     dev = fabs(x - z);
     if (dev < 0.05)
     {
	if (z < 0.0)
	   run = (0.05 + z) * (0.05 + z) * 0.15 * d_arr[k];
	else 
	   run = (-0.05 + z) * (-0.05 + z) * 0.15 * d_arr[k]; 
     }
     else
	run =  d_arr[k] * x * x;

     sum += run;
  }

  rtn = (float) sum;
  return(rtn);
}
/*********************************/
/* end of f20()                  */
/*********************************/


/* *********************************** */
/* Colville Function F21 Code          */
/* *********************************** */
#ifdef ANSI
float  
f21( double *Params,
    int     Nb_params)
#else
float 
f21(Params, Nb_params)
  double *Params;
  int     Nb_params;
#endif
{
  double sum, x1, x2, x3, x4;
  float rtn;

  sum = 0.0;  
  x1 = Params[0];
  x2 = Params[1];
  x3 = Params[2];
  x4 = Params[3];

   sum = 100 * (x2 - (x1 * x1)) * (x2 - (x1 * x1)) + (1 - x1) * (1 - x1) +
          90 * (x4 - (x3 * x3)) * (x4 - (x3 * x3)) + (1 - x3) * (1 - x3) +
        10.1 * ((x2 - 1) * (x2 - 1)  + (x4 - 1) * (x4 - 1)) +
        19.8 * (x2 - 1) * (x4 - 1);

  rtn = (float) sum;
  return(rtn);
}
/*****************************************************/
/* end of f21()                                      */
/*****************************************************/
	
   
/* ************************************************ */
/* Function to Average Nb_params number of numbers  */
/* ************************************************ */
#ifdef ANSI
float  
f22( double *Params,
     int    Nb_params)
#else
float 
f22(Params, Nb_params)
  double *Params;
  int    Nb_params;
#endif
{
  int    k;
  float  rtn;
  double sum;

  sum = 0.0;
  for (k = 0; k < Nb_params; k++)
     sum += Params[k];

  sum /= Nb_params;

  rtn = (float)(sum);
  return(rtn);
}
/*********************************/
/* end of f22()                  */
/*********************************/

/* ********************************************************* */
/* CSU n-degray function to be used with raw_eval            */
/*     This function isolates a single term from f98 by      */
/*     performing degray's n times on a parameter set.       */
/*  NOTE: THIS IS NOT A TEST SUITE FUNCTION!  It is just     */
/*        to use to get an idea of what the transformations  */
/*        do to one of the independent f98 terms.            */
/* ********************************************************* */
#ifdef ANSI
float  
f97( GENE_DATA  *buff,
     double     *Params,
     int        Nb_params,
     int        cf,
     int        num_dgry)
#else
float 
f97(buff, Params, Nb_params, cf, num_dgry)
  GENE_DATA  *buff;
  double     *Params;
  int        Nb_params;
  int        cf;
  int        num_dgry;
#endif
{

  float      (*eval)();
  int       buf_toggle=0;
  GENE_DATA  buf[2][800];

  int        i;

  f_fx(cf, &(eval));

  if (num_dgry>0)
   {
      /* Generate a degray term */
      degray(buff, buf[buf_toggle], Nb_params, Work_params);
      for (i=1; i<num_dgry; i++)
       {
	  degray(buf[buf_toggle], buf[!buf_toggle], Nb_params, Work_params);
	  buf_toggle=!buf_toggle;
       }
      gen_trnsfrm(buf[buf_toggle], Params, Nb_params, Work_params);
   }
  else
    gen_trnsfrm(buff, Params, Nb_params, Work_params);

  return ((*eval)(Params, Nb_params));       

}
/*****************************************************/
/* end of f97()                                      */
/*****************************************************/

/* ************************************* */
/* CSU Gray-BCD Equivalent Function Code */
/* ************************************* */
#ifdef ANSI
float  
f98( GENE_DATA  *buff,
     double     *Params,
     int        Nb_params,
     int        cf,
     int        terms)
#else
float 
f98(buff, Params, Nb_params, cf, terms)
  GENE_DATA  *buff;
  double     *Params;
  int        Nb_params;
  int        cf;
  int        terms;
#endif
{

  int        prms, psn;
  float      sum;
  float      (*eval)();
  int       buf_toggle=0;
  GENE_DATA  buf[2][800];

  int        i,j;

  f_fx(cf, &(eval));
  prms = Nb_params / terms;
  psn = 0;

  gen_trnsfrm(buff, Params, prms, &Work_params[psn]);
  sum = (*eval)(Params, prms);          /* first terms */

  for (i=1; i<(terms); i++)            /* generate inner terms */
   {
     psn += prms;
      /* Generate a degray term */
      degray(&buff[Work_params[psn].posn], buf[buf_toggle], prms, 
	     &Work_params[psn]);
      for (j=1; j<i; j++)
       {
	degray(buf[buf_toggle], buf[!buf_toggle], prms, &Work_params[psn]);
	buf_toggle=!buf_toggle;
       }
      gen_trnsfrm(buf[buf_toggle], Params, prms, &Work_params[psn]);
      
      sum += (*eval)(Params, prms);       

   }

  return(sum);
}
/*****************************************************/
/* end of f98()                                      */
/*****************************************************/

/* ********************************************** */
/* CSU Gray-BCD Close Approximation Function Code */
/* ********************************************** */
#ifdef ANSI
float  
f99( GENE_DATA  *buff,
     double     *Params,
     int        Nb_params,
     int        cf)
#else
float 
f99(buff, Params, Nb_params, cf)
  GENE_DATA  *buff;
  double     *Params;
  int        Nb_params;
  int        cf;
#endif
{

  int        prms, psn;
  float      sum;
  float      (*eval)();
  GENE_DATA  gbuf[800], dbuf[800], ddbuf[800];


  f_fx(cf, &(eval));
  prms = Nb_params / 4.0;
  psn = prms;

  sum = (*eval)(Params, prms);          /* first set of function terms */

  
  gray(&buff[Work_params[psn].posn], gbuf, prms, &Work_params[psn]);
  sym_trnsfrm(gbuf, Params, prms, &Work_params[psn]);
  sum += (*eval)(Params, prms);          /* second set of function terms */
  psn += prms;

  degray(&buff[Work_params[psn].posn], dbuf, prms, &Work_params[psn]);
  sym_trnsfrm(dbuf, Params, prms, &Work_params[psn]);
  sum += (*eval)(Params, prms);          /* third set of function terms */
  psn += prms;


  degray(&buff[Work_params[psn].posn], dbuf, prms, &Work_params[psn]);
  degray(dbuf, ddbuf, prms, &Work_params[psn]);
  sym_trnsfrm(ddbuf, Params, prms, &Work_params[psn]);
  sum += (*eval)(Params, prms);          /* fourth set of function terms */

  return(sum);
}
/*****************************************************/
/* end of f99()                                      */
/*****************************************************/

/* ******************************* */
/*          SHG Code               */
/* ******************************* */
#ifdef ANSI
float
f100( double *Params,
    int     Nb_params)
#else
float
f100(Params, Nb_params)
  double *Params;
  int     Nb_params;
#endif
{
  float rtn;
  int status, igen=1,icalc=1,params=Nb_params;
  double shg,dummy;
/*double  xvec[params];
for (unsigned i = 0; i<params ; i++)
        xvec[i] = Params[i];
*/
shgval_(&myproc,&icalc,&params,&shg,Params);
//result = -1.0*shg;
//cout << "C alignshg done" <<endl;
rtn = (float) (9.070365-shg);
return(rtn);
}
/*********************************/
/* end of f100()                   */
/*********************************/
	
	

