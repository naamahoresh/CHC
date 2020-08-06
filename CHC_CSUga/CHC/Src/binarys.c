/**************    Binary Generic Functions ***************/


#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include "gene.h"
#include "binarys.h"


#ifdef ANSI
void
allocate_bin_params (double **prb_params,
                     struct p_list **Work_params,
                     int Nb_params,
		     int Length,
		     char **gray_wrk)

#else
void allocate_bin_params (prb_params, Work_params, Nb_params, Length, gray_wrk)
double **prb_params;
struct p_list **Work_params;
int Nb_params;
int Length;
char **gray_wrk;
#endif
{

  *prb_params = (double *) calloc(Nb_params, sizeof(double));
  *gray_wrk = (char *) calloc(Length, sizeof(char));
  *Work_params = (struct p_list *) malloc (Nb_params * sizeof(struct p_list));
}
/*****************************************************/
/* end of allocate_bin_params()                       */
/*****************************************************/


/* *********************************************************** */
/*  ctod (Char TO Double)                                       */
/*    takes in a string of binary characters.  It then         */
/*    then picks off the next weight of "length" number        */
/*    of bits and coverts it to an integer using Horner's      */
/*    Scheme. That integer is then converted to some real      */
/*    discrete value using the step value.  For example,       */
/*    if the step value is 0.25, then the possible values      */
/*    generated are 0.0, 0.25, 0.50, 0.75, etc.                */
/*    Also, note that we decode negative integers by           */
/*    treating the binary characters as sign magnitude         */
/*    strings.                                                 */
/*                                                             */
/*    A typical call to ctod                                   */
/*                                                             */
/*   g = ctod(&w[ndx*bits_per_param], bits_per_param, scale);  */
/*                                                             */
/* Author:  Darrell Whitley                                    */
/************************************************************* */
/*                                                             */
/*  Changed 7/4/90                                             */
/*  Now accepts a parameter to determine if it is a SIGNED     */
/*  number that is being converted or if the number is not     */
/*  SIGNED.                                                    */
/*                                                             */
/*    A typical call to ctod                                   */
/*                                                             */
/*   g = ctod(&w[ndx*bits_per_param], bits_per_param, scale, sign);  */
/*                                                             */
/************************************************************* */
#ifdef ANSI
double
ctod ( char   *string,
       int    length,
       double step_value,
       int    sign)
#else
double 
ctod (string, length, step_value, sign)
    char *string;
    int length;
    double step_value;
    int sign;
#endif
{
    int i;
    double sum;
    int sign_val;


    sign_val = 0;
    sum = 0;
    if ( ((length < 2) && sign) || (length < 1) )
        return(sum);

    if (sign)
    {
       if (*string == '1')
          sign_val = 1;
       string++;
       length--;
    }

    for (i=1; i <= length; i++)
    {
        sum += sum;
        if (*string == '1')
            sum++;
        string++;
    }

    if (sign_val)
       sum = (-1) * sum;
    
    /* to facilitate new mapping where 100000 and 000000 are not the same. */
    if (sign_val)
       sum = sum - 1;
    /*********************************************************************/

    sum = sum * step_value;

    return (sum);
}
/*****************************************************/
/* end of ctod()                                     */
/*****************************************************/


/*****************************************************************/
/*              Transform is set up to decode the bit            */
/*              string into the correct real valued parameters   */
/*              to be used in the transformations.               */
/*****************************************************************/


#ifdef ANSI
void 
transform( GENE_DATA       *buff,
              double          *Vals,
              int             Nb_params,
              struct p_list   *Work_params)
#else
void 
transform(buff, Vals, Nb_params, Work_params)
  GENE_DATA     *buff;
  double        *Vals;
  int           Nb_params;
  struct p_list *Work_params;
#endif
{
 
  int i;
נ
  /* decode the bit patterns for the parameters */
  for (i = 0; i < Nb_params; i++)
     Vals[i] = ctod(&buff[Work_params[i].posn], Work_params[i].lngth,
                    Work_params[i].scale, Work_params[i].sign);


}
/******************************************/
/* end of transform()                  */
/******************************************/


/* ******************  Low Order Bit by parameter ***************** */

#ifdef ANSI
int
ham_dist( GENE_DATA *bf1,
          GENE_DATA *bf2,
          int Nb_params,
          struct p_list *Work_params)
#else
int ham_dist(bf1, bf2, Nb_params, Work_params)
  GENE_DATA *bf1, *bf2;
  int Nb_params;
  struct p_list *Work_params;
#endif
{

  int i;
  double diffnt;
  double delta1;
  double delta2;


  diffnt = 0.0;
  for (i = 0; i < Nb_params; i++)
  {
     delta1 = ctod(&bf1[Work_params[i].posn], Work_params[i].lngth,
                   NO_CHANGE, UNSIGNED);
     delta2 = ctod(&bf2[Work_params[i].posn], Work_params[i].lngth,
                   NO_CHANGE, UNSIGNED);
     diffnt += fabs(delta1 - delta2);
     if (diffnt > HAM_LIMIT)
        return(APART);
  }

  return(SAME);
}
/*****************************************************/
/* end of ham_dist()                                 */
/*****************************************************/



/* ******************  Hamming distance up to a limit ***************** */

#ifdef ANSI
int
hamming_diff( GENE_DATA *bf1,
              GENE_DATA *bf2,
              int       length)
#else
int hamming_diff(bf1, bf2, length)
  GENE_DATA *bf1, *bf2;
  int       length;
#endif
{

  int k, diff;

  diff = 0;
  for (k = 0; k < length; k++)
  {
     if (bf1[k] != bf2[k])
	diff++;
     if (diff > HAM_LIMIT)
        return(APART);
  }

  return(SAME);
}
/*****************************************************/
/* end of hamming_diff()                             */
/*****************************************************/


/* ******************  Hamming distance ***************** */

#ifdef ANSI
int
hamming_distance(GENE_DATA *bf1,
                 GENE_DATA *bf2,
                 int       length)
#else
int 
hamming_distance(bf1, bf2, length)
     GENE_DATA *bf1, *bf2;
     int       length;
#endif
{

  int k, diff;

  diff = 0;
  for (k = 0; k < length; k++)
  {
     if (bf1[k] != bf2[k])
	diff++;
  }

  return(diff);
}
/*****************************************************/
/* end of hamming_distance()                             */
/*****************************************************/


#ifdef ANSI
void 
gray (char *inbuf,
      char *outbuf,
      int  Nb_params,
      struct p_list *Work_params)

#else
void gray (inbuf, outbuf, Nb_params, Work_params)
       char *inbuf;
       char *outbuf;
       int  Nb_params;
       struct p_list *Work_params;
#endif
{
  int i, k, m;
  int lngth;
  char last;
  
  k = 0;
  for (m=0; m<Nb_params; m++)
  {
     lngth = Work_params[m].lngth;
     last = '0';
     for (i=0; i<lngth; i++, k++)
     {
        if (inbuf[k] != last)
	   outbuf[k] = '1';
        else
	   outbuf[k] = '0';

        last = inbuf[k];
     }
  }
}
/*****************************************************/
/* end of gray( )                                    */
/*****************************************************/

#ifdef ANSI
void 
degray (char *inbuf,
        char *outbuf,
        int  Nb_params,
	struct p_list *Work_params)

#else
void degray (inbuf, outbuf, Nb_params, Work_params)
       char *inbuf;
       char *outbuf;
       int  Nb_params;
       struct p_list *Work_params;
#endif
{
  int i, k, m;
  int lngth;
  char last;

   k = 0;
   for (m=0; m<Nb_params; m++)
   {
      lngth = Work_params[m].lngth;
      last = '0';
      for (i=0; i<lngth; i++, k++)
      {
         if (inbuf[k] == '1')  
         {
 	    if (last == '0')
	       last = '1';
            else
	       last = '0';
         }
         outbuf[k] = last;
         last = outbuf[k];
      }
   }
}
/*****************************************************/
/* end of degray( )                                  */
/*****************************************************/

