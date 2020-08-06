#include <stdio.h>
#include <math.h>



#include "gene.h"
#include "ga_random.h"
#include "op_hux.h"

/***************************************************************************
 * FUNCTION: allocate_hux_params
 *
 * DESCRIBE:  Allocates space for a difference table
 *            (defined in op_hux_params.h)
 *
 * INPUT PARAMETERS: 
 *
 * RETURN VALUE: 
 ****************************************************************************/
#ifdef ANSI
int *
allocate_hux_params (int   string_length)

#else
int *
allocate_hux_params (string_length)
                     int   string_length;
#endif
{

  int  *arry;

  if (!(arry = (int *) calloc(string_length,sizeof(int))))
     fatal_error("allocate hux params: Malloc failure for difference array.\n");

  return (arry);
}

/*****************************************************/
/* nint() isn't on all architectures, so this is     */
/* our version, our_nint.  It rounds up at .5        */
/*****************************************************/
#ifdef ANSI
int our_nint(double q)
#else
int our_nint(q)
double q;
#endif
{
double z;
z = floor(q);
if (q - z >= .5) 
      return(z+1);
else 
      return (z);
}

/*****************************************************/
/* This routine performs the HUX crossover as defined
   by L. Eschelman.                                  */
/*****************************************************/
#ifdef ANSI
void 
hux( char    *mom,
     char    *dad,
     int     *dif_ary,
     int     length)
#else
void 
hux(mom, dad, dif_ary, length)
    char    *mom;
    char    *dad;
    int     *dif_ary;
    int     length;
#endif
{
   int  k;
   int  loop, postn, indx;
   char temp;
   int  diff_cnt = 0;

   for (k=0; k<length; k++)
   {
      if (mom[k] != dad[k])
      {
	 dif_ary[diff_cnt] = k;
	 diff_cnt++;
      }
   }
   loop = our_nint((double) diff_cnt/2);
   for (k=0; k<loop; k++)
   {
      postn = (int)(fracrand() * diff_cnt);
      while(dif_ary[postn] < 0)
	 postn = (int)(fracrand() * diff_cnt);

      indx = dif_ary[postn];
      temp = mom[indx];
      mom[indx] = dad[indx];
      dad[indx] = temp;
      dif_ary[postn] = -1;
   }

}
