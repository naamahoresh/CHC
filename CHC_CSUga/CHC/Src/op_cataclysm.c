#include <stdio.h>

#include "gene.h"
#include "ga_random.h"
#include "op_cataclysm.h"


/****************************************************************************
 * FUNCTION:  cataclysm
 *
 * DESCRIBE:  This routine performs cataclysmic mutation by using 
 *            some string as a pattern and muatating some percentage 
 *            of the bits in the pattern randomly to form new strings
 *            with which to seed the population.
 *
 *            This routine expects the seed string to be in position
 *            seed_idx and then initializes the population from that
 *            index up to the (size + index) limit.  It also expects 
 *            pct_mut to be a value to be between 0 and 0.99.
 *            A partial seeding can be done by passing a different size 
 *            or length for the population or strings, depending on 
 *            the desired affect.
 * 
 *
 * INPUT PARAMETERS:  pointer to the gene pool
 *                    size - number of consecutive strings to be seeded
 *                    length of the string
 *                    index of string to be used as a pattern
 *                    percent of mutation to apply to pattern
 *                    
 * RETURN VALUE:  None.  
 *                However, the entire pool is reseeded for continued search.
 *                
 ****************************************************************************/
#ifdef ANSI
void 
cataclysm( POOLPTR    pool,
	   int        size,
	   int        length,
	   int        seed_idx,
	   float      pct_mut)
#else
void 
cataclysm(pool, size, length, seed_idx, pct_mut)
          POOLPTR    pool;
	  int        size;
	  int        length;
	  int        seed_idx;
	  float      pct_mut;
#endif
{

   int k, m;
   int mutate, to_seed;

   to_seed = seed_idx + 1;
   for (m=0; m<size; m++)
   {
      for (k=0; k<length; k++)
      {
	 mutate = (int) (fracrand() + pct_mut);
	 if (mutate)
	 {
	    if (pool->data[seed_idx].string[k] == '1')
	       pool->data[to_seed].string[k] = '0';
            else
	       pool->data[to_seed].string[k] = '1';
	 }
	 else
	    pool->data[to_seed].string[k] = pool->data[seed_idx].string[k];
      }
      to_seed++;
   }

}


/****************************************************************************
 * FUNCTION:  parm_cataclysm
 *
 * DESCRIBE:  This routine performs cataclysmic mutation by using 
 *            some string as a pattern and muatating some percentage 
 *            of the paramters in the pattern string, which are randomly 
 *            chosen, to form new strings with which to seed the population.
 *            The parameters are mutated by zeroing out all bits in the 
 *            parameter and then probabalistically turning a small number 
 *            of them on.
 *
 *            This routine expects the seed string to be in position
 *            seed_idx and then initializes the population from that
 *            index up to the (index + size) limit.  It also expects 
 *            pct_mut to be a value to be between 0 and 0.99.
 *            A partial seeding can be done by passing a different size 
 *            or length for the population or strings, depending on 
 *            the desired affect.
 * 
 *
 * INPUT PARAMETERS:  pointer to the gene pool
 *                    size - number of consecutive strings to be seeded
 *                    number of parameters in the string
 *                    number of bits in any one parameter 
 *                    length of the string
 *                    index of string to be used as a pattern
 *                    percent of mutation to apply to pattern
 *                    
 * RETURN VALUE:  None.  
 *                However, the entire pool is reseeded for continued search.
 *                
 ****************************************************************************/
#ifdef ANSI
void 
parm_cataclysm( POOLPTR    pool,
	        int        size,
                int        nb_parms,
                int        sparse,
                int        parm_lngth,
	        int        seed_idx,
	        float      pct_mut)
#else
void 
parm_cataclysm(pool, size, nb_parms, sparse, parm_lngth, seed_idx, pct_mut)
          POOLPTR    pool;
	  int        size;
          int        nb_parms;
          int        sparse;
          int        parm_lngth;
	  int        seed_idx;
	  float      pct_mut;
#endif
{

   int    j, k, m;
   int    mutate, to_seed, idx;
   float  seed_prob;


   seed_prob = (float)sparse / (float)parm_lngth;
   to_seed = seed_idx + 1;
   for (m=0; m<size; m++)
   {
      idx = 0;
      for (k=0; k<nb_parms; k++)
      {
	 mutate = (int) (fracrand() + pct_mut);
	 if (mutate)                        /* if parameter is to be mutated */
	 {
            for (j = idx ; j < (parm_lngth + idx); j++)        
	    {
	       if (fracrand() <=  seed_prob)
	          pool->data[to_seed].string[j] = '1';
	       else
	          pool->data[to_seed].string[j] = '0';
	    }
	 }
	 else
	 {
	    for (j = idx; j < (parm_lngth + idx); j++)
	       pool->data[to_seed].string[j] = pool->data[seed_idx].string[j];
	 }
	 idx += parm_lngth;
      }
      to_seed++;
   }
}
