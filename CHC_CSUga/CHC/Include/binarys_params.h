  /* **************************************************************** */
  /*          Parameters needed by Binary Transform routines          */
  /* **************************************************************** */

#ifndef _BINARYS_PARAMS_H_
#define _BINARYS_PARAMS_H_


  int Nb_params;                    /* number of parameters in the run */

  struct p_list *Work_params;       /* Parameters describing the genes */

  double *Params;                   /* Work space for gene translation */

  char *Grey_gene;                  /* Work Space for de-greying a gene */


#endif
