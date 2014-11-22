#include "par.h"

/* ============================================================ *
 * Returns a copy of the array par.				*
 * ============================================================ */
par_t *copy_par_t(const par_t *par, int npar, error **err)
{
   int i;
   par_t *new;

   if (par==NULL) return NULL;

   new = malloc_err(sizeof(par_t)*npar, err);
   forwardError(*err, __LINE__, NULL);
   for (i=0; i<npar; i++) {
      new[i] = par[i];
   }

   return new;
}

/* ============================================================ *
 * Creates and returns an array par filled with the par_t types *
 * corresponding to the array of strings spar.			*
 * ============================================================ */
void spar_to_par(par_t **par, int npar, const char *spar[], error **err)
{
   int i, j;
   par_t *ploc;

   ploc = malloc_err(sizeof(par_t)*npar, err);      forwardError(*err, __LINE__,);

   for (i=0; i<npar; i++) {
      STRING2ENUM(ploc[i], spar[i], par_t, spar_t, j, Npar_t, err);
   }
   *par = ploc;
}
