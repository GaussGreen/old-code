/* =================================================================================
   FILENAME :     SrtInitCorrMatrix.c

   PURPOSE:       The Function to call to set up a Correlation Matrix in Grfn
   =================================================================================
 */

#include "SrtAccess.h"

/*  ----------------------------------------------------------------------
    Static Declaration:
                very important:stores the complete correlation term structure
        ----------------------------------------------------------------------
 */

static SrtCorrLstPtr __corr_list = NULL;

/* ------------------------------------------------------------------------
   Buils te Full Correlation List for any Grfn Calculation
   ----------------------------------------------------------------------- */
Err SrtInitCorrelationMatrix(int ndates, int ncorr, double **correl,
                             double *dates, String **und_names) {
  Err err = NULL;

  err = srt_f_init_Corr_TermStruct(ndates, ncorr, correl, dates, und_names,
                                   &__corr_list);

  return err;
}

/* -----------------------------------------------------------------------
   Returns the adress of the static used to store the correlation list.
   ----------------------------------------------------------------------- */
SrtCorrLstPtr srt_f_GetTheCorrelationList(void) { return __corr_list; }

/* -----------------------------------------------------------------------
   Creates the STATIC SrtCorrLstPtr __corr_list  , and gives it a name...
   ----------------------------------------------------------------------- */
Err create_correlation_list(String corr_list_name) {
  Err err;

  if (err = srt_f_corrlstcreate(&__corr_list, corr_list_name)) {
    __corr_list = NULL;
    return (err);
  }
  return NULL;
}

/* -----------------------------------------------------------------------
   Destroys the STATIC SrtCorrLstPtr __corr_list
   ----------------------------------------------------------------------- */
Err destroy_correlation_list() {
  Err err;

  if (err = srt_f_corrlstdelete(&__corr_list))
    return (err);

  __corr_list = NULL;
  return NULL;
}

/* ----------------------------------------------------------------------- */
