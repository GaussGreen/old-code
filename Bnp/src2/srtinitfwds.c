/* ==========================================================================
   FILENAME: SrtInitFwds.C

   PURPOSE: Initialise a forward object stores it in the Curves list
   ========================================================================== */

#include "SrtAccess.h"
#include "grf_h_all.h"
#include "srt_h_all.h"
#include "srt_h_fwdcurve.h"
#include "srt_h_repocurve.h"

char *SrtInitDividends(long today, char *dvdName, char *ccy, int nRows,
                       int nCols, double **dvdVals) /* [Date][Value] */
{
  Err err = NULL;

  /* Call the initialisation function */
  err = srt_f_init_DividendCurve(today, ccy, dvdVals, nCols, nRows, dvdName);

  /* Return a success message or an error */
  return err;

} /* END char* SrtInitDividends(...) */

char *SrtInitRepos(long today, char *repo_name, char *ccy, int nrows, int ncols,
                   double **repo_vals) /*[Date][Value] */
{
  Err err = NULL;

  err = srt_f_init_EQRepoCurve(today, ccy, repo_vals, ncols, nrows, repo_name);
  return err;
}
