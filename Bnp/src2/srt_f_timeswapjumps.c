/**************************************************************************************************
 *																								  *
 *		Function: srt_f_timeswapjumps
 **
 *																								  *
 *		Calculates the PV of the fixed leg of a  range swap with minimum
 *daily coupon and		  * an approximation of the PV when there is a
 *minimum PERIOD coupon.						  *
 *																								  *
 *		Note that the daily PV is always bigger than the period PV. The
 *two prices computed here  * are upper bounds to the real PV.
 **
 *																								  *
 *		If the TimeSwap flag is set to YES the PV of the fixed leg of a
 *regular TimeSwap		  * is computed.
 **
 *																								  *
 *		see the FIRST Document : Time Swap and Range Swap with a Minimum
 *Coupon					  * for details.
 **
 *																								  *
 *		Author: Ezra Nahum
 ** Last modified Jan 2000.
 **
 *																								  *
 ***************************************************************************************************/

#include "OPFNCTNS.H"
#include "UTALLHDR.H"
#include "math.h"
#include "num_h_gamma.h"
#include "srt_h_all.h"
#include "swp_h_cms.h"
#include <srt_h_resetable.h"
#include <swp_h_all.h"
#include <swp_h_cmsopt.h"
#include <swp_h_external_fct.h"
#include <swp_h_vol.h"

#define NMAXSTRIKES 10

/*Free Memory functions*/

Err free_resetableperiod(double **paramperiod, double *periodmat,
                         double *DRStsStart, double *DRSts, double *isfriday,
                         double nfixdays, double *driftinterp);

Err free_resetablegeneral_ts(double **paramatrix, int n_resetdates) {
  if (paramatrix)
    free_dmatrix(paramatrix, 1, n_resetdates, 1, 5);
  paramatrix = NULL;

  return NULL;
}

/* Function which calculates the value of a Fra */

Err get_fra(long Fradate, int spotlag, char *cRefRateCode,
            char *szYieldCurveName, double *dFra);

/* Function which calibrates the parameters of the jump model to a caplet smile
corresponding to a specific matirity date */

Err calibtosmile(long today, long dDate, char *cMarketId, char *szVolCurveName,
                 char *cRefRateCode, char *szYieldCurveName, double *param,
                 double *dFra, long *ia);

/* General calibration function which calibrates the parameters to the smiles of
the resetdates of the resetable as well as the smile of the end date of the
product (going from the longer maturity to the shorter one)*/

Err resetablecalib(long today, int n_resetdates, double *resetdates,
                   char *cMarketId, char *szVolCurveName, char *cRefRateCode,
                   char *szYieldCurveName, double *paraminitguess,
                   double **paramatrix);

/* function which linearly interpolates the parameters obtained from the smiles
of the beginning and end of the period   , for a maturity in between*/

Err get_param(double Fradate, int indexStart, double *resetdates,
              double **paramatrix, double *paramfra);

Err get_param_first_leg(double Fradate, double *resetdates, double **paramatrix,
                        double *paramfra) {
  Err err = NULL;
  double datespread;

  /*linear interpolation of the parameters*/

  datespread = (resetdates[2] - resetdates[1]);

  paramfra[1] = (paramatrix[2][1] - paramatrix[1][1]) *
                    (Fradate - resetdates[1]) / datespread +
                paramatrix[1][1];

  paramfra[2] = (paramatrix[2][2] - paramatrix[1][2]) *
                    (Fradate - resetdates[1]) / datespread +
                paramatrix[1][2];

  paramfra[3] = paramatrix[1][3];

  paramfra[4] = (paramatrix[2][4] - paramatrix[1][4]) *
                    (Fradate - resetdates[1]) / datespread +
                paramatrix[1][4];

  paramfra[5] = paramatrix[1][5];
  return err;
}

/* gets the number of fixing days within a period */

Err get_NumFixingDays(int indexStart, double *resetdates, double *nfixdays);

/* Gets the maturity date of each of the fixing days within a period and assigns
the value 1 in the vector isfriday if the day is a friday or a 0 otherwise */

Err get_periodmatandisfriday(int indexStart, double *resetdates,
                             double nfixdays, double *periodmat,
                             double *isfriday);

/* Gets the DRS value for each fixing day within a period */

Err get_DRSts(long today, long Paydate, int nfixdays, double *periodmat,
              char *cRefRateCode, char *szYieldCurveName, double *DRSts);

/* Using the Merton diffusion with fixed parameters  , computes the value of a
DRs corresponding to the solution of the SDE it satisfies*/

Err DRSCalc(double *wn1n2, int nfixdays, double resetmat, double **periodparam,
            double *drift, double *DRSts);

/* Main function*/

Err srt_f_timeswapjumps(int n_resetdates, double *resetdates, double *initparam,
                        char *FloatCoupon, double margin, double upperbarrier,
                        double lowerbarrier, char *cRefRateCode,
                        char *cMarketId, char *szYieldCurveName,
                        char *szVolCurveName, double *pv) {

  Err err = NULL;
  int perindex, j;
  double nfixdays;
  double **paramatrix;
  double *periodmat;
  double *isfriday;
  double *DRSts, *DRStsStart;
  double **paramperiod;
  double width;
  double prix;
  double eps;
  double ntotaldays;
  double drift;
  double *driftinterp;
  double *maturities;
  double cvg;
  double dFradrift;
  int spotlag;
  double callspread = 0.001;
  /* double *discfact;*/
  long today;
  SrtBasisCode float_basis;
  SrtCompounding float_compounding;
  SrtCrvPtr yldcrv;

  yldcrv = lookup_curve(szYieldCurveName);
  today = get_today_from_curve(yldcrv);
  spotlag = get_spotlag_from_curve(yldcrv);
  width = upperbarrier - lowerbarrier;

  /*Get details for the index Fra*/
  if (err = swp_f_get_ref_rate_details(cRefRateCode, &float_basis,
                                       &float_compounding)) {
    smessage("Error in swp_f_get_ref_rate_details");
  }

  /*memory allocation*/
  paramatrix = dmatrix(1, n_resetdates, 1, 5);

  /* calibrate the parameters to the caplets smiles of the resetdates and end
   * date */

  err =
      resetablecalib(today, n_resetdates, resetdates, cMarketId, szVolCurveName,
                     cRefRateCode, szYieldCurveName, initparam, paramatrix);
  if (err) {
    free_resetablegeneral_ts(paramatrix, n_resetdates);
    return err;
  }

  /* Big loop begins over the number of periods in the resetable */

  for (perindex = 1; perindex < n_resetdates; perindex++) {

    err = get_NumFixingDays(perindex, resetdates, &nfixdays);
    if (err) {
      free_resetablegeneral_ts(paramatrix, n_resetdates);
      return err;
    }

    ntotaldays = resetdates[perindex + 1] - resetdates[perindex];
    pv[perindex] = 0;

    /*memory allocation*/
    periodmat = dvector(1, (long)nfixdays);
    DRStsStart = dvector(1, (long)nfixdays);
    isfriday = dvector(1, (long)nfixdays);
    paramperiod = dmatrix(1, (long)nfixdays, 1, 5);
    DRSts = dvector(1, (long)nfixdays);
    driftinterp = dvector(1, (long)nfixdays);
    maturities = dvector(1, (long)nfixdays);

    err = get_periodmatandisfriday(perindex, resetdates, nfixdays, periodmat,
                                   isfriday);
    if (err) {
      free_resetableperiod(paramperiod, periodmat, DRStsStart, DRSts, isfriday,
                           nfixdays, driftinterp);
      free_resetablegeneral_ts(paramatrix, n_resetdates);
      return err;
    }

    err = get_DRSts(today, (long)resetdates[perindex + 1], (int)nfixdays,
                    periodmat, cRefRateCode, szYieldCurveName, DRStsStart);
    if (err) {
      free_resetableperiod(paramperiod, periodmat, DRStsStart, DRSts, isfriday,
                           nfixdays, driftinterp);
      free_resetablegeneral_ts(paramatrix, n_resetdates);
      return err;
    }

    cvg = coverage((long)resetdates[perindex + 1],
                   add_unit((long)resetdates[perindex + 1],
                            (int)(12 / float_compounding), SRT_MONTH,
                            MODIFIED_SUCCEEDING),
                   float_basis);
    /* here  , we get the parameters corresponding to each fixing day within the
     period and put them in
     paramperiod which is a matrix with nfixdays rows and 5 columns (for the 5
     parameters)*/

    err = get_fra((long)resetdates[perindex + 1], spotlag, cRefRateCode,
                  szYieldCurveName, &dFradrift);

    drift = paramatrix[perindex][1] * cvg * dFradrift / (1 + cvg * dFradrift);

    for (j = 1; j <= nfixdays; j++) {
      err = get_param(periodmat[j], perindex, resetdates, paramatrix,
                      paramperiod[j]);

      driftinterp[j] = drift * (periodmat[j] - resetdates[perindex]) /
                       (resetdates[perindex + 1] - resetdates[perindex]);

      if (err) {
        free_resetableperiod(paramperiod, periodmat, DRStsStart, DRSts,
                             isfriday, nfixdays, driftinterp);
        free_resetablegeneral_ts(paramatrix, n_resetdates);
        return err;
      }
    }

    /* maximum values for the number of positive and negative jumps are computed
    for the resetdate (start date) of the period in question */
    for (j = 1; j <= nfixdays; j++) {
      maturities[j] = (periodmat[j] - today) / 365.0;
      DRSts[j] = DRStsStart[j] * exp(driftinterp[j] * maturities[j]);
    }

    eps = upperbarrier - DRSts[1];

    prix = PrixTimeSwap(eps, DRSts, (int)nfixdays, maturities, isfriday,
                        paramperiod, width, callspread);

    if (err) {
      free_resetableperiod(paramperiod, periodmat, DRStsStart, DRSts, isfriday,
                           nfixdays, driftinterp);
      free_resetablegeneral_ts(paramatrix, n_resetdates);
      return err;
    }

    if (strcmp(FloatCoupon, "YES") == 0) {
      pv[perindex] += (DRSts[1] + margin) * (-prix) / ntotaldays;
    }

    else {
      pv[perindex] += (margin) * (-prix) / ntotaldays;
    }

    free_resetableperiod(paramperiod, periodmat, DRStsStart, DRSts, isfriday,
                         nfixdays, driftinterp);
    free_dvector(maturities, 1, (long)nfixdays);
  }

  free_resetablegeneral_ts(paramatrix, n_resetdates);

  return err;
}
