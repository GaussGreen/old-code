/******************************************************************************
  SYSTEM    :SRT
  SUB-SYSTEM    :GRFN
  MODULE NAME   :GRF_f_MIDAT_CLSDFRM.C

  AUTHOR        :Eric Auld
  CREATED   :april 1994
  VERSION       :1.00
  DESCRIPTION   :


  AMENDMENTS
    Reference   :
    Author      : E AULD
    Date        : 21 july 94
    Description : changed function so that bond strike is assumed
        to be payed at first coupon date after today  , if today is
        not a coupon date.

******************************************************************************/

#include "grf_h_public.h"
#include "srt_h_all.h"
#include "srt_h_grfclsdfrm.h"
#include "stdlib.h"

static void grfn_mk_midat_str(String opt_str, double bond_strike,
                              SrtReceiverType rec_pay, SRT_Boolean pv) {
  String pvstring;
  if (pv == SRT_YES)
    pvstring = "PV[0]";
  else
    pvstring = "0.0";

  if (rec_pay == SRT_RECEIVER) {
    sprintf(opt_str,
            "max(pvrng(0  ,1  ,a[2  ,i])-%.12lf*df(now  ,a[2  ,i])  ,%s)",
            bond_strike, pvstring);
  } else {
    sprintf(opt_str,
            "max(-pvrng(0  ,1  ,a[2  ,i])+%.12lf*df(now  ,a[2  ,i])  ,%s)",
            bond_strike, pvstring);
  }
}

/*****************************************************************************
   FUNCTION     : grfn_mk_midat_GrfnCells
   DESCRIPTION      : generate the Grfn Tableau describing
            a midat swaption. American swaption
            is considered a type of midat swaption.

   AMENDMENTS       :
    Reference   :
    Author          :
    Date            :
    Description     :

******************************************************************************/

static Err
grfn_mk_midat_GrfnCells(long num_exercise_dates, Date *exercise_dates,
                        double *exercise_premiums, SrtReceiverType rec_pay,
                        Date clcn_date, long *num_eventdates, Date **eventdates,
                        long *nrows, long *ncols, GrfnCell ***sprdsht) {
  Err err = NULL;
  int i;

  /*** only include exercise dates after clcn_date ***/
  if (num_exercise_dates <= 0)
    return serror("no valid exercise dates");

  /*** allocate space ***/
  *ncols = 1;
  *nrows = *num_eventdates = num_exercise_dates;
  *eventdates = (Date *)srt_calloc(num_exercise_dates, sizeof(Date));
  *sprdsht = GrfnCellmatrix(*nrows, *ncols, GRFN_DEF_ARGBUFSZ);

  /*** set grfn event dates ***/
  memcpy(*eventdates, exercise_dates, sizeof(Date) * num_exercise_dates);

  /*** enter payoff string for each date ***/
  for (i = 0; i < num_exercise_dates - 1; i++) {
    grfn_mk_midat_str((*sprdsht)[i][0].sval, exercise_premiums[i], rec_pay,
                      SRT_YES);
    (*sprdsht)[i][0].type = GRFNSCELL;
    if (err)
      return err;
  }
  grfn_mk_midat_str((*sprdsht)[i][0].sval, exercise_premiums[i], rec_pay,
                    SRT_NO);
  (*sprdsht)[i][0].type = GRFNSCELL;
  if (err)
    return err;

  return err;
}

/*****************************************************************************
   FUNCTION     : grfn_midat_clsdfrm
   DESCRIPTION      : PRICE A MIDAT SWAPTION USING GRFN.

   AMENDMENTS       :
    Reference   :
    Author          :
    Date            :
    Description     :

******************************************************************************/
Err grfn_midat_clsdfrm(long num_exercise_dates,    /* len of next two arrays */
                       Date *exercise_dates,       /* dates when option
                                     can be exercise_d */
                       Date *exercise_start_dates, /* dates when bonds start  ,
                                 >= correponding exercise
                                 dates. */
                       double *exercise_premiums,  /* amounts that must be
                                     payed to exercise_ options */
                       long num_prod_dates, Date *prod_dates, double *prod_cfs,

                       SrtReceiverType rec_pay,
                       SrtUndPtr und,           /* name of underlying to use */
                       SrtGrfnParam *grfnparam, /* model and implementation
                           details */
                       double *answer           /* value returned here*/
) {
  Err err = NULL;
  long ncols = 0, nrows = 0, num_eventdates = 0;
  Date *eventdates = NULL;
  GrfnCell **sprdsht = NULL;
  double **aux;
  long *auxlen;
  long auxwidth = 0;
  int i, inc;
  SrtIOStruct *iolist; /* list of requests for prices */
  Date today;

  today = get_today_from_underlying(und);

  inc = (long)(grfnparam->end_of_day_flg ? 1 : 0);

  err = srt_f_IOstructcreate(&iolist, "");
  /* Create a list with a request for the price */

  /*** create auxiliary vars ***/
  if (num_prod_dates < 1) {
    *answer = 0;
    return NULL;
  }

  auxwidth = 3;
  auxlen = srt_calloc(auxwidth, sizeof(long));

  /* [MN. Bug Fix] Remove aux ranges corresponding to (start) exercise dates
                   before today - Use end-of-day flag
  */

  while ((num_exercise_dates > 0) && (exercise_dates[0] < (today + inc))) {

    num_exercise_dates--; /* Decrement number of prod dates */

    exercise_start_dates++; /* Move to next exercise start date */

    exercise_dates++;

    exercise_premiums++;
  }

  aux = (double **)srt_calloc(auxwidth, sizeof(double *));

  aux[0] = (double *)srt_calloc(num_prod_dates, sizeof(double));
  aux[1] = (double *)srt_calloc(num_prod_dates, sizeof(double));
  aux[2] = (double *)srt_calloc(num_exercise_dates, sizeof(double));

  for (i = 0; i < num_prod_dates; i++) {
    aux[0][i] = (double)prod_dates[i];
  }
  memcpy(aux[1], prod_cfs, num_prod_dates * sizeof(double));
  auxlen[0] = auxlen[1] = num_prod_dates;
  for (i = 0; i < num_exercise_dates; i++) {
    aux[2][i] = (double)exercise_start_dates[i];
  }
  auxlen[2] = num_exercise_dates;

  /* CONSTRUCT A GRFN SPREADSHEET THAT DESCRIBES THE MIDAT SWAPTION */
  err = grfn_mk_midat_GrfnCells(
      num_exercise_dates, exercise_dates, exercise_premiums, rec_pay, today,
      &num_eventdates, &eventdates, &nrows, &ncols, &sprdsht);

  if (!err)
    err =
        srt_f_grfn(und, grfnparam, num_eventdates, &eventdates, &nrows, &ncols,
                   &sprdsht, 0, 0, auxwidth, auxlen, aux, iolist, 0, 0);

  err = srt_f_IOstructgetpremiumval(*iolist, answer);

  if (eventdates)
    srt_free(eventdates);
  if (sprdsht) {
    grfn_free_GrfnCellmatrix(sprdsht, nrows, ncols);
  }
  /** free auxiliary vars **/
  srt_free(aux[0]);
  srt_free(aux[1]);
  srt_free(aux[2]);
  srt_free(aux);
  srt_free(auxlen);

  err = srt_f_IOstructfree(&iolist);

  return err;
}
