/* ===========================================================================

         FILENAME:     swp_h_swap_utils.h

     PURPOSE:      A few useful functions when dealing with a swap structure

   ===========================================================================
 */
#ifndef SWP_H_SWAP_UTILS_H
#define SWP_H_SWAP_UTILS_H

/* ------------- Definition of a Double List: Dlist
 * --------------------------------- */
typedef struct {
  double *d;
  int len;
} Dlist;

/* ------------- Functions operating on a Dlist
 * --------------------------------- */

Err new_Dlist(int len, Dlist *t);

int free_inDlist(Dlist *dl);

/* -------------------------------------------------------------------- */

void Dlist_oper(Dlist l, double x, char oper);

Err Dlist_func(Dlist l1, Dlist l2, char oper, Dlist *l3);

Err Dlist_copy(Dlist l1, Dlist *l2);

/* -------------------------------------------------------------------- */

Err time_list(DateList dl, Date today, Dlist *t);

Err maturity_list(DateList fixing_date, Date today, Dlist *t);

Err cpn_list(BasisCode bc, SrtCompounding compd, DateList start_date,
             DateList end_date, DateList pay_date, double strike, double ini,
             double fin, StructType type, Dlist *cpn);

Err df_list(DateList pay_date, String ycname, Dlist *x);

Err fixing_list(DateList start_date, int spot_lag, DateList *fixing);

Err cvg_list(DateList start_date, DateList end_date, BasisCode bc, Dlist *x);

Err const_list(DateList dl, double constant, Dlist *l);

Err fwd_cash_list(DateList start_date, DateList end_date,
                  SrtBasisCode basiscode, String ycname, Dlist *fwd_cash);

Err spread_list(DateList start_date, DateList end_date, String ref_rate_code,
                Dlist *spread);

/* --------------------------------------------------------------------- */

typedef enum CMRateType { SWAP_RATE, TREASURY_RATE } CMRateType;

Err fwd_treas_list(DateList dl, SrtCurvePtr cmt_crv, Dlist *fwd_treas);

Err fwd_swap_list(DateList dl, SrtCurvePtr cmt_crv, Dlist *fwd_swap);

Err swap_vol_list(DateList dl, Dlist swap_rate, SrtCurvePtr cmt_crv,
                  Dlist *swap_vol);

Err treas_vol_list(DateList dl, Dlist swap_rate, Dlist treas_rate,
                   SrtCurvePtr cmt_crv, Dlist *treas_vol);

Err conv_fwd_list(long today, DateList dl, Dlist fwd_rate, Dlist vol,
                  SrtCurvePtr cmt_crv, CMRateType rate_type,
                  Dlist *conv_fwd_treas);

#endif
