#ifndef __SMM_H__
#define __SMM_H__

typedef struct _smm_swap {
  // fixed leg:
  int nfix;        // number of coupons
  long *fix_d;     // dates
  double *fix_t;   // times
  double *fix_cvg; // coverages
  double *fix_df;  // discount factors

  // floating	leg:
  int nflt;        // number of coupons
  long *flt_d;     // dates
  double *flt_t;   // times
  double *flt_cvg; // coverages
  double *flt_df;  // discount factors
  double *flt_spr; // spreads

  // fix  , start & end:
  long fxg_date, start_date, end_date;
  double fxg_time, start_time, end_time;
  double df_fix, df_start, df_end;

  double swap0, lvl0;     // swap & level values today
  double lvlc0;           // level cash at forward
  double cms_fwd;         // CMS forward rate
  double swpn_strikes[3]; // Swaption strikes (for smile calib)
  double swpn_values[3];  // Swaption values (for smile calib)

  double beta, fwd1, fwd2, sigbeta1, sigbeta2, pi; // BMM parameters

} smm_swap;

typedef struct _smm_params {
  double nstd, maxerror; // for marginals cumulative calculation
  int cum_npts;          // number of cumulative discretization points asked
  SrtMCSamType mc_type;  // type of the MC method used
  long npaths;           // number of MC paths
  double smile_nstd; // number of standard deviations to get strikes for smile
                     // calibration
  double vol_mult;   // vol mult factor for tests
  int do_match;      // 1: match market (weighted MC)  , 0: original simulation
  double sv_min;     // relative minimum of a singular value
  int vol_hyp; // bit0 : converging / sliding on sigma-beta  , bit1: c/s on sabr
               // params
  int use_fee; // 1: match IV  , 0: no fees
  int otc_type; // 0: IV  , 1: call  , 2: put
} smm_params;

typedef struct _smm_swaptions_desc {
  long *cf_idx;  // for each instrument - index of the first cashflow
  int *ncf;      // for each instrument - number of cash flows
  double *cf;    // all cash flows
  long *dfs_idx; // for each cash flow - index of the corresponding df

  SrtReceiverType pay_rec;

  // test CMS:
  double *K;       // strikes (for each swaption)
  long *cffix_idx; // index of the first fixed leg cashflow (for each swaption)

} smm_swaptions_desc;

typedef struct _cif_fund_cpn {
  double cpn;
  long df_idx;
} cif_fund_cpn;

typedef struct _cif_exo_cpn {
  // for CMS calculation:
  int nfix, nflt;
  double *fix_cpn, *flt_cpn;
  long *fix_df_idx, *flt_df_idx;

  long start, end;

  // call string:
  double wcst, wfwd;
  int nstrikes;
  double *strikes, *weights;

  // pay info:
  double cvg;
  long df_idx;
} cif_exo_cpn;

typedef struct _cif_desc {
  int fund_ncpn, exo_ncpn;
  cif_fund_cpn *fund_cpn;
  cif_exo_cpn *exo_cpn;

  int *ex_fund_idx, *ex_exo_idx;
} cif_desc;

Err smm_swaptions_init_product(
    SProductDesc *g, // The product description structure to initialise
    char *yc,        // Yield curve name
    char *freq_str,  // Underlying swaps fixed freq
    char *bas_str,   // Underlying swaps fixed basis
    char *ref_rate,  // Underlying swaps floating refrate
    char *pay_rec,   // Receiver / payer
    int is_cms,      // 1: cms options  , 0: swaptions
    long start,      // Common swaps start
    int nends,       // Number of underlying swaps
    long *ends,      // Underlying swaps ends
    int *nK,         // Number of strikes per swap
    double **K);     // Strikes (per swap)

void smm_swaptions_free_product(SProductDesc *g);

Err smm_otc_price(char *yc,       // Market yield curve
                  char *vc,       // Market vol curve
                  char *freq_str, // Fixed freq of swaps used in SMM model
                  char *bas_str,  // Fixed basis of swaps used in SMM model
                  char *ref_rate, // Floating refrate of swaps used in SMM model
                  char *corr_cube,    // Swaps correlation cube name
                  smm_params *params, // Numerical parameters
                  SProductDesc *g,    // Product to price
                  int iex,            // Exercise date index
                  double *pv,         // PV (per instrument)
                  double *stddev);    // Standard deviation (per instrument)

Err smm_dfs_init_product(SProductDesc *g, char *yc, long ex_date, int ndates,
                         long *dates);
void smm_dfs_free_product(SProductDesc *g);

#endif // #ifndef __SMM_H__