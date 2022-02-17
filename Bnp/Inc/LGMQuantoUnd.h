#ifndef __LGM_QUANTO_UND_H
#define __LGM_QUANTO_UND_H

//	Structures and functions for the underlying and its term structures

//	Underlying term structs
typedef struct {
  char name[256];
  long today;
  char dom_yc[256];
  char for_yc[256];
  char dom_vc[256];
  char for_vc[256];
  char dom_ref[256]; //	Reference rate used for getting the vol
  char dom_swap_freq[256];
  char dom_swap_basis[256];
  char for_ref[256]; //	Reference rate used for getting the vol
  char for_swap_freq[256];
  char for_swap_basis[256];

  long sigma_n;
  double *sigma_date;
  double *sigma_time;
  double *dom_sigma;  //	Domestic Vol
  double *for_sigma;  //	Foreign Vol
  double *fx_sigma;   //	FX Vol
  double *domfor_rho; //	Dom/For Correl
  double *quanto_rho; //	FX/For Correl

  double dom_lambda;
  double for_lambda;

  //	Calibration instrument data
  int has_inst_data;
  cpd_calib_inst_data dom_inst_data;
  cpd_calib_inst_data for_inst_data;

  //	Value of the fwd receiver CCF
  int has_fwd_iv;
  int nb_fwdiv;
  double *exercise_date;
  double *market_fwdiv;
  double *model_fwdiv;
  double *extra_fees;

} lgmQto_und, *LGMQTO_UND;

/*	Fill underlying structure from a predefined underlying */
Err fill_lgmQto_und(char *und3dfx, char *dom_vc, char *for_vc, char *dom_ref,
                    char *dom_swap_freq, char *dom_swap_basis, char *for_ref,
                    char *for_swap_freq, char *for_swap_basis, LGMQTO_UND und);

/*	Fill underlying structure from calibration instruments */
/*
Err calib_lgmQto_und(
long		today  ,
//	EOD Flag
int			eod_flag  ,
//	0: I  , 1: E
char		*yc  , //	yc char		*vc  ,
//	vc (only if calib) char		*ref  ,
//	ref rate (only if calib) char		*swap_freq  ,
//	swap freq (only if calib) char		*swap_basis  ,
//	swap basis (only if calib) int			lam_ts  ,
//	0: use unique lambda  , 1: use ts double		lambda  ,
//	lambda if unique int			tsnlam  ,
//	number of lambdas if ts double		*tslamtime  ,
//	lambda times i.e. (date - today) / 365 if ts
double		*tslam  ,
//	corresponding lambdas if ts double		alpha  ,
//	alpha
double		gamma  , //	gamma
double		rho  , //	rho
//	Calib params
int			force_atm  ,
//	force atm calib double		max_std_long  , double
max_std_short  , int			fix_lambda  ,
//	0: calib lambda to cap  , 1: fix lambda calib
                                                                                        //			to diagonal
int			one_f_equi  ,
//	1F equivalent flag:
                                                                                        //			if set to 1  , then 2F lambda will calibrate
                                                                                        //			to the cap priced within calibrated 1F
                                                                                        //			with the given lambda
int			skip_last  ,
//	If 1  , the last option is disregarded
                                                                                        //			and the forward volatility is flat from option
                                                                                        //			n-1
//	End of calib params
CTSQUANTO_STR		ctsquanto  ,
//	structure Err			(*get_cash_vol)(
//	function to get IR cash vol from the markets char	*vol_curve_name
, double	start_date  , double	end_date  , double	cash_strike  ,
                                int		zero  ,
                                char	*ref_rate_name  ,
                                double	*vol  ,
                                double	*power)  ,
LGMQTO_UND		und  ,
long		spread_vol_n  ,
double		*spread_vol_time  ,
double		*spread_vol_floor  ,
double		*spread_vol_cap  ,
int			cvg_sv  ,
int			is_corr);

*/

void copy_lgmQto_und(LGMQTO_UND src, LGMQTO_UND dest);

Err free_lgmQto_und(LGMQTO_UND und);

/*	Arguments to the adi function */
/*	----------------------------- */

typedef struct {
  int nstp;
  double *time;
  double *date;
  int nstpx;
  double *sig_time;
  double *domsig;
  double *forsig;
  double *fxsig;
  int nb_sig;
  double *domforrho;
  double *quantorho;
  double domlambda;
  double forlambda;
  void **void_prm;
  int *is_event;
  double *dom_ifr;
  double *for_ifr;
  char dom_yc[256];
  char for_yc[256];
} lgmQto_adi_arg, *LGMQTO_ADI_ARG;

Err lgmQto_fill_adi_arg(
    LGMQTO_UND und, int nb_event_dates, double *event_time, void *prm,
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Required number of steps */
    int req_stp, int req_stpx, LGMQTO_ADI_ARG adi_arg);

Err lgmQto_free_adi_arg(LGMQTO_ADI_ARG adi_arg);

#endif