#ifndef _SRTGRFNAUTOCAL_H_
#define _SRTGRFNAUTOCAL_H_

typedef struct {
  long num_stp;
  long num_paths;
  long num_stpx;
  int do_pecs;
} fx3datc_param, *FX3DATC_PARAM;

void fx3datc_set_default_param(FX3DATC_PARAM param);

Err GrfnAutocalCaller(
    /*	Today */
    long today, double spot_fx,
    /*	Get Cash Vol function */
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Domestic market */
    char *dom_und_name, char *dom_ccy,
    char *dom_yc_name,        /*	Name of the yield curve */
    char *dom_vol_curve_name, /*	Name of the market vol curve */
    char *dom_ref_rate_name,  /*	Name of the reference rate */
    char *dom_instr_freq,     /*	Frequency and basis of instruments */
    char *dom_instr_basis,
    /*	Domestic Model */
    double dom_lambda, /*	Lambda: may NOT be changed in the process */
    /*	Foreign market */
    char *for_und_name, char *for_ccy,
    char *for_yc_name,        /*	Name of the yield curve */
    char *for_vol_curve_name, /*	Name of the market vol curve */
    char *for_ref_rate_name,  /*	Name of the reference rate */
    char *for_instr_freq,     /*	Frequency and basis of instruments */
    char *for_instr_basis,
    /*	Domestic Model */
    double for_lambda, /*	Lambda: may NOT be changed in the process */
    /*	Fx market */
    long *fx_mkt_vol_date, /*	Option maturity dates */
    double *fx_mkt_vol,    /*	Option BS vol */
    int num_fx_mkt_vol,    /*	Number of options */
    /*	Fx model */
    double *corr_times, double *correl_dom_for, double *correl_dom_fx,
    double *correl_for_fx, long corr_n_times,
    /*	Structure */
    /*	If ex_date is NULL      ,
    exercise dates will be generated 2bd before start */
    int num_ex_dates, /*	Exercise dates      ,
                                                  all supposed to be on or
                     after today */
    long *ex_date,    /*	Supposed to be sorted */
    int *cal_date,    /*	1: use ex_date as calibration date      , 0: don't */
    /*	If calibration dates are overwritten */
    int owr_cal_dates,    /*	1: use overwritten dates      , 0: use ex_date
                           */
    int num_owr_dates,    /*	Number of overwritten dates */
    long *owr_date,       /*	Overwritten dates */
    char **dom_end_tenor, /*	Tenors of the underlying instruments or "DIAG"
                           */
    char **for_end_tenor, long end_date, /*	End date for diagonal */
    double *dom_strike,                  /*	Domestic strikes 0: ATM */
    double *for_strike,                  /*	Foreign strikes 0: ATM */
    /*	Calibration parameters */
    CPD_DIAG_CALIB_PARAM calib_param,
    /*	GRFN Tableau */
    long tableauRows, long tableauCols, char ***tableauStrings,
    int **tableauMask, long auxWidth, long *auxLen, double **aux,
    /* End of Day Flags */
    int is_end_of_day_fixing, int is_end_of_day_payment,
    /*	Pricing method */
    int mth, /*	0: Tree      , 1: Monte-Carlo */
    /*	Pricing parameters */
    FX3DATC_PARAM pricing_param,
    /*	Calibration output */
    int *num_sig, /*	Answer */
    double **sig_time, double **dom_sig, double **for_sig, int *num_fx_vol,
    double **fx_vol_time, double **fx_vol,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA
        dom_inst_data, /*	NULL = don't save calibration instrument data */
    CPD_CALIB_INST_DATA
        for_inst_data, /*	NULL = don't save calibration instrument data */
    /*	Pricing output */
    int *num_prod, double **prod_val,
    /*	Whether to keep the data */
    int keep_calib_data, int keep_und,
    /*	Quanto or Hybrid ?	*/
    int quantoflag);

#endif