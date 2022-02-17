/************************************************************************/
/* MODULE NAME : 	SRT_H_2F_COR_CAL.H				*/
/* Version 		1.00						*/
/* Author		Olivier VAN EYSEREN				*/
/************************************************************************/

double  srt_f_hist_cor_fit(double *rate_date, double *hist_func, double *weight,
			int fit_lambda_1, double lambda_1_guess,  int ndata,
                        int iter_numb,
			double converg_date,
			double hist_start_date, double hist_end_date);

double display_corr_struct( String s );

double tf_corr_value(double vol_ratio, double lambda_1, double lambda_2, 
			double rho,
			double rate_date1, double rate_date2,
			double hist_start, double hist_end, double fwd_date);

double tf_corr_func(double vol_ratio, double lambda_1, double lambda_2, 
			double rho,
			double rate_date, 
			double hist_start, double hist_end, double fwd_date);
