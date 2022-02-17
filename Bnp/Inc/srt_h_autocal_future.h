#ifndef LGMFUTURE_H_
#define LGMFUTURE_H_

Err srt_f_price_future_lgm (
			long			  today,
			long              *fut_dates, 
			long              num_fut_dates,
			char              *yc_name,
			char              *ref_rate_code,
			char              *calibrated_underlying,
			int				  result_type,
			double            *fut_price);

Err srt_f_autocal_future (
			long              *fut_dates, 
			long              num_fut_dates,
			Err			(*GetVol)(double dStart, double dEnd, double dStrike, 
						double dForward, double dSpread, double *pdBsVol),
							  /* function that returns the volatility of the reference swptns */
		    String            logNormStr,
			char              *mdl_name, 
			char              *yc_name,
			char              *ref_rate_code,
			double            tau, 
			double            alpha, 
			double            gamma, 
			double            rho,
			long              numpaths,
			double            *fut_price);

Err srt_f_autocal_future_input_instr (
			long              *fut_dates, 
			long              num_fut_dates,
			long			  ninstruments,
			char			  **type_str,
			long			  *start_dates,
			long			  *end_dates,
			char			  **compd_str,
			char 			  **basis_str,
			char			  **recpay_str,
			char			  **refratecode,
			char              *mdl_name, 
			char              *yc_name,
			char			  *volcurvename,
			char              *ref_rate_code,
			double            tau, 
			double            alpha, 
			double            gamma, 
			double            rho,
			double            *fut_price,
			double			  *instr_price,
			int				  bumpinst);

Err srt_f_lgmcashfuture (
			Ddate             start_date,
			Ddate             end_date,
			SrtBasisCode      basis, 
			char              *und_name,
			double            * future_price);

Err srt_f_cheycashfuture (
			Ddate             start_date, 
			Ddate             end_date,
			char             *basis,
			char             *und_name,
			long              numpaths,
			double           *future_price);
#endif