#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "math.h"

static double Z_Func(double x, double t1, double t2)
{
	return ((exp(x * t2) - exp(x * t1))/ x);
}

static double Phi_Func(double x, double T, double s, double t)
{
double result;

	result = (exp (-x * (T - t)) - exp (-x * (T - s))) / x;

	return result;
}

static double Etha_Func(double x, double T, double s, double t)
{
double result;

	result = (t - s - Phi_Func (x, T, s, t)) / x;

	return result;
}

static double B_Func(double x, double T, double s, double t)
{
double result;

	result = - (t * t - s * s) / 2.0 
		+ 1.0 / x * ((t - 1.0 / x) * exp (-x * (T - t)) 
		- (s - 1.0 / x) * exp(-x * (T - s)));

	return result;
}

static double Psi_Func(double x, double y, double T, double s, double t)
{
double result;

	result = 1.0 / (x * y) 
		* (t - s 
		- Phi_Func (x, T, s, t) 
		- Phi_Func (y, T, s, t) 
		+ Phi_Func (x + y, T, s, t));

	return result;
}

static double Psi2_Func(double x, double y, double Tx, double Ty, double s, double t)
{
double result;

	result = 1.0 / (x * y) 
		* (t - s 
		- Phi_Func(x, Tx, s, t) 
		- Phi_Func(y, Ty, s, t) 
		+ exp (-x * (Tx - Ty)) * Phi_Func (x + y, Ty, s, t));

	return result;
}

Err fill_Betamc_init(
					double	*date,
					double	*time,
					long	nb_dates,
					double	*sig_dates, 
					long	nb_sig_dates,
					double	*sig_curve_dom,
					double	lda_dom,
					double	*sig_curve_for,
					double	lda_for,
					double	*sig_curve_fx,
					char	*dom_yc,
					char	*for_yc,
					double	*dom_ifr,
					double	*dom_std,
					double	*dom_phi,					
					double	*for_ifr,					
					double	*for_std,
					double	*for_phi,
					double	*fx_std
					)
{
double	sig_dom2, sig_for2, sig_fx2;
double	var_dom, var_for, var_fx;
double	phi_dom, phi_for;
double	y_dom, y_for;
double	T1, T2, start_date, end_date, start_mat, end_mat;

int		i, k;
long	StartIndex, EndIndex;
Err		err = NULL;

	dom_phi[0] = 0;
	phi_dom = 0;
	for_phi[0] = 0;
	phi_for = 0;
		
	start_date = date[0];
	start_mat = time[0];
	StartIndex = Get_Index (start_mat, sig_dates, nb_sig_dates);
	
	for (k=0; k<nb_dates-1; k++)
	{
		end_date = date[k+1];					
		end_mat = time[k+1];
		EndIndex = Get_Index (end_mat, sig_dates, nb_sig_dates);	

		var_dom = 0.0;		
		var_for = 0.0;
		var_fx = 0.0;
		
		dom_ifr[k] = swp_f_zr (start_date, end_date, dom_yc);
		for_ifr[k] = swp_f_zr (start_date, end_date, for_yc);
				
		for (i = StartIndex; i < EndIndex + 1; i++)
		{
			if (i > StartIndex)
			{
				T1 = sig_dates[i-1];
			}
			else
			{
				/* First part */
				T1 = start_mat;
			}
			

			if (i == EndIndex || StartIndex == EndIndex)
			{
				/* Last part */
				T2 = end_mat;
			}
			else
			{
				T2 = sig_dates[i];
			}

			sig_dom2 = sig_curve_dom[i];
			sig_dom2 *= sig_dom2;			
			sig_for2 = sig_curve_for[i];							
			sig_for2 *= sig_for2; 						
			sig_fx2 = sig_curve_fx[i];		
			sig_fx2 *= sig_fx2;
						
			y_dom = Z_Func(2 * lda_dom, T1, T2);
			y_for = Z_Func(2 * lda_for, T1, T2);
			
			var_dom += sig_dom2 * y_dom;							
			var_for += sig_for2 * y_for;
			var_fx += sig_fx2 * (T2 - T1);
		}
		
		phi_dom += var_dom;
		phi_for += var_for;

		/* Standard deviation */		
		dom_std[k+1] = exp(-lda_dom * end_mat) * sqrt(var_dom);		
		for_std[k+1] = exp(-lda_for * end_mat) * sqrt(var_for);		
		fx_std[k+1] = sqrt(var_fx);

		/* Phi */
		dom_phi[k+1] = phi_dom * exp(-2 * lda_dom * end_mat);
		for_phi[k+1] = phi_for * exp(-2 * lda_for * end_mat);

		start_date = end_date;
		start_mat = end_mat;
		StartIndex = EndIndex;
	}
	
	dom_ifr[k] = swp_f_zr (date[k], date[k]+1, dom_yc);
	for_ifr[k] = swp_f_zr (date[k], date[k]+1, for_yc);
	
	return err;
}