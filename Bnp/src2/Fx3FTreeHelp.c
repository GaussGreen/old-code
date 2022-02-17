
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "math.h"

/* calculate both L and P in same function to save time */

static double Phi_Func(double x, double T, double s, double t)
{
double result;

	result = (exp(-x*(T-t)) - exp(-x*(T-s))) / x;

	return result;
}

static double Etha_Func(double x, double T, double s, double t)
{
double result;

	if (fabs(x) < 1.0e-8)
	{
		result = (t-s) * (T-(t+s)/2);
	}
	else
	{
		result = (t-s-Phi_Func(x, T, s, t)) / x;
	}

	return result;
}

static double B_Func(double x, double T, double s, double t)
{
double result;

	result = -(t*t-s*s)/2 + 1/x * ((t-1/x)*exp(-x*(T-t)) - (s-1/x)*exp(-x*(T-s)));

	return result;
}

static double Psi_Func(double x, double y, double T, double s, double t)
{
double result;

	if ((fabs(x) > 1.0e-8) && (fabs(y) > 1.0e-8))
	{
		/* both x and y are different from 0 */
		result = 1/(x*y) * (t-s - Phi_Func(x, T, s, t) - Phi_Func(y, T, s, t) + Phi_Func(x+y, T, s, t));
	}
	else
	{
		if (fabs(x) > 1.0e-8)
		{
			/* x different from 0, y equal to 0 */
			result = 1/x * (T * (t-s - Phi_Func(x, T, s, t)) + B_Func(x, T, s, t));
		}
		else
		{
			if (fabs(y) > 1.0e-8)
			{
				/* y different from 0, x equal to 0 */
				result = 1/y * (T * (t-s - Phi_Func(y, T, s, t)) + B_Func(y, T, s, t));
			}
			else
			{
				/* both equal to 0 */
				result = (pow(T-s, 3) - pow(T-t, 3)) / 3;
			}
		}
	}

	return result;
}


void fill_fwd_var(	
					long		nstp,
					double		*time,
					double		*date,
					double		*sig_dom,		
					double		*sig_for, 
					double		*sig_fx,
					double		dom_lam, 
					double		for_lam, 
					double		corr_dom_for,
					double		corr_dom_fx,
					double		corr_for_fx,
					char		*dom_yc,
					char		*for_yc,
					/*	To be allocated by caller, filled by the function */
					double		*dom_ifr,
					double		*dom_fwd,
					double		*dom_var,
					double		*for_ifr,
					double		*for_fwd,
					double		*for_var,
					double		*fx_fwd,
					double		*fx_var)
{
	int		i, j;
	double	a, b, c;					
	double	vol_dom, vol_for, vol_fx;
	double	var_dom, var_for, var_fx;
	double	fwd_dom, fwd_for;
	double	prev_mean, prev_rd, prev_rf;
	double	t, prev_t, t1, t2;
	double	aux1, aux2, L, P;

	dom_ifr[0] = swp_f_zr (date[0], date[1], dom_yc);
	for_ifr[0] = swp_f_zr (date[0], date[1], for_yc);

	fx_fwd[0] = 0.0;
	
	prev_mean = 0.0;
	prev_rd = dom_ifr[0];
	prev_rf = for_ifr[0];
	prev_t = 0.0;

	for (i=1; i<nstp ; i++)
	{
		if (i < nstp - 1)
		{
			dom_ifr[i] = swp_f_zr (date[i], date[i+1], dom_yc);
			for_ifr[i] = swp_f_zr (date[i], date[i+1], for_yc);
		}
		else
		{
			dom_ifr[i] = swp_f_zr (date[i], date[i]+1, dom_yc);
			for_ifr[i] = swp_f_zr (date[i], date[i]+1, for_yc);
		}
		
		t = time[i];		
		var_fx = 0.0;		
		fwd_dom = 0.0;
		var_dom = 0.0;
		fwd_for = 0.0;
		var_for = 0.0;

		t1 = 0.0;
		for (j=0; j<i ; j++)
		{
			t2 = time[j+1];
			
			/* on ]j, j+1], we use the volatility at time j+1 */
			vol_dom = sig_dom[j+1];
			vol_for = sig_for[j+1];
			vol_fx = sig_fx[j+1];

			/* coefficient of sig_fx 2 */
			a = t2 - t1;
			/* coefficient of sig_fx   */
			b = -2.0 * corr_for_fx * vol_for * Etha_Func(for_lam, t, t1, t2) 
				+ 2.0 * corr_dom_fx * vol_dom * Etha_Func(dom_lam, t, t1, t2);
			/* constant coefficient    */		   	
			c = vol_for * vol_for * Psi_Func(for_lam, for_lam, t, t1, t2)
				+ vol_dom * vol_dom * Psi_Func(dom_lam, dom_lam, t, t1, t2)
				- 2.0 * corr_dom_for * vol_dom * vol_for * Psi_Func(for_lam, dom_lam, t, t1, t2);
			
			var_fx += a * vol_fx * vol_fx + b * vol_fx + c;

			/* domestic fwd and var */
			aux1 = exp(-dom_lam * (t - t1));
			aux2 = exp(-dom_lam * (t - t2));

			L = 1.0 / (dom_lam * dom_lam) * (aux2 * (1.0 - 0.5 * aux2) - aux1 * (1.0 - 0.5 * aux1));
			P = 1.0 / (2.0 * dom_lam) * (aux2 * aux2 - aux1 * aux1);
			
			fwd_dom += vol_dom * vol_dom * L;
			var_dom += vol_dom * vol_dom * P;

			/* foreign fwd and var */
			aux1 = exp(-for_lam *(t-t1));
			aux2 = exp(-for_lam *(t-t2));

			L = 1.0 / (for_lam * for_lam) * (aux2 * (1.0 - 0.5 * aux2) - aux1 * (1.0 - 0.5 * aux1));
			P = 1.0/ (2.0 * for_lam) * (aux2 * aux2 - aux1 * aux1);
	
			fwd_for += vol_for * vol_for * L;
			var_for += vol_for * vol_for * P;

			t1 = t2;
		}
		dom_fwd[i] = fwd_dom;
		dom_var[i] = var_dom;

		for_fwd[i] = fwd_for;
		for_var[i] = var_for;

		fx_var[i] = var_fx;		
		fx_fwd[i] = fx_fwd[i-1] + (prev_rd - prev_rf - vol_fx * vol_fx / 2.0) * (t - prev_t);
	
		/* Compute prev */
		prev_rd = fwd_dom + dom_ifr[i];
		prev_rf = fwd_for + for_ifr[i];
		prev_t = t;
	}
}

void fill_fwd_var_corr(	
					long		nstp,
					double		*time,
					double		*date,
					double		*sig_dom,		
					double		*sig_for, 
					double		*sig_fx,
					double		dom_lam, 
					double		for_lam, 
					double		*correl_dom_for,
					double		*correl_dom_fx,
					double		*correl_for_fx,
					char		*dom_yc,
					char		*for_yc,
					/*	To be allocated by caller, filled by the function */
					double		*dom_ifr,
					double		*dom_fwd,
					double		*dom_var,
					double		*for_ifr,
					double		*for_fwd,
					double		*for_var,
					double		*fx_fwd,
					double		*fx_var)
{
	int		i, j;
	double	a, b, c;					
	double	vol_dom, vol_for, vol_fx;
	double	var_dom, var_for, var_fx;
	double	fwd_dom, fwd_for;
	double	corr_dom_for, corr_dom_fx, corr_for_fx;
	double	prev_mean, prev_rd, prev_rf;
	double	t, prev_t, t1, t2;
	double	aux1, aux2, L, P;

	dom_ifr[0] = swp_f_zr (date[0], date[1], dom_yc);
	for_ifr[0] = swp_f_zr (date[0], date[1], for_yc);

	fx_fwd[0] = 0.0;
	
	prev_mean = 0.0;
	prev_rd = dom_ifr[0];
	prev_rf = for_ifr[0];
	prev_t = 0.0;

	for (i=1; i<nstp ; i++)
	{
		if (i < nstp - 1)
		{
			dom_ifr[i] = swp_f_zr (date[i], date[i+1], dom_yc);
			for_ifr[i] = swp_f_zr (date[i], date[i+1], for_yc);
		}
		else
		{
			dom_ifr[i] = swp_f_zr (date[i], date[i]+1, dom_yc);
			for_ifr[i] = swp_f_zr (date[i], date[i]+1, for_yc);
		}
		
		t = time[i];		
		var_fx = 0.0;		
		fwd_dom = 0.0;
		var_dom = 0.0;
		fwd_for = 0.0;
		var_for = 0.0;

		t1 = 0.0;
		for (j=0; j<i ; j++)
		{
			t2 = time[j+1];
			
			/* on ]j, j+1], we use the volatility at time j+1 */
			vol_dom = sig_dom[j+1];
			vol_for = sig_for[j+1];
			vol_fx = sig_fx[j+1];
			corr_dom_for = correl_dom_for[j+1];
			corr_dom_fx = correl_dom_fx[j+1];
			corr_for_fx = correl_for_fx[j+1];

			/* coefficient of sig_fx 2 */
			a = t2 - t1;
			/* coefficient of sig_fx   */
			b = -2.0 * corr_for_fx * vol_for * Etha_Func(for_lam, t, t1, t2) 
				+ 2.0 * corr_dom_fx * vol_dom * Etha_Func(dom_lam, t, t1, t2);
			/* constant coefficient    */		   	
			c = vol_for * vol_for * Psi_Func(for_lam, for_lam, t, t1, t2)
				+ vol_dom * vol_dom * Psi_Func(dom_lam, dom_lam, t, t1, t2)
				- 2.0 * corr_dom_for * vol_dom * vol_for * Psi_Func(for_lam, dom_lam, t, t1, t2);
			
			var_fx += a * vol_fx * vol_fx + b * vol_fx + c;

			/* domestic fwd and var */
			aux1 = exp(-dom_lam * (t - t1));
			aux2 = exp(-dom_lam * (t - t2));

			L = 1.0 / (dom_lam * dom_lam) * (aux2 * (1.0 - 0.5 * aux2) - aux1 * (1.0 - 0.5 * aux1));
			P = 1.0 / (2.0 * dom_lam) * (aux2 * aux2 - aux1 * aux1);
			
			fwd_dom += vol_dom * vol_dom * L;
			var_dom += vol_dom * vol_dom * P;

			/* foreign fwd and var */
			aux1 = exp(-for_lam *(t-t1));
			aux2 = exp(-for_lam *(t-t2));

			L = 1.0 / (for_lam * for_lam) * (aux2 * (1.0 - 0.5 * aux2) - aux1 * (1.0 - 0.5 * aux1));
			P = 1.0/ (2.0 * for_lam) * (aux2 * aux2 - aux1 * aux1);
	
			fwd_for += vol_for * vol_for * L;
			var_for += vol_for * vol_for * P;

			t1 = t2;
		}
		dom_fwd[i] = fwd_dom;
		dom_var[i] = var_dom;

		for_fwd[i] = fwd_for;
		for_var[i] = var_for;

		fx_var[i] = var_fx;		
		fx_fwd[i] = fx_fwd[i-1] + (prev_rd - prev_rf - vol_fx * vol_fx / 2.0) * (t - prev_t);
	
		/* Compute prev */
		prev_rd = fwd_dom + dom_ifr[i];
		prev_rf = fwd_for + for_ifr[i];
		prev_t = t;
	}
}