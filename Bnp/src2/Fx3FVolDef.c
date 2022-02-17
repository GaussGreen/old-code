/* ==========================================================================
   FILE_NAME:	Fx3FVolDef.c

   PURPOSE:		Defines the long term Fx volatility using the 3F model
				5 strategies
				1.)	Calibrate all
				2.)	Calibrate all, optimise fvol stationarity in slr
				3.) Calibrate partial, extrapolate with slr
				4.) Calibrate partial, best fit rest with slr
				5.) Optimise fvol stationarity under bid/offer constraint

   DATE:		08/07/00
   
   AUTHOR:		A.S.
   ========================================================================== */

#include "srt_h_all.h"
#include "SrtAccess.h"
#include "srt_h_allFx3F.h"
#include "num_h_levenberg.h"
#include "opfnctns.h"
#include "math.h"
#include "old_cpd_calib_wrapper.h"

/* Static data */

/* Name of the Initial Underlyings */
char			*DOM_UND = "Domestic_Init";
char			*FOR_UND = "Foreign_Init";

/* Just useful numbers */
double			ONE;
double			SMALL;
double			BIG;

/* Implied volatility term structure */
int				num_opt_;
double			*opt_exe_;
double			*opt_mat_;
double			*ivol_;
double			*bid_;
double			*offer_;

/* Model parameters and bounds */
double			*fvol_;
double			*sig_dom_;
double			sig_dom_min_;
double			sig_dom_max_;
double			*sig_for_;
double			sig_for_min_;
double			sig_for_max_;
double			*lam_dom_;
double			lam_dom_min_;
double			lam_dom_max_;
double			*lam_for_;
double			lam_for_min_;
double			lam_for_max_;
double			*rho_dom_for_;
double			rho_dom_for_min_;
double			rho_dom_for_max_;
double			*rho_dom_fx_;
double			rho_dom_fx_min_;
double			rho_dom_fx_max_;
double			*rho_for_fx_;
double			rho_for_fx_min_;
double			rho_for_fx_max_;
int				force_eq_sig_;
int				force_eq_lam_;
int				ext_;
double			exta_;
double			extb_;

/* Optimisation parameters */

/* Equations */

/* Total number of equations */
int				ndata_;
/* Needs recalibration */
int				do_calib_;
/* Number of instruments for (partial) calibration */
int				num_calib_;

/* First part: bounds */
int				do_bnd_;
int				ndata_bnd_first_;
int				ndata_bnd_num_;
double			**ndata_bnd_val_;
double			*ndata_bnd_min_;
double			*ndata_bnd_max_;

/* Second part: implied volatility */
int				do_ivol_;
int				ndata_ivol_first_;
int				ndata_ivol_first_opt_;
int				ndata_ivol_num_;

/* Third part: implied volatility bounds */
int				do_ivol_bnd_;
int				ndata_ivol_bnd_first_;
int				ndata_ivol_bnd_first_opt_;
int				ndata_ivol_bnd_num_;

/* Forth part: forward volatility stationarity */
int				do_fvol_;
int				ndata_fvol_first_;
int				ndata_fvol_first_opt_;
int				ndata_fvol_num_;

/* The parameters */

int				ma_;
double			**a_;

/* Levenberg evaluation function */
static Err funcs1 (		double	x,
						double	a[],
						double	*yfit,
						int		ma)
{
	int i, j;
	double *temp_fvol;
	double temp_ivol;
	double max, min, maxmv;
	Err err = NULL;

	/*	Parameter set */

	for (i=0; i<ma; i++)
	{
		(*(a_[i])) = a[i+1];
	}

	if (force_eq_sig_)
	{
		*sig_for_ = *sig_dom_;
	}

	if (force_eq_lam_)
	{
		*lam_for_ = *lam_dom_;
	}

	/*	Equation number */
	j = DTOL (x);

	/*	Bounds */
	if (do_bnd_ && j >= ndata_bnd_first_ && j < ndata_bnd_first_ + ndata_bnd_num_)
	{
		/*	Number of bnd equation */
		j -= ndata_bnd_first_;

		if (*(ndata_bnd_val_[j]) < ndata_bnd_min_[j])
		{
			*yfit = ndata_bnd_min_[j] - *(ndata_bnd_val_[j]);
		}
		else
		if (*(ndata_bnd_val_[j]) > ndata_bnd_max_[j])
		{
			*yfit = *(ndata_bnd_val_[j]) - ndata_bnd_max_[j];
		}
		else
		{
			*yfit = 0.0;
		}

		*yfit *= 100.0;
	}

	/*	Implied volatility */
	if (do_ivol_ && j >= ndata_ivol_first_ && j < ndata_ivol_first_ + ndata_ivol_num_)
	{
		/*	Number of ivol equation */
		j -= ndata_ivol_first_;

		if (do_calib_ && num_calib_)
		{
			/*	Calibration */
			err = Fx3DtsCalibration (	opt_exe_,
										opt_mat_,
										ivol_,
										num_calib_,
										&ONE,
										1,
										sig_dom_,
										*lam_dom_,
										sig_for_,
										*lam_for_,
										*rho_dom_for_,
										*rho_dom_fx_,
										*rho_for_fx_,
										&temp_fvol);

			if (err)
			{
				/*	Couldn't calibrate: return very BIG value */
				*yfit = 1000.0;
				return NULL;
			}

			for (i=0; i<num_calib_; i++)
			{
				fvol_[i] = temp_fvol[i];
			}
			/*	If calibration is partial, extrapolate */
			for (i=num_calib_; i<num_opt_; i++)
			{
				if (ext_)
				{
					fvol_[i] = fvol_[num_calib_-1] * exp (-exta_ * (opt_mat_[i] - opt_mat_[num_calib_-1]))
						+ extb_ * (1.0 - exp (-exta_ * (opt_mat_[i] - opt_mat_[num_calib_-1])));
				}
				else
				{
					fvol_[i] = fvol_[i-1];
				}
			}

			free (temp_fvol);
		}

		/*	Implied volatility valuation */

		err = Fx3DtsImpliedVol (	opt_exe_[ndata_ivol_first_opt_+j],
									0.0,
									opt_mat_[ndata_ivol_first_opt_+j],
									&ONE,									
									1,
									sig_dom_,
									*lam_dom_,
									sig_for_,
									*lam_for_,
									opt_mat_,
									fvol_,
									num_opt_,
									*rho_dom_for_,
									*rho_dom_fx_, 
									*rho_for_fx_,
									yfit);
		
		if (err)
		{
			return err;
		}

		*yfit -= ivol_[ndata_ivol_first_opt_+j];
	}

	/*	Implied volatility bounds */
	if (do_ivol_bnd_ && j >= ndata_ivol_bnd_first_ && j < ndata_ivol_bnd_first_ + ndata_ivol_bnd_num_)
	{
		/*	Number of ivol bnd equation */
		j -= ndata_ivol_bnd_first_;

		if (do_calib_ && num_calib_)
		{
			/*	Calibration */
			err = Fx3DtsCalibration (	
										opt_exe_,
										opt_mat_,
										ivol_,
										num_calib_,
										&ONE,
										1,
										sig_dom_,										
										*lam_dom_,										
										sig_for_,										
										*lam_for_,
										*rho_dom_for_,
										*rho_dom_fx_,
										*rho_for_fx_,
										&temp_fvol);

			if (err)
			{
				/*	Couldn't calibrate: return very BIG value */
				*yfit = 1000.0;
				return NULL;
			}

			for (i=0; i<num_calib_; i++)
			{
				fvol_[i] = temp_fvol[i];
			}
			/*	If calibration is partial, extrapolate */
			for (i=num_calib_; i<num_opt_; i++)
			{
				if (ext_)
				{
					fvol_[i] = fvol_[num_calib_-1] * exp (-exta_ * (opt_mat_[i] - opt_mat_[num_calib_-1]))
						+ extb_ * (1.0 - exp (-exta_ * (opt_mat_[i] - opt_mat_[num_calib_-1])));
				}
				else
				{
					fvol_[i] = fvol_[i-1];
				}
			}

			free (temp_fvol);
		}

		/*	Implied volatility valuation */

		err = Fx3DtsImpliedVol (	opt_exe_[ndata_ivol_bnd_first_opt_+j],
									0.0,
									opt_mat_[ndata_ivol_bnd_first_opt_+j],
									&ONE,
									1,
									sig_dom_,
									*lam_dom_,
									sig_for_,
									*lam_for_,
									opt_mat_,
									fvol_,
									num_opt_,
									*rho_dom_for_,
									*rho_dom_fx_,
									*rho_for_fx_,
									&temp_ivol);
		
		if (err)
		{
			return err;
		}

		if (temp_ivol < bid_[ndata_ivol_bnd_first_opt_+j])
		{
			*yfit = bid_[ndata_ivol_bnd_first_opt_+j] - temp_ivol;
		}
		else
		if (temp_ivol > offer_[ndata_ivol_bnd_first_opt_+j])
		{
			*yfit = temp_ivol - offer_[ndata_ivol_bnd_first_opt_+j];
		}
		else
		{
			*yfit = 0.0;
		}

		*yfit *= 100.0;
	}

	/*	Forward volatility stationarity */
	if (do_fvol_ && j >= ndata_fvol_first_ && j < ndata_fvol_first_ + ndata_fvol_num_)
	{
		/*	Number of fvol equation */
		j -= ndata_fvol_first_;

		if (do_calib_ && num_calib_)
		{
			/*	Calibration */
			err = Fx3DtsCalibration (	
										opt_exe_,
										opt_mat_,
										ivol_,
										num_calib_,
										&ONE,
										1,
										sig_dom_,
										*lam_dom_,
										sig_for_,
										*lam_for_,
										*rho_dom_for_,
										*rho_dom_fx_,
										*rho_for_fx_,
										&temp_fvol);

			if (err)
			{
				/*	Couldn't calibrate: return very BIG value */
				*yfit = 1000.0;
				return NULL;
			}

			for (i=0; i<num_calib_; i++)
			{
				fvol_[i] = temp_fvol[i];
			}
			/*	If calibration is partial, extrapolate */
			for (i=num_calib_; i<num_opt_; i++)
			{
				if (ext_)
				{
					fvol_[i] = fvol_[num_calib_-1] * exp (-exta_ * (opt_mat_[i] - opt_mat_[num_calib_-1]))
						+ extb_ * (1.0 - exp (-exta_ * (opt_mat_[i] - opt_mat_[num_calib_-1])));
				}
				else
				{
					fvol_[i] = fvol_[i-1];
				}
			}

			free (temp_fvol);
		}

		max = SMALL;
		min = BIG;
		maxmv = 0.0;
		
		for (i=ndata_fvol_first_opt_; i<num_opt_; i++)
		{
			if (i > ndata_fvol_first_opt_ && fabs (fvol_[i] - fvol_[i-1]) > maxmv)
			{	
				maxmv = fabs (fvol_[i] - fvol_[i-1]);
			}

			if (fvol_[i] > max)
			{
				max = fvol_[i];
			}

			if (fvol_[i] < min)
			{
				min = fvol_[i];
			}
		} 
		
		*yfit = max - min + maxmv;
	}

	return NULL;
}

/* Levenberg derivative function */
#define DER_SHIFT 0.0001
static Err funcs2 (		double	x,
						double	a[],
						double	*yfit,
						double	dyda[],
						int		ma)
{
	int i;
	double y;
	Err err;

	err = funcs1 (x, a, yfit, ma);
	if (err)
	{
		return err;
	}

	for (i=0; i<ma; i++)
	{
		a[i+1] += DER_SHIFT;
		err = funcs1 (x, a, &y, ma);
		if (err)
		{
			return err;
		}
		dyda[i+1] = (y - *yfit) / DER_SHIFT;
		a[i+1] -= DER_SHIFT;
	}

	return NULL;
}

/* Main function */
/* Defines the long term Fx volatility using the 3F model
				5 strategies
				1.)	Calibrate all
				2.)	Calibrate all, optimise fvol stationarity in slr
				3.) Calibrate partial, extrapolate with slr
				4.) Calibrate partial, best fit rest with slr
				5.) Optimise fvol stationarity under bid/offer constraint */
Err Fx3DDefImpVol(		/*	Total number of options */
						int		num_opt,
						/*	Number of options to calibrate to */
						int		num_calib,
						/*	Implied volatility */
						double	*opt_mat,
						double	*ivol,
						double	*bid,
						double	*offer,
						/*	Model parameters and bounds */
						double	*fvol,
						double	*sig_dom,
						double	sig_dom_min,
						double	sig_dom_max,
						double	*sig_for,
						double	sig_for_min,
						double	sig_for_max,
						double	*lam_dom,
						double	lam_dom_min,
						double	lam_dom_max,
						double	*lam_for,
						double	lam_for_min,
						double	lam_for_max,
						double	*rho_dom_for,
						double	rho_dom_for_min,
						double	rho_dom_for_max,
						double	*rho_dom_fx,
						double	rho_dom_fx_min,
						double	rho_dom_fx_max,
						double	*rho_for_fx,
						double	rho_for_fx_min,
						double	rho_for_fx_max,
						int		ext,
						double	exta,
						double	extb,
						/*	Number of iterations */
						int		num_iter,
						/*	What to optimise on */
						int		opt_sigma,
						int		opt_lambda,
						int		opt_correl,
						int		force_eq_sig,
						int		force_eq_lam,
						/*	Method 1 to 5 */
						int		method)
{
	double *data = NULL, *target = NULL, *weight = NULL, *parm = NULL;	
	double *temp_fvol = NULL;
	double chisq;
	int i, j;
	Err err = NULL;

	/*	A. Setup global pointers and numbers */

	ndata_bnd_val_ = NULL;
	ndata_bnd_min_ = NULL;
	ndata_bnd_max_ = NULL;
	a_ = NULL;

	ONE = 1.0;
	SMALL = 0.0001;
	BIG = 10.0;

	/*	B. Setup global data */

	num_opt_ = num_opt;
	opt_mat_ = opt_mat;
	opt_exe_ = opt_mat;
	ivol_ = ivol;
	bid_ = bid;
	offer_ = offer;
	fvol_ = fvol;
	
	sig_dom_ = sig_dom;
	sig_dom_min_ = sig_dom_min;
	sig_dom_max_ = sig_dom_max;
	sig_for_ = sig_for;
	sig_for_min_ = sig_for_min;
	sig_for_max_ = sig_for_max;
	lam_dom_ = lam_dom;
	lam_dom_min_ = lam_dom_min;
	lam_dom_max_ = lam_dom_max;
	lam_for_ = lam_for;
	lam_for_min_ = lam_for_min;
	lam_for_max_ = lam_for_max;
	rho_dom_for_ = rho_dom_for;
	rho_dom_for_min_ = rho_dom_for_min;
	rho_dom_for_max_ = rho_dom_for_max;
	rho_dom_fx_ = rho_dom_fx;
	rho_dom_fx_min_ = rho_dom_fx_min;
	rho_dom_fx_max_ = rho_dom_fx_max;
	rho_for_fx_ = rho_for_fx;
	rho_for_fx_min_ = rho_for_fx_min;
	rho_for_fx_max_ = rho_for_fx_max;
	force_eq_sig_ = force_eq_sig;
	force_eq_lam_ = force_eq_lam;
	ext_ = ext;
	exta_ = exta;
	extb_ = extb;

	/*	C. Setup calibration parameters */
	
	if (method == 1 || method == 2)
	{
		/*	Full calibration */
		do_calib_ = 1;
		num_calib_ = num_opt;
	}
	else
	if (method == 3 || method == 4)
	{
		/*	Partial calibration */
		do_calib_ = 1;
		num_calib_ = num_calib;
		if (num_calib_ < 1)
		{
			return "At least one option must be calibrated for methods 3 and 4";
		}
	}
	else
	{
		/*	No calibration */
		do_calib_ = 0;
		num_calib_ = 0;
	}

	/*	D. Setup bounds */
	
	j = 0;

	if (method == 1 || method == 3)
	{
		/*	No optimisation on slr */
		do_bnd_ = 0;
		ndata_bnd_first_ = 0;
		ndata_bnd_num_ = 0;
	}
	else
	{
		ndata_bnd_first_ = 0;
		ndata_bnd_num_ = 0;

		if (opt_sigma) 
		{
			if (force_eq_sig)
			{
				ndata_bnd_num_ += 1;
			}
			else
			{
				ndata_bnd_num_ += 2;
			}
		}

		if (opt_lambda)
		{
			if (force_eq_lam)
			{
				ndata_bnd_num_ += 1;
			}
			else
			{
				ndata_bnd_num_ += 2;
			}
		}
		
		if (opt_correl) 
		{
			ndata_bnd_num_ += 3;
		}

		if (ndata_bnd_num_ >= 1)
		{
			do_bnd_ = 1;

			ndata_bnd_val_ = (double**) calloc (ndata_bnd_num_, sizeof (double*));
			ndata_bnd_min_ = (double*) calloc (ndata_bnd_num_, sizeof (double));
			ndata_bnd_max_ = (double*) calloc (ndata_bnd_num_, sizeof (double));

			i = 0;
			if (opt_sigma)
			{
				ndata_bnd_val_[i] = sig_dom_;
				ndata_bnd_min_[i] = sig_dom_min_;
				ndata_bnd_max_[i] = sig_dom_max_;
				
				i += 1;

				if (!force_eq_sig)
				{
					ndata_bnd_val_[i] = sig_for_;
					ndata_bnd_min_[i] = sig_for_min_;
					ndata_bnd_max_[i] = sig_for_max_;
					
					i += 1;
				}
			}

			if (opt_lambda)
			{
				ndata_bnd_val_[i] = lam_dom_;
				ndata_bnd_min_[i] = lam_dom_min_;
				ndata_bnd_max_[i] = lam_dom_max_;
				
				i += 1;

				if (!force_eq_lam)
				{
					ndata_bnd_val_[i] = lam_for_;
					ndata_bnd_min_[i] = lam_for_min_;
					ndata_bnd_max_[i] = lam_for_max_;
					
					i += 1;
				}
			}

			if (opt_correl)
			{
				ndata_bnd_val_[i] = rho_dom_for_;
				ndata_bnd_min_[i] = rho_dom_for_min_;
				ndata_bnd_max_[i] = rho_dom_for_max_;
				
				ndata_bnd_val_[i+1] = rho_dom_fx_;
				ndata_bnd_min_[i+1] = rho_dom_fx_min_;
				ndata_bnd_max_[i+1] = rho_dom_fx_max_;

				ndata_bnd_val_[i+2] = rho_for_fx_;
				ndata_bnd_min_[i+2] = rho_for_fx_min_;
				ndata_bnd_max_[i+2] = rho_for_fx_max_;
			}	
			
		}
		else
		{
			/*	No optimisation on slr */
			do_bnd_ = 0;
			ndata_bnd_first_ = 0;
			ndata_bnd_num_ = 0;
		}
	}

	j += ndata_bnd_num_;

	/*	E. Setup implied vol */
	
	if (method == 4 && ndata_bnd_num_ >= 1)
	{
		do_ivol_ = 1;
		ndata_ivol_first_ = j;
		ndata_ivol_first_opt_ = num_calib_;
		ndata_ivol_num_ = num_opt_ - num_calib_;
	}
	else
	{
		/*	No optimisation on implied vol */
		do_ivol_ = 0;
		ndata_ivol_first_ = 0;
		ndata_ivol_first_opt_ = 0;
		ndata_ivol_num_ = 0;
	}

	j += ndata_ivol_num_;

	/*	F. Setup implied vol bounds */
	
	if (method == 5)
	{
		do_ivol_bnd_ = 1;
		ndata_ivol_bnd_first_ = j;
		ndata_ivol_bnd_first_opt_ = 0;
		ndata_ivol_bnd_num_ = num_opt - ndata_ivol_bnd_first_opt_;
	}
	else
	{
		/*	No bounds on implied vol */
		do_ivol_bnd_ = 0;
		ndata_ivol_bnd_first_ = 0;
		ndata_ivol_bnd_first_opt_ = 0;
		ndata_ivol_bnd_num_ = 0;
	}

	j += ndata_ivol_bnd_num_;

	/*	G. Setup fvol stationarity */
	
	if (method == 2 || method == 5)
	{
		do_fvol_ = 1;
		ndata_fvol_first_ = j;
		ndata_fvol_first_opt_ = 0;
		ndata_fvol_num_ = 1;
	}
	else
	{
		do_fvol_ = 0;
		ndata_fvol_first_ = 0;
		ndata_fvol_first_opt_ = 0;
		ndata_fvol_num_ = 0;
	}

	j += ndata_fvol_num_;
	
	ndata_ = j;

	/*	H.	Setup parameters */

	if (method == 1 || method == 3)
	{
		/*	No optimisation */
		ma_ = 0;
	}
	else
	if (method == 2 || method == 4)
	{
		/*	Optimise in slr only */
		ma_ = ndata_bnd_num_;
	}
	else
	{
		/*	Optimise in slr and fvol */
		ma_ = ndata_bnd_num_ + num_opt_;
	}
	
	if (ma_ >= 1)
	{
		a_ = (double**) calloc (ma_, sizeof (double*));
		
		if (ndata_bnd_num_ >= 1)
		{
			for (i=0; i<ndata_bnd_num_; i++)
			{
				a_[i] = ndata_bnd_val_[i];
			}
		}

		if (ma_ > ndata_bnd_num_)
		{
			for (i=0; i<num_opt_; i++)
			{
				a_[ndata_bnd_num_+i] = &(fvol_[i]);
			}
		}
	}

	/*	I. Setup local data */	
	
	data = (double*) calloc (ndata_, sizeof (double));
	for (i=0; i<ndata_; i++)
	{
		data[i] = 1.0 * i;
	}

	target = (double*) calloc (ndata_, sizeof (double));
	for (i=0; i<ndata_; i++)
	{
		target[i] = 0.0;
	}

	weight = (double*) calloc (ndata_, sizeof (double));
	for (i=0; i<ndata_; i++)
	{
		weight[i] = 1.0;
	}

	parm = (double*) calloc (ma_, sizeof (double));

	/*	J. Calibrate fvol */

	if (method == 1 || method == 2 || method == 5)
	{
		/*	Full calibration */

		err = Fx3DtsCalibration (	opt_exe_,
									opt_mat_,
									ivol_,
									num_opt_,
									&ONE,
									1,
									sig_dom_,
									*lam_dom_,
									sig_for_,
									*lam_for_,
									*rho_dom_for_,
									*rho_dom_fx_,
									*rho_for_fx_,
									&temp_fvol);

		if (err)
		{
			goto FREE_RETURN;
		}

		for (i=0; i<num_opt_; i++)
		{
			fvol_[i] = temp_fvol[i];
		}

		if (temp_fvol) free (temp_fvol);
	}
	else
	{
		/*	Partial calibration */

		err = Fx3DtsCalibration (	
									opt_exe_,
									opt_mat_,
									ivol_,
									num_calib_,
									&ONE,
									1,
									sig_dom_,
									*lam_dom_,
									sig_for_,
									*lam_for_,
									*rho_dom_for_,
									*rho_dom_fx_,
									*rho_for_fx_,
									&temp_fvol);

		if (err)
		{
			goto FREE_RETURN;
		}

		for (i=0; i<num_calib_; i++)
		{
			fvol_[i] = temp_fvol[i];
		}
		
		/*	Calibration is partial, extrapolate */
		for (i=num_calib_; i<num_opt_; i++)
		{
			if (ext_)
			{
				fvol_[i] = fvol_[num_calib_-1] * exp (-exta_ * (opt_mat_[i] - opt_mat_[num_calib_-1]))
					+ extb_ * (1.0 - exp (-exta_ * (opt_mat_[i] - opt_mat_[num_calib_-1])));
			}
			else
			{
				fvol_[i] = fvol_[i-1];
			}
		}

		if (temp_fvol) free (temp_fvol);
	}

	/*	K. Call Levenberg */
	
	for (i=0; i<ma_; i++)
	{
		parm[i] = *(a_[i]);
	}

	if (ndata_ >= 1 && ma_ >= 1)
	{
		err = levenberg_marquardt (	data-1,
									target-1,
									weight-1,
									ndata_,
									parm-1,
									ma_,
									num_iter,
									funcs2,
									&chisq);
	}

	if (force_eq_sig_)
	{
		*sig_for_ = *sig_dom_;
	}

	if (force_eq_lam_)
	{
		*lam_for_ = *lam_dom_;
	}
	
	if (!err)
	{
		for (i=0; i<ma_; i++)
		{
			*(a_[i]) = parm[i];
		}

		for (i=0; i<num_opt_; i++)
		{
			err = Fx3DtsImpliedVol (	opt_exe_[i],
										0.0,
										opt_mat_[i],
										&ONE,
										1,
										sig_dom_,
										*lam_dom_,
										sig_for_,
										*lam_for_,
										opt_mat_,
										fvol_,
										num_opt_,
										*rho_dom_for_,
										*rho_dom_fx_,
										*rho_for_fx_,
										&(ivol_[i]));
		}
	}

FREE_RETURN:

	if (a_) free (a_);
	if (ndata_bnd_val_) free (ndata_bnd_val_);
	if (ndata_bnd_min_) free (ndata_bnd_min_);
	if (ndata_bnd_max_) free (ndata_bnd_max_);
	if (data) free (data);
	if (target) free (target);
	if (weight) free (weight);
	if (parm) free (parm);

	return err;
}

static Err Get_Calib_LGM_Underlying(								
										SrtUndPtr   undptr,
										char		*und_name,
										double		tau,
										long		today,
										long		spot_date,
										double		maturity_date,
										char		*yc,
										char		*ref_rate,
										char		*basis,
										char		*cpd,
										char		*ccy,
										char		*call_freq,
										long		last_opt_date,
										Err			(*pfGetVol)(double dStart, double dEnd, double dStrike, 
																double dForward, double dSpread, double *pdBsVol),
										char		*log_nor,
										SrtGrfnParam       *grfnparam
									)

{
/* params of calboot */
double		**vol_curve			= NULL,
			**tau_curve			= NULL,			
			*undPrice			= NULL,
			*strike				= NULL,
			*bondStrike			= NULL;
			
char		**paramStrings		= NULL,
			**valueStrings		= NULL,			
			**refRateStr		= NULL;

SrtReceiverType	*recPayStr		= NULL;
			
		
int			nInput, i,
			*type				= NULL;

long		call_start, call_end;
double		vol, Spread;
double		maturity_date_adj;
SrtDateList	dateList;
SwapDP      *sdp				= NULL;

SwapDP		Swap;
Err			err = NULL;

		/* initialization */
		dateList.len = 0;

		/* Check if we interpolate or if we extrapolate */
		if (maturity_date < last_opt_date)
		{
			maturity_date_adj = last_opt_date;
		}
		else
		{
			maturity_date_adj = maturity_date;
		}

		if (maturity_date_adj < (today + 365.0 / 2.0))
		{
			err = "Maturity date is smaller than 6m from today, we cannot calibrate";
			goto FREE_RETURN;
		}


		/* 1) Initialisation of the Underlyings */

		/* allocation */
		vol_curve = dmatrix(0, 2, 0, 0);
		tau_curve = dmatrix(0, 1, 0, 0);

		if (!vol_curve || !tau_curve )
		{
			err = "Memory allocation failure (1) in Get_Calib_LGM_Underlying";
			goto FREE_RETURN;
		}


		vol_curve[0][0] = tau_curve[0][0] = today + 100.0;
		vol_curve[1][0] = 0.01;
		tau_curve[1][0] = tau;

		err = SrtInitIRUnd	(und_name, yc, "LGM1F",
							1, 3, vol_curve,
							1, 2, tau_curve,
							0, 0, 0, 0, 0, 0,
							0.0,0,0,NULL);
		if (err)
			goto FREE_RETURN;

		

		/* 2) Construct the calibration schedule */

		call_end = spot_date;
		nInput = 0;

		/* look for the first maturity in number of year wich is above maturity_date */
		while (call_end < maturity_date_adj)
		{
			nInput += 1;
			call_end = add_unit (spot_date, nInput, SRT_YEAR, MODIFIED_SUCCEEDING);
		}

		if (nInput == 1)
		{
			call_start = add_unit (spot_date, 6, SRT_MONTH, MODIFIED_SUCCEEDING);	
			call_end = add_unit (call_start, 6, SRT_MONTH, MODIFIED_SUCCEEDING);
		}
		else
		{
			call_start = spot_date;

			/* get the schedule */
			err = swp_f_initSwapDP(	call_start, 
									call_end,
									call_freq, 
                  					basis,
									&Swap);
			if (err)
				goto FREE_RETURN;

			dateList = SwapDP_to_DateList(&Swap, NOTBROKEN);

			/*
			if (!dateList)
			{
				err = "Cannot construct schedule";
				goto FREE_RETURN;
			}
			*/

			nInput = dateList.len - 2;
		}

		
		sdp = calloc (nInput, sizeof(SwapDP));

		recPayStr = calloc(nInput, sizeof(SrtReceiverType));
		refRateStr = calloc(nInput, sizeof(char *));
		type = calloc(nInput, sizeof(int));
		strike = calloc(nInput, sizeof(double));
		bondStrike = calloc(nInput, sizeof(double));
		undPrice = calloc(nInput, sizeof(double));

		if (!sdp || !recPayStr || !type || !strike || !bondStrike || !undPrice) 
		{
			err = "Memory allocation failure (2) in Get_Calib_LGM_Underlying";
			goto FREE_RETURN;
		}

		for (i=0; i<nInput; i++)
		{
			if (nInput > 1)
			{
				type[i] = SWAPTION;

				err = swp_f_initSwapDP(dateList.date[i+1], call_end, cpd, basis, &(sdp[i]));				
				
				if (err)
					goto FREE_RETURN;
				
				err = swp_f_ForwardRate(dateList.date[i+1], call_end, cpd, basis, yc, ref_rate, &(strike[i]));
				
				if (err)
					goto FREE_RETURN;
				
				err = swp_f_ForwardRate(dateList.date[i+1], call_end, cpd, basis, yc, "CASH", &Spread);
				
				if (err)
					goto FREE_RETURN;

				Spread -= strike[i];
				Spread *= -1;

				if (err)
					goto FREE_RETURN;

				err = pfGetVol(dateList.date[i+1], call_end, strike[i], strike[i], Spread, &vol);
				if (err)
					goto FREE_RETURN;

				err = swp_f_Swaption(dateList.date[i+1], call_end, cpd, basis,
										vol,
										strike[i],
										"PAY",
										ref_rate,
										yc,
										"PREMIUM",
										log_nor,
										&(undPrice[i]));
				if (err)
					goto FREE_RETURN;
			}
			else
			{
				type[i] = CAPFLOOR;

				err = swp_f_initSwapDP(call_start, call_end, cpd, basis, &(sdp[i]));
				if (err)
					goto FREE_RETURN;
				err = swp_f_ForwardRate(call_start, call_end, cpd, basis, yc, ref_rate, &(strike[i]));
				if (err)
					goto FREE_RETURN;
				err = pfGetVol(call_start, call_end, strike[i], strike[i], Spread, &vol);
				if (err)
					goto FREE_RETURN;
				err = swp_f_CapFloor(
									call_start,
									call_end,
									strike[i],
									pfGetVol,
									"CAP",									
									ref_rate,
									yc,
									"PREMIUM",
									log_nor,
									&(undPrice[i]));
				if (err)
					goto FREE_RETURN;				
			}
			
			recPayStr[i] = SRT_PAYER;
			refRateStr[i] = ref_rate;
		}

		undptr = lookup_und(und_name);

		err =  srt_f_bootstrap_calibrate(
										undptr,
										grfnparam,
										sdp,
										strike,
										bondStrike,
										type,
										recPayStr,
										refRateStr,
										undPrice,
										undPrice,
										tau,
										0,
										0,
										0,
										&nInput
										);

		if (err)
			goto FREE_RETURN;

FREE_RETURN:

		if (vol_curve)
		{
			free_dmatrix(vol_curve, 0, 2, 0, 0);
		}

		if (tau_curve)
		{
			free_dmatrix(tau_curve, 0, 1, 0, 0);
		}
				
		if (nInput > 1 && dateList.len > 0)
			swp_f_free_in_DateList(dateList);		

		if (type)
		{
			free(type);
		}

		if (recPayStr)
		{
			free(recPayStr);
		}

		if (refRateStr)
		{
			free(refRateStr);
		}
	
		if (undPrice)
		{
			free (undPrice);
		}

		if (strike)
		{
			free (strike);
		}

		if (bondStrike)
		{
			free (bondStrike);
		}	

		if (sdp)
		{
			free (sdp);
		}
	
	return err;
	
}

Err Fx3DDefImpVol2(		/*	Total number of options */
						int		num_opt,

						/*	Number of options to calibrate to */
						int		num_calib,

						/*	Implied volatility */
						double	*opt_exe,
						double	*opt_mat,
						double	*ivol,

						/*	Market parameters */		
						long	today,
						long	spot_date_dom,
						long	spot_date_for,
						long	maturity_date,
						char    *dom_yc,
						char	*vol_dom,
						char	*ref_rate_dom,
						char	*basis_dom,
						char	*cpd_dom,
						char	*ccy_dom,					
						char	*for_yc,
						char	*vol_for,
						char	*ref_rate_for, 						
						char	*basis_for,
						char	*cpd_for,
						char	*ccy_for,

						/* Model parameters */						
						double	lam_dom,
						double	lam_for,
						double	rho_dom_for,
						double	rho_dom_fx,
						double	rho_for_fx,

						/* Extrapolation parameters */
						double	vol_lim,
						double	vol_lam,

						/* Calibration parameters */
						char	*call_freq,
						double	dom_vol_shift,
						double	for_vol_shift,

						/* Vol Function */
						Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
														char	*vol_curve_name,	
														double	start_date, 
														double	end_date,
														double	cash_strike,
														int		zero,
														char	*ref_rate_name,
														double	*vol,
														double	*power),

						/* Answer */
						double	*vol
					)
{
SrtGrfnParam    grfnparam;
double			*fx_vol_curve = NULL;
double			opt_maturity, exe_maturity;
long			opt_date_lim, maturity_calib;
long			exercise_date;

long			sigma_n_dom, sigma_n_for;
long			nb_merge_dates, spot_date;
double			*sigma_date_dom = NULL,
				*sigma_dom		= NULL,								
				*sigma_date_for	= NULL,
				*sigma_for		= NULL,
				*merge_dates	= NULL,
				*sig_dom		= NULL,
				*sig_for		= NULL;

long			*ex_date		= NULL;

SrtCompounding	call_comp, dom_comp, for_comp;

int				i, call_nb, nb_exe, call_nb_max;
long			cur_date, cur_date_adj;

Err				err = NULL;

	
	err =  srt_f_set_default_GrfnParams(&grfnparam);

	if (err)
	{
		goto FREE_RETURN;
	}

	exercise_date = add_unit (maturity_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING); 
	
	opt_maturity = (maturity_date - today) / 365.0;
	exe_maturity = (exercise_date - today) / 365.0;

	spot_date = max(spot_date_dom, spot_date_for);

	err = interp_compounding(call_freq, &call_comp);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = interp_compounding(cpd_dom, &dom_comp);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = interp_compounding(cpd_for, &for_comp);
	if (err)
	{
		goto FREE_RETURN;
	}

	call_nb_max = 12 / (min(min(call_comp, dom_comp), for_comp));

	opt_date_lim = add_unit (spot_date, 2 * call_nb_max, SRT_MONTH, MODIFIED_SUCCEEDING);

	if (maturity_date < opt_date_lim)
	{
		maturity_calib = opt_date_lim;
	}
	else
	{
		maturity_calib = maturity_date;
	}

	/* calibration dates */

	call_nb = 12 / call_comp;

	cur_date = maturity_calib;
	cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

	nb_exe = 0;

	while (cur_date_adj > spot_date)
	{
		cur_date = add_unit(cur_date, -call_nb, SRT_MONTH, NO_BUSDAY_CONVENTION);
		cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

		nb_exe ++;
	}

	nb_exe --;

	ex_date = calloc(nb_exe, sizeof(long));

	if (!ex_date)
	{
		err = "Memory Allocation faillure in FX3DVolDef2";
		goto FREE_RETURN;
	}


	cur_date = maturity_calib;

	for (i=nb_exe-1; i>=0; i--)
	{
		cur_date = add_unit(cur_date, -call_nb, SRT_MONTH, NO_BUSDAY_CONVENTION);
		cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

		ex_date[i] = add_unit(cur_date_adj, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
	}


	err = old_cpd_calib_diagonal_wrapper(
								dom_yc,
								vol_dom,
								ref_rate_dom,
								get_cash_vol,
								dom_vol_shift,
								1,
								nb_exe,
								ex_date,
								maturity_calib,
								NULL,
								NULL,
								0,
								1.0,
								1.0,
								cpd_dom,
								basis_dom,
								1,
								0,
								0,
								0,
								0,
								0,
								0,
								NULL,
								&lam_dom,
								1,
								0.0,
								0.0,
								0.0,
								&sigma_n_dom,
								&sigma_date_dom,
								&sigma_dom,
								NULL);
										
	if (err)
	{
		goto FREE_RETURN;
	}

	err = old_cpd_calib_diagonal_wrapper(
								for_yc,
								vol_for,
								ref_rate_for,
								get_cash_vol,
								for_vol_shift,
								1,
								nb_exe,
								ex_date,
								maturity_calib,
								NULL,
								NULL,
								0,
								1.0,
								1.0,
								cpd_for,
								basis_for,
								1,
								0,
								0,
								0,
								0,
								0,
								0,
								NULL,
								&lam_for,
								1,
								0.0,
								0.0,
								0.0,
								&sigma_n_for,
								&sigma_date_for,
								&sigma_for,
								NULL);
										
	if (err)
	{
		goto FREE_RETURN;
	}

	err = merge_rates_ts (
							sigma_date_dom, 
							sigma_dom,
							sigma_n_dom,
							sigma_date_for, 
							sigma_for,
							sigma_n_for,
							&merge_dates,
							&sig_dom,
							&sig_for,
							&nb_merge_dates);

	if (err)
	{
		goto FREE_RETURN;
	}

	err = Fx3DtsCalibration(
							opt_exe,
							opt_mat,
							ivol,
							num_calib,
							merge_dates,
							nb_merge_dates,
							sig_dom,
							lam_dom,
							sig_for,
							lam_for,
							rho_dom_for,
							rho_dom_fx,
							rho_for_fx,
							&fx_vol_curve);

	
	if (err)
	{
		goto FREE_RETURN;
	}

	err = Fx3DtsImpliedVolExtrapol(	opt_maturity,
									0,
									exe_maturity,
									merge_dates,
									nb_merge_dates,
									sig_dom,
									lam_dom,
									sig_for,
									lam_for,
									opt_exe,
									fx_vol_curve,
									num_calib, 
									vol_lim,
									vol_lam,
									rho_dom_for,
									rho_dom_fx,
									rho_for_fx,
									vol);

FREE_RETURN:

	if (fx_vol_curve)
	{
		free (fx_vol_curve);
	}
	
	if (sigma_date_dom)
	{
		free (sigma_date_dom);
	}

	if (sigma_dom)
	{
		free (sigma_dom);
	}

	if (sigma_date_for)
	{
		free (sigma_date_for);
	}

	if (sigma_for)
	{
		free (sigma_for);
	}

	if (sig_dom)
	{
		free (sig_dom);
	}

	if (sig_for)
	{
		free (sig_for);
	}

	if (merge_dates)
	{
		free (merge_dates);
	}

	if (ex_date)
	{
		free (ex_date);
	}	

	return err;
}

Err Fx3DDefImpVol2_corr(		/*	Total number of options */
						int		num_opt,

						/*	Number of options to calibrate to */
						int		num_calib,

						/*	Implied volatility */
						double	*opt_exe,
						double	*opt_mat,
						double	*ivol,

						/*	Market parameters */		
						long	today,
						long	spot_date_dom,
						long	spot_date_for,
						long	maturity_date,
						char    *dom_yc,
						char	*vol_dom,
						char	*ref_rate_dom,
						char	*basis_dom,
						char	*cpd_dom,
						char	*ccy_dom,					
						char	*for_yc,
						char	*vol_for,
						char	*ref_rate_for, 						
						char	*basis_for,
						char	*cpd_for,
						char	*ccy_for,

						/* Model parameters */						
						double	lam_dom,
						double	lam_for,
						double	*rho_times,
						double	*rho_dom_for_ts,
						double	*rho_dom_fx_ts,
						double	*rho_for_fx_ts,
						long	nb_rho,

						/* Extrapolation parameters */
						double	vol_lim,
						double	vol_lam,

						/* Calibration parameters */
						char	*call_freq,
						double	dom_vol_shift,
						double	for_vol_shift,

						/* Vol Function */
						Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
														char	*vol_curve_name,	
														double	start_date, 
														double	end_date,
														double	cash_strike,
														int		zero,
														char	*ref_rate_name,
														double	*vol,
														double	*power),

						/* Answer */
						double	*vol
					)
{
SrtGrfnParam    grfnparam;
double			*fx_vol_curve = NULL;
double			opt_maturity, exe_maturity;
long			opt_date_lim, maturity_calib;
long			exercise_date;

long			sigma_n_dom, sigma_n_for;
long			nb_merge_dates, spot_date;
double			*sigma_date_dom = NULL,
				*sigma_dom		= NULL,								
				*sigma_date_for	= NULL,
				*sigma_for		= NULL,
				*merge_dates	= NULL,
				*sig_dom		= NULL,
				*sig_for		= NULL;

long			*ex_date		= NULL;

SrtCompounding	call_comp, dom_comp, for_comp;

int				i, call_nb, nb_exe, call_nb_max;
long			cur_date, cur_date_adj;

Err				err = NULL;

	
	err =  srt_f_set_default_GrfnParams(&grfnparam);

	if (err)
	{
		goto FREE_RETURN;
	}

	exercise_date = add_unit (maturity_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING); 
	
	opt_maturity = (maturity_date - today) / 365.0;
	exe_maturity = (exercise_date - today) / 365.0;

	spot_date = max(spot_date_dom, spot_date_for);

	err = interp_compounding(call_freq, &call_comp);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = interp_compounding(cpd_dom, &dom_comp);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = interp_compounding(cpd_for, &for_comp);
	if (err)
	{
		goto FREE_RETURN;
	}

	call_nb_max = 12 / (min(min(call_comp, dom_comp), for_comp));

	opt_date_lim = add_unit (spot_date, 2 * call_nb_max, SRT_MONTH, MODIFIED_SUCCEEDING);

	if (maturity_date < opt_date_lim)
	{
		maturity_calib = opt_date_lim;
	}
	else
	{
		maturity_calib = maturity_date;
	}

	/* calibration dates */

	call_nb = 12 / call_comp;

	cur_date = maturity_calib;
	cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

	nb_exe = 0;

	while (cur_date_adj > spot_date)
	{
		cur_date = add_unit(cur_date, -call_nb, SRT_MONTH, NO_BUSDAY_CONVENTION);
		cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

		nb_exe ++;
	}

	nb_exe --;

	ex_date = calloc(nb_exe, sizeof(long));

	if (!ex_date)
	{
		err = "Memory Allocation faillure in FX3DVolDef2";
		goto FREE_RETURN;
	}


	cur_date = maturity_calib;

	for (i=nb_exe-1; i>=0; i--)
	{
		cur_date = add_unit(cur_date, -call_nb, SRT_MONTH, NO_BUSDAY_CONVENTION);
		cur_date_adj = add_unit(cur_date, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

		ex_date[i] = add_unit(cur_date_adj, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
	}


	err = old_cpd_calib_diagonal_wrapper(
								dom_yc,
								vol_dom,
								ref_rate_dom,
								get_cash_vol,
								dom_vol_shift,
								1,
								nb_exe,
								ex_date,
								maturity_calib,
								NULL,
								NULL,
								0,
								1.0,
								1.0,
								cpd_dom,
								basis_dom,
								1,
								0,
								0,
								0,
								0,
								0,
								0,
								NULL,
								&lam_dom,
								1,
								0.0,
								0.0,
								0.0,
								&sigma_n_dom,
								&sigma_date_dom,
								&sigma_dom,
								NULL);
										
	if (err)
	{
		goto FREE_RETURN;
	}

	err = old_cpd_calib_diagonal_wrapper(
								for_yc,
								vol_for,
								ref_rate_for,
								get_cash_vol,
								for_vol_shift,
								1,
								nb_exe,
								ex_date,
								maturity_calib,
								NULL,
								NULL,
								0,
								1.0,
								1.0,
								cpd_for,
								basis_for,
								1,
								0,
								0,
								0,
								0,
								0,
								0,
								NULL,
								&lam_for,
								1,
								0.0,
								0.0,
								0.0,
								&sigma_n_for,
								&sigma_date_for,
								&sigma_for,
								NULL);
										
	if (err)
	{
		goto FREE_RETURN;
	}

	err = merge_rates_ts (
							sigma_date_dom, 
							sigma_dom,
							sigma_n_dom,
							sigma_date_for, 
							sigma_for,
							sigma_n_for,
							&merge_dates,
							&sig_dom,
							&sig_for,
							&nb_merge_dates);

	if (err)
	{
		goto FREE_RETURN;
	}

	err = Fx3DtsCalibration_corr(
							opt_exe,
							opt_mat,
							ivol,
							num_calib,
							merge_dates,
							nb_merge_dates,
							sig_dom,
							lam_dom,
							sig_for,
							lam_for,
							rho_times,
							rho_dom_for_ts,
							rho_dom_fx_ts,
							rho_for_fx_ts,
							nb_rho,
							&fx_vol_curve);

	
	if (err)
	{
		goto FREE_RETURN;
	}

	err = Fx3DtsImpliedVolExtrapol_corr(	opt_maturity,
											0,
											exe_maturity,
											merge_dates,
											nb_merge_dates,
											sig_dom, lam_dom,
											sig_for, lam_for,
											opt_exe,
											fx_vol_curve,
											num_calib, 
											vol_lim,
											vol_lam,
											rho_times,
											rho_dom_for_ts,
											rho_dom_fx_ts,
											rho_for_fx_ts,
											nb_rho,
											vol);

FREE_RETURN:

	if (fx_vol_curve)
	{
		free (fx_vol_curve);
	}
	
	if (sigma_date_dom)
	{
		free (sigma_date_dom);
	}

	if (sigma_dom)
	{
		free (sigma_dom);
	}

	if (sigma_date_for)
	{
		free (sigma_date_for);
	}

	if (sigma_for)
	{
		free (sigma_for);
	}

	if (sig_dom)
	{
		free (sig_dom);
	}

	if (sig_for)
	{
		free (sig_for);
	}

	if (merge_dates)
	{
		free (merge_dates);
	}

	if (ex_date)
	{
		free (ex_date);
	}	

	return err;
}

Err Fx3DDefImpVolBeta(		/*	Total number of options */
						int		num_opt,

						/*	Number of options to calibrate to */
						int		num_calib,

						/*	Implied volatility */
						double	*opt_exe,
						double	*opt_mat,
						double	*ivol,

						/*	Market parameters */		
						long	today,
						long	spot_date_dom,
						long	spot_date_for,
						long	maturity_date,
						double	spot_fx,
						char    *dom_yc,
						char	*vol_dom,
						char	*ref_rate_dom,
						char	*basis_dom,
						char	*cpd_dom,
						char	*ccy_dom,					
						char	*for_yc,
						char	*vol_for,
						char	*ref_rate_for, 						
						char	*basis_for,
						char	*cpd_for,
						char	*ccy_for,

						/* Model parameters */						
						double	lam_dom,
						double	lam_for,
						double	rho_dom_for,
						double	rho_dom_fx,
						double	rho_for_fx,
						double	beta,

						/* Extrapolation parameters */
						double	vol_lim,
						double	vol_lam,

						/* Calibration parameters */
						char	*call_freq,
						double	dom_vol_shift,
						double	for_vol_shift,

						long	nbSteps,
						double	disc_dt,
						double	fx_dt,
						long	nbIterMax,

						/* Vol Function */
						Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
														char	*vol_curve_name,	
														double	start_date, 
														double	end_date,
														double	cash_strike,
														int		zero,
														char	*ref_rate_name,
														double	*vol,
														double	*power),

						/* Answer */
						double	*vol
					)
{
SrtGrfnParam    grfnparam;
double			*fx_vol_curve = NULL;
double			opt_maturity, exe_maturity;
long			opt_date_lim, maturity_calib;
long			exercise_date;
double			fwd, price, df;

long			sigma_n_dom, sigma_n_for;
long			i, nb_merge_dates, nb_fx_vol;
double			*sigma_date_dom = NULL,
				*sigma_dom		= NULL,								
				*sigma_date_for	= NULL,
				*sigma_for		= NULL,
				*merge_dates	= NULL,
				*sig_dom		= NULL,
				*sig_for		= NULL,
				*beta_tab		= NULL,
				*new_vol_mat	= NULL;


Err				err = NULL;

	
	err =  srt_f_set_default_GrfnParams(&grfnparam);

	if (err)
		goto FREE_RETURN;

	exercise_date = add_unit (maturity_date, -2, SRT_BDAY, MODIFIED_SUCCEEDING); 
	
	opt_maturity = (maturity_date - today) / 365.0;
	exe_maturity = (exercise_date - today) / 365.0;

	opt_date_lim = add_unit (spot_date_dom, 6, SRT_MONTH, MODIFIED_SUCCEEDING) + 1;

	if (maturity_date < opt_date_lim)
	{
		maturity_calib = opt_date_lim;
	}
	else
	{
		maturity_calib = maturity_date;
	}

	spot_fx *= swp_f_df (today, spot_date_dom, dom_yc) / swp_f_df (today, spot_date_for, for_yc);
	df = swp_f_df(today, exercise_date, dom_yc);
	fwd = spot_fx * swp_f_df(today, exercise_date, for_yc) / df;


	err = old_cpd_calib_diagonal_wrapper(
								dom_yc,
								vol_dom,
								ref_rate_dom,
								get_cash_vol,
								dom_vol_shift,
								1,
								0,
								NULL,
								maturity_calib,
								NULL,
								NULL,
								0,
								1.0,
								1.0,
								cpd_dom,
								basis_dom,
								1,
								0,
								0,
								0,
								0,
								0,
								0,
								NULL,
								&lam_dom,
								1,
								0.0,
								0.0,
								0.0,
								&sigma_n_dom,
								&sigma_date_dom,
								&sigma_dom,
								NULL);
										
	if (err)
		goto FREE_RETURN;
	
	opt_date_lim = add_unit (spot_date_for, 6, SRT_MONTH, MODIFIED_SUCCEEDING) + 1;

	if (maturity_date < opt_date_lim)
	{
		maturity_calib = opt_date_lim;
	}
	else
	{
		maturity_calib = maturity_date;
	}

	err = old_cpd_calib_diagonal_wrapper(
								for_yc,
								vol_for,
								ref_rate_for,
								get_cash_vol,
								for_vol_shift,
								1,
								0,
								NULL,
								maturity_calib,
								NULL,
								NULL,
								0,
								1.0,
								1.0,
								cpd_for,
								basis_for,
								1,
								0,
								0,
								0,
								0,
								0,
								0,
								NULL,
								&lam_for,
								1,
								0.0,
								0.0,
								0.0,
								&sigma_n_for,
								&sigma_date_for,
								&sigma_for,
								NULL);
										
	if (err)
		goto FREE_RETURN;

	err = merge_rates_ts (
							sigma_date_dom, 
							sigma_dom,
							sigma_n_dom,
							sigma_date_for, 
							sigma_for,
							sigma_n_for,
							&merge_dates,
							&sig_dom,
							&sig_for,
							&nb_merge_dates);

	if (err)
		goto FREE_RETURN;

	beta_tab = dvector(0, num_calib - 1);

	for (i=0; i<num_calib; i++)
	{
		beta_tab[i] = beta;
	}

	err = Fx3DBetatsCalibration(
								today,
								opt_exe,
								opt_mat,
								ivol,
								num_calib,
								merge_dates,
								nb_merge_dates,
								sig_dom,
								lam_dom,
								sig_for,
								lam_for,
								beta_tab,
								spot_fx,
								rho_dom_for,
								rho_dom_fx,
								rho_for_fx,
								dom_yc,
								for_yc,
								&fx_vol_curve,
								disc_dt,
								fx_dt,
								nbIterMax);

	
	if (err)
		goto FREE_RETURN;

	/* Now we need to fill the new fx beta vol term structure */
	/* we add a point every year */


	nb_fx_vol = num_calib + (int) ((exe_maturity - opt_exe[num_calib - 1]) / 0.5) + 1;

	if (fabs(opt_exe[num_calib - 1] + (nb_fx_vol - num_calib -1) * 0.5 - exe_maturity) < 0.01)
	{
		nb_fx_vol -= 1;
	}


	fx_vol_curve = realloc(fx_vol_curve, nb_fx_vol * sizeof(double));	
	new_vol_mat = dvector(0, nb_fx_vol - 1);

	if (!new_vol_mat || !fx_vol_curve)
		goto FREE_RETURN;

	for (i=0; i<num_calib; i++)
	{
		new_vol_mat[i] = opt_exe[i];
	}

	for (i=num_calib; i<nb_fx_vol-1; i++)
	{
		new_vol_mat[i] = new_vol_mat[i-1] + 0.5;

		fx_vol_curve[i] = vol_lim + (fx_vol_curve[num_calib - 1] - vol_lim) 
							* exp(-vol_lam * (new_vol_mat[i] - new_vol_mat[num_calib - 1]));
	}

	new_vol_mat[nb_fx_vol - 1] = exe_maturity;
	fx_vol_curve[nb_fx_vol - 1] = vol_lim + (fx_vol_curve[num_calib - 1] - vol_lim) 
							* exp(-vol_lam * (exe_maturity - new_vol_mat[num_calib - 1]));


	err = Fx3DBetatsTreeFxOptions(	
								today,
								exercise_date,
								&fwd,
								1,
								merge_dates,
								nb_merge_dates,
								sig_dom,
								lam_dom,
								sig_for,
								lam_for,
								new_vol_mat,
								nb_fx_vol,
								fx_vol_curve,
								0.0,
								beta,
								spot_fx,
								rho_dom_for,
								rho_dom_fx,
								rho_for_fx,
								dom_yc,
								for_yc,												
								&price,
								nbSteps
								);
	if (err)
		goto FREE_RETURN;

	err = srt_f_optimpvol(		
								price,
								fwd,
								fwd,
								exe_maturity,								
								df,
								SRT_CALL,
								SRT_LOGNORMAL,
								vol);	

FREE_RETURN:

	if (fx_vol_curve)
		free (fx_vol_curve);
	
	if (sigma_date_dom)
	{
		free (sigma_date_dom);
	}

	if (sigma_dom)
	{
		free (sigma_dom);
	}

	if (sigma_date_for)
	{
		free (sigma_date_for);
	}

	if (sigma_for)
	{
		free (sigma_for);
	}

	if (sig_dom)
	{
		free (sig_dom);
	}

	if (sig_for)
	{
		free (sig_for);
	}

	if (merge_dates)
	{
		free (merge_dates);
	}

	if (beta_tab)
	{
		free_dvector(beta_tab, 0, nb_merge_dates);
	}

	if (new_vol_mat)
	{
		free_dvector(new_vol_mat, 0, nb_fx_vol-1);
	}

	return err;
}
