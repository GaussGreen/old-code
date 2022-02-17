

/******************************************************************************/
/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_optblkvoltobetavol(...)                                   
*                                                                       
* PURPOSE      	: A Quick vol transformation to go from a BS vol to a BS beta vol   
*                                                                      
* DESCRIPTION  	: The volatility transformation is based on a small noise expansion
                  of the density for both BS model and its equivalent BS beta                       
				  (see Pat Hagan's paper on the subject)
*                                                                   
*		                                    
* PARAMETERS                                                           
*	INPUT	    : fwd_price	    - forward underlying price    	 
*              	: strike      	- strike price		      	  
*              	: bsvol         - annual bs volatility   	      	  
*              	: mat         	- initial time, in years  	      	   	      
*				: beta          - dF = a*F^{beta}dW
*
*                                                      
* RETURNS      	: betavol       - equivalent BS Beta vol	      	 
*                                                                       
*******************************************************************************/

/*******************************************************************************
**                      Include Files
*******************************************************************************/


#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

static	double	static_forward;							
static	double	static_strike;
static	double	static_mat;
static	double	static_beta;
static  double  static_bsvol;
static  double  static_alpha;
static  double  static_rho;

/* ------------------------------------------------------------------------------ */

static Err num_f_optblkvoltobetavol(
		double            forward,							
		double            strike, 
		double            bsvol, 
		double            mat,  
		double            beta,
		double           *betavol)
{
double   middle;
double   distance;
double   space_correction;
double   time_correction;
double   variance;
Err      err                = NULL;

/* A few security checks */
	if ((mat <= 0.0 ) || (bsvol == 0.0))
	{
		*betavol = 0.0;
		return NULL;
	}

	if ((beta > 1.0 ) || (beta < 0.0 ))
		return serror("Beta has to be between 0.0 and 1.0");

/* The point around which the expansion is made : 0.5*(forward + strike) */
	middle = 0.5 * ( forward + strike );

/* The distance between the forward and the strike : 0.5 * ( forward- strike ) */
	distance = forward - strike;

/* The variance of the distribution */
	variance = bsvol * bsvol * mat;
	
/* The leading order term for corection in space (diffusion term) */
	space_correction = (1.0 - beta)*(2.0 + beta)/24.0*(distance/middle)*(distance/middle);

/* The leading order term for corection due to time */
	time_correction = (1.0 - beta)*(1.0 - beta)/24.0 * variance;
	
/* The equivalent BlackScholes BEta vol (constant elasticity) */
	*betavol = bsvol * pow(middle, 1.0 - beta ) 
			* ( 1.0 - (space_correction + time_correction) / (1.0 + space_correction + 3.0 * time_correction));

/* Return a success message */
	return NULL;
}


/* ------------------------------------------------------------------------------- */

static Err rtsafe_func(double x, double *fx, double *dfx)
{
  double dx = 1.0e-6;
  double vol_tmp, vol_tmp2;
  Err  err;

	err = srt_f_optbetavoltoblkvol(
		static_forward,							
		static_strike, 
		x, 
		static_mat,  
		static_beta,
		&vol_tmp);

	if(err) return err;

	err = srt_f_optbetavoltoblkvol(
		static_forward,							
		static_strike, 
		x + dx, 
		static_mat,  
		static_beta,
		&vol_tmp2);
	
	if(err) return err;

	 *fx = vol_tmp - static_bsvol;
	 *dfx = (vol_tmp2 - vol_tmp) / dx;

	 return err;
 
}

/* ---------------------------------------------------------------------------------- */
Err srt_f_optblkvoltobetavol(
		double            forward,							
		double            strike, 
		double            bsvol, 
		double            mat,  
		double            beta,
		double           *betavol)
{
double   guess_vol_down;
double   guess_vol_up;
double   guess_vol;
double   direc = 1;

double   vol_log;
double   vol_log_up;
double   vol_log_down;

double   stop_eps = 1e-8;
double   vol_bump = 0.04;

int		 setup = 0;
int      niter =0;
Err      err = NULL;

/* A few security checks */
	if ((mat <= 0.0 ) || (bsvol == 0.0))
	{
		*betavol = 0.0;
		return NULL;
	}

	if ((beta > 1.0 ) || (beta < 0.0 ))
		return serror("Beta has to be between 0.0 and 1.0");

	static_forward = forward;							
	static_strike = strike;
	static_mat = mat;
	static_beta = beta;
	static_bsvol = bsvol;

	err = num_f_optblkvoltobetavol(
		forward,							
		strike, 
		bsvol, 
		mat,  
		beta,
		&guess_vol);
	if (err) return err;

	err = srt_f_optbetavoltoblkvol(
		forward,							
		strike, 
		guess_vol, 
		mat,  
		beta,
		&vol_log);

	if (fabs(bsvol - vol_log) < stop_eps)
	{
		*betavol = guess_vol;
		return NULL;
	}

	if (vol_log > bsvol) direc = -1;

	while ((setup == 0) && (niter <=5))
	{
		if ((bsvol - vol_bump <= 0.0) && (direc <0.0)) vol_bump = 0.9* bsvol; 

		err = num_f_optblkvoltobetavol(
			forward,							
			strike, 
			bsvol + direc * vol_bump, 
			mat,  
			beta,
			&guess_vol_down);
		if (err) return err;

		err = srt_f_optbetavoltoblkvol(
			forward,							
			strike, 
			guess_vol_down, 
			mat,  
			beta,
			&vol_log_down);

		if ((vol_log_down - bsvol)*(vol_log - bsvol) < 0.0)
		{
			if (direc < 0 )
			{
				vol_log_up = vol_log;
				guess_vol_up = guess_vol;
				setup = 1;
			}
			else
			{
				vol_log_up = vol_log_down;
				guess_vol_up = guess_vol_down;
				vol_log_down = vol_log;
				guess_vol_down = guess_vol;
				setup = 1;
			}

		}
		
		niter ++;
		vol_bump *= 1.5;
	}

	if (setup == 0)
	{
		err = "sorry, couldn't be able to initialise";
		return err;
	}
	
	err = rtsafe(rtsafe_func,guess_vol_down,guess_vol_up,stop_eps,10, betavol);
	if (err) return err;

/* Return a success message */
	return NULL;
}

static Err rtsafe_func2(double x, double *fx, double *dfx)
{
  double dx = 1.0e-6;
  double vol_tmp, vol_tmp2;
  Err  err;


	err = srt_f_optbetastochvoltoblkvol(
				static_forward,							
				static_strike, 
				x,
				static_alpha,
				static_rho,
				static_mat,  
				static_beta,
				&vol_tmp);
	
	if(err) return err;

	err = srt_f_optbetastochvoltoblkvol(
				static_forward,							
				static_strike, 
				x + dx,
				static_alpha,
				static_rho,
				static_mat,  
				static_beta,
				&vol_tmp2);

	if(err) return err;

	*fx = vol_tmp - static_bsvol;
	*dfx = (vol_tmp2 - vol_tmp) / dx;

	return err;

}

/* ---------------------------------------------------------------------------------- */
Err srt_f_optblkvolATMtobetavolStochVol(
		double            forward,							
		double            strike, 
		double            bsvol, 
		double            mat, 
		double            alpha,
		double            beta,
		double            rho,
		double           *betavol)
{
double   guess_vol_down;
double   guess_vol_up;
double   guess_vol;
double   direc = 1;

double   vol_log;
double   vol_log_up;
double   vol_log_down;

double   stop_eps = 1e-8;
double   vol_bump = 0.03;

int		 setup = 0;
int      niter =0;
Err      err = NULL;

/* A few security checks */
	if ((mat <= 0.0 ) || (bsvol == 0.0))
	{
		*betavol = 0.0;
		return NULL;
	}

	if ((beta > 1.0 ) || (beta < 0.0 ))
		return serror("Beta has to be between 0.0 and 1.0");

	static_forward = forward;							
	static_strike = strike;
	static_mat = mat;
	static_beta = beta;
	static_alpha = alpha;
	static_rho = rho;
	static_bsvol = bsvol;

	err = num_f_optblkvoltobetavol(
		forward,							
		strike, 
		bsvol, 
		mat,  
		beta,
		&guess_vol);
	if (err) return err;

	err = srt_f_optbetastochvoltoblkvol(
				forward,							
				strike, 
				guess_vol,
				alpha,
				rho,
				mat,  
				beta,
				&vol_log);

	if (fabs(bsvol - vol_log) < stop_eps)
	{
		*betavol = guess_vol;
		return NULL;
	}

	if (vol_log > bsvol) direc = -1;

	while ((setup == 0) && (niter <=10))
	{
		if ((bsvol - vol_bump <= 0.0) && (direc <0.0)) vol_bump = 0.9* bsvol; 

		err = num_f_optblkvoltobetavol(
			forward,							
			strike, 
			bsvol + direc * vol_bump, 
			mat,  
			beta,
			&guess_vol_down);
		if (err) return err;

		err = srt_f_optbetastochvoltoblkvol(
				forward,							
				strike, 
				guess_vol_down,
				alpha,
				rho,
				mat,  
				beta,
				&vol_log_down);

		if ((vol_log_down - bsvol)*(vol_log - bsvol) < 0.0)
		{
			if (direc < 0 )
			{
				vol_log_up = vol_log;
				guess_vol_up = guess_vol;
				setup = 1;
			}
			else
			{
				vol_log_up = vol_log_down;
				guess_vol_up = guess_vol_down;
				vol_log_down = vol_log;
				guess_vol_down = guess_vol;
				setup = 1;
			}

		}
		
		niter ++;
		vol_bump *= 1.5;
	}

	if (setup == 0)
	{
		err = "sorry, couldn't be able to initialise";
		return err;
	}
	
	err = rtsafe(rtsafe_func2,guess_vol_down,guess_vol_up,stop_eps,40, betavol);
	if (err) return err;

/* Return a success message */
	return NULL;
}