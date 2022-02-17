/* ==========================================================================
   FILE_NAME:	FxOMModel.c

   PURPOSE:		Modelling of the spot FX vol by taking in consideration a stochastic 
				spread between the two interest rates - domestic and foreign
   ========================================================================== */

#include "utallhdr.h"
#include "math.h"

static double L_Func(double lambda, double t, double T)
{
double result;

	result = (1. - exp(-lambda*(T-t))) / lambda;

	return (result);
}

static double K_func(double lambda, double T1, double T2, double T)
{
double result;

	result = (exp(-lambda*(T-T2)) - exp(-lambda*(T-T1))) / lambda; 

	return (result);
}

static double CalcVar(double vol, double spreadvol, double rho, double lambda,
					  double T1, double T2, double T)
{
double var;

	var = vol*vol*(T2-T1);
	var += spreadvol*spreadvol*
				(T2-T1 
				+ K_func(2.*lambda, T1, T2, T)
				- 2. * K_func(lambda, T1, T2, T) )
			/ (lambda*lambda);

	var += 2.*rho*vol*spreadvol*
				(T2-T1 - K_func(lambda, T1, T2, T))
			/ lambda;

	return var;
}


/* Get the implied volatility of a FXoption after having calibrated the model */
Err GetImpVol(	double optmat,
				double * mats,
				double * fxvols,
				long nbropts,
				double spreadvol,
				double rho,
				double lambda,
				double * impvol)
{
int i;

	/* we do not allow extrapolation so far */
	/* (will have to modify the program to allow that) */
	if (optmat > mats[nbropts]) 
		return ("Cannot compute the implied volatility - no extrapolation allowed");

	(*impvol) = 0.;
	for (i=1 ; (i <= nbropts) && (mats[i] < optmat) ; i++)
	{

		(*impvol) += CalcVar(fxvols[i], spreadvol, rho, lambda, 
							mats[i-1], mats[i], optmat);
	}

	/* add the last part of the cumulative variance */
	/* need the previous index */
	(*impvol) += CalcVar(fxvols[i], spreadvol, rho, lambda, 
							mats[i-1], optmat, optmat);

	(*impvol) = sqrt((*impvol) / optmat);

	return NULL;
}



Err CalibFxOMModel(	double * impvols,
					double * optmats,
					long nbropts,
					double spreadvol,
					double rho,
					double lambda,
					double ** fxvols)
{
int idxopt, idxcumvar;
double cumvar;
double a,b,c, delta;


	/* memory allocation */
	(*fxvols) = (double *) srt_malloc((nbropts+1)*sizeof(double));
	if (!(*fxvols))
	{
		return ("Memory allocation error");
	}

	/* loop on the number of options */
	for (idxopt = 1 ; idxopt <= nbropts ; idxopt++)
	{
		/* Get the cumulated variance till the previous maturity */
		cumvar = 0.;
		for (idxcumvar = 1 ; idxcumvar < idxopt ; idxcumvar++)
		{
			cumvar += CalcVar((*fxvols)[idxcumvar], spreadvol, rho, lambda, 
								optmats[idxcumvar-1], optmats[idxcumvar], optmats[idxopt]);
		}

		/* find the spot fx volatility fxvols[idxopt] to match */
		/* the implied volatility of the idxopt option */
		a = optmats[idxopt] - optmats[idxopt-1];
		b = 2.*rho*spreadvol*(a-K_func(lambda,optmats[idxopt-1],optmats[idxopt],optmats[idxopt]))/lambda;
		
		c = spreadvol * spreadvol * (a 
									+ K_func(2.*lambda, optmats[idxopt-1],optmats[idxopt],optmats[idxopt])
									- 2. * K_func(lambda, optmats[idxopt-1],optmats[idxopt],optmats[idxopt]) );
		c /= (lambda*lambda);
		c += cumvar - impvols[idxopt]*impvols[idxopt]*optmats[idxopt];

		/* just solve the second order equation */
		delta = b*b - 4*a*c;
		if (delta <= 0)
		{
			free((*fxvols));
			return serror("Cannot find solution to match the %d option vol", idxopt);
		}
		
		/* if there is two positive solutions we are taking the smallest one */
		/* it would be easier to match the next one */
		if ((b + sqrt(delta)) < 0)
		{
			(*fxvols)[idxopt] = (-b - sqrt(delta)) / (2.*a);
		}
		else
		{
			(*fxvols)[idxopt] = (-b + sqrt(delta)) / (2.*a);
		}
	}

	return NULL;
}

/* Main Function */
Err ImpFxOpt(double optmat,
			 double * calibimpvols, double * caliboptmats, long nbrcalibopts, 
			 double spreadvol, double rho, double lambda,
			 double * impvol)
{
Err err = NULL;
double * calibfxvols = NULL; /* calibrated fx spot vols */
	
	/* Calibration of the model */
	err = CalibFxOMModel(calibimpvols,
						 caliboptmats,
						 nbrcalibopts,
						 spreadvol,	
						 rho,
						 lambda,					
						 &calibfxvols);

	if (err)
		return err;

	/* Get the implied vol */
	err = GetImpVol(optmat,
					caliboptmats,
					calibfxvols,
					nbrcalibopts,
					spreadvol,
					rho,
					lambda,
					impvol);
	
	return err;
}