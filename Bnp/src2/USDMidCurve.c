/* ------------------------------------------------------------------------------------------------------------------------------------------------- */
/* USDMidCurve.c */
/* Contains functions: */
/*		*/


/* Include files */
#include "opFnctns.h"
#include "utTypes.h"
#include "math.h"
#include "USDMidCurve.h"
#include "swp_h_spread.h"
#include "swp_h_all.h"
#include "utallhdr.h"
#include "srt_h_all.h"
#include "swp_h_ccy_param.h"



char* EurodollarOptionCaller(  // For fwds and vols
						long lToday,
						char* szYieldCurveName,
						char* szVolCurve,
						char* szRefRate,
						char* (*getCashVol)(	char *szVolCurve,
													double	start_date, 
													double	end_date,
													double	cash_strike,
													int		zero,
													char	*ref_rate_name,
													double	*vol,
													double	*power),
						char* (*getDF)(char *szYieldCurve, double dStart, double dEnd, double *dDF),
						char* (*getSpread)(long start_date,long end_date,const char *szRefRate, double *dSpread),
						long lOptionExpiry,
						double dStrike,
						long lFutureStart,
						double dFuturePrice,
						int iConvexityModel,
						SrtDiffusionType enumDiffType,
						SrtCallPutType enumCallPut, // Call = long
						SrtGreekType enumGreek,
						double dDifVol,
						double dJumpAvg,
						double dJumpVol,
						double dJumpInt,
						long lIterNo, 
						double* out_dPrice,
						double* out_dVol
						)
{

// Variable declaration
	long lFutureEnd;
	SrtCurvePtr srtCurvePtr;
	SrtCcyParam *pCcyParam;
	SrtBasisCode    basis;
	SrtCompounding  compounding;
	
//  Model arguments
	double dCashStrike;
	double dCashFwd;
	double dOptionExpiry;
	double dOptionExpiryDf;

	char* err;

// Initialize return values
	*out_dPrice = 0.0;
	*out_dVol = 0.0;


// Get the yield curve pointer
	if (  !(srtCurvePtr = swp_f_lookup_curve( szYieldCurveName ))  )
		return "Could not find yield curve";

// Get the market conventions
	if ( err = swp_f_get_CcyParam_from_CcyStr( srtCurvePtr->curve_ccy, &pCcyParam )  )
		return err;

// Get the compounding from the refrate
	if ( err = swp_f_get_ref_rate_details( szRefRate, &basis, &compounding)  )
		return err;

// Calculate the end date of the fwd
	lFutureEnd = add_unit( lFutureStart, 12/compounding, SRT_MONTH, pCcyParam->swap_bus_day_conv );

	dCashStrike = (100 - dStrike);

// Fwd, option maturity and df
	dOptionExpiry = ( lOptionExpiry - lToday ) / 365.0;
	dCashFwd=(100 - dFuturePrice);

	if ( err = (*getDF)( szYieldCurveName, (double) lToday, (double) lOptionExpiry, &dOptionExpiryDf ) )
			return err;

// Calculate the price
	return EurodollarOption( 
						lToday,
						dOptionExpiry,
						dCashStrike,
						dCashFwd,
						enumDiffType,
						enumCallPut,
						enumGreek,
						dDifVol,
						dJumpAvg,
						dJumpVol,
						dJumpInt,
						lIterNo,
						dOptionExpiryDf,
						out_dPrice,
						out_dVol
						);
}	
// A MidCurve option wrapper
// ---------------------------------------------------------------------------------- //


// ---------------------------------------------------------------------------------- //

char* EurodollarOption( 
						long lToday,
						double dOptionExpiry,
						double dCashStrike,
						double dCashFwd,
						SrtDiffusionType enumDiffType,
						SrtCallPutType enumCallPut, // Call = long
						SrtGreekType enumGreek,
						double dDifVol,
						double dJumpAvg,
						double dJumpVol,
						double dJumpInt,
						long lIterNo,
						double dOptionExpiryDf,
						double* out_dPrice,
						double* out_dVol
						)
{

// Variable Declarations
	long n;
	double n_f=1;

	double dFwd_n;
	double dVol_n;
	double dPremium;
	double bs;
	
	char * err;

	
// Initialize the return quantities	
	*out_dVol = -1.0;
	*out_dPrice = -1.0;

// If the option date is negative, return 0 value
	if ( dOptionExpiry < 0.0 )
	{
		*out_dPrice = 0.0;
		*out_dVol = 0.0;
		return 0;
	}

// Reverse the call put so that it is on the yield
	if ( enumCallPut == SRT_CALL )
	{
		enumCallPut = SRT_PUT;
	}
	else
	{
		enumCallPut = SRT_CALL;
	}

// Calculate the price
	
	dPremium = 0.0;
	for (n = 0; n <= lIterNo; n++)
	{
		if (n == 0) n_f = 1;
		else n_f *= n;

		dVol_n = pow(dDifVol,2) + n * pow(dJumpVol,2)/dOptionExpiry;
		dVol_n = pow(dVol_n,0.5);

		dFwd_n= dCashFwd * exp(dJumpInt * (1-dJumpAvg) * dOptionExpiry) * pow(dJumpAvg,n);
		
		bs = srt_f_optblksch(dFwd_n, dCashStrike, dVol_n, dOptionExpiry, dOptionExpiryDf, 
										enumCallPut, enumGreek );

		dPremium += 1/n_f * exp(-dJumpInt * dOptionExpiry) * pow(dJumpInt * dOptionExpiry, n) * bs;

	}

	*out_dPrice = dPremium;
	
// Calculate Implied Vol
	if (err = srt_f_optimpvol(dPremium, dCashFwd, dCashStrike, dOptionExpiry, dOptionExpiryDf, enumCallPut, enumDiffType, out_dVol))
		return  err;
	

	
// Check that the the price is greater than the intrinsic
//	switch ( enumGreek )
//	{
//		case SRT_PREMIUM:
//			*out_dPrice = *out_dPrice < dIntrinsic ? dIntrinsic : *out_dPrice;
//			break;
//		case DELTA:
//			*out_dPrice = - (*out_dPrice);
//			break;
//		default:
//			break;
//	}


// Return the price
	return 0;
}	
// A Eurodollar option pricer


/* ---------------------------------------------------------------------------------------------------------------- 

  Functions for the function

									f(t) = ( a + bt + gt^2 )exp{-ct} + d

  ---------------------------------------------------------------------------------------------------------------- */

 /* ---------------------------------------------------------------------------------------------------------------- 

  Functions for the function

									f(t) = ( a + bt + gt^2 )exp{-ct} + d

  ---------------------------------------------------------------------------------------------------------------- */


struct PL_parameters
{
	double a, b, c, d, g;
};


struct PL_input
{
	double inital, asymptotic, max, time, width;
};


char* PL_errorFunction( double index , double in_par[], double* out_y, double* out_dyda, int in_nPar )
{

}

char* PL_convertInput( struct PL_input* input, struct PL_parameters* out_p )
{
/* Variable declaration */
	double data[4], target[4], weight[4], param[6];
	double chisq;
	char* err = 0;

/* Set the parameters that we can */
	out_p->d = input->asymptotic;
	out_p->a = input->inital - input->asymptotic;


/* Initial guess */
	out_p->b = out_p->a / ( 2.0 * input->width - input->time ); 
	out_p->c = out_p->b / ( out_p->b * input->time + out_p->a ); 
	out_p->g = 0.0; 

	data[1] = out_p->b;
	data[2] = out_p->c;
	data[3] = out_p->g;

/* weights */
	weight[1] = 1;
	weight[2] = 1;
	weight[3] = 1;

/* Solve */
/*	if ( err = levenberg_marquardt(
										data,   // From [1] to [ndata] 
										target, // From [1] to [ndata] 
										weight, // From [1] to [ndata] 
										3,
										param,  // From [1] to [nparam] 
										long nparam,
										long niter,
										Err (*funcs)(double, double[], double*, double[], int),
										chisq)	)
						return err;
*/
/* Transfer the solution */
	out_p->b = data[1];
	out_p->c = data[2];
	out_p->g = data[3];

/* return */
	return 0;
}



