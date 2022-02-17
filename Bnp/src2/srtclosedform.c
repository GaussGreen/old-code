/**********************************************************************
 *      Name: SrtClosedForm.c                                         * 
 *  Function: Entry point to srt_f_grfn_clsdfrm with raw data         *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 15/11/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects sUndPtr and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 15/11/95 FOS     Created for SORT5-GRFN3 port to NT based on the   *
 *                  gr_clsdfrm_addin of 2020                          *
 **********************************************************************/
#include "srt_h_all.h"
#include "grf_h_public.h"
#include "srt_h_grfclsdfrm.h"
#include "SrtAccess.h" 
#include "opfnctns.h"
	  
char *SrtClosedForm(
					long    start, 
					long    nfp_or_end, 
					char    *cpdStr,
				    char    *basisStr, 
					char    *recPayStr, 
					char    *typeStr,
                    double  strike, 
					double  bondStrike,
					char    *refRateCodeStr,
                    int     numParams, 
					char    **paramStrings, 
					char    **valueStrings,
                    char    *sUndPtrName, 
					double  *price, 
					double  **sigma_vega, 
					long    *num_sig, 
					double  **tau_vega, 
					long    *num_tau,
					char    *greekStr)
{
Err                err = NULL;
int                status = 0;
SrtUndPtr          sUndPtr;
SrtGrfnParam       grfnparam;
SwapDP             sdp;
StructType         type;
SrtReceiverType    rec_pay;
SrtGreekType       greek;

 /* Initialise and Overwrite defaults parameters with user defined parameters */
	if (err = srt_f_set_GrfnParams(numParams,paramStrings,
								  valueStrings,&grfnparam))
	{
	  return err;
	}

 /* Gets the sUndPtrerlying through its name */
	if (sUndPtrName)
	{
		sUndPtr = lookup_und(sUndPtrName);
	}
	else
	{
		return serror("Heinous error: No sUndPtrerlying passed to SrtClosedForm");
	} 

/* Checks the sUndPtrerlying type */
	if (sUndPtr)
	{
	   if (!ISUNDTYPE(sUndPtr,INTEREST_RATE_UND))
	   {
          return serror("sUndPtrerlying must be of type Interest Rate");
	   }
	}
	else
	{
	   return serror("Could not find sUndPtrerlying in market list");
	}

/* Populates the SwapDp structire with Spot Lag, Start Date,...*/
	sdp.spot_lag = get_spotlag_from_underlying(sUndPtr);

    if (err = swp_f_initSwapDP(start,nfp_or_end,cpdStr,basisStr,&sdp))
	{
	   return err;
	}
	

/* Interprets the Rec Pay string into a type  */
    if (err = interp_rec_pay(recPayStr,&rec_pay))
	{
	   return err;
	}

/* Interprets the Product type (SWAPTION, CAPFLOOR,...) */
	if (err = interp_struct(typeStr,&type))
	{
	   return err;
	}

/* Interprets the Greek type if it exists, or set it to default */
	greek = PREMIUM;
	
	if (strcmp(greekStr, "") != 0)
	{
		if (err = interp_greeks(greekStr, &greek))
		{
			return err;
		}
	}
		
/* Call the relevant function */
	if (greek == PREMIUM)
	{
		*num_sig = 0;
		*num_tau = 0;

		(*sigma_vega) = (double *) srt_malloc( (*num_sig) * sizeof(double));
		(*tau_vega) = (double *) srt_malloc( (*num_tau) * sizeof(double)); 

		err = srt_f_grfn_clsdfrm(sUndPtr,&grfnparam,&sdp,strike,bondStrike,rec_pay,type,refRateCodeStr, 
			price);
	}
	else if (greek == VEGA)
	{
		err = srt_f_grfn_clsdfrm_vega(sUndPtr,&grfnparam,&sdp,strike,bondStrike,rec_pay,type,
			refRateCodeStr, price,sigma_vega,num_sig,tau_vega,num_tau);
	}
	else
	{
		err = "Greek should be premium or vega ";
	}

/* Return a error message (or success) */
	return err;
}

/* --------------------------------------------------------------------------- */

Err   SrtFutureClosedForm(
					long       lFutureDate,
					long       lStartDate, 
					long       lNfpOrEnd, 
					char     *szComp,
				    char     *szBasis, 
					char     *szRecPay, 
					char     *szType,
                    double     dStrike, 
					double     dBondStrike,
					char     *szRefRateCode,
                    int        iNumParams, 
					char   **pszParamStrings, 
					char   **pszValueStrings,
                    char     *szsUndPtrName, 
					double   *pdFuturePrice)
{
Err                err = NULL;
int                status = 0;
SrtUndPtr          sUndPtr;
SrtGrfnParam       grfnparam;
SwapDP             sdp;
StructType         type;
SrtReceiverType    rec_pay;

/* Checks the Future Date is before the Start Date */
	if (lFutureDate > lStartDate)
		return serror("Start Date before Future Date in SrtFutureClosedForm");

/* Initialise and Overwrite Grfn defaults parameters with user defined parameters */
	if (err = srt_f_set_GrfnParams(iNumParams,pszParamStrings,
			pszValueStrings,&grfnparam))
	{
		return err;
	}

/* Gets the sUndPtrerlying through its name */
	if (szsUndPtrName)
	{
		sUndPtr = lookup_und(szsUndPtrName);
	}
	else
	{
		return serror("No Underlying passed to SrtFutureClosedForm");
	} 

/* Checks the sUndPtrerlying type */
	if (sUndPtr)
	{
		if (!ISUNDTYPE(sUndPtr,INTEREST_RATE_UND))
		{
			return serror("Underlying must be of type Interest Rate");
		}
	}
	else
	{
		return serror("Could not find sUndPtrerlying in market list");
	}

/* Populates the SwapDp structire with Spot Lag, Start Date,...*/
	sdp.spot_lag = get_spotlag_from_underlying(sUndPtr);
	if (err = swp_f_initSwapDP(lStartDate,lNfpOrEnd,szComp,szBasis,&sdp))
	{
	   return err;
	}
	

/* Interprets the Rec Pay string into a type  */
    if (err = interp_rec_pay(szRecPay,&rec_pay))
	{
	   return err;
	}

/* Interprets the Product type (SWAPTION, CAPFLOOR,...) */
	if (err = interp_struct(szType,&type))
	{
	   return err;
	}

/* Call to the main function */
	err = srt_f_grfn_future_closedform(
				lFutureDate,
				sUndPtr,
				&grfnparam,
				&sdp,
				dStrike,
				dBondStrike,
				rec_pay,
				type,
				szRefRateCode,
				pdFuturePrice);
	if(err)
	{
		return err;
	}

/* Return a success message */
	return NULL;

}/* END SrtFutureClosedForm(...) */


/* ------------------------------------------------------------------------------ */

Err   SrtFutureVolatility(
					long       lFutureDate,
					long       lStartDate, 
					long       lNfpOrEnd, 
					char     *szComp,
				    char     *szBasis, 
					char     *szRecPay, 
					char     *szType,
                    double     dStrike, 
					double     dBondStrike,
					char     *szRefRateCode,
                    int        iNumParams, 
					char   **pszParamStrings, 
					char   **pszValueStrings,
                    char     *szsUndPtrName, 
					char    *szVolType,
					double   *pdFutureVolatility)
{
Err                err = NULL;
int                status = 0;
SrtUndPtr          sUndPtr;
SrtGrfnParam       grfnparam;
SwapDP             sdp;
StructType         type;
SrtReceiverType    eRecPay;
double             dFuturePrice;
String            szYcName;
long               lFixingDate;
long               lSpotLag;
double             dDfFutureDate;
double             dLevel;
double             dForward;
SrtDiffusionType   eVolType;
SrtCallPutType     eCallPut;

/* Checks the Future Date is before the Start Date */
	if (lFutureDate > lStartDate)
		return serror("Start Date before Future Date in SrtFutureClosedForm");

/* Initialise and Overwrite Grfn defaults parameters with user defined parameters */
	if (err = srt_f_set_GrfnParams(iNumParams,pszParamStrings,
			pszValueStrings,&grfnparam))
	{
		return err;
	}

/* Gets the sUndPtrerlying through its name */
	if (szsUndPtrName)
	{
		sUndPtr = lookup_und(szsUndPtrName);
	}
	else
	{
		return serror("No Underlying passed to SrtFutureClosedForm");
	} 

/* Checks the sUndPtrerlying type */
	if (sUndPtr)
	{
		if (!ISUNDTYPE(sUndPtr,INTEREST_RATE_UND))
		{
			return serror("Underlying must be of type Interest Rate");
		}
	}
	else
	{
		return serror("Could not find sUndPtrerlying in market list");
	}

/* Populates the SwapDp structire with Spot Lag, Start Date,...*/
	sdp.spot_lag = get_spotlag_from_underlying(sUndPtr);
	if (err = swp_f_initSwapDP(lStartDate,lNfpOrEnd,szComp,szBasis,&sdp))
	{
	   return err;
	}
	

/* Interprets the Rec Pay string into a type  */
    if (err = interp_rec_pay(szRecPay,&eRecPay))
	{
	   return err;
	}

/* Interprets the Product type (SWAPTION, CAPFLOOR,...) */
	if (err = interp_struct(szType,&type))
	{
	   return err;
	}

/* Call to the main function to compute the future price (using Grfn) */
	err = srt_f_grfn_future_closedform(
				lFutureDate,
				sUndPtr,
				&grfnparam,
				&sdp,
				dStrike,
				dBondStrike,
				eRecPay,
				type,
				szRefRateCode,
				&dFuturePrice);
	if(err)
	{
		return err;
	}

/* Gets th Yield Curve through the name attached to the sUndPtrerlying */
	szYcName = get_discname_from_underlying(sUndPtr) ;

/* Computes the forward rate for the relevant period (to compute an implied volatility) */
	err = swp_f_ForwardRate (lStartDate, 
							 lNfpOrEnd, 
							 szComp, 
							 szBasis, 
							 szYcName, 
							 szRefRateCode, 
							 &dForward);
	if (err)
	{
		return serror ("Error in SrtFutureVolatility (srt_f_ForwardRate): %s", err);
	}

/* Compute the discount factor from Today until Future Date */
	dDfFutureDate = swp_f_df(0.0,lFutureDate,szYcName);

/* Computes the Level Payment for the relevant period (to compute an implied volatility) */
	err = swp_f_LevelPayment (lStartDate,
							  lNfpOrEnd, 
							  szComp, 
							  szBasis, 
							  szYcName, 
							  szRefRateCode, 
							  &dLevel);
	if (err)
	{
		return serror ("Error in SrtFutureVolatility (srt_f_LevelPayment): %s", err);
	}
	dLevel /= dDfFutureDate, 
	
/* Gets the Fixing Date (Spot Lag "BD" before the Start Date)  */
	lSpotLag  = get_spotlag_from_underlying(sUndPtr);
    lFixingDate = add_unit (lStartDate, - lSpotLag, SRT_BDAY, MODIFIED_SUCCEEDING);
	
/*  CALL or PUT ? */
	if (eRecPay == SRT_PAYER)
		eCallPut = SRT_CALL;
	else if (eRecPay == SRT_RECEIVER)
		eCallPut = SRT_PUT;

/*	Lognormal / Normal */
	err = interp_diffusion_type (szVolType, &eVolType);
	if (err)
	{
		return serror ("Error in SrtFutureVolatility (interp_diffusion_type): %s", err);
	}

/*  Compute the Implied Volatility for this Future Price */
	err = srt_f_optimpvol (dFuturePrice , dForward, 
		dStrike, 
		(lFixingDate - lFutureDate) * YEARS_IN_DAY, 
		dLevel , 
		eCallPut, 
		eVolType, 
		pdFutureVolatility);

	if (err)
	{
		return serror ("Error in SrtFutureVolatility (srt_f_optimpvol): %s", err);
	}

	
/* Return a success message */
	return NULL;

}/* END SrtFutureClosedForm(...) */


/* ------------------------------------------------------------------------------------------ */