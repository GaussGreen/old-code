
#include "math.h"
#include "opfnctns.h"
#include "GenericMidatAutocal.h"
#include "GenericMidatApplications.h"
#include "DiagCalibDLM.h"
#include "swp_h_vol.h"

Err	GMA_Midat_allocate_params(	int				iNbShortVolShift,
								double			*dShortVolShiftTime,
								double			*dShortVolShift,
								GMA_MIDATPARAMS	sParams)
{
Err	err = NULL;

	sParams->iNbShortVolShift = iNbShortVolShift;

	if (dShortVolShift && dShortVolShiftTime)
	{
		sParams->dShortVolShiftTime = calloc(iNbShortVolShift, sizeof(double));
		sParams->dShortVolShift = calloc(iNbShortVolShift, sizeof(double));

		if (!dShortVolShiftTime || !dShortVolShift)
		{
			err = "Memory allocation faillure in GMA_Midat_allocate_params";
			return err;
		}

		memcpy(sParams->dShortVolShiftTime, dShortVolShiftTime, iNbShortVolShift * sizeof(double));
		memcpy(sParams->dShortVolShift, dShortVolShift, iNbShortVolShift * sizeof(double));
	}

	sParams->dVolShift = 0.0;
	sParams->iVolShiftType = 0;

	return NULL;
}

void GMA_MIDAT_set_default_params(GMA_MIDATPARAMS	sParams)
{
	sParams->iCalibATM = 1;
	sParams->iCalibFwdVol = 0;
	sParams->cCorrelName[0] = '\0';

	sParams->iNbShortVolShift = 0;
	sParams->dShortVolShift = NULL;
	sParams->dShortVolShiftTime = NULL;

	sParams->iNewAlgo = 1;
}

void GMA_MIDAT_free_params(GMA_MIDATPARAMS	sParams)
{
	if (sParams)
	{
		if (sParams->dShortVolShift) free (sParams->dShortVolShift);
		if (sParams->dShortVolShiftTime) free (sParams->dShortVolShiftTime);
	}
}

Err GMA_MidatAutocal(		/* Cash Flows description */
						long	lToday,
						char	*cYieldCurve,
						char	*cVolCurve,

						/*		Generic Info */
						long	lStartDate,
						long	lTheoEndDate,

						/* Fixed Leg */
						int		iNbFixedCoupon,
						double	*dCoupon,
						long	*lCouponStartDate,
						long	*lCouponEndDate,
						long	*lCouponPayDate,
						char	*cCouponBasis,

						/* Floating Leg */
						int		iNbFixedFunding,
						char	*cRefRate,
						double	*dFundingMargin,
						long	*lFundingStartDate,
						long	*lFundingEndDate,
						long	*lFundingPayDate,
						char	*cFundingBasis,

						/* Exercise Dates */
						int		iNbExercise,
						long	*lExerciseDate,
						long	*lSettlementDate,
						double	*lExerciseFee,

						int		iPayRec,

						/* Model parameters */						
						GENMIDAT_MODEL			sModel,

						/* Parameters */
						GENMIDAT_AUTOCALPARAMS	sAutocalParams,
						GENMIDAT_CALIBPARAMS	sCalibParams,
						GENMIDAT_PDEPAMS		sPDEParams,

						/* Extra Parameters */
						GMA_MIDATPARAMS			sMidatParams,
				
						/* Outputs */						
						double	*dIV,
						double	*dCall,
						GENMIDAT_AUTOCALINFO	sInfos)
{
double	*dFras				= NULL,
		*dLongLevel			= NULL,
		*dShortLevel		= NULL,
		*dExeTimes			= NULL,

		*dLongForward		= NULL,
		*dShortForward		= NULL,
		*dLongOptions		= NULL,
		*dShortOptions		= NULL,

		*dATMLongForward	= NULL,
		*dATMShortForward	= NULL,
		*dATMLongOptions	= NULL,
		*dATMShortOptions	= NULL;

int		i, j, i0;
double			dDfStart, dDfEnd, dSpread, dCoverage;
double			dShortForwardi, dLongForwardi, dStrike, dShortVolatility, dLongVolatility, dCorrel;
long			lCorrelStart, lCorrelEnd;
double			dCorrelStrike, dLongWeight, dShortWeight, dPartialVolatility;
SrtBasisCode	sFundBasis, sCouponBasis;
int		iIndexCoupon;
double	power, power2;
SrtCallPutType	sOptType;
double	dVolShift;


Err		err = NULL;

	/* Initialisation */
	*dIV = 0.0;
	*dCall = 0.0;

	if (iPayRec > 0.5)
	{
		sOptType = SRT_PUT;
	}
	else
	{
		sOptType = SRT_CALL;
	}

	if (sMidatParams->iCalibFwdVol)
	{
		if (sModel->iNbFactor > 1)
		{
			err = "2 Factor doesn't support Fwd Vol Calibration for now";
			goto FREE_RETURN;
		}

		sCalibParams->iShortIsCap = 0;
	}

	/* Calculate the Short Forward ATM */
	dFras = calloc(iNbFixedFunding, sizeof(double));
	dExeTimes = calloc(iNbFixedFunding, sizeof(double));

	dShortLevel = calloc(iNbFixedCoupon, sizeof(double));
	dLongLevel = calloc(iNbFixedCoupon, sizeof(double));
	dATMShortForward = calloc(iNbFixedCoupon, sizeof(double));
	dATMLongForward = calloc(iNbFixedCoupon, sizeof(double));
	dShortForward = calloc(iNbFixedCoupon, sizeof(double));
	dLongForward = calloc(iNbFixedCoupon, sizeof(double));

	dATMShortOptions = calloc(iNbFixedCoupon, sizeof(double));
	dATMLongOptions = calloc(iNbFixedCoupon, sizeof(double));
	dShortOptions = calloc(iNbFixedCoupon, sizeof(double));
	dLongOptions = calloc(iNbFixedCoupon, sizeof(double));

	if (!dFras || !dExeTimes || !dShortLevel || !dLongLevel || !dATMShortForward || !dATMLongForward ||
		!dShortForward || !dLongForward || !dATMShortOptions || !dATMLongOptions || !dShortOptions || !dLongOptions)
	{
		err = "Memory allocation faillure in GMA_MidatAutocal";
		goto FREE_RETURN;
	}

	err = interp_basis (cFundingBasis, &sFundBasis);
	
	if (err) goto FREE_RETURN;	

	iIndexCoupon = 0;

	for (i=0; i<iNbFixedFunding; i++)
	{		
		dDfStart = swp_f_df(lToday, lFundingStartDate[i], cYieldCurve);
		dDfEnd = swp_f_df(lToday, lFundingPayDate[i], cYieldCurve);
		dSpread = swp_f_spread(lToday, lFundingPayDate[i], cRefRate);
		dCoverage = coverage(lFundingStartDate[i], lFundingEndDate[i], sFundBasis);

		dFras[i] = (dDfStart - dDfEnd) / (dCoverage * dDfEnd) + dSpread;
		
		dATMShortForward[iIndexCoupon] += dFras[i] * dCoverage * dDfEnd;
		dShortForward[iIndexCoupon] -= (dFras[i] + dFundingMargin[i]) * dCoverage * dDfEnd;
		
		if (i < iNbFixedFunding - 1 && lCouponPayDate[iIndexCoupon] < lFundingPayDate[i+1] - 5)
		{
			for (j=0; j<=iIndexCoupon; j++)
			{
				dATMLongForward[j] += dATMShortForward[iIndexCoupon];
				dLongForward[j] += dShortForward[iIndexCoupon];
			}

			iIndexCoupon++;
		}
	}

	/* Last Coupon */
	for (j=0; j<=iIndexCoupon; j++)
	{
		dATMLongForward[j] += dATMShortForward[iIndexCoupon];
		dLongForward[j] += dShortForward[iIndexCoupon];
	}

	err = interp_basis (cCouponBasis, &sCouponBasis);
	
	if (err) goto FREE_RETURN;

	/* Calculate the Short and LongForward ATM */		

	for (i=0; i<iNbFixedCoupon; i++)
	{
		dDfEnd = swp_f_df(lToday, lCouponPayDate[i], cYieldCurve);
		dCoverage = coverage(lCouponStartDate[i], lCouponEndDate[i], sCouponBasis);
		dShortLevel[i] = dCoverage * dDfEnd;
		
		dShortForward[i] += dShortLevel[i] * dCoupon[i];

		for (j=0; j<=i; j++)
		{
			dLongLevel[j] += dShortLevel[i];
			dLongForward[j] += dShortLevel[i] * dCoupon[i];
		}
	}

	*dIV = -dLongForward[0] * iPayRec * 1.0;

	/* Option Calculations */
	/* skip non calibration dates */
	i0 = 0;
	while (i0 < iNbFixedCoupon && lCouponStartDate[i0] < lSettlementDate[0] - 5)
	{
		i0++;
	}

	for (i=i0; i<iNbFixedCoupon; i++)
	{
		dExeTimes[i] = (lExerciseDate[i-i0] - lToday) * YEARS_IN_DAY;

		/* Calculate Forward */
		dShortForwardi = dATMShortForward[i] / dShortLevel[i];
		dLongForwardi = dATMLongForward[i] / dLongLevel[i];		

		/* SHORT */				
		if (sMidatParams->dShortVolShift && sMidatParams->dShortVolShiftTime)
		{
			/* shift the volatility */
			constant_or_linear_interpolation_dlm(	sMidatParams->dShortVolShiftTime,
													sMidatParams->dShortVolShift,
													sMidatParams->iNbShortVolShift,
													dExeTimes[i],
													&dVolShift);
		}
		else
		{
			dVolShift = 0.0;
		}
		
		if (sMidatParams->iCalibATM)
		{
			err = swp_f_vol(cVolCurve, lCouponStartDate[i], lCouponEndDate[i], dShortForwardi, &dShortVolatility, &power);
			
			if (err) goto FREE_RETURN;

			dShortVolatility += dVolShift * dShortVolatility;
			dShortVolatility += sMidatParams->dVolShift * dShortVolatility;

			if (power > 0.5)
			{
				dATMShortOptions[i] = srt_f_optblksch(	dShortForwardi,
														dShortForwardi,
														dShortVolatility,
														dExeTimes[i],
														dShortLevel[i],
														sOptType,
														PREMIUM);

				if (sCalibParams->iShortIsCap == 0)
				{
					/* Convert into Normal Vol */
					err = srt_f_optimpvol(	dATMShortOptions[i],
											dShortForwardi,
											dShortForwardi,
											dExeTimes[i],
											dShortLevel[i],
											sOptType,
											SRT_NORMAL,
											&dShortVolatility);

					if (err) goto FREE_RETURN;
				}
			}
			else
			{
				dATMShortOptions[i] = srt_f_optblknrm(	dShortForwardi,
														dShortForwardi,
														dShortVolatility,
														dExeTimes[i],
														dShortLevel[i],
														sOptType,
														PREMIUM);
			}			
		}
		else
		{
			dStrike = (dShortForward[i] + dATMShortForward[i]) / dShortLevel[i];

			err = swp_f_vol(cVolCurve, lCouponStartDate[i], lCouponEndDate[i], dStrike, &dShortVolatility, &power);

			if (err) goto FREE_RETURN;

			dShortVolatility += dVolShift * dShortVolatility;
			dShortVolatility += sMidatParams->dVolShift * dShortVolatility;

			if (power > 0.5)
			{
				dShortOptions[i] = srt_f_optblksch(	dShortForwardi,
													dStrike,
													dShortVolatility,
													dExeTimes[i],
													dShortLevel[i],
													sOptType,
													PREMIUM);
			}
			else
			{
				dShortOptions[i] = srt_f_optblknrm(	dShortForwardi,
													dStrike,
													dShortVolatility,
													dExeTimes[i],
													dShortLevel[i],
													sOptType,
													PREMIUM);
			}
		}

		/* LONG */						
		if (sMidatParams->iCalibATM)
		{
			err = swp_f_vol(cVolCurve, lCouponStartDate[i], lTheoEndDate, dLongForwardi, &dLongVolatility, &power);
			if (err) goto FREE_RETURN;

			dLongVolatility += sMidatParams->dVolShift * dLongVolatility;

			if (power > 0.5)
			{
				dATMLongOptions[i] = srt_f_optblksch(	dLongForwardi,
														dLongForwardi,
														dLongVolatility,
														dExeTimes[i],
														dLongLevel[i],
														sOptType,
														PREMIUM);

				if (sCalibParams->iShortIsCap == 0)
				{
					/* Convert into Normal Vol */
					err = srt_f_optimpvol(	dATMLongOptions[i],
											dLongForwardi,
											dLongForwardi,
											dExeTimes[i],
											dLongLevel[i],
											sOptType,
											SRT_NORMAL,
											&dLongVolatility);

					if (err) goto FREE_RETURN;
				}
			}
			else
			{
				dATMLongOptions[i] = srt_f_optblknrm(	dLongForwardi,
														dLongForwardi,
														dLongVolatility,
														dExeTimes[i],
														dLongLevel[i],
														sOptType,
														PREMIUM);
			}
		}

		/* Adjust for forward Volatility if needed */
		if (sCalibParams->iShortIsCap == 0 && sMidatParams->iCalibATM)		
		{
			if (i < iNbFixedCoupon - 1)
			{
				/* Get the correlation */
				dCorrelStrike = dExeTimes[i] / 100.0;
				lCorrelStart = lToday + lCouponEndDate[i] - lCouponStartDate[i];
				lCorrelEnd = lCorrelStart + lTheoEndDate - lCouponStartDate[i];

				err = swp_f_vol(sMidatParams->cCorrelName, lCorrelStart, lCorrelEnd, dCorrelStrike, &dCorrel, &power2);
				if (err) goto FREE_RETURN;

				dLongWeight = dLongLevel[i] / (dLongLevel[i] - dShortLevel[i]);
				dShortWeight = dShortLevel[i] / (dLongLevel[i] - dShortLevel[i]);

				dPartialVolatility = dLongWeight * dLongWeight * dLongVolatility * dLongVolatility
									+ dShortWeight * dShortWeight * dShortVolatility * dShortVolatility
									- 2.0 * dCorrel * dLongWeight * dShortWeight * dLongVolatility * dShortVolatility;

				if (dPartialVolatility < 1.0E-08)
				{
					err = serror("Correlation %.2f between %.2f y %.2f y and %.2f y %.2f is too low",
							dCorrel, dExeTimes[i], (lCouponEndDate[i] - lCouponStartDate[i]) * DAYS_IN_YEAR, dExeTimes[i], (lTheoEndDate - lCouponStartDate[i]) * DAYS_IN_YEAR);
				}

				dPartialVolatility = sqrt(dPartialVolatility);
				
				dATMShortOptions[i] = srt_f_optblknrm(	dATMLongForward[i+1] / dLongLevel[i+1],
														dATMLongForward[i+1] / dLongLevel[i+1],
														dPartialVolatility,
														dExeTimes[i],
														dLongLevel[i+1],
														sOptType,
														PREMIUM);
			}
			else
			{
				dShortOptions[i] = 0.0;
			}
		}

		dStrike = (dLongForward[i] + dATMLongForward[i]) / dLongLevel[i];

		err = swp_f_vol(cVolCurve, lCouponStartDate[i], lTheoEndDate, dStrike, &dLongVolatility, &power);
		if (err) goto FREE_RETURN;

		dLongVolatility += sMidatParams->dVolShift * dLongVolatility;

		if (power > 0.5)
		{
			dLongOptions[i] = srt_f_optblksch(	dLongForwardi,
													dStrike,
													dLongVolatility,
													dExeTimes[i],
													dLongLevel[i],
													sOptType,
													PREMIUM);
		}
		else
		{
			dLongOptions[i] = srt_f_optblknrm(	dLongForwardi,
													dStrike,
													dLongVolatility,
													dExeTimes[i],
													dLongLevel[i],
													sOptType,
													PREMIUM);
		}
	}

	/* Set the numeraire */
	sModel->dNumeraire = dLongLevel[i0];

	for (i=0; i<iNbFixedCoupon; i++)
	{
		dATMLongForward[i] = 0.0;
		dLongForward[i] /= sModel->dNumeraire;
		dLongForward[i] *= iPayRec * 1.0;
	}
		
	if (sMidatParams->iCalibATM)
	{
		sCalibParams->iCalibLambda = 1;
		
		/* Initialialise the Model */
		err = genmidat_init_model(	sModel->dAlpha,
									sModel->dGamma,
									sModel->dRho,
									&(dExeTimes[i0]),
									&(dATMLongForward[i0]),
									NULL,
									NULL,
									NULL,
									NULL,
									sModel->dStartBeta2,
									sModel->dNumeraire,
									sModel);

		if (err) goto FREE_RETURN;

		/* First ATM Calibration */
		err = generic_midat_calib(	iNbExercise,
									&(dATMLongOptions[i0]),
									&(dATMShortOptions[i0]),
									sModel,
									sCalibParams,
									NULL);

		if (err) goto FREE_RETURN;

		sCalibParams->iCalibLambda = 0;
	}

	/* Call Midat Autocal */
	/* Initialialise the Model */
	err = genmidat_init_model(	sModel->dAlpha,
								sModel->dGamma,
								sModel->dRho,
								&(dExeTimes[i0]),
								&(dLongForward[i0]),
								NULL,
								NULL,
								NULL,
								NULL,
								sModel->dStartBeta2,
								sModel->dNumeraire,
								sModel);

	if (err) goto FREE_RETURN;

	err = GenericMidatAutocal(	iNbExercise,
								&(dLongOptions[i0]),
								&(dShortOptions[i0]),
								NULL,
								sModel,
								sAutocalParams,
								sCalibParams,
								sPDEParams,
								dCall,
								sInfos);

	if (err) goto FREE_RETURN;

	

FREE_RETURN:

	if (dFras) free (dFras);	
	if (dExeTimes) free (dExeTimes);

	if (dShortLevel) free (dShortLevel);
	if (dLongLevel) free (dLongLevel);
	if (dATMShortForward) free (dATMShortForward);
	if (dATMLongForward) free (dATMLongForward);
	if (dShortForward) free (dShortForward);
	if (dLongForward) free (dLongForward);

	if (dATMShortOptions) free (dATMShortOptions);
	if (dATMLongOptions) free (dATMLongOptions);
	if (dShortOptions) free (dShortOptions);
	if (dLongOptions) free (dLongOptions);

	return err;
}

Err GMA_MidatAutocalNew(		/* Cash Flows description */
						long	lToday,
						char	*cYieldCurve,
						char	*cVolCurve,
						int		iEODFixFlag,
						int		iEODPayFlag,
						int		iEODExeFlag,

						/*		Generic Info */
						long	lStartDate,
						long	lTheoEndDate,
						
						/* Floating Leg */
						int		iNbFunding,
						char	*cRefRate,
						long	*lFundingFixDate,
						long	*lFundingStartDate,
						long	*lFundingPayDate,						
						double	*dFundingCoverage,						
						double	*dFundingMargin,
						double	*dPastFixings,

						/* Fixed Leg */
						int		iNbCoupon,						
						long	*lCouponStartDate,
						long	*lCouponPayDate,
						double	*dCouponCoverage,
						double	*dCoupon,

						/* Exercise Dates */
						int		iNbExercise,
						long	*lExerciseDate,
						long	*lSettlementDate,
						double	*lExerciseFee,

						double	dPayRec,

						/* Model parameters */						
						GENMIDAT_MODEL			sModel,

						/* Parameters */
						GENMIDAT_AUTOCALPARAMS	sAutocalParams,
						GENMIDAT_CALIBPARAMS	sCalibParams,
						GENMIDAT_PDEPAMS		sPDEParams,

						/* Extra Parameters */
						GMA_MIDATPARAMS			sMidatParams,
				
						/* Outputs */						
						double	*dIV,
						double	*dCall,
						GENMIDAT_AUTOCALINFO	sInfos)
{
Err		err = NULL;

double	*dFras			= NULL,
		*dExoticPV		= NULL,
		*dFundingPV		= NULL,
		*dFlatFundingPV	= NULL;

double	*dExeTimes			= NULL,
		*dFwdIV				= NULL,
		*dATMFwdIV			= NULL,
		*dLongLevel			= NULL,
		*dLongForward		= NULL,
		*dATMShortOptions	= NULL,
		*dATMLongOptions	= NULL,
		*dATMLongVol		= NULL,
		*dLongOptions		= NULL;

int		iFundStartIdx, iCouponStartIdx, iExeStartIdx;
int		iFundIdx, iCouponIdx;
double	dMarketIV, dDfStart, dDfEnd, dSpread;
double	dLongForwardi, dLongLeveli, dFlatLongForward, dLongStrike, dLongATMVoli, dLongStrikeVol;
double	dShortForwardi, dShortLeveli, dShortATMVoli;
long	lCorrelStart, lCorrelEnd;
double	dCorrelStrike, dLongWeight, dShortWeight, dCorrel, dPartialVolatility;
double	dPower;
int		i, j;

SrtCallPutType	sOptType;

	/* Initialisation */
	*dIV = 0.0;
	*dCall = 0.0;

	if (dPayRec > 0.5)
	{
		sOptType = SRT_PUT;
	}
	else
	{
		sOptType = SRT_CALL;
	}

	if (sMidatParams->iCalibFwdVol)
	{
		if (sModel->iNbFactor > 1)
		{
			err = "2 Factor doesn't support Fwd Vol Calibration for now";
			goto FREE_RETURN;
		}

		sCalibParams->iShortIsCap = 0;
	}

	/* Memory Allocation */
	dFras = calloc(iNbFunding, sizeof(double));
	dExoticPV = calloc(iNbCoupon, sizeof(double));
	dFundingPV = calloc(iNbFunding, sizeof(double));
	dFlatFundingPV = calloc(iNbFunding, sizeof(double));

	dExeTimes = calloc(iNbExercise, sizeof(double));
	dFwdIV = calloc(iNbExercise, sizeof(double));
	dATMFwdIV = calloc(iNbExercise, sizeof(double));
	dLongLevel = calloc(iNbExercise, sizeof(double));
	dLongForward = calloc(iNbExercise, sizeof(double));
	dATMShortOptions = calloc(iNbExercise, sizeof(double));
	dATMLongOptions = calloc(iNbExercise, sizeof(double));
	dATMLongVol = calloc(iNbExercise, sizeof(double));
	dLongOptions = calloc(iNbExercise, sizeof(double));

	if (!dFras || !dExoticPV || !dFundingPV || !dFlatFundingPV ||
		!dExeTimes || !dFwdIV || !dATMFwdIV || !dLongLevel || !dLongForward ||
		!dATMShortOptions || !dATMLongOptions || !dATMLongVol || !dLongOptions)
	{
		err = "Memory allocation faillure in GMA_MidatAutocalNew";
		goto FREE_RETURN;
	}

	/* Cash-Flows computation */
	/* ********************** */

	dMarketIV = 0.0;

	/* Funding */	
	iFundStartIdx = 0;

	while (iFundStartIdx < iNbFunding && lFundingFixDate[iFundStartIdx] < lToday + iEODFixFlag)
	{
		if (lFundingPayDate[iFundStartIdx] < lToday + iEODPayFlag)
		{
			dFundingPV[iFundStartIdx] = 0.0;
		}
		else
		{
			dDfEnd = swp_f_df(lToday, lFundingPayDate[iFundStartIdx], cYieldCurve);
			dFundingPV[iFundStartIdx] = (dPastFixings[iFundStartIdx] + dFundingMargin[iFundStartIdx]) * dFundingCoverage[iFundStartIdx] * dDfEnd;
		}

		dMarketIV += dFundingPV[iFundStartIdx];

		iFundStartIdx++;
	}

	for (i=iFundStartIdx; i<iNbFunding; i++)
	{	
		dDfStart = swp_f_df(lToday, lFundingStartDate[i], cYieldCurve);
		dDfEnd = swp_f_df(lToday, lFundingPayDate[i], cYieldCurve);
		dSpread = swp_f_spread(lToday, lFundingPayDate[i], cRefRate);

		dFlatFundingPV[i] = dDfStart - (1.0 - dSpread) * dDfEnd;
		dFundingPV[i] = dFlatFundingPV[i] + dFundingMargin[i] * dFundingCoverage[i] * dDfEnd;		

		dMarketIV += dFundingPV[i];
	}

	/* Exotic */
	iCouponStartIdx = 0;

	while (iCouponStartIdx < iNbCoupon && lCouponPayDate[iCouponStartIdx] < lToday + iEODPayFlag)
	{
		dExoticPV[iCouponStartIdx] = 0.0;
		iCouponStartIdx++;
	}

	for (i=iCouponStartIdx; i<iNbCoupon; i++)
	{		
		dDfEnd = swp_f_df(lToday, lCouponPayDate[i], cYieldCurve);
		dExoticPV[i] = dCoupon[i] * dCouponCoverage[i] * dDfEnd;
		dMarketIV -= dExoticPV[i];
	}

	*dIV = dPayRec * dMarketIV;

	/* Calibration Instruments Calculation */
	/* *********************************** */

	iExeStartIdx = 0;

	while (iExeStartIdx < iNbExercise && lExerciseDate[iExeStartIdx] < lToday + iEODExeFlag)
	{
		iExeStartIdx++;
	}

	for (i=iExeStartIdx; i<iNbExercise; i++)
	{
		dExeTimes[i] = (lExerciseDate[i] - lToday) * YEARS_IN_DAY;

		/* Find the Funding Idx */
		iFundIdx = iFundStartIdx;

		while (iFundIdx < iNbFunding && lFundingStartDate[iFundIdx] < lExerciseDate[i])
		{
			iFundIdx++;
		}

		if (iFundIdx == iNbFunding)
		{
			err = serror("Exercise %d doesn't call any funding", i+1);
			goto FREE_RETURN;
		}

		/* Find the Exotic Idx */
		iCouponIdx = iCouponStartIdx;

		while (iCouponIdx < iNbCoupon && lCouponStartDate[iCouponIdx] < lExerciseDate[i])
		{
			iCouponIdx++;
		}		

		/* Compute the Long Forward */
		dLongForwardi = 0.0;
		dLongLeveli = 0.0;

		/* Add the Broken Period */
		if (iCouponIdx > 0)
		{
			dLongLeveli += (lCouponPayDate[iCouponIdx-1] - lSettlementDate[i]) / (1.0 * (lCouponPayDate[iCouponIdx-1] - lCouponStartDate[iCouponIdx-1])) * dExoticPV[iCouponIdx-1] / dCoupon[iCouponIdx-1];
			dLongForwardi += (lCouponPayDate[iCouponIdx-1] - lSettlementDate[i]) / (1.0 * (lCouponPayDate[iCouponIdx-1] - lCouponStartDate[iCouponIdx-1])) * dExoticPV[iCouponIdx-1];
		}

		for (j=iCouponIdx; j<iNbCoupon; j++)
		{
			dLongLeveli += dExoticPV[j] / dCoupon[j];
			dLongForwardi += dExoticPV[j];
		}

		dLongLevel[i] = dLongLeveli;
		
		/* Substract Funding */
		dFlatLongForward = 0.0;

		for (j=iFundIdx; j<iNbFunding; j++)
		{
			dLongForwardi -= dFundingPV[j];
			dFlatLongForward += dFlatFundingPV[j];
		}

		/* Add the Fee */
		dLongForwardi = dLongForwardi * dPayRec - lExerciseFee[i] * swp_f_df(lToday, lSettlementDate[i], cYieldCurve);

		dFwdIV[i] = dLongForwardi;

		dLongStrike = (dLongForwardi + dFlatLongForward) / dLongLeveli;
		dLongForward[i] = dFlatLongForward / dLongLeveli;

		err = swp_f_vol(cVolCurve, lSettlementDate[i], lTheoEndDate, dLongForward[i], &dLongATMVoli, &dPower);
		if (err) goto FREE_RETURN;

		err = swp_f_vol(cVolCurve, lSettlementDate[i], lTheoEndDate, dLongStrike, &dLongStrikeVol, &dPower);
		if (err) goto FREE_RETURN;

		if (dPower > 0.5)
		{
			dLongOptions[i] = srt_f_optblksch(	dLongForward[i],
												dLongStrike,
												dLongStrikeVol,
												dExeTimes[i],
												dLongLevel[i],
												sOptType,
												PREMIUM);

			dATMLongOptions[i] = srt_f_optblksch(	dLongForward[i],
													dLongForward[i],
													dLongATMVoli,
													dExeTimes[i],
													dLongLevel[i],
													sOptType,
													PREMIUM);

			/* Convert into Normal Vol */
			err = srt_f_optimpvol(	dATMLongOptions[i],
									dLongForward[i],
									dLongForward[i],
									dExeTimes[i],
									dLongLevel[i],
									sOptType,
									SRT_NORMAL,
									&dLongATMVoli);
		}
		else
		{
			dLongOptions[i] = srt_f_optblknrm(	dLongForward[i],
												dLongStrike,
												dLongStrikeVol,
												dExeTimes[i],
												dLongLevel[i],
												sOptType,
												PREMIUM);

			dATMLongOptions[i] = srt_f_optblknrm(	dLongForward[i],
													dLongForward[i],
													dLongATMVoli,
													dExeTimes[i],
													dLongLevel[i],
													sOptType,
													PREMIUM);
		}

		dATMLongVol[i] = dLongATMVoli;

		if (i > iExeStartIdx)
		{
			/* Compute the previous ATM Short Option */
			dShortLeveli = dLongLevel[i-1] - dLongLevel[i];
			dShortForwardi = (dLongLevel[i-1] * dLongForward[i-1] - dLongLevel[i] * dLongForward[i]) / dShortLeveli;

			err = swp_f_vol(cVolCurve, lSettlementDate[i-1], lSettlementDate[i], dShortForwardi, &dShortATMVoli, &dPower);
			if (err) goto FREE_RETURN;

			if (dPower > 0.5)
			{				
				dATMShortOptions[i-1] = srt_f_optblksch(dShortForwardi,
														dShortForwardi,
														dShortATMVoli,
														dExeTimes[i-1],
														dShortLeveli,
														sOptType,
														PREMIUM);

				/* Convert into Normal Vol */
				err = srt_f_optimpvol(	dATMShortOptions[i-1],
										dShortForwardi,
										dShortForwardi,
										dExeTimes[i-1],
										dShortLeveli,
										sOptType,
										SRT_NORMAL,
										&dShortATMVoli);
			}
			else
			if (sCalibParams->iShortIsCap)
			{				
				dATMShortOptions[i-1] = srt_f_optblknrm(dShortForwardi,
														dShortForwardi,
														dShortATMVoli,
														dExeTimes[i-1],
														dShortLeveli,
														sOptType,
														PREMIUM);
			}
			
			if (!sCalibParams->iShortIsCap)
			{
				dCorrelStrike = dExeTimes[i-1] / 100.0;
				lCorrelStart = lToday + lSettlementDate[i] - lSettlementDate[i-1];
				lCorrelEnd = lCorrelStart + lTheoEndDate - lSettlementDate[i-1];

				err = swp_f_vol(sMidatParams->cCorrelName, lCorrelStart, lCorrelEnd, dCorrelStrike, &dCorrel, &dPower);
				if (err) goto FREE_RETURN;

				dLongWeight = dLongLevel[i-1] / dLongLevel[i];
				dShortWeight = dShortLeveli / dLongLevel[i];

				dPartialVolatility = dLongWeight * dLongWeight * dATMLongVol[i-1] * dATMLongVol[i-1]
									+ dShortWeight * dShortWeight * dShortATMVoli * dShortATMVoli
									- 2.0 * dCorrel * dLongWeight * dShortWeight * dATMLongVol[i-1] * dShortATMVoli;

				if (dPartialVolatility < 1.0E-08)
				{
					err = serror("Correlation %.2f between %.2f y %.2f y and %.2f y %.2f is too low",
							dCorrel, dExeTimes[i-1], (lSettlementDate[i] - lSettlementDate[i-1]) * DAYS_IN_YEAR, dExeTimes[i-1], (lTheoEndDate - lSettlementDate[i-1]) * DAYS_IN_YEAR);
				}

				dPartialVolatility = sqrt(dPartialVolatility);
				
				dATMShortOptions[i-1] = srt_f_optblknrm(dLongForward[i],
														dLongForward[i],
														dPartialVolatility,
														dExeTimes[i-1],
														dLongLevel[i],
														sOptType,
														PREMIUM);
			}

		}
	}

	/* Set the numeraire */
	sModel->dNumeraire = dLongLevel[iExeStartIdx];

	for (i=iExeStartIdx; i<iNbExercise; i++)
	{
		dATMFwdIV[i] = 0.0;
		dFwdIV[i] /= sModel->dNumeraire;
	}

	sModel->iNbVols = iNbExercise - iExeStartIdx;
		
	if (sMidatParams->iCalibATM)
	{
		sCalibParams->iCalibLambda = 1;
		
		/* Initialialise the Model */
		err = genmidat_init_model(	sModel->dAlpha,
									sModel->dGamma,
									sModel->dRho,
									&(dExeTimes[iExeStartIdx]),
									&(dATMFwdIV[iExeStartIdx]),
									NULL,
									NULL,
									NULL,
									NULL,
									sModel->dStartBeta2,
									sModel->dNumeraire,
									sModel);

		if (err) goto FREE_RETURN;

		/* First ATM Calibration */
		err = generic_midat_calib(	iNbExercise - iExeStartIdx,
									&(dATMLongOptions[iExeStartIdx]),
									&(dATMShortOptions[iExeStartIdx]),
									sModel,
									sCalibParams,
									NULL);

		if (err) goto FREE_RETURN;

		sCalibParams->iCalibLambda = 0;
	}

	/* Call Midat Autocal */
	/* Initialialise the Model */
	err = genmidat_init_model(	sModel->dAlpha,
								sModel->dGamma,
								sModel->dRho,
								&(dExeTimes[iExeStartIdx]),
								&(dFwdIV[iExeStartIdx]),
								NULL,
								NULL,
								NULL,
								NULL,
								sModel->dStartBeta2,
								sModel->dNumeraire,
								sModel);

	if (err) goto FREE_RETURN;

	err = GenericMidatAutocal(	iNbExercise - iExeStartIdx,
								&(dLongOptions[iExeStartIdx]),
								NULL,
								NULL,
								sModel,
								sAutocalParams,
								sCalibParams,
								sPDEParams,
								dCall,
								sInfos);

	if (err) goto FREE_RETURN;	

FREE_RETURN:

	if (dFras) free(dFras);
	if (dExoticPV) free(dExoticPV);
	if (dFundingPV) free(dFundingPV);
	if (dFlatFundingPV) free(dFlatFundingPV);

	if (dExeTimes) free(dExeTimes);
	if (dFwdIV) free(dFwdIV);
	if (dATMFwdIV) free(dATMFwdIV);
	if (dLongLevel) free(dLongLevel);
	if (dLongForward) free(dLongForward);
	if (dATMShortOptions) free(dATMShortOptions);
	if (dATMLongOptions) free(dATMLongOptions);
	if (dATMLongVol) free(dATMLongVol);
	if (dLongOptions) free(dLongOptions);

	return err;
}
						

