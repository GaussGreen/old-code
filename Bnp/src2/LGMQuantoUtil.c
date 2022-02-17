/****************************************************************************/
/* Purpose */
/****************************************************************************/

/****************************************************************************/
/* EXTERNAL ROUTINES USED */
/****************************************************************************/
/*

*/
/****************************************************************************/
/* Headers */
/****************************************************************************/
#include "SRT_H_ALL.H>
#include "srt_h_lgmtypes.h"
#include "srt_h_lgmprotos.h"
#include "math.h"

/****************************************************************************/
/* -- Step 5 of the quanto autocal: Calibration of the Fx gamma function -- */
/****************************************************************************/

LGMErr LGMExtractAlphaFromTS(LGM_TS *tsPtrA, LGM_TS *tsPtrB,
							 /* Outputs two alpha term structure */
							 double **AlphaCcyA, double **AlphaCcyB)
{
	long	NumZCcyA, NumZCcyB, i;

	NumZCcyA = tsPtrA->numZ; 
	*AlphaCcyA = (double *) srt_malloc(NumZCcyA * sizeof(double));
	if (*AlphaCcyA == NULL) return "Allocation Error in Quanto Callable";

	NumZCcyB = tsPtrB->numZ; 
	*AlphaCcyB = (double *) srt_malloc(NumZCcyB * sizeof(double));
	if (*AlphaCcyB == NULL) 
	{
		free(*AlphaCcyA);
		return "Allocation Error in Quanto Callable";
	}
	(*AlphaCcyA)[0] = .0; (*AlphaCcyB)[0] = .0;

	for (i=1; i<NumZCcyA; i++)
		(*AlphaCcyA)[i] = sqrt(max((tsPtrA->zeta[i] - tsPtrA->zeta[i-1])
					/ ((double) (tsPtrA->zdate[i] - tsPtrA->zdate[i-1])), 0));

	for (i=1; i<NumZCcyB; i++)
		(*AlphaCcyB)[i] = sqrt(max((tsPtrB->zeta[i] - tsPtrB->zeta[i-1])
					/ ((double) (tsPtrB->zdate[i] - tsPtrB->zdate[i-1])), 0));

	return NULL;
}

LGMErr LGMExpandAlphaUAlphaV(long lNumU, Date *datesU, 
						    double *dAlphaU, /* AlphaU[0..NumU - 1] */
						    long lNumV, Date *datesV, 
						    double *dAlphaV, /* AlphaV[0..NumV - 1] */
						    long *lNum, Date **dates,
							double **dNewAlphaU, double **dNewAlphaV)
{
	long lNumMax, i, j, k;
	
	lNumMax = lNumU + lNumV;
	*dNewAlphaU = (double *) srt_malloc(lNumMax * sizeof(double));
	*dNewAlphaV = (double *) srt_malloc(lNumMax * sizeof(double));
	*dates = (Date *) srt_malloc(lNumMax * sizeof(Date));

	i = 0; j = 0; 
	for (k = 0; k < lNumMax || i < lNumU || j < lNumV; k++)
	{
		(*dNewAlphaU)[k] = dAlphaU[i]; 
		(*dNewAlphaV)[k] = dAlphaV[j];
		if (datesU[i] > datesV[j]) (*dates)[k] = datesV[j++];
		else
		{
			if (datesU[i] < datesV[j]) (*dates)[k] = datesU[i++];
			else if (datesU[i++] == datesV[j]) (*dates)[k] = datesV[j++];
		}
	}

	if (i == lNumU)
		for (i = j; i < lNumV; i++)
		{
			(*dNewAlphaV)[k] = dAlphaV[i];
			(*dates)[k++] = datesV[i++];
		}
	else 
	if (j == lNumV)
		for (j = i; j < lNumU; j++)
		{
			(*dNewAlphaV)[k] = dAlphaU[j];
			(*dates)[k++] = datesU[j++];
		}

	*lNum = k;
	*dNewAlphaU = (double*) realloc(*dNewAlphaU, (*lNum)*sizeof(double));
	*dNewAlphaV = (double*) realloc(*dNewAlphaV, (*lNum)*sizeof(double));
	*dates = (Date *) realloc(*dates, (*lNum)*sizeof(Date));
	return NULL;
}

double LGMIntAlpha(long		lNum, 
				   Date		*dates, 
				   double	*dAlpha, 
				   Date		lStartDate,
				   Date		lEndDate)
{
	long i;
	double dIntegral = .0;

	if (dates[0] > lEndDate) return .0;
	
	for (i = 0; i < lNum && lStartDate >= dates[i]; i++) ;
	
	dIntegral = dAlpha[i] * (dates[i] - lStartDate);
	
	for (;i < lNum && dates[i] < lEndDate; i++)
		dIntegral += dAlpha[i] * (dates[i] - dates[i-1]);

	return dIntegral + dAlpha[i] * (lEndDate - dates[i-1]);
}

double LGMIntAlphaUAlphaV(long		lNum, 
						  Date		*dates, 
						  double	*dAlphaU, 
						  double	*dAlphaV, 
						  Date		lStartDate,
						  Date		lEndDate)
{
	long i;
	double dIntegral = .0;

	if (dates[0] > lEndDate) return .0;
	
	for (i = 0; i < lNum && lStartDate >= dates[i]; i++) ;
	
	dIntegral = dAlphaU[i] * dAlphaV[i] * (dates[i] - lStartDate);
	
	for (;i < lNum && dates[i] < lEndDate; i++)
		dIntegral += dAlphaU[i] * dAlphaV[i] * (dates[i] - dates[i-1]);

	return dIntegral + dAlphaU[i] * dAlphaV[i] * (lEndDate - dates[i-1]);
}

/* Main function of the calibration of Fx gamma: */

LGMErr LGMCalibGammaFromTSnFx ( LGM_TS		*tsPtrA,
								LGM_TS		*tsPtrB,
								long		lNumFx,
								Date		*datesFx,
								double		*VolFx,
								double		dCorrelAB,
								double		dCorrelAFx,
								double		dCorrelBFx,
								/* outputs linked to lNumFx and datesFx*/
								double		**dgamma)
{
	double	*AlphaA = NULL, *AlphaB = NULL, 
			*AlphaCcyA = NULL, *AlphaCcyB = NULL;
	double	dInt2Fx = .0, dIntAFx =.0, dIntBFx = .0, 
			dIntAB = .0, dIntA = .0, dIntB = .0;
	double	dBi, dCi, dHAi, dHBi, dksiA, dksiB, dDeltatFx; 
	Date	*dDatesAB;
	long	i, lNumAB;
	LGMErr	error;

	*dgamma = NULL;
	*dgamma = (double *) srt_malloc(lNumFx * sizeof(double));
	if (*dgamma ==NULL) return "Error in CalibGammaFromTSnFx";

	if(error = LGMExtractAlphaFromTS(tsPtrA, tsPtrB, &AlphaA, &AlphaB))
		return error;
	if(error = LGMExpandAlphaUAlphaV(tsPtrA->numZ, tsPtrA->zdate, AlphaA,
									 tsPtrB->numZ, tsPtrB->zdate, AlphaB,
									 &lNumAB, &dDatesAB, &AlphaCcyA, &AlphaCcyB))
		return error;

	/* the Fx calibration is just a recursion if you rearrange a the terms */
	(*dgamma)[0] = .0;
	for (i = 1; i < lNumFx; i++)
	{
		dDeltatFx = datesFx[i] - datesFx[i-1];

		dHAi= LGMGFromTS(datesFx[i], tsPtrA);
		dHBi = LGMGFromTS(datesFx[i], tsPtrB);
		dksiA = LGMZetaFromTS(datesFx[i], 0, tsPtrA);
		dksiB = LGMZetaFromTS(datesFx[i], 0, tsPtrB);

		dCi = dInt2Fx + 2.0 * dCorrelBFx * dHBi * dIntBFx 
			  - 2.0 * dCorrelAFx * dHAi * dIntAFx 
			  + dHAi * dHAi * dksiA + dHBi * dHBi * dksiB
			  - 2.0 * dCorrelAB * dHAi * dHBi * dIntAB 
			  - VolFx[i] * VolFx[i] * dDeltatFx;

		dBi = 2 * (dCorrelBFx * dHBi * dIntB - dCorrelAFx * dHAi * dIntA);
		
		if(dBi * dBi >= 4.0 * dDeltatFx * dCi && dCi <= .0)
			(*dgamma)[i] = (sqrt(dBi * dBi - 4.0 * dDeltatFx * dCi) - dBi) / 2.0 / dDeltatFx;
		else
			(*dgamma)[i] = .0;

		dInt2Fx += (*dgamma)[i] * (*dgamma)[i] * dDeltatFx;

		dIntA = LGMIntAlpha(tsPtrA->numZ, tsPtrA->zdate, AlphaA, 
							datesFx[i-1], datesFx[i]);

		dIntB = LGMIntAlpha(tsPtrB->numZ, tsPtrB->zdate, AlphaB, 
							datesFx[i-1], datesFx[i]);

		dIntAFx += (*dgamma)[i] * dIntA;
		dIntBFx += (*dgamma)[i] * dIntB;
		dIntAB += LGMIntAlphaUAlphaV(lNumAB, dDatesAB, 
									 AlphaCcyA, AlphaCcyB, 
									 datesFx[i-1], datesFx[i]);
	}

	if(AlphaA != NULL) srt_free(AlphaA);
	if(AlphaB != NULL) srt_free(AlphaB);
	if(dDatesAB != NULL) srt_free(dDatesAB); 
	if(AlphaCcyA != NULL) srt_free(AlphaCcyA);
	if(AlphaCcyB != NULL) srt_free(AlphaCcyB);

	return error;
}

/****************************************************************************/
/* ------------ Step 6 of the quanto autocal: Matching moments,
				from 3 gaussian variables to 2 gaussian variables --------- */
/****************************************************************************/
/*

LGMErr LGMComputeLambdaForMatchM(LGM_TS		*tsPtrB,
								 long		lNumEx,
								 Date		*datesEx,	
								 long		lNumCashFlow,
								 Date		*PayDatesB,
								 Date		*StartDatesB,
								 Date		tNow,
								 char		*cYCName,
								 double		*db,
								 double		*dbReduce,
								 double		*dStrikes,
								 double		**dLambda)
{
	long	i, k;
	double	dSum1, dSum2, dCashFlow;
	*dLambda = (double *) srt_malloc(lNumEx * sizeof(double));
	if (*dLambda == NULL) return "Allocation Error in Quanto Callable";
	
	for (i = 1; i < lNumEx; i++)
	{
		dSum1 = dSum2 = .0;
		
		for(k = 1; k < lNumCashFlow; k++)
		{
			if(StartDatesB[k]>datesEx[i])
			{
				dCashFlow = swp_f_df(tNow, PayDatesB[k], cYCName)*db[k];
				dSum1 += dCashFlow;
				dSum2 += dCashFlow*LGMGFromTS(PayDatesB[k], tsPtrB);
			}
		}
	}
}
*/

/****************************************************************************/
/* ------------ Step 7 of the quanto autocal: Matching moments,
				from 3 gaussian variables to 2 gaussian variables --------- */
/****************************************************************************/

LGMErr LGMReducAndMatchMoments (LGM_TS		*tsPtrA,
								LGM_TS		*tsPtrB,
								long		lNumEx,
								Date		*datesEx,
								double		*dLambda,
								long		lNumFx,
								Date		*datesFx,
								double		*dGamma,
								double		dCorrelAB,
								double		dCorrelAFx,
								double		dCorrelBFx,
								double		**dSki12,
								double		**dSki22)
{
	double		*dGammaFxA, *dAlphaAFx, *dGammaFxB, *dAlphaBFx, 
				*AlphaA, *AlphaB;
	Date		*dDatesAFx, *dDatesBFx;
	long		lNumAFx, lNumBFx, i;
	double		dIntGB = .0, dIntB = .0, dIntA = .0, 
				dIntGA = .0, dIntG2 = .0, dIntAB = .0,
				m = 1.0, dDeltaLambda, dksiB;
	LGMErr		error;

	if(error = LGMExtractAlphaFromTS(tsPtrA, tsPtrB, &AlphaA, &AlphaB))
		return error;

	if(error = LGMExpandAlphaUAlphaV(tsPtrA->numZ, tsPtrA->zdate, AlphaA,
									 lNumFx, datesFx, dGamma,
									 &lNumAFx, &dDatesAFx, &dAlphaAFx, &dGammaFxA))
		return error;

	if(error = LGMExpandAlphaUAlphaV(tsPtrB->numZ, tsPtrB->zdate, AlphaB,
									 lNumFx, datesFx, dGamma,
									 &lNumBFx, &dDatesBFx, &dAlphaBFx, &dGammaFxB))
		return error;

	for (i = 2; i < lNumEx; i++)
	{
		dDeltaLambda = dLambda[i] - dLambda[i-1];
		dIntGA += LGMIntAlphaUAlphaV(lNumAFx, dDatesAFx, dGammaFxA, dAlphaAFx,
									 datesEx[i-1], datesEx[i]);

		dIntGB += LGMIntAlphaUAlphaV(lNumBFx, dDatesBFx, dGammaFxB, dAlphaBFx,
									 datesEx[i-1], datesEx[i]);

		dIntG2 += LGMIntAlpha(lNumFx, datesFx, dGamma, 
							 datesEx[i-1], datesEx[i]);

		dksiB = LGMZetaFromTS(datesEx[i], 0, tsPtrB);

		m *= (1 + dDeltaLambda * dCorrelBFx * dIntGB
			   + dLambda[i] * dDeltaLambda * dksiB) / 
			(dIntG2 + 2 * dLambda[i] * dCorrelBFx * dIntGB
			   + dLambda[i] * dLambda[i] * dksiB);

		(*dSki12)[i] = (dCorrelAFx * dIntGA + dLambda[i] * dCorrelAB * dIntAB) / m;

		(*dSki22)[i] = (dIntG2 + 2.0 * dLambda[i] * dCorrelBFx * dIntGB
					+ dLambda[i] * dLambda[i] * dksiB) / m / m;
	}
	return NULL;
}
