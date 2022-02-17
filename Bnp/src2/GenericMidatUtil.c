
#include "GenericMidatUtil.h"
#include "math.h"
#include "opfnctns.h"

Err	genmidat_alloc_model(	int				iNbExe,
							int				iNbFactor,
							GENMIDAT_MODEL	sModel)
{
Err	err = NULL;

	sModel->iNbFactor = iNbFactor;

	/* Memory allocation */
	sModel->iNbVols = iNbExe;
	sModel->dForward = calloc(sModel->iNbVols, sizeof(double));
	sModel->dTimes = calloc(sModel->iNbVols, sizeof(double));
	sModel->dSigma = calloc(sModel->iNbVols, sizeof(double));
	sModel->dBeta = calloc(sModel->iNbVols, sizeof(double));
	sModel->dLambda = calloc(sModel->iNbVols, sizeof(double));

	if (!sModel->dTimes || !sModel->dTimes || !sModel->dSigma || !sModel->dBeta || !sModel->dLambda)
	{
		err = "Memory allocation faillure in genmidat_alloc_model";
		return err;
	}

	if (iNbFactor == 2)
	{		
		sModel->dSigma2 = calloc(sModel->iNbVols, sizeof(double));
		sModel->dBeta2 = calloc(sModel->iNbVols, sizeof(double));
		sModel->dLambda2 = calloc(sModel->iNbVols, sizeof(double));
		sModel->dCorrel = calloc(sModel->iNbVols, sizeof(double));

		if (!sModel->dSigma2 || !sModel->dBeta2 || !sModel->dLambda2 || !sModel->dCorrel)
		{
			err = "Memory allocation faillure (2) in genmidat_alloc_model";
			return err;
		}
	}

	/* Set the default */
	sModel->dNumeraire = 1.0;
	sModel->dStartBeta = 2.0;
	sModel->dStartBeta2 = 1.8;

	return err;
}

void genmidat_init_NULL_model(GENMIDAT_MODEL	sModel)
{
	sModel->iNbVols = 0;
	sModel->dTimes = NULL;
	sModel->dForward = NULL;
	sModel->dBeta = NULL;
	sModel->dLambda = NULL;
	sModel->dSigma = NULL;
	sModel->dLambda2 = NULL;
	sModel->dBeta2 = NULL;
	sModel->dSigma2 = NULL;
	sModel->dCorrel = NULL;
}

Err	genmidat_init_model(double			dAlpha,
						double			dGamma,
						double			dRho,
						double			*dTimes,
						double			*dForward,
						double			*dSigma,						
						double			*dBeta,
						double			*dSigma2,
						double			*dCorrel,
						double			dStartBeta2,
						double			dNumeraire,
						GENMIDAT_MODEL	sModel)
{
Err	err = NULL;

	sModel->dAlpha = dAlpha;
	sModel->dGamma = dGamma;
	sModel->dRho = dRho;
	sModel->dAlpha2 = sModel->dAlpha * sModel->dAlpha;
	sModel->dRhoAlpha = sModel->dRho * sModel->dAlpha;

	sModel->dNumeraire = dNumeraire;

	if (dTimes)
	{
		memcpy(sModel->dTimes, dTimes, sModel->iNbVols * sizeof(double));		
	}

	if (dForward)
	{
		memcpy(sModel->dForward, dForward, sModel->iNbVols * sizeof(double));		
	}

	if (dSigma)
	{
		memcpy(sModel->dSigma, dSigma, sModel->iNbVols * sizeof(double));		
	}

	if (dBeta)
	{
		memcpy(sModel->dBeta, dBeta, sModel->iNbVols * sizeof(double));

		genmidat_convert_beta_into_lambda(sModel);
	}

	if (sModel->iNbFactor == 2)
	{
		if (dSigma2)
		{
			memcpy(sModel->dSigma2, dSigma2, sModel->iNbVols * sizeof(double));		
		}

		if (dCorrel)
		{
			memcpy(sModel->dCorrel, dCorrel, sModel->iNbVols * sizeof(double));
		}

		sModel->dStartBeta2 = dStartBeta2;
	}

	if (sModel->iNbFactor == 2 && dSigma2 && dCorrel && dBeta)
	{
		genmidat_fill_2factor(sModel);
	}
	
	return NULL;
}

void genmidat_free_model(GENMIDAT_MODEL	sModel)
{
	if (sModel)
	{
		if (sModel->dTimes) free (sModel->dTimes);
		if (sModel->dForward) free (sModel->dForward);
		if (sModel->dSigma) free (sModel->dSigma);
		if (sModel->dBeta) free (sModel->dBeta);
		if (sModel->dLambda) free (sModel->dLambda);

		if (sModel->dBeta2) free (sModel->dBeta2);
		if (sModel->dLambda2) free (sModel->dLambda2);
		if (sModel->dSigma2) free (sModel->dSigma2);
		if (sModel->dCorrel) free (sModel->dCorrel);
	}
}

Err	genmidat_copy_model(GENMIDAT_MODEL	sSourceModel,
						GENMIDAT_MODEL	sDestModel)
{
Err	err = NULL;

	err = genmidat_alloc_model(	sSourceModel->iNbVols,
								sSourceModel->iNbFactor,
								sDestModel);

	if (err) return err;

	err = genmidat_init_model(
							sSourceModel->dAlpha,
							sSourceModel->dGamma,
							sSourceModel->dRho,
							sSourceModel->dTimes,
							sSourceModel->dForward,
							sSourceModel->dSigma,
							sSourceModel->dBeta,
							sSourceModel->dSigma2,
							sSourceModel->dCorrel,
							sSourceModel->dStartBeta2,
							sSourceModel->dNumeraire,							
							sDestModel);

	if (err) return err;

	return err;
}

void genmidat_convert_beta_into_lambda(GENMIDAT_MODEL sModel)
{
int		i;
double	dSumLam;

	if (0)
	{
		if (0)
		{
			for (i=0; i<sModel->iNbVols; i++)
			{
				if (i > 0)
				{
					sModel->dLambda[i] = -log(sModel->dBeta[i] / sModel->dBeta[i-1]) / (sModel->dTimes[i] - sModel->dTimes[i-1]);
				}
				else
				{
					sModel->dLambda[i] = -log(sModel->dBeta[i]) / sModel->dTimes[0];			
				}
				
				if (sModel->iNbFactor == 2)
				{
					/* Update Beta2 and Lambda2 */
					sModel->dLambda2[i] = sModel->dLambda[i] + sModel->dGamma;
					sModel->dBeta2[i] = sModel->dBeta[i] * exp(-sModel->dGamma * sModel->dTimes[i]);
				}
			}
		}
		else
		{
			
			dSumLam = 0.0;

			for (i=0; i<sModel->iNbVols; i++)
			{
				if (i > 0)
				{
					if (sModel->dBeta[i-1] - sModel->dBeta[i] > 1.0E-08)
					{
						sModel->dLambda[i] = (-log(sModel->dBeta[i-1] - sModel->dBeta[i]) -dSumLam) / (sModel->dTimes[i] - sModel->dTimes[i-1]);
					}
					else
					{
						sModel->dLambda[i] = 100.0;
					}

					dSumLam += sModel->dLambda[i] * (sModel->dTimes[i] - sModel->dTimes[i-1]);
				}
				else
				{
					if (sModel->dBeta[i] >= 0.9999999999)
					{
						sModel->dLambda[i] = 0.0;
					}
					else
					{
						sModel->dLambda[i] = -log(1.0 - sModel->dBeta[i]) / sModel->dTimes[i];
					}

					dSumLam += sModel->dLambda[i] * sModel->dTimes[i];
				}
			}
		}
	}
}

void genmidat_convert_lambda_into_beta(GENMIDAT_MODEL sModel)
{
int		i;
double	sum_lam;

	sum_lam = 0.0;

	for (i=0; i<sModel->iNbVols; i++)
	{
		if (i > 0)
		{
			sum_lam += sModel->dLambda[i] * (sModel->dTimes[i] - sModel->dTimes[i-1]);
		}
		else
		{
			sum_lam += sModel->dLambda[i] * sModel->dTimes[i];	
		}

		sModel->dBeta[i] = exp(-sum_lam);

		if (sModel->iNbFactor == 2)
		{
			/* Update Beta2 and Lambda2 */
			sModel->dLambda2[i] = sModel->dLambda[i] + sModel->dGamma;
			sModel->dBeta2[i] = sModel->dBeta[i] * exp(-sModel->dGamma * sModel->dTimes[i]);
		}
	}
}

Err genmidat_fill_2factor(GENMIDAT_MODEL sModel)
{
Err		err = NULL;
int		i;
double	exp_fact;
	
	sModel->dBeta2[0] = sModel->dStartBeta2 + (sModel->dBeta[0] - sModel->dStartBeta) * exp(-sModel->dGamma * sModel->dTimes[0]);
	sModel->dLambda2[0] = 0.0;

	for (i=1; i<sModel->iNbVols; i++)
	{
		exp_fact = exp(sModel->dGamma * sModel->dTimes[i]);
		sModel->dBeta2[i] = sModel->dBeta2[i-1] + (sModel->dBeta[i] - sModel->dBeta[i-1]) / exp_fact;
		sModel->dLambda2[i] = sModel->dLambda[i] + sModel->dGamma;		
	}

	return err;
}

double	genmidat_Swap_correlation(	double	dTime1,
								   double	dTime2,
								   double	dAlpha,
								   double	dGamma,
								   double	dRho)
{
double	dVol1, dVol2;
double	rho;

	dVol1 = dAlpha * exp(-dGamma * dTime1);
	dVol2 = dAlpha * exp(-dGamma * dTime2);

	rho = 1.0 + dRho * (dVol1 + dVol2) + dVol1 * dVol2;

	dVol1 = 1.0 + 2.0 * dRho * dVol1 + dVol1 * dVol1;
	dVol2 = 1.0 + 2.0 * dRho * dVol2 + dVol2 * dVol2;

	rho /= sqrt(dVol1 * dVol2);

	return rho;
}

Err	genmidat_free_und_struct(SrtUndPtr pUndDesc)
{	
genmidat_model	*sModel = NULL;
	
	sModel = (genmidat_model *) (pUndDesc->spec_desc);
	genmidat_free_model(sModel);	
	free(pUndDesc->spec_desc);
	free(pUndDesc);
	pUndDesc = NULL;
	return NULL;
}

Err SrtInitGenericMidatUnd(	char	*undName,	/* und name */
							int		iNbFactor,
							double	dNumeraire,
							double	dAlpha,
							double	dGamma,
							double	dRho,
							int		iNbVols,
							double	*dTimes,
							double	*dForward,
							double	*dSigma,
							double	*dBeta,
							double	*dSigma2,
							double	*dCorrel,
							double	dStartBeta2)
{	
	SrtUndPtr		pUndPtr = NULL;
	SrtUndListPtr	und_list;

	genmidat_model	*sModel = NULL;

	Err				err = NULL;		
											
	/* Allocate the structure */
	sModel = calloc(1, sizeof(genmidat_model));
	pUndPtr = (SrtUndPtr) calloc(1, sizeof(SrtUndDesc));

	if (!sModel || !pUndPtr)
	{
		err = "Memory allocation faillure in SrtInitGenericMidatUnd";
		goto FREE_RETURN;
	}
				
	pUndPtr->underl_type = GENMIDAT_UND;
	strcpy(pUndPtr->underl_name, undName);
	strupper(pUndPtr ->underl_name);
	strip_white_space(pUndPtr->underl_name);
	strcpy(pUndPtr->underl_lbl, "GENMIDAT_UND");			
	pUndPtr->underl_ccy = "USD";
	
	/* Fill the model */
	err = genmidat_alloc_model(	iNbVols,
								iNbFactor,
								sModel);

	if (err) goto FREE_RETURN;
	
	err = genmidat_init_model(	dAlpha,
								dGamma,
								dRho,
								dTimes,
								dForward,
								dSigma,
								dBeta,
								dSigma2,
								dCorrel,
								dStartBeta2,
								dNumeraire,
								sModel);

	if (err) goto FREE_RETURN;
	
	pUndPtr->spec_desc = sModel;

	/* Put the underlying into the depot */

	und_list = get_underlying_list();

	err = srt_f_lstins(		und_list,
							pUndPtr ->underl_name,
							0.0,
							OBJ_PTR_UND,
							(void*)pUndPtr,
							&genmidat_free_und_struct,
							&(pUndPtr ->underl_ticker));

FREE_RETURN:

	if (err)
	{
		genmidat_free_model(sModel);
		if (sModel) free (sModel);
		if (pUndPtr) free (pUndPtr);
	}

	return err;
}

Err	genmidat_get_model(	char			*cUndName,
						GENMIDAT_MODEL	sModel)
{
	SrtUndPtr		und;
	GENMIDAT_MODEL	sSavedModel;
	Err				err = NULL;

	und = lookup_und (cUndName);

	if (!und)
	{
		err = serror("Couldn't fin underlying named %s", cUndName);
		goto FREE_RETURN;
	}

	if (get_underlying_type (und) != GENMIDAT_UND)
	{
		return serror("Underlying %s is not of type GENMIDAT", cUndName);
		goto FREE_RETURN;
	}

	sSavedModel = (GENMIDAT_MODEL) (und->spec_desc);

	err = genmidat_copy_model(	sSavedModel,
								sModel);

	if (err) goto FREE_RETURN;

FREE_RETURN:

	return err;
}

Err	genmidat_implied_volatility_and_price_model(	int				iIntegStartIndex,
													int				iIntegEndIndex,
													int				iUndIndex,
													int				iNbPeriod,
													GENMIDAT_MODEL	sModel,
													double			*dVolatility,
													double			*dPrice)
{
Err	err = NULL;
int		i;
double	varX, varY, covar;
double	dStartTime, dEndTime;
double	dVol;
double	dBeta1, dBeta2;
int		iUndIndex2;
double	dForward;

	/* Checks */
	if (iIntegStartIndex > sModel->iNbVols || iIntegEndIndex > sModel->iNbVols
		|| iIntegStartIndex >= iIntegEndIndex || iUndIndex > sModel->iNbVols - 1
		|| iUndIndex < iIntegEndIndex)
	{
		err = "Check Indexes in genmidat_implied_volatility_model";
		return err;
	}

	/* Special case for last one */
	if (iUndIndex + iNbPeriod == sModel->iNbVols)
	{
		iNbPeriod = 0;
	}		

	varX = 0.0;
	varY = 0.0;
	covar = 0.0;

	/* Get the total variances */
	if (iIntegStartIndex < 0)
	{
		/* first period */
		varX += sModel->dSigma[0] * sModel->dSigma[0] * sModel->dTimes[0];

		if (sModel->iNbFactor == 2)
		{
			varY += sModel->dSigma2[0] * sModel->dSigma2[0] * sModel->dTimes[0];
			covar += sModel->dCorrel[0] * sModel->dSigma[0] * sModel->dSigma2[0] * sModel->dTimes[0];
		}

		dStartTime = 0.0;
		iIntegStartIndex = 0;
	}
	else
	{
		dStartTime = sModel->dTimes[iIntegStartIndex];
	}	

	for (i=iIntegStartIndex; i<iIntegEndIndex; i++)
	{
		varX += sModel->dSigma[i+1] * sModel->dSigma[i+1] * (sModel->dTimes[i+1] - sModel->dTimes[i]);

		if (sModel->iNbFactor == 2)
		{
			varY += sModel->dSigma2[i+1] * sModel->dSigma2[i+1] * (sModel->dTimes[i+1] - sModel->dTimes[i]);
			covar += sModel->dCorrel[i+1] * sModel->dSigma[i+1] * sModel->dSigma2[i+1] * (sModel->dTimes[i+1] - sModel->dTimes[i]);
		}

	}

	dEndTime = sModel->dTimes[iIntegEndIndex];
		
	if (iNbPeriod == 0)
	{		
		dBeta1 = sModel->dBeta[iUndIndex];
		
		if (sModel->iNbFactor == 1)
		{	
			dBeta2 = 0.0;
		}
		else
		{
			dBeta2 = sModel->dBeta2[iUndIndex];
		}

		dForward = sModel->dForward[iUndIndex];
	}
	else
	{
		iUndIndex2 = iUndIndex + iNbPeriod;

		if (iUndIndex2 > sModel->iNbVols - 1)
		{			
			err = "Invalid Index in genmidat_implied_volatility_model";
			return err;
		}

		dBeta1 = sModel->dBeta[iUndIndex] - sModel->dBeta[iUndIndex2];

		if (sModel->iNbFactor == 1)
		{	
			dBeta2 = 0.0;
		}
		else
		{
			dBeta2 = sModel->dBeta2[iUndIndex] - sModel->dBeta2[iUndIndex2];
		}

		dForward = sModel->dForward[iUndIndex] - sModel->dForward[iUndIndex2];
	}

	dVol = dBeta1 * dBeta1 * varX + dBeta2 * dBeta2 * varY + 2.0 * dBeta1 * dBeta2 * covar;
		
	if (dVol > 0.0)
	{
		dVol = sqrt(dVol / (dEndTime - dStartTime));
	}
	else
	{
		err = "Negative variance in genmidat_implied_volatility";
		return err;
	}

	*dVolatility = dVol;

	/* Calculation of the price */
	*dPrice = srt_f_optblknrm(	dForward,
								0.0,
								*dVolatility,
								dEndTime - dStartTime,
								sModel->dNumeraire,
								SRT_CALL,
								PREMIUM);

	return err;
}

Err	genmidat_implied_volatility_and_price_und(	char			*cUndName,
												int				iIntegStartIndex,
												int				iIntegEndIndex,
												int				iUndIndex,
												int				iNbPeriod,
												double			*dVolatility,
												double			*dPrice)
		{
GENMIDAT_MODEL	sModel;
Err				err = NULL;

	/* Allocation */
	sModel = calloc(1, sizeof(genmidat_model));

	if (!sModel)
	{
		err = "Memory allocation faillure in genmidat_implied_volatility_und";
		goto FREE_RETURN;
	}

	/* Get the Model */
	err = genmidat_get_model(	cUndName,
								sModel);

	if (err) goto FREE_RETURN;

	/* Defaults */
	if (iIntegEndIndex == -1)
	{
		iIntegEndIndex = iUndIndex;
	}	

	err = genmidat_implied_volatility_and_price_model(
											iIntegStartIndex,
											iIntegEndIndex,
											iUndIndex,
											iNbPeriod,
											sModel,
											dVolatility,
											dPrice);

	if (err) goto FREE_RETURN;

FREE_RETURN:

	if (sModel)
	{
		genmidat_free_model(sModel);
		free(sModel);
	}

	return err;
}

Err	genmidat_implied_correlation_model(	int				iIntegStartIndex,
										int				iIntegEndIndex,
										int				iUndIndex1,
										int				iNbPeriod1,
										int				iUndIndex2,
										int				iNbPeriod2,
										GENMIDAT_MODEL	sModel,
										double			*dCorrelation)
{
Err	err = NULL;
int		i;
double	varX, varY, covar;
double	dStartTime, dEndTime;
double	dVol1, dVol2, dCovar;
double	dBetaX1, dBetaY1, dBetaX2, dBetaY2;
int		iUndIndex12, iUndIndex22;

	/* Checks */
	if (iIntegStartIndex > sModel->iNbVols || iIntegEndIndex > sModel->iNbVols || iIntegStartIndex >= iIntegEndIndex 
		|| iUndIndex1 > sModel->iNbVols - 1 || iUndIndex1 < iIntegEndIndex || iUndIndex2 > sModel->iNbVols - 1 || iUndIndex2 < iIntegEndIndex)
	{
		err = "Check Indexes in genmidat_implied_correlation";
		return err;
	}

	/* Special case for last one */
	if (iUndIndex1 + iNbPeriod1 == sModel->iNbVols)
	{
		iNbPeriod1 = 0;
	}

	if (iUndIndex2 + iNbPeriod2 == sModel->iNbVols)
	{
		iNbPeriod2 = 0;
	}

	varX = 0.0;
	varY = 0.0;
	covar = 0.0;

	/* Get the total variances */
	if (iIntegStartIndex < 0)
	{
		/* first period */
		varX += sModel->dSigma[0] * sModel->dSigma[0] * sModel->dTimes[0];

		if (sModel->iNbFactor == 2)
		{
			varY += sModel->dSigma2[0] * sModel->dSigma2[0] * sModel->dTimes[0];
			covar += sModel->dCorrel[0] * sModel->dSigma[0] * sModel->dSigma2[0] * sModel->dTimes[0];
		}

		dStartTime = 0.0;
		iIntegStartIndex = 0;
	}
	else
	{
		dStartTime = sModel->dTimes[iIntegStartIndex];
	}	

	for (i=iIntegStartIndex; i<iIntegEndIndex; i++)
	{
		varX += sModel->dSigma[i+1] * sModel->dSigma[i+1] * (sModel->dTimes[i+1] - sModel->dTimes[i]);

		if (sModel->iNbFactor == 2)
		{
			varY += sModel->dSigma2[i+1] * sModel->dSigma2[i+1] * (sModel->dTimes[i+1] - sModel->dTimes[i]);
			covar += sModel->dCorrel[i+1] * sModel->dSigma[i+1] * sModel->dSigma2[i+1] * (sModel->dTimes[i+1] - sModel->dTimes[i]);
		}

	}

	dEndTime = sModel->dTimes[iIntegEndIndex];
		
	if (iNbPeriod1 == 0)
	{		
		dBetaX1 = sModel->dBeta[iUndIndex1];
		
		if (sModel->iNbFactor == 1)
		{	
			dBetaY1 = 0.0;
		}
		else
		{
			dBetaY1 = sModel->dBeta2[iUndIndex1];
		}		
	}
	else
	{
		iUndIndex12 = iUndIndex1 + iNbPeriod1;

		if (iUndIndex12 > sModel->iNbVols - 1)
		{			
			err = "Invalid Index in genmidat_implied_volatility_model";
			return err;
		}

		dBetaX1 = sModel->dBeta[iUndIndex1] - sModel->dBeta[iUndIndex12];

		if (sModel->iNbFactor == 1)
		{	
			dBetaY1 = 0.0;
		}
		else
		{
			dBetaY1 = sModel->dBeta2[iUndIndex1] - sModel->dBeta2[iUndIndex12];
		}
	}

	if (iNbPeriod2 == 0)
	{		
		dBetaX2 = sModel->dBeta[iUndIndex2];
		
		if (sModel->iNbFactor == 1)
		{	
			dBetaY2 = 0.0;
		}
		else
		{
			dBetaY2 = sModel->dBeta2[iUndIndex2];
		}		
	}
	else
	{
		iUndIndex22 = iUndIndex2 + iNbPeriod2;

		if (iUndIndex22 > sModel->iNbVols - 1)
		{			
			err = "Invalid Index in genmidat_implied_volatility_model";
			return err;
		}

		dBetaX2 = sModel->dBeta[iUndIndex2] - sModel->dBeta[iUndIndex22];

		if (sModel->iNbFactor == 1)
		{	
			dBetaY2 = 0.0;
		}
		else
		{
			dBetaY2 = sModel->dBeta2[iUndIndex2] - sModel->dBeta2[iUndIndex22];
		}
	}

	dVol1 = dBetaX1 * dBetaX1 * varX + dBetaY1 * dBetaY1 * varY + 2.0 * dBetaX1 * dBetaY1 * covar;
	dVol2 = dBetaX2 * dBetaX2 * varX + dBetaY2 * dBetaY2 * varY + 2.0 * dBetaX2 * dBetaY2 * covar;
	dCovar = dBetaX1 * dBetaX2 * varX 
			+ (dBetaX1 * dBetaY2 + dBetaY1 * dBetaX2) * covar
			+ dBetaY1 * dBetaY2 * varY;			

	*dCorrelation = dCovar / sqrt(dVol1 * dVol2);	

	return err;
}

Err	genmidat_implied_correlation_und(	char			*cUndName,
									int				iIntegStartIndex,
									int				iIntegEndIndex,
									int				iUndIndex1,
									int				iNbPeriod1,
									int				iUndIndex2,
									int				iNbPeriod2,
									double			*dCorrelation)
{
GENMIDAT_MODEL	sModel;
Err				err = NULL;

	/* Allocation */
	sModel = calloc(1, sizeof(genmidat_model));

	if (!sModel)
	{
		err = "Memory allocation faillure in genmidat_implied_volatility_und";
		goto FREE_RETURN;
	}

	/* Get the Model */
	err = genmidat_get_model(	cUndName,
								sModel);

	if (err) goto FREE_RETURN;

	/* Defaults */
	if (iIntegEndIndex == -1)
	{
		iIntegEndIndex = min(iUndIndex1, iUndIndex2);
	}	

	err = genmidat_implied_correlation_model(	iIntegStartIndex,
												iIntegEndIndex,
												iUndIndex1,
												iNbPeriod1,
												iUndIndex2,
												iNbPeriod2,
												sModel,
												dCorrelation);

	if (err) goto FREE_RETURN;

FREE_RETURN:

	if (sModel)
	{
		genmidat_free_model(sModel);
		free(sModel);
	}

	return err;
}