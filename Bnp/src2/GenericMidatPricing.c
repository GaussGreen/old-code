
#include "math.h"
#include "opfnctns.h"
#include "GenericMidatPricing.h"
#include "Fx3FUtils.h"


void genmidat_set_pdeparams_default(GENMIDAT_PDEPAMS sParams)
{
	sParams->iPricingMethod = 0;

	sParams->iDiscType = 0;
	
	sParams->lNbTime = 75;
	sParams->lNbX = 101;

	sParams->dTheta = 0.5;
	sParams->dNbStd = 5.0;

	sParams->iJumpQuadra = 1;

	sParams->iDiffMethod = 0;
	sParams->iAdaptGrid = 0;

	sParams->lNbPaths = 20000;	
}

void GenericMidat_OneFactor_ChangeBetaGrid(	long				lStart,
											long				lEnd,
											long				lStartProd,
											long				lEndProd,
																
											double				*dOldGridX,
											double				dOldBeta,
											double				*dNewGridX,
											double				dNewBeta,
											
											GENMIDAT_MODEL		sModel,
											GENMIDAT_PDEPAMS	sParams,

											double				**dOldValues,
											double				**dNewValues)
{
Err	err = NULL;
double	dOldX;
int		oldI;
double	dFirstX, dx, dx2;
double	dRatio, dProba, dProbaDown, dProbaMid, dProbaUp;
int		i, j;

	dRatio = dOldBeta / dNewBeta;
	dx = (dOldGridX[lStart + 1] - dOldGridX[lStart]);
	dx2 = dx * dx;
	dFirstX = dOldGridX[lStart];

	if (sParams->iJumpQuadra == 0)
	{
		/* Linear Interpolation */
		for (i=lStart; i<=lEnd; i++)
		{
			dOldX = dRatio * dNewGridX[i];
			oldI = max(min(lStart + (int) ((dOldX - dFirstX) / dx), lEnd - 1), lStart);
			dProba = (dOldX - dOldGridX[oldI]) / (dOldGridX[oldI+1] - dOldGridX[oldI]);

			for (j=lStartProd; j<=lEndProd; j++)
			{
				dNewValues[i][j] = (1.0 - dProba) * dOldValues[oldI][j] + dProba * dOldValues[oldI+1][j];
			}			
		}
	}
	else
	{
		/* Quadratic Interpolation */
		for (i=lStart; i<=lEnd; i++)
		{
			dOldX = dRatio * dNewGridX[i];
			oldI = max(min(lStart + (int) ((dOldX - dFirstX) / dx + 0.5), lEnd - 1), lStart + 1);

			dProbaDown = (dOldX * (dOldX - dOldGridX[oldI] - dOldGridX[oldI+1]) + dOldGridX[oldI] * dOldGridX[oldI+1]) / (2.0 * dx2);
			dProbaMid = (dOldX * (dOldX - dOldGridX[oldI-1] - dOldGridX[oldI+1]) + dOldGridX[oldI-1] * dOldGridX[oldI+1]) / (-dx2);
			dProbaUp = (dOldX * (dOldX - dOldGridX[oldI-1] - dOldGridX[oldI]) + dOldGridX[oldI-1] * dOldGridX[oldI]) / (2.0 * dx2);

			for (j=lStartProd; j<=lEndProd; j++)
			{
				dNewValues[i][j] = dProbaDown * dOldValues[oldI-1][j]
									+ dProbaMid * dOldValues[oldI][j]
									+ dProbaUp * dOldValues[oldI+1][j];
			}
		}

	}
}

void GenericMidat_TwoFactor_ChangeBetaGrid(	/* Old Grids */
											long				lOldStartX,
											long				lOldEndX,
											double				*dOldGridX,

											long				lOldStartZ,
											long				lOldEndZ,
											double				*dOldGridZ,
											
											/* New Grids */
											long				lNewStartX,
											long				lNewEndX,
											double				*dNewGridX,

											long				lNewStartZ,
											long				lNewEndZ,
											double				*dNewGridZ,

											/* Model parameters */
											int					iIndex,
											GENMIDAT_MODEL		sModel,
											GENMIDAT_PDEPAMS	sParams,

											/* Product Values */
											long				lStartProd,
											long				lEndProd,
											double				***dTempValues,
											double				***dOldValues,
											double				***dNewValues)
{
Err	err = NULL;
double	dOldX, dOldZ;
int		oldI, oldZ;
double	dFirstX, dFirstZ, dx, dx2, dz, dz2;
double	dRatio1, dRatio2, dCoefX, dCoefZ, dProba, dProbaDown, dProbaMid, dProbaUp;
int		i, j, k;

	dRatio1 = sModel->dBeta[iIndex] / sModel->dBeta[iIndex-1];
	dRatio2 = sModel->dBeta2[iIndex] / sModel->dBeta2[iIndex-1];

	dCoefX = 1.0 / sqrt(1.0 - sModel->dCorrel[iIndex] * sModel->dCorrel[iIndex])
			* sModel->dBeta2[iIndex] / sModel->dBeta[iIndex-1]
		* (	sModel->dCorrel[iIndex-1] * sModel->dSigma2[iIndex-1] / sModel->dSigma[iIndex-1]
			- sModel->dCorrel[iIndex] * sModel->dSigma2[iIndex] / sModel->dSigma[iIndex]);

	dCoefZ = dRatio2 * sqrt(1.0 - sModel->dCorrel[iIndex-1] * sModel->dCorrel[iIndex-1])
		/ sqrt(1.0 - sModel->dCorrel[iIndex] * sModel->dCorrel[iIndex]);

	dx = (dOldGridX[lOldStartX + 1] - dOldGridX[lOldStartX]);
	dx2 = dx * dx;
	dFirstX = dOldGridX[lOldStartX];

	dz = (dOldGridZ[lOldStartZ + 1] - dOldGridZ[lOldStartZ]);
	dz2 = dz * dz;
	dFirstZ = dOldGridZ[lOldStartZ];

	if (sParams->iJumpQuadra == 0)
	{
		/* Linear Interpolation */

		/* First in X */
		for (i=lNewStartX; i<=lNewEndX; i++)
		{
			dOldX = dRatio1 * dNewGridX[i];
			oldI = max(min(lOldStartX + (int) ((dOldX - dFirstX) / dx), lOldEndX - 1), lOldStartX);

			dProba = (dOldX - dOldGridX[oldI]) / (dOldGridX[oldI+1] - dOldGridX[oldI]);

			for (j=lNewStartZ; j<=lNewEndZ; j++)
			{
				for (k=lStartProd; k<=lEndProd; k++)
				{
					dTempValues[i][j][k] = (1.0 - dProba) * dOldValues[oldI][j][k]
										+ dProba * dOldValues[oldI+1][j][k];
				}
			}
		}

		/* Then in Z */
		for (i=lNewStartX; i<=lNewEndX; i++)
		{
			for (j=lNewStartZ; j<=lNewEndZ; j++)
			{
				dOldZ = dCoefX * dNewGridX[i] + dCoefZ * dNewGridZ[j];
				oldZ = max(min(lOldStartZ + (int) ((dOldZ - dFirstZ) / dz), lOldEndZ - 1), lOldStartZ);

				dProba = (dOldZ - dOldGridZ[oldZ]) / (dOldGridZ[oldZ+1] - dOldGridZ[oldZ]);
				
				for (k=lStartProd; k<=lEndProd; k++)
				{
					dNewValues[i][j][k] = (1.0 - dProba) * dTempValues[i][oldZ][k]
										+ dProba * dTempValues[i][oldZ+1][k];
				}
			}
		}
	}
	else
	{
		/* Quadratic Interpolation */

		/* First in X */
		for (i=lNewStartX; i<=lNewEndX; i++)
		{
			dOldX = dRatio1 * dNewGridX[i];
			oldI = max(min(lOldStartX + (int) ((dOldX - dFirstX) / dx + 0.5), lOldEndX - 1), lOldStartX + 1);

			dProbaDown = (dOldX * (dOldX - dOldGridX[oldI] - dOldGridX[oldI+1]) + dOldGridX[oldI] * dOldGridX[oldI+1]) / (2.0 * dx2);
			dProbaMid = (dOldX * (dOldX - dOldGridX[oldI-1] - dOldGridX[oldI+1]) + dOldGridX[oldI-1] * dOldGridX[oldI+1]) / (-dx2);
			dProbaUp = (dOldX * (dOldX - dOldGridX[oldI-1] - dOldGridX[oldI]) + dOldGridX[oldI-1] * dOldGridX[oldI]) / (2.0 * dx2);

			for (j=lNewStartZ; j<=lNewEndZ; j++)
			{
				for (k=lStartProd; k<=lEndProd; k++)
				{
					dTempValues[i][j][k] = dProbaDown * dOldValues[oldI-1][j][k]
										+ dProbaMid * dOldValues[oldI][j][k]
										+ dProbaUp * dOldValues[oldI+1][j][k];
				}
			}
		}

		/* Then in Z */
		for (i=lNewStartX; i<=lNewEndX; i++)
		{
			for (j=lNewStartZ; j<=lNewEndZ; j++)
			{
				dOldZ = dCoefX * dNewGridX[i] + dCoefZ * dNewGridZ[j];
				oldZ = max(min(lOldStartZ + (int) ((dOldZ - dFirstZ) / dz + 0.5), lOldEndZ - 1), lOldStartZ + 1);

				dProbaDown = (dOldZ * (dOldZ - dOldGridZ[oldZ] - dOldGridZ[oldZ+1]) + dOldGridZ[oldZ] * dOldGridZ[oldZ+1]) / (2.0 * dz2);
				dProbaMid = (dOldZ * (dOldZ - dOldGridZ[oldZ-1] - dOldGridZ[oldZ+1]) + dOldGridZ[oldZ-1] * dOldGridZ[oldZ+1]) / (-dz2);
				dProbaUp = (dOldZ * (dOldZ - dOldGridZ[oldZ-1] - dOldGridZ[oldZ]) + dOldGridZ[oldZ-1] * dOldGridZ[oldZ]) / (2.0 * dz2);

				for (k=lStartProd; k<=lEndProd; k++)
				{
					dNewValues[i][j][k] = dProbaDown * dTempValues[i][oldZ-1][k]
										+ dProbaMid * dTempValues[i][oldZ][k]
										+ dProbaUp * dTempValues[i][oldZ+1][k];
				}
			}
		}
	}
}

Err	GenericMidat_ConvertModelLinear(	int				iNbVols,
										double			*dTimes,
										double			*dSigmaTarget,
										double			dSigma0,
										double			*dCoefs,										
										double			*dSigma)
{
int		i;
double	dLastTime, dLastSigma, dMinVol;
double	a, b, c, delta;
Err		err = NULL;

	dLastTime = 0.0;
	
	dLastSigma = dSigma0;
	dLastTime = 0.0;
	dMinVol = 0.0001;

	a = 1.0 / 3.0;

	for (i=0; i<iNbVols; i++)
	{
		b = dLastSigma;
		c = dLastSigma * dLastSigma - dSigmaTarget[i] * dSigmaTarget[i];
		delta = b * b - 4.0 * a * c;

		if (delta < 0.0)
		{
			/* Put the best we can do: */
			dCoefs[i] = (dMinVol - dLastSigma) / (dTimes[i] - dLastTime);			
		}
		else
		{
			dCoefs[i] = (-b + sqrt(delta)) / (2.0 * a * (dTimes[i] - dLastTime));
			dCoefs[i] = max(dCoefs[i], (dMinVol - dLastSigma) / (dTimes[i] - dLastTime));
		}

		dLastSigma +=  dCoefs[i] * (dTimes[i] - dLastTime);
		dSigma[i] = dLastSigma;
		dLastTime = dTimes[i];
	}

	return err;
}

static Err GenericMidat_VarModelLinear(	int				iNbVols,
										double			*dTimes,
										double			*dSigmaTarget,
										double			dSigma0,
										double			*dCoefs,										
										double			*dSigma,
										double			*dResult)
{
Err	err = NULL;
int	i;

	err = GenericMidat_ConvertModelLinear(	iNbVols,
											dTimes,
											dSigmaTarget,
											dSigma0,
											dCoefs,	
											dSigma);

	if (err) return err;

	*dResult = 0.0;

	for (i=0; i<iNbVols; i++)
	{
		*dResult += dCoefs[i] * dCoefs[i];
	}

	*dResult = sqrt(*dResult);

	return err;
}

Err	GenericMidat_OptimalConvertModelLinear(	int				iNbVols,
											double			*dTimes,
											double			*dSigmaTarget,
											double			*dSigma0,
											double			*dCoefs,										
											double			*dSigma)
{
Err	err = NULL;
int		i;
double	dLastTime;
double	sig1, sig2, sig;
double	var1, var2, var;
double	step;
int		optim;

	dLastTime = 0.0;
	sig1 = 0.0;

	/* find a first guess */
	for (i=0; i<iNbVols; i++)
	{
		sig1 += dSigmaTarget[i] * dSigmaTarget[i] * (dTimes[i] - dLastTime);
		dLastTime = dTimes[i];
	}
	
	sig1 = sqrt(sig1 / dTimes[iNbVols-1]);	
	sig2 = sig1;
	
	err = GenericMidat_VarModelLinear(	iNbVols,
										dTimes,
										dSigmaTarget,
										sig2,
										dCoefs,	
										dSigma,
										&var2);

	i = 1;

	while (err && i < 4)
	{
		sig2 *= 2.0;

		err = GenericMidat_VarModelLinear(	iNbVols,
										dTimes,
										dSigmaTarget,
										sig2,
										dCoefs,	
										dSigma,
										&var2);

		i++;
	}

	if (err)
	{
		/* wrong direction */
		i = 1;
		sig2 = sig1;

		while (err && i < 4)
		{
			sig2 /= 2.0;

			err = GenericMidat_VarModelLinear(	iNbVols,
											dTimes,
											dSigmaTarget,
											sig2,
											dCoefs,	
											dSigma,
											&var2);

			i++;
		}
	}

	if (err) return err;

	/* We have a first point, look for the optimal */
	sig1 = sig2;
	var1 = var2;

	step = sig1 / 20.0;

	/* Find the direction */
	sig2 = sig1 + step;

	err = GenericMidat_VarModelLinear(	iNbVols,
										dTimes,
										dSigmaTarget,
										sig2,
										dCoefs,	
										dSigma,
										&var2);

	if (err) return err;

	if (var2 > var1)
	{
		step *= -1.0;
		sig2 = sig1;
		sig = sig1;
		var = var1;
	}
	else
	{
		sig = sig2;
		var = var2;	
	}

	i = 1;
	optim = 2;

	while (optim && i < 50)
	{
		sig2 += step;

		err = GenericMidat_VarModelLinear(	iNbVols,
										dTimes,
										dSigmaTarget,
										sig2,
										dCoefs,	
										dSigma,
										&var2);

		if (err)
		{
			optim = 0;
			err = NULL;
		}
		else
		{
			if (var2 < var)
			{
				sig = sig2;
				var = var2;
				optim = 2;
			}
			else
			{
				optim --;
			}
		}

		i++;
	}

	*dSigma0 = sig;

	return err;
}

/* solve (exp(lam *t2) -exp(lam*t1)) / lam = target */
static Err solve_integ_exponential(	double	t1,
									 double	t2,
									 double	target,
									 double	*lam)
{
Err	err = NULL;
double	lam1, lam2;
double	price1, price2, delta;
double	lammax;
int		i;

	lammax = 10.0 / t2;

	/* Check for the null lambda */
	if (fabs(target - (t2 - t1)) < 1.0E-08)
	{
		*lam = 0.0;
		return err;
	}

	/* First guess */
	lam1 = (target - (t2 - t1)) / (0.5 * (t2 - t1) * (t2 + t1));
	lam1 = min(lam1, lammax);

	price1 = (exp(lam1 * t2) - exp(lam1 * t1)) / lam1;

	if (price1 > target)
	{
		lam2 = 0.5 * lam1;
	}
	else
	{
		lam2 = 2.0 * lam1;
	}

	i = 1;

	price2 = 1.0E10;

	while (i < 20 && fabs(price2 - target) > 1.0E-10)
	{
		price2 = (exp(lam2 * t2) - exp(lam2 * t1)) / lam2;

		delta = (price2 - price1) / (lam2 - lam1);

		lam1 = lam2;
		price1 = price2;

		lam2 += (target - price2) / delta;
		lam2 = min(lam2, lammax);

		if (fabs(lam2 - lam1) < 1.0E-08)
		{
			break;
		}

		i++;
	}

	if (i == 20)
	{		
		err = "Cannot find a solution in solve_integ_exponential";
	}
	else
	{
		*lam = lam2;
	}	

	return err;
}

Err	GenericMidat_ConvertModelLambda(	int				iNbVols,
										double			*dTimes,
										double			*dSigmaTarget,
										double			*dLambda,
										double			*dSigma0,
										double			*dSigma)
{
int		i;
double	dLastTime, dLastSigma, dTarget;
double	dVar;
Err		err = NULL;
double	dMinVol = 0.0001;
double	dMaxLam = 3.0;

	dLastTime = 0.0;
	dVar = 0.0;

	/* find a good sigma0 */
	for (i=0; i<iNbVols; i++)
	{
		dVar += dSigmaTarget[i] * dSigmaTarget[i] * (dTimes[i] - dLastTime);
		dLastTime = dTimes[i];
	}
	
	*dSigma0 = sqrt(dVar / dTimes[iNbVols-1]);

	err = GenericMidat_OptimalConvertModelLinear(	iNbVols,
													dTimes,
													dSigmaTarget,
													dSigma0,
													dLambda,													
													dSigma);

	if (err) return err;

	dLastSigma = *dSigma0;
	dLastTime = 0.0;

	for (i=0; i<iNbVols; i++)
	{
		dTarget = dSigmaTarget[i] * dSigmaTarget[i] * (dTimes[i] - dLastTime) 
			/ dLastSigma / dLastSigma;

		err = solve_integ_exponential(	0,
										dTimes[i] - dLastTime,
										dTarget,
										&(dLambda[i]));

		if (err)
		{
			/* Apply a floor to avoid negative vols */
			dLambda[i] = log(dMinVol / dLastSigma) / (dTimes[i] - dLastTime);
		}
		else
		{		
			dLambda[i] /= 2.0;
			dLambda[i] = max(dLambda[i], log(dMinVol / dLastSigma) / (dTimes[i] - dLastTime));
		}

		dLambda[i] = max(-dMaxLam, dLambda[i]);
		dLambda[i] = min(dMaxLam, dLambda[i]);

		dLastSigma *= exp(dLambda[i] * (dTimes[i] - dLastTime));
		dLastSigma = max(dLastSigma, 0.0001);
		dSigma[i] = dLastSigma;
		dLastTime = dTimes[i];
	}

	return err;
}

Err	GenericMidat_OneFactor_PDE(	/* Time Discretisation */
								long				lNbTime,
								double				*dTimes,

								/* Model */
								GENMIDAT_MODEL		sModel,

								/*	Product data */
								void				**sFuncParmTab, 
								int					*iEvalEvt,

								/*	Payoff Function */
								Err (*payoff_func)(/* Event */
													double	dTime,

													/* Parameters */
													void	*sFuncParm,																										
													
													/* Gride data	*/
													long			lStartIndex,
													long			lEndIndex,											
													double			*dXt,
													GENMIDAT_RECONS	sReconsParams,
																			
													/* Vector of results to be updated */
													int		iNbProd,
													double	**dProdVal),

								/* PDE Parameters */
								GENMIDAT_PDEPAMS	sParams,

								/* Result */
								int					iNbProd,
								double				*dProdVal)
{
Err	err = NULL;

double	*dXtGrid		= NULL,
		*dNextXtGrid	= NULL,
		*dXtGridTemp,
		*dXtDrift		= NULL,
		*dXtVar			= NULL,
		*dXtVarGrid		= NULL,
		*dShortR		= NULL,
		*dXtTotStd		= NULL,
		*dBeta			= NULL,
		*dLambda		= NULL,
		*dSigma			= NULL;

double	dSigma0;

double	**dNextValues	= NULL,
		**dValues		= NULL,
		**dValuesTemp;

GENMIDAT_RECONS	sReconsParams = NULL;

double	dt, dx;

int		i, j, iIndexMid;
int		step;
int		iFirstEvent;
double	lambda;

CNPDE_TEMP	*pde = NULL;

	
	sParams->lNbX = ((int) (sParams->lNbX / 2.0)) * 2 + 1;

	/* Memory allocations */
	dXtGrid = calloc(sParams->lNbX, sizeof(double));
	dNextXtGrid = calloc(sParams->lNbX, sizeof(double));
	dXtDrift = calloc(sParams->lNbX, sizeof(double));	
	dXtVarGrid = calloc(sParams->lNbX, sizeof(double));
	dShortR = calloc(sParams->lNbX, sizeof(double));

	dXtVar = calloc(lNbTime, sizeof(double));
	dXtTotStd = calloc(lNbTime, sizeof(double));
	dBeta = calloc(lNbTime+1, sizeof(double));

	dNextValues = dmatrix(0, sParams->lNbX - 1, 0, iNbProd - 1);
	dValues = dmatrix(0, sParams->lNbX - 1, 0, iNbProd - 1);

	sReconsParams = calloc(1, sizeof(genmidat_recons));
	
	if (!dXtGrid || !dNextXtGrid || !dXtDrift || !dXtVar || !dXtVarGrid || !dShortR || !dXtTotStd ||
		!dNextValues || !dValues || !dBeta || !sReconsParams)
	{
		err = "Memory allocation faillure in GenericMidat_OneFactor_PDE";
		goto FREE_RETURN;
	}

	/* Fill the dynamics constant */
	sReconsParams->iDiffMethod = sParams->iDiffMethod;
	dXtTotStd[0] = 0.0;

	for (i=1; i<lNbTime; i++)
	{
		dt = dTimes[i] - dTimes[i-1];		
		dXtVar[i] = sModel->dSigma[Get_Index(dTimes[i], sModel->dTimes, sModel->iNbVols)];
		dXtVar[i] *= dXtVar[i] * dt;				
		dXtTotStd[i] = sqrt(dXtTotStd[i-1] * dXtTotStd[i-1] + dXtVar[i]);		
	}

	if (sParams->iDiffMethod == 1)
	{
		/* We diffuse Beta * Xt */
		dBeta[0] = sModel->dBeta[0];
		dXtTotStd[0] = 0.0;
		j = 0;

		for (i=1; i<lNbTime; i++)
		{
			if (j < sModel->iNbVols-1 && dTimes[i] > sModel->dTimes[j])
			{
				j++;
			}

			if (j > 0)
			{
				dBeta[i] = sModel->dBeta[j-1];
			}
			else
			{
				dBeta[i] = sModel->dBeta[0];
			}

			dXtTotStd[i] *= dBeta[i];

			if (!sParams->iAdaptGrid)
			{
				dXtTotStd[i] = max(dXtTotStd[i], dXtTotStd[i-1]);
			}

			dXtVar[i] *= dBeta[i] * dBeta[i];
		}

		if (j < sModel->iNbVols-1 && dTimes[lNbTime-1] > sModel->dTimes[j])
		{
			j++;
		}

		dBeta[lNbTime] = sModel->dBeta[j];
	}

	if (sParams->iDiffMethod == 2)
	{
		/* We diffuse Yt = Xt / sigma(t) */

		dLambda = calloc(sModel->iNbVols, sizeof(double));
		dSigma = calloc(sModel->iNbVols, sizeof(double));

		if (!dLambda || !dSigma)
		{
			err = "Memory allocation faillure in GenericMidat_OneFactor_PDE";
			goto FREE_RETURN;
		}

		/* We need to convert the model */
		err = GenericMidat_ConvertModelLambda(	sModel->iNbVols,
												sModel->dTimes,
												sModel->dSigma,
												dLambda,
												&dSigma0,
												dSigma);

		if (err) goto FREE_RETURN;

		for (i=1; i<lNbTime; i++)
		{
			dt = dTimes[i] - dTimes[i-1];		
			dXtVar[i] = dt;
			dXtTotStd[i] /= dSigma[Get_Index(dTimes[i], sModel->dTimes, sModel->iNbVols)];
			dXtTotStd[i] = max(dXtTotStd[i], dXtTotStd[i-1]);
		}
	}

	for (i=0; i<sParams->lNbX; i++)
	{
		dXtDrift[i] = 0.0;
		dShortR[i] = 0.0;
	}

	/* Create the Grid */
	dx = sParams->dNbStd * dXtTotStd[lNbTime-1] / ((sParams->lNbX - 1.0) * 0.5);
	iIndexMid = (int) (sParams->lNbX / 2.0);

	dXtGrid[0] = -sParams->dNbStd * dXtTotStd[lNbTime-1];

	for (i=1; i<sParams->lNbX; i++)
	{
		dXtGrid[i] = dXtGrid[i-1] + dx;
	}
	
	/* Last Payoff */
	if (iEvalEvt[lNbTime-1])
	{
		/* Fill Recons Param */
		switch (sParams->iDiffMethod)
		{
			case 0:
			{
				sReconsParams->dCoefX1 = 1.0;
				break;
			}

			case 1:				
			{
				sReconsParams->dCoefX1 = 1.0 / dBeta[lNbTime-1];
				break;
			}

			case 2:
			{
				sReconsParams->dCoefX1 = dSigma[Get_Index(dTimes[lNbTime-1], sModel->dTimes, sModel->iNbVols)];
				break;
			}

			default:
			{
				err = "PDE Diff Method not recognised";
				goto FREE_RETURN;
			}
		}

		err = payoff_func(	dTimes[lNbTime-1],
							sFuncParmTab[lNbTime-1],
							0,
							sParams->lNbX-1,
							dXtGrid,
							sReconsParams,
							iNbProd,
							dNextValues);

		if (err) goto FREE_RETURN;

		iFirstEvent = 2;
	}
	else
	{
		iFirstEvent = 0;
	}
	
	/* Initialise the PDE */
	pde = calloc(1, sizeof(CNPDE_TEMP));

	if (!pde)
	{
		err = "Memory allocation faillure (2) in GenericMidat_OneFactor_PDE";
		goto FREE_RETURN;
	}

	num_f_pde_init(	pde,
					sParams->lNbX,
					iNbProd);
	
	for (step=lNbTime-2; step>=0; step--)
	{
		/* Fill the local volatility */
		for (i=0; i<sParams->lNbX; i++)
		{
			dXtVarGrid[i] = dXtVar[step+1];
		}

		if (sParams->iDiffMethod == 2)
		{
			/* Fill the drift as well */
			lambda = dLambda[Get_Index(dTimes[step+1], sModel->dTimes, sModel->iNbVols)];
			lambda *= dTimes[step+1] - dTimes[step];

			for (i=0; i<sParams->lNbX; i++)
			{
				dXtDrift[i] = -lambda * dXtGrid[i];
			}
		}

		num_f_pde_one_step_backward(	pde,
										sParams->lNbX,
										dXtGrid,
										iNbProd,
										dNextValues,
										dXtDrift,
										dXtVarGrid,
										dShortR,
										sParams->dTheta,
										dValues);

		/*	Eval payoff */
		if (iEvalEvt[step])
		{
			/* Fill Recons Param */
			switch (sParams->iDiffMethod)
			{
				case 0:
				{
					sReconsParams->dCoefX1 = 1.0;
					break;
				}

				case 1:				
				{
					if (iFirstEvent == 0)
					{
						sReconsParams->dCoefX1 = 1.0 / dBeta[step];
					}
					else
					{
						sReconsParams->dCoefX1 = 1.0 / dBeta[step+1];
					}

					break;
				}

				case 2:
				{
					sReconsParams->dCoefX1 = dSigma[Get_Index(dTimes[step], sModel->dTimes, sModel->iNbVols)];
					break;
				}

				default:
				{
					err = "PDE Diff Method not recognised";
					goto FREE_RETURN;
				}
			}
			
			err = payoff_func(	dTimes[step],
								sFuncParmTab[step],
								0,
								sParams->lNbX-1,
								dXtGrid,
								sReconsParams,
								iNbProd,
								dValues);

			if (err) goto FREE_RETURN;
			
			iFirstEvent ++;
		}

		if (sParams->iDiffMethod == 1 && step > 0 && (fabs(dBeta[step+1] - dBeta[step]) > 1.0E-04))
		{
			if (iFirstEvent == 1)
			{
				iFirstEvent++;
			}
			else
			{
				if (sParams->iAdaptGrid)
				{
					/* Create the new grid */
					dx = sParams->dNbStd * dXtTotStd[step] / ((sParams->lNbX - 1.0) * 0.5);				

					dNextXtGrid[0] = -sParams->dNbStd * dXtTotStd[step];

					for (i=1; i<sParams->lNbX; i++)
					{
						dNextXtGrid[i] = dNextXtGrid[i-1] + dx;
					}

					dXtGridTemp = dNextXtGrid;
				}
				else
				{
					dXtGridTemp = dXtGrid;
				}

				GenericMidat_OneFactor_ChangeBetaGrid(	0,
														sParams->lNbX-1,
														0,
														iNbProd-1,
														dXtGrid,
														dBeta[step+1],
														dXtGridTemp,
														dBeta[step],
														sModel,
														sParams,
														dValues,
														dNextValues);

				if (sParams->iAdaptGrid)
				{
					/* Swap Grids */
					dXtGridTemp = dNextXtGrid;		
					dNextXtGrid = dXtGrid;
					dXtGrid = dXtGridTemp;
				}

				/* Swap */
				dValuesTemp = dNextValues;		
				dNextValues = dValues;
				dValues = dValuesTemp;
			}			
		}

		/* Reverse Values */
		dValuesTemp = dNextValues;		
		dNextValues = dValues;
		dValues = dValuesTemp;
	}

	/* copy the result					*/
	for (i=0; i<iNbProd; i++)
	{
		dProdVal[i] = dNextValues[iIndexMid][i];
	}
	
FREE_RETURN:

	if (dXtGrid) free (dXtGrid);
	if (dNextXtGrid) free (dNextXtGrid);
	if (dXtDrift) free (dXtDrift);
	if (dXtVar) free (dXtVar);
	if (dXtVarGrid) free (dXtVarGrid);
	if (dShortR) free (dShortR);
	if (dXtTotStd) free (dXtTotStd);
	if (dBeta) free(dBeta);

	if (dNextValues) free_dmatrix(dNextValues, 0, sParams->lNbX - 1, 0, iNbProd - 1);
	if (dValues) free_dmatrix(dValues, 0, sParams->lNbX - 1, 0, iNbProd - 1);

	if (sReconsParams) free (sReconsParams);
	if (dLambda) free (dLambda);
	if (dSigma) free (dSigma);

	if (pde)
	{
		num_f_pde_free(pde, sParams->lNbX, iNbProd);
		free (pde);
	}

	return err;
}

/* New Value in (Xt, Zt) = Old Value in (Xt, dCoefX * Xt + dCoefZ * Zt) */
Err	GenericMidat_ValueJump_PDE(	long				lStartX,
								long				lEndX,
								long				lStartZ,
								long				lEndZ,
								long				lStartProd,
								long				lEndProd,
													
								double				*dGridX,
								double				*dGridZ,
								double				*dNextGridX,
								double				*dNextGridZ,
													
								int					iIndex,
								
								GENMIDAT_MODEL		sModel,
								GENMIDAT_PDEPAMS	sParams,
								GENMIDAT_RECONS		sReconsParams,

								double				***dTempValues,
								double				***dOldValues,
								double				***dNewValues)
{
Err		err = NULL;
double	dCoefX, dCoefZ;
double	dNewX, dNewZ, dNewZi;
int		newI, newJ;
double	dFirstX, dx, dx2, dFirstZ, dz, dz2;
double	dProba, dProbaDown, dProbaMid, dProbaUp;
int		i, j, k;

	if (fabs(sReconsParams->dCoefX1 - sReconsParams->dLastCoefX1) < 1.0E-10 &&
		fabs(sReconsParams->dCoefX2 - sReconsParams->dLastCoefX2) < 1.0E-10)
	{
		/* The reconstruction is only in Z */
		dCoefX = (sReconsParams->dLastCoefY1 - sReconsParams->dCoefY1) / sReconsParams->dCoefY2;
		dCoefZ = sReconsParams->dLastCoefY2 / sReconsParams->dCoefY2;

		dz = dNextGridZ[lStartZ + 1] - dNextGridZ[lStartZ];
		dz2 = dz * dz;
		dFirstZ = dNextGridZ[lStartZ];

		if (sParams->iJumpQuadra == 0)
		{
			/* Linear Interpolation */
			for (i=lStartX; i<=lEndX; i++)
			{
				dNewZi = dCoefX * dGridX[i];

				for (j=lStartZ; j<=lEndZ; j++)
				{
					dNewZ = dNewZi + dCoefZ * dGridZ[j];

					newJ = max(min(lStartZ + (int) ((dNewZ - dFirstZ) / dz), lEndZ - 1), lStartZ);
					dProba = (dNewZ - dNextGridZ[newJ]) / (dNextGridZ[newJ+1] - dNextGridZ[newJ]);

					for (k=lStartProd; k<=lEndProd; k++)
					{
						dNewValues[i][j][k] = (1.0 - dProba) * dOldValues[i][newJ][k] + dProba * dOldValues[i][newJ+1][k];
					}
				}
			}
		}
		else
		{
			/* Quadratic Interpolation */
			for (i=lStartX; i<=lEndX; i++)
			{
				dNewZi = dCoefX * dGridX[i];

				for (j=lStartZ; j<=lEndZ; j++)
				{
					dNewZ = dNewZi + dCoefZ * dGridZ[j];

					newJ = max(min(lStartZ + (int) ((dNewZ - dFirstZ) / dz + 0.5), lEndZ - 1), lStartZ + 1);

					dProbaDown = (dNewZ * (dNewZ - dNextGridZ[newJ] - dNextGridZ[newJ+1]) + dNextGridZ[newJ] * dNextGridZ[newJ+1]) / (2.0 * dz2);
					dProbaMid = (dNewZ * (dNewZ - dNextGridZ[newJ-1] - dNextGridZ[newJ+1]) + dNextGridZ[newJ-1] * dNextGridZ[newJ+1]) / (-dz2);
					dProbaUp = (dNewZ * (dNewZ - dNextGridZ[newJ-1] - dNextGridZ[newJ]) + dNextGridZ[newJ-1] * dNextGridZ[newJ]) / (2.0 * dz2);

					for (k=lStartProd; k<=lEndProd; k++)
					{
						dNewValues[i][j][k] = dProbaDown * dOldValues[i][newJ-1][k]
											+ dProbaMid * dOldValues[i][newJ][k]
											+ dProbaUp * dOldValues[i][newJ+1][k];
					}
				}
			}
		}
	}
	else
	{
		/* Reconstruction in X as well but no coef in Z */
		/* We just do the step in X and call the code again */

		dx = dNextGridX[lStartX + 1] - dNextGridX[lStartX];
		dx2 = dx * dx;
		dFirstX = dNextGridX[lStartX];

		dCoefX = sReconsParams->dLastCoefX1 / sReconsParams->dCoefX1;

		if (sParams->iJumpQuadra == 0)
		{
			for (i=lStartX; i<=lEndX; i++)
			{
				dNewX = dCoefX * dGridX[i];
				newI = max(min(lStartX + (int) ((dNewX - dFirstX) / dx), lEndX - 1), lStartX);				
				dProba = (dNewX - dNextGridX[newI]) / (dNextGridX[newI+1] - dNextGridX[newI]);

				for (j=lStartZ; j<=lEndZ; j++)
				{
					for (k=lStartProd; k<=lEndProd; k++)
					{
						dTempValues[i][j][k] = (1.0 - dProba) * dOldValues[newI][j][k] + dProba * dOldValues[newI+1][j][k];
					}
				}
			}
		}
		else
		{
			/* Quadratic Interpolation */
			for (i=lStartX; i<=lEndX; i++)
			{
				dNewX = dCoefX * dGridX[i];
				newI = max(min(lStartX + (int) ((dNewX - dFirstX) / dx + 0.5), lEndX - 1), lStartX + 1);

				dProbaDown = (dNewX * (dNewX - dNextGridX[newI] - dNextGridX[newI+1]) + dNextGridX[newI] * dNextGridX[newI+1]) / (2.0 * dx2);
				dProbaMid = (dNewX * (dNewX - dNextGridX[newI-1] - dNextGridX[newI+1]) + dNextGridX[newI-1] * dNextGridX[newI+1]) / (-dx2);
				dProbaUp = (dNewX * (dNewX - dNextGridX[newI-1] - dNextGridX[newI]) + dNextGridX[newI-1] * dNextGridX[newI]) / (2.0 * dx2);

				for (j=lStartZ; j<=lEndZ; j++)
				{
					for (k=lStartProd; k<=lEndProd; k++)
					{
						dTempValues[i][j][k] = dProbaDown * dOldValues[newI-1][j][k]
											+ dProbaMid * dOldValues[newI][j][k]
											+ dProbaUp * dOldValues[newI+1][j][k];
					}
				}
			}
		}

		/* Call the function to interpolate in Z */
		sReconsParams->dCoefX1 = 1.0;
		sReconsParams->dCoefX2 = 0.0;
		sReconsParams->dLastCoefX1 = 1.0;
		sReconsParams->dLastCoefX2 = 0.0;

		sReconsParams->dCoefY1 *= dCoefX;

		err = GenericMidat_ValueJump_PDE(	lStartX,
											lEndX,
											lStartZ,
											lEndZ,
											lStartProd,
											lEndProd,
											dGridX,
											dGridZ,
											dNextGridX,
											dNextGridZ,
											iIndex,
											sModel,
											sParams,
											sReconsParams,
											dTempValues,
											dTempValues,
											dNewValues);
	}

	return err;
}


Err	GenericMidat_TwoFactor_PDE(	/* Time Discretisation */
								long				lNbTime,
								double				*dTimes,

								/* Model */
								GENMIDAT_MODEL		sModel,

								/*	Product data */
								void				**sFuncParmTab, 
								int					*iEvalEvt,

								/*	Payoff Function */
								Err (*payoff_func)(/* Event */
													double	dTime,

													/* Parameters */
													void	*sFuncParm,																										
													
													/* Gride data	*/
													long	lStartIndex1,
													long	lEndIndex1,
													long	lStartIndex2,
													long	lEndIndex2,

													double	*dXt,
													double	*dZt,
													GENMIDAT_RECONS	sReconsParams,
																			
													/* Vector of results to be updated */
													int		iNbProd,
													double	***dProdVal),

								/* PDE Parameters */
								GENMIDAT_PDEPAMS	sParams,

								/* Result */
								int					iNbProd,
								double				*dProdVal)
{
Err	err = NULL;

double	*dXtGrid		= NULL,
		*dNextXtGrid	= NULL,
		*dXtGridTemp,
		*dZtGrid		= NULL,
		*dNextZtGrid	= NULL,
		*dZtGridTemp,
		*dXtVar			= NULL,
		*dZtVar			= NULL,
		*dCorrel		= NULL,
		**dXtDrift		= NULL,
		**dZtDrift		= NULL,
		**dXtVarGrid	= NULL,
		**dZtVarGrid	= NULL,
		**dShortR		= NULL,
		*dXtTotStd		= NULL,
		*dZtTotStd		= NULL,		
		*dLambda1		= NULL,
		*dSigma1		= NULL,
		*dLambda2		= NULL,
		*dSigma2		= NULL;

double	dSigma01, dSigma02;

double	*dBeta1			= NULL,
		*dBeta2			= NULL;

double	***dNextValues	= NULL,
		***dValues		= NULL,
		***dTempValues	= NULL,
		***dValuesTemp;

int		*ValueJump		= NULL;

GENMIDAT_RECONS	sReconsParams = NULL;

double	dt, dx, dz;

int		i, j, iIndexXMid, iIndexZMid, index;
int		step;
int		iFirstEvent;
int		iNbX, iNbZ, iNextNbX, iNextNbZ;
int		lx, ux, lz, uz;
double	lambda1, lambda2;

CNPDE_TEMP_2D_ADI	*pde = NULL;

	/* Memory allocations */
	dXtVar = calloc(lNbTime, sizeof(double));
	dZtVar = calloc(lNbTime, sizeof(double));
	dCorrel	= calloc(lNbTime, sizeof(double));
	dXtTotStd = calloc(lNbTime, sizeof(double));
	dZtTotStd = calloc(lNbTime, sizeof(double));
	ValueJump = calloc(lNbTime - 1, sizeof(int));

	dBeta1 = calloc(lNbTime+1, sizeof(double));
	dBeta2 = calloc(lNbTime+1, sizeof(double));

	sReconsParams = calloc(1, sizeof(genmidat_recons));

	if (!dXtVar || !dZtVar || !dCorrel || !dXtTotStd || !dZtTotStd || !ValueJump ||
		!dBeta1 || !dBeta2 || !sReconsParams)
	{
		err = "Memory allocation faillure in GenericMidat_TwoFactor_PDE";
		goto FREE_RETURN;
	}

	/* Fill the dynamics constant */
	sReconsParams->iDiffMethod = sParams->iDiffMethod;

	dXtTotStd[0] = 0.0;
	dZtTotStd[0] = 0.0;
	dCorrel[0] = sModel->dCorrel[0];

	for (i=1; i<lNbTime; i++)
	{
		dt = dTimes[i] - dTimes[i-1];
		index = Get_Index(dTimes[i], sModel->dTimes, sModel->iNbVols);

		dXtVar[i] = sModel->dSigma[index];
		dXtVar[i] *= dXtVar[i] * dt;
		dXtTotStd[i] = sqrt(dXtTotStd[i-1] * dXtTotStd[i-1] + dXtVar[i]);

		dZtVar[i] = sModel->dSigma2[index];
		dZtVar[i] *= dZtVar[i] * dt;				
		dZtTotStd[i] = sqrt(dZtTotStd[i-1] * dZtTotStd[i-1] + dZtVar[i]);

		dCorrel[i] = sModel->dCorrel[index];

		if (fabs(dCorrel[i] - dCorrel[i-1]) > 1.0E-10)
		{
			ValueJump[i-1] = index + 1;
		}
	}

	if (sParams->iDiffMethod == 1)
	{
		dBeta1[0] = sModel->dBeta[0];
		dBeta2[0] = sModel->dBeta2[0];

		dXtTotStd[0] = 0.0;
		dZtTotStd[0] = 0.0;

		j = 0;

		for (i=1; i<lNbTime; i++)
		{
			if (j < sModel->iNbVols-1 && dTimes[i] > sModel->dTimes[j])
			{
				j++;
			}

			if (j > 0)
			{
				dBeta1[i] = sModel->dBeta[j-1];
				dBeta2[i] = sModel->dBeta2[j-1];
			}
			else
			{
				dBeta1[i] = sModel->dBeta[0];
				dBeta2[i] = sModel->dBeta2[0];
			}

			dXtTotStd[i] *= dBeta1[i];
			dZtTotStd[i] *= dBeta2[i];

			if (!sParams->iAdaptGrid)
			{
				dXtTotStd[i] = max(dXtTotStd[i], dXtTotStd[i-1]);
				dZtTotStd[i] = max(dZtTotStd[i], dZtTotStd[i-1]);
			}

			dXtVar[i] *= dBeta1[i] * dBeta1[i];
			dZtVar[i] *= dBeta2[i] * dBeta2[i];
		}

		if (j < sModel->iNbVols-1 && dTimes[lNbTime-1] > sModel->dTimes[j])
		{
			j++;
		}

		dBeta1[lNbTime] = sModel->dBeta[j];
		dBeta2[lNbTime] = sModel->dBeta2[j];
	}

	if (sParams->iDiffMethod == 2)
	{
		/* We diffuse Yt = Xt / sigma(t) */

		dLambda1 = calloc(sModel->iNbVols, sizeof(double));
		dSigma1 = calloc(sModel->iNbVols, sizeof(double));
		dLambda2 = calloc(sModel->iNbVols, sizeof(double));
		dSigma2 = calloc(sModel->iNbVols, sizeof(double));

		if (!dLambda1 || !dSigma1 || !dLambda2 || !dSigma2)
		{
			err = "Memory allocation faillure in GenericMidat_TwoFactor_PDE";
			goto FREE_RETURN;
		}

		/* We need to convert the model */
		err = GenericMidat_ConvertModelLambda(	sModel->iNbVols,
												sModel->dTimes,
												sModel->dSigma,
												dLambda1,
												&dSigma01,
												dSigma1);

		if (err) goto FREE_RETURN;

		err = GenericMidat_ConvertModelLambda(	sModel->iNbVols,
												sModel->dTimes,
												sModel->dSigma2,
												dLambda2,
												&dSigma02,
												dSigma2);

		if (err) goto FREE_RETURN;

		for (i=1; i<lNbTime; i++)
		{
			dt = dTimes[i] - dTimes[i-1];		
			dXtVar[i] = dt;
			dZtVar[i] = dt;

			index = Get_Index(dTimes[i], sModel->dTimes, sModel->iNbVols);

			if (index > 0)
			{
				dXtTotStd[i] /= 0.5 * (sModel->dSigma[index] + sModel->dSigma[index-1]);
				dZtTotStd[i] /= 0.5 * (sModel->dSigma2[index] + sModel->dSigma2[index-1]);
			}
			else
			{
				dXtTotStd[i] /= 0.5 * (sModel->dSigma[index] + dSigma01);
				dZtTotStd[i] /= 0.5 * (sModel->dSigma2[index] + dSigma02);
			}

			dXtTotStd[i] = max(dXtTotStd[i], dXtTotStd[i-1]);
			dZtTotStd[i] = max(dZtTotStd[i], dZtTotStd[i-1]);
		}
	}

	/* Recalculates the number of points */
	iNbX = (int) (sParams->lNbX * sqrt(dXtTotStd[lNbTime-1] / dZtTotStd[lNbTime-1]) + 0.5);
	iNbZ = (int) (sParams->lNbX * sqrt(dZtTotStd[lNbTime-1] / dXtTotStd[lNbTime-1]) + 0.5);

	iNbX = max(((int) (iNbX / 2.0)) * 2 + 1, 3);
	iNbZ = max(((int) (iNbZ / 2.0)) * 2 + 1, 3);

	/* Allocates memory */
	dXtGrid = calloc(iNbX, sizeof(double));
	dNextXtGrid = calloc(iNbX, sizeof(double));
	dZtGrid = calloc(iNbZ, sizeof(double));
	dNextZtGrid = calloc(iNbZ, sizeof(double));	

	dXtDrift = dmatrix(0, iNbX-1, 0, iNbZ-1);
	dZtDrift = dmatrix(0, iNbX-1, 0, iNbZ-1);
	dXtVarGrid = dmatrix(0, iNbX-1, 0, iNbZ-1);
	dZtVarGrid = dmatrix(0, iNbX-1, 0, iNbZ-1);
	dShortR = dmatrix(0, iNbX-1, 0, iNbZ-1);
		
	dNextValues = f3tensor(0, iNbX-1, 0, iNbZ-1, 0, iNbProd - 1);
	dValues = f3tensor(0, iNbX-1, 0, iNbZ-1, 0, iNbProd - 1);	
	dTempValues = f3tensor(0, iNbX-1, 0, iNbZ-1, 0, iNbProd - 1);

	if (!dXtGrid || !dNextXtGrid || !dZtGrid || !dNextZtGrid ||
		!dXtDrift || !dZtDrift || !dXtVarGrid || !dZtVarGrid 
		|| !dShortR || !dNextValues || !dValues || !dTempValues)
	{
		err = "Memory allocation faillure (2) in GenericMidat_TwoFactor_PDE";
		goto FREE_RETURN;
	}

	for (i=0; i<iNbX; i++)
	{
		for (j=0; j<iNbZ; j++)
		{
			dXtDrift[i][j] = 0.0;
			dZtDrift[i][j] = 0.0;
			dShortR[i][j] = 0.0;
		}
	}

	/* Create the Two Grid */
	dx = sParams->dNbStd * dXtTotStd[lNbTime-1] / ((iNbX - 1.0) * 0.5);
	iIndexXMid = (int) (iNbX / 2.0);

	dXtGrid[0] = -sParams->dNbStd * dXtTotStd[lNbTime-1];

	for (i=1; i<iNbX; i++)
	{
		dXtGrid[i] = dXtGrid[i-1] + dx;
	}

	dz = sParams->dNbStd * dZtTotStd[lNbTime-1] / ((iNbZ - 1.0) * 0.5);
	iIndexZMid = (int) (iNbZ / 2.0);

	dZtGrid[0] = -sParams->dNbStd * dZtTotStd[lNbTime-1];

	for (i=1; i<iNbZ; i++)
	{
		dZtGrid[i] = dZtGrid[i-1] + dz;
	}

	lx = 0;
	ux = iNbX - 1;
	lz = 0;
	uz = iNbZ - 1;

	/* Last Payoff */
	if (iEvalEvt[lNbTime-1])
	{
		/* Fill Recons Param */
		switch (sParams->iDiffMethod)
		{
			case 0:
			{
				index = Get_Index(dTimes[lNbTime-1], sModel->dTimes, sModel->iNbVols);

				sReconsParams->dCoefX1 = 1.0;
				sReconsParams->dCoefX2 = 0.0;
				sReconsParams->dCoefY1 = sModel->dCorrel[index]	* sModel->dSigma2[index] / sModel->dSigma[index];
				sReconsParams->dCoefY2 = sqrt(1.0 - sModel->dCorrel[index] * sModel->dCorrel[index]);

				break;
			}

			case 1:				
			{
				/*
				sReconsParams->dCoefX1 = 1.0 / dBeta[lNbTime-1];
				*/
				break;
			}

			case 2:
			{
				index = Get_Index(dTimes[lNbTime-1], sModel->dTimes, sModel->iNbVols);

				sReconsParams->dCoefX1 = dSigma1[index];
				sReconsParams->dCoefX2 = 0.0;
				sReconsParams->dCoefY1 = sModel->dCorrel[index]	* sModel->dSigma2[index] / sModel->dSigma[index] * dSigma1[index];
				sReconsParams->dCoefY2 = sqrt(1.0 - sModel->dCorrel[index] * sModel->dCorrel[index]) * dSigma2[index];

				break;
			}

			default:
			{
				err = "PDE Diff Method not recognised";
				goto FREE_RETURN;
			}
		}

		err = payoff_func(	dTimes[lNbTime-1],
							sFuncParmTab[lNbTime-1],
							lx,
							ux,
							lz,
							uz,
							dXtGrid,
							dZtGrid,
							sReconsParams,
							iNbProd,
							dNextValues);

		if (err) goto FREE_RETURN;

		iFirstEvent = 2;
	}
	else
	{
		iFirstEvent = 0;
	}

	/* Initialise the PDE */
	pde = calloc(1, sizeof(CNPDE_TEMP_2D_ADI));

	if (!pde)
	{
		err = "Memory allocation faillure (3) in GenericMidat_TwoFactor_PDE";
		goto FREE_RETURN;
	}

	num_f_pde_init_2d_adi(	pde,
							iNbX,
							iNbZ,
							iNbProd);
	
	for (step=lNbTime-2; step>=0; step--)
	{
		/* Fill the local volatility */
		for (i=0; i<iNbX; i++)
		{
			for (j=0; j<iNbZ; j++)
			{
				dXtVarGrid[i][j] = dXtVar[step+1];
				dZtVarGrid[i][j] = dZtVar[step+1];
			}
		}

		if (sParams->iDiffMethod == 2)
		{
			/* Fill the drift as well */
			index = Get_Index(dTimes[step+1], sModel->dTimes, sModel->iNbVols);
			lambda1 = dLambda1[index];
			lambda1 *= dTimes[step+1] - dTimes[step];
			lambda2 = dLambda2[index];
			lambda2 *= dTimes[step+1] - dTimes[step];

			for (i=0; i<iNbX; i++)
			{
				for (j=0; j<iNbZ; j++)
				{
					dXtDrift[i][j] = -lambda1 * dXtGrid[i];
					dZtDrift[i][j] = -lambda2 * dZtGrid[j];
				}
			}
		}

		num_f_pde_one_step_backward_2f_adi(	pde,
											iNbX,
											dXtGrid,
											iNbZ,
											dZtGrid,
											0,
											iNbProd - 1,
											dNextValues,
											dXtDrift,
											dZtDrift,
											dXtVarGrid,
											dZtVarGrid,
											dShortR,
											dValues,
											lx,
											ux,
											lz,
											uz);

		/*	Jump Valuation */
		if (ValueJump[step] > 0 && (sParams->iDiffMethod == 0 || sParams->iDiffMethod == 2))
		{
			/* Fill Recons Param */
			switch (sParams->iDiffMethod)
			{
				case 0:
				{
					index = ValueJump[step] - 1;

					sReconsParams->dCoefX1 = 1.0;
					sReconsParams->dCoefX2 = 0.0;
					sReconsParams->dCoefY1 = sModel->dCorrel[index]	* sModel->dSigma2[index] / sModel->dSigma[index];
					sReconsParams->dCoefY2 = sqrt(1.0 - sModel->dCorrel[index] * sModel->dCorrel[index]);					

					sReconsParams->dLastCoefX1 = 1.0;
					sReconsParams->dLastCoefX2 = 0.0;
					sReconsParams->dLastCoefY1 = sModel->dCorrel[index-1]	* sModel->dSigma2[index-1] / sModel->dSigma[index-1];
					sReconsParams->dLastCoefY2 = sqrt(1.0 - sModel->dCorrel[index-1] * sModel->dCorrel[index-1]);

					break;
				}

				case 1:				
				{
					/*
					sReconsParams->dCoefX1 = 1.0 / dBeta[lNbTime-1];
					*/
					break;
				}

				case 2:
				{
					index = ValueJump[step] - 1;

					sReconsParams->dCoefX1 = dSigma1[index];
					sReconsParams->dCoefX2 = 0.0;
					sReconsParams->dCoefY1 = sModel->dCorrel[index]	* sModel->dSigma2[index] / sModel->dSigma[index] * dSigma1[index];
					sReconsParams->dCoefY2 = sqrt(1.0 - sModel->dCorrel[index] * sModel->dCorrel[index]) * dSigma2[index];

					sReconsParams->dLastCoefX1 = dSigma1[index];
					sReconsParams->dLastCoefX2 = 0.0;
					sReconsParams->dLastCoefY1 = sModel->dCorrel[index-1] * sModel->dSigma2[index-1] / sModel->dSigma[index-1] * dSigma1[index];
					sReconsParams->dLastCoefY2 = sqrt(1.0 - sModel->dCorrel[index-1] * sModel->dCorrel[index-1]) * dSigma2[index];

					break;
				}

				default:
				{
					err = "PDE Diff Method not recognised";
					goto FREE_RETURN;
				}
			}

			if (sParams->iAdaptGrid)
			{
				/* Create the new grid */
				dx = sParams->dNbStd * dXtTotStd[step] / ((iNbX - 1.0) * 0.5);
				dz = sParams->dNbStd * dZtTotStd[step] / ((iNbZ - 1.0) * 0.5);

				dNextXtGrid[0] = -sParams->dNbStd * dXtTotStd[step];

				for (i=1; i<iNbX; i++)
				{
					dNextXtGrid[i] = dNextXtGrid[i-1] + dx;
				}

				dXtGridTemp = dNextXtGrid;

				dNextZtGrid[0] = -sParams->dNbStd * dZtTotStd[step];

				for (i=1; i<iNbZ; i++)
				{
					dNextZtGrid[i] = dNextZtGrid[i-1] + dz;
				}

				dZtGridTemp = dNextZtGrid;
			}
			else
			{
				dXtGridTemp = dXtGrid;
				dZtGridTemp = dZtGrid;

				iNextNbX = iNbX;
				iNextNbZ = iNbZ;
			}

			err = GenericMidat_ValueJump_PDE(	lx,
												ux,
												lz,
												uz,
												0,
												iNbProd-1,
												dXtGridTemp,
												dZtGridTemp,
												dXtGrid,
												dZtGrid,
												ValueJump[step] - 1,
												sModel,
												sParams,
												sReconsParams,
												dTempValues,
												dValues,
												dNextValues); 

			/* Swap */
			dValuesTemp = dNextValues;
			dNextValues = dValues;
			dValues = dValuesTemp;

			if (sParams->iAdaptGrid)
			{
				/* Swap Grids */
				dXtGridTemp = dNextXtGrid;		
				dNextXtGrid = dXtGrid;
				dXtGrid = dXtGridTemp;

				dZtGridTemp = dNextZtGrid;		
				dNextZtGrid = dZtGrid;
				dZtGrid = dZtGridTemp;
			}
		}

		/*	Eval payoff */
		if (iEvalEvt[step])
		{
			/* Fill Recons Param */
			switch (sParams->iDiffMethod)
			{
				case 0:
				{
					index = Get_Index(dTimes[step], sModel->dTimes, sModel->iNbVols);

					sReconsParams->dCoefX1 = 1.0;
					sReconsParams->dCoefX2 = 0.0;
					sReconsParams->dCoefY1 = sModel->dCorrel[index]	* sModel->dSigma2[index] / sModel->dSigma[index];
					sReconsParams->dCoefY2 = sqrt(1.0 - sModel->dCorrel[index] * sModel->dCorrel[index]);

					break;
				}

				case 1:				
				{
					/*
					sReconsParams->dCoefX1 = 1.0 / dBeta[lNbTime-1];
					*/
					break;
				}

				case 2:
				{
					index = Get_Index(dTimes[step], sModel->dTimes, sModel->iNbVols);

					sReconsParams->dCoefX1 = dSigma1[index];
					sReconsParams->dCoefX2 = 0.0;
					sReconsParams->dCoefY1 = sModel->dCorrel[index]	* sModel->dSigma2[index] / sModel->dSigma[index] * dSigma1[index];
					sReconsParams->dCoefY2 = sqrt(1.0 - sModel->dCorrel[index] * sModel->dCorrel[index]) * dSigma2[index];

					break;
				}

				default:
				{
					err = "PDE Diff Method not recognised";
					goto FREE_RETURN;
				}
			}

			err = payoff_func(	dTimes[step],
								sFuncParmTab[step],
								lx,
								ux,
								lz,
								uz,
								dXtGrid,
								dZtGrid,
								sReconsParams,
								iNbProd,
								dValues);

			if (err) goto FREE_RETURN;

			iFirstEvent ++;
		}

		if (sParams->iDiffMethod == 1 && step > 0 && (fabs(dBeta1[step+1] - dBeta1[step]) > 1.0E-04))
		{
			if (iFirstEvent == 1)
			{
				iFirstEvent++;
			}
			else
			{
				if (sParams->iAdaptGrid)
				{
					/* Create the new grid */
					dx = sParams->dNbStd * dXtTotStd[step] / ((iNbX - 1.0) * 0.5);
					dz = sParams->dNbStd * dZtTotStd[step] / ((iNbZ - 1.0) * 0.5);

					iNextNbX = iNbX;
					iNextNbZ = iNbZ;

					dNextXtGrid[0] = -sParams->dNbStd * dXtTotStd[step];

					for (i=1; i<iNbX; i++)
					{
						dNextXtGrid[i] = dNextXtGrid[i-1] + dx;
					}

					dXtGridTemp = dNextXtGrid;

					dNextZtGrid[0] = -sParams->dNbStd * dZtTotStd[step];

					for (i=1; i<iNbZ; i++)
					{
						dNextZtGrid[i] = dNextZtGrid[i-1] + dz;
					}

					dZtGridTemp = dNextZtGrid;
				}
				else
				{
					dXtGridTemp = dXtGrid;
					dZtGridTemp = dZtGrid;

					iNextNbX = iNbX;
					iNextNbZ = iNbZ;
				}
			
				GenericMidat_TwoFactor_ChangeBetaGrid(	0,
														iNbX-1,
														dXtGrid,

														0,
														iNbZ-1,
														dZtGrid,

														0,
														iNextNbX-1,
														dXtGridTemp,

														0,
														iNextNbZ-1,
														dZtGridTemp,

														ValueJump[step] - 2,

														sModel,
														sParams,
														0,
														iNbProd-1,
														dTempValues,
														dValues,
														dNextValues);

				if (sParams->iAdaptGrid)
				{
					/* Swap Grids */
					dXtGridTemp = dNextXtGrid;		
					dNextXtGrid = dXtGrid;
					dXtGrid = dXtGridTemp;

					dZtGridTemp = dNextZtGrid;		
					dNextZtGrid = dZtGrid;
					dZtGrid = dZtGridTemp;
				}

				/* Swap */
				dValuesTemp = dNextValues;		
				dNextValues = dValues;
				dValues = dValuesTemp;
			}			
		}

		/* Reverse Values */
		dValuesTemp = dNextValues;		
		dNextValues = dValues;
		dValues = dValuesTemp;
	}

	/* copy the result					*/
	for (i=0; i<iNbProd; i++)
	{
		dProdVal[i] = dNextValues[iIndexXMid][iIndexZMid][i];
	}
	
FREE_RETURN:

	if (dXtVar) free (dXtVar);
	if (dZtVar) free (dZtVar);
	if (dCorrel) free(dCorrel);
	if (dXtTotStd) free (dXtTotStd);
	if (dZtTotStd) free (dZtTotStd);
	if (ValueJump) free (ValueJump);

	if (dBeta1) free(dBeta1);
	if (dBeta2) free(dBeta2);
	
	if (dXtGrid) free (dXtGrid);
	if (dNextXtGrid) free (dNextXtGrid);
	if (dZtGrid) free (dZtGrid);
	if (dNextZtGrid) free (dNextZtGrid);

	if (sReconsParams) free (sReconsParams);
	if (dLambda1) free (dLambda1);
	if (dSigma1) free (dSigma1);
	if (dLambda2) free (dLambda2);
	if (dSigma2) free (dSigma2);

	if (dXtDrift) free_dmatrix(dXtDrift, 0, iNbX-1, 0, iNbZ-1);
	if (dZtDrift) free_dmatrix(dZtDrift, 0, iNbX-1, 0, iNbZ-1);
	if (dXtVarGrid) free_dmatrix(dXtVarGrid, 0, iNbX-1, 0, iNbZ-1);
	if (dZtVarGrid) free_dmatrix(dZtVarGrid, 0, iNbX-1, 0, iNbZ-1);
	if (dShortR) free_dmatrix(dShortR, 0, iNbX-1, 0, iNbZ-1);
	
	if (dNextValues) free_f3tensor(dNextValues, 0, iNbX-1, 0, iNbZ-1, 0, iNbProd - 1);
	if (dValues) free_f3tensor(dValues, 0, iNbX-1, 0, iNbZ-1, 0, iNbProd - 1);
	if (dTempValues) free_f3tensor(dTempValues, 0, iNbX-1, 0, iNbZ-1, 0, iNbProd - 1);	

	if (pde)
	{
		num_f_pde_free_2d_adi(pde, iNbX, iNbZ, iNbProd);
		free (pde);
	}

	return err;
}


Err	GenericMidat_OneFactor_MC(	/* Time Discretisation */
								int					iNbEvent,
								double				*dTimes,

								/* Model */
								GENMIDAT_MODEL		sModel,

								/*	Product data */
								void				**sFuncParmTab, 

								/*	Payoff Function */
								Err (*payoff_func)(/* Event */
													double	dTime,

													/* Parameters */
													void			*sFuncParm,																								
													
													/* Gride data	*/
													long			lNumPaths,
													double			*dXt,
													GENMIDAT_RECONS	sReconsParams,
																			
													/* Vector of results to be updated */
													int				iNbProd,
													double			**dProdVal),

								/* MC Parameters */
								GENMIDAT_PDEPAMS	sParams,

								/* for Optimisation of exercise boundary */
								int					iDoOptim,
								int					*iOptimise,
								MCEBPARAMS			sMCEBParams,

								/* Result */
								int					iNbProd,
								double				**dProdVal)
{
Err	err = NULL;

double	*dXt			= NULL,
		*dInitGauss		= NULL,
		*dBrownian		= NULL,
		**dPayoff		= NULL,
		**dResEvt		= NULL,
		***dSaveValues	= NULL;

GENMIDAT_RECONS	sReconsParams = NULL;

long	iNumPaths, rand;
double	step, prob;
long	seed = -123456789;
double	dLastTime, dStdX;
int		iLastIndex, iIndex;
double	dCoef;

long	i, j, k;

	iNumPaths = 2 * ((long) (sParams->lNbPaths / 2)) + 1;

	dXt = calloc(iNumPaths, sizeof(double));
	dInitGauss = calloc(iNumPaths, sizeof(double));
	dBrownian = calloc(iNumPaths, sizeof(double));
	dPayoff = dmatrix(0, iNumPaths-1, 0, iNbProd - 1);	
	dResEvt = dmatrix(0, iNumPaths-1, 0, iNbProd - 1);

	sReconsParams = calloc(1, sizeof(genmidat_recons));

	if (!dXt || !dInitGauss || !dBrownian || !dPayoff || !dResEvt || !sReconsParams)
	{
		err = "Memory allocation (2) failure in GenericMidat_OneFactor_MC";
		goto FREE_RETURN;
	}

	if (iDoOptim)
	{
		dSaveValues = f3tensor(0, iNbEvent - 1, 0, sMCEBParams->iNbIndex, 0, iNumPaths - 1);		

		if (!dSaveValues)
		{
			err = "Memory allocation (2) failure in GenericMidat_OneFactor_MC";
			goto FREE_RETURN;
		}
	}

	/* Gauss initialisation */
	iNumPaths -= 1;
	iNumPaths /= 2;	
	step = 0.5 / (iNumPaths + 1);
	prob = step;

	/* Generation of the fractiles of the gaussian */
	for (i=0; i<iNumPaths; i++)
	{
		dInitGauss[i] = inv_cumnorm_fast(prob);
		dInitGauss[iNumPaths + i + 1] = -dInitGauss[i];
		prob += step;
	}

	dInitGauss[iNumPaths] = 0.0;	
	iNumPaths *= 2;
	iNumPaths += 1;

	/* Initialisation */
	for (i=0; i<iNumPaths; i++)
	{
		dXt[i] = 0.0;
	}	

	dLastTime = 0.0;
	iLastIndex = -1;	

	for (j=0; j<iNbEvent; j++)
	{
		/* Calculates StdX */
		dStdX = 0.0;
		iIndex = Get_Index(dTimes[j], sModel->dTimes, sModel->iNbVols);

		for (k=iLastIndex+1; k<=iIndex; k++)
		{
			dStdX += sModel->dSigma[k] * sModel->dSigma[k] * (sModel->dTimes[k] - dLastTime);
			dLastTime = sModel->dTimes[k];
		}

		dStdX = sqrt(dStdX);

		/* We diffuse Beta * Xt */
		dStdX *= sModel->dBeta[iIndex];

		/* Generate the random numbers */		
		for (i=0; i<iNumPaths-1; i++)
		{
			/* rand = random_int(nbPaths-1-i, &seed) + i; */
			rand = i + (int) ((iNumPaths-i) * uniform (&seed));
			dBrownian[i] = dInitGauss[rand];
			dInitGauss[rand] = dInitGauss[i];
			dInitGauss[i] = dBrownian[i];
		}

		dBrownian[iNumPaths-1] = dInitGauss[iNumPaths-1];

		/* Update State variabes */

		if (j > 0)
		{
			dCoef = sModel->dBeta[iIndex] / sModel->dBeta[iLastIndex];
		}
		else
		{
			dCoef = 0.0;
		}

		for (i=0; i<iNumPaths; i++)
		{
			dXt[i] = dCoef * dXt[i] + dStdX * dBrownian[i];
		}

		/* Payoff Evaluation */	

		if (sFuncParmTab[j])
		{
			sReconsParams->dCoefX1 = 1.0 / sModel->dBeta[iIndex];

			err = payoff_func(	dTimes[j],
								sFuncParmTab[j],
								iNumPaths,
								dXt,
								sReconsParams,
								iNbProd,
								dResEvt);

			if (err) goto FREE_RETURN;
		}

		/* Update the sum */
		for (i=0; i<iNumPaths; i++)
		{
			for (k=0; k<iNbProd; k++)
			{
				dPayoff[i][k] += dResEvt[i][k];
			}
		}

		if (iDoOptim)
		{
			for (i=0; i<iNumPaths; i++)
			{
				for (k=0; k<sMCEBParams->iNbIndex; k++)
				{
					dSaveValues[j][k][i] = dResEvt[i][sMCEBParams->iColBound+k];
				}

				if (iOptimise[j])
				{
					dSaveValues[j][sMCEBParams->iNbIndex][i] = dResEvt[i][sMCEBParams->iColPay];
				}
				else
				{
					dSaveValues[j][sMCEBParams->iNbIndex][i] = 0.0;
				}
			}
		}

		iLastIndex = iIndex;
	}
	
	for (k=0; k<iNbProd; k++)
	{
		dProdVal[k][0] = 0.0;
		dProdVal[k][1] = 0.0;
	}

	for (i=0; i<iNumPaths; i++)
	{
		for (k=0; k<iNbProd; k++)
		{
			dProdVal[k][0] += dPayoff[i][k] / iNumPaths;
			dProdVal[k][1] += dPayoff[i][k] * dPayoff[i][k] / iNumPaths;
		}
	}

	for (k=0; k<iNbProd; k++)
	{	
		dProdVal[k][1] = (dProdVal[k][1] - dProdVal[k][0] * dProdVal[k][0]) / iNumPaths;

		if (dProdVal[k][1] > 0.0)
		{
			dProdVal[k][1] = sqrt(dProdVal[k][1]);
		}
		else
		{
			dProdVal[k][1] = 0.0;
		}
	}

	if (iDoOptim)
	{
		if (sMCEBParams->iFindBestOptim)
		{
			err = find_best_optim_dates(	dSaveValues,
											iNbEvent,
											sMCEBParams->iFindBestOptim,
											iNumPaths,
											iOptimise,
											sMCEBParams,
											&(dProdVal[iNbProd][0]),
											&(dProdVal[iNbProd][1]));
		}
		else
		{
			err = find_and_optimise_boundary(	dSaveValues,
												iNbEvent,
												iNumPaths,
												iOptimise,
												sMCEBParams,
												&(dProdVal[iNbProd][0]),
												&(dProdVal[iNbProd][1]));
		}

		if (err) goto FREE_RETURN;		
	}

FREE_RETURN:

	if (dXt) free (dXt);
	if (dInitGauss) free (dInitGauss);
	if (dBrownian) free (dBrownian);
	if (dPayoff) free_dmatrix(dPayoff, 0, iNumPaths-1, 0, iNbProd - 1);
	if (dResEvt) free_dmatrix(dResEvt, 0, iNumPaths-1, 0, iNbProd - 1);

	if (sReconsParams) free(sReconsParams);

	if (dSaveValues) free_f3tensor(dSaveValues, 0, iNbEvent - 1, 0, sMCEBParams->iNbIndex, 0, iNumPaths - 1);	

	return err;

}	

Err	GenericMidat_TwoFactor_MC(	/* Time Discretisation */
								int					iNbEvent,
								double				*dTimes,

								/* Model */
								GENMIDAT_MODEL		sModel,

								/*	Product data */
								void				**sFuncParmTab, 

								/*	Payoff Function */
								Err (*payoff_func)(/* Event */
													double	dTime,

													/* Parameters */
													void			*sFuncParm,																								
													
													/* Gride data	*/
													long			lNumPaths,
													double			*dXt,
													double			*dYt,
													GENMIDAT_RECONS	sReconsParams,
																			
													/* Vector of results to be updated */
													int				iNbProd,
													double			**dProdVal),

								/* MC Parameters */
								GENMIDAT_PDEPAMS	sParams,

								/* for Optimisation of exercise boundary */
								int					iDoOptim,
								int					*iOptimise,
								MCEBPARAMS			sMCEBParams,

								/* Result */
								int					iNbProd,
								double				**dProdVal)
{
Err	err = NULL;

double	*dXt			= NULL,
		*dYt			= NULL,
		*dInitGauss		= NULL,
		**dBrownian		= NULL,
		**dPayoff		= NULL,
		**dResEvt		= NULL,
		***dSaveValues	= NULL;

GENMIDAT_RECONS	sReconsParams = NULL;

long	iNumPaths, rand;
double	step, prob;
long	seed = -123456789;
double	dLastTime, dVarX, dStdX, dVarY, dStdY, dCovar, dCorrel, dCorrel2;
int		iLastIndex, iIndex;
double	dCoef11, dCoef12, dCoef21, dCoef22;
double	dLastCoef11, dLastCoef12, dLastCoef21, dLastCoef22;
double	dDelta, dTransCoef11, dTransCoef12, dTransCoef21, dTransCoef22;
double	temp;
long	i, j, k;

	iNumPaths = 2 * ((long) (sParams->lNbPaths / 2)) + 1;

	dXt = calloc(iNumPaths, sizeof(double));
	dYt = calloc(iNumPaths, sizeof(double));
	dInitGauss = calloc(iNumPaths, sizeof(double));
	dBrownian = dmatrix(0, 1, 0, iNumPaths-1);
	dPayoff = dmatrix(0, iNumPaths-1, 0, iNbProd - 1);	
	dResEvt = dmatrix(0, iNumPaths-1, 0, iNbProd - 1);

	sReconsParams = calloc(1, sizeof(genmidat_recons));

	if (!dXt || !dYt || !dInitGauss || !dBrownian || !dPayoff || !dResEvt || !sReconsParams)
	{
		err = "Memory allocation (2) failure in GenericMidat_OneFactor_MC";
		goto FREE_RETURN;
	}

	if (iDoOptim)
	{
		dSaveValues = f3tensor(0, iNbEvent - 1, 0, sMCEBParams->iNbIndex, 0, iNumPaths - 1);
		
		if (!dSaveValues)
		{
			err = "Memory allocation (2) failure in GenericMidat_OneFactor_MC";
			goto FREE_RETURN;
		}
	}

	/* Gauss initialisation */
	iNumPaths -= 1;
	iNumPaths /= 2;	
	step = 0.5 / (iNumPaths + 1);
	prob = step;

	/* Generation of the fractiles of the gaussian */
	for (i=0; i<iNumPaths; i++)
	{
		dInitGauss[i] = inv_cumnorm_fast(prob);
		dInitGauss[iNumPaths + i + 1] = -dInitGauss[i];
		prob += step;
	}

	dInitGauss[iNumPaths] = 0.0;	
	iNumPaths *= 2;
	iNumPaths += 1;

	/* Initialisation */
	for (i=0; i<iNumPaths; i++)
	{
		dXt[i] = 0.0;
		dYt[i] = 0.0;
	}

	sReconsParams->dCoefX1 = 1.0;

	dLastTime = 0.0;
	iLastIndex = -1;	

	for (j=0; j<iNbEvent; j++)
	{
		/* Calculates StdX */
		dVarX = 0.0;
		dVarY = 0.0;
		dCovar = 0.0;

		iIndex = Get_Index(dTimes[j], sModel->dTimes, sModel->iNbVols);

		for (k=iLastIndex+1; k<=iIndex; k++)
		{
			dVarX += sModel->dSigma[k] * sModel->dSigma[k] * (sModel->dTimes[k] - dLastTime);
			dVarY += sModel->dSigma2[k] * sModel->dSigma2[k] * (sModel->dTimes[k] - dLastTime);
			dCovar += sModel->dCorrel[j] * sModel->dSigma[k] * sModel->dSigma2[k] * (sModel->dTimes[k] - dLastTime);

			dLastTime = sModel->dTimes[k];
		}		

		/* We diffuse Beta1 * Xt + Beta2 * Yt */
		dCoef11 = sModel->dBeta[iIndex];
		dCoef12 = sModel->dBeta2[iIndex];

		dStdX = dCoef11 * dCoef11 * dVarX
				+ dCoef12 * dCoef12 * dVarY
				+ 2.0 * dCoef11 * dCoef12 * dCovar;

		dStdX = sqrt(dStdX);

		/* we diffuse (Beta1 - NextBeta1) * Xt + (Beta2 - NextBeta2) * Yt */
		if (j < sModel->iNbVols - 1)
		{
			dCoef21 = sModel->dBeta[iIndex] - sModel->dBeta[iIndex+1];
			dCoef22 = sModel->dBeta2[iIndex] - sModel->dBeta2[iIndex+1];

			dStdY = dCoef21 * dCoef21 * dVarX
				+ dCoef22 * dCoef22 * dVarY
				+ 2.0 * dCoef21 * dCoef22 * dCovar;

			dStdY = sqrt(dStdY);
		}
		else
		{
			dCoef21 = dCoef11;
			dCoef22 = dCoef12;

			dStdY = dStdX;
		}		

		dCorrel = dCoef11 * dCoef21 * dVarX + dCoef12 * dCoef22 * dVarY
				+ (dCoef11 * dCoef22 + dCoef12 * dCoef21) * dCovar;

		dCorrel /= dStdX * dStdY;
		dCorrel2 = sqrt(1.0 - dCorrel * dCorrel);

		/* Generate the random numbers */		
		for (k=0; k<2; k++)
		{
			for (i=0; i<iNumPaths-1; i++)
			{
				/* rand = random_int(nbPaths-1-i, &seed) + i; */
				rand = i + (int) ((iNumPaths-i) * uniform (&seed));
				dBrownian[k][i] = dInitGauss[rand];
				dInitGauss[rand] = dInitGauss[i];
				dInitGauss[i] = dBrownian[k][i];
			}

			dBrownian[k][iNumPaths-1] = dInitGauss[iNumPaths-1];
		}

		/* Update State variabes */

		if (j > 0)
		{
			dDelta = dLastCoef11 * dLastCoef22 - dLastCoef21 * dLastCoef12;

			if (fabs(dDelta) < 1.0E-10)
			{
				err = "Cannot perform MC";
				goto FREE_RETURN;
			}

			dTransCoef11 = 1.0 / dDelta * (dCoef11 * dLastCoef22 - dCoef12 * dLastCoef21);
			dTransCoef12 = 1.0 / dDelta * (-dCoef11 * dLastCoef12 + dCoef12 * dLastCoef11);
			dTransCoef21 = 1.0 / dDelta * (dCoef21 * dLastCoef22 - dCoef22 * dLastCoef21);
			dTransCoef22 = 1.0 / dDelta * (-dCoef21 * dLastCoef12 + dCoef22 * dLastCoef11);
		}
		else
		{
			dTransCoef11 = 0.0;
			dTransCoef12 = 0.0;
			dTransCoef21 = 0.0;
			dTransCoef22 = 0.0;
		}

		for (i=0; i<iNumPaths; i++)
		{
			temp = dTransCoef11 * dXt[i] + dTransCoef12 * dYt[i] + dStdX * dBrownian[0][i];
			dYt[i] = dTransCoef21 * dXt[i] + dTransCoef22 * dYt[i] + dStdY * (dCorrel * dBrownian[0][i] + dCorrel2 * dBrownian[1][i]);
			dXt[i] = temp;
		}

		/* Payoff Evaluation */		
		if (sFuncParmTab[j])
		{
			err = payoff_func(	dTimes[j],
								sFuncParmTab[j],
								iNumPaths,
								dXt,
								dYt,
								sReconsParams,
								iNbProd,
								dResEvt);

			if (err) goto FREE_RETURN;
		}

		/* Update the sum */
		for (i=0; i<iNumPaths; i++)
		{
			for (k=0; k<iNbProd; k++)
			{
				dPayoff[i][k] += dResEvt[i][k];
			}
		}

		if (iDoOptim)
		{
			for (i=0; i<iNumPaths; i++)
			{
				for (k=0; k<sMCEBParams->iNbIndex; k++)
				{
					dSaveValues[j][k][i] = dResEvt[i][sMCEBParams->iColBound+k];
				}

				if (iOptimise[j])
				{
					dSaveValues[j][sMCEBParams->iNbIndex][i] = dResEvt[i][sMCEBParams->iColPay];
				}
				else
				{
					dSaveValues[j][sMCEBParams->iNbIndex][i] = 0.0;
				}
			}
		}

		iLastIndex = iIndex;

		dLastCoef11 = dCoef11;
		dLastCoef12 = dCoef12;
		dLastCoef21 = dCoef21;
		dLastCoef22 = dCoef22;
	}
	
	for (k=0; k<iNbProd; k++)
	{
		dProdVal[k][0] = 0.0;
		dProdVal[k][1] = 0.0;
	}

	for (i=0; i<iNumPaths; i++)
	{
		for (k=0; k<iNbProd; k++)
		{
			dProdVal[k][0] += dPayoff[i][k] / iNumPaths;
			dProdVal[k][1] += dPayoff[i][k] * dPayoff[i][k] / iNumPaths;
		}
	}

	for (k=0; k<iNbProd; k++)
	{	
		dProdVal[k][1] = (dProdVal[k][1] - dProdVal[k][0] * dProdVal[k][0]) / iNumPaths;

		if (dProdVal[k][1] > 0.0)
		{
			dProdVal[k][1] = sqrt(dProdVal[k][1]);
		}
		else
		{
			dProdVal[k][1] = 0.0;
		}
	}

	if (iDoOptim)
	{
		if (sMCEBParams->iFindBestOptim)
		{
			err = find_best_optim_dates(	dSaveValues,
											iNbEvent,
											sMCEBParams->iFindBestOptim,
											iNumPaths,
											iOptimise,
											sMCEBParams,
											&(dProdVal[iNbProd][0]),
											&(dProdVal[iNbProd][1]));
		}
		else
		{
			err = find_and_optimise_boundary(	dSaveValues,
												iNbEvent,
												iNumPaths,
												iOptimise,
												sMCEBParams,
												&(dProdVal[iNbProd][0]),
												&(dProdVal[iNbProd][1]));
		}

		if (err) goto FREE_RETURN;		
	}

FREE_RETURN:

	if (dXt) free (dXt);
	if (dInitGauss) free (dInitGauss);
	if (dBrownian) free_dmatrix(dBrownian, 0, 1, 0, iNumPaths-1);
	if (dPayoff) free_dmatrix(dPayoff, 0, iNumPaths-1, 0, iNbProd - 1);
	if (dResEvt) free_dmatrix(dResEvt, 0, iNumPaths-1, 0, iNbProd - 1);

	if (sReconsParams) free(sReconsParams);

	if (dSaveValues) free_f3tensor(dSaveValues, 0, iNbEvent - 1, 0, sMCEBParams->iNbIndex, 0, iNumPaths - 1);
	
	return err;

}	
