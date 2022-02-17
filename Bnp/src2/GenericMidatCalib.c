
#include "GenericMidatCalib.h"
#include "math.h"
#include "opfnctns.h"
#include "DiagCalibGen.h"

void genmidat_set_calibparams_default(GENMIDAT_CALIBPARAMS params)
{	
	params->iCalibLambda = 1;
	params->iShortIsCap = 1;
	params->dMinTime = 1.0;
	params->iSkipLast = 0;

	
	params->iContinueOnFail = 1;
	params->iInterpFailedVol = 1;
	params->iInterpVolType = 0;

	params->dBetaStep = 0.1;
	params->iMaxBetaIterations = 75;
	params->iMaxBadIter = 3;
	params->iMaxReallyBadIter = 2;
	params->iMaxSameDir = 4;

	params->dMinChangeVol = 0.01;
	params->dMaxChangeVol = 100.0;

	params->dLastBeta2 = 0.0;
	params->dLastBeta2Tolerance = 0.001;	
}

Err	genmidat_alloc_calibinfos(	int					iNbExe,
								GENMIDAT_CALIBINFOS	sCalibInfos)
{
Err	err = NULL;

	sCalibInfos->iNbExe = iNbExe;
	sCalibInfos->dLongVols = calloc(sCalibInfos->iNbExe, sizeof(double));
	sCalibInfos->dShortVols = calloc(sCalibInfos->iNbExe, sizeof(double));

	if (!sCalibInfos->dLongVols || !sCalibInfos->dShortVols)
	{
		err = "Memory allocation faillure in genmidat_alloc_calibinfos";
		return err;
	}

	return err;
}

void genmidat_free_calibinfos(	GENMIDAT_CALIBINFOS	sCalibInfos)
{
	if (sCalibInfos)
	{
		if (sCalibInfos->dLongVols) free (sCalibInfos->dLongVols);
		if (sCalibInfos->dShortVols) free (sCalibInfos->dShortVols);
	}
}

/* Interpolation in Sigma / (T - t) */
static Err generic_midat_interpolate_vol(int	iNbExe,
										 double	*dMaturity,
										 double	*dVols,
										 int	*iCalibVols,
										 int	iInterpType)
{
int		i, index_down, index_up;
double	dLastMaturity;
Err		err = NULL;

	if (iNbExe == 1 && iCalibVols[0] == 0)
	{
		err = "Cannot interpolate volatility with only one point";
		return err;
	}
	
	dLastMaturity = dMaturity[iNbExe-1] + (dMaturity[iNbExe-1] - dMaturity[iNbExe-2]);

	for (i=0; i<iNbExe; i++)
	{
		if (!iCalibVols[i])
		{
			/* We have a problem !!! */
			/* Find the indexes */
			index_down = i-1;
			while (index_down >= 0 && !iCalibVols[index_down]) index_down--;
			index_up = i+1;
			while (index_up < iNbExe && !iCalibVols[index_up]) index_up++;

			if (index_down < 0)
			{
				/* no down limit ! */				
				index_down = index_up;
				index_up = index_down + 1;
				while (index_up < iNbExe && !iCalibVols[index_up]) index_up++;

				if (index_up >= iNbExe)
				{
					err = "Cannot interpolate volatility, not enough points";
					return err;
				}
			}

			if (index_up >= iNbExe)
			{
				/* no up limit */
				index_up = index_down;
				index_down = index_down - 1;
				while (index_down >= 0 && !iCalibVols[index_down]) index_down--;

				if (index_down < 0)
				{
					err = "Cannot interpolate volatility, not enough points";
					return err;
				}
			}

			/* If everything's OK, let' interpolate */
			switch (iInterpType)
			{
				case 0:
				{
					/* Interpolation on Vol */
					dVols[i] = dVols[index_down]
						+	(dVols[index_up] - dVols[index_down] ) / (dMaturity[index_up] - dMaturity[index_down])
							* (dMaturity[i] - dMaturity[index_down]);
					break;
				}

				case 1:
				{
					/* Interpolation on Vol / (TEnd - TStar) */							
					dVols[i] = dVols[index_down] / (dLastMaturity - dMaturity[index_down])
						+	(dVols[index_up] / (dLastMaturity - dMaturity[index_up]) - 
							dVols[index_down] / (dLastMaturity - dMaturity[index_down])) / (dMaturity[index_up] - dMaturity[index_down])
							* (dMaturity[i] - dMaturity[index_down]);

					dVols[i] *= dLastMaturity - dMaturity[i];
				}
			}

			/* Check that vol > 0 */
			if (dVols[i] < 0.0)
			{
				/* It means the extrapolation failed */
				if (i < index_down)
				{
					/* we remove the volatility on index_down */
					iCalibVols[index_down] = 0;
					/* and we start again */
					i--;
				}
				else
				{
					/* we set equal to Vol up */
					dVols[i] = dVols[index_up];
				}
			}
		}
	}

	for (i=0; i<iNbExe; i++)
	{
		if (!iCalibVols[i]) iCalibVols[i] = 1;
	}

	return err;
}

Err generic_midat_calib(	/* Input */
							int						iNbExe,
							double					*dLongOption,
							double					*dShortOption,

							/* Model */
							GENMIDAT_MODEL			sModel,

							/* Parameters */
							GENMIDAT_CALIBPARAMS	sParams,

						   /* Optional Output */
							GENMIDAT_CALIBINFOS		sCalibInfos)
{
Err	err = NULL;
double	*dLongVols		= NULL,
		*dShortVols		= NULL;

int		*iCalibLong		= NULL,
		*iCalibShort	= NULL;

double	dShiftFwd, dShortFwd;
int	i;
	
	/* Memory allocation */
	dLongVols = calloc(iNbExe, sizeof(double));
	dShortVols = calloc(iNbExe, sizeof(double));
	iCalibLong = calloc(iNbExe, sizeof(int));
	iCalibShort = calloc(iNbExe, sizeof(int));
	
	if (!dLongVols || !dShortVols || !iCalibLong || !iCalibShort)
	{
		err = "Memory allocation faillure in generic_midat_calib_one_factor";
		goto FREE_RETURN;
	}

	/* Conversion of Long price into Vols */	
	for (i=0; i<iNbExe; i++)
	{
		if (dLongOption[i] / sModel->dNumeraire >= sModel->dForward[i] && fabs(dLongOption[i]) > 1.0E-10)
		{
			dShiftFwd = max(2.0 * fabs(sModel->dForward[i]), 1.0);

			err = srt_f_optimpvol(	dLongOption[i],
									sModel->dForward[i] + dShiftFwd,
									dShiftFwd,
									sModel->dTimes[i],
									sModel->dNumeraire,
									SRT_CALL,
									SRT_NORMAL,
									&(dLongVols[i]));

			if (err)
			{
				iCalibLong[i] = 0;
				err = NULL;
			}
			else
			{
				iCalibLong[i] = 1;
			}
		}
		else
		{
			if (sParams->iInterpFailedVol && 
				(fabs(dLongOption[i]) < 1.0E-10 ||
				(dLongOption[i] / sModel->dNumeraire - sModel->dForward[i] > -dLongOption[i] / 1000.0)))
			{
				iCalibLong[i] = 0;
			}
			else
			{
				err = serror("Long Option %d: Option < Intrinseq or cannot solve for volatility !" , i+1);
				goto FREE_RETURN;
			}
		}
	}

	if (sParams->iInterpFailedVol)
	{
		err = generic_midat_interpolate_vol(iNbExe,
											sModel->dTimes,
											dLongVols,
											iCalibLong,
											sParams->iInterpVolType);

		if (err) goto FREE_RETURN;
	}

	/* Conversion of Short prices into Vols */
	if (sParams->iCalibLambda)
	{
		if (sParams->iShortIsCap)
		{
			for (i=0; i<iNbExe-1; i++)
			{
				dShortFwd = sModel->dForward[i] - sModel->dForward[i+1];

				if (dShortOption[i] / sModel->dNumeraire >= dShortFwd
					 && fabs(dShortOption[i]) > 1.0E-10)
				{
					
					dShiftFwd = max(2.0 * fabs(dShortFwd), 1.0);
					

					err = srt_f_optimpvol(	dShortOption[i],
											dShortFwd + dShiftFwd,
											dShiftFwd,
											sModel->dTimes[i],
											sModel->dNumeraire,
											SRT_CALL,
											SRT_NORMAL,
											&(dShortVols[i]));

					if (err || fabs((dShortVols[i] - (dShortFwd + dShiftFwd) * 1.0e-9)) < 1.0e-15) // To avoid the imput of misleading information
					{
						iCalibShort[i] = 0;
						err = NULL;
					}
					else
					{
						iCalibShort[i] = 1;
					}
				}
				else
				{
					if (sParams->iInterpFailedVol && 
						(fabs(dShortOption[i]) < 1.0E-10 ||
						(dShortOption[i] / sModel->dNumeraire - dShortFwd > -dShortOption[i] / 1000.0)))
					{
						iCalibShort[i] = 0;
					}
					else
					{
						err = serror("Short Option %d: Option < Intrinseq or cannot solve for vol !" , i+1);
						goto FREE_RETURN;
					}
				}
			}
		}
		else
		{
			/* Short are partial long options !!! */			
			for (i=0; i<iNbExe-1; i++)
			{
				if (dShortOption[i] / sModel->dNumeraire >= sModel->dForward[i+1]
					&& fabs(dShortOption[i]) > 1.0E-10)
				{
					dShiftFwd = max(2.0 * fabs(sModel->dForward[i+1]), 1.0);

					err = srt_f_optimpvol(	dShortOption[i],
											sModel->dForward[i+1] + dShiftFwd,
											dShiftFwd,
											sModel->dTimes[i],
											sModel->dNumeraire,
											SRT_CALL,
											SRT_NORMAL,
											&(dShortVols[i]));

					if (err)
					{
						iCalibShort[i] = 0;
						err = NULL;
					}
					else
					{
						iCalibShort[i] = 1;
					}
				}
				else
				{
					if (sParams->iInterpFailedVol && 
						(fabs(dShortOption[i]) < 1.0E-10 ||
						(dShortOption[i] / sModel->dNumeraire - sModel->dForward[i+1] > -dShortOption[i] / 1000.0)))
					{
						iCalibShort[i] = 0;
					}
					else
					{
						err = serror("Short Option %d: Option < Intrinseq or cannot solve for vol !" , i+1);
						goto FREE_RETURN;
					}
				}
			}
		}

		if (sParams->iInterpFailedVol)
		{
			err = generic_midat_interpolate_vol(iNbExe-1,
												sModel->dTimes,
												dShortVols,
												iCalibShort,
												sParams->iInterpVolType);

			if (err) goto FREE_RETURN;
		}
	}	

	if (sModel->iNbFactor == 1)
	{
		err = generic_midat_calib_one_factor(	iNbExe,
												iCalibLong,
												dLongVols,
												iCalibShort,
												dShortVols,
												sModel,
												sParams,
												sCalibInfos);
	}
	else if (sModel->iNbFactor == 2)
	{
		err = generic_midat_calib_two_factor_total(	iNbExe,
													iCalibLong,
													dLongVols,
													iCalibShort,
													dShortVols,
													sModel,
													sParams,
													sCalibInfos);
	}
	else
	{
		err = "No more than 2 Factor implemented yet !";
		goto FREE_RETURN;
	}		

	if (err) goto FREE_RETURN;

FREE_RETURN:

	if (dLongVols) free(dLongVols);
	if (dShortVols) free(dShortVols);

	if (iCalibLong) free(iCalibLong);
	if (iCalibShort) free(iCalibShort);

	return err;
}

Err generic_midat_find_const_lambda(double	dT0,
									double	dT1,
									double	dRatio,
									double	*dLambda)
{
double	lam1, lam2, res1, res2;
double	target, diff;
int		i;

	target = 1.0 - dRatio;

	lam1 = 1.0;
	res1 = exp(-lam1 * dT1) / (1.0 - exp(-lam1 * dT0));
	
	if (fabs(res1 - target) < 0.0001)
	{
		*dLambda = res1;
		return NULL;
	}

	if (res1 > target)
	{
		lam2 = lam1 * 1.2;
	}
	else
	{
		lam2 = lam1 * 0.8;
	}

	i = 1;
	res2 = 1.0E10;

	while (i < 20 && fabs(res2 - target) > 0.0001)
	{
		res2 = exp(-lam2 * dT1) / (1.0 - exp(-lam2 * dT0));

		diff = (res2 - res1) / (lam2 - lam1);

		lam1 = lam2;
		res1 = res2;

		lam2 += (target - res2) / diff;		

		i++;
	}

	*dLambda = lam2;	

	return NULL;
}
Err generic_midat_calib_one_factor(	/* Input */
									int						iNbExe,
									int						*iCalibLong,
									double					*dLongVols,
									int						*iCalibShort,
									double					*dShortVols,

									/* Model */									
									GENMIDAT_MODEL			sModel,

									/* Parameters */
									GENMIDAT_CALIBPARAMS	sParams,

								   /* Optional Infos */
									GENMIDAT_CALIBINFOS		sCalibInfos)
{
Err	err = NULL;
double	ratio;
int	i;
	
	/* No check for the moment ! */
	
	/* First, Lambda Calibration */
	if (sParams->iCalibLambda && iNbExe > 1)
	{				
		/* Special case for the first Lambda */
		if (sParams->iShortIsCap)
		{
			ratio = 1.0 - dShortVols[0] / dLongVols[0];
		}
		else
		{
			ratio = dShortVols[0] / dLongVols[0];
		}

		/*
		err = generic_midat_find_const_lambda(	sModel->dTimes[0],
												sModel->dTimes[1],
												ratio,
												&(sModel->dLambda[0]));		
		
		sModel->dBeta[0] = 1.0 - exp(-sModel->dLambda[0] * sModel->dTimes[0]);				
		*/

		sModel->dBeta[0] = 1.0;

		for (i=1; i<iNbExe; i++)
		{
			if (sParams->iShortIsCap)
			{
				ratio = 1.0 - dShortVols[i-1] / dLongVols[i-1];
			}
			else
			{
				ratio = dShortVols[i-1] / dLongVols[i-1];
			}

			if (ratio > 0.0)
			{
				sModel->dBeta[i] = sModel->dBeta[i-1] * ratio;
			}
			else
			{
				err = serror("Cannot calibrate lambda on option %d", i);
				goto FREE_RETURN;
			}
		}
	}
	else if (sParams->iCalibLambda)
	{
		sModel->dBeta[0] = 1.0;
	}

	genmidat_convert_beta_into_lambda(sModel);

	/* Then volatility calibration */

	/* First Volatility */
	sModel->dSigma[0] = dLongVols[0] / sModel->dBeta[0];

	/* Then the others */
	for (i=1; i<iNbExe; i++)
	{
		sModel->dSigma[i] = dLongVols[i] * dLongVols[i] * sModel->dTimes[i] 
			- pow(sModel->dBeta[i] / sModel->dBeta[i-1] * dLongVols[i-1], 2.0) * sModel->dTimes[i-1];

		if (sModel->dSigma[i] > 0.0)
		{
			sModel->dSigma[i] = sqrt(sModel->dSigma[i] / (sModel->dTimes[i] - sModel->dTimes[i-1])) / sModel->dBeta[i];
		}
		else
		{
			if (sParams->iContinueOnFail)
			{
				smessage("Cannot calibrate sigma on option %d", i + 1);
				sModel->dSigma[i] = 0.0001;
			}
			else
			{				
				err = serror("Cannot calibrate sigma on option %d", i + 1);
				goto FREE_RETURN;
			}
		}

		/* Floor and Cap */
		sModel->dSigma[i] = max(sModel->dSigma[i], sParams->dMinChangeVol * sModel->dSigma[i-1]);
		sModel->dSigma[i] = min(sModel->dSigma[i], sParams->dMaxChangeVol * sModel->dSigma[i-1]);
	}			
	
FREE_RETURN:

	return err;
}

static double generic_midat_get_last_beta2(GENMIDAT_MODEL	sModel)
{
double	dLastBeta2, dLastMat;

	if (sModel->iNbVols > 1)
	{
		dLastMat = 2.0 * sModel->dTimes[sModel->iNbVols-1] - sModel->dTimes[sModel->iNbVols-2];
		dLastBeta2 = sModel->dBeta2[sModel->iNbVols-1] - exp(-sModel->dGamma * dLastMat) * sModel->dBeta[sModel->iNbVols-1];
	}
	else
	{
		dLastBeta2 = 0.0;
	}

	return dLastBeta2;
}


Err generic_midat_calib_two_factor_total(	/* Input */
									int						iNbExe,
									int						*iCalibLong,
									double					*dLongVols,
									int						*iCalibShort,
									double					*dShortVols,

									/* Model */									
									GENMIDAT_MODEL			sModel,

									/* Parameters */
									GENMIDAT_CALIBPARAMS	sParams,

								   /* Optional Infos */
									GENMIDAT_CALIBINFOS		sCalibInfos)
{
Err	err = NULL;
int		j, iter, lastiter, l;
GENMIDAT_CALIBINFOS	sUsedCalibInfos = NULL;
double	**dResIter		= NULL,
		**dResIterBeta	= NULL;

double	dStartBeta2, dLastSigma2, dSigma2Max, dLastBeta2;
int		dLastOpt, iIndexMax, iOptMax;
int		iNbMax, iNbBadIter, iNbReallyBadIter;
int		iSavedForceCalib;
int		iNbSameDir;

	/* Allocations */
	if (!sCalibInfos)
	{
		/* Allocate one */
		sUsedCalibInfos = calloc(1, sizeof(genmidat_calibinfos));

		if (!sUsedCalibInfos)
		{
			err = "Memory allocation faillure in generic_midat_calib_two_factor_total";
			goto FREE_RETURN;
		}
	}
	else
	{
		sUsedCalibInfos = sCalibInfos;
	}

	dResIter = dmatrix(0, 3, 0, sParams->iMaxBetaIterations);
	dResIterBeta = dmatrix(0, sParams->iMaxBetaIterations, 0, 1);

	if (!dResIter || !dResIterBeta)
	{
		err = "Memory allocation faillure in generic_midat_calib_two_factor_total";
		goto FREE_RETURN;
	}

	iNbMax = 0;
	iNbSameDir = 0;

	/* Switch off the Fail flag */
	iSavedForceCalib = sParams->iContinueOnFail;

	if (sParams->iCalibLambda)
	{
		sParams->iContinueOnFail = 0;
	}

	/* First Steps */
	err = generic_midat_calib_two_factor_new(	iNbExe,
												iCalibLong,
												dLongVols,
												iCalibShort,
												dShortVols,
												sModel,
												sParams,
												sUsedCalibInfos);

	
	if (!sParams->iCalibLambda) goto FREE_RETURN;

	if (!err)
	{
		dLastBeta2  = generic_midat_get_last_beta2(sModel);

		if (fabs(dLastBeta2 - sParams->dLastBeta2) < sParams->dLastBeta2Tolerance)
		{
			goto FREE_RETURN;
		}
	}
	else
	{
		dLastBeta2 = 0.0;
	}
	
	iter = 0;
	lastiter = 0;
	iNbBadIter = 0;
	iNbReallyBadIter = 0;

	/* First Stage: reach the last option */
	dResIter[0][iter] = sModel->dStartBeta2;	
	dResIter[1][iter] = sUsedCalibInfos->iFailedOpt;	
	dResIter[2][iter] = sUsedCalibInfos->dFailedSigma2;
	dResIter[3][iter] = iter;

	iter++;

	if (!err)
	{
		dResIterBeta[0][0] = sModel->dStartBeta2;
		dResIterBeta[0][1] = dLastBeta2;
		lastiter++;
	}

	dLastOpt = sUsedCalibInfos->iFailedOpt;
	dLastSigma2 = sUsedCalibInfos->dFailedSigma2;

	if (dLastOpt > iNbExe - 1.5)
	{
		iNbMax++;
	}

	iIndexMax = 0;
	iOptMax = dLastOpt;
	dSigma2Max = dLastSigma2;

	/* Next guess */
	dStartBeta2 = sModel->dStartBeta2 - sParams->dBetaStep;

	iNbSameDir = -1;	

	while (iter < sParams->iMaxBetaIterations)
	{
		/* Second Step */
		iter++;

		sModel->dStartBeta2 = dStartBeta2;

		err = generic_midat_calib_two_factor_new(	iNbExe,
													iCalibLong,
													dLongVols,
													iCalibShort,
													dShortVols,
													sModel,
													sParams,
													sUsedCalibInfos);

		if (!err)
		{
			lastiter++;

			dLastBeta2  = generic_midat_get_last_beta2(sModel);

			if (fabs(dLastBeta2 - sParams->dLastBeta2) < sParams->dLastBeta2Tolerance)
			{
				goto FREE_RETURN;
			}

			/* Save Res */
			l = 0;
			while (l < lastiter-1 && dResIterBeta[l][0] < dStartBeta2)
			{
				l++;
			}

			if (l < lastiter-1)
			{
				for (j=lastiter-2; j>=l; j--)
				{
					dResIterBeta[j+1][0] = dResIterBeta[j][0];
					dResIterBeta[j+1][1] = dResIterBeta[j][1];
				}
			}			

			dResIterBeta[l][0] = dStartBeta2;
			dResIterBeta[l][1] = dLastBeta2;
		}

		dLastOpt = sUsedCalibInfos->iFailedOpt;
		dLastSigma2 = sUsedCalibInfos->dFailedSigma2;

		if (dLastOpt > iNbExe - 0.5)
		{
			iNbMax++;
		}

		/* Save Res */
		l = 0;
		while (l < iter-1 && dResIter[0][l] < dStartBeta2)
		{
			l++;
		}

		if (l < iter-1)
		{
			for (j=iter-2; j>=l; j--)
			{
				dResIter[0][j+1] = dResIter[0][j];
				dResIter[1][j+1] = dResIter[1][j];
				dResIter[2][j+1] = dResIter[2][j];
				dResIter[3][j+1] = dResIter[3][j];
			}
		}

		if (l <= iIndexMax)
		{
			iIndexMax++;
		}

		dResIter[0][l] = dStartBeta2;
		dResIter[1][l] = dLastOpt;		
		dResIter[2][l] = dLastSigma2;
		dResIter[3][l] = iter;
		
		if (dLastOpt > iOptMax || (dLastOpt == iOptMax && dLastSigma2 > dSigma2Max))
		{
			if (dLastOpt > iOptMax)
			{
				/* True improvement */
				iNbSameDir = 0;
				iNbBadIter = 0;
				iNbReallyBadIter = 0;
			}

			iOptMax = dLastOpt;
			dSigma2Max = dLastSigma2;
			iIndexMax = l;						
		}
		else
		{
			iNbBadIter++;
		}

		/* Find Next Guess */
		if (err || lastiter == 1)
		{
			if (fabs(iNbSameDir) > sParams->iMaxSameDir)
			{
				/* we change of direction */
				if (iNbSameDir > 0)
				{
					dStartBeta2 = dResIter[0][0] - (dResIter[0][1] - dResIter[0][0]);
					iNbSameDir--;
				}
				else
				{
					dStartBeta2 = dResIter[0][iter-1] + (dResIter[0][iter-1] - dResIter[0][iter-2]);
					iNbSameDir++;
				}
			}
			else
			{
				if (iIndexMax > 0 && iIndexMax < iter - 1)
				{
					/* We are in the middle */
					/* Find the best second guess */				
					if (dResIter[1][iIndexMax-1] > dResIter[1][iIndexMax+1])
					{
						if (iNbBadIter <= sParams->iMaxBadIter)
						{
							dStartBeta2 = 0.5 * (dResIter[0][iIndexMax-1] + dResIter[0][iIndexMax]);
							
							iNbSameDir = min(iNbSameDir - 1, -1);												
						}
						else
						{
							dStartBeta2 = 0.5 * (dResIter[0][iIndexMax] + dResIter[0][iIndexMax+1]);

							iNbSameDir = max(iNbSameDir + 1, 1);
							
							iNbBadIter -= 2;

							if (iNbReallyBadIter > sParams->iMaxReallyBadIter)
							{
								/* look for the last good on the left */
								l = iIndexMax - 1;
								while (l >= 0 && dResIter[1][iIndexMax-1] == dResIter[1][l])
								{
									l--;
								}

								if (l < 0)
								{
									dStartBeta2 = dResIter[0][0] - (dResIter[0][1] - dResIter[0][0]);
								}
								else
								{
									dStartBeta2 = 0.5 * (dResIter[0][l] + dResIter[0][l+1]);
								}

								iNbReallyBadIter -= 2;
							}
							else
							{
								iNbReallyBadIter++;
							}
						}
					}
					else if (dResIter[1][iIndexMax-1] < dResIter[1][iIndexMax+1])
					{
						if (iNbBadIter <= sParams->iMaxBadIter)
						{
							dStartBeta2 = 0.5 * (dResIter[0][iIndexMax] + dResIter[0][iIndexMax+1]);

							iNbSameDir = max(iNbSameDir + 1, 1);							
						}
						else
						{
							dStartBeta2 = 0.5 * (dResIter[0][iIndexMax-1] + dResIter[0][iIndexMax]);

							iNbSameDir = min(iNbSameDir - 1, -1);

							iNbBadIter -= 2;

							if (iNbReallyBadIter > sParams->iMaxReallyBadIter)
							{
								/* look for the last good on the left */
								l = iIndexMax + 1;
								while (l < iter && dResIter[1][iIndexMax+1] == dResIter[1][l])
								{
									l++;
								}

								if (l == iter)
								{
									dStartBeta2 = dResIter[0][iter-1] + (dResIter[0][iter-1] - dResIter[0][iter-2]);
								}
								else
								{
									dStartBeta2 = 0.5 * (dResIter[0][l] + dResIter[0][l-1]);
								}

								iNbReallyBadIter -= 2;
							}
							else
							{
								iNbReallyBadIter++;
							}
						}
					}
					else if (dResIter[2][iIndexMax-1] > dResIter[2][iIndexMax+1])
					{
						if (iNbBadIter <= sParams->iMaxBadIter)
						{
							dStartBeta2 = 0.5 * (dResIter[0][iIndexMax-1] + dResIter[0][iIndexMax]);

							iNbSameDir = min(iNbSameDir - 1, -1);
						}
						else
						{
							dStartBeta2 = 0.5 * (dResIter[0][iIndexMax] + dResIter[0][iIndexMax+1]);

							iNbSameDir = max(iNbSameDir + 1, 1);						
							
							iNbBadIter -= 2;

							if (iNbReallyBadIter > sParams->iMaxReallyBadIter)
							{
								/* look for the last good on the left */
								l = iIndexMax - 1;
								while (l >= 0 && dResIter[1][iIndexMax-1] == dResIter[1][l])
								{
									l--;
								}

								if (l < 0)
								{
									dStartBeta2 = dResIter[0][0] - (dResIter[0][1] - dResIter[0][0]);
								}
								else
								{
									dStartBeta2 = 0.5 * (dResIter[0][l] + dResIter[0][l+1]);
								}

								iNbReallyBadIter -= 2;
							}
							else
							{
								iNbReallyBadIter++;
							}
						}
					}
					else
					{
						if (iNbBadIter <= sParams->iMaxBadIter)
						{
							dStartBeta2 = 0.5 * (dResIter[0][iIndexMax] + dResIter[0][iIndexMax+1]);

							iNbSameDir = max(iNbSameDir + 1, 1);
						}
						else
						{
							dStartBeta2 = 0.5 * (dResIter[0][iIndexMax-1] + dResIter[0][iIndexMax]);

							iNbSameDir = min(iNbSameDir - 1, -1);

							iNbBadIter -= 2;
							
							if (iNbReallyBadIter > sParams->iMaxReallyBadIter)
							{
								/* look for the last good on the left */
								l = iIndexMax + 1;
								while (l < iter && dResIter[1][iIndexMax+1] == dResIter[1][l])
								{
									l++;
								}

								if (l == iter)
								{
									dStartBeta2 = dResIter[0][iter-1] + (dResIter[0][iter-1] - dResIter[0][iter-2]);
								}
								else
								{
									dStartBeta2 = 0.5 * (dResIter[0][l] + dResIter[0][l-1]);
								}

								iNbReallyBadIter -= 2;
							}
							else
							{
								iNbReallyBadIter++;
							}
						}
					}
				}
				else if (iIndexMax == 0)
				{
					if (iNbBadIter <= sParams->iMaxBadIter)
					{
						dStartBeta2 = dResIter[0][iIndexMax] - (dResIter[0][iIndexMax+1] - dResIter[0][iIndexMax]);

						iNbSameDir = min(iNbSameDir - 1, -1);
					}
					else
					{
						dStartBeta2 = 0.5 * (dResIter[0][iIndexMax] + dResIter[0][iIndexMax+1]);

						iNbSameDir = max(iNbSameDir + 1, 1);
						iNbBadIter -= 2;
						iNbReallyBadIter++;
					}
				}
				else
				{
					if (iNbBadIter <= sParams->iMaxBadIter)
					{
						dStartBeta2 = dResIter[0][iIndexMax] + (dResIter[0][iIndexMax] - dResIter[0][iIndexMax-1]);

						iNbSameDir = max(iNbSameDir + 1, 1);
					}
					else
					{
						dStartBeta2 = 0.5 * (dResIter[0][iIndexMax-1] + dResIter[0][iIndexMax]);

						iNbSameDir = min(iNbSameDir - 1, -1);
						iNbBadIter -= 2;
						iNbReallyBadIter++;
					}
				}		
			}
		}
		else
		{
			/* We need to solve for Last Beta */
			dStartBeta2 = solve_for_next_coef_gen(	dResIterBeta,
													lastiter,
													sParams->dLastBeta2,
													sParams->dLastBeta2Tolerance,
													2);

			/* Check Limits */
		}
	}

	if (err && iSavedForceCalib)
	{
		/* do a last calibration  */
		sParams->iContinueOnFail = 1;

		err = generic_midat_calib_two_factor_new(	iNbExe,
													iCalibLong,
													dLongVols,
													iCalibShort,
													dShortVols,
													sModel,
													sParams,
													sUsedCalibInfos);

		if (err) goto FREE_RETURN;
	}

FREE_RETURN:

	/* Restore Flag */
	sParams->iContinueOnFail = iSavedForceCalib;

	if (!sCalibInfos && sUsedCalibInfos) free (sUsedCalibInfos);

	if (dResIter) free_dmatrix(dResIter, 0, 3, 0, sParams->iMaxBetaIterations);

	if (dResIterBeta) free_dmatrix(dResIterBeta, 0, sParams->iMaxBetaIterations, 0, 1);

	return err;
}

Err generic_midat_calib_two_factor_new(	/* Input */
									int						iNbExe,
									int						*iCalibLong,
									double					*dLongVols,
									int						*iCalibShort,
									double					*dShortVols,

									/* Model */									
									GENMIDAT_MODEL			sModel,

									/* Parameters */
									GENMIDAT_CALIBPARAMS	sParams,

								   /* Optional Infos */
									GENMIDAT_CALIBINFOS		sCalibInfos)
{
Err	err = NULL;
double	dSigma2;
double	dPhi0, dPhi1, dPhi2;
double	dVarX, dVarY, dCoVar;
double	dFact0, dFact1, dFact2;
double	dSumLam, dLastSumLam;
double	dExpSumLam ,dLastExpSumLam;
double	dt;
int	i;
double	dStartBeta, dStartBeta2;
double	dDeltaCalib;
	
	/* Initialisation */
	dVarX = 0.0;
	dCoVar = 0.0;
	dVarY = 0.0;

	dPhi0 = 0.0;
	dPhi1 = 0.0;
	dPhi2 = 0.0;

	if (iNbExe > 1)
	{
		dDeltaCalib = sModel->dTimes[1] - sModel->dTimes[0];
	}
	else
	{
		dDeltaCalib = sModel->dTimes[0];
	}

	dStartBeta = sModel->dStartBeta;
	dStartBeta2 = sModel->dStartBeta2;
	
	/* First, Lambda Calibration */
	if (sParams->iCalibLambda)
	{				
		/* Special case for the first Lambda */		
		sModel->dBeta[0] = 1.0;
		sModel->dBeta2[0] = dStartBeta2 + (sModel->dBeta[0] - dStartBeta) * exp(-sModel->dGamma * dDeltaCalib);
	}		
	
	sModel->dLambda[0] = -log(dStartBeta - sModel->dBeta[0]) / sModel->dTimes[0];

	dLastSumLam = 0.0;
	dLastExpSumLam = 1.0;
	dSumLam = sModel->dLambda[0] * sModel->dTimes[0];
	dExpSumLam = exp(-dSumLam);

	/* Then the others */
	for (i=0; i<iNbExe; i++)
	{	
		if (i > 0)
		{
			dt = sModel->dTimes[i] - sModel->dTimes[i-1];

					/* Calculate the Factors */		
			dFact0 = (1.0 / dExpSumLam / dExpSumLam 
					- 1.0 / dLastExpSumLam / dLastExpSumLam) / (2.0 * sModel->dLambda[i]);
			dFact1 = (exp(sModel->dGamma * sModel->dTimes[i]) / dExpSumLam / dExpSumLam
					- exp(sModel->dGamma * sModel->dTimes[i-1]) / dLastExpSumLam / dLastExpSumLam) / (2.0 * sModel->dLambda[i] + sModel->dGamma);
			dFact2 = (exp(2.0 * sModel->dGamma * sModel->dTimes[i]) / dExpSumLam / dExpSumLam
					- exp(2.0 * sModel->dGamma * sModel->dTimes[i-1]) / dLastExpSumLam / dLastExpSumLam) / (2.0 * (sModel->dLambda[i] + sModel->dGamma));
		}
		else
		{
			dt = sModel->dTimes[0];

			if (fabs(sModel->dLambda[0]) < 1.0E-10)
			{
				dFact0 = dt;
			}
			else
			{
				dFact0 = (exp(2.0 * sModel->dLambda[0] * sModel->dTimes[0]) - 1.0) / (2.0 * sModel->dLambda[0]);
			}

			dFact1 = (exp((2.0 * sModel->dLambda[0] + sModel->dGamma) * sModel->dTimes[0]) - 1.0) / (2.0 * sModel->dLambda[0] + sModel->dGamma);
			dFact2 = (exp(2.0 * (sModel->dLambda[0] + sModel->dGamma) * sModel->dTimes[0]) - 1.0) / (2.0 * (sModel->dLambda[0] + sModel->dGamma));
		}

		/* Calculate the Volatility */
		dSigma2 = dLongVols[i] * dLongVols[i] * sModel->dTimes[i] 
			- sModel->dBeta[i] * sModel->dBeta[i] * dVarX
			- 2.0 * sModel->dBeta[i] * sModel->dBeta2[i] * dCoVar
			- sModel->dBeta2[i] * sModel->dBeta2[i] * dVarY;
		
		dSigma2 /= sModel->dBeta[i] * sModel->dBeta[i] * dFact0 
			+ 2.0 * sModel->dBeta[i] * sModel->dBeta2[i] * sModel->dRhoAlpha * dFact1
			+ sModel->dBeta2[i] * sModel->dBeta2[i] * sModel->dAlpha2 * dFact2;		

		if (dSigma2 < 0.0)
		{
			if (sParams->iContinueOnFail)
			{
				dSigma2 = 0.0001 * 0.0001 * dt / dFact0;
				smessage("Cannot calibrate sigma on option %d", i + 1);
				err = NULL;
			}
			else
			{
				err = serror("Cannot calibrate sigma on option %d", i + 1);
				goto FREE_RETURN;
			}
		}

		/* Floor and Cap */
		if (i > 0)
		{
			dSigma2 = max(dSigma2, (sParams->dMinChangeVol * sModel->dSigma[i-1]) * (sParams->dMinChangeVol * sModel->dSigma[i-1]) * dt / dFact0);
			dSigma2 = min(dSigma2, (sParams->dMaxChangeVol * sModel->dSigma[i-1]) * (sParams->dMaxChangeVol * sModel->dSigma[i-1]) * dt / dFact0);
		}

		/* Update variables */		
		sModel->dSigma[i] = dSigma2 * dFact0;
		sModel->dCorrel[i] = sModel->dRhoAlpha * dSigma2 * dFact1;
		sModel->dSigma2[i] = sModel->dAlpha2 * dSigma2 * dFact2;

		dVarX += sModel->dSigma[i];
		dCoVar += sModel->dCorrel[i];
		dVarY += sModel->dSigma2[i];
				
		if (i < iNbExe - 1 && sParams->iCalibLambda)
		{						
			/* Find the next Lambda */
			dPhi0 = dVarX;
			dPhi1 = exp(-sModel->dGamma * (sModel->dTimes[i+1] - sModel->dTimes[i])) * dPhi1
				+ exp(-sModel->dGamma * sModel->dTimes[i+1]) * sModel->dCorrel[i];

			dPhi2 = exp(-2.0 * sModel->dGamma * (sModel->dTimes[i+1] - sModel->dTimes[i])) * dPhi2
				+ exp(-2.0 * sModel->dGamma * sModel->dTimes[i+1]) * sModel->dSigma2[i];

			sModel->dBeta[i+1] = dShortVols[i] * dShortVols[i] * sModel->dTimes[i]
				/ (dPhi0 + 2.0 * dPhi1 + dPhi2);
			sModel->dBeta[i+1] = sqrt(sModel->dBeta[i+1]);

			sModel->dLambda[i+1] = (-log(sModel->dBeta[i+1]) - dSumLam) / (sModel->dTimes[i+1] - sModel->dTimes[i]);			
			sModel->dBeta[i+1] = sModel->dBeta[i] - sModel->dBeta[i+1];			
		}

		if (i < iNbExe - 1)
		{
			sModel->dLambda2[i+1] = sModel->dLambda[i+1] + sModel->dGamma;
			sModel->dBeta2[i+1] = sModel->dBeta2[i] + (sModel->dBeta[i+1] - sModel->dBeta[i]) * exp(-sModel->dGamma * sModel->dTimes[i+1]);
			dLastSumLam = dSumLam;
			dLastExpSumLam = dExpSumLam;
			dSumLam = dLastSumLam + sModel->dLambda[i+1] * (sModel->dTimes[i+1] - sModel->dTimes[i]);
			dExpSumLam = exp(-dSumLam);
		}

		sModel->dSigma[i] = sqrt(sModel->dSigma[i] / dt);
		sModel->dSigma2[i] = sqrt(sModel->dSigma2[i] / dt);
		sModel->dCorrel[i] /= sModel->dSigma[i] * sModel->dSigma2[i] * dt;
	}

	/*
	genmidat_convert_beta_into_lambda(sModel);
	*/
	
FREE_RETURN:

	if (sCalibInfos)
	{
		sCalibInfos->iFailedOpt = i;
		sCalibInfos->dFailedSigma2 = dSigma2;
	}

	return err;
}

static double	generic_midat_LocalVolInteg(double			t1,
											double			t2,
											double			T1,
											double			T2,
											GENMIDAT_MODEL	sModel)
{
double	res;

	res = t2 - t1;

	res += sModel->dAlpha * sModel->dRho
		* (exp(-sModel->dGamma * (T1 - t2)) + exp(-sModel->dGamma * (T2 - t2)) 
		- exp(-sModel->dGamma * (T1 - t1)) - exp(-sModel->dGamma * (T2 - t1))) / sModel->dGamma;

	res += sModel->dAlpha * sModel->dAlpha
		* (exp(-sModel->dGamma * (T1 + T2 - 2.0 * t2)) - exp(-sModel->dGamma * (T1 + T2 - 2.0 * t1))) 
		/ (2.0 * sModel->dGamma);

	return res;
}

/* All the Functions for calibration */
/* ********************************* */

Err	GenMidatCalib_GetTarget(void				*Inst_,
							void				*Params_,
							void				*Model_,
							CALIBGEN_PARAMS		CalibConsts,
							double				*target)
{
GenMidatCalib_INST		Inst	= (GenMidatCalib_INST) Inst_;
	
	*target = Inst->sLastInst->dShortVol;	

	return NULL;
}

Err	GenMidatCalib_GetFirstGuess(void				*Model_,
								void				*Params_,
								int					index_param,
								double				target,
								double				*param1)
{
GenMidatCalib_MODEL		Model	= (GenMidatCalib_MODEL) Model_;

	*param1 = Model->sModel->dLambda[index_param];

	return NULL;
}

Err	GenMidatCalib_GetSecondGuess(void				*Model_,
								void				*Params_,
								int					index_param,
								double				param1,
								double				price1,
								double				target,
								double				*param2)
{	
	*param2 = param1 * 1.1;	

	return NULL;
}

Err	GenMidatCalib_GetLimitAndLastParam(	void				*Model_,
										CALIBGEN_PARAMS	CalibConsts,
										void				*Params_,
										int					index_param,
										double				*last_param,
										double				*limit_down,
										double				*limit_up)
{
	*limit_down = 0.01;
	*limit_up = 1.0e10;

	return NULL;
}

Err	GenMidatCalib_BumpParam(void				*Model_,
							void				*Params_,
							int					index_param,
							double				param1)
{
GenMidatCalib_MODEL		Model	= (GenMidatCalib_MODEL) Model_;

	Model->dLambda = param1;

	return NULL;
}

Err	GenMidatCalib_SetParam(	void				*Model_,
							CALIBGEN_PARAMS		CalibConsts,
							void				*Params_,
							int					index_param,
							double				param1)
{
GenMidatCalib_MODEL		Model	= (GenMidatCalib_MODEL) Model_;
GenMidatCalib_PARAMS	Params	= (GenMidatCalib_PARAMS) Params_;
GenMidatCalib_INST		Inst	= Params->AllInst[index_param];

Err	err = NULL;

	if (!CalibConsts->do_calib)
	{	
		param1 = Model->sModel->dLambda[index_param];
	}

	/* We set the Lambda */
	Model->dLambda = param1;
	Model->sModel->dLambda[index_param] = param1;
	Model->sModel->dLambda2[index_param] = param1 + Model->sModel->dGamma;

	/* Update Beta */
	Model->dBeta = Model->dLastBeta * exp(-Model->dLambda * Inst->dDt);
	Model->sModel->dBeta[index_param] = Model->dBeta;
	Model->sModel->dBeta2[index_param] = Model->dBeta * Inst->dExpFact1;	
			
	/* Calibrate the volatility */
	
	/* Update Factors */
	Inst->dFact0 = (1.0 / (Model->dBeta * Model->dBeta) - 1.0 / (Model->dLastBeta * Model->dLastBeta))
		/ (2.0 * Model->dLambda);

	Inst->dFact1 = (1.0 / (Inst->dExpFact1 * Model->dBeta * Model->dBeta) - 1.0 / (Inst->dLastExpFact1 * Model->dLastBeta * Model->dLastBeta))
		/ (2.0 * Model->dLambda + Model->sModel->dGamma);

	Inst->dFact2 = (1.0 / (Inst->dExpFact2 * Model->dBeta * Model->dBeta) - 1.0 / (Inst->dLastExpFact2 * Model->dLastBeta * Model->dLastBeta))
		/ (2.0 * (Model->dLambda + Model->sModel->dGamma));

	/* Calibrate the volatility */
	
	Model->dSigma2 = Inst->dLongVariance / (Model->dBeta * Model->dBeta)
		- Model->dVarX 
		- 2.0 * Model->dCoVar * Inst->dExpFact1 
		- Model->dVarY * Inst->dExpFact2;

	Model->dSigma2 /= Inst->dFact0 
		+ 2.0 * Model->sModel->dRhoAlpha * Inst->dExpFact1 * Inst->dFact1
		+ Model->sModel->dAlpha2 * Inst->dExpFact2 * Inst->dFact2;

	if (Model->dSigma2 < 0.0)
	{
		err = serror("Cannot calibrate sigma on option %d", index_param + 1);
		return err;
	}

	/* Save Local volatility and update variance */
	Model->sModel->dSigma[index_param] = Model->dSigma2 * Inst->dFact0;
	Model->sModel->dCorrel[index_param] = Model->sModel->dRhoAlpha * Model->dSigma2 * Inst->dFact1;
	Model->sModel->dSigma2[index_param] = Model->sModel->dAlpha2 * Model->dSigma2 * Inst->dFact2;
	
	Model->dVarX = Model->dVarX + Model->sModel->dSigma[index_param];
	Model->dCoVar = Model->dCoVar + Model->sModel->dCorrel[index_param];
	Model->dVarY = Model->dVarY + Model->sModel->dSigma2[index_param];

	Model->sModel->dSigma[index_param] = sqrt(Model->sModel->dSigma[index_param] / Inst->dDt);
	Model->sModel->dSigma2[index_param] = sqrt(Model->sModel->dSigma2[index_param] / Inst->dDt);
	Model->sModel->dCorrel[index_param] /= Model->sModel->dSigma[index_param] * Model->sModel->dSigma2[index_param] * Inst->dDt;

	/* Shift Beta */
	Model->dLastBeta = Model->dBeta;

	return NULL;
}

Err	GenMidatCalib_ExtrapolParam(void				*Model_,
								void				*Params_,
								int					index_param)
{
	return NULL;
}

Err	GenMidatCalib_UpdateConstsAfterParam(	void				*Inst_,
											void				*InstConst_,
											void				*Params_,
											void				*Model_,
											CALIBGEN_PARAMS		CalibConsts)
{
GenMidatCalib_PARAMS	Params	= (GenMidatCalib_PARAMS) Params_;
GenMidatCalib_MODEL		Model	= (GenMidatCalib_MODEL) Model_;
GenMidatCalib_INST		Inst	= (GenMidatCalib_INST) Inst_;

Err	err = NULL;

	/* Update first params */
	if (Inst->iIndex == 1)
	{
		/* Reinitialise all the cons*/
		
		Model->dLastBeta = 1.0;
		Model->dVarX = 0.0;
		Model->dVarY = 0.0;
		Model->dCoVar = 0.0;		
		
		err = GenMidatCalib_SetParam(	Model_,
										CalibConsts,
										Params_,
										0,
										Model->dLambda);
	}

	/* Update Beta */
	if (CalibConsts->do_calib)
	{
		Model->dBeta = Model->dLastBeta * exp(-Model->dLambda * Inst->dDt);	
	}	
			
	return err;
}

Err	GenMidatCalib_PriceInst(void				*Inst_,
							void				*InstConst_,
							void				*Params_,
							void				*Model_,
							double				*InstPrice)
{
GenMidatCalib_PARAMS	Params	= (GenMidatCalib_PARAMS) Params_;
GenMidatCalib_MODEL		Model	= (GenMidatCalib_MODEL) Model_;
GenMidatCalib_INST		NxtInst	= (GenMidatCalib_INST) Inst_;
GenMidatCalib_INST		Inst	= NxtInst->sLastInst;

double	vol;
double	diffbet1, diffbet2;

	diffbet1 = Model->dBeta - Model->dLastBeta;
	diffbet2 = Model->dBeta * NxtInst->dExpFact1  - Model->dLastBeta * Inst->dExpFact1; 

	vol = diffbet1 * diffbet1 * Model->dVarX;
	vol += 2.0 * diffbet1 * diffbet2 * Model->dCoVar;
	vol += diffbet2 * diffbet2 * Model->dVarY;
	
	*InstPrice = sqrt(vol / Inst->dTime);

	return NULL;
}

Err generic_midat_calib_two_factor(	/* Input */
									int						iNbExe,
									int						*iCalibLong,
									double					*dLongVols,
									int						*iCalibShort,
									double					*dShortVols,

									/* Model */									
									GENMIDAT_MODEL			sModel,

									/* Parameters */
									GENMIDAT_CALIBPARAMS	sParams,

								   /* Optional Infos */
									GENMIDAT_CALIBINFOS		sCalibInfos)
{
Err		err = NULL;
double	ratio;
int		i;
CALIBGEN_PARAMS			CalibParams		= NULL;
CALIBFUNCTIONS			AllFunctions	= NULL;
GenMidatCalib_INST		*AllInst		= NULL;
GenMidatCalib_MODEL		CalibModel		= NULL;
GenMidatCalib_PARAMS	Params			= NULL;
			
	/* Memory allocations */
	CalibParams = calloc(1, sizeof(CALIBGEN_Params));	
	AllFunctions = calloc(1, sizeof(CalibFunctions));
	AllInst = calloc(iNbExe, sizeof(GenMidatCalib_INST));
	CalibModel = calloc(1, sizeof(genmidat_calib_model));
	Params = calloc(1, sizeof(genmidat_calib_params));

	if (!CalibParams || !AllFunctions || !AllInst || !CalibModel || !Params)
	{
		err = "Memory allocation faillure in generic_midat_calib_two_factor_total";
		goto FREE_RETURN;
	}

	/* First, Guess */
	if (sParams->iCalibLambda && iNbExe > 1)
	{				
		/* Special case for the first Lambda */
		ratio = 1.0 - dShortVols[0] / dLongVols[0];

		if (ratio > 0.0)
		{
			sModel->dBeta[0] = pow(ratio, sModel->dTimes[0] / (sModel->dTimes[1] - sModel->dTimes[0]));
			sModel->dBeta[1] = pow(sModel->dBeta[0], sModel->dTimes[1] / sModel->dTimes[0]);
		}
		else
		{
			err = serror("Cannot calibrate lambda on option %d", 1);
			goto FREE_RETURN;
		}

		for (i=2; i<iNbExe; i++)
		{
			ratio = 1.0 - dShortVols[i-1] / dLongVols[i-1];

			if (ratio > 0.0)
			{
				sModel->dBeta[i] = sModel->dBeta[i-1] * ratio;
			}
			else
			{
				err = serror("Cannot calibrate lambda on option %d", i);
				goto FREE_RETURN;
			}
		}
	}

	genmidat_convert_beta_into_lambda(sModel);

	/* Fill the structures */

	/* Initialise Newton params */
	err = Initialise_CalibParams(	1,
									0.00001,
									7,
									1,
									sParams->iCalibLambda,
									1,
									0,
									0.0,
									CalibParams);

	if (err) goto FREE_RETURN;

	/* Fill Instruments */
	for (i=0; i<iNbExe; i++)
	{
		AllInst[i] = calloc(1, sizeof(genmidat_calib_inst));

		if (!AllInst[i])
		{
			err = "Memory allocation faillure in generic_midat_calib_two_factor_total";
			goto FREE_RETURN;
		}

		AllInst[i]->iIndex = i;
		AllInst[i]->dTime = sModel->dTimes[i];
		
		if (i > 0)
		{
			AllInst[i]->dDt = sModel->dTimes[i] - sModel->dTimes[i-1];
			AllInst[i]->sLastInst = AllInst[i-1];
			AllInst[i]->dLastExpFact1 = AllInst[i-1]->dLastExpFact1;
			AllInst[i]->dLastExpFact2 = AllInst[i-1]->dLastExpFact2;
		}
		else
		{
			AllInst[i]->dDt = sModel->dTimes[i];
			AllInst[i]->sLastInst = NULL;
			AllInst[i]->dLastExpFact1 = 1.0;
			AllInst[i]->dLastExpFact2 = 1.0;
		}

		AllInst[i]->dExpFact1 = exp(-sModel->dGamma * sModel->dTimes[i]);
		AllInst[i]->dExpFact2 = AllInst[i]->dExpFact1 * AllInst[i]->dExpFact1;
		
		AllInst[i]->dLongVariance = dLongVols[i] * dLongVols[i] * sModel->dTimes[i];
		AllInst[i]->dShortVol = dShortVols[i];
	}

	/* Fill the Model */
	CalibModel->sModel = sModel;

	if (sParams->iCalibLambda)
	{
		CalibModel->dLastBeta = 1.0;
	}
	else
	{
		CalibModel->dLastBeta = sModel->dBeta[0];
	}

	CalibModel->dVarX = 0.0;
	CalibModel->dVarY = 0.0;
	CalibModel->dCoVar = 0.0;

	/* Fill the functions */
	AllFunctions->GetTarget = GenMidatCalib_GetTarget;
	AllFunctions->BumpParam = GenMidatCalib_BumpParam;
	AllFunctions->ExtrapolParam = GenMidatCalib_ExtrapolParam;
	AllFunctions->GetFirstGuess = GenMidatCalib_GetFirstGuess;
	AllFunctions->GetLimitAndLastParam = GenMidatCalib_GetLimitAndLastParam;
	AllFunctions->GetSecondGuess = GenMidatCalib_GetSecondGuess;
	AllFunctions->PriceInst = GenMidatCalib_PriceInst;
	AllFunctions->SetParam = GenMidatCalib_SetParam;
	AllFunctions->UpdateConstsAfterParam = GenMidatCalib_UpdateConstsAfterParam;

	/* Fill Params */
	Params->AllInst = AllInst;

	/* Launch the Newton */	
	err = CalibrateParamTS(	1,
							iNbExe-1,
							AllInst,
							AllInst,
							Params,
							CalibModel,
							CalibParams,
							AllFunctions);

	if (err) goto FREE_RETURN;	
	
FREE_RETURN:

	if (CalibParams)
	{
		Free_CalibParams(CalibParams);
		free(CalibParams);
	}	

	if (AllFunctions) free(AllFunctions);

	if (AllInst)
	{
		for (i=0; i<iNbExe; i++)
		{
			if (AllInst[i]) free (AllInst[i]);
		}

		free (AllInst);
	}

	if (CalibModel) free (CalibModel);
	if (Params) free (Params);

	return err;
}

