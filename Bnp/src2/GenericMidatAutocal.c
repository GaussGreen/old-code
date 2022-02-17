
#include "math.h"
#include "opfnctns.h"
#include "GenericMidatAutocal.h"
#include "Fx3FUtils.h"

void genmidat_set_autocalparams_default(GENMIDAT_AUTOCALPARAMS sParams)
{
	sParams->iOneTimeCallIndex = 0;
	sParams->iCalcExeBound = 0;
	sParams->iOneTimeShortLong = 1;

	sParams->lFindBestOptim = 0;
}

void genmidat_init_autocalinfos(GENMIDAT_AUTOCALINFO sInfos)
{
	sInfos->sCalibInfos = NULL;
	
	sInfos->iNbExe = 0;
	sInfos->dExeBound = NULL;
	sInfos->dCoefLin = NULL;
}

Err genmidat_allocate_autocalinfos(	int						iNbExe,
									GENMIDAT_AUTOCALPARAMS	sParams,
									GENMIDAT_AUTOCALINFO	sInfos)
{
Err err = NULL;

	sInfos->iCalcExeBound = sParams->iCalcExeBound;	

	if (sInfos->iCalcExeBound)
	{
		sInfos->iNbExe = iNbExe;
		sInfos->dExeBound = calloc(iNbExe, sizeof(double));
		sInfos->dCoefLin = dmatrix(0, iNbExe-1, 0, 1);

		if (!sInfos->dExeBound || !sInfos->dCoefLin)
		{
			err = "Memory allocation faillure in genmidat_allocate_autocalinfos";
			return err;
		}
	}
	else
	{
		sInfos->iNbExe = 0;
		sInfos->dExeBound = NULL;
		sInfos->dCoefLin = NULL;
	}

	return err;
}

void genmidat_free_autocalinfos(GENMIDAT_AUTOCALINFO	sInfos)
{
	if (sInfos)
	{
		if (sInfos->dExeBound) free (sInfos->dExeBound);
		if (sInfos->dCoefLin) free_dmatrix(sInfos->dCoefLin, 0, sInfos->iNbExe-1, 0, 1);
	}
}

Err GenericMidat_OneFactorPayoff_PDE(/* Event */
									double	dTime,

									/* Parameters */
									void	*sFuncParm,																										
									
									/* Gride data	*/
									long	lStartIndex,
									long	lEndIndex,											
									double	*dXt,
									GENMIDAT_RECONS	sReconsParams,

									/* Vector of results to be updated */
									int		iNbProd,
									double	**dProdVal)
{
GENMIDAT_PAYOFFPARAMS sPayoffParam;
int		i, called;
double	dForward;
double	dCoef;

	sPayoffParam = (GENMIDAT_PAYOFFPARAMS) sFuncParm;
	
	dCoef = sPayoffParam->dBeta1 * sReconsParams->dCoefX1;	

	if (!sPayoffParam->sInfos->iCalcExeBound)
	{
		for (i=lStartIndex; i<=lEndIndex; i++)
		{
			/* Reconstruction */
			dForward = sPayoffParam->dForward + dCoef * dXt[i];

			if (dForward > dProdVal[i][0])
			{
				dProdVal[i][0] = dForward;
			}
		}
	}
	else
	{
		called = 0;

		for (i=lStartIndex; i<=lEndIndex; i++)
		{
			/* Reconstruction */
			dForward = sPayoffParam->dForward + dCoef * dXt[i];

			if (dForward > dProdVal[i][0])
			{
				dProdVal[i][0] = dForward;

				if (!called)
				{
					sPayoffParam->sInfos->dExeBound[sPayoffParam->iIndex] = dForward;
					called = 1;
				}
			}
		}
	}	

	return NULL;
}

Err GenericMidat_TwoFactorPayoff_PDE(/* Event */
									double	dTime,

									/* Parameters */
									void	*sFuncParm,																										
									
									/* Gride data	*/
									long	lStartX,
									long	lEndX,
									long	lStartZ,
									long	lEndZ,

									double	*dXt,
									double	*dZt,

									GENMIDAT_RECONS	sReconsParams,
															
									/* Vector of results to be updated */
									int		iNbProd,
									double	***dProdVal)
{
GENMIDAT_PAYOFFPARAMS sPayoffParam;
int		i, j;
double	dForwardi, dForward;
double	coeftotX, coeftotZ;

	sPayoffParam = (GENMIDAT_PAYOFFPARAMS) sFuncParm;

	coeftotX = sPayoffParam->dBeta1 * sReconsParams->dCoefX1 + sPayoffParam->dBeta2 * sReconsParams->dCoefY1;
	coeftotZ = sPayoffParam->dBeta1 * sReconsParams->dCoefX2 + sPayoffParam->dBeta2 * sReconsParams->dCoefY2;

	for (i=lStartX; i<=lEndX; i++)
	{
		dForwardi = sPayoffParam->dForward + coeftotX * dXt[i];

		for (j=lStartZ; j<=lEndZ; j++)
		{
			/* Reconstruction */
			dForward = dForwardi + coeftotZ * dZt[j];

			if (dForward > dProdVal[i][j][0])
			{
				dProdVal[i][j][0] = dForward;
			}
		}
	}

	return NULL;
}

Err GenericMidat_OneFactorPayoff_MC(/* Event */
									double	dTime,

									/* Parameters */
									void	*sFuncParm,																										
									
									/* Gride data	*/
									long	iNumPaths,
									double	*dXt,
									GENMIDAT_RECONS	sReconsParams,

									/* Vector of results to be updated */
									int		iNbProd,
									double	**dProdVal)
{
GENMIDAT_PAYOFFPARAMS sPayoffParam;
int		i;
double	dForward;
double	dCoef;

	sPayoffParam = (GENMIDAT_PAYOFFPARAMS) sFuncParm;
	
	dCoef = sPayoffParam->dBeta1 * sReconsParams->dCoefX1;	
			
	for (i=0; i<iNumPaths; i++)
	{
		/* Reconstruction */
		dForward = sPayoffParam->dForward + dCoef * dXt[i];

		dProdVal[i][0] = dForward;
	}	

	return NULL;
}

Err GenericMidat_TwoFactorPayoff_MC(/* Event */
									double	dTime,

									/* Parameters */
									void	*sFuncParm,																										
									
									/* Gride data	*/
									long	iNumPaths,
									double	*dXt,
									double	*dYt,
									GENMIDAT_RECONS	sReconsParams,

									/* Vector of results to be updated */
									int		iNbProd,
									double	**dProdVal)
{
GENMIDAT_PAYOFFPARAMS sPayoffParam;
int		i;

	sPayoffParam = (GENMIDAT_PAYOFFPARAMS) sFuncParm;		
			
	for (i=0; i<iNumPaths; i++)
	{
		/* Reconstruction */		
		dProdVal[i][0] = sPayoffParam->dForward + dXt[i];
		dProdVal[i][1] = sPayoffParam->dForward2 + dYt[i];
	}	

	return NULL;
}

Err	GenericMidatAutocal(int						iNbExe,
						double					*dLongOptions,
						double					*dShortOptions,
						int						*iIsCallDate,

						/* Model */
						GENMIDAT_MODEL			sModel,

						/* Parameters */
						GENMIDAT_AUTOCALPARAMS	sAutocalParams,
						GENMIDAT_CALIBPARAMS	sCalibParams,
						GENMIDAT_PDEPAMS		sPDEParams,

						/* Outputs */
						double					*MultiCall,
						GENMIDAT_AUTOCALINFO	sInfos)
{
Err	err = NULL;
double	*dTimes			= NULL,
		**dMCProdVal	= NULL;

int		*iOptimise	= NULL;

MCEBPARAMS sMCEBParams = NULL;

int		i, k, iIndexEvent, iNbEvent, iNbStep;

GENMIDAT_PAYOFFPARAMS	*sPayoffParams	= NULL;
int						*iEvalEvent		= NULL;

	/* FIRST STEP: Calibration */
	/* *********************** */

	err = generic_midat_calib(	iNbExe,
								dLongOptions,
								dShortOptions,
								sModel,
								sCalibParams,
								sInfos->sCalibInfos);

	if (err) goto FREE_RETURN;

	/* SECOND STEP: PDE parameters */
	/* *************************** */

	iNbStep = iNbExe;
	dTimes = calloc(iNbExe, sizeof(double));

	if (!dTimes)
	{
		err = "Memory allocation faillure in GenericMidatAutocal";
		goto FREE_RETURN;
	}

	memcpy(dTimes, sModel->dTimes, sModel->iNbVols * sizeof(double));

	if (sPDEParams->iPricingMethod == 0)
	{
		/* PDE */
		err = fill_time_vector (	&dTimes, 
									&iNbStep, 
									0,
									NULL,
									0, 
									NULL, 
									sPDEParams->lNbTime);

		if (err) goto FREE_RETURN;
	}

	sPayoffParams = calloc(iNbStep, sizeof(GENMIDAT_PAYOFFPARAMS));
	iEvalEvent = calloc(iNbStep, sizeof(int));

	if (!sPayoffParams || !iEvalEvent)
	{
		err = "Memory allocation faillure (2) in GenericMidatAutocal";
		goto FREE_RETURN;
	}

	iIndexEvent = iNbExe - 1;
	iNbEvent = 0;
	
	for (i=iNbStep-1; i>=0; i--)
	{
		if (iIndexEvent >= 0 && fabs (dTimes[i] - sModel->dTimes[iIndexEvent]) < 1.0e-04)
		{
			if ((sAutocalParams->iOneTimeCallIndex == 0 || iIndexEvent ==  sAutocalParams->iOneTimeCallIndex - 1)
				&& (!iIsCallDate || iIsCallDate[iIndexEvent]))
			{
				/* Fill the Event */
				iEvalEvent[i] = 1;
				iNbEvent++;

				sPayoffParams[i] = calloc(1, sizeof(genmidat_payoffparams));

				if (!sPayoffParams[i])
				{
					err = "Memory allocation faillure (3) in GenericMidatAutocal";
					goto FREE_RETURN;
				}

				sPayoffParams[i]->iIndex = iIndexEvent;
				sPayoffParams[i]->iEventIndex = iNbEvent;

				sPayoffParams[i]->dForward = sModel->dForward[iIndexEvent];				
				sPayoffParams[i]->dBeta1 = sModel->dBeta[iIndexEvent];								

				if (iIndexEvent < sModel->iNbVols - 1)
				{
					sPayoffParams[i]->dForward2 = sModel->dForward[iIndexEvent] - sModel->dForward[iIndexEvent+1];
				}
				else
				{
					sPayoffParams[i]->dForward2 = sModel->dForward[iIndexEvent];				
				}

				if (sModel->iNbFactor == 2)
				{
					sPayoffParams[i]->dBeta2 = sModel->dBeta2[iIndexEvent];
				}

				if (sAutocalParams->iOneTimeCallIndex > 0 && sAutocalParams->iOneTimeShortLong < 0.5 && iIndexEvent < iNbExe - 1)
				{
					sPayoffParams[i]->dForward -= sModel->dForward[iIndexEvent+1];
					sPayoffParams[i]->dBeta1 -= sModel->dBeta[iIndexEvent+1];

					if (sModel->iNbFactor == 2)
					{
						sPayoffParams[i]->dBeta2 -= sModel->dBeta2[iIndexEvent+1];
					}
				}

				sPayoffParams[i]->sPdeParams = sPDEParams;

				sPayoffParams[i]->sInfos = sInfos;
			}
			else
			{
				iEvalEvent[i] = 0;
				sPayoffParams[i] = NULL;
			}

			iIndexEvent--;
		}
		else
		{
			iEvalEvent[i] = 0;
			sPayoffParams[i] = NULL;
		}
	}	


	/* THIRD STEP: Call the Pricer */
	/* *************************** */

	switch (sPDEParams->iPricingMethod)
	{
		case 0:
		{
			/* PDE */
			if (sModel->iNbFactor == 1)
			{
				err = GenericMidat_OneFactor_PDE(	iNbStep,
													dTimes,
													sModel,
													sPayoffParams,
													iEvalEvent,
													GenericMidat_OneFactorPayoff_PDE,
													sPDEParams,
													1,
													MultiCall);
			}
			else if (sModel->iNbFactor == 2)
			{
				err = GenericMidat_TwoFactor_PDE(	iNbStep,
													dTimes,
													sModel,
													sPayoffParams,
													iEvalEvent,
													GenericMidat_TwoFactorPayoff_PDE,
													sPDEParams,
													1,
													MultiCall);
			}

			if (err) goto FREE_RETURN;

			break;
		}

		case 1:
		{
			iOptimise = calloc(iNbStep, sizeof(int));
			sMCEBParams = calloc(1, sizeof(MCEBParams));
			dMCProdVal = dmatrix(0, sModel->iNbFactor + iNbStep, 0, 2 + sModel->iNbFactor);

			if (!iOptimise || !sMCEBParams || !dMCProdVal)
			{
				err = "Memory allocation faillure in Generic MidatAutocal";
				goto FREE_RETURN;
			}

			mceb_set_default_params(sMCEBParams);

			/* Set MCEB Params */
			sMCEBParams->iCallCurrent = 1;
			sMCEBParams->iIsKO = 0;

			if (sAutocalParams->iOneTimeCallIndex == 0 || sAutocalParams->iOneTimeShortLong || sModel->iNbFactor == 1)
			{
				sMCEBParams->iColPay = 0;
			}
			else
			{
				sMCEBParams->iColPay = 1;
			}

			sMCEBParams->iColBound = 0;
			sMCEBParams->iMultiIndex = 1;
			sMCEBParams->iNbIndex = sModel->iNbFactor;
			sMCEBParams->iRemoveLastOnLast = 1;
			sMCEBParams->iDoInfos = 0;
			sMCEBParams->iFindBestOptim = sAutocalParams->lFindBestOptim;
			
			err = mceb_allocate_params(sMCEBParams, iNbStep);
			if (err) goto FREE_RETURN;

			/* Check that the Exe Bounds are allocated */
			if (!sAutocalParams->iCalcExeBound)
			{
				sAutocalParams->iCalcExeBound = 1;
				err = genmidat_allocate_autocalinfos(	iNbStep,
														sAutocalParams,
														sInfos);

				if (err) goto FREE_RETURN;
			}

			for (i=0; i<iNbStep; i++)
			{
				if (iEvalEvent[i])
				{
					iOptimise[i] = 1;
				}
				else
				{
					iOptimise[i] = 0;
				}
			}

			/* MC */
			if (sModel->iNbFactor == 1)
			{
				err = GenericMidat_OneFactor_MC(	iNbStep,
													dTimes,
													sModel,
													sPayoffParams,
													GenericMidat_OneFactorPayoff_MC,
													sPDEParams,
													1,
													iOptimise,
													sMCEBParams,
													sModel->iNbFactor,
													dMCProdVal);
			}
			else if (sModel->iNbFactor == 2)
			{
				err = GenericMidat_TwoFactor_MC(	iNbStep,
													dTimes,
													sModel,
													sPayoffParams,
													GenericMidat_TwoFactorPayoff_MC,
													sPDEParams,
													1,
													iOptimise,
													sMCEBParams,
													sModel->iNbFactor,
													dMCProdVal);				
			}

			if (err) goto FREE_RETURN;

			/* Recopy Barrier / CoefLin for the moment */
			for (i=0; i<iNbStep; i++)
			{
				dMCProdVal[i][2] = sMCEBParams->dBarrier[i];

				for (k=0; k<sMCEBParams->iNbIndex; k++)
				{
					dMCProdVal[i][3+k] = sMCEBParams->dCoefLin[i][k+1];
				}
			}

			break;
		}

		default:
		{
			err = "Pricing Method not recognised in Generic Midat Autocal";
			goto FREE_RETURN;
		}
	}

	/* FOURTH STEP: Discount and Return */
	/* ******************************** */

	switch (sPDEParams->iPricingMethod)
	{
		case 0:
		{
			*MultiCall *= sModel->dNumeraire;
			break;
		}

		case 1:
		{
			*MultiCall = dMCProdVal[sModel->iNbFactor][0] * sModel->dNumeraire;

			for (i=0; i<iNbStep; i++)
			{
				sInfos->dExeBound[i] = dMCProdVal[i][2] * sModel->iNbVols;
				sInfos->dCoefLin[i][0] = dMCProdVal[i][3];

				if (sModel->iNbFactor == 2)
				{
					sInfos->dCoefLin[i][1] = dMCProdVal[i][4];
				}
				else
				{
					sInfos->dCoefLin[i][1] = 0.0;
				}
			}

			break;
		}
	}


FREE_RETURN:

	if (dTimes) free (dTimes);

	if (sPayoffParams)
	{
		for (i=0; i<iNbStep; i++)
		{
			if (sPayoffParams[i]) free (sPayoffParams[i]);
		}

		free (sPayoffParams);
	}

	if (iEvalEvent) free (iEvalEvent);

	if (iOptimise) free(iOptimise);
	if (sMCEBParams)
	{
		mceb_free_params(sMCEBParams);
		free(sMCEBParams);
	}

	if (dMCProdVal) free_dmatrix(dMCProdVal, 0, sModel->iNbFactor + iNbStep, 0, 2 + sModel->iNbFactor);

	return err;
}




