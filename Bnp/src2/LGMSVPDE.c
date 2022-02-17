/* ==========================================================================
   FILE_NAME:	LGMSVPDE.c

   PURPOSE:		PDE implementation of the 1 Factor LGM model with stochastic vol.
	
   DATE:		02/25/03
   
   AUTHOR:		P.A. and L.C
   ========================================================================== */

#include "LGMSVPDE.h"
#include "LGMSVGrfn.h"
#include "LGMSVUtil.h"

/* For swaping two variables */
#define LGMSVPDESwap(Type,A,B) {Type C; C=A; A=B; B=C;}

#define	LGMSV_MINPSIVOL		0.0001
#define	LGMSV_MINVTVOL		0.0001
#define	LGMSV_MINORTHOVOL	0.0001

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEExpectations														  

	Calculation of expectations for the definition of the grids 		
		in f(t,Tstar)-f(0,Tstar) - u(t)*dRho/dAlpha*(Vt-1)
		in Zt = Vt 
		in Psit=Phit*exp(-2*dLambdaX*(Tstar-t))

	The Model is under QTstar

	df(t,Tstar) = u(t) * sqrt(Vt) * dW1(t)
	dPsit = u(t)^2*Vt dt
	dVt = -dLambdaV*(Vt-1)dt + dAlpha * sqrt(Vt)*dW2(t)
	<dW1(t),dW2(t)> = dRho * dt

	where 
		u(t) is piecewise constant

	Orthogonalisation

		Xt = f(t,Tstar)-f(0,Tstar) - u(t)*dRho/dAlpha*(Vt-1)
			is such that <dXt,dW2(t)>=0

		Zt = Vt 
	

   ------------------------------------------------------------------------------------------------ */
static Err LGMSVPDEExpectations_UtPieceWise(
					/* Inputs */
					 int	iNbTime,			
					 double	*dTime,
					
					 /* For f(t,Tstar) diffusion and Xt*/
					 double dLambdaX,
					 double dSigma,
					 double dTstar,
					 
					 /* definition of u(t) piecewise constant */
					 double *SigTime,
					 double	*Sig,
					 int	iNbSigTime,

					 /* For Vt diffusion with TS on same dates as the volatility */
					 double *dLambdaV,
					 double	*dLvlV,
					 double *dAlphaV,
					 double	*dRho,

					 /* Numerical parameters */
					 double	dMinTime,

					 /* Outputs */
					 double *Sigt,			/* Value of u(t) = sig(t) for all t in dtime */
					 double	*dLambdaVt,		/* Value of lameps(t) for all t in dtime */
					 double	*dLvlVt,		/* Value of mean-revert level */
					 double	*dAlphaVt,		/* Value of alphaeps(t) for all t in dtime */
					 double	*dRhot,			/* Value of rhoeps(t) for all t in dtime */
					 
					 double *PsitMean,
					 double *PsitStd,
					 
					 double *XtMean,
					 double *XtStd,

					 double *ZtMean,
					 double *ZtStd,
					 double *ImpliedRhoTimesStdVt)
{	
/* Declaration of locals variables */
int iNumTime;
double dDt;

double			*lambda			= NULL;

LGMSVSolFunc1	*FuncVarft		= NULL,
				*FuncVt2Mean	= NULL,
				*FuncVtPsitMean	= NULL,
				*FuncPsit2Mean	= NULL,
				*FuncVtftMean	= NULL;

double	dVarft,
		dNewVarft,
		dVt2Mean,
		dNewVt2Mean,
		dPsitMean,
		dNewPsitMean,
		dVtPsitMean,
		dNewVtPsitMean,
		dPsit2Mean,
		dNewPsit2Mean,
		dVtftMean,
		dNewVtftMean;

double dRatio, dDtRatio;
double dTemp, dTempVtFtMean;

double dFinalTime, dLastTime, dNewTime, dLastPeriodTime, dLocalMinTime;
int	   iVolIndex, iNbPeriod, iPeriodIndex, iLocalIndexTime;


double du2;
Err	   err = NULL;

	/* All allocation errors */
	
	lambda = calloc(100, sizeof(double));
	FuncVarft = calloc(1, sizeof(LGMSVSolFunc1));
	FuncVt2Mean = calloc(1, sizeof(LGMSVSolFunc1));
	FuncVtPsitMean = calloc(1, sizeof(LGMSVSolFunc1));
	FuncPsit2Mean = calloc(1, sizeof(LGMSVSolFunc1));
	FuncVtftMean = calloc(1, sizeof(LGMSVSolFunc1));

	if (!lambda || !FuncVarft || !FuncVt2Mean || 
		!FuncVtPsitMean || !FuncPsit2Mean || !FuncVtftMean)
	{
		err = "Memory allocation faillure in LGMSVPDEExpectations_UtPieceWise";
		goto FREE_RETURN;
	}

	/* Initialisation of the LGMSVSolFunc */
		
	/* Definition of the function Var(f(t,Tstar)-f(0,Tstar)) */
	/* dE(f^2)/dt = ut^2 */
	FuncVarft->dLambda = 0.0;
	FuncVarft->iNbFunction = 1;
	
	FuncVarft->bIsft1[0] = 1;

	/* Definition of the function E(Vt^2) */
	/* dE(Vt^2)/dt = (dAlphaV*dAlphaV+2*dLvlt)-2*dLambdaVt*E(Vt^2) */	
	/*FuncVt2Mean->dLambda = 2*dLambdaV;*/
	FuncVt2Mean->iNbFunction = 1;

	FuncVt2Mean->bIsft1[0] = 1;
	/*FuncVt2Mean->a[0] = dAlphaV*dAlphaV+2*dLvlV;*/

	/* Definition of the function E(Psit) */
	/* dE(Psit)/dt = ut^2 = dE(f^2)/dt */
	/* => E(Psit) = Var(f(t,Tstar)-f(0,Tstar)) */	

	/* Definition of the function E(VtPsit) */
	/* dE(VtPsit)/dt = dLvlV*E(Psit) + ut^2 *E(Vt^2) - dLambdaV*E(VtPsit) */
	/*FuncVtPsitMean->dLambda = dLambdaV;*/
	FuncVtPsitMean->iNbFunction = 2;

	FuncVtPsitMean->bIsft1[0] = 0;
	FuncVtPsitMean->pft[0] = FuncVt2Mean;
	
	FuncVtPsitMean->bIsft1[1] = 0;
	/*FuncVtPsitMean->a[1] = dLvlV;*/
	FuncVtPsitMean->pft[1] = FuncVarft;
	
	/* Definition of the function E(Psit^2) */
	/* dE(Psit^2)/dt = 2 * ut^2 * E(VtPsit) */
	FuncPsit2Mean->dLambda = 0;
	FuncPsit2Mean->iNbFunction = 1;

	FuncPsit2Mean->bIsft1[0] = 0;
	FuncPsit2Mean->pft[0] = FuncVtPsitMean;

	/* Definition of the function E(Vt*(f(t,Tstar)-f(0,Tstar))) */
	/* dE(Vt*ft)/dt = dAlphaV*dRho*ut-dLambdaV*E(Vt*ft) */
	/*FuncVtftMean->dLambda = dLambdaV;*/
	FuncVtftMean->iNbFunction = 1;

	FuncVtftMean->bIsft1[0] = 1;

	/* Filling the Sigt array 
	   piecewise interpolation from the Sig array on SigTime 
	   u(t) is right continuous									*/		
	
	/* Init Value at t=0 */
	Sigt[0] = Sig[0];
	dLambdaVt[0] = dLambdaV[0];
	dLvlVt[0] = dLvlV[0];
	dAlphaVt[0] = dAlphaV[0];
	dRhot[0] = dRho[0];

	dVarft = 0.0;
	dVt2Mean = 1.0;
	dPsitMean = 0.0;
	dVtPsitMean = 0.0;
	dPsit2Mean = 0.0;
	dVtftMean = 0.0;
	
	PsitMean[0] = 0.0;
	PsitStd[0] = 0.0;
	ZtMean[0] = 1.0;
	ZtStd[0] = 0.0;
	XtMean[0] = 0.0;
	XtStd[0] = 0.0;
	ImpliedRhoTimesStdVt[0] = 0.0;

	dLastTime = 0.0;
	dFinalTime = dTime[iNbTime - 1];
	iVolIndex = 0;
	iNumTime = 1;

	while (iNumTime < iNbTime)
	{
		if (iVolIndex == iNbSigTime || SigTime[iVolIndex] > dFinalTime - 1.0E-10)
		{			
			dNewTime = dFinalTime;

			if (iVolIndex == iNbSigTime) iVolIndex--;
		}
		else
		{
			dNewTime = SigTime[iVolIndex];
		}

		/* Nb Period and Dt */

		/* Check what is the min time in the period */
		dLocalMinTime = dTime[iNumTime] - dTime[iNumTime-1];

		iLocalIndexTime = iNumTime;

		while (dTime[iLocalIndexTime] < dNewTime - 1.0E-12 && dLocalMinTime > dMinTime)
		{
			iLocalIndexTime++;

			dLocalMinTime = min(dLocalMinTime, dTime[iLocalIndexTime] - dTime[iLocalIndexTime-1]);
		}

		dLocalMinTime = max(dLocalMinTime, dMinTime);

		iNbPeriod = max((int) ((dNewTime - dLastTime) / dLocalMinTime), 1);
		dDt = (dNewTime - dLastTime) / (iNbPeriod * 1.0);

		/* Constant calculations */

		/* Calculation of u^2 applied between t-1 and t*/
		du2 = Sig[iVolIndex] * Sig[iVolIndex];

		/* Calcuatio of  sig * rho / alpha */
		dRatio = Sig[iVolIndex] * dRho[iVolIndex] / dAlphaV[iVolIndex];
				
		/* Function constant */
		FuncVarft->a[0] = du2; 		
		FuncVt2Mean->dLambda = 2.0 * dLambdaV[iVolIndex];
		FuncVt2Mean->a[0] = dAlphaV[iVolIndex] * dAlphaV[iVolIndex] + 2.0 * dLvlV[iVolIndex];		
		FuncVtPsitMean->dLambda = dLambdaV[iVolIndex];
		FuncVtPsitMean->a[0] = du2;
		FuncVtPsitMean->a[1] = dLvlV[iVolIndex];				
		FuncPsit2Mean->a[0] = 2.0 * du2;				
		FuncVtftMean->dLambda = dLambdaV[iVolIndex];
		FuncVtftMean->a[0] = dAlphaV[iVolIndex] * dRho[iVolIndex] * Sig[iVolIndex];

		for (iPeriodIndex=0; iPeriodIndex<iNbPeriod; iPeriodIndex++)
		{
			if (iPeriodIndex == iNbPeriod - 1)
			{
				dDt = dNewTime - dLastTime;
			}

			dLastPeriodTime = dLastTime + dDt;

			/* Update value at t-1 on all func */
			FuncVarft->dXt1 = dVarft;
			FuncVt2Mean->dXt1 = dVt2Mean;
			FuncVtPsitMean->dXt1 = dVtPsitMean;
			FuncPsit2Mean->dXt1 = dPsit2Mean;
			FuncVtftMean->dXt1 = dVtftMean;

			/* Calculation at time t */

			/* Var(f(t,Tstar)-f(0,Tstar)) at t knowing t-1 */
			dNewVarft = LGMSVFuncValue1(FuncVarft, dDt, lambda, 0);

			/* E(Vt^2) at t knowing t-1 */
			dNewVt2Mean = LGMSVFuncValue1(FuncVt2Mean, dDt, lambda, 0);

			/* E(Psit) at t knowing t-1 */
			dNewPsitMean = dNewVarft;
			
			/* E(VtPsit) at t knowing t-1 */
			dNewVtPsitMean = LGMSVFuncValue1(FuncVtPsitMean, dDt, lambda, 0);

			/* E(Psit^2) at t knowing t-1 */
			dNewPsit2Mean = LGMSVFuncValue1(FuncPsit2Mean, dDt, lambda, 0);

			/* E(Vt*(f(t,Tstar)-f(0,Tstar))) at t knowing t-1 */
			dNewVtftMean = LGMSVFuncValue1(FuncVtftMean, dDt, lambda, 0);

			/* Fill the intermediary points */
			while (iNumTime < iNbTime && dTime[iNumTime] < dLastPeriodTime + 1.0E-10)
			{
				Sigt[iNumTime] = Sig[iVolIndex];
				dLambdaVt[iNumTime] = dLambdaV[iVolIndex];
				dLvlVt[iNumTime] = dLvlV[iVolIndex];
				dAlphaVt[iNumTime] = dAlphaV[iVolIndex];
				dRhot[iNumTime] = dRho[iVolIndex];

				dDtRatio = (dTime[iNumTime] - dLastTime) / dDt;

				/* Mean of Phi at t knowing t-1 */
				PsitMean[iNumTime] = (dNewVarft - dVarft) * dDtRatio + dVarft;

				dTemp = (dNewPsit2Mean - dPsit2Mean) * dDtRatio + dPsit2Mean - PsitMean[iNumTime] * PsitMean[iNumTime];

				if (dTemp > 0.0)
				{
					PsitStd[iNumTime] = sqrt(dTemp);
				}
				else
				{
					/* numerical problem... */
					PsitStd[iNumTime] = LGMSV_MINPSIVOL * sqrt(dDt);
				}

				/* Calculation of Vt Mean */	
				ZtMean[iNumTime] = 1;

				/* Calculation of Vt Std */

				dTemp = (dNewVt2Mean - dVt2Mean) * dDtRatio + dVt2Mean;

				if (dTemp > 1.0)
				{
					ZtStd[iNumTime] = sqrt(dTemp - 1.0);
				}
				else
				{
					ZtStd[iNumTime] = LGMSV_MINVTVOL * sqrt(dDt);
				}
				
				/* Mean of f(t,Tstar)-f(0,Tstar)- u(t)*dRho/dAlpha*(Vt-1) at t knowing t-1 */
				XtMean[iNumTime] = 0;		

				/* Std of Yt = f(t,Tstar)-f(0,Tstar) - u(t)*dRho/dAlpha*(Vt-1) at t knowing t-1 */
				/* We want that YtStd(t) <= YtStd(T) we use for the grid */
				/* => YtStd(t) should increase */

				dTempVtFtMean = (dNewVtftMean - dVtftMean) * dDtRatio + dVtftMean;

				dTemp = PsitMean[iNumTime] + pow(dRatio, 2) * (dTemp - 1.0)
					- 2.0 * dRatio * dTempVtFtMean;

				if (dTemp > 0.0)
				{
					XtStd[iNumTime] = max(XtStd[iNumTime-1], sqrt(dTemp));
				}
				else
				{
					if (iNumTime > 0)
					{
						XtStd[iNumTime] = XtStd[iNumTime-1];
					}
					else
					{
						XtStd[iNumTime] = LGMSV_MINORTHOVOL * sqrt(dDt);
					}
				}

				/* Calculation of the implied rho between Vt and ft times Std of Vt */
				ImpliedRhoTimesStdVt[iNumTime] = dTempVtFtMean / sqrt(PsitMean[iNumTime]);

				iNumTime++;
			}

			/* Update */
			dVarft = dNewVarft;
			dVt2Mean = dNewVt2Mean;
			dVtPsitMean = dNewVtPsitMean;
			dPsit2Mean = dNewPsit2Mean;
			dVtftMean = dNewVtftMean;

			dLastTime = dLastPeriodTime;
		}

		iVolIndex++;
		
	}

FREE_RETURN:

	if (lambda) free(lambda);
	if (FuncVarft) free(FuncVarft);
	if (FuncVt2Mean) free(FuncVt2Mean);
	if (FuncVtPsitMean) free(FuncVtPsitMean);
	if (FuncPsit2Mean) free(FuncPsit2Mean);
	if (FuncVtftMean) free(FuncVtftMean);

	return err;
}

static Err LGMSVPDEExpectations_UtPieceWise_OldVersion(
					/* Inputs */
					 int	iNbTime,			
					 double	*dTime,
					
					 /* For f(t,Tstar) diffusion and Xt*/
					 double dLambdaX,
					 double dSigma,
					 double dTstar,
					 
					 /* definition of u(t) piecewise constant */
					 double *SigTime,
					 double	*Sig,
					 int	iNbSigTime,

					 /* For Vt diffusion with TS on same dates as the volatility */
					 double *dLambdaV,
					 double	*dLvlV,
					 double *dAlphaV,
					 double	*dRho,
					 
					 /* Outputs */
					 double *Sigt,			/* Value of u(t) = sig(t) for all t in dtime */
					 double	*dLambdaVt,		/* Value of lameps(t) for all t in dtime */
					 double	*dLvlVt,		/* Value of mean-revert level */
					 double	*dAlphaVt,		/* Value of alphaeps(t) for all t in dtime */
					 double	*dRhot,			/* Value of rhoeps(t) for all t in dtime */
					 
					 double *PsitMean,
					 double *PsitStd,
					 
					 double *XtMean,
					 double *XtStd,

					 double *ZtMean,
					 double *ZtStd,
					 double *ImpliedRhoTimesStdVt)
{	
/* Declaration of locals variables */
int iNumTime, iNumSigTime;
double dDt;

double			*lambda			= NULL;

LGMSVSolFunc1	*FuncVarft		= NULL,
				*FuncVt2Mean	= NULL,
				*FuncVtPsitMean	= NULL,
				*FuncPsit2Mean	= NULL,
				*FuncVtftMean	= NULL;

double dVarft, dVt2Mean, dPsitMean;
double dVtPsitMean, dPsit2Mean, dVtftMean;
double dRatio;
double dTemp;

double du2;
Err	   err = NULL;

	/* All allocation errors */
	
	lambda = calloc(100, sizeof(double));
	FuncVarft = calloc(1, sizeof(LGMSVSolFunc1));
	FuncVt2Mean = calloc(1, sizeof(LGMSVSolFunc1));
	FuncVtPsitMean = calloc(1, sizeof(LGMSVSolFunc1));
	FuncPsit2Mean = calloc(1, sizeof(LGMSVSolFunc1));
	FuncVtftMean = calloc(1, sizeof(LGMSVSolFunc1));

	if (!lambda || !FuncVarft || !FuncVt2Mean || 
		!FuncVtPsitMean || !FuncPsit2Mean || !FuncVtftMean)
	{
		err = "Memory allocation faillure in LGMSVPDEExpectations_UtPieceWise";
		goto FREE_RETURN;
	}

	/* Initialisation of the LGMSVSolFunc */
		
	/* Definition of the function Var(f(t,Tstar)-f(0,Tstar)) */
	/* dE(f^2)/dt = ut^2 */
	FuncVarft->dLambda = 0.0;
	FuncVarft->iNbFunction = 1;
	
	FuncVarft->bIsft1[0] = 1;

	/* Definition of the function E(Vt^2) */
	/* dE(Vt^2)/dt = (dAlphaV*dAlphaV+2*dLvlt)-2*dLambdaVt*E(Vt^2) */	
	/*FuncVt2Mean->dLambda = 2*dLambdaV;*/
	FuncVt2Mean->iNbFunction = 1;

	FuncVt2Mean->bIsft1[0] = 1;
	/*FuncVt2Mean->a[0] = dAlphaV*dAlphaV+2*dLvlV;*/

	/* Definition of the function E(Psit) */
	/* dE(Psit)/dt = ut^2 = dE(f^2)/dt */
	/* => E(Psit) = Var(f(t,Tstar)-f(0,Tstar)) */
	

	/* Definition of the function E(VtPsit) */
	/* dE(VtPsit)/dt = dLvlV*E(Psit) + ut^2 *E(Vt^2) - dLambdaV*E(VtPsit) */
	/*FuncVtPsitMean->dLambda = dLambdaV;*/
	FuncVtPsitMean->iNbFunction = 2;

	FuncVtPsitMean->bIsft1[0] = 0;
	FuncVtPsitMean->pft[0] = FuncVt2Mean;
	
	FuncVtPsitMean->bIsft1[1] = 0;
	/*FuncVtPsitMean->a[1] = dLvlV;*/
	FuncVtPsitMean->pft[1] = FuncVarft;
	
	/* Definition of the function E(Psit^2) */
	/* dE(Psit^2)/dt = 2 * ut^2 * E(VtPsit) */
	FuncPsit2Mean->dLambda = 0;
	FuncPsit2Mean->iNbFunction = 1;

	FuncPsit2Mean->bIsft1[0] = 0;
	FuncPsit2Mean->pft[0] = FuncVtPsitMean;

	/* Definition of the function E(Vt*(f(t,Tstar)-f(0,Tstar))) */
	/* dE(Vt*ft)/dt = dAlphaV*dRho*ut-dLambdaV*E(Vt*ft) */
	/*FuncVtftMean->dLambda = dLambdaV;*/
	FuncVtftMean->iNbFunction = 1;

	FuncVtftMean->bIsft1[0] = 1;

	/* Filling the Sigt array 
	   piecewise interpolation from the Sig array on SigTime 
	   u(t) is right continuous									*/
	
	iNumSigTime = 0;

	for (iNumTime = 0; iNumTime<iNbTime; iNumTime++)
	{
		if ((dTime[iNumTime] > SigTime[iNumSigTime] + 1.0E-09)
			&& (iNumSigTime<iNbSigTime-1))
		{
			iNumSigTime++;
		}
		
		Sigt[iNumTime] = Sig[iNumSigTime];
		dLambdaVt[iNumTime] = dLambdaV[iNumSigTime];
		dLvlVt[iNumTime] = dLvlV[iNumSigTime];
		dAlphaVt[iNumTime] = dAlphaV[iNumSigTime];
		dRhot[iNumTime] = dRho[iNumSigTime];
	}
	
	/* Init Value at t=0 */
	dVarft = 0.0;
	dVt2Mean = 1.0;
	dPsitMean = 0.0;
	dVtPsitMean = 0.0;
	dPsit2Mean = 0.0;
	dVtftMean = 0.0;
	
	PsitMean[0] = 0.0;
	PsitStd[0] = 0.0;

	ZtMean[0] = 1.0;
	ZtStd[0] = 0.0;

	XtMean[0] = 0.0;
	XtStd[0] = 0.0;

	ImpliedRhoTimesStdVt[0] = 0.0;
	
	/* Calculation of all the mean and variances */
	for (iNumTime = 1;iNumTime<iNbTime;iNumTime++)
	{
		/* Delta time between t and t-1 */
		dDt = dTime[iNumTime] - dTime[iNumTime-1];
		
		/* Calculation of u^2 applied between t-1 and t*/
		du2 = Sigt[iNumTime] * Sigt[iNumTime];

		/* Calcuatio of  sig * rho / alpha */
		dRatio = Sigt[iNumTime] * dRhot[iNumTime] / dAlphaVt[iNumTime];

		/* Update value at t-1 on all func */
		FuncVarft->dXt1 = dVarft;
		FuncVarft->a[0] = du2;
 
		FuncVt2Mean->dXt1 = dVt2Mean;
		FuncVt2Mean->dLambda = 2.0 * dLambdaVt[iNumTime];
		FuncVt2Mean->a[0] = dAlphaVt[iNumTime] * dAlphaVt[iNumTime] + 2.0 * dLvlVt[iNumTime];

		FuncVtPsitMean->dXt1 = dVtPsitMean;
		FuncVtPsitMean->dLambda = dLambdaVt[iNumTime];
		FuncVtPsitMean->a[0] = du2;
		FuncVtPsitMean->a[1] = dLvlVt[iNumTime];
		
		FuncPsit2Mean->dXt1 = dPsit2Mean;
		FuncPsit2Mean->a[0] = 2.0 * du2;
		
		FuncVtftMean->dXt1 = dVtftMean;
		FuncVtftMean->dLambda = dLambdaVt[iNumTime];
		FuncVtftMean->a[0] = dAlphaVt[iNumTime] * dRhot[iNumTime] * Sigt[iNumTime];

		/* Calculation at time t */

		/* Var(f(t,Tstar)-f(0,Tstar)) at t knowing t-1 */
		dVarft = LGMSVFuncValue1(FuncVarft, dDt, lambda, 0);

		/* E(Vt^2) at t knowing t-1 */
		dVt2Mean = LGMSVFuncValue1(FuncVt2Mean, dDt, lambda, 0);

		/* E(Psit) at t knowing t-1 */
		dPsitMean = dVarft;
		
		/* E(VtPsit) at t knowing t-1 */
		dVtPsitMean = LGMSVFuncValue1(FuncVtPsitMean, dDt, lambda, 0);

		/* E(Psit^2) at t knowing t-1 */
		dPsit2Mean = LGMSVFuncValue1(FuncPsit2Mean, dDt, lambda, 0);

		/* E(Vt*(f(t,Tstar)-f(0,Tstar))) at t knowing t-1 */
		dVtftMean = LGMSVFuncValue1(FuncVtftMean, dDt, lambda, 0);

		/* Mean of Phi at t knowing t-1 */
		PsitMean[iNumTime] = dPsitMean;

		dTemp = dPsit2Mean - dPsitMean * dPsitMean;

		if (dTemp > 0.0)
		{
			PsitStd[iNumTime] = sqrt(dTemp);
		}
		else
		{
			/* numerical problem... */
			PsitStd[iNumTime] = LGMSV_MINPSIVOL * sqrt(dDt);
		}
		
		/* Calculation of Vt Mean */	
		ZtMean[iNumTime] = 1;

		/* Calculation of Vt Std */
		if (dVt2Mean > 1.0)
		{
			ZtStd[iNumTime] = sqrt(dVt2Mean - 1.0);
		}
		else
		{
			PsitStd[iNumTime] = LGMSV_MINVTVOL * sqrt(dDt);
		}
		
		/* Mean of f(t,Tstar)-f(0,Tstar)- u(t)*dRho/dAlpha*(Vt-1) at t knowing t-1 */
		XtMean[iNumTime] = 0;		

		/* Std of Yt = f(t,Tstar)-f(0,Tstar) - u(t)*dRho/dAlpha*(Vt-1) at t knowing t-1 */
		/* We want that YtStd(t) <= YtStd(T) we use for the grid */
		/* => YtStd(t) should increase */
		dTemp = dVarft + pow(dRatio, 2) * (dVt2Mean - 1.0) - 2.0 * dRatio * dVtftMean;

		if (dTemp > 0.0)
		{
			XtStd[iNumTime] = max(XtStd[iNumTime-1], sqrt(dTemp));
		}
		else
		{
			if (iNumTime > 0)
			{
				XtStd[iNumTime] = XtStd[iNumTime-1];
			}
			else
			{
				XtStd[iNumTime] = LGMSV_MINORTHOVOL * sqrt(dDt);
			}
		}

		/* Calculation of the implied rho between Vt and ft times Std of Vt */
		ImpliedRhoTimesStdVt[iNumTime] = dVtftMean / sqrt(dPsitMean);
	}

FREE_RETURN:

	if (lambda) free(lambda);
	if (FuncVarft) free(FuncVarft);
	if (FuncVt2Mean) free(FuncVt2Mean);
	if (FuncVtPsitMean) free(FuncVtPsitMean);
	if (FuncPsit2Mean) free(FuncPsit2Mean);
	if (FuncVtftMean) free(FuncVtftMean);

	return err;
}


/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEGridX_UtPieceWise


   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEGridX_UtPieceWise(/* Inputs */
								int iNb,
								double dStd,
								double iNbSigmaMaxLeft,
								double iNbSigmaLeft,
								double iNbSigmaRight,
								double iNbSigmaMaxRight,
								double dPercent,
								int IsEqSpacedExtremePoint,
									 
							 /* Outputs */
								double *Grid,
								int	*pIndex0)
{
/* Declaration of locals variables */
int iNumStep, iNb1,iNbLeft,iNbRight;
double dStep,dStepLeft; /*,dStepRight;*/
double dMin;

double dMax,dCumulStep;

	if (iNbSigmaMaxLeft<iNbSigmaLeft)
		iNbSigmaMaxLeft = iNbSigmaLeft;
	if (iNbSigmaMaxRight<iNbSigmaRight)
		iNbSigmaMaxRight = iNbSigmaRight;

	if ((iNbSigmaMaxLeft == iNbSigmaLeft)&&(iNbSigmaMaxRight==iNbSigmaRight))
		dPercent = 0;

	/* Calculation of the step in X grid */
	if (iNb<2)
	{
		Grid[0]=0;
		*pIndex0 = 0;
	}
	else
	{
		iNb1 = iNb - (int) (dPercent*iNb);
		
		dStep = ((iNbSigmaLeft+iNbSigmaRight)*dStd)/(iNb1-1);

		/* Calculation of the min value of X */
		dMin =  - iNbSigmaLeft*dStd;

		/* Shift the grid so that X=0.0 is in the grid */
		if ((dMin <=0)&&(iNb>=2))
		{
			*pIndex0 =(int)((0-dMin)/dStep+1e-08);
			dMin = 0-(*pIndex0)*dStep;
		}
		else
		{
			dMin=0;
			*pIndex0 = 0;
		}
		
		if (iNbSigmaMaxRight == iNbSigmaRight)
		{
			iNbLeft = (int) (dPercent*iNb/1.0);
		}
		else
		{
			if (iNbSigmaMaxLeft == iNbSigmaLeft)
			{
				iNbLeft = 0;
			}
			else
			{
				iNbLeft = (int) (dPercent*iNb/2.0);
			}
		}

		if (iNbLeft == 1)
		{
			if (dMin>-iNbSigmaMaxLeft*dStd)
			{
				Grid[0] = -iNbSigmaMaxLeft*dStd;
			}
			else
			{
				iNbLeft = 0;
			}
		}
		else if (iNbLeft>0)
		{
		
			dStepLeft = (dMin-(-iNbSigmaMaxLeft)*dStd)/(iNbLeft);
			for (iNumStep=0; iNumStep<iNbLeft; iNumStep++ )
				Grid[iNumStep]=(-iNbSigmaMaxLeft)*dStd+iNumStep*dStepLeft;	
		}

		*pIndex0 = *pIndex0 + iNbLeft;

		/* Filling the center grid in X */
		/* X = 0 is in the grid  */
		Grid[iNbLeft] = dMin;
		for (iNumStep=iNbLeft+1; iNumStep<iNb1+iNbLeft; iNumStep++ )
			Grid[iNumStep]=Grid[iNumStep-1]+dStep;

		iNbRight = iNb -(iNbLeft +iNb1);
		
		if (iNbRight == 1)
		{
			if (Grid[iNb1+iNbLeft-1]<iNbSigmaMaxRight*dStd)
			{
				Grid[iNb1+iNbLeft] = iNbSigmaMaxRight*dStd;
			}
			else
			{
				Grid[iNb1+iNbLeft] = Grid[iNb1+iNbLeft-1]+dStep;
			}
		}
		else if (iNbRight>0)
		{
			/* iNbRight is > than 1 */
			/*dStepRight = (iNbSigmaMaxRight*dStd-Grid[iNb1+iNbLeft-1])/(iNbRight-1.0);
			for (iNumStep=iNb1+iNbLeft; iNumStep<iNb; iNumStep++ )
				Grid[iNumStep]=Grid[iNb1+iNbLeft-1]+(iNumStep-(iNb1+iNbLeft)+1)*dStepRight;	
			*/

			dMax = iNbSigmaMaxRight*dStd;

			dCumulStep = ((dMax - Grid[iNb1+iNbLeft-1])-iNbRight*dStep)*2.0/(iNbRight*(iNbRight+1.0));
			if ((dCumulStep <0) || (IsEqSpacedExtremePoint == LGMSV_TRUE) )
			{
				dStep = (dMax - Grid[iNb1+iNbLeft-1])/(iNbRight-1.0);
				for (iNumStep=iNb1+iNbLeft; iNumStep<iNb; iNumStep++ )
					Grid[iNumStep]=Grid[iNumStep-1]+dStep;
			}
			else
			{
				for (iNumStep=iNb1+iNbLeft; iNumStep<iNb; iNumStep++ )
					Grid[iNumStep]=Grid[iNumStep-1]+dStep+(iNumStep-(iNb1+iNbLeft)+1)*dCumulStep;
			}


		}
	}
	/* Force to 0 */
	Grid[*pIndex0] = 0.0;

}


/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEGridZ_UtPieceWise


   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEGridZ_UtPieceWise(/* Inputs */
								int iNb,
								double dMean,
								double dStd,
								double iNbSigmaLeft,
								double iNbSigmaRight,
								double iNbSigmaMaxRight,
								double dPercent,
								int IsEqSpacedExtremePoint,
									 
							 /* Outputs */
								double *Grid,
								int	*pIndex0)
{
/* Declaration of locals variables */
int iNumStep, iNb1, iNb2;
double dStep;
double dMin, dMin1;
double dMax, dCumulStep;

	if (iNbSigmaMaxRight == iNbSigmaRight)
		dPercent = 0;

	/* Calculation of the step in X grid */
	if (iNb<2)
	{
		Grid[0]=dMean;
		*pIndex0 = 0;
	}
	else
	{
		iNb1 = iNb - (int) (dPercent*iNb);
		/* Force to odd Nb */
		iNb1 = ((int)(iNb1/2.0))*2+1;

		if (iNbSigmaLeft*dStd > dMean)
		{
			dMin = dMean*exp(-iNbSigmaLeft*dStd);
		}
		else
		{
			dMin = dMean-iNbSigmaLeft*dStd;
		}

		dStep = (iNbSigmaRight*dStd+dMean-dMin) / (iNb1-1.0);
		dMin1 = dMin;

		if ((dMin1 <=dMean)&&(iNb>=2))
		{
			*pIndex0 =(int)((dMean-dMin1)/dStep+1e-08);
			dMin1 = dMean-(*pIndex0)*dStep;
		}
		else
		{
			dMin1=dMean;
			*pIndex0 = 0;
		}

		if ((fabs(dMin1-dMin)>0) && (iNbSigmaLeft*dStd > dMean))
		{
			*pIndex0 = *pIndex0 + 1;
			Grid[0] = dMin;
			for (iNumStep=1; iNumStep<iNb1; iNumStep++ )
				Grid[iNumStep]=dMin1+(iNumStep-1)*dStep;
		}
		else
		{
			for (iNumStep=0; iNumStep<iNb1; iNumStep++ )
				Grid[iNumStep]=dMin1+(iNumStep)*dStep;
		}
		
		
		/* Filling the grid in Z */

		iNb2 = iNb - iNb1;
		dMax = dMean+iNbSigmaMaxRight*dStd;
		if (iNb2>0)
		{
			dCumulStep = ((dMax - Grid[iNb1-1])-iNb2*dStep)*2.0/(iNb2*(iNb2+1.0));
			if ((dCumulStep <0) || (IsEqSpacedExtremePoint == LGMSV_TRUE) )
			{
				dStep = (dMax - Grid[iNb1-1])/(iNb2-1.0);
				for (iNumStep=iNb1; iNumStep<iNb; iNumStep++ )
					Grid[iNumStep]=Grid[iNumStep-1]+dStep;
			}
			else
			{
				for (iNumStep=iNb1; iNumStep<iNb; iNumStep++ )
					Grid[iNumStep]=Grid[iNumStep-1]+dStep+(iNumStep-iNb1+1)*dCumulStep;
			}
		}


	}

}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEGridZ_UtPieceWise


   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEGridZ_UtPieceWise_New(/* Inputs */
								int iNb,
								double dMean,
								double dStd,
								double iNbSigmaLeft,
								double iNbSigmaRight,
								double iNbSigmaMaxRight,
								double dPercent,
								int IsEqSpacedExtremePoint,
									 
							 /* Outputs */
								double *Grid,
								int	*pIndex0)
{
/* Declaration of locals variables */
int iNumStep, iNb1, iNb2;
double dStep;
double dMin;
double dMax, dCumulStep;
int iNb1_r, iNb1_g;
double dNbRightSym;

	if (iNbSigmaMaxRight == iNbSigmaRight)
		dPercent = 0;

	/* Calculation of the step in X grid */
	if (iNb<2)
	{
		Grid[0]=dMean;
		*pIndex0 = 0;
	}
	else
	{
		if (iNbSigmaLeft*dStd > dMean)
		{
			dMin = dMean*exp(-iNbSigmaLeft*dStd);
		}
		else
		{
			dMin = dMean-iNbSigmaLeft*dStd;
		}
		
		dNbRightSym = (dMean-dMin)/dStd;
		if (dNbRightSym>iNbSigmaRight)
				dNbRightSym = iNbSigmaRight;

		
		iNb1 = iNb - (int) (dPercent*iNb);
		/* Force to odd Nb */	
		iNb1 = ((int)(iNb1/2.0))*2+1;
		
		iNb1_g = (int) (iNb1*(dNbRightSym+iNbSigmaLeft)/(iNbSigmaLeft+iNbSigmaRight));
		/* Force to odd Nb */
		iNb1_g = ((int)(iNb1_g/2.0))*2+1;

		dStep = (dNbRightSym*dStd+dMean-dMin) / (iNb1_g-1.0);
		*pIndex0 = (int)(iNb1_g/2.0);
		
		for (iNumStep=0; iNumStep<iNb1_g; iNumStep++ )
			Grid[iNumStep]=dMin+(iNumStep)*dStep;

		iNb1_r = iNb1-iNb1_g;
		if (iNb1_r >= 1.0 )
		{
			dStep = (iNbSigmaRight*dStd+dMean- Grid[iNb1_g-1])/ (iNb1_r);
		}

		for (iNumStep=iNb1_g; iNumStep<iNb1; iNumStep++ )
				Grid[iNumStep]=Grid[iNumStep-1]+dStep;
		
		/* Filling the grid in Z */

		iNb2 = iNb - iNb1;
		dMax = dMean+iNbSigmaMaxRight*dStd;
		if (iNb2>0)
		{
			dCumulStep = ((dMax - Grid[iNb1-1])-iNb2*dStep)*2.0/(iNb2*(iNb2+1.0));
			if ((dCumulStep <0) || (IsEqSpacedExtremePoint == LGMSV_TRUE) )
			{
				dStep = (dMax - Grid[iNb1-1])/(iNb2-1.0);
				for (iNumStep=iNb1; iNumStep<iNb; iNumStep++ )
					Grid[iNumStep]=Grid[iNumStep-1]+dStep;
			}
			else
			{
				for (iNumStep=iNb1; iNumStep<iNb; iNumStep++ )
					Grid[iNumStep]=Grid[iNumStep-1]+dStep+(iNumStep-iNb1+1)*dCumulStep;
			}
		}


	}

}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEMakePsiGrid														  

	Compute the dynamic grid on Psi

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEMakePsiGrid_UtPieceWise(/* Inputs */
								int iNbPsi,
								double dPsiMean,
								double dPsiStd,
								double iNbSigmaPsiGrid,

								/* Parameter */
								int IsExtremePoints,
								double iNbSigmaExtremePoints,
								int IsPsiBoundExp,
								
								/* Outputs */
								double *GridPsi,
								int	   *pIndexPsiMean)
{
/* Declaration of locals variables */
double  dLogPsiStd, dLogPsiMean;
double dPsiStep, dPsiMin, dNewPsiMin;
int iNumStep;

double PsiFloor = 0.0;	
double dPsiMinExtreme, dPsiMaxExtreme;
		
	if (iNbPsi == 1)
	{
		GridPsi[0] = dPsiMean;
		*pIndexPsiMean = 0;
	}
	else
	{
		switch (IsPsiBoundExp)
		{
		/*--------------------
			PsiBoundExp Grid
		 ----------------------*/
		case LGMSV_TRUE:
			{
				if (dPsiMean ==0)
				{
					/* Filling the grid in Psi with constant value */
					for (iNumStep=0; iNumStep<iNbPsi; iNumStep++ )
						GridPsi[iNumStep] = dPsiMean;
				}
				else		
				{
					/* Calculation of mean and Std of log Psi*/
					dLogPsiStd = sqrt(log(1+dPsiStd*dPsiStd/dPsiMean/dPsiMean));
					dLogPsiMean = log(dPsiMean)-dLogPsiStd*dLogPsiStd/2;
					
					switch (IsExtremePoints)
					{
					case LGMSV_TRUE:
						{
							if (iNbPsi == 3)
							{
								GridPsi[0]=dPsiMean*exp(-iNbSigmaExtremePoints*dLogPsiStd);
								GridPsi[1]=dPsiMean;
								GridPsi[2]=dPsiMean*exp(iNbSigmaExtremePoints*dLogPsiStd);
								*pIndexPsiMean = 1;
							
							}
							else
							{
								/* Calculation of the step in Psi grid */
								dPsiStep = dPsiMean*(exp(iNbSigmaPsiGrid*dLogPsiStd)-
												exp(-iNbSigmaPsiGrid*dLogPsiStd))/(iNbPsi-3);

								/* Calculation of the min value of Psi */
								dPsiMin = dPsiMean*exp(-iNbSigmaPsiGrid*dLogPsiStd);

								/* Shift the grid so that dPsiMean is in the grid */
								*pIndexPsiMean =(int) floor((dPsiMean-dPsiMin)/dPsiStep);
								dNewPsiMin = dPsiMean-(*pIndexPsiMean)*dPsiStep;

								*pIndexPsiMean += 1;
							
								/* Filling the grid in Psi */
								for (iNumStep=1; iNumStep<iNbPsi-1; iNumStep++ )
									GridPsi[iNumStep]=dNewPsiMin+(iNumStep-1)*dPsiStep;
								GridPsi[0] = dPsiMean*exp(-iNbSigmaExtremePoints*dLogPsiStd);
								GridPsi[iNbPsi-1] = dPsiMean*exp(iNbSigmaExtremePoints*dLogPsiStd);
							}
							break;
						}
					
					default :
						{
							/* Calculation of the step in Psi grid */
							dPsiStep = dPsiMean*(exp(iNbSigmaPsiGrid*dLogPsiStd)-
												exp(-iNbSigmaPsiGrid*dLogPsiStd))/(iNbPsi-1);

							/* Calculation of the min value of LogPsi */
							dPsiMin = dPsiMean*exp(-iNbSigmaPsiGrid*dLogPsiStd);

							/* Shift the grid so that log(dPsiMean) is in the grid */
							*pIndexPsiMean =(int) floor((dPsiMean-dPsiMin)/dPsiStep);
							dNewPsiMin = dPsiMean-(*pIndexPsiMean)*dPsiStep;
							
							/* Filling the grid in Psi */
							for (iNumStep=0; iNumStep<iNbPsi; iNumStep++ )
								GridPsi[iNumStep]=dNewPsiMin+(iNumStep-1)*dPsiStep;

							break;

						}
					}	
				}
				break;
			}
		/*-------------
			Unif Grid
		 --------------*/
		default :
			{
				switch (IsExtremePoints)
				{
				case LGMSV_TRUE:
					{
						if (iNbPsi == 3)
						{
							/* Psi Should be positive */
							GridPsi[0]=max(dPsiMean-iNbSigmaExtremePoints*dPsiStd,PsiFloor);
							
							GridPsi[1]=dPsiMean;
							GridPsi[2]=dPsiMean+iNbSigmaExtremePoints*dPsiStd;
							*pIndexPsiMean = 1;
						}
						else
						{
							/* Calculation of the Psi Min Extreme point */
							dPsiMinExtreme = max(dPsiMean-iNbSigmaExtremePoints*dPsiStd,PsiFloor);
							dPsiMaxExtreme = dPsiMean+iNbSigmaExtremePoints*dPsiStd;
							
							/* Calculation of the step in Psi grid */
							dPsiStep = (dPsiMean+iNbSigmaPsiGrid*dPsiStd-max(dPsiMean-iNbSigmaPsiGrid*dPsiStd,PsiFloor))
								/(iNbPsi-3);

							/* Calculation of the min value of Psi */
							dPsiMin = max(dPsiMean - iNbSigmaPsiGrid*dPsiStd,PsiFloor);

							/* Shift the grid so that dPsiMean is in the grid */
							*pIndexPsiMean =(int) floor((dPsiMean-dPsiMin)/dPsiStep);
							dNewPsiMin = dPsiMean-(*pIndexPsiMean)*dPsiStep;

							
							/* Filling the grid in Psi */
							if (dNewPsiMin != dPsiMinExtreme)
							{
								*pIndexPsiMean += 1;
								
								for (iNumStep=1; iNumStep<iNbPsi-1; iNumStep++ )
									GridPsi[iNumStep]=dNewPsiMin+(iNumStep-1)*dPsiStep;
							
								GridPsi[0] = dPsiMinExtreme;
								GridPsi[iNbPsi-1] = dPsiMaxExtreme;
							}
							else
							{	
								for (iNumStep=0; iNumStep<iNbPsi-1; iNumStep++ )
									GridPsi[iNumStep]=dNewPsiMin+iNumStep*dPsiStep;
							
								GridPsi[iNbPsi-1] = dPsiMaxExtreme;
							}

						}
						break;
					}
				default:
					{
						/* Calculation of the step in Psi grid */
						dPsiStep = (dPsiMean+iNbSigmaPsiGrid*dPsiStd-max(dPsiMean-iNbSigmaPsiGrid*dPsiStd,PsiFloor))
							/(iNbPsi-1);

						/* Calculation of the min value of Psi */
						dPsiMin = max(dPsiMean - iNbSigmaPsiGrid*dPsiStd,PsiFloor);

						/* Shift the grid so that dPsiMean is in the grid */
						*pIndexPsiMean =(int) floor((dPsiMean-dPsiMin)/dPsiStep);
						dPsiMin = dPsiMean-(*pIndexPsiMean)*dPsiStep;
							
						/* Filling the grid in Psi */
						for (iNumStep=0; iNumStep<iNbPsi; iNumStep++ )
							GridPsi[iNumStep]=dPsiMin+iNumStep*dPsiStep;
						break;
					}
				}
				break;
			}
		}	
	}
}


/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDECalculationPayoffMatrixfromPsit												  

	

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDECalculationValueMatrixfromPsit_UtPieceWise(
											 /* Inputs */
											 int	iIndexMinZ,
											 int	iIndexMaxZ,
											 int	iIndexMinX,
											 int	iIndexMaxX,
											 
											 int	iNbPsi_t,
											 int	iNbPsi_t1,
											 int	iProductNb,
											 
											 double	*GridftTstar, 
											 double	*GridZCenter,
											 
											 double	*GridPsi_t,
											 double	*GridPsi_t1,


											 /* Martingale Mesure QTstar information */
											 double	dTstar,		/* In time */

											 /* 4 dimensions : Zt, Xt , Phit, Product	*/
											 double ****ValueTensor,

											 double dSig_t1,
											 double dLambdaX,
											 double dAlphaV,
											 double	dRho,

											 double t, /*time(t) */
											 double t1, /*time(t+1) */

											 int TypeInterpolationPhi,
											 int IsExtrapolationFlat,

											 double **InterpQuadraCoef,

											 /* Output */
											 /* 4 dimensions : Zt, Xt , Phit, Product	*/
											 double ****ValueTensorNew)
{
/* Declaration of locals variables */
double dDt, dVt, dPsi_t1, dPsi_t1_sq;
int iNumZ, iNumX,iNumPsi, iNumProduct;
int iIndexPsi;

double dAlphaRho_Sig_t1, dSig_t12Dt;

int iIndex1, iIndex2, iIndex3;
double	***ValueTensorNewZ, **ValueTensorNewZX;
double	***ValueTensorZ, **ValueTensorZX, *ValueTensorZX1, *ValueTensorZX2, *ValueTensorZX3;

	/* variation of time between t and t+1 */
	dDt = t1-t;

	/* Calculation of constant */
	dAlphaRho_Sig_t1 = dAlphaV * dRho; /*/dSig_t1; */
	dSig_t12Dt = dSig_t1 * dSig_t1 * dDt;
	
	for (iNumZ=iIndexMinZ;iNumZ<=iIndexMaxZ;iNumZ++)
	{
		ValueTensorNewZ = ValueTensorNew[iNumZ];
		ValueTensorZ = ValueTensor[iNumZ];

		/* Calculation of the Vt from  the couple (ft,Zt) */
		/* dVt = max((GridZCenter[iNumZ]+dZtMean)+dAlpha*dRho/dSigt*GridftTstar[iNumX],dVFloor); */
		dVt = dSig_t12Dt * GridZCenter[iNumZ];

		for (iNumX=iIndexMinX;iNumX<=iIndexMaxX;iNumX++)
		{
			ValueTensorNewZX = ValueTensorNewZ[iNumX];
			ValueTensorZX = ValueTensorZ[iNumX];

			for (iNumPsi = 0; iNumPsi<iNbPsi_t; iNumPsi++)
			{
				iIndexPsi = iNumPsi;
		
				/* Calculation of the Psi(t+1)    */
				/*  = dPsit+dSigt2*dDt*dVt; */
				/* Because u(t) is continuous to the right, from t->t+dt */
				/* the applied vol is u(t+dt), the Vt used is also the Vt with u(t+dt) */
				dPsi_t1 = GridPsi_t[iNumPsi] + dVt;
				dPsi_t1_sq = dPsi_t1 * dPsi_t1;

				/* Find the Index such that GridPsi_t1[index]<dPsi_t1 and GridPsi_t1[index+1]>= dPsi_t1 */
				while ((iIndexPsi < iNbPsi_t1) && (GridPsi_t1[iIndexPsi] < dPsi_t1))
				{
					iIndexPsi++;
				}
				iIndexPsi--;

				if (IsExtrapolationFlat == LGMSV_TRUE)
				{
					/* Interpolate value */
					if ((iIndexPsi > 1) && (iIndexPsi < iNbPsi_t1-3) && (TypeInterpolationPhi == LGMSV_MIXQUADRALINEAR))
					{
						/* Quadratic Interpolation  */
						if (dPsi_t1 - GridPsi_t1[iIndexPsi] > 0.5 * (GridPsi_t1[iIndexPsi+1] - GridPsi_t1[iIndexPsi]))
						{						
							iIndex1 = iIndexPsi;
						}					
						else
						{
							iIndex1 = iIndexPsi - 1;
						}					

						iIndex2 = iIndex1 + 1;
						iIndex3 = iIndex2 + 1;
						
						ValueTensorZX1 = ValueTensorZX[iIndex1];
						ValueTensorZX2 = ValueTensorZX[iIndex2];
						ValueTensorZX3 = ValueTensorZX[iIndex3];

						for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
						{
							ValueTensorNewZX[iNumPsi][iNumProduct] = 
								ValueTensorZX1[iNumProduct]*InterpQuadraCoef[iIndex1][0]*
								(dPsi_t1_sq-InterpQuadraCoef[iIndex1][1]*dPsi_t1+InterpQuadraCoef[iIndex1][2])+
								ValueTensorZX2[iNumProduct]*InterpQuadraCoef[iIndex1][3]*
								(dPsi_t1_sq-InterpQuadraCoef[iIndex1][4]*dPsi_t1+InterpQuadraCoef[iIndex1][5])+
								ValueTensorZX3[iNumProduct]*InterpQuadraCoef[iIndex1][6]*
								(dPsi_t1_sq-InterpQuadraCoef[iIndex1][7]*dPsi_t1+InterpQuadraCoef[iIndex1][8]);					
						}					
					}
					else
					if (iIndexPsi == -1)
					{
						/* Extrapolation flat */						
						for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
						{
							ValueTensorNewZX[iNumPsi][iNumProduct] = 
									ValueTensorZX[0][iNumProduct];
						}
					}	
					else
					if (iIndexPsi == iNbPsi_t1-1)
					{
						/* Extrapolation flat */						
						for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
						{
							ValueTensorNewZX[iNumPsi][iNumProduct] = 
								ValueTensorZX[iNbPsi_t1-1][iNumProduct];
						}
					}
					else
					{
						/* Linear Interpolation */						
						for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
						{
							ValueTensorNewZX[iNumPsi][iNumProduct] = (dPsi_t1-GridPsi_t1[iIndexPsi])/(GridPsi_t1[iIndexPsi+1]-GridPsi_t1[iIndexPsi])*
								(ValueTensorZX[iIndexPsi+1][iNumProduct]-ValueTensorZX[iIndexPsi][iNumProduct])
								+ ValueTensorZX[iIndexPsi][iNumProduct];
						}
					}
				}
				else
				{
					/* Interpolate value */
					if ((iIndexPsi>1)&&(iIndexPsi<iNbPsi_t1-3)&&(TypeInterpolationPhi==LGMSV_MIXQUADRALINEAR))
					{
						/* Quadratic Interpolation  */
						if (dPsi_t1-GridPsi_t1[iIndexPsi]>0.5*(GridPsi_t1[iIndexPsi+1]-GridPsi_t1[iIndexPsi]))
						{
							iIndex1 = iIndexPsi;
							iIndex2 = iIndexPsi+1;
							iIndex3 = iIndexPsi+2;
						}
						else
						{
							iIndex1 = iIndexPsi-1;
							iIndex2 = iIndexPsi;
							iIndex3 = iIndexPsi+1;						
						}
						
						for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
						{
							ValueTensorNew[iNumZ][iNumX][iNumPsi][iNumProduct] = 
								ValueTensor[iNumZ][iNumX][iIndex1][iNumProduct] * InterpQuadraCoef[iIndex1][0] *
								(dPsi_t1_sq - InterpQuadraCoef[iIndex1][1] * dPsi_t1 + InterpQuadraCoef[iIndex1][2]) +
								ValueTensor[iNumZ][iNumX][iIndex2][iNumProduct] * InterpQuadraCoef[iIndex1][3] *
								(dPsi_t1_sq - InterpQuadraCoef[iIndex1][4] * dPsi_t1 + InterpQuadraCoef[iIndex1][5]) +
								ValueTensor[iNumZ][iNumX][iIndex3][iNumProduct] * InterpQuadraCoef[iIndex1][6] *
								(dPsi_t1_sq - InterpQuadraCoef[iIndex1][7] * dPsi_t1 + InterpQuadraCoef[iIndex1][8]);
						}
					}
					else
					{
						/* Interpolate value */
						if (iIndexPsi == -1)
						{
							iIndexPsi = 0;		
						}	
						else
						if (iIndexPsi == iNbPsi_t1-1)
						{
							/* Extrapolation flat */
							iIndexPsi= max(iNbPsi_t1 - 2, 0);
						}

						if (iNbPsi_t1 == 1)
						{
							/* Extrapolation flat */						
							for (iNumProduct=0; iNumProduct<iProductNb; iNumProduct++)
							{
								ValueTensorNew[iNumZ][iNumX][iNumPsi][iNumProduct] = 
										ValueTensor[0][iNumZ][iNumX][iNumProduct];
							}
						}
						else
						{
							/* Linear Interpolation */
							for (iNumProduct=0; iNumProduct<iProductNb; iNumProduct++)
							{
								ValueTensorNew[iNumZ][iNumX][iNumPsi][iNumProduct] = (dPsi_t1-GridPsi_t1[iIndexPsi])/(GridPsi_t1[iIndexPsi+1]-GridPsi_t1[iIndexPsi])*
									(ValueTensor[iNumZ][iNumX][iIndexPsi+1][iNumProduct]-ValueTensor[iNumZ][iNumX][iIndexPsi][iNumProduct])
									+ ValueTensor[iNumZ][iNumX][iIndexPsi][iNumProduct];
							}
						}
					}
				}
			}
		}
	}
}


/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDECalculationRDriftVar_UtPieceWise												  

	

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDECalculationRDriftVar_UtPieceWise(/* Inputs */
										int iIndexMinX,
										int	iIndexMaxX,
										int	iIndexMinZ,
										int	iIndexMaxZ,
										
										double	*GridX,
										double	*GridZ,
									
										double dSig_t1,

										double dLambdaV,
										double dLvlV,
										double dAlphaV,
										double dRho,

										double t, /*time(t) */
										double t1, /*time(t+1) */
																				
										/* Outputs */
										double **DXDrift,
										double **DZDrift,
										double **DXVar,
										double **DZVar,
										double **r)
{
/* Declaration of locals variables */
double dDt;
double dVt, dVt_m1;
int iNumX, iNumZ;
double dConstXDrift, dConstXVar;
double dConstZDrift, dConstZVar;
double dRatio;

	/* variation of time between t and t+1 */
	dDt = t1-t;

	/* Constants */
	dConstXDrift = dLambdaV*dRho*dSig_t1/dAlphaV*dDt;
	dConstXVar = dSig_t1*dSig_t1*(1-dRho*dRho)*dDt;
	dConstZDrift = -dLambdaV*dDt;
	dConstZVar = dAlphaV*dAlphaV*dDt;
	
	if (fabs(dLambdaV) < 1.0E-08)
	{
		dRatio = 0.0;
	}
	else
	{
		dRatio = dLvlV / dLambdaV;
	}

	/* Calculaton drift and variance */
	for (iNumZ=iIndexMinZ;iNumZ<=iIndexMaxZ;iNumZ++)
	{
		/* Calculation of Vt */
		dVt = GridZ[iNumZ];	
		/* Calculation of Vt-1 */
		dVt_m1 = dVt - dRatio;

		r[iNumZ][iIndexMinX] = 0;
						
		/* dt*Drifts under QTstar and dt*Variances under QTstar */
			
		/* E[Xt+1 - Xt / Ft] and V[Xt+1 - Xt / Ft]*/
		DXDrift[iNumZ][iIndexMinX] = dConstXDrift*dVt_m1;
		DXVar[iNumZ][iIndexMinX] = dConstXVar*dVt;
		
		/* E[Zt+1 - Zt / Ft] and V[Zt+1 - Zt / Ft]*/
		DZDrift[iNumZ][iIndexMinX] = dConstZDrift*dVt_m1;
		DZVar[iNumZ][iIndexMinX] = dConstZVar*dVt;
				
		for (iNumX=iIndexMinX+1;iNumX<=iIndexMaxX;iNumX++)
		{	
			/* Discount rate * dt */
			r[iNumZ][iNumX] = 0;
						
			/* dt*Drifts under QTstar and dt*Variances under QTstar */
			
			/* E[Xt+1 - Xt / Ft] and V[Xt+1 - Xt / Ft]*/
			DXDrift[iNumZ][iNumX] = DXDrift[iNumZ][iIndexMinX];
			DXVar[iNumZ][iNumX] = DXVar[iNumZ][iIndexMinX];
			
			/* E[Zt+1 - Zt / Ft] and V[Zt+1 - Zt / Ft]*/
			DZDrift[iNumZ][iNumX] = DZDrift[iNumZ][iIndexMinX];
			DZVar[iNumZ][iNumX] = DZVar[iNumZ][iIndexMinX];
		}
	}
}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEJumpUtTreatment												  

	

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEJumpUtTreatment(/* Inputs */
									int iIndexMinX,
									int	iIndexMaxX,
									int	iIndexMinZ,
									int	iIndexMaxZ,
									int iIndexMinProduct,
									int iIndexMaxProduct,
									int iNbPsi,
									
									double	*GridX,
									double	*GridZ,
								
									double dSigt,
									double dSig_t1,									
									double dAlphaVt,
									double dAlphaV_t1,
									double dRhot,
									double dRho_t1,

									int TypeInterpolation,
									int IsExtrapolationFlat,
									
									double ****ValueGridZtXtPsitProduct_t1,
									double **InterpQuadraCoef,

									/* Outputs */
									double ****ValueGridZtXtPsitProduct_t)
{
/* Declaration of locals variables */
double dX;
int iNumX, iNumZ, iNumPsi, iNumProduct;
int iIndexXplus, iIndex1, iIndex2, iIndex3;
double dX_t1;
	
	for (iNumZ=iIndexMinZ;iNumZ<=iIndexMaxZ;iNumZ++)
	{	
		/* Calculation of the (Xt+)-Xt */
		dX=(GridZ[iNumZ]-1) * (dRhot / dAlphaVt * dSigt - dRho_t1 / dAlphaV_t1 * dSig_t1);
			
		for (iNumX=iIndexMinX;iNumX<=iIndexMaxX;iNumX++)
		{			
			iIndexXplus = iIndexMinX-1;
			
			dX_t1 = GridX[iNumX]+dX;
			
			/* Find the index such that GridX[iIndexXplus]>= dX_t1 */
			while ((iIndexXplus<iIndexMaxX)&&(GridX[iIndexXplus+1]<dX_t1))
				iIndexXplus ++;
			
			if (IsExtrapolationFlat == LGMSV_TRUE)
			{
				/* Interpolate value */
				if (iIndexXplus == iIndexMinX-1)
				{
					/* Extrapolation flat */
					for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
						for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)	
							ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct] =
								ValueGridZtXtPsitProduct_t1[iNumZ][iIndexMinX][iNumPsi][iNumProduct];
				}
				else 
				if (iIndexXplus == iIndexMaxX)
				{
					/* Extrapolation flat */
					for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
						for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)	
							ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct]=
								ValueGridZtXtPsitProduct_t1[iNumZ][iIndexMaxX][iNumPsi][iNumProduct];
				}
				else if ((iIndexXplus>1+iIndexMinX)&&(iIndexXplus<iIndexMaxX-2)
					&&(TypeInterpolation==LGMSV_MIXQUADRALINEAR))
				{
					/* Quadratic Interpolation  */
					if (dX_t1-GridX[iIndexXplus]>0.5*(GridX[iIndexXplus+1]-GridX[iIndexXplus]))
					{
						iIndex1 = iIndexXplus;
						iIndex2 = iIndexXplus+1;
						iIndex3 = iIndexXplus+2;
						for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
							for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)
							{
								ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct] = 
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex1][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][0]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][1]*dX_t1+InterpQuadraCoef[iIndex1][2])+
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex2][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][3]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][4]*dX_t1+InterpQuadraCoef[iIndex1][5])+
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex3][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][6]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][7]*dX_t1+InterpQuadraCoef[iIndex1][8]);					
							}
					}
					else
					{
						iIndex1 = iIndexXplus-1;
						iIndex2 = iIndexXplus;
						iIndex3 = iIndexXplus+1;

						for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
							for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)
							{
								ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct] = 
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex1][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][0]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][1]*dX_t1+InterpQuadraCoef[iIndex1][2])+
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex2][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][3]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][4]*dX_t1+InterpQuadraCoef[iIndex1][5])+
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex3][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][6]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][7]*dX_t1+InterpQuadraCoef[iIndex1][8]);					
							}
					}
				}
				else
				{
					/* Linear Interpolation */
					for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
						for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)
						{
							ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct] = 
								(dX_t1-GridX[iIndexXplus])/(GridX[iIndexXplus+1]-GridX[iIndexXplus])*
								(ValueGridZtXtPsitProduct_t1[iNumZ][iIndexXplus+1][iNumPsi][iNumProduct]-
								ValueGridZtXtPsitProduct_t1[iNumZ][iIndexXplus][iNumPsi][iNumProduct])
								+ ValueGridZtXtPsitProduct_t1[iNumZ][iIndexXplus][iNumPsi][iNumProduct];
						}

				}
			}
			else
			{
				/* Interpolate value */
				if ((iIndexXplus>1+iIndexMinX)&&(iIndexXplus<iIndexMaxX-2)
					&&(TypeInterpolation==LGMSV_MIXQUADRALINEAR))
				{
					/* Quadratic Interpolation  */
					if (dX_t1-GridX[iIndexXplus]>0.5*(GridX[iIndexXplus+1]-GridX[iIndexXplus]))
					{
						iIndex1 = iIndexXplus;
						iIndex2 = iIndexXplus+1;
						iIndex3 = iIndexXplus+2;
						for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
							for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)
							{
								ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct] = 
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex1][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][0]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][1]*dX_t1+InterpQuadraCoef[iIndex1][2])+
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex2][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][3]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][4]*dX_t1+InterpQuadraCoef[iIndex1][5])+
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex3][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][6]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][7]*dX_t1+InterpQuadraCoef[iIndex1][8]);					
							}
					}
					else
					{
						iIndex1 = iIndexXplus-1;
						iIndex2 = iIndexXplus;
						iIndex3 = iIndexXplus+1;

						for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
							for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)
							{
								ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct] = 
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex1][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][0]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][1]*dX_t1+InterpQuadraCoef[iIndex1][2])+
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex2][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][3]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][4]*dX_t1+InterpQuadraCoef[iIndex1][5])+
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndex3][iNumPsi][iNumProduct]*InterpQuadraCoef[iIndex1][6]*
									(dX_t1*dX_t1-InterpQuadraCoef[iIndex1][7]*dX_t1+InterpQuadraCoef[iIndex1][8]);					
							}	
					}
				}
				else
				{
					
					/* Interpolate value */
					if (iIndexXplus==iIndexMinX-1)
					{
						iIndexXplus = iIndexMinX;		
					}	
					else
					if (iIndexXplus==iIndexMaxX)
					{
						/* Extrapolation flat */
						iIndexXplus= max(iIndexMaxX-1,0);
					
					}

					if (iIndexMaxX == iIndexMinX)
						/* Extrapolation flat */
							for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
								for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)
									ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct] = 
											ValueGridZtXtPsitProduct_t1[iNumZ][iIndexMinX][iNumPsi][iNumProduct];
					else
					{

						/* Linear Interpolation */
						for (iNumPsi=0;iNumPsi<iNbPsi;iNumPsi++)
							for (iNumProduct=iIndexMinProduct;iNumProduct<=iIndexMaxProduct;iNumProduct++)
							{
								ValueGridZtXtPsitProduct_t[iNumZ][iNumX][iNumPsi][iNumProduct] = 
									(dX_t1-GridX[iIndexXplus])/(GridX[iIndexXplus+1]-GridX[iIndexXplus])*
									(ValueGridZtXtPsitProduct_t1[iNumZ][iIndexXplus+1][iNumPsi][iNumProduct]-
									ValueGridZtXtPsitProduct_t1[iNumZ][iIndexXplus][iNumPsi][iNumProduct])
									+ ValueGridZtXtPsitProduct_t1[iNumZ][iIndexXplus][iNumPsi][iNumProduct];
							}
					}
				}
			}
		}
	}
}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVCalculationNbSigMax
																			
   ------------------------------------------------------------------------------------------------ */
void LGMSVCalculationNbSigMax( 
					/* Inputs */
					int	iNbTime,			
			
					double *StdVt,
					double *ImpliedRhoTimesStdVt,

					/* Parameters */
					LGMSVParam *Params,
					 
					 /* Outputs */
					double	*GridXNbSigmaMaxRight,
					double	*GridZNbSigmaMaxRight,
					double	*GridXNbSigmaMaxLeft,
					double	*GridZNbSigmaMaxLeft)
{
/* Declaration of locals variables */
int iNumTime;
double dStdVt, dRhoTimesStdVt;

double	dC0_XLeft = -3.45;
double	dC1_XLeft = -1.56;
double	dC2_XLeft = 0.57;

double	dC0_XRight = 4.12;
double	dC1_XRight = 1.74;
double	dC2_XRight = 0.86;

double	dC0_ZRight = 4.41;
double	dC1_ZRight = 3.49;
double	dC2_ZRight = 0.0;

	for (iNumTime = 0; iNumTime<iNbTime; iNumTime++)
	{
		/* Calculation of the alphaEq*sqrt(T) and Rho*alphaEq*sqrt(T) */
		dStdVt = StdVt[iNumTime];
		dRhoTimesStdVt = ImpliedRhoTimesStdVt[iNumTime];

		/* Gestion */
		if (Params->IsNbSigmaMaxXGridLeftAuto == LGMSV_TRUE)
		{
			GridXNbSigmaMaxLeft[iNumTime] = -(dC0_XLeft+dC1_XLeft*dStdVt+dC2_XLeft*dRhoTimesStdVt);
			if (GridXNbSigmaMaxLeft[iNumTime]<Params->iNbSigmaXGridLeft)
				GridXNbSigmaMaxLeft[iNumTime] = Params->iNbSigmaXGridLeft;
		}
		else
		{
			GridXNbSigmaMaxLeft[iNumTime] = Params->iNbSigmaMaxXGridLeft;
		}
		
		if (Params->IsNbSigmaMaxXGridRightAuto == LGMSV_TRUE)
		{
			GridXNbSigmaMaxRight[iNumTime] = (dC0_XRight+dC1_XRight*dStdVt+dC2_XRight*dRhoTimesStdVt);
			if (GridXNbSigmaMaxRight[iNumTime]<Params->iNbSigmaXGridRight)
				GridXNbSigmaMaxRight[iNumTime] = Params->iNbSigmaXGridRight;
		}
		else
		{
			GridXNbSigmaMaxRight[iNumTime] = Params->iNbSigmaMaxXGridRight;
		}
		
		if (Params->IsNbSigmaMaxZGridRightAuto == LGMSV_TRUE)
		{
			GridZNbSigmaMaxRight[iNumTime] = (dC0_ZRight+dC1_ZRight*dStdVt+dC2_ZRight*dRhoTimesStdVt);
			if (GridZNbSigmaMaxRight[iNumTime]<Params->iNbSigmaZGridRight)
				GridZNbSigmaMaxRight[iNumTime] = Params->iNbSigmaZGridRight;
		}
		else
		{
			GridZNbSigmaMaxRight[iNumTime] = Params->iNbSigmaMaxZGridRight;
		}

		GridZNbSigmaMaxLeft[iNumTime] = Params->iNbSigmaZGridLeft;
	}
}
/* ----------------------------------------------------------------------------------------------- 
	LGMSVSaveTensor
	
	  Save the tensor for analysis
																			
   ------------------------------------------------------------------------------------------------ */

void LGMSVSaveTensor_UtPieceWise( /* Inputs */
					 int iNbPhi,
					 int iNbX,
					 int iNbEps,
					 
					 double *GridPhi,
					 double *GridX,
					 double *GridZ,
					 double ****ValueTensor_t1)
{	
int  iNumPhi, iNumX, iNumEps;
FILE *fid;
	
	fid = fopen("C:\\toto.txt","wb");
	fwrite(&iNbPhi,1,sizeof(long),fid);
	fwrite(&iNbX,1,sizeof(long),fid);
	fwrite(&iNbEps,1,sizeof(long),fid);

	/* Save the grids */
	for (iNumPhi = 0; iNumPhi<iNbPhi; iNumPhi++)
		fwrite(&GridPhi[iNumPhi],1,sizeof(double),fid);
	for (iNumX = 0; iNumX<iNbX; iNumX++)
		fwrite(&GridX[iNumX],1,sizeof(double),fid);
	for (iNumEps = 0; iNumEps<iNbEps; iNumEps++)
		fwrite(&GridZ[iNumEps],1,sizeof(double),fid);


	for (iNumPhi = 0; iNumPhi<iNbPhi; iNumPhi++)
			for (iNumX = 0; iNumX<iNbX; iNumX++)
				for (iNumEps = 0; iNumEps<iNbEps; iNumEps++)
						fwrite(&ValueTensor_t1[iNumPhi][iNumX][iNumEps][0],1,sizeof(double),fid);
	fclose(fid);

}
/* ----------------------------------------------------------------------------------------------- 
	Main function 
	------------- 


	Modelisation															
	-------------														
																							
	The Model is under QTstar

	d(f(t,Tstar)-f(0,Tstar)) = u(t) * sqrt(Vt) * dW1(t)
	dPsit = u(t)^2*Vt dt
	dVt = (dLvlV-dLambdaV*Vt) * dt + dAlpha * sqrt(Vt)*dW2(t)
	<dW1(t),dW2(t)> = dRho * dt

	where 
		u(t) is piecewise constant
		Psit=Phit*exp(-2*dLambdaX*(Tstar-t))

	Orthogonalisation
	-----------------

	  Zt = Vt ;
	  Xt = f(t,Tstar)-f(0,Tstar) - u(t)*dRho/dAlpha*(Vt-1)
			is such that <dXt,dW2(t)>=0


	Reconstruction
	--------------


																			
   ------------------------------------------------------------------------------------------------ */

Err	 lgmSV_adi_UtPieceWise(	
					/*	Time Information  */
					int			iNbTime,					
					double		*dTime,
					double		*dDate,

					/*	Space Discretisation	*/
					int			iNbPsi,
					int			iNbX,
					int			iNbZ,
						
					/*	Model data Information	*/
					double		dLambdaX,
					double		dSigma,	/* No More Use */
					
					double		*SigTime,
					double		*Sig,
					int			iNbSigTime,
					
					double		*dLambdaV,
					double		*dLvlV,
					double		*dAlphaV,
					double		*dRho,

					/* Parameters */
					LGMSVPARAM	Params,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,
					
					/*	Market data */
					double		*Ifr,					
					char		*cYieldCurve,
					
					/*	Payoff function */
					Err (*payoff_func)(/* Time Information */
										double	dDate,
										double	dTime,
										void	*func_parm,
										
										/* Market data	*/										
										void	*cYieldCurve,
										
										/*	Model data Information	*/
										double	dLambdaX,
										
										/* Martingale Mesure QTstar information */
										double	dTstar,		/* In time */

										double  dAlphaV,
										double  dRho,
										double  dSigt,
										
										/* Grid data Information	*/
										int		iIndPsitMin,
										int		iIndPsitMax,
										int		iIndXtMin,
										int		iIndXtMax,
										int		iIndZtMin,
										int		iIndZtMax,
										
										double	*GridPsi,
										double	*GridX,
										double	*GridZ,
						
										/* Tensor of results to be updated		*/
										/* 4 dimensions : Psit,Xt,Zt,Product	*/
										int		*iNbProduct,
										double	****PayoffTensor),
					/*	Result */
					int			iNbProduct,
					int			iMaxNbProduct,
					double		*dProductArrayPv)
{
/* Declaration of locals variables */
double			*Sigt		= NULL,
				*dLambdaVt	= NULL,
				*dLvlVt		= NULL,
				*dAlphaVt	= NULL,
				*dRhot		= NULL;	

double			*GridX = NULL;
double			*GridZ = NULL;

double			*GridXNbSigmaMaxRight = NULL;
double			*GridZNbSigmaMaxRight = NULL;
double			*GridXNbSigmaMaxLeft = NULL;
double			*GridZNbSigmaMaxLeft = NULL;

double			*GridPsi_t = NULL;
double			*GridPsi_t1 = NULL;

double			*PsiMean = NULL;
double			*PsiStd = NULL;
double			*XtMean = NULL;
double			*XtStd = NULL;
double			*ZtMean = NULL;
double			*ZtStd = NULL;
double			*ImpliedRhoTimesStdVt = NULL;

double			**DXDrift = NULL;
double			**DZDrift = NULL;
double			**DXVar = NULL;
double			**DZVar = NULL;
double			**r = NULL;

double			***ValueGridXtZtProduct = NULL;
double			****ValueTensor_t = NULL;
double			****ValueTensor_t1 = NULL;

Err				err = NULL;
clock_t			t1, t2;
double			time_psi = 0.0,
				time_pde = 0.0,
				time_jmp = 0.0,
				time_pay = 0.0;

CNPDE_TEMP_2D_LGMSV_ADI	pdestr, *pde = &pdestr; 

int iNumTime, iNumPsi, iNumX, iNumZ, iNumProduct;
int iIndexZ0, iIndexX0, iIndexPsi;

double dTstar = Params->Tstar;
double dB0Tstar;

/* Variables for computation optimisation */
int iNbPsiGrid_t, iNbPsiGrid_t1;
int iIndexMinX, iIndexMaxX, iIndexMinZ, iIndexMaxZ;

double a1,a2,a3;
double **InterpQuadraCoef = NULL;
double **InterpQuadraCoefGridX = NULL;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */

	/* For computational time calculation				 */
	t1 = clock();

	/* grid should be at least a 3x3 in Xt , Zt			*/
	iNbX = max(iNbX,3);
	iNbZ = max(iNbZ,3);

	/* All Number of points in the grid must be odd		*/
	iNbX = ((int) (iNbX/2))*2 +1;
	iNbZ = ((int) (iNbZ/2))*2 +1;
	iNbPsi = ((int) (iNbPsi/2))*2 +1;
	
	/*  Memory allocations								 */
	Sigt = dvector(0, iNbTime-1);
	dLambdaVt = dvector(0, iNbTime-1);
	dLvlVt = dvector(0, iNbTime-1);
	dAlphaVt = dvector(0, iNbTime-1);
	dRhot = dvector(0, iNbTime-1);

	GridX = dvector(0, iNbX-1);
	GridZ = dvector(0, iNbZ-1);

	GridXNbSigmaMaxRight = dvector(0, iNbTime-1);
	GridZNbSigmaMaxRight = dvector(0, iNbTime-1);
	GridXNbSigmaMaxLeft = dvector(0, iNbTime-1);
	GridZNbSigmaMaxLeft = dvector(0, iNbTime-1);
	
	GridPsi_t = dvector(0, iNbPsi-1);	/* Grid of Psi at t   */
	GridPsi_t1 = dvector(0, iNbPsi-1);	/* Grid of Psi at t+1 */

	PsiMean = dvector(0, iNbTime-1);
	PsiStd = dvector(0, iNbTime-1);
	
	XtMean = dvector(0, iNbTime-1);
	XtStd = dvector(0, iNbTime-1);

	ZtMean = dvector(0, iNbTime-1);
	ZtStd = dvector(0, iNbTime-1);

	ImpliedRhoTimesStdVt = dvector(0, iNbTime-1);

	DXDrift = dmatrix(0, iNbZ-1, 0, iNbX-1);
	DZDrift = dmatrix(0, iNbZ-1, 0, iNbX-1);
	DXVar = dmatrix(0, iNbZ-1, 0, iNbX-1);
	DZVar = dmatrix(0, iNbZ-1, 0, iNbX-1);
	r = dmatrix(0, iNbZ-1, 0, iNbX-1);
	
	ValueGridXtZtProduct =f3tensor(	0, iNbZ-1,
									0, iNbX-1,
									0, iMaxNbProduct-1);
		
	ValueTensor_t =f4tensor(0, iNbZ-1,
							0, iNbX-1,
							0, iNbPsi-1,
							0, iMaxNbProduct-1);
	
	ValueTensor_t1 =f4tensor(	0, iNbZ-1,
								0, iNbX-1,
								0, iNbPsi-1,
								0, iMaxNbProduct-1);

	InterpQuadraCoef = dmatrix(0, iNbPsi-1, 0, 8);
	InterpQuadraCoefGridX = dmatrix(0, iNbX-1, 0, 8);
	

	/* Backward PDE */

	/* Init of the Backward PDE : Allocation of arrays    */	
	num_f_pde_init_2d_LGMSV_adi(pde,								
								iNbPsi,
								iNbZ,
								iNbX,
								iMaxNbProduct); 

	if (!pde)
	{
		err = "Memory allocation error (2) in LGMSVPDE";
		goto FREE_RETURN;
	}

	/* Gestion of allocation errors */
	if (!pde || !Sigt || !dLambdaVt || !dLvlVt || !dAlphaVt || !dRhot || !GridX || !GridZ || !GridPsi_t || !GridPsi_t1 || !PsiMean || !PsiStd 
		 || !XtMean || !XtStd || !ZtMean || !ZtStd|| !DXDrift || !DZDrift || !DXVar || !DZVar || !r 
		 || !ValueGridXtZtProduct || !ValueTensor_t || !ValueTensor_t1 || !InterpQuadraCoef || 
		 !InterpQuadraCoefGridX || !GridXNbSigmaMaxRight || !GridZNbSigmaMaxRight || 
		 !GridXNbSigmaMaxLeft || !GridZNbSigmaMaxLeft)
	{
		err = "Memory allocation error (1) in LGMSVPDE";
		goto FREE_RETURN;
	}

	/* Calculation of B(0,Tstar) */
	dB0Tstar  = swp_f_df (dDate[0], dDate[0]+dTstar*365.0, cYieldCurve);

	/* Expectations calculations */
	if (Params->iUseOldExpect)
	{
		err = LGMSVPDEExpectations_UtPieceWise_OldVersion(
								/* Inputs */
								 iNbTime,			
								 dTime,
								
								 /* For f(t,Tstar) diffusion and Xt*/
								 dLambdaX,
								 dSigma,
								 dTstar,
								 
								 /* definition of u(t) piecewise constant */
								 SigTime,
								 Sig,
								 iNbSigTime,

								 /* For Vt diffusion */
								 dLambdaV,
								 dLvlV,
								 dAlphaV,
								 dRho,
								 
								 /* Outputs */
								 Sigt,	/* Value of u(t) = sig(t) for all t in dtime */
								 dLambdaVt,
								 dLvlVt,
								 dAlphaVt,
								 dRhot,

								 PsiMean,
								 PsiStd,
								 XtMean,	
								 XtStd,		
								 ZtMean,
								 ZtStd,
								 ImpliedRhoTimesStdVt);
	}
	else
	{
		err = LGMSVPDEExpectations_UtPieceWise(
									/* Inputs */
									 iNbTime,			
									 dTime,
									
									 /* For f(t,Tstar) diffusion and Xt*/
									 dLambdaX,
									 dSigma,
									 dTstar,
									 
									 /* definition of u(t) piecewise constant */
									 SigTime,
									 Sig,
									 iNbSigTime,

									 /* For Vt diffusion */
									 dLambdaV,
									 dLvlV,
									 dAlphaV,
									 dRho,

									 /* Numerical */
									 Params->dMultiIntegMinTime,
									 
									 /* Outputs */
									 Sigt,	/* Value of u(t) = sig(t) for all t in dtime */
									 dLambdaVt,
									 dLvlVt,
									 dAlphaVt,
									 dRhot,

									 PsiMean,
									 PsiStd,
									 XtMean,	
									 XtStd,		
									 ZtMean,
									 ZtStd,
									 ImpliedRhoTimesStdVt);
	}

	if (err)
	{
		goto FREE_RETURN;
	}

	/* Verifying that ZtStd and XtStd >0 at time T */
	if ((ZtStd[iNbTime-1]==0)||(XtStd[iNbTime-1]==0))
	{
		err = "error in LGMSVPDE : Vt or ft is not stochastic";
		goto FREE_RETURN;
	}

	/* Calculation of the NbSigmaMax Right and left for X and Z Grid */

	/*if (dRho < 0.0)
	{
		Params->IsNbSigmaMaxXGridLeftAuto = LGMSV_FALSE;
	}*/

	LGMSVCalculationNbSigMax( 
					/* Inputs */
					iNbTime,			
			
					ZtStd,
					ImpliedRhoTimesStdVt,

					/* Parameters */
					Params,
					 
					 /* Outputs */
					GridXNbSigmaMaxRight,
					GridZNbSigmaMaxRight,
					GridXNbSigmaMaxLeft,
					GridZNbSigmaMaxLeft);
	
	/* Making the grid on X */	
	LGMSVPDEGridX_UtPieceWise(/* Inputs */
								iNbX,
								XtStd[iNbTime-1],
								GridXNbSigmaMaxLeft[iNbTime-1],
								Params->iNbSigmaXGridLeft,
								Params->iNbSigmaXGridRight,
								GridXNbSigmaMaxRight[iNbTime-1],
								Params->dPercentExtremePointXGrid,
								Params->IsEqSpacedExtremePointXGrid,
									 
							 /* Outputs */
								GridX,
								&iIndexX0);
					
	/* Making the grid on Z */			
	LGMSVPDEGridZ_UtPieceWise(/* Inputs */
								iNbZ,
								ZtMean[iNbTime-1],
								ZtStd[iNbTime-1],
								Params->iNbSigmaZGridLeft,
								Params->iNbSigmaZGridRight,
								GridZNbSigmaMaxRight[iNbTime-1],
								Params->dPercentExtremePointZGrid,
								Params->IsEqSpacedExtremePointZGrid,
									 
							 /* Outputs */
								GridZ,
								&iIndexZ0);

	
	/* Final Grid of Psi , Grid of Psi at time T */
	LGMSVPDEMakePsiGrid_UtPieceWise(/* Inputs */
						iNbPsi,
						PsiMean[iNbTime-1],
						PsiStd[iNbTime-1],
						Params->iNbSigmaPhiGrid,

						/* Parameter */
						Params->IsExtremePoints,
						Params->iNbSigmaExtremePoints,
						Params->IsPhiBoundExp,
								
						/* Outputs */
						GridPsi_t1,
						&iIndexPsi);

	
	/* Final Payoff valuation   */
	err = payoff_func(/* Time Information */
						dDate[iNbTime-1],
						dTime[iNbTime-1],
						func_parm_tab[iNbTime-1],
										
						/* Market data	*/										
						cYieldCurve,
										
						/*	Model data Information	*/
						dLambdaX,
						
						/* Martingale Mesure QTstar information */
						dTstar,		/* In time */
						dAlphaVt[iNbTime-1],
						dRhot[iNbTime-1],
						Sigt[iNbTime-1],
										
						/* Grid data Information	*/						
						0,
						iNbPsi-1,
						0,
						iNbX-1,
						0,
						iNbZ-1,

						GridPsi_t1,
						GridX,
						GridZ,
						
						/* Tensor of results to be updated		*/
						/* 4 dimensions : Psit,Xt,Zt,Product	*/
						&iNbProduct,
						ValueTensor_t1);

	/* Calculation of Coefficients for quadratic interpolation in GridX */
	/* Used for Jump treatment between t+ -> t- */
	for (iNumX=1; iNumX<iNbX-3; iNumX++)
	{
		a1=GridX[iNumX];
		a2=GridX[iNumX+1];
		a3=GridX[iNumX+2];

		InterpQuadraCoefGridX[iNumX][0]=1/(a1-a2)/(a1-a3);
		InterpQuadraCoefGridX[iNumX][1]=a2+a3;
		InterpQuadraCoefGridX[iNumX][2]=a2*a3;

		InterpQuadraCoefGridX[iNumX][3]=1/(a2-a3)/(a2-a1);
		InterpQuadraCoefGridX[iNumX][4]=a3+a1;
		InterpQuadraCoefGridX[iNumX][5]=a3*a1;

		InterpQuadraCoefGridX[iNumX][6]=1/(a3-a1)/(a3-a2);
		InterpQuadraCoefGridX[iNumX][7]=a1+a2;
		InterpQuadraCoefGridX[iNumX][8]=a1*a2;		
	}
	
	if (Params->VerifExpectation>0)
	{
		
		if (Params->VerifExpectation == 1 )
		{
			/* Pricing the E(Vt) */			
			for (iNumZ = 0; iNumZ<iNbZ; iNumZ++)
				for (iNumX = 0; iNumX<iNbX; iNumX++)
					for (iNumPsi = 0; iNumPsi<iNbPsi; iNumPsi++)
						for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
						ValueTensor_t1[iNumZ][iNumX][iNumPsi][iNumProduct] *=GridZ[iNumZ];
		}
		if (Params->VerifExpectation == 2 )
		{
			/* Pricing the E(Psit) */
			for (iNumZ = 0; iNumZ<iNbZ; iNumZ++)
				for (iNumX = 0; iNumX<iNbX; iNumX++)
					for (iNumPsi = 0; iNumPsi<iNbPsi; iNumPsi++)
						for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
						ValueTensor_t1[iNumZ][iNumX][iNumPsi][iNumProduct] *=GridPsi_t1[iNumPsi];
		}
		if (Params->VerifExpectation == 3 )
		{
			/* Pricing the E(Vt^2) */
			for (iNumZ = 0; iNumZ<iNbZ; iNumZ++)
				for (iNumX = 0; iNumX<iNbX; iNumX++)
					for (iNumPsi = 0; iNumPsi<iNbPsi; iNumPsi++)
						for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
						ValueTensor_t1[iNumZ][iNumX][iNumPsi][iNumProduct] *=GridZ[iNumZ]*GridZ[iNumZ];
		}
		if (Params->VerifExpectation == 4 )
		{
			/* Pricing the E(ft) */
			for (iNumZ = 0; iNumZ<iNbZ; iNumZ++)
				for (iNumX = 0; iNumX<iNbX; iNumX++)
					for (iNumPsi = 0; iNumPsi<iNbPsi; iNumPsi++)
						for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
						ValueTensor_t1[iNumZ][iNumX][iNumPsi][iNumProduct] *=
							GridX[iNumX]+dRhot[iNbTime-1]*Sigt[iNbTime-1]/dAlphaVt[iNbTime-1]*(GridZ[iNumZ]-1.0);
		}
		if (Params->VerifExpectation == 5 )
		{
			/* Pricing the E(ft^2) */
			for (iNumZ = 0; iNumZ<iNbZ; iNumZ++)
				for (iNumX = 0; iNumX<iNbX; iNumX++)
					for (iNumPsi = 0; iNumPsi<iNbPsi; iNumPsi++)
						for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
						ValueTensor_t1[iNumZ][iNumX][iNumPsi][iNumProduct] *=
							pow(GridX[iNumX]+dRhot[iNbTime-1]*Sigt[iNbTime-1]/dAlphaVt[iNbTime-1]*(GridZ[iNumZ]-1.0),2);
		}
		if (Params->VerifExpectation == 6 )
		{
			/* Pricing the E(Zt^2) */
			for (iNumZ = 0; iNumZ<iNbZ; iNumZ++)
				for (iNumX = 0; iNumX<iNbX; iNumX++)
					for (iNumPsi = 0; iNumPsi<iNbPsi; iNumPsi++)
						for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
						ValueTensor_t1[iNumZ][iNumX][iNumPsi][iNumProduct] *=GridZ[iNumZ]*GridZ[iNumZ];
		}
		if (Params->VerifExpectation == 7 )
		{
			/* Pricing the V<0 */
			for (iNumZ = 0; iNumZ<iNbZ; iNumZ++)
				for (iNumX = 0; iNumX<iNbX; iNumX++)
					for (iNumPsi = 0; iNumPsi<iNbPsi; iNumPsi++)
						for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
						ValueTensor_t1[iNumZ][iNumX][iNumPsi][iNumProduct] *=(GridZ[iNumZ]<-1e-8);
		}
	}

	/* Init Nb Psi at t+1 */
	iNbPsiGrid_t1 = iNbPsi;

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	t2 = clock();
	
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Phase 2 -convolution, NbSigXRightMax: %f ", GridXNbSigmaMaxRight[iNbTime-1]);
	smessage ("Phase 2 -convolution, NbSigXLeftMax: %f ", GridXNbSigmaMaxLeft[iNbTime-1]);
	smessage ("Phase 2 -convolution, NbSigZRightMax: %f ", GridZNbSigmaMaxRight[iNbTime-1]);

	time_pay = t2 - t1;

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	/* Initial Index value */
	iIndexMinX = 0;
	iIndexMaxX = iNbX-1;
	iIndexMinZ = 0;
	iIndexMaxZ = iNbZ-1;

if (Params->ModeExport!=LGMSV_TRUE)
{	
	/* For each time t backward */
	for (iNumTime=iNbTime-2; iNumTime>=0; iNumTime--)
	{
		/* Update Parameters for calculation optimisation */
		
		/* Number of Psi in the grid at time t */
		iNbPsiGrid_t = iNbPsi-(Params->GridPhiOptim == LGMSV_TRUE)*(int)((1-iNumTime/((double)(iNbTime-1)))
			*(iNbPsi-1) + 1.0e-08);
		iNbPsiGrid_t = ((int) (iNbPsiGrid_t/2))*2 +1;
		
		
		if (Params->GridXADIOptim == LGMSV_TRUE)
		{
			while ((iIndexMinX <iIndexX0-1)&&(GridX[iIndexMinX+1]<-GridXNbSigmaMaxLeft[iNumTime+1]*XtStd[iNumTime+1]))
				iIndexMinX ++;

			while ((iIndexMaxX >iIndexX0+1)&&(GridX[iIndexMaxX-1]>GridXNbSigmaMaxRight[iNumTime+1]*XtStd[iNumTime+1]))
				iIndexMaxX --;
		}

		if (Params->GridZADIOptim == LGMSV_TRUE)
		{
			while ((iIndexMinZ <iIndexZ0-1)&&(GridZ[iIndexMinZ+1]<1.0-GridZNbSigmaMaxLeft[iNumTime+1]*ZtStd[iNumTime+1]))
				iIndexMinZ ++;

			while ((iIndexMaxZ >iIndexZ0+1)&&(GridZ[iIndexMaxZ-1]>1.0+GridZNbSigmaMaxRight[iNumTime+1]*ZtStd[iNumTime+1]))
				iIndexMaxZ --;
		}

		if (Params->SaveTensor == LGMSV_TRUE)
		{
			/* Save the Tensor */
			LGMSVSaveTensor_UtPieceWise( /* Inputs */
					iNbPsi,
					iNbX,
					iNbZ,
					 
					GridPsi_t1,
					GridX,
					GridZ,
					ValueTensor_t1);
		}
		
		/* Pre Calculation for code optimisation */
		for (iNumPsi=1; iNumPsi<iNbPsiGrid_t1-3; iNumPsi++)
		{
			a1=GridPsi_t1[iNumPsi];
			a2=GridPsi_t1[iNumPsi+1];
			a3=GridPsi_t1[iNumPsi+2];

			InterpQuadraCoef[iNumPsi][0]=1/(a1-a2)/(a1-a3);
			InterpQuadraCoef[iNumPsi][1]=a2+a3;
			InterpQuadraCoef[iNumPsi][2]=a2*a3;

			InterpQuadraCoef[iNumPsi][3]=1/(a2-a3)/(a2-a1);
			InterpQuadraCoef[iNumPsi][4]=a3+a1;
			InterpQuadraCoef[iNumPsi][5]=a3*a1;

			InterpQuadraCoef[iNumPsi][6]=1/(a3-a1)/(a3-a2);
			InterpQuadraCoef[iNumPsi][7]=a1+a2;
			InterpQuadraCoef[iNumPsi][8]=a1*a2;		
		}

		/* Make the grid in Psit at time t */
		LGMSVPDEMakePsiGrid_UtPieceWise(/* Inputs */
						iNbPsiGrid_t,
						PsiMean[iNumTime],
						PsiStd[iNumTime],
						Params->iNbSigmaPhiGrid,

						/* Parameter */
						Params->IsExtremePoints,
						Params->iNbSigmaExtremePoints,
						Params->IsPhiBoundExp,
								
						/* Outputs */
						GridPsi_t,
						&iIndexPsi);

		
		/* Calculation of the drifts and variances under Qt+1 */
		LGMSVPDECalculationRDriftVar_UtPieceWise(/* Inputs */
										iIndexMinX,
										iIndexMaxX,
										iIndexMinZ,
										iIndexMaxZ,
										
										GridX,
										GridZ,
									
										Sigt[iNumTime+1],

										dLambdaVt[iNumTime+1],
										dLvlVt[iNumTime+1],
										dAlphaVt[iNumTime+1],
										dRhot[iNumTime+1],

										dTime[iNumTime], /*time(t) */
										dTime[iNumTime+1], /*time(t+1) */
																				
										/* Outputs */
										DXDrift,
										DZDrift,
										DXVar,
										DZVar,
										r);
		
		/* Calculation of the value at t+1 on the grid (X,Z)			 */
		/* at the deduced forward Psi t+1 from Psit						 */
		/* Interpolation from the value in the grid (X,Z,Psi) at t+1	 */

		t1 = clock();

		LGMSVPDECalculationValueMatrixfromPsit_UtPieceWise(
									 /* Inputs */
									 iIndexMinZ,
									 iIndexMaxZ,
									 iIndexMinX,
									 iIndexMaxX,										 
									 
									 iNbPsiGrid_t,
									 iNbPsiGrid_t1,
									 iNbProduct,
									 
									 GridX, 
									 GridZ,
									 GridPsi_t,
									 GridPsi_t1,

									 /* Martingale Mesure QTstar information */
									 dTstar,	/* In time */

									 /* 4 dimensions : Zt, Xt , Phit, Product	*/
									 ValueTensor_t1,

									 Sigt[iNumTime+1],
									 dLambdaX,
									 dAlphaVt[iNumTime+1],
									 dRhot[iNumTime+1],

									 dTime[iNumTime], /*time(t) */
									 dTime[iNumTime+1], /*time(t+1) */

									 Params->TypeInterpolationPhi,
									 Params->IsExtrapolationFlat,

									 InterpQuadraCoef,

									 /* Output */
									 /* 4 dimensions : Zt, Xt , Phit, Product	*/
									 ValueTensor_t);

		t2 = clock();
		time_psi += (double) (t2 - t1);

		/* One Step Backward 2D using Crank-Nicolson on grid (X,Z) */	
		/* From t+1 to t										   */
		/* Launch The Crank-Nicolson					*/
		/* Calculation of the payoff tensor at time t	*/
		/* for the correspondant Psi(t)					*/

		num_f_pde_one_step_backward_2f_save_LGMSV_adi(	pde,														
														iNbZ,
														GridZ,
														iNbX,
														GridX,
														0,
														iNbProduct-1,
														ValueTensor_t,
														DZDrift,
														DXDrift,
														DZVar,
														DXVar,
														r,
														ValueTensor_t1,
														0,
														iNbPsiGrid_t-1,
														iIndexMinZ,
														iIndexMaxZ,
														iIndexMinX,
														iIndexMaxX);

		t1 = clock();
		time_pde += (double) (t1 - t2);

		/* Treatment if jump of u(t) between t+ -> t- */
		if ((fabs(Sigt[iNumTime+1]-Sigt[iNumTime]) > 1e-8) ||
			(fabs(dAlphaVt[iNumTime+1]-dAlphaVt[iNumTime]) > 1e-8) ||
			(fabs(dRhot[iNumTime+1]-dRhot[iNumTime]) > 1e-8))
		{
			LGMSVPDEJumpUtTreatment(/* Inputs */
										iIndexMinX,
										iIndexMaxX,
										iIndexMinZ,
										iIndexMaxZ,
										0,
										iNbProduct-1,
										iNbPsiGrid_t,
										
										GridX,
										GridZ,
									
										Sigt[iNumTime],
										Sigt[iNumTime+1],
										dAlphaVt[iNumTime],
										dAlphaVt[iNumTime+1],
										dRhot[iNumTime],
										dRhot[iNumTime+1],

										LGMSV_MIXQUADRALINEAR, /*LGMSV_LINEAR, /*LGMSV_MIXQUADRALINEAR, */
										LGMSV_TRUE,
										
										ValueTensor_t1,
										InterpQuadraCoefGridX,

										/* Outputs */
										ValueTensor_t);

			LGMSVPDESwap(double ****,ValueTensor_t,ValueTensor_t1);
		}

		t2 = clock();
		time_jmp += (double) (t2 - t1);

		/*	Case of evaluation events */
		if (EvalEvent[iNumTime])
		{
			t1 = clock();

			/* Modification of the Payoff at t */			
			err = payoff_func(/* Time Information */
								dDate[iNumTime],
								dTime[iNumTime],
								func_parm_tab[iNumTime],
												
								/* Market data	*/										
								cYieldCurve,
												
								/*	Model data Information	*/
								dLambdaX,
								
								/* Martingale Mesure QTstar information */
								dTstar,		/* In time */
								dAlphaVt[iNumTime],
								dRhot[iNumTime],
								Sigt[iNumTime],
												
								/* Grid data Information	*/
								0,
								iNbPsiGrid_t-1,
								iIndexMinX,
								iIndexMaxX,
								iIndexMinZ,
								iIndexMaxZ,
												
								GridPsi_t,
								GridX,
								GridZ,
								
								/* Tensor of results to be updated		*/
								/* 4 dimensions : Psit,ft,Epst,Product	*/
								&iNbProduct,
								ValueTensor_t1);
			if (err)
			{
					goto FREE_RETURN;
			}

			t2 = clock();

			time_pay += t2 - t1;
		}

		/* Swaping, t+1 becomes t */
		//LGMSVPDESwap(double ****,ValueTensor_t,ValueTensor_t1);
		LGMSVPDESwap(double *,GridPsi_t,GridPsi_t1);
		LGMSVPDESwap(int ,iNbPsiGrid_t,iNbPsiGrid_t1);
	}
	/* End for each t */

	/* Multiplication of the results by B(0,Tstar) to obtain the PV */
	for (iNumProduct = 0; iNumProduct<iNbProduct; iNumProduct++)
			ValueTensor_t1[iIndexZ0][iIndexX0][iIndexPsi][iNumProduct] *= dB0Tstar;

	/* Copy of the result */
	/* The grids are made such as 
		center on Psi = 0, X=0, Z=1  at t=0 */
	memcpy(dProductArrayPv,ValueTensor_t1[iIndexZ0][iIndexX0][iIndexPsi],iNbProduct*sizeof(double));

	/* Convolution time display */
	t1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (t1 - t2) / CLOCKS_PER_SEC);

	smessage ("time psi %.2f", (double) time_psi / CLOCKS_PER_SEC);
	smessage ("time pde %.2f", (double) time_pde / CLOCKS_PER_SEC);
	smessage ("time jump %.2f", (double) time_jmp / CLOCKS_PER_SEC);
	smessage ("time payoff %.2f", (double) time_pay / CLOCKS_PER_SEC);
}
else
{
	dProductArrayPv[0]=	PsiStd[iNbTime-1]; 
	dProductArrayPv[1]= PsiMean[iNbTime-1]; 
}


FREE_RETURN:

	/* free all the arrays */
	if (pde) num_f_pde_free_2d_LGMSV_adi(pde, iNumPsi, iNbX, iNbZ, iNbProduct); 
	
	if (Sigt)
		free_dvector(Sigt, 0, iNbTime-1);

	if (dLambdaVt)
		free_dvector(dLambdaVt, 0, iNbTime-1);

	if (dLvlVt)
		free_dvector(dLvlVt, 0, iNbTime-1);

	if (dAlphaVt)
		free_dvector(dAlphaVt, 0, iNbTime-1);

	if (dRhot)
		free_dvector(dRhot, 0, iNbTime-1);
	
	if (GridX)
		free_dvector(GridX, 0, iNbX-1);

	if (GridZ)
		free_dvector(GridZ, 0, iNbZ-1);

	if (GridXNbSigmaMaxRight)
		free_dvector(GridXNbSigmaMaxRight, 0, iNbTime-1);

	if (GridZNbSigmaMaxRight)
		free_dvector(GridZNbSigmaMaxRight, 0, iNbTime-1);

	if (GridXNbSigmaMaxLeft)
		free_dvector(GridXNbSigmaMaxLeft, 0, iNbTime-1);

	if (GridZNbSigmaMaxLeft)
		free_dvector(GridZNbSigmaMaxLeft, 0, iNbTime-1);

	if (GridPsi_t)
		free_dvector(GridPsi_t, 0, iNbPsi-1);
	if (GridPsi_t1)
		free_dvector(GridPsi_t1, 0, iNbPsi-1);

	if (PsiMean)
		free_dvector(PsiMean, 0, iNbTime-1);
	if (PsiStd)
		free_dvector(PsiStd, 0, iNbTime-1);

	if (XtMean)
		free_dvector(XtMean, 0, iNbTime-1);
	if (XtStd)
		free_dvector(XtStd, 0, iNbTime-1);
	
	if (ZtMean)
		free_dvector(ZtMean, 0, iNbTime-1);
	if (ZtStd)
		free_dvector(ZtStd, 0, iNbTime-1);
	if (ImpliedRhoTimesStdVt)
		free_dvector(ImpliedRhoTimesStdVt, 0, iNbTime-1);

	if (DXDrift)
		free_dmatrix(DXDrift, 0, iNbX-1, 0, iNbZ-1);
	if (DZDrift)
		free_dmatrix(DZDrift, 0, iNbX-1, 0, iNbZ-1);
	if (DXVar)
		free_dmatrix(DXVar, 0, iNbX-1, 0, iNbZ-1);
	if (DZVar)
		free_dmatrix(DZVar, 0, iNbX-1, 0, iNbZ-1);
	if (r)
		free_dmatrix(r, 0, iNbX-1, 0, iNbZ-1);
	
	if (ValueGridXtZtProduct)
		free_f3tensor(ValueGridXtZtProduct, 0, iNbZ-1, 0, iNbX-1, 0, iNbProduct-1);
		
	if (ValueTensor_t)
		free_f4tensor(ValueTensor_t, 0, iNbZ-1, 0, iNbX-1, 0, iNbPsi-1, 0, iNbProduct-1);
	
	if (ValueTensor_t1)
		free_f4tensor(ValueTensor_t1, 0, iNbZ-1, 0, iNbX-1, 0, iNbPsi-1, 0, iNbProduct-1);
			
	
	if (InterpQuadraCoef)
		free_dmatrix(InterpQuadraCoef, 0, iNbPsi-1, 0, 8);

	if (InterpQuadraCoefGridX)
		free_dmatrix(InterpQuadraCoefGridX, 0, iNbX-1, 0, 8);




	/* Return the error message */
	return err;

}

Err Fill_lgmSV_defaultParam(LGMSVPARAM	Params)
{

	Params->IsQTstarModel = LGMSV_IsQTstarModel;

	// X Grid
	Params->IsNbSigmaMaxXGridRightAuto = LGMSV_IsNbSigmaMaxXGridRightAuto;
	Params->iNbSigmaMaxXGridRight = LGMSV_iNbSigmaMaxXGridRight;
	Params->iNbSigmaXGridRight = LGMSV_iNbSigmaXGridRight;
	
	Params->IsNbSigmaMaxXGridLeftAuto = LGMSV_IsNbSigmaMaxXGridLeftAuto;
	Params->iNbSigmaMaxXGridLeft = LGMSV_iNbSigmaMaxXGridLeft;
	Params->iNbSigmaXGridLeft = LGMSV_iNbSigmaXGridLeft;
	
	Params->dPercentExtremePointXGrid = LGMSV_dPercentExtremePointXGrid;
	Params->IsEqSpacedExtremePointXGrid = LGMSV_IsEqSpacedExtremePointXGrid;

	// Use for Old Version
	Params->IsNbSigmaXGridRightAuto  = LGMSV_IsNbSigmaXGridRightAuto;

	// Z Grid
	Params->IsNbSigmaMaxZGridRightAuto = LGMSV_IsNbSigmaMaxZGridRightAuto;
	Params->iNbSigmaMaxZGridRight = LGMSV_iNbSigmaMaxZGridRight;
	Params->iNbSigmaZGridRight = LGMSV_iNbSigmaZGridRight;
	Params->iNbSigmaZGridLeft = LGMSV_iNbSigmaZGridLeft;
	
	Params->dPercentExtremePointZGrid = LGMSV_dPercentExtremePointZGrid;
	Params->IsEqSpacedExtremePointZGrid = LGMSV_IsEqSpacedExtremePointZGrid;

	// Phi Grid
	Params->TypeInterpolationPhi = LGMSV_TypeInterpolationPhi;
	Params->IsExtrapolationFlat = LGMSV_IsExtrapolationFlat;

	Params->iNbSigmaPhiGrid = LGMSV_iNbSigmaPhiGrid;
	Params->IsExtremePoints = LGMSV_IsExtremePoints;
	Params->iNbSigmaExtremePoints = LGMSV_iNbSigmaExtremePoints;
	Params->IsPhiBoundExp = LGMSV_IsPhiBoundExp;

	// Jump
	Params->TypeInterpolationJump = LGMSV_TypeInterpolationJump;
	Params->IsExtrapolationFlatJump = LGMSV_IsExtrapolationFlatJump;

	Params->VerifExpectation = LGMSV_VerifExpectation;

	// Debug
	Params->ModeExport = LGMSV_ModeExport;
	Params->SaveTensor = LGMSV_SaveTensor;

	// Use For Old Version
	Params->VFloor = LGMSV_VFloor;

	/* For Calculation Optimisation */
	Params->GridXADIOptim = LGMSV_GridXADIOptim;
	Params->GridZADIOptim = LGMSV_GridZADIOptim;
	Params->GridPhiOptim = LGMSV_GridPhiOptim;

	/* For MC */
	Params->UseBalsamGen = LGMSV_TRUE;
	Params->UseReverseMC = 2;

	Params->Tstar = LGMSV_Tstar;
	Params->UseNewTStarMC = LGMSV_FALSE;
	Params->iSchemeOrder = 1;

	Params->dMultiIntegMinTime = 1.0 / 4.0;
	Params->iUseOldExpect = LGMSV_FALSE;

	return NULL;
}