/* ==========================================================================
   FILE_NAME:	LGMSVPDE_Old.c

   PURPOSE:		PDE implementation of the 1 Factor LGM model with stochastic vol.
	
   DATE:		11/28/01
   
   AUTHOR:		P.A.
   ========================================================================== */

#include "LGMSVPDE.h"
#include "LGMSVGrfn.h"
#include "LGMSVUtil.h"

/* For swaping two variables */
#define LGMSVPDESwap(Type,A,B) {Type C; C=A; A=B; B=C;}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEExpectations														  

	Calculation of expectations for the definition of the grids 
		in Xt (constant grid center in XTMean =E(XT / X0) and +- iNbSigmaXGrid XTStd)
		in Epst (constant grid center in EpsTMean =E(EpsT / Eps0) and +- iNbSigmaEpsGrid EpsStd)
		in Phit (variable grid center at each t on PhitMean(t) 

   ------------------------------------------------------------------------------------------------ */
static Err LGMSVPDEExpectations(
					/* Inputs */
					 int	iNbTime,			
					 double	*dTime,
					
					 /* For X diffusion */
					 double dLambdaX,
					 double dSigma,
					 double dTstar,
					 
					 double *SigTime,
					 double	*Sig,
					 int	iNbSigTime,

					 /* For Eps diffusion */
					 double dLambdaEps,
					 double dAlpha,
					 double	dRho,
					 
					 /* Outputs */
					 double *Sigt,			/* Value of sig(t) for all t in dtime */
					 
					 double *PhitMean,
					 double *PhitStd,
					 
					 double *ftMean,
					 double *ftStd,

					 double *ZtMean,
					 double *ZtStd)
{	
/* Declaration of locals variables */
int iNumTime, iNumSigTime;
double dDt;

double lambda[100];

LGMSVSolFunc FuncDerBeta, FuncDerBeta2, FuncVarft, FuncVt2Mean, FuncPhitMean;
LGMSVSolFunc FuncVtPhitMean, FuncPhit2Mean, FuncVtftMean;

double dDerBeta, dDerBeta2, dVarft, dVt2Mean, dPhitMean;
double dVtPhitMean, dPhit2Mean, dVtftMean;

double dg2;
Err	   err = NULL;

	memset(lambda, 0, 100 * sizeof(double));

	/* Initialisation of the LGMSVSolFunc */
	
	/* Definition of the function (dBeta/DT(t,Tstar)) */
	FuncDerBeta.dLambda = -dLambdaX;
	FuncDerBeta.bIsft1 = 1;
	FuncDerBeta.a =0;
	FuncDerBeta.bIsgt1 = 1;
	FuncDerBeta.b = 0;
	FuncDerBeta.bIsht1 = 1;
	FuncDerBeta.c =0;
	
	/* Definition of the function (dBeta/DT(t,Tstar))^2 */
	FuncDerBeta2.dLambda = -2*dLambdaX;
	FuncDerBeta2.bIsft1 = 1;
	FuncDerBeta2.a =0;
	FuncDerBeta2.bIsgt1 = 1;
	FuncDerBeta2.b = 0;
	FuncDerBeta2.bIsht1 = 1;
	FuncDerBeta2.c =0;

	
	/* Definition of the function Var(f(t,Tstar)-f(0,Tstar)) */
	FuncVarft.dLambda = 0;
	FuncVarft.bIsft1 = 0;
	FuncVarft.pft = &FuncDerBeta2;
	FuncVarft.bIsgt1 = 1;
	FuncVarft.b = 0;
	FuncVarft.bIsht1 = 1;
	FuncVarft.c =0;

	/* Definition of the function E(Vt^2) */
	FuncVt2Mean.dLambda = 2*dLambdaEps;
	FuncVt2Mean.bIsft1 = 1;
	FuncVt2Mean.a = dAlpha*dAlpha+2*dLambdaEps;
	FuncVt2Mean.bIsgt1 = 1;
	FuncVt2Mean.b = 0;
	FuncVt2Mean.bIsht1 = 1;
	FuncVt2Mean.c =0;

	/* Definition of the function E(Phit) */
	FuncPhitMean.dLambda = 2*dLambdaX;
	FuncPhitMean.bIsft1 = 1;
	FuncPhitMean.bIsgt1 = 1;
	FuncPhitMean.b = 0;
	FuncPhitMean.bIsht1 = 1;
	FuncPhitMean.c =0;

	/* Definition of the function E(VtPhit) */
	FuncVtPhitMean.dLambda = 2*dLambdaX+dLambdaEps;
	FuncVtPhitMean.bIsft1 = 0;
	FuncVtPhitMean.pft = &FuncVt2Mean;
	FuncVtPhitMean.bIsgt1 = 0;
	FuncVtPhitMean.b = dLambdaEps;
	FuncVtPhitMean.pgt = &FuncPhitMean;
	FuncVtPhitMean.bIsht1 = 1;
	FuncVtPhitMean.c =0;


	/* Definition of the function E(Phit^2) */
	FuncPhit2Mean.dLambda = 4*dLambdaX;
	FuncPhit2Mean.bIsft1 = 0;
	FuncPhit2Mean.pft = &FuncVtPhitMean;
	FuncPhit2Mean.bIsgt1 = 1;
	FuncPhit2Mean.b = 0;
	FuncPhit2Mean.bIsht1 = 1;
	FuncPhit2Mean.c =0;


	/* Definition of the function E(Vt*(f(t,Tstar)-f(0,Tstar))) */
	FuncVtftMean.dLambda = dLambdaEps;
	FuncVtftMean.bIsft1 = 0;
	FuncVtftMean.pft = &FuncDerBeta;
	FuncVtftMean.bIsgt1 = 1;
	FuncVtftMean.b = 0;
	FuncVtftMean.bIsht1 = 1;
	FuncVtftMean.c =0;

	/* Filling the Sigt array 
	   piecewise interpolation from the Sig array on SigTime */
	
	iNumSigTime = 0;
	for (iNumTime = 0; iNumTime<iNbTime; iNumTime++)
	{
		if ((dTime[iNumTime]>SigTime[iNumSigTime]) && (iNumSigTime<iNbSigTime-1))
			iNumSigTime++;
		
		Sigt[iNumTime] = Sig[iNumSigTime];	
	}

	
	/* Init Value at t=0 */
	dDerBeta = exp(-dLambdaX*dTstar);
	dDerBeta2 = exp(-2*dLambdaX*dTstar);
	dVarft = 0;
	dVt2Mean = 1;
	dPhitMean = 0;
	dVtPhitMean = 0;
	dPhit2Mean = 0;
	dVtftMean = 0;
	
	PhitMean[0] = dPhitMean;
	PhitStd[0] = sqrt(dPhit2Mean-dPhitMean*dPhitMean);

	ZtMean[0] = 1;
	ZtStd[0] = 0;

	ftMean[0] = 0;
	ftStd[0] = 0;
	
	/* Calculation of all the mean and variances */
	for (iNumTime = 1;iNumTime<iNbTime;iNumTime++)
	{
		/* Delta time between t and t-1 */
		dDt = dTime[iNumTime]-dTime[iNumTime-1];
		
		/* Calculation of g(t-1)^2 */
		dg2 = Sigt[iNumTime-1]*Sigt[iNumTime-1];

		/* Update value at t-1 on all func */
		FuncDerBeta.dXt1 = dDerBeta; /*exp(-dLambdaX*(dTstar-dTime[iNumTime-1]));*/

		FuncDerBeta2.dXt1 = dDerBeta2; /*exp(-2*dLambdaX*(dTstar-dTime[iNumTime-1]));*/
		
		FuncVarft.dXt1 = dVarft;
		FuncVarft.a=dg2;

		FuncVt2Mean.dXt1 = dVt2Mean;

		FuncPhitMean.dXt1=dPhitMean;
		FuncPhitMean.a=dg2;

		FuncVtPhitMean.dXt1=dVtPhitMean;
		FuncVtPhitMean.a=dg2;

		FuncPhit2Mean.dXt1=dPhit2Mean;
		FuncPhit2Mean.a=2*dg2;

		FuncVtftMean.dXt1=dVtftMean;
		FuncVtftMean.a = dAlpha*dRho*Sigt[iNumTime-1];


		/* Calculation at time t */

		/* (dBeta/DT(t,Tstar)) at t knowing t-1 */
		dDerBeta = LGMSVFuncValue(FuncDerBeta,dDt,lambda,0);

		/* (dBeta/DT(t,Tstar))^2 at t knowing t-1 */
		dDerBeta2 = LGMSVFuncValue(FuncDerBeta2,dDt,lambda,0);

		/* Var(f(t,Tstar)-f(0,Tstar)) at t knowing t-1 */
		dVarft = LGMSVFuncValue(FuncVarft,dDt,lambda,0);

		/* E(Vt^2) at t knowing t-1 */
		dVt2Mean = LGMSVFuncValue(FuncVt2Mean,dDt,lambda,0);

		/* E(Phit) at t knowing t-1 */
		dPhitMean = LGMSVFuncValue(FuncPhitMean,dDt,lambda,0);
		
		/* E(VtPhit) at t knowing t-1 */
		dVtPhitMean = LGMSVFuncValue(FuncVtPhitMean,dDt,lambda,0);

		/* E(Phit^2) at t knowing t-1 */
		dPhit2Mean = LGMSVFuncValue(FuncPhit2Mean,dDt,lambda,0);

		/* E(Vt*(f(t,Tstar)-f(0,Tstar))) at t knowing t-1 */
		dVtftMean = LGMSVFuncValue(FuncVtftMean,dDt,lambda,0);

		/* Mean of Phi at t knowing t-1 */
		PhitMean[iNumTime] = dPhitMean;
		PhitStd[iNumTime] = sqrt(dPhit2Mean-dPhitMean*dPhitMean);
		
		/* Mean of f(t,Tstar)-f(0,Tstar) at t knowing t-1 */
		ftMean[iNumTime] = 0;

		/* Std of f(t,Tstar)-f(0,Tstar) at t knowing t-1 */
		ftStd[iNumTime] = sqrt(dVarft);

		/* Calculation of Zt Mean */
		ZtMean[iNumTime] = 1-dAlpha*dRho/Sigt[iNumTime]/dDerBeta*ftMean[iNumTime];

		/* Calculation of Zt Mean */
		ZtStd[iNumTime] = sqrt(dVt2Mean-1 + pow(dAlpha*dRho/(Sigt[iNumTime]*dDerBeta),2)*dVarft
		-2*dAlpha*dRho/(Sigt[iNumTime]*dDerBeta)*dVtftMean);

		/* We want that ZtStd increases */
		ZtStd[iNumTime] = max (ZtStd[iNumTime], ZtStd[iNumTime-1]);

		/* Error Check */
		if ((dPhit2Mean-dPhitMean*dPhitMean<0)||(dVarft <0) ||(dVt2Mean-1 + pow(dAlpha*dRho/(Sigt[iNumTime]*dDerBeta),2)*dVarft
			-2*dAlpha*dRho/(Sigt[iNumTime]*dDerBeta)*dVtftMean<0))
		{
			err ="Err LGMSVPDE : Numerical Pb in Std Calculation, Tau could be too big";
		}
	}

	return err;
}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEGridUnif

	compute a uniform (equally spaced) grid from dMean-iNbSigmaLeft*dStd to dMean+iNbSigmaRight*dStd
	and shift the grid so that 0 is in the grid

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEGridUnif(/* Inputs */
								int iNb,
								double dMean,
								double dStd,
								double iNbSigmaLeft,
								double iNbSigmaRight,
									 
							 /* Outputs */
								double *Grid,
								int	*pIndex0)
{
/* Declaration of locals variables */
int iNumStep;
double dStep;
double dMin;

	/* Calculation of the step in X grid */
	if (iNb<2)
	{
		Grid[0]=0;
		*pIndex0 = 0;
	}
	else
	{
		dStep = ((iNbSigmaLeft+iNbSigmaRight)*dStd)/(iNb-1);

		/* Calculation of the min value of X */
		dMin = dMean - iNbSigmaLeft*dStd;

		/* Shift the grid so that X=0 is in the grid */
		if ((dMin <=0)&&(iNb>=2))
		{
			*pIndex0 =(int)(-dMin/dStep+1e-08);
			dMin = -(*pIndex0)*dStep;
		}
		else
		{
			dMin=0;
			*pIndex0 = 0;
		}
		
		/* Filling the grid in X */
		/* X = 0 is in the grid  */
		for (iNumStep=0; iNumStep<iNb; iNumStep++ )
			Grid[iNumStep]=dMin+iNumStep*dStep;
	}

}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEMakePhiGrid														  

	Compute the dynamic grid on Phi

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEMakePhiGrid(/* Inputs */
								int iNbPhi,
								double dPhiMean,
								double dPhiStd,
								double iNbSigmaPhiGrid,

								/* Parameter */
								int IsExtremePoints,
								double iNbSigmaExtremePoints,
								int IsPhiBoundExp,
								
								/* Outputs */
								double *GridPhi,
								int	   *pIndexPhiMean)
{
/* Declaration of locals variables */
double  dLogPhiStd, dLogPhiMean;
double dPhiStep, dPhiMin, dNewPhiMin;
int iNumStep;

double PhiFloor = 0.0;	
double dPhiMinExtreme, dPhiMaxExtreme;
		
	if (iNbPhi == 1)
	{
		GridPhi[0] = dPhiMean;
		*pIndexPhiMean = 0;
	}
	else
	{
		switch (IsPhiBoundExp)
		{
		/*--------------------
			PhiBoundExp Grid
		 ----------------------*/
		case LGMSV_TRUE:
			{
				if (dPhiMean ==0)
				{
					/* Filling the grid in Phi with constant value */
					for (iNumStep=0; iNumStep<iNbPhi; iNumStep++ )
						GridPhi[iNumStep] = dPhiMean;
				}
				else		
				{
					/* Calculation of mean and Std of log Phi*/
					dLogPhiStd = sqrt(log(1+dPhiStd*dPhiStd/dPhiMean/dPhiMean));
					dLogPhiMean = log(dPhiMean)-dLogPhiStd*dLogPhiStd/2;
					
					switch (IsExtremePoints)
					{
					case LGMSV_TRUE:
						{
							if (iNbPhi == 3)
							{
								GridPhi[0]=dPhiMean*exp(-iNbSigmaExtremePoints*dLogPhiStd);
								GridPhi[1]=dPhiMean;
								GridPhi[2]=dPhiMean*exp(iNbSigmaExtremePoints*dLogPhiStd);
								*pIndexPhiMean = 1;
							
							}
							else
							{
								/* Calculation of the step in Phi grid */
								dPhiStep = dPhiMean*(exp(iNbSigmaPhiGrid*dLogPhiStd)-
												exp(-iNbSigmaPhiGrid*dLogPhiStd))/(iNbPhi-3);

								/* Calculation of the min value of Phi */
								dPhiMin = dPhiMean*exp(-iNbSigmaPhiGrid*dLogPhiStd);

								/* Shift the grid so that dPhiMean is in the grid */
								*pIndexPhiMean =(int) floor((dPhiMean-dPhiMin)/dPhiStep);
								dNewPhiMin = dPhiMean-(*pIndexPhiMean)*dPhiStep;

								*pIndexPhiMean += 1;
							
								/* Filling the grid in Phi */
								for (iNumStep=1; iNumStep<iNbPhi-1; iNumStep++ )
									GridPhi[iNumStep]=dNewPhiMin+(iNumStep-1)*dPhiStep;
								GridPhi[0] = dPhiMean*exp(-iNbSigmaExtremePoints*dLogPhiStd);
								GridPhi[iNbPhi-1] = dPhiMean*exp(iNbSigmaExtremePoints*dLogPhiStd);
							}
							break;
						}
					
					default :
						{
							/* Calculation of the step in Phi grid */
							dPhiStep = dPhiMean*(exp(iNbSigmaPhiGrid*dLogPhiStd)-
												exp(-iNbSigmaPhiGrid*dLogPhiStd))/(iNbPhi-1);

							/* Calculation of the min value of LogPhi */
							dPhiMin = dPhiMean*exp(-iNbSigmaPhiGrid*dLogPhiStd);

							/* Shift the grid so that log(dPhiMean) is in the grid */
							*pIndexPhiMean =(int) floor((dPhiMean-dPhiMin)/dPhiStep);
							dNewPhiMin = dPhiMean-(*pIndexPhiMean)*dPhiStep;
							
							/* Filling the grid in Phi */
							for (iNumStep=0; iNumStep<iNbPhi; iNumStep++ )
								GridPhi[iNumStep]=dNewPhiMin+(iNumStep-1)*dPhiStep;

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
						if (iNbPhi == 3)
						{
							/* Phi Should be positive */
							GridPhi[0]=max(dPhiMean-iNbSigmaExtremePoints*dPhiStd,PhiFloor);
							
							GridPhi[1]=dPhiMean;
							GridPhi[2]=dPhiMean+iNbSigmaExtremePoints*dPhiStd;
							*pIndexPhiMean = 1;
						}
						else
						{
							/* Calculation of the Phi Min Extreme point */
							dPhiMinExtreme = max(dPhiMean-iNbSigmaExtremePoints*dPhiStd,PhiFloor);
							dPhiMaxExtreme = dPhiMean+iNbSigmaExtremePoints*dPhiStd;
							
							/* Calculation of the step in Phi grid */
							dPhiStep = (dPhiMean+iNbSigmaPhiGrid*dPhiStd-max(dPhiMean-iNbSigmaPhiGrid*dPhiStd,PhiFloor))
								/(iNbPhi-3);

							/* Calculation of the min value of Phi */
							dPhiMin = max(dPhiMean - iNbSigmaPhiGrid*dPhiStd,PhiFloor);

							/* Shift the grid so that dPhiMean is in the grid */
							*pIndexPhiMean =(int) floor((dPhiMean-dPhiMin)/dPhiStep);
							dNewPhiMin = dPhiMean-(*pIndexPhiMean)*dPhiStep;

							
							/* Filling the grid in Phi */
							if (dNewPhiMin != dPhiMinExtreme)
							{
								*pIndexPhiMean += 1;
								
								for (iNumStep=1; iNumStep<iNbPhi-1; iNumStep++ )
									GridPhi[iNumStep]=dNewPhiMin+(iNumStep-1)*dPhiStep;
							
								GridPhi[0] = dPhiMinExtreme;
								GridPhi[iNbPhi-1] = dPhiMaxExtreme;
							}
							else
							{	
								for (iNumStep=0; iNumStep<iNbPhi-1; iNumStep++ )
									GridPhi[iNumStep]=dNewPhiMin+iNumStep*dPhiStep;
							
								GridPhi[iNbPhi-1] = dPhiMaxExtreme;
							}

						}
						break;
					}
				default:
					{
						/* Calculation of the step in Phi grid */
						dPhiStep = (dPhiMean+iNbSigmaPhiGrid*dPhiStd-max(dPhiMean-iNbSigmaPhiGrid*dPhiStd,PhiFloor))
							/(iNbPhi-1);

						/* Calculation of the min value of Phi */
						dPhiMin = max(dPhiMean - iNbSigmaPhiGrid*dPhiStd,PhiFloor);

						/* Shift the grid so that dPhiMean is in the grid */
						*pIndexPhiMean =(int) floor((dPhiMean-dPhiMin)/dPhiStep);
						dPhiMin = dPhiMean-(*pIndexPhiMean)*dPhiStep;
							
						/* Filling the grid in Phi */
						for (iNumStep=0; iNumStep<iNbPhi; iNumStep++ )
							GridPhi[iNumStep]=dPhiMin+iNumStep*dPhiStep;
						break;
					}
				}
				break;
			}
		}	
	}
}


/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEInterpolate													  

	Interpolation

   ------------------------------------------------------------------------------------------------ */
static double LGMSVPDEInterpolate(/* Inputs */
								   double x0,
								   int iNbfx,
								   double *x,
								   double ****fx,
								   int iNumX,
								   int iNumEps,
								   int iNumProduct,
								   int InterpolationType,
								   int IsExtrapolationFlat,
								   int iIndex)
{
/* declaration of locals variables */
double fx0 ;
int iIndex1, iIndex2, iIndex3;

	/* Extrapolation flat */
	if (((x0<x[0]) || (x0>x[iNbfx-1])) && (IsExtrapolationFlat == 1))
	{
		if (x0>x[iNbfx-1])
			fx0 = fx[iNbfx-1][iNumX][iNumEps][iNumProduct];
		else
			fx0 = fx[0][iNumX][iNumEps][iNumProduct];
	}
	else
	{
		/* Interpolation Treatment */
		switch (InterpolationType)
		{
			case LGMSV_LINEAR:
			{
				/* linear interpolation */
				if (iNbfx == 1) 
					fx0 = fx[0][iNumX][iNumEps][iNumProduct];
				else
				{
					if (x[0]==x[1])	
						fx0 = fx[0][iNumX][iNumEps][iNumProduct];
					else
					{
						if (iIndex ==0)
							/* Linear Extrapolation */
							fx0 = (x0-x[iIndex])/(x[iIndex+1]-x[iIndex])*(fx[iIndex+1][iNumX][iNumEps][iNumProduct]-
								fx[iIndex][iNumX][iNumEps][iNumProduct]) + fx[iIndex][iNumX][iNumEps][iNumProduct];
						else
							fx0 = (x0-x[iIndex-1])/(x[iIndex]-x[iIndex-1])*(fx[iIndex][iNumX][iNumEps][iNumProduct]-
								fx[iIndex-1][iNumX][iNumEps][iNumProduct]) + fx[iIndex-1][iNumX][iNumEps][iNumProduct];
					}
				}
				break;
			}
			
			case  LGMSV_QUADRATIC:
			{
				/* Quadratic interpolation by default */
				if (iNbfx == 1) 
					fx0 = fx[0][iNumX][iNumEps][iNumProduct];
				else
				{
					if (x[0]==x[1])
						fx0 = fx[0][iNumX][iNumEps][iNumProduct];
					else
					{
						/* Calculation of the 3 nearest points */
							
						if ((iIndex == 0) || (iIndex ==1))
						{
							iIndex1 = 0;
							iIndex2 = 1;
							iIndex3 = 2;
						}
						else
						{
							if (iIndex == iNbfx-1)
							{
								iIndex1 = iNbfx-3;
								iIndex2 = iNbfx-2;
								iIndex3 = iNbfx-1;
							}
							else
							{
								if (fabs(x0-x[iIndex-2])>fabs(x0-x[iIndex+1]))
								{
									iIndex1 = iIndex-1;
									iIndex2 = iIndex;
									iIndex3 = iIndex+1;
								}
								else
								{
									iIndex1 = iIndex-2;
									iIndex2 = iIndex-1;
									iIndex3 = iIndex;
								}
							}
						}

						/* Interpolation */
						fx0 = fx[iIndex1][iNumX][iNumEps][iNumProduct]/(x[iIndex1]-x[iIndex2])/(x[iIndex1]-x[iIndex3])*
								(x0*x0-(x[iIndex2]+x[iIndex3])*x0+x[iIndex2]*x[iIndex3])+
							fx[iIndex2][iNumX][iNumEps][iNumProduct]/(x[iIndex2]-x[iIndex3])/(x[iIndex2]-x[iIndex1])*
								(x0*x0-(x[iIndex3]+x[iIndex1])*x0+x[iIndex3]*x[iIndex1])+
							fx[iIndex3][iNumX][iNumEps][iNumProduct]/(x[iIndex3]-x[iIndex1])/(x[iIndex3]-x[iIndex2])*
								(x0*x0-(x[iIndex1]+x[iIndex2])*x0+x[iIndex1]*x[iIndex2]);
					}
				}			
				break;
			}
			
			default:
			{
				/*  */
				if (iNbfx == 1) 
					fx0 = fx[0][iNumX][iNumEps][iNumProduct];
				else
				{
					if (x[0]==x[1])
						fx0 = fx[0][iNumX][iNumEps][iNumProduct];
					else
					{
						if ((iNbfx>5) && (x0>x[2]) && (x0<x[iNbfx-3]))
						{
							/* Interpolation quadratic */
							
							/* Calculation of the 3 nearest points */
								
							if ((iIndex == 0) || (iIndex ==1))
							{
								iIndex1 = 0;
								iIndex2 = 1;
								iIndex3 = 2;
							}
							else
							{
								if (iIndex == iNbfx-1)
								{
									iIndex1 = iNbfx-3;
									iIndex2 = iNbfx-2;
									iIndex3 = iNbfx-1;
								}
								else
								{
									if (fabs(x0-x[iIndex-2])>fabs(x0-x[iIndex+1]))
									{
										iIndex1 = iIndex-1;
										iIndex2 = iIndex;
										iIndex3 = iIndex+1;
									}
									else
									{
										iIndex1 = iIndex-2;
										iIndex2 = iIndex-1;
										iIndex3 = iIndex;
									}
								}
							}

							/* Interpolation */
							fx0 = fx[iIndex1][iNumX][iNumEps][iNumProduct]/(x[iIndex1]-x[iIndex2])/(x[iIndex1]-x[iIndex3])*
									(x0*x0-(x[iIndex2]+x[iIndex3])*x0+x[iIndex2]*x[iIndex3])+
								fx[iIndex2][iNumX][iNumEps][iNumProduct]/(x[iIndex2]-x[iIndex3])/(x[iIndex2]-x[iIndex1])*
									(x0*x0-(x[iIndex3]+x[iIndex1])*x0+x[iIndex3]*x[iIndex1])+
								fx[iIndex3][iNumX][iNumEps][iNumProduct]/(x[iIndex3]-x[iIndex1])/(x[iIndex3]-x[iIndex2])*
									(x0*x0-(x[iIndex1]+x[iIndex2])*x0+x[iIndex1]*x[iIndex2]);
						}
						else
						{
							/* linear interpolation */
							if (iIndex ==0)
								/* Linear Extrapolation */
								fx0 = (x0-x[iIndex])/(x[iIndex+1]-x[iIndex])*(fx[iIndex+1][iNumX][iNumEps][iNumProduct]-
									fx[iIndex][iNumX][iNumEps][iNumProduct]) + fx[iIndex][iNumX][iNumEps][iNumProduct];
							else
								fx0 = (x0-x[iIndex-1])/(x[iIndex]-x[iIndex-1])*(fx[iIndex][iNumX][iNumEps][iNumProduct]-
									fx[iIndex-1][iNumX][iNumEps][iNumProduct]) + fx[iIndex-1][iNumX][iNumEps][iNumProduct];

						}
					}
				}			
				break;
			}
		}

	}
	if (fx0<0)
	{
		fx0 = fx0;
	}

	return fx0;
}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDECalculationPayoffMatrixfromPhit												  

	

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDECalculationValueMatrixfromPhit(
											 /* Inputs */
											 int	iIndexMinX,
											 int	iIndexMaxX,
											 int	iIndexMinZ,
											 int	iIndexMaxZ,

											 int	iNbPhi,
											 int	iProductNb,
											 double	*GridftTstar, 
											 double	*GridZCenter,
											 double	*GridPhi_t1,

											 double dVFloor,

											 /* Martingale Mesure QTstar information */
											 double	dTstar,		/* In time */

											 /* 4 dimensions : Phit,Xt,Zt,Product	*/
											 double ****ValueTensor,

											 double dZtMean,

											 double dPhit,
											 double dSigma,
											 double dSigt,
											 double dLambdaX,
											 double dAlpha,
											 double	dRho,

											 double t, /*time(t) */
											 double t1, /*time(t+1) */

											 int TypeInterpolationPhi,
											 int IsExtrapolationFlat,

											 double **InterpQuadraCoef,

											 /* Output */
											 /* 3 dimensions : Xt,Zt,Product	*/
											 double ***ValueGridXtZtProduct)
{
/* Declaration of locals variables */
double dDt, dVt, dPhi_t1;
int iNumZ, iNumX, iNumProduct;
int iIndexPhi;

double dDerBeta, dBetaDerBeta;
double dAlphaRho_SigtDerBeta, dSigma2Sigt2Dt, dPhit_1m2LxDt;

int iIndex1, iIndex2, iIndex3;

	
	/* variation of time between t and t+1 */
	dDt = t1-t;

	/* Calculation of beta(t,Tstar)*D/DTstar ( beta(t,Tstar)) and D/DTstar ( beta(t,Tstar)) */
	dDerBeta =  exp(-dLambdaX*(dTstar-t));
	dBetaDerBeta = dDerBeta * (1-dDerBeta)/dLambdaX;

	/* Calculation of constant */
	dAlphaRho_SigtDerBeta = dAlpha*dRho/dSigt/dDerBeta;
	dSigma2Sigt2Dt = dSigma*dSigt*dSigma*dSigt*dDt;
	dPhit_1m2LxDt=dPhit*(1-2*dLambdaX*dDt);

	
	for (iNumX=iIndexMinX;iNumX<=iIndexMaxX;iNumX++)	
	{
		iIndexPhi=-1;
		for (iNumZ=iIndexMinZ;iNumZ<=iIndexMaxZ;iNumZ++)
		{
			/* Calculation of Xt from the grid of f(t,Tstar) */
			/*dXt = (GridftTstar[iNumX] - dBetaDerBeta*dPhit)/dDerBeta;*/
			
			/* Calculation of the Vt from  the couple (ft,Zt) */
			/* dVt = max((GridZCenter[iNumZ]+dZtMean)+dAlpha*dRho/dSigt/dDerBeta*GridftTstar[iNumX],dVFloor); */
			dVt = max((GridZCenter[iNumZ]+dZtMean)+dAlphaRho_SigtDerBeta*GridftTstar[iNumX],dVFloor);
			
			/* Calculation of the Phi(t+1) */
			/* dPhit+(dSigma2Sigt2*dVt-2*dLambdaX*dPhit)*dDt; */
			dPhi_t1=dPhit_1m2LxDt+dSigma2Sigt2Dt*dVt;

			/* Find the Index such that GridPhi_t1[index]<dPhi_t1 and GridPhi_t1[index+1]>= dPhi_t1 */
			while ((iIndexPhi <iNbPhi-1) && (GridPhi_t1[iIndexPhi+1]<dPhi_t1))
				iIndexPhi++;

			/* Interpolate value */
			if (iIndexPhi==-1)
				/* Extrapolation flat */
				for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
					ValueGridXtZtProduct[iNumX][iNumZ][iNumProduct] = 
						ValueTensor[0][iNumX][iNumZ][iNumProduct];
			else
			if (iIndexPhi==iNbPhi-1)
				/* Extrapolation flat */
				for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
					ValueGridXtZtProduct[iNumX][iNumZ][iNumProduct] = 
						ValueTensor[iNbPhi-1][iNumX][iNumZ][iNumProduct];
			else
			if ((iIndexPhi>1)&&(iIndexPhi<iNbPhi-3)&&(TypeInterpolationPhi==LGMSV_MIXQUADRALINEAR))
			{
				/* Quadratic Interpolation  */
				if (dPhi_t1-GridPhi_t1[iIndexPhi]>0.5*(GridPhi_t1[iIndexPhi+1]-GridPhi_t1[iIndexPhi]))
				{
					/*
					iIndex1 = iIndexPhi;
					iIndex2 = iIndexPhi+1;
					iIndex3 = iIndexPhi+2;

					for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
					{
						ValueGridXtZtProduct[iNumX][iNumZ][iNumProduct] = 
							ValueTensor[iIndex1][iNumX][iNumZ][iNumProduct]/(GridPhi_t1[iIndex1]-GridPhi_t1[iIndex2])
							/(GridPhi_t1[iIndex1]-GridPhi_t1[iIndex3])*(dPhi_t1*dPhi_t1-(GridPhi_t1[iIndex2]+GridPhi_t1[iIndex3])*dPhi_t1
							+GridPhi_t1[iIndex2]*GridPhi_t1[iIndex3])+
							ValueTensor[iIndex2][iNumX][iNumZ][iNumProduct]/(GridPhi_t1[iIndex2]-GridPhi_t1[iIndex3])
							/(GridPhi_t1[iIndex2]-GridPhi_t1[iIndex1])*(dPhi_t1*dPhi_t1-(GridPhi_t1[iIndex3]+GridPhi_t1[iIndex1])*dPhi_t1
							+GridPhi_t1[iIndex3]*GridPhi_t1[iIndex1])+
							ValueTensor[iIndex3][iNumX][iNumZ][iNumProduct]/(GridPhi_t1[iIndex3]-GridPhi_t1[iIndex1])
							/(GridPhi_t1[iIndex3]-GridPhi_t1[iIndex2])*(dPhi_t1*dPhi_t1-(GridPhi_t1[iIndex1]+GridPhi_t1[iIndex2])*dPhi_t1
							+GridPhi_t1[iIndex1]*GridPhi_t1[iIndex2]);					
					}
					*/
					iIndex1 = iIndexPhi;
					iIndex2 = iIndexPhi+1;
					iIndex3 = iIndexPhi+2;
					for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
					{
						ValueGridXtZtProduct[iNumX][iNumZ][iNumProduct] = 
							ValueTensor[iIndex1][iNumX][iNumZ][iNumProduct]*InterpQuadraCoef[iIndex1][0]*
							(dPhi_t1*dPhi_t1-InterpQuadraCoef[iIndex1][1]*dPhi_t1+InterpQuadraCoef[iIndex1][2])+
							ValueTensor[iIndex2][iNumX][iNumZ][iNumProduct]*InterpQuadraCoef[iIndex1][3]*
							(dPhi_t1*dPhi_t1-InterpQuadraCoef[iIndex1][4]*dPhi_t1+InterpQuadraCoef[iIndex1][5])+
							ValueTensor[iIndex3][iNumX][iNumZ][iNumProduct]*InterpQuadraCoef[iIndex1][6]*
							(dPhi_t1*dPhi_t1-InterpQuadraCoef[iIndex1][7]*dPhi_t1+InterpQuadraCoef[iIndex1][8]);					
					}
				}
				else
				{
					iIndex1 = iIndexPhi-1;
					iIndex2 = iIndexPhi;
					iIndex3 = iIndexPhi+1;

					for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
					{
						ValueGridXtZtProduct[iNumX][iNumZ][iNumProduct] = 
							ValueTensor[iIndex1][iNumX][iNumZ][iNumProduct]*InterpQuadraCoef[iIndex1][0]*
							(dPhi_t1*dPhi_t1-InterpQuadraCoef[iIndex1][1]*dPhi_t1+InterpQuadraCoef[iIndex1][2])+
							ValueTensor[iIndex2][iNumX][iNumZ][iNumProduct]*InterpQuadraCoef[iIndex1][3]*
							(dPhi_t1*dPhi_t1-InterpQuadraCoef[iIndex1][4]*dPhi_t1+InterpQuadraCoef[iIndex1][5])+
							ValueTensor[iIndex3][iNumX][iNumZ][iNumProduct]*InterpQuadraCoef[iIndex1][6]*
							(dPhi_t1*dPhi_t1-InterpQuadraCoef[iIndex1][7]*dPhi_t1+InterpQuadraCoef[iIndex1][8]);					
					}
				}
			}
			else
			{
				/* Linear Interpolation */
				for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
				{
					ValueGridXtZtProduct[iNumX][iNumZ][iNumProduct] = (dPhi_t1-GridPhi_t1[iIndexPhi])/(GridPhi_t1[iIndexPhi+1]-GridPhi_t1[iIndexPhi])*
					 (ValueTensor[iIndexPhi+1][iNumX][iNumZ][iNumProduct]-ValueTensor[iIndexPhi][iNumX][iNumZ][iNumProduct])
					 + ValueTensor[iIndexPhi][iNumX][iNumZ][iNumProduct];
				}

			}

		}
	}


}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDECalculationRDriftVar												  

	

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDECalculationRDriftVar(/* Inputs */
										int iIndexMinX,
										int	iIndexMaxX,
										int	iIndexMinZ,
										int	iIndexMaxZ,
										
										double	*GridftTstar,
										double	*GridZCenter,
									
										double dSigma,
										double dSigt,
										double dSig_t1,
										double dXtMean,
										double dZtMean,
										double dPhitMean,
										double dLambdaX,
										double dLambdaEps,
										double dAlpha,
										double dRho,
										double dIfrt,
										double t, /*time(t) */
										double t1, /*time(t+1) */

										double dVFloor,

										/* Martingale Mesure QTstar information */
										double	dTstar,		/* In time */
																				
										/* Outputs */
										double **DfDrift,
										double **DZDrift,
										double **DfVar,
										double **DZVar,
										double **r)
{
/* Declaration of locals variables */
double dDt;
double dVt, dZt;
int iNumX, iNumZ;
double dDerBeta, dBetaDerBeta;

double dAlphaRho_DerBeta, dAlphaRho_SigtDerBeta;
double dDerBeta2Sigt2Dt, dAlpha2_1mRho2Dt, DZDriftCoef;

	/* variation of time between t and t+1 */
	dDt = t1-t;

	/* Calculation of beta(t,Tstar)*D/DTstar ( beta(t,Tstar)) and D/DTstar ( beta(t,Tstar)) */
	dDerBeta = exp(-dLambdaX*(dTstar-t));
	dBetaDerBeta = dDerBeta * (1-dDerBeta)/dLambdaX;

	/* Calculation of constant */
	dAlphaRho_DerBeta = dAlpha*dRho/dDerBeta;
	dAlphaRho_SigtDerBeta = dAlphaRho_DerBeta/dSigt;
	dDerBeta2Sigt2Dt = dDerBeta*dDerBeta*dSigt*dSigt*dDt;
	
	DZDriftCoef  = dAlphaRho_DerBeta*((1/dSig_t1-1/dSigt)-dLambdaX/dSigt*dDt);
	dAlpha2_1mRho2Dt = dAlpha*dAlpha*(1-dRho*dRho)*dDt;
	
	for (iNumX=iIndexMinX;iNumX<=iIndexMaxX;iNumX++)
	{
		/* Calculation of Xt from the grid of f(t,Tstar) */
		/*dXt = (GridftTstar[iNumX]- dBetaDerBeta*dPhit)/dDerBeta; */
		
		for (iNumZ=iIndexMinZ;iNumZ<=iIndexMaxZ;iNumZ++)
		{
			dZt = GridZCenter[iNumZ]+dZtMean;
						
			/* Calculation of the Vt from  the couple (ft,Zt) */
			/* dVt = max(dZt+dAlpha*dRho/dSigt/dDerBeta*GridftTstar[iNumX],dVFloor); */
			dVt = max(dZt+dAlphaRho_SigtDerBeta*GridftTstar[iNumX],dVFloor);
			
			/* Discount rate * dt */
			r[iNumX][iNumZ] = 0;/**(dXt+dIfrt)*dDt;*/
						
			/* dt*Drifts under Qt+1 and dt*Variances under Qt+1 */
			
			/* Et+1 [Xt+1 - Xt / Ft] and Vt+1 [Xt+1 - Xt / Ft]*/
			/*	DfDrift[iNumX][iNumZ] = 0;
				DfVar[iNumX][iNumZ] = dDerBeta*dDerBeta*dSigt*dSigt*dVt*dDt; */
			DfDrift[iNumX][iNumZ] = 0;
			DfVar[iNumX][iNumZ] = dDerBeta2Sigt2Dt*dVt;
			
			/* Et+1 [Zt+1 - Zt / Ft] and Vt+1 [Zt+1 - Zt / Ft]*/
			/*	DZDrift[iNumX][iNumZ] = -dLambdaEps*(dVt-1)*dDt-dAlpha*dRho/dDerBeta*GridftTstar[iNumX]*
					((1/dSig_t1-1/dSigt)-dLambdaX/dSigt*dDt);
				DZVar[iNumX][iNumZ] = dAlpha*dAlpha*dVt*(1-dRho*dRho)*dDt;			*/
			DZDrift[iNumX][iNumZ] = -dLambdaEps*(dVt-1)*dDt-DZDriftCoef*GridftTstar[iNumX];

			DZVar[iNumX][iNumZ] = dAlpha2_1mRho2Dt*dVt;
		}
	}
}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVSaveTensor
	
	  Save the tensor for analysis
																			
   ------------------------------------------------------------------------------------------------ */

void LGMSVSaveTensor( /* Inputs */
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
																				
	Xt= rt - f(0,t)															
	dXt =  [Phit-lambda_X*Xt]dt + sigma(Xt,t,omega) dWt						
	dPhit = [sigma*sigma - 2*lambda_X*Phit] dt								
																				
	sigma : time, spot dependant and stochastic of the following form		
		sigma(Xt,t,omega)  = sig(t)*sigma* Epst							

		Where																	
			The stochastic part of the sigma : Epst follows :					
				dEpst = -lambda_Eps*[Epst-1]dt + alpha * Epst * dZt						
				(mean reverting to 1, log normal vol)									
																			
   ------------------------------------------------------------------------------------------------ */

Err	 lgmSV_adi(	
					/*	Time Information  */
					int			iNbTime,					
					double		*dTime,
					double		*dDate,

					/*	Space Discretisation	*/
					int			iNbPhi,
					int			iNbX,
					int			iNbEps,
						
					/*	Model data Information	*/
					double		dLambdaX,
					double		dSigma,
					
					double		*SigTime,
					double		*Sig,
					int			iNbSigTime,
					
					double		dLambdaEps,
					double		dAlpha,
					double		dRho,

					/* Parameters */
					LGMSVParam Params,

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
										
										/* Grid data Information	*/
										int		iIndPhitMin,
										int		iIndPhitMax,
										int		iIndXtMin,
										int		iIndXtMax,
										int		iIndEpstMin,
										int		iIndEpstMax,
										
										double	*GridPhi,
										double	*GridftTstar,
						
										/* Tensor of results to be updated		*/
										/* 4 dimensions : Phit,Xt,Epst,Product	*/
										int		iNbProduct,
										double	****PayoffTensor),
					/*	Result */
					int			iNbProduct, 
					double		*dProductArrayPv)
{

/* Declaration of locals variables */
double			*Sigt = NULL;

double			*GridftTstar = NULL;
double			*GridZCenter = NULL;

double			*GridPhi_t = NULL;
double			*GridPhi_t1 = NULL;

double			*PhiMean = NULL;
double			*PhiStd = NULL;
double			*ftMean = NULL;
double			*ftStd = NULL;
double			*ZtMean = NULL;
double			*ZtStd = NULL;

double			**DfDrift = NULL;
double			**DZDrift = NULL;
double			**DfVar = NULL;
double			**DZVar = NULL;
double			**r = NULL;

double			***ValueGridftZtProduct = NULL;
double			****ValueTensor_t = NULL;
double			****ValueTensor_t1 = NULL;

Err				err = NULL;
clock_t			t1, t2;

double			dLambdaXGrid, dLambdaEpsGrid;

CNPDE_TEMP_2D_ADI	pdestr, *pde = &pdestr; 
/*CNPDE_TEMP_2D_LOD	pdestr, *pde = &pdestr;*/

int iNumTime, iNumPhi, iNumX, iNumEps, iNumProduct;
int iIndexZ0, iIndexf0, iIndexPhi;

double dTstar = Params.Tstar;
double dB0Tstar;

/* Variables for computation optimisation */
int iNbPhiGrid_t, iNbPhiGrid_t1;
int iIndexMinX, iIndexMaxX, iIndexMinZ, iIndexMaxZ;

double a1,a2,a3;
double **InterpQuadraCoef = NULL;


	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */

	/* For computational time calculation				 */
	t1 = clock();

	/* Auto size the Right ZGridCenter if applied */
	/* depend on density of X at time T			  */
	/* fct(dRho, Equivalent Alpha at T)			  */
	/* with alpha and lambdaEps on Eps !!!		  */
	/* Empirical coefficients					  */
	if (Params.IsNbSigmaXGridRightAuto == LGMSV_TRUE)
	{
		if (fabs(dLambdaEps)<1E-06)
		{
			Params.iNbSigmaXGridRight = 7.5+25*max(dRho,0)*dAlpha/2.0;
		}
		else
		{
			Params.iNbSigmaXGridRight = 7.5+25*max(dRho,0)*dAlpha/2.0*
				sqrt((1-exp(-dLambdaEps*dTime[iNbTime-1]))/(dLambdaEps*dTime[iNbTime-1]));
		}
	}

	/* grid should be at least a 3x3 in X, Eps			*/
	iNbX = max(iNbX,3);
	iNbEps = max(iNbEps,3);

	/* All Number of points in the grid must be odd		*/
	iNbX = ((int) (iNbX/2))*2 +1;
	iNbEps = ((int) (iNbEps/2))*2 +1;
	iNbPhi = ((int) (iNbPhi/2))*2 +1;
	
	/*  Memory allocations								 */
	Sigt = dvector(0, iNbTime-1);

	GridftTstar = dvector(0, iNbX-1);
	GridZCenter = dvector(0, iNbEps-1);
	
	GridPhi_t = dvector(0, iNbPhi-1);	/* Grid of Phi at t   */
	GridPhi_t1 = dvector(0, iNbPhi-1);	/* Grid of Phi at t+1 */

	PhiMean = dvector(0, iNbTime-1);
	PhiStd = dvector(0, iNbTime-1);
	
	ftMean = dvector(0, iNbTime-1);
	ftStd = dvector(0, iNbTime-1);

	ZtMean = dvector(0, iNbTime-1);
	ZtStd = dvector(0, iNbTime-1);

	DfDrift = dmatrix(0, iNbX-1, 0, iNbEps-1);
	DZDrift = dmatrix(0, iNbX-1, 0, iNbEps-1);
	DfVar = dmatrix(0, iNbX-1, 0, iNbEps-1);
	DZVar = dmatrix(0, iNbX-1, 0, iNbEps-1);
	r = dmatrix(0, iNbX-1, 0, iNbEps-1);
	
	ValueGridftZtProduct =f3tensor(0, iNbX-1,
									  0, iNbEps-1,
									  0, iNbProduct-1);
		
	ValueTensor_t =f4tensor(0, iNbPhi-1,
							0, iNbX-1,
							0, iNbEps-1,
							0, iNbProduct-1);
	
	ValueTensor_t1 =f4tensor(0, iNbPhi-1,
							 0, iNbX-1,
							 0, iNbEps-1, 
							 0, iNbProduct-1);

	InterpQuadraCoef = dmatrix(0, iNbPhi-1, 0, 8);

	/* Backward PDE */

	/* Init of the Backward PDE : Allocation of arrays    */
	num_f_pde_init_2d_adi(pde,
						iNbX,
						iNbEps,
						iNbProduct); 
	/*num_f_pde_init_2d_lod(pde,
						iNbX,
						iNbEps,
						iNbProduct);*/

	if (!pde)
	{
		err = "Memory allocation error (2) in LGMSVPDE";
		goto FREE_RETURN;
	}

	/* Gestion of allocation errors */
	if (!pde ||  !Sigt || !GridftTstar || !GridZCenter || !GridPhi_t || !GridPhi_t1 || !PhiMean || !PhiStd 
		 || !ftMean || !ftStd || !ZtMean || !ZtStd|| !DfDrift || !DZDrift || !DfVar || !DZVar || !r 
		 || !ValueGridftZtProduct || !ValueTensor_t || !ValueTensor_t1)
	{
		err = "Memory allocation error (1) in LGMSVPDE";
		goto FREE_RETURN;
	}

	/* Constant	*/
	dLambdaXGrid = dLambdaX;			/* Use of LambdaX for the grid */
	dLambdaEpsGrid = dLambdaEps;		/* Use of LambdaEps for the grid */ 

	/* Calculation of B(0,Tstar) */
	dB0Tstar  = swp_f_df (dDate[0], dDate[0]+dTstar*365.0, cYieldCurve);

	
	/* Expectations calculations */
	err = LGMSVPDEExpectations(
								/* Inputs */
								 iNbTime,			
								 dTime,
								
								 /* For X diffusion */
								 dLambdaXGrid,
								 dSigma,
								 dTstar,
								 
								 SigTime,
								 Sig,
								 iNbSigTime,

								 /* For Eps diffusion */
								 dLambdaEpsGrid,
								 dAlpha,
								 dRho,
								 
								 /* Outputs */
								 Sigt,

								 PhiMean,
								 PhiStd,
								 ftMean,
								 ftStd,
								 ZtMean,
								 ZtStd);

	if (err)
	{
		goto FREE_RETURN;
	}


	/* Verifying that ZtStd and ftStd >0 at time T */
	if ((ZtStd[iNbTime-1]==0)||(ftStd[iNbTime-1]==0))
	{
		err = "error in LGMSVPDE : Zt or ft is not stochastic";
		goto FREE_RETURN;
	}

	/* Making the grid on fcenter and ZCenter */	
								
	LGMSVPDEGridUnif(/* Inputs */
						iNbEps,
						0,
						ZtStd[iNbTime-1],

						Params.iNbSigmaZGridLeft,
						Params.iNbSigmaZGridRight,
									 
						/* Outputs */
						GridZCenter,
						&iIndexZ0); 

	LGMSVPDEGridUnif(/* Inputs */
						iNbX,
						0,
						ftStd[iNbTime-1],
						Params.iNbSigmaXGridLeft,
						Params.iNbSigmaXGridRight,
									 
						/* Outputs */
						GridftTstar,
						&iIndexf0); 


	/* Final Grid of Phi , Grid of phi at time T */
	LGMSVPDEMakePhiGrid(/* Inputs */
						iNbPhi,
						PhiMean[iNbTime-1],
						PhiStd[iNbTime-1],
						Params.iNbSigmaPhiGrid,

						/* Parameter */
						Params.IsExtremePoints,
						Params.iNbSigmaExtremePoints,
						Params.IsPhiBoundExp,
								
						/* Outputs */
						GridPhi_t1,
						&iIndexPhi);

	
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
										
						/* Grid data Information	*/
						0,
						iNbPhi-1,
						0,
						iNbX-1,
						0,
						iNbEps-1,
										
						GridPhi_t1,
						GridftTstar,
						
						/* Tensor of results to be updated		*/
						/* 4 dimensions : Phit,Xt,Epst,Product	*/
						iNbProduct,
						ValueTensor_t1);
	
	if (Params.VerifExpectation==LGMSV_TRUE)
	{
		for (iNumPhi = 0; iNumPhi<iNbPhi; iNumPhi++)
			for (iNumX = 0; iNumX<iNbX; iNumX++)
				for (iNumEps = 0; iNumEps<iNbEps; iNumEps++)
					for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
					ValueTensor_t1[iNumPhi][iNumX][iNumEps][iNumProduct]=GridZCenter[iNumEps]*GridZCenter[iNumEps];
	}

	/* Init Nb Phi at t+1 */
	iNbPhiGrid_t1 = iNbPhi;


	/* End of initialisation treatment */
	
	/* Initialisation time display */
	t2 = clock();
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Phase 2 -convolution, NbSigXGridRight: %f ", Params.iNbSigmaXGridRight);

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

if (Params.ModeExport!=LGMSV_TRUE)
{	
	/* For each time t backward */
	for (iNumTime=iNbTime-2; iNumTime>=0; iNumTime--)
	{
		/* Update Parameters for calculation optimisation */
		
		/* Number of Phi in the grid at time t */
		iNbPhiGrid_t = iNbPhi-(Params.GridPhiOptim == LGMSV_TRUE)*(int)((1-iNumTime/((double)(iNbTime-1)))
			*(iNbPhi-1) + 1.0e-08);
		iNbPhiGrid_t = ((int) (iNbPhiGrid_t/2))*2 +1;
		
		
		if (Params.GridXADIOptim == LGMSV_TRUE)
		{
			iIndexMinX = (int) ((-Params.iNbSigmaXGridLeft*ftStd[iNumTime+1]-GridftTstar[0])
				/(GridftTstar[1]-GridftTstar[0]) + 1.0e-08);
			iIndexMaxX = max((int) ((Params.iNbSigmaXGridRight*ftStd[iNumTime+1]-GridftTstar[0])
				/(GridftTstar[1]-GridftTstar[0]) + 1.0e-08),min(iIndexMinX+2,iNbX-1));
		}
		else
		{
			iIndexMinX = 0;
			iIndexMaxX = iNbX-1;
		}

		if (Params.GridZADIOptim == LGMSV_TRUE)
		{
			iIndexMinZ = (int) ((-Params.iNbSigmaZGridLeft*ZtStd[iNumTime+1]-GridZCenter[0])
				/(GridZCenter[1]-GridZCenter[0]) + 1.0e-08);
			iIndexMaxZ = max((int) ((Params.iNbSigmaZGridRight*ZtStd[iNumTime+1]-GridZCenter[0])
				/(GridZCenter[1]-GridZCenter[0]) + 1.0e-08),min(iIndexMinZ+2,iNbEps-1));
		}
		else
		{	
			iIndexMinZ = 0;
			iIndexMaxZ = iNbEps-1;
		}

		if (Params.SaveTensor == LGMSV_TRUE)
		{
			/* Save the Tensor */
			LGMSVSaveTensor( /* Inputs */
					iNbPhi,
					iNbX,
					iNbEps,
					 
					GridPhi_t1,
					GridftTstar,
					GridZCenter,
					ValueTensor_t1);
		}
		
		/* Pre Calculation for code optimisation */
		for (iNumPhi=1; iNumPhi<iNbPhiGrid_t1-3; iNumPhi++)
		{
			a1=GridPhi_t1[iNumPhi];
			a2=GridPhi_t1[iNumPhi+1];
			a3=GridPhi_t1[iNumPhi+2];

			InterpQuadraCoef[iNumPhi][0]=1/(a1-a2)/(a1-a3);
			InterpQuadraCoef[iNumPhi][1]=a2+a3;
			InterpQuadraCoef[iNumPhi][2]=a2*a3;

			InterpQuadraCoef[iNumPhi][3]=1/(a2-a3)/(a2-a1);
			InterpQuadraCoef[iNumPhi][4]=a3+a1;
			InterpQuadraCoef[iNumPhi][5]=a3*a1;

			InterpQuadraCoef[iNumPhi][6]=1/(a3-a1)/(a3-a2);
			InterpQuadraCoef[iNumPhi][7]=a1+a2;
			InterpQuadraCoef[iNumPhi][8]=a1*a2;		
		}

		/* Make the grid in Phit at time t */
		LGMSVPDEMakePhiGrid(/* Inputs */
						iNbPhiGrid_t,
						PhiMean[iNumTime],
						PhiStd[iNumTime],
						Params.iNbSigmaPhiGrid,

						/* Parameter */
						Params.IsExtremePoints,
						Params.iNbSigmaExtremePoints,
						Params.IsPhiBoundExp,
								
						/* Outputs */
						GridPhi_t,
						&iIndexPhi);

		/* Calculation of the drifts and variances under Qt+1 */
		LGMSVPDECalculationRDriftVar(/* Inputs */
									iIndexMinX,
									iIndexMaxX,
									iIndexMinZ,
									iIndexMaxZ,

									GridftTstar, 
									GridZCenter,

									dSigma,
									Sigt[iNumTime],
									Sigt[iNumTime+1],
									ftMean[iNumTime],
									ZtMean[iNumTime],
									PhiMean[iNumTime],
									dLambdaX,
									dLambdaEps,
									dAlpha,
									dRho,
									Ifr[iNumTime],
									dTime[iNumTime], /*time(t) */
									dTime[iNumTime+1], /*time(t+1) */

									Params.VFloor,

									dTstar,										
									
									/* Outputs */
									DfDrift,
									DZDrift,
									DfVar,
									DZVar,
									r);

		/* For each Phit in the grid at time t*/
		for (iNumPhi=0;iNumPhi<iNbPhiGrid_t;iNumPhi++)
		{		
				
			/* Calculation of the Value at t+1 on the grid (XCenter,ZCenter)			 */
			/* at the deduced forwards Phi t+1 from Phit								 */
			/* Interpolation From the value in the grid (XCenter,ZCenter,Phi) at t+1	 */
			LGMSVPDECalculationValueMatrixfromPhit(
											 /* Inputs */
											 iIndexMinX,
											 iIndexMaxX,
											 iIndexMinZ,
											 iIndexMaxZ,

											 iNbPhiGrid_t1,
											 iNbProduct,
											 GridftTstar, 
											 GridZCenter,
											 GridPhi_t1,

											 Params.VFloor,

											 dTstar,

											 ValueTensor_t1,

											 ZtMean[iNumTime],
											 
											 GridPhi_t[iNumPhi],
											 dSigma,
											 Sigt[iNumTime],
											 dLambdaX,
											 dAlpha,
											 dRho,

											 dTime[iNumTime], /*time(t) */
											 dTime[iNumTime+1], /*time(t+1) */

											 Params.TypeInterpolationPhi,
											 Params.IsExtrapolationFlat,

											 InterpQuadraCoef,

											 /* Output */
											 ValueGridftZtProduct);

			/* One Step Backward 2D using Crank-Nicolson on grid (Xcenter,Zcenter) */	
			/* From t+1 to t											           */

			
			/* Launch The Crank-Nicolson					*/
			/* Calculation of the payoff tensor at time t	*/
			/* for the correspondant phi(t)					*/
			num_f_pde_one_step_backward_2f_adi(	pde,
												iNbX,
												GridftTstar,
												iNbEps,
												GridZCenter,
												0,
												iNbProduct-1,
												ValueGridftZtProduct,
												DfDrift,
												DZDrift,
												DfVar,
												DZVar,
												r,
												ValueTensor_t[iNumPhi],
												iIndexMinX,
												iIndexMaxX,
												iIndexMinZ,
												iIndexMaxZ); 

			/*num_f_pde_one_step_backward_2f_lod(	pde,
												iNbX,
												GridXCenter,
												iNbEps,
												GridZCenter,
												iNbProduct,
												ValueGridXtZtProduct,
												DXDrift,
												DZDrift,
												DXVar,
												DZVar,
												r,
												0.5,
												ValueTensor_t[iNumPhi],
												0,
												iNbX-1,
												0,
												iNbEps-1);*/

		}

		/* End for each Phit */

		/*	Case of evaluation events */
		if (EvalEvent[iNumTime])
		{
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
												
								/* Grid data Information	*/
								0,
								iNbPhiGrid_t-1,
								iIndexMinX,
								iIndexMaxX,
								iIndexMinZ,
								iIndexMaxZ,
												
								GridPhi_t,
								GridftTstar,
								
								/* Tensor of results to be updated		*/
								/* 4 dimensions : Phit,Xt,Epst,Product	*/
								iNbProduct,
								ValueTensor_t);
			if (err)
			{
					goto FREE_RETURN;
			}
		}

		/* Swaping, t+1 becomes t */
		LGMSVPDESwap(double ****,ValueTensor_t,ValueTensor_t1);
		LGMSVPDESwap(double *,GridPhi_t,GridPhi_t1);
		LGMSVPDESwap(int ,iNbPhiGrid_t,iNbPhiGrid_t1);

	}
	/* End for each t */

	/* Multiplication of the results by B(0,Tstar) to obtain the PV */
	for (iNumProduct = 0; iNumProduct<iNbProduct; iNumProduct++)
			ValueTensor_t1[iIndexPhi][iIndexf0][iIndexZ0][iNumProduct] *= dB0Tstar;

	/* Copy of the result */
	/* The grids are made such as 
		center on Phi = 0, X=0,  at t=0 */
	memcpy(dProductArrayPv,ValueTensor_t1[iIndexPhi][iIndexf0][iIndexZ0],iNbProduct*sizeof(double));

	/* Convolution time display */
	t1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (t1 - t2) / CLOCKS_PER_SEC);

}
else
{
	dProductArrayPv[0]=	PhiStd[iNbTime-1]; /*dZTStd; /*PhiStd[iNbTime-1];*/
	dProductArrayPv[1]=PhiMean[iNbTime-1]; /*dXTStd; */
}


FREE_RETURN:

	/* free all the arrays */
	if (pde) num_f_pde_free_2d_adi(pde, iNbX, iNbEps, iNbProduct); 
	/*if (pde) num_f_pde_free_2d_lod(pde, iNbX, iNbEps, iNbProduct); */
	
	if (Sigt)
		free_dvector(Sigt, 0, iNbTime-1);

	if (GridftTstar)
		free_dvector(GridftTstar, 0, iNbX-1);
	if (GridZCenter)
		free_dvector(GridZCenter, 0, iNbEps-1);

	if (GridPhi_t)
		free_dvector(GridPhi_t, 0, iNbPhi-1);
	if (GridPhi_t1)
		free_dvector(GridPhi_t1, 0, iNbPhi-1);

	if (PhiMean)
		free_dvector(PhiMean, 0, iNbTime-1);
	if (PhiStd)
		free_dvector(PhiStd, 0, iNbTime-1);

	if (ftMean)
		free_dvector(ftMean, 0, iNbTime-1);
	if (ftStd)
		free_dvector(ftStd, 0, iNbTime-1);
	
	if (ZtMean)
		free_dvector(ZtMean, 0, iNbTime-1);
	if (ZtStd)
		free_dvector(ZtStd, 0, iNbTime-1);

	if (DfDrift)
		free_dmatrix(DfDrift, 0, iNbX-1, 0, iNbEps-1);
	if (DZDrift)
		free_dmatrix(DZDrift, 0, iNbX-1, 0, iNbEps-1);
	if (DfVar)
		free_dmatrix(DfVar, 0, iNbX-1, 0, iNbEps-1);
	if (DZVar)
		free_dmatrix(DZVar, 0, iNbX-1, 0, iNbEps-1);
	if (r)
		free_dmatrix(r, 0, iNbX-1, 0, iNbEps-1);
	
	if (ValueGridftZtProduct)
		free_f3tensor(ValueGridftZtProduct, 0, iNbX-1, 0, iNbEps-1, 0, iNbProduct-1);
		
	if (ValueTensor_t)
		free_f4tensor(ValueTensor_t, 0, iNbPhi-1, 0, iNbX-1, 0, iNbEps-1, 0, iNbProduct-1);
	
	if (ValueTensor_t1)
		free_f4tensor(ValueTensor_t1,0, iNbPhi-1, 0, iNbX-1, 0, iNbEps-1,  0, iNbProduct-1);
			
	
	if (InterpQuadraCoef)
		free_dmatrix(InterpQuadraCoef, 0, iNbPhi-1, 0, 8);


	/* Return the error message */
	return err;

}


/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/

/************************************************************
 *						MODEL UNDER QBeta					*
 ************************************************************/

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEExpectations														  

	Calculation of expectations for the definition of the grids 
		in Xt (constant grid center in XTMean =E(XT / X0) and +- iNbSigmaXGrid XTStd)
		in Epst (constant grid center in EpsTMean =E(EpsT / Eps0) and +- iNbSigmaEpsGrid EpsStd)
		in Phit (variable grid center at each t on PhitMean(t) 

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEExpectations_QBeta(
					/* Inputs */
					 int	iNbTime,			
					 double	*dTime,
					
					 /* For X diffusion */
					 double dLambdaX,
					 double dSigma,
					 
					 double *SigTime,
					 double	*Sig,
					 int	iNbSigTime,

					 /* For Eps diffusion */
					 double dLambdaEps,
					 double dAlpha,
					 double	dRho,
					 
					 /* Outputs */
					 double *Sigt,			/* Value of sig(t) for all t in dtime */
					 
					 double *PhitMean,
					 double *PhitStd,
					 
					 double *XtMean,
					 double *XtStd,

					 double *ZtMean,
					 double *ZtStd)
{	
/* Declaration of locals variables */
int iNumTime, iNumSigTime;
double dDt;

double lambda[100];

LGMSVSolFunc FuncPhitMean, FuncPhit2Mean, FuncVtPhitMean, FuncVt2Mean;
LGMSVSolFunc FuncXtMean, FuncXt2Mean, FuncXtPhitMean, FuncXtVtMean;

double dPhitMean, dPhit2Mean, dVtPhitMean, dVt2Mean;
double dXtMean, dXt2Mean, dXtPhitMean, dXtVtMean;

double dg2;


	/* Initialisation of the LGMSVSolFunc */

	memset(lambda, 0, 100 * sizeof(double));
	
	/* Definition of the function E(Phit) */
	FuncPhitMean.dLambda = 2*dLambdaX;
	FuncPhitMean.bIsft1 = 1;
	FuncPhitMean.bIsgt1 = 1;
	FuncPhitMean.b = 0;
	FuncPhitMean.bIsht1 = 1;
	FuncPhitMean.c =0;

	/* Definition of the function E(Vt^2) */
	FuncVt2Mean.dLambda = 2*dLambdaEps;
	FuncVt2Mean.bIsft1 = 1;
	FuncVt2Mean.a = dAlpha*dAlpha+2*dLambdaEps;
	FuncVt2Mean.bIsgt1 = 1;
	FuncVt2Mean.b = 0;
	FuncVt2Mean.bIsht1 = 1;
	FuncVt2Mean.c =0;


	/* Definition of the function E(VtPhit) */
	FuncVtPhitMean.dLambda = 2*dLambdaX+dLambdaEps;
	FuncVtPhitMean.bIsft1 = 0;
	FuncVtPhitMean.pft = &FuncVt2Mean;
	FuncVtPhitMean.bIsgt1 = 0;
	FuncVtPhitMean.b = dLambdaEps;
	FuncVtPhitMean.pgt = &FuncPhitMean;
	FuncVtPhitMean.bIsht1 = 1;
	FuncVtPhitMean.c =0;


	/* Definition of the function E(Phit^2) */
	FuncPhit2Mean.dLambda = 4*dLambdaX;
	FuncPhit2Mean.bIsft1 = 0;
	FuncPhit2Mean.pft = &FuncVtPhitMean;
	FuncPhit2Mean.bIsgt1 = 1;
	FuncPhit2Mean.b = 0;
	FuncPhit2Mean.bIsht1 = 1;
	FuncPhit2Mean.c =0;

	/* Definition of the function E(Xt) */
	FuncXtMean.dLambda = dLambdaX;
	FuncXtMean.bIsft1 = 0;
	FuncXtMean.a = 1;
	FuncXtMean.pft = &FuncPhitMean;
	FuncXtMean.bIsgt1 = 1;
	FuncXtMean.b = 0;
	FuncXtMean.bIsht1 = 1;
	FuncXtMean.c =0;

	/* Definition of the function E(Xt^2) */
	FuncXt2Mean.dLambda = 2*dLambdaX;
	FuncXt2Mean.bIsft1 = 0;
	FuncXt2Mean.a = 2;
	FuncXt2Mean.pft = &FuncXtPhitMean;
	FuncXt2Mean.bIsgt1 = 1;
	FuncXt2Mean.bIsht1 = 1;
	FuncXt2Mean.c =0;

	/* Definition of the function E(Xt*Phit) */
	FuncXtPhitMean.dLambda = 3*dLambdaX;
	FuncXtPhitMean.bIsft1 = 0;
	FuncXtPhitMean.a = 1;
	FuncXtPhitMean.pft = &FuncPhit2Mean;
	FuncXtPhitMean.bIsgt1 = 0;
	FuncXtPhitMean.pgt = &FuncXtVtMean;
	FuncXtPhitMean.bIsht1 = 1;
	FuncXtPhitMean.c =0;

	/* Definition of the function E(Xt*Vt) */
	FuncXtVtMean.dLambda = dLambdaX+dLambdaEps;
	FuncXtVtMean.bIsft1 = 0;
	FuncXtVtMean.a = 1;
	FuncXtVtMean.pft = &FuncVtPhitMean;
	FuncXtVtMean.bIsgt1 = 0;
	FuncXtVtMean.b = dLambdaEps;
	FuncXtVtMean.pgt = &FuncXtMean;
	FuncXtVtMean.bIsht1 = 1;

	/* Filling the Sigt array 
	   piecewise interpolation from the Sig array on SigTime */
	
	iNumSigTime = 0;
	for (iNumTime = 0; iNumTime<iNbTime; iNumTime++)
	{
		if ((dTime[iNumTime]>SigTime[iNumSigTime]) && (iNumSigTime<iNbSigTime-1))
			iNumSigTime++;
		
		Sigt[iNumTime] = Sig[iNumSigTime];	
	}

	
	/* Init Value at t=0 */
	dPhitMean = 0;
	dVt2Mean = 1;
	dVtPhitMean = 0;
	dPhit2Mean = 0;

	dXtMean = 0;
	dXt2Mean = 0;
	dXtPhitMean = 0;
	dXtVtMean = 0;
	
	PhitMean[0] = dPhitMean;
	PhitStd[0] = sqrt(dPhit2Mean-dPhitMean*dPhitMean);

	ZtMean[0] = 1;
	ZtStd[0] = 0;

	XtMean[0] = 0;
	XtStd[0] = 0;
	
	/* Calculation of all the mean and variances */
	for (iNumTime = 1;iNumTime<iNbTime;iNumTime++)
	{
		/* Delta time between t and t-1 */
		dDt = dTime[iNumTime]-dTime[iNumTime-1];
		
		/* Calculation of g(t-1)^2 */
		dg2 = Sigt[iNumTime-1]*Sigt[iNumTime-1];

		/* Update value at t-1 on all func */
		FuncPhitMean.dXt1=dPhitMean;
		FuncPhitMean.a=dg2;

		FuncVt2Mean.dXt1 = dVt2Mean;

		FuncVtPhitMean.dXt1=dVtPhitMean;
		FuncVtPhitMean.a=dg2;

		FuncPhit2Mean.dXt1=dPhit2Mean;
		FuncPhit2Mean.a=2*dg2;

		FuncXtMean.dXt1=dXtMean;

		FuncXt2Mean.dXt1=dXt2Mean;
		FuncXt2Mean.b=dg2;

		FuncXtPhitMean.dXt1=dXtPhitMean;
		FuncXtPhitMean.b=dg2;

		FuncXtVtMean.dXt1=dXtVtMean;
		FuncXtVtMean.c =dAlpha*dRho*Sigt[iNumTime-1];

		/* Calculation at time t */

		/* E(Phit) at t knowing t-1 */
		dPhitMean = LGMSVFuncValue(FuncPhitMean,dDt,lambda,0);

		/* E(Vt^2) at t knowing t-1 */
		dVt2Mean = LGMSVFuncValue(FuncVt2Mean,dDt,lambda,0);

		/* E(VtPhit) at t knowing t-1 */
		dVtPhitMean = LGMSVFuncValue(FuncVtPhitMean,dDt,lambda,0);

		/* E(Phit^2) at t knowing t-1 */
		dPhit2Mean = LGMSVFuncValue(FuncPhit2Mean,dDt,lambda,0);

		/* E(Xt) at t knowing t-1 */
		dXtMean = LGMSVFuncValue(FuncXtMean,dDt,lambda,0);

		/* E(Xt^2) at t knowing t-1 */
		dXt2Mean = LGMSVFuncValue(FuncXt2Mean,dDt,lambda,0);

		/* E(Xt*Phit) at t knowing t-1 */
		dXtPhitMean = LGMSVFuncValue(FuncXtPhitMean,dDt,lambda,0);

		/* E(Xt*Vt) at t knowing t-1 */
		dXtVtMean = LGMSVFuncValue(FuncXtVtMean,dDt,lambda,0);


		/* Mean and std of Phi at t knowing t-1 */
		PhitMean[iNumTime] = dPhitMean;
		PhitStd[iNumTime] = sqrt(dPhit2Mean-dPhitMean*dPhitMean);

		/* Mean and std of Xt at t knowing t-1 */
		XtMean[iNumTime] = dXtMean;
		XtStd[iNumTime] = sqrt(dXt2Mean-dXtMean*dXtMean);

		/* Mean and std of Zt at t knowing t-1 */
		ZtMean[iNumTime] = 1-dAlpha*dRho/Sigt[iNumTime]*XtMean[iNumTime];
		ZtStd[iNumTime] = sqrt(dVt2Mean + pow(dAlpha*dRho/Sigt[iNumTime],2)*dXt2Mean
		-2*dAlpha*dRho/Sigt[iNumTime]*dXtVtMean-ZtMean[iNumTime]*ZtMean[iNumTime]);
	}
}


/* ----------------------------------------------------------------------------------------------- 
	

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEGridUnif_QBeta(/* Inputs */
								int iNb,
								double dMean,
								double dStd,
								double iNbSigmaLeft,
								double iNbSigmaRight,
									 
							 /* Outputs */
								double *Grid,
								int	*pIndex0)
{
/* Declaration of locals variables */
int iNumStep;
double dStep;
double dMin;

	/* Calculation of the step in X grid */
	if (iNb<2)
	{
		Grid[0]=0;
		*pIndex0 = 0;
	}
	else
	{
		dStep = ((iNbSigmaLeft+iNbSigmaRight)*dStd)/(iNb-1);

		/* Calculation of the min value of X */
		dMin = dMean - iNbSigmaLeft*dStd;

		/* Shift the grid so that X=0 is in the grid */
		if ((dMin <=0)&&(iNb>=2))
		{
			*pIndex0 =(int)(-dMin/dStep+1e-08);
			dMin = -(*pIndex0)*dStep;
		}
		else
		{
			dMin=0;
			*pIndex0 = 0;
		}
		
		/* Filling the grid in X */
		/* X = 0 is in the grid  */
		for (iNumStep=0; iNumStep<iNb; iNumStep++ )
			Grid[iNumStep]=dMin+iNumStep*dStep;
	}

}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEMakePhiGrid														  

	Compute the dynamic grid on Phi

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDEMakePhiGrid2_QBeta(/* Inputs */
								int iNbPhi,
								double dPhiMean,
								double dPhiStd,
								double iNbSigmaPhiGrid,

								/* Parameter */
								int IsExtremePoints,
								double iNbSigmaExtremePoints,
								int IsLogPhiGrid,
								
								/* Outputs */
								double *GridPhi,
								int	   *pIndexPhiMean)
{
/* Declaration of locals variables */
double dLogPhiStep, dLogPhiStd, dLogPhiMean, dLogPhiMin, dNewLogPhiMin;
double dPhiStep, dPhiMin, dNewPhiMin;
int iNumStep;

double PhiFloor = 0.0;	
double dPhiMinExtreme, dPhiMaxExtreme;
		
	if (iNbPhi == 1)
	{
		GridPhi[0] = dPhiMean;
		*pIndexPhiMean = 0;
	}
	else
	{
		switch (IsLogPhiGrid)
		{
		/*-------------
			Log Grid
		 --------------*/
		case LGMSV_TRUE:
			{
				if (dPhiMean ==0)
				{
					/* Filling the grid in Phi with constant value */
					for (iNumStep=0; iNumStep<iNbPhi; iNumStep++ )
						GridPhi[iNumStep] = dPhiMean;
				}
				else		
				{
					/* Calculation of mean and Std of log Phi*/
					dLogPhiStd = sqrt(log(1+dPhiStd*dPhiStd/dPhiMean/dPhiMean));
					dLogPhiMean = log(dPhiMean)-dLogPhiStd*dLogPhiStd/2;
					
					switch (IsExtremePoints)
					{
					case LGMSV_TRUE:
						{
							if (iNbPhi == 3)
							{
								GridPhi[0]=dPhiMean*exp(-iNbSigmaExtremePoints*dLogPhiStd);
								GridPhi[1]=dPhiMean;
								GridPhi[2]=dPhiMean*exp(iNbSigmaExtremePoints*dLogPhiStd);
								*pIndexPhiMean = 1;
							
							}
							else
							{
								/* Calculation of the step in LogPhi grid */
								dLogPhiStep = dPhiMean*(exp(iNbSigmaPhiGrid*dLogPhiStd)-
												exp(-iNbSigmaPhiGrid*dLogPhiStd))/(iNbPhi-3);

								/* Calculation of the min value of LogPhi */
								dLogPhiMin = dPhiMean*exp(-iNbSigmaPhiGrid*dLogPhiStd);

								/* Shift the grid so that dPhiMean is in the grid */
								*pIndexPhiMean =(int) floor((dPhiMean-dLogPhiMin)/dLogPhiStep);
								dNewLogPhiMin = dPhiMean-(*pIndexPhiMean)*dLogPhiStep;

								*pIndexPhiMean += 1;
							
								/* Filling the grid in Phi */
								for (iNumStep=1; iNumStep<iNbPhi-1; iNumStep++ )
									GridPhi[iNumStep]=dNewLogPhiMin+(iNumStep-1)*dLogPhiStep;
								GridPhi[0] = dPhiMean*exp(-iNbSigmaExtremePoints*dLogPhiStd);
								GridPhi[iNbPhi-1] = dPhiMean*exp(iNbSigmaExtremePoints*dLogPhiStd);
							}
							break;
						}
					
					default :
						{
							/* Calculation of the step in LogPhi grid */
							dLogPhiStep = (exp(dLogPhiMean+iNbSigmaPhiGrid*dLogPhiStd)-
												exp(dLogPhiMean-iNbSigmaPhiGrid*dLogPhiStd))/(iNbPhi-1);

							/* Calculation of the min value of LogPhi */
							dLogPhiMin = exp(dLogPhiMean-iNbSigmaPhiGrid*dLogPhiStd);

							/* Shift the grid so that log(dPhiMean) is in the grid */
							*pIndexPhiMean =(int) floor((dPhiMean-dLogPhiMin)/dLogPhiStep);
							dNewLogPhiMin = dPhiMean-(*pIndexPhiMean)*dLogPhiStep;
							
							/* Filling the grid in Phi */
							for (iNumStep=0; iNumStep<iNbPhi; iNumStep++ )
								GridPhi[iNumStep]=dNewLogPhiMin+(iNumStep-1)*dLogPhiStep;

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
						if (iNbPhi == 3)
						{
							/* Phi Should be positive */
							GridPhi[0]=max(dPhiMean-iNbSigmaExtremePoints*dPhiStd,PhiFloor);
							
							GridPhi[1]=dPhiMean;
							GridPhi[2]=dPhiMean+iNbSigmaExtremePoints*dPhiStd;
							*pIndexPhiMean = 1;
						}
						else
						{
							/* Calculation of the Phi Min Extreme point */
							dPhiMinExtreme = max(dPhiMean-iNbSigmaExtremePoints*dPhiStd,PhiFloor);
							dPhiMaxExtreme = dPhiMean+iNbSigmaExtremePoints*dPhiStd;
							
							/* Calculation of the step in Phi grid */
							dPhiStep = (dPhiMean+iNbSigmaPhiGrid*dPhiStd-max(dPhiMean-iNbSigmaPhiGrid*dPhiStd,PhiFloor))
								/(iNbPhi-3);

							/* Calculation of the min value of Phi */
							dPhiMin = max(dPhiMean - iNbSigmaPhiGrid*dPhiStd,PhiFloor);

							/* Shift the grid so that dPhiMean is in the grid */
							*pIndexPhiMean =(int) floor((dPhiMean-dPhiMin)/dPhiStep);
							dNewPhiMin = dPhiMean-(*pIndexPhiMean)*dPhiStep;

							
							/* Filling the grid in Phi */
							if (dNewPhiMin != dPhiMinExtreme)
							{
								*pIndexPhiMean += 1;
								
								for (iNumStep=1; iNumStep<iNbPhi-1; iNumStep++ )
									GridPhi[iNumStep]=dNewPhiMin+(iNumStep-1)*dPhiStep;
							
								GridPhi[0] = dPhiMinExtreme;
								GridPhi[iNbPhi-1] = dPhiMaxExtreme;
							}
							else
							{	
								for (iNumStep=0; iNumStep<iNbPhi-1; iNumStep++ )
									GridPhi[iNumStep]=dNewPhiMin+iNumStep*dPhiStep;
							
								GridPhi[iNbPhi-1] = dPhiMaxExtreme;
							}

						}
						break;
					}
				default:
					{
						/* Calculation of the step in Phi grid */
						dPhiStep = (dPhiMean+iNbSigmaPhiGrid*dPhiStd-max(dPhiMean-iNbSigmaPhiGrid*dPhiStd,PhiFloor))
							/(iNbPhi-1);

						/* Calculation of the min value of Phi */
						dPhiMin = max(dPhiMean - iNbSigmaPhiGrid*dPhiStd,PhiFloor);

						/* Shift the grid so that dPhiMean is in the grid */
						*pIndexPhiMean =(int) floor((dPhiMean-dPhiMin)/dPhiStep);
						dPhiMin = dPhiMean-(*pIndexPhiMean)*dPhiStep;
							
						/* Filling the grid in Phi */
						for (iNumStep=0; iNumStep<iNbPhi; iNumStep++ )
							GridPhi[iNumStep]=dPhiMin+iNumStep*dPhiStep;
						break;
					}
				}
				break;
			}
		}	
	}
}


/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDEInterpolate													  

	Interpolation

   ------------------------------------------------------------------------------------------------ */
static double LGMSVPDEInterpolate_QBeta(/* Inputs */
								   double x0,
								   int iNbfx,
								   double *x,
								   double ****fx,
								   int iNumX,
								   int iNumEps,
								   int iNumProduct,
								   int InterpolationType,
								   int IsExtrapolationFlat,
								   int iIndex)
{
/* declaration of locals variables */
double fx0 ;
int iIndex1, iIndex2, iIndex3;

	/* Extrapolation flat */
	if (((x0<x[0]) || (x0>x[iNbfx-1])) && (IsExtrapolationFlat == 1))
	{
		if (x0>x[iNbfx-1])
			fx0 = fx[iNbfx-1][iNumX][iNumEps][iNumProduct];
		else
			fx0 = fx[0][iNumX][iNumEps][iNumProduct];
	}
	else
	{
		/* Interpolation Treatment */
		switch (InterpolationType)
		{
			case LGMSV_LINEAR:
			{
				/* linear interpolation */
				if (iNbfx == 1) 
					fx0 = fx[0][iNumX][iNumEps][iNumProduct];
				else
				{
					if (x[0]==x[1])	
						fx0 = fx[0][iNumX][iNumEps][iNumProduct];
					else
					{
						if (iIndex ==0)
							/* Linear Extrapolation */
							fx0 = (x0-x[iIndex])/(x[iIndex+1]-x[iIndex])*(fx[iIndex+1][iNumX][iNumEps][iNumProduct]-
								fx[iIndex][iNumX][iNumEps][iNumProduct]) + fx[iIndex][iNumX][iNumEps][iNumProduct];
						else
							fx0 = (x0-x[iIndex-1])/(x[iIndex]-x[iIndex-1])*(fx[iIndex][iNumX][iNumEps][iNumProduct]-
								fx[iIndex-1][iNumX][iNumEps][iNumProduct]) + fx[iIndex-1][iNumX][iNumEps][iNumProduct];
					}
				}
				break;
			}
			
			case  LGMSV_QUADRATIC:
			{
				/* Quadratic interpolation by default */
				if (iNbfx == 1) 
					fx0 = fx[0][iNumX][iNumEps][iNumProduct];
				else
				{
					if (x[0]==x[1])
						fx0 = fx[0][iNumX][iNumEps][iNumProduct];
					else
					{
						/* Calculation of the 3 nearest points */
							
						if ((iIndex == 0) || (iIndex ==1))
						{
							iIndex1 = 0;
							iIndex2 = 1;
							iIndex3 = 2;
						}
						else
						{
							if (iIndex == iNbfx-1)
							{
								iIndex1 = iNbfx-3;
								iIndex2 = iNbfx-2;
								iIndex3 = iNbfx-1;
							}
							else
							{
								if (fabs(x0-x[iIndex-2])>fabs(x0-x[iIndex+1]))
								{
									iIndex1 = iIndex-1;
									iIndex2 = iIndex;
									iIndex3 = iIndex+1;
								}
								else
								{
									iIndex1 = iIndex-2;
									iIndex2 = iIndex-1;
									iIndex3 = iIndex;
								}
							}
						}

						/* Interpolation */
						fx0 = fx[iIndex1][iNumX][iNumEps][iNumProduct]/(x[iIndex1]-x[iIndex2])/(x[iIndex1]-x[iIndex3])*
								(x0*x0-(x[iIndex2]+x[iIndex3])*x0+x[iIndex2]*x[iIndex3])+
							fx[iIndex2][iNumX][iNumEps][iNumProduct]/(x[iIndex2]-x[iIndex3])/(x[iIndex2]-x[iIndex1])*
								(x0*x0-(x[iIndex3]+x[iIndex1])*x0+x[iIndex3]*x[iIndex1])+
							fx[iIndex3][iNumX][iNumEps][iNumProduct]/(x[iIndex3]-x[iIndex1])/(x[iIndex3]-x[iIndex2])*
								(x0*x0-(x[iIndex1]+x[iIndex2])*x0+x[iIndex1]*x[iIndex2]);
					}
				}			
				break;
			}
			
			default:
			{
				/*  */
				if (iNbfx == 1) 
					fx0 = fx[0][iNumX][iNumEps][iNumProduct];
				else
				{
					if (x[0]==x[1])
						fx0 = fx[0][iNumX][iNumEps][iNumProduct];
					else
					{
						if ((iNbfx>5) && (x0>x[2]) && (x0<x[iNbfx-3]))
						{
							/* Interpolation quadratic */
							
							/* Calculation of the 3 nearest points */
								
							if ((iIndex == 0) || (iIndex ==1))
							{
								iIndex1 = 0;
								iIndex2 = 1;
								iIndex3 = 2;
							}
							else
							{
								if (iIndex == iNbfx-1)
								{
									iIndex1 = iNbfx-3;
									iIndex2 = iNbfx-2;
									iIndex3 = iNbfx-1;
								}
								else
								{
									if (fabs(x0-x[iIndex-2])>fabs(x0-x[iIndex+1]))
									{
										iIndex1 = iIndex-1;
										iIndex2 = iIndex;
										iIndex3 = iIndex+1;
									}
									else
									{
										iIndex1 = iIndex-2;
										iIndex2 = iIndex-1;
										iIndex3 = iIndex;
									}
								}
							}

							/* Interpolation */
							fx0 = fx[iIndex1][iNumX][iNumEps][iNumProduct]/(x[iIndex1]-x[iIndex2])/(x[iIndex1]-x[iIndex3])*
									(x0*x0-(x[iIndex2]+x[iIndex3])*x0+x[iIndex2]*x[iIndex3])+
								fx[iIndex2][iNumX][iNumEps][iNumProduct]/(x[iIndex2]-x[iIndex3])/(x[iIndex2]-x[iIndex1])*
									(x0*x0-(x[iIndex3]+x[iIndex1])*x0+x[iIndex3]*x[iIndex1])+
								fx[iIndex3][iNumX][iNumEps][iNumProduct]/(x[iIndex3]-x[iIndex1])/(x[iIndex3]-x[iIndex2])*
									(x0*x0-(x[iIndex1]+x[iIndex2])*x0+x[iIndex1]*x[iIndex2]);
						}
						else
						{
							/* linear interpolation */
							if (iIndex ==0)
								/* Linear Extrapolation */
								fx0 = (x0-x[iIndex])/(x[iIndex+1]-x[iIndex])*(fx[iIndex+1][iNumX][iNumEps][iNumProduct]-
									fx[iIndex][iNumX][iNumEps][iNumProduct]) + fx[iIndex][iNumX][iNumEps][iNumProduct];
							else
								fx0 = (x0-x[iIndex-1])/(x[iIndex]-x[iIndex-1])*(fx[iIndex][iNumX][iNumEps][iNumProduct]-
									fx[iIndex-1][iNumX][iNumEps][iNumProduct]) + fx[iIndex-1][iNumX][iNumEps][iNumProduct];

						}
					}
				}			
				break;
			}
		}

	}
	if (fx0<0)
	{
		fx0 = fx0;
	}

	return fx0;
}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDECalculationPayoffMatrixfromPhit												  

	

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDECalculationValueMatrixfromPhit_QBeta(
											 /* Inputs */
											 int	iIndexMinX,
											 int	iIndexMaxX,
											 int	iIndexMinZ,
											 int	iIndexMaxZ,

											 int	iNbPhi,
											 int	iProductNb,
											 double	*GridXCenter, 
											 double	*GridZCenter,
											 double	*GridPhi_t1,

											 double dVFloor,

											 /* 4 dimensions : Phit,Xt,Zt,Product	*/
											 double ****ValueTensor,

											 double dXtMean,
											 double dZtMean,

											 double dPhit,
											 double dSigma,
											 double dSigt,
											 double dLambdaX,
											 double dAlpha,
											 double	dRho,

											 double t, /*time(t) */
											 double t1, /*time(t+1) */

											 int TypeInterpolationPhi,
											 int IsExtrapolationFlat,

											 /* Output */
											 /* 3 dimensions : Xt,Zt,Product	*/
											 double ***ValueGridXtZtProduct)
{
/* Declaration of locals variables */
double dDt, dVt, dPhi_t1, dXt, dZt;
int iNumZ, iNumX, iNumProduct;
int iIndexPhi;
	
	/* variation of time between t and t+1 */
	dDt = t1-t;
	
	for (iNumZ=iIndexMinZ;iNumZ<=iIndexMaxZ;iNumZ++)
	{
		for (iNumX=iIndexMinX;iNumX<=iIndexMaxX;iNumX++)
		{
			/* Calculation of Xt from the grid of Xcenter */
			dXt = GridXCenter[iNumX]+dXtMean;

			/* Calculation of Zt from the grid of Zcenter */
			dZt = GridZCenter[iNumX]+dZtMean;
			
			/* Calculation of the Vt from  the couple (Xt,Zt) */
			dVt = max(dZt+dAlpha*dRho/dSigt*dXt,dVFloor);
			
			/* Calculation of the Phi(t+1) */
			dPhi_t1=dPhit+(dSigma*dSigt*dSigma*dSigt*dVt-2*dLambdaX*dPhit)*dDt;

			/* Find the Index such that GridPhi_t1[index-1]< dPhi_t1 and GridPhi_t1[index]>= dPhi_t1 */
			iIndexPhi=0;
			while ((iIndexPhi <iNbPhi-1) && (GridPhi_t1[iIndexPhi]<dPhi_t1))
				iIndexPhi++;

			/* Interpolate value */
			
			for (iNumProduct=0;iNumProduct<iProductNb;iNumProduct++)
			{
				ValueGridXtZtProduct[iNumX][iNumZ][iNumProduct] = LGMSVPDEInterpolate_QBeta(/* Inputs */
									   dPhi_t1,
									   iNbPhi,
									   GridPhi_t1,
									   ValueTensor,
									   iNumX,
									   iNumZ,
									   iNumProduct,
									   TypeInterpolationPhi,
									   IsExtrapolationFlat,
									   iIndexPhi);

			}

		}
	}


}

/* ----------------------------------------------------------------------------------------------- 
	LGMSVPDECalculationRDriftVar												  

	

   ------------------------------------------------------------------------------------------------ */
static void LGMSVPDECalculationRDriftVar_QBeta(/* Inputs */
										int iIndexMinX,
										int	iIndexMaxX,
										int	iIndexMinZ,
										int	iIndexMaxZ,
										
										double	*GridXCenter,
										double	*GridZCenter,
										double dPhit,
										double dSigma,
										double dSigt,
										double dSig_t1,
										double dXtMean,
										double dZtMean,
										double dPhitMean,
										double dLambdaX,
										double dLambdaEps,
										double dAlpha,
										double dRho,
										double dIfrt,
										double t, /*time(t) */
										double t1, /*time(t+1) */

										double dVFloor,
																				
										/* Outputs */
										double **DXDrift,
										double **DZDrift,
										double **DXVar,
										double **DZVar,
										double **r)
{
/* Declaration of locals variables */
double dDt;
double dXt, dVt, dZt;
int iNumX, iNumZ;

	/* variation of time between t and t+1 */
	dDt = t1-t;
	
	for (iNumX=iIndexMinX;iNumX<=iIndexMaxX;iNumX++)
	{
		/* Calculation of Xt from the grid of XCenter */
		dXt = GridXCenter[iNumX]+dXtMean; 
		
		for (iNumZ=iIndexMinZ;iNumZ<=iIndexMaxZ;iNumZ++)
		{
			dZt = GridZCenter[iNumZ]+dZtMean;
						
			/* Calculation of the Vt from  the couple (Xt,Zt) */
			dVt = max(dZt+dAlpha*dRho/dSigt*dXt,dVFloor);
			
			/* Discount rate * dt */
			r[iNumX][iNumZ] =(dXt+dIfrt)*dDt;
						
			/* dt*Drifts under Qt+1 and dt*Variances under Qt+1 */
			
			/* Et+1 [Xt+1 - Xt / Ft] and Vt+1 [Xt+1 - Xt / Ft]*/
			DXDrift[iNumX][iNumZ] = ((dPhit-dPhitMean)-dLambdaX*(dXt-dXtMean))*dDt;
			DXVar[iNumX][iNumZ] = dSigt*dSigt*dVt*dDt;
			
			/* Et+1 [Zt+1 - Zt / Ft] and Vt+1 [Zt+1 - Zt / Ft]*/
			DZDrift[iNumX][iNumZ] = -dLambdaEps*(dVt-1)*dDt-dAlpha*dRho*(
				(1/dSig_t1-1/dSigt)*(dXt-dXtMean)+((dPhit-dPhitMean)-dLambdaX*(dXt-dXtMean))/dSigt*dDt);

			DZVar[iNumX][iNumZ] = dAlpha*dAlpha*dVt*(1-dRho*dRho)*dDt;
		}
	}
}

void LGMSVSaveTensor_QBeta( /* Inputs */
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
																				
	Xt= rt - f(0,t)															
	dXt =  [Phit-lambda_X*Xt]dt + sigma(Xt,t,omega) dWt						
	dPhit = [sigma*sigma - 2*lambda_X*Phit] dt								
																				
	sigma : time, spot dependant and stochastic of the following form		
		sigma(Xt,t,omega)  = sig(t)*sigma* Epst							

		Where																	
			The stochastic part of the sigma : Epst follows :					
				dEpst = -lambda_Eps*[Epst-1]dt + alpha * Epst * dZt						
				(mean reverting to 1, log normal vol)									
																			
   ------------------------------------------------------------------------------------------------ */

Err	 lgmSV_adi_QBeta(	
					/*	Time Information  */
					int			iNbTime,					
					double		*dTime,
					double		*dDate,

					/*	Space Discretisation	*/
					int			iNbPhi,
					int			iNbX,
					int			iNbEps,
						
					/*	Model data Information	*/
					double		dLambdaX,
					double		dSigma,
					
					double		*SigTime,
					double		*Sig,
					int			iNbSigTime,
					
					double		dLambdaEps,
					double		dAlpha,
					double		dRho,

					/* Parameters */
					LGMSVParam Params,

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
										
										/* Grid data Information	*/
										int		iIndPhitMin,
										int		iIndPhitMax,
										int		iIndXtMin,
										int		iIndXtMax,
										int		iIndEpstMin,
										int		iIndEpstMax,
										
										double	*GridPhi,
										double	*GridX,
						
										/* Tensor of results to be updated		*/
										/* 4 dimensions : Phit,Xt,Epst,Product	*/
										int		iNbProduct,
										double	****PayoffTensor),
					/*	Result */
					int			iNbProduct, 
					double		*dProductArrayPv)
{

/* Declaration of locals variables */
double			*Sigt = NULL;

double			*GridX = NULL;
double			*GridXCenter = NULL;
double			*GridZCenter = NULL;

double			*GridPhi_t = NULL;
double			*GridPhi_t1 = NULL;

double			*PhiMean = NULL;
double			*PhiStd = NULL;
double			*XtMean = NULL;
double			*XtStd = NULL;
double			*ZtMean = NULL;
double			*ZtStd = NULL;

double			**DXDrift = NULL;
double			**DZDrift = NULL;
double			**DXVar = NULL;
double			**DZVar = NULL;
double			**r = NULL;

double			***ValueGridXtZtProduct = NULL;
double			****ValueTensor_t = NULL;
double			****ValueTensor_t1 = NULL;

Err				err = NULL;
clock_t			t1;

double			dLambdaXGrid, dLambdaEpsGrid;

CNPDE_TEMP_2D_ADI	pdestr, *pde = &pdestr; 
/*CNPDE_TEMP_2D_LOD	pdestr, *pde = &pdestr;*/

int iNumTime, iNumPhi, iNumX, iNumEps, iNumProduct;
int iIndexZ0, iIndexX0, iIndexPhi;

int iNbPhiGrid_t, iNbPhiGrid_t1;
int iIndexMinX, iIndexMaxX, iIndexMinZ, iIndexMaxZ;


	/* Inputs Verification */	

	/* grid should be at least a 3x3 in X, Eps */
	iNbX = max(iNbX,3);
	iNbEps = max(iNbEps,3);

	/* All Number of points in the grid must be odd */
	iNbX = ((int) (iNbX/2))*2 +1;
	iNbEps = ((int) (iNbEps/2))*2 +1;
	iNbPhi = ((int) (iNbPhi/2))*2 +1;

	/* For computational time calculation */
	t1 = clock();
	
	/*  Memory allocations */
	Sigt = dvector(0, iNbTime-1);

	GridX = dvector(0, iNbX-1);
	GridXCenter = dvector(0, iNbX-1);
	GridZCenter = dvector(0, iNbEps-1);
	
	GridPhi_t = dvector(0, iNbPhi-1);	/* Grid of Phi at t   */
	GridPhi_t1 = dvector(0, iNbPhi-1);	/* Grid of Phi at t+1 */

	PhiMean = dvector(0, iNbTime-1);
	PhiStd = dvector(0, iNbTime-1);
	
	XtMean = dvector(0, iNbTime-1);
	XtStd = dvector(0, iNbTime-1);

	ZtMean = dvector(0, iNbTime-1);
	ZtStd = dvector(0, iNbTime-1);

	DXDrift = dmatrix(0, iNbX-1, 0, iNbEps-1);
	DZDrift = dmatrix(0, iNbX-1, 0, iNbEps-1);
	DXVar = dmatrix(0, iNbX-1, 0, iNbEps-1);
	DZVar = dmatrix(0, iNbX-1, 0, iNbEps-1);
	r = dmatrix(0, iNbX-1, 0, iNbEps-1);
	
	ValueGridXtZtProduct =f3tensor(0, iNbX-1,
									  0, iNbEps-1,
									  0, iNbProduct-1);
		
	ValueTensor_t =f4tensor(0, iNbPhi-1,
							0, iNbX-1,
							0, iNbEps-1,
							0, iNbProduct-1);
	
	ValueTensor_t1 =f4tensor(0, iNbPhi-1,
							 0, iNbX-1,
							 0, iNbEps-1, 
							 0, iNbProduct-1);

	/* Gestion of allocation errors */
	if (!pde || !GridX || !GridXCenter ||  !Sigt || !GridZCenter || !GridPhi_t || !GridPhi_t1 || !PhiMean || !PhiStd 
		 || !XtMean || !XtStd || !ZtMean || !ZtStd|| !DXDrift || !DZDrift || !DXVar || !DZVar || !r || !ValueGridXtZtProduct
		 || !ValueTensor_t || !ValueTensor_t1)
	{
		err = "Memory allocation error (1) in LGMSVPDE";
		goto FREE_RETURN;
	}

	/* Constant	*/
	dLambdaXGrid = dLambdaX;			/* Use of LambdaX for the grid */
	dLambdaEpsGrid = dLambdaEps;		/* Use of LambdaEps for the grid */ 


	
	/* Expectations calculations */
	LGMSVPDEExpectations_QBeta(
					/* Inputs */
					 iNbTime,			
					 dTime,
					
					 /* For X diffusion */
					 dLambdaXGrid,
					 dSigma,
					 
					 SigTime,
					 Sig,
					 iNbSigTime,

					 /* For Eps diffusion */
					 dLambdaEpsGrid,
					 dAlpha,
					 dRho,
					 
					 /* Outputs */
					 Sigt,

					 PhiMean,
					 PhiStd,
					 XtMean,
					 XtStd,
					 ZtMean,
					 ZtStd);


	/* Making the grid on Xcenter and ZCenter */	
								
	LGMSVPDEGridUnif_QBeta(/* Inputs */
						iNbEps,
						0,
						ZtStd[iNbTime-1],

						Params.iNbSigmaZGridLeft,
						Params.iNbSigmaZGridRight,
									 
						/* Outputs */
						GridZCenter,
						&iIndexZ0); 

	LGMSVPDEGridUnif_QBeta(/* Inputs */
						iNbX,
						0,
						XtStd[iNbTime-1],
						Params.iNbSigmaXGridLeft,
						Params.iNbSigmaXGridRight,
									 
						/* Outputs */
						GridXCenter,
						&iIndexX0); 

	/* Backward PDE */

	/* Init of the Backward PDE : Allocation of arrays    */
	num_f_pde_init_2d_adi(pde,
						iNbX,
						iNbEps,
						iNbProduct); 
	/*num_f_pde_init_2d_lod(pde,
						iNbX,
						iNbEps,
						iNbProduct);*/

	if (!pde)
	{
		err = "Memory allocation error (2) in LGMSVPDE";
		goto FREE_RETURN;
	}

	/* Final Grid of Phi , Grid of phi at time T */
	LGMSVPDEMakePhiGrid2_QBeta(/* Inputs */
						iNbPhi,
						PhiMean[iNbTime-1],
						PhiStd[iNbTime-1],
						Params.iNbSigmaPhiGrid,

						/* Parameter */
						Params.IsExtremePoints,
						Params.iNbSigmaExtremePoints,
						Params.IsPhiBoundExp,
								
						/* Outputs */
						GridPhi_t1,
						&iIndexPhi);

	
	/* Final Payoff valuation   */
	for (iNumX = 0; iNumX<iNbX; iNumX++)
		GridX[iNumX]=GridXCenter[iNumX]+XtMean[iNbTime-1]; 

	err = payoff_func(/* Time Information */
						dDate[iNbTime-1],
						dTime[iNbTime-1],
						func_parm_tab[iNbTime-1],
										
						/* Market data	*/										
						cYieldCurve,
										
						/*	Model data Information	*/
						dLambdaX,			
										
						/* Grid data Information	*/
						0,
						iNbPhi-1,
						0,
						iNbX-1,
						0,
						iNbEps-1,
										
						GridPhi_t1,
						GridX,
						
						/* Tensor of results to be updated		*/
						/* 4 dimensions : Phit,Xt,Epst,Product	*/
						iNbProduct,
						ValueTensor_t1);
	
	if (Params.VerifExpectation==LGMSV_TRUE)
	{
		for (iNumPhi = 0; iNumPhi<iNbPhi; iNumPhi++)
			for (iNumX = 0; iNumX<iNbX; iNumX++)
				for (iNumEps = 0; iNumEps<iNbEps; iNumEps++)
					for (iNumProduct = 0; iNumProduct<1; iNumProduct++)
					ValueTensor_t1[iNumPhi][iNumX][iNumEps][iNumProduct]=GridZCenter[iNumEps]*GridZCenter[iNumEps];
	}

	/* Init Nb Phi at t+1 */
	iNbPhiGrid_t1 = iNbPhi;


if (Params.ModeExport!=LGMSV_TRUE)
{	
	/* For each time t backward */
	for (iNumTime=iNbTime-2; iNumTime>=0; iNumTime--)
	{
		/* Update Parameters for calculation optimisation */
		
		/* Number of Phi in the grid at time t */
		iNbPhiGrid_t = iNbPhi-(Params.GridPhiOptim == LGMSV_TRUE)*(int)((1-iNumTime/((double)(iNbTime-1)))
			*(iNbPhi-1) + 1.0e-08);
		iNbPhiGrid_t = ((int) (iNbPhiGrid_t/2))*2 +1;
		
		
		if (Params.GridXADIOptim == LGMSV_TRUE)
		{
			iIndexMinX = (int) ((-Params.iNbSigmaXGridLeft*XtStd[iNumTime+1]-GridXCenter[0])
				/(GridXCenter[1]-GridXCenter[0]) + 1.0e-08);
			iIndexMaxX = (int) ((Params.iNbSigmaXGridRight*XtStd[iNumTime+1]-GridXCenter[0])
				/(GridXCenter[1]-GridXCenter[0]) + 1.0e-08);
		}
		else
		{
			iIndexMinX = 0;
			iIndexMaxX = iNbX-1;
		}

		if (Params.GridZADIOptim == LGMSV_TRUE)
		{
			iIndexMinZ = (int) ((-Params.iNbSigmaZGridLeft*ZtStd[iNumTime+1]-GridZCenter[0])
				/(GridZCenter[1]-GridZCenter[0]) + 1.0e-08);
			iIndexMaxZ = (int) ((Params.iNbSigmaZGridRight*ZtStd[iNumTime+1]-GridZCenter[0])
				/(GridZCenter[1]-GridZCenter[0]) + 1.0e-08);
		}
		else
		{	
			iIndexMinZ = 0;
			iIndexMaxZ = iNbEps-1;
		}

		if (Params.SaveTensor == LGMSV_TRUE)
		{
			/* Save the Tensor */
			LGMSVSaveTensor_QBeta( /* Inputs */
					iNbPhi,
					iNbX,
					iNbEps,
					 
					GridPhi_t1,
					GridXCenter,
					GridZCenter,
					ValueTensor_t1);
		}
			
		/* Make the grid in Phit at time t */
		LGMSVPDEMakePhiGrid2_QBeta(/* Inputs */
						iNbPhiGrid_t,
						PhiMean[iNumTime],
						PhiStd[iNumTime],
						Params.iNbSigmaPhiGrid,

						/* Parameter */
						Params.IsExtremePoints,
						Params.iNbSigmaExtremePoints,
						Params.IsPhiBoundExp,
								
						/* Outputs */
						GridPhi_t,
						&iIndexPhi);

		/* For each Phit in the grid at time t*/
		for (iNumPhi=0;iNumPhi<iNbPhiGrid_t;iNumPhi++)
		{		
				
			/* Calculation of the Value at t+1 on the grid (XCenter,ZCenter)			 */
			/* at the deduced forwards Phi t+1 from Phit								 */
			/* Interpolation From the value in the grid (XCenter,ZCenter,Phi) at t+1	 */
			LGMSVPDECalculationValueMatrixfromPhit_QBeta(
											 /* Inputs */
											 iIndexMinX,
											 iIndexMaxX,
											 iIndexMinZ,
											 iIndexMaxZ,

											 iNbPhiGrid_t1,
											 iNbProduct,
											 GridXCenter, 
											 GridZCenter,
											 GridPhi_t1,

											 Params.VFloor,

											 ValueTensor_t1,
											 
											 XtMean[iNumTime],
											 ZtMean[iNumTime],
											 
											 GridPhi_t[iNumPhi],
											 dSigma,
											 Sigt[iNumTime],
											 dLambdaX,
											 dAlpha,
											 dRho,

											 dTime[iNumTime], /*time(t) */
											 dTime[iNumTime+1], /*time(t+1) */

											 Params.TypeInterpolationPhi,
											 Params.IsExtrapolationFlat,

											 /* Output */
											 ValueGridXtZtProduct);

			/* One Step Backward 2D using Crank-Nicolson on grid (Xcenter,Zcenter) */	
			/* From t+1 to t											           */

			/* Calculation of the drifts and variances under Qt+1 */
			LGMSVPDECalculationRDriftVar_QBeta(/* Inputs */
										iIndexMinX,
										iIndexMaxX,
										iIndexMinZ,
										iIndexMaxZ,

										GridXCenter, 
										GridZCenter,
										GridPhi_t[iNumPhi],
										dSigma,
										Sigt[iNumTime],
										Sigt[iNumTime+1],
										XtMean[iNumTime],
										ZtMean[iNumTime],
										PhiMean[iNumTime],
										dLambdaX,
										dLambdaEps,
										dAlpha,
										dRho,
										Ifr[iNumTime],
										dTime[iNumTime], /*time(t) */
										dTime[iNumTime+1], /*time(t+1) */

										Params.VFloor,

										/* Outputs */
										DXDrift,
										DZDrift,
										DXVar,
										DZVar,
										r);

			/* Launch The Crank-Nicolson					*/
			/* Calculation of the payoff tensor at time t	*/
			/* for the correspondant phi(t)					*/
			num_f_pde_one_step_backward_2f_adi(	pde,
												iNbX,
												GridXCenter,
												iNbEps,
												GridZCenter,
												0,
												iNbProduct-1,
												ValueGridXtZtProduct,
												DXDrift,
												DZDrift,
												DXVar,
												DZVar,
												r,
												ValueTensor_t[iNumPhi],
												iIndexMinX,
												iIndexMaxX,
												iIndexMinZ,
												iIndexMaxZ); 

			/*num_f_pde_one_step_backward_2f_lod(	pde,
												iNbX,
												GridXCenter,
												iNbEps,
												GridZCenter,
												iNbProduct,
												ValueGridXtZtProduct,
												DXDrift,
												DZDrift,
												DXVar,
												DZVar,
												r,
												0.5,
												ValueTensor_t[iNumPhi],
												0,
												iNbX-1,
												0,
												iNbEps-1);*/

		}

		/* End for each Phit */

		/*	Case of evaluation events */
		if (EvalEvent[iNumTime])
		{
			/* Calculation of the grid in X from the gridXCenter */
			for (iNumX = 0; iNumX<iNbX; iNumX++)
				GridX[iNumX]=GridXCenter[iNumX]+XtMean[iNumTime]; 
			
			
			/* Modification of the Payoff at t */			
			err = payoff_func(/* Time Information */
								dDate[iNumTime],
								dTime[iNumTime],
								func_parm_tab[iNumTime],
												
								/* Market data	*/										
								cYieldCurve,
												
								/*	Model data Information	*/
								dLambdaX,			
												
								/* Grid data Information	*/
								0,
								iNbPhiGrid_t-1,
								iIndexMinX,
								iIndexMaxX,
								iIndexMinZ,
								iIndexMaxZ,
												
								GridPhi_t,
								GridX,
								
								/* Tensor of results to be updated		*/
								/* 4 dimensions : Phit,Xt,Epst,Product	*/
								iNbProduct,
								ValueTensor_t);
			
			if (err)
			{
					goto FREE_RETURN;
			}
		}

		/* Swaping, t+1 becomes t */
		LGMSVPDESwap(double ****,ValueTensor_t,ValueTensor_t1);
		LGMSVPDESwap(double *,GridPhi_t,GridPhi_t1);
		LGMSVPDESwap(int ,iNbPhiGrid_t,iNbPhiGrid_t1);

	}
	/* End for each t */

	/* Copy of the result */
	/* The grids are made such as 
		center on Phi = 0, X=0,  at t=0 */
	memcpy(dProductArrayPv,ValueTensor_t1[iIndexPhi][iIndexX0][iIndexZ0],iNbProduct*sizeof(double));
}
else
{
	dProductArrayPv[0]=	PhiStd[iNbTime-1]; /*dZTStd; /*PhiStd[iNbTime-1];*/
	dProductArrayPv[1]=PhiMean[iNbTime-1]; /*dXTStd; */
}


FREE_RETURN:

	/* free all the arrays */
	if (pde) num_f_pde_free_2d_adi(pde, iNbX, iNbEps, iNbProduct); 
	/*if (pde) num_f_pde_free_2d_lod(pde, iNbX, iNbEps, iNbProduct); */
	
	if (Sigt)
		free_dvector(Sigt, 0, iNbTime-1);

	if (GridX)
		free_dvector(GridX, 0, iNbX-1);
	if (GridXCenter)
		free_dvector(GridXCenter, 0, iNbX-1);
	if (GridZCenter)
		free_dvector(GridZCenter, 0, iNbEps-1);

	if (GridPhi_t)
		free_dvector(GridPhi_t, 0, iNbPhi-1);
	if (GridPhi_t1)
		free_dvector(GridPhi_t1, 0, iNbPhi-1);

	if (PhiMean)
		free_dvector(PhiMean, 0, iNbTime-1);
	if (PhiStd)
		free_dvector(PhiStd, 0, iNbTime-1);

	if (XtMean)
		free_dvector(XtMean, 0, iNbTime-1);
	if (XtStd)
		free_dvector(XtStd, 0, iNbTime-1);
	
	if (ZtMean)
		free_dvector(ZtMean, 0, iNbTime-1);
	if (ZtStd)
		free_dvector(ZtStd, 0, iNbTime-1);

	if (DXDrift)
		free_dmatrix(DXDrift, 0, iNbX-1, 0, iNbEps-1);
	if (DZDrift)
		free_dmatrix(DZDrift, 0, iNbX-1, 0, iNbEps-1);
	if (DXVar)
		free_dmatrix(DXVar, 0, iNbX-1, 0, iNbEps-1);
	if (DZVar)
		free_dmatrix(DZVar, 0, iNbX-1, 0, iNbEps-1);
	if (r)
		free_dmatrix(r, 0, iNbX-1, 0, iNbEps-1);
	
	if (ValueGridXtZtProduct)
		free_f3tensor(ValueGridXtZtProduct, 0, iNbX-1, 0, iNbEps-1, 0, iNbProduct-1);
		
	if (ValueTensor_t)
		free_f4tensor(ValueTensor_t, 0, iNbPhi-1, 0, iNbX-1, 0, iNbEps-1, 0, iNbProduct-1);
	
	if (ValueTensor_t1)
		free_f4tensor(ValueTensor_t1,0, iNbPhi-1, 0, iNbX-1, 0, iNbEps-1,  0, iNbProduct-1);
			
	
	/* Return the error message */
	return err;

}


/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/

