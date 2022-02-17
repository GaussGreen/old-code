/* ==========================================================================
   FILE_NAME:	LGMSVClosedForm2.c

   PURPOSE:		
	
   DATE:		05/01/02
   
   AUTHOR:		L.C.
   ========================================================================== */

#include "LGMSVPDE.h"
#include "opfnctns.h"
#include "Fx3FUtils.h"
#include "FxSabrAdi.h"
#include "Fx3FCalib.h"
#include "lgmsvclosedform.h"


#define	PI	3.141592653589793238462643383279502884197169399375105820974944592308
#define PI2	6.28318530717958647692528676655900576839433879875021164194988918462
#define ONE_MONTH 0.083333333
#define MAX_CPN			600
#define	LGM_VEGA_SHIFT	0.025

static void prod_comp(double	re1,
					  double	im1,
					  double	re2,
					  double	im2,
					  double	*re3,
					  double	*im3)
{
	*re3 = re1 * re2 - im1 * im2;
	*im3 = re1 * im2 + im1 * re2;
}

void LGMSVMomentInit2(
					/* Inputs */					 
					double			dLambdaX,					 
					double			dAlpha,
					double			dLambdaEps,
					double			dRho,
										
					LGMSVSolFunc	*FuncPhi,
					LGMSVSolFunc	*FuncPhi2,
					LGMSVSolFunc	*FuncV2,
					LGMSVSolFunc	*FuncPhiV)
{	
	FuncPhi->dLambda = 2.0 * dLambdaX;
	FuncPhi->bIsft1 = 1;	
	FuncPhi->bIsgt1 = 1;
	FuncPhi->b = 0;
	FuncPhi->bIsht1 = 1;
	FuncPhi->c =0;	
	
	FuncPhi2->dLambda = 4.0 * dLambdaX;
	FuncPhi2->bIsft1 = 0;
	FuncPhi2->pft = FuncPhiV;
	FuncPhi2->bIsgt1 = 1;
	FuncPhi2->b = 0.0;
	FuncPhi2->bIsht1 = 1;
	FuncPhi2->c = 0.0;	

	FuncV2->dLambda = 2.0 * dLambdaEps;
	FuncV2->bIsft1 = 1;
	FuncV2->a = 2.0 * dLambdaEps + dAlpha * dAlpha;
	FuncV2->bIsgt1 = 1;
	FuncV2->b = 0.0;
	FuncV2->bIsht1 = 1;
	FuncV2->c = 0.0;	

	FuncPhiV->dLambda = 2.0 * dLambdaX + dLambdaEps;
	FuncPhiV->bIsft1 = 0;
	FuncPhiV->pft = FuncPhi;
	FuncPhiV->a = dLambdaEps;
	FuncPhiV->bIsgt1 = 0;
	FuncPhiV->pgt = FuncV2;
	FuncPhiV->bIsht1 = 1;
	FuncPhiV->c = 0.0;
}

void LGMSVMomentCalculation2(

					double			Sig2,
					double			dt,
					
					double			InitPhi,
					double			InitPhi2,
					double			InitV2,
					double			InitPhiV,
					
					LGMSVSolFunc	*FuncPhi,
					LGMSVSolFunc	*FuncPhi2,
					LGMSVSolFunc	*FuncV2,
					LGMSVSolFunc	*FuncPhiV,
					double			*lambda,
					
					double			*ResPhi,
					double			*ResPhi2,
					double			*ResV2,
					double			*ResPhiV)
{
	/* Initialisation */
	FuncPhi->a = Sig2;
	FuncPhi->dXt1 = InitPhi;

	FuncPhi2->a = 2.0 * Sig2;
	FuncPhi2->dXt1 = InitPhi2;

	FuncV2->dXt1 = InitV2;

	FuncPhiV->b = Sig2;
	FuncPhiV->dXt1 = InitPhiV;

	/* Calculation */
	LGMSVFuncValue2(FuncPhi, dt, lambda, 0, ResPhi);
	LGMSVFuncValue2(FuncPhi2, dt, lambda, 0, ResPhi2);
	LGMSVFuncValue2(FuncV2, dt, lambda, 0, ResV2);
	LGMSVFuncValue2(FuncPhiV, dt, lambda, 0, ResPhiV);
}



void	LGMSVSolveODELogTSEvalFunc(double	Coef1Re,
								   double	Coef2Re,
								   double	Coef2Im,
								   double	Coef3Re,
								   double	Coef3Im,
								   double	expt,
								   double	DRe,
								   double	DIm,
								   double	*ResRe,
								   double	*ResIm)
{
	prod_comp(Coef1Re * DRe + Coef2Re, Coef1Re * DIm + Coef2Im * expt, DRe, DIm, ResRe, ResIm);
	expt *= expt;
	*ResRe += Coef3Re * expt;
	*ResIm += Coef3Im * expt;
}


void LGMSVOptionPricePrecalc(	int		iNbPhi,
								int		iNbft,								

								double	dLambdaX,
								double	dTStar,
								int		endi,
								double	expt2init,
								double	*Coef1ReT,
								double	*Coef2ReT,
								double	*Coef2ImT,
								double	*Coef3ReT,
								double	*Coef3ImT,
								double	*CoefAT,
								double	*CoefexpT,
								int		*nbpT,
								
								double	dPhitMean,									
								int		SaveFile,
								double	***Density,
								double	**DensitySpeq,
								double	**PayOff,

								double	dExTime,
								int		iNbCoupon,								
								double	*CouponTime,
								double	*Coupon,
																
								double	dPhiFreqStep,
								double	dPhiStep,
								int		iIndexPhiMean,
								double	dftFreqStep,
								double	dftStep,
								int		iIndexft0,
								
								double	*Price)
{
int		iNumPhi, iNumft;
double	dIntegral;
clock_t	time1, time2, time3;

	time1 = clock();

	
	LGMSVCalculateDensityPrecalc(iNbPhi, iNbft, dPhiFreqStep, dftFreqStep, endi, expt2init,
								Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT, CoefexpT, nbpT,
								dPhitMean, SaveFile, Density, DensitySpeq);

	time2 = clock();

	LGMSVFillPayoff(
					iNbPhi,
					iNbft,
					dLambdaX,
					dPhitMean,
					iIndexPhiMean,
					dPhiStep,
					iIndexft0,
					dftStep,
					dExTime,
					dTStar,					
					iNbCoupon,
					CouponTime,
					Coupon,
					PayOff);

	/* Calculation of the Integral of the payoff times the density */
	dIntegral = 0;
	for (iNumPhi=1;iNumPhi<=iNbPhi;iNumPhi++)
		for (iNumft=1;iNumft<=iNbft;iNumft++)
			/* Calculation of the payoff */
				dIntegral += Density[1][iNumPhi][iNumft]*PayOff[iNumPhi][iNumft];
	dIntegral *= dPhiStep*dftStep;

	time3 = clock();

	/* Return the result */
	*Price = dIntegral;
}		
	

void LGMSVFindFreq(	int		iNbPhi,
					int		iNbft,
					
					double	AlphaEq,
					double	Rho,

					int		endi,
					double	expt2init,
					double	*Coef1ReT,
					double	*Coef2ReT,
					double	*Coef2ImT,
					double	*Coef3ReT,
					double	*Coef3ImT,
					double	*CoefAT,
					double	*CoefexpT,
					int		*nbpT,

					double	dPhitMean,
					double	dPhitStd,
					double	dftStd,

					double	iNbSigmaPhiGridLeft,
					double	iNbSigmaPhiGridRight,
					double	iNbSigmaftLeft,
					double	iNbSigmaftRight,
					double	iLimitPhi,
					double	iLimitft,
					int		iPriorityFreqPhi,
					int		iPriorityFreqFt,
					
					/* Outputs */
					double	*dPhiFreqStep,
					double	*dPhiStep,
					int		*iIndexPhiMean,
					double	*dftFreqStep,
					double	*dftStep,
					int		*iIndexft0)
{
double	dPhiMin, dPhiMax, dftMin;
double	dPhiFreqMin, dFtFreqMin;
double	LogTFRe, LogTFIm;
int		i, NbMinMean;


	/* Find Limit on Phi */
	if (iNbSigmaPhiGridLeft == 0)
	{
		iNbSigmaPhiGridLeft = 5.0;
	}
	if (iNbSigmaPhiGridRight == 0)
	{
		iNbSigmaPhiGridRight = 5.0 + 5.0 * AlphaEq;
	}

	dPhiMin = max(dPhitMean - iNbSigmaPhiGridLeft * dPhitStd, 0.0);
	dPhiMax = dPhitMean + iNbSigmaPhiGridRight * dPhitStd;
	*dPhiStep = (dPhiMax - dPhiMin) / iNbPhi;
	*iIndexPhiMean = max((int) ((dPhitMean - dPhiMin) / *dPhiStep + 0.5), 1) + 1;
	*dPhiStep = (dPhitMean - dPhiMin) / (*iIndexPhiMean - 1);
	
	if (iPriorityFreqPhi)
	{
		dPhiFreqMin = 1.0 / *dPhiStep * PI;
		
		LGMSVCalculateDensityTPPoint(dPhiFreqMin, 0.0, endi, expt2init, Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT,
									CoefexpT, nbpT, dPhitMean, &LogTFRe, &LogTFIm);

		i = 0;
		while (LogTFRe > iLimitPhi && i++ < 200)
		{
			dPhiFreqMin *= 1.2;
			LGMSVCalculateDensityTPPoint(dPhiFreqMin, 0.0, endi, expt2init, Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT,
									CoefexpT, nbpT, dPhitMean, &LogTFRe, &LogTFIm);			
		}

		dPhiFreqMin /= PI;
		*dPhiStep = 1.0 / dPhiFreqMin;
		*iIndexPhiMean = max((int) ((dPhitMean - dPhiMin) / *dPhiStep + 0.5), 1) + 1;
		*dPhiStep = (dPhitMean - dPhiMin) / (*iIndexPhiMean - 1);
	}
		
	*dPhiFreqStep = 1.0 / *dPhiStep / iNbPhi;
	NbMinMean = (int) ((dPhitMean - dPhiMin) / *dPhiStep) + 2;
	if (NbMinMean < iNbPhi)
	{
		*iIndexPhiMean = NbMinMean;
	}

	/* Find Limit on ft */
	if (iNbSigmaftLeft == 0)
	{
		iNbSigmaftLeft = 5.0 + 20 * AlphaEq - 10.0 * AlphaEq * Rho;
	}
	if (iNbSigmaftRight == 0)
	{
		iNbSigmaftRight = 5.0 + 20 * AlphaEq + 10.0 * AlphaEq * Rho;
	}

	*dftStep = (iNbSigmaftLeft + iNbSigmaftRight) * dftStd / iNbft;
	dftMin = -iNbSigmaftLeft * dftStd;
	*iIndexft0 = max((int) (-dftMin / *dftStep + 0.5), 1) + 1;
	*dftStep =  -dftMin / (*iIndexft0 - 1);
	
	if (iPriorityFreqFt)
	{
		dFtFreqMin = 1.0 / *dftStep * PI;
		
		LGMSVCalculateDensityTPPoint(0.0, dFtFreqMin, endi, expt2init, Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT,
									CoefexpT, nbpT, dPhitMean, &LogTFRe, &LogTFIm);

		i = 0;
		while (LogTFRe > iLimitft && i++ < 200)
		{
			dFtFreqMin *= 1.2;
			LGMSVCalculateDensityTPPoint(0.0, dFtFreqMin, endi, expt2init, Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT,
									CoefexpT, nbpT, dPhitMean, &LogTFRe, &LogTFIm);			
		}

		dFtFreqMin /= PI;
		*dftStep = 1.0 / dFtFreqMin;
		*iIndexft0 = max((int) (-dftMin / *dftStep + 0.5), 1) + 1;
		*dftStep = -dftMin / (*iIndexft0 - 1);
	}
		
	*dftFreqStep = 1.0 / *dftStep / iNbft;
	NbMinMean = (int) (-dftMin / *dftStep) + 2;
	if (NbMinMean < iNbft)
	{
		*iIndexft0 = NbMinMean;
	}
}

void LGMSVCalculateDensityTPPoint(	double	dPhifreq,
									double	dftFreq,
									
									int		endi,
									double	expt2init,
									double	*Coef1ReT,
									double	*Coef2ReT,
									double	*Coef2ImT,
									double	*Coef3ReT,
									double	*Coef3ImT,
									double	*CoefAT,
									double	*CoefexpT,
									int		*nbpT,

									double	PhiMean,

									/* Outputs */
									double	*LogTFRe,
									double	*LogTFIm)
{
double	ARe, AIm, DRe, DIm, DReTemp, DImTemp;
double	K1Re, K1Im, K2Re, K2Im, K3Re, K3Im, K4Re, K4Im;
double	CoefA, Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im;
double	expt2, coefexp;
int		nbp;

int		i, j;	

	/* Initialisation */	
	ARe = 0.0;
	AIm = -PhiMean * dPhifreq;
	DRe = 0.0;
	DIm = 0.0;
	expt2 = expt2init;
	
	for (i=endi; i>=0; i--)
	{					
		/* Solve the ODE between t1 and t2 */												
		coefexp = CoefexpT[i];
		Coef1Re = Coef1ReT[i];
		Coef2Re = Coef2ReT[i];
		Coef2Im = Coef2ImT[i] * dftFreq;
		Coef3Re = Coef3ReT[i] * dftFreq * dftFreq;
		Coef3Im = Coef3ImT[i] * dPhifreq;
		CoefA = CoefAT[i];
		nbp = nbpT[i];

		for (j=0; j<nbp; j++)
		{
			DReTemp = DRe;
			DImTemp = DIm;								
			
			LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
										expt2, DRe, DIm, &K1Re, &K1Im);
			expt2 *= coefexp;
			LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
										expt2, DRe + 0.5 * K1Re, DIm + 0.5 * K1Im,  &K2Re, &K2Im);
			LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
										expt2, DRe + 0.5 * K2Re, DIm + 0.5 * K2Im,  &K3Re, &K3Im);
			expt2 *= coefexp;
			LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
										expt2, DRe + K3Re, DIm + K3Im,  &K4Re, &K4Im);

			DRe += (K1Re + K4Re) / 6.0 + (K2Re + K3Re) / 3.0;
			DIm += (K1Im + K4Im) / 6.0 + (K2Im + K3Im) / 3.0;

			ARe += CoefA * (DRe + DReTemp);
			AIm += CoefA * (DIm + DImTemp);					
		}
	}

	*LogTFRe = ARe + DRe;
	*LogTFIm = AIm + DIm;
}

void LGMSVCalculateDensityPrecalc(	int		iNbPhi,
									int		iNbft,							
									double	dPhiFreqStep,
									double	dftFreqStep,
									
									int		endi,
									double	expt2init,
									double	*Coef1ReT,
									double	*Coef2ReT,
									double	*Coef2ImT,
									double	*Coef3ReT,
									double	*Coef3ImT,
									double	*CoefAT,
									double	*CoefexpT,
									int		*nbpT,

									double	PhiMean,									
									int		SaveFile,

									/* Outputs */
									double	***Density,
									double	**DensitySpeq)
{
double	ARe, AIm, DRe, DIm, DReTemp, DImTemp;
double	K1Re, K1Im, K2Re, K2Im, K3Re, K3Im, K4Re, K4Im;
double	dftFreq, dPhifreq, dExpReal, dImag, PhiTrans, PhiTransStep, JumpPhi, JumpPhiTrans;
double	CoefA, Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im;
double	dCoef, expt2, coefexp;
int		iNumPhiFreq, iNumftFreq, nbFtHalf, nbPhiHalf, nbp;

int		i, j;

	
	dCoef = 2.0 * dftFreqStep * dPhiFreqStep;

	dftFreqStep *= PI2;
	dPhiFreqStep *= PI2;

	nbFtHalf = iNbft / 2 + 1;
	nbPhiHalf = iNbPhi / 2;

	PhiTransStep = dPhiFreqStep * PhiMean;
	JumpPhi = iNbPhi * dPhiFreqStep;
	JumpPhiTrans = iNbPhi * PhiTransStep;
		
	dPhifreq = 0.0;
	PhiTrans = 0.0;

	for (iNumPhiFreq=1; iNumPhiFreq<=iNbPhi; iNumPhiFreq++)
	{					
		dftFreq = 0.0;		

		for (iNumftFreq=1; iNumftFreq<=nbFtHalf; iNumftFreq++)
		{
			/* Initialisation */	
			ARe = 0.0;
			AIm = -PhiTrans;
			DRe = 0.0;
			DIm = 0.0;
			expt2 = expt2init;
			
			for (i=endi; i>=0; i--)
			{					
				/* Solve the ODE between t1 and t2 */												
				coefexp = CoefexpT[i];
				Coef1Re = Coef1ReT[i];
				Coef2Re = Coef2ReT[i];
				Coef2Im = Coef2ImT[i] * dftFreq;
				Coef3Re = Coef3ReT[i] * dftFreq * dftFreq;
				Coef3Im = Coef3ImT[i] * dPhifreq;
				CoefA = CoefAT[i];
				nbp = nbpT[i];

				for (j=0; j<nbp; j++)
				{
					DReTemp = DRe;
					DImTemp = DIm;								
					
					LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
												expt2, DRe, DIm, &K1Re, &K1Im);
					expt2 *= coefexp;
					LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
												expt2, DRe + 0.5 * K1Re, DIm + 0.5 * K1Im,  &K2Re, &K2Im);
					LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
												expt2, DRe + 0.5 * K2Re, DIm + 0.5 * K2Im,  &K3Re, &K3Im);
					expt2 *= coefexp;
					LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
												expt2, DRe + K3Re, DIm + K3Im,  &K4Re, &K4Im);

					DRe += (K1Re + K4Re) / 6.0 + (K2Re + K3Re) / 3.0;
					DIm += (K1Im + K4Im) / 6.0 + (K2Im + K3Im) / 3.0;

					ARe += CoefA * (DRe + DReTemp);
					AIm += CoefA * (DIm + DImTemp);					
				}
			}

			if (iNumftFreq < nbFtHalf)
			{
				dExpReal = dCoef * exp(ARe + DRe);
				dImag = fmod(AIm + DIm, PI2);

				/* Real Part */
				Density[1][iNumPhiFreq][2*iNumftFreq-1] = dExpReal * cos(dImag);		
				/* Imaginary Part */
				Density[1][iNumPhiFreq][2*iNumftFreq] = dExpReal * sin(dImag);
			}
			else
			{
				DensitySpeq[1][2*iNumPhiFreq-1] = dCoef * exp(ARe + DRe);
				DensitySpeq[1][2*iNumPhiFreq] = 0.0;
				/*
				Density[1][iNumPhiFreq][2*iNumftFreq-1] = 0.0;
				Density[1][iNumPhiFreq][2*iNumftFreq] = 0.0;
				*/
			}

			dftFreq += dftFreqStep;
		}

		dPhifreq += dPhiFreqStep;
		PhiTrans += PhiTransStep;

		if (iNumPhiFreq == nbPhiHalf)
		{
			dPhifreq -= JumpPhi;
			PhiTrans -= JumpPhiTrans;
		}
	}

	/* Save the TF density */
	if (SaveFile)
	{
		LGMSVSaveTFDensity( iNbPhi,
							iNbft,					 
							Density);
	}

	/* Calculation of the density by iFFT */
	rlft3(Density,DensitySpeq,1,iNbPhi,iNbft,-1);
}

//LGMSVOptionPricePrecalc(iNbPhi, iNbft, );

Err LGMSVCalculateTFDensity(	int		iNbPhi,
								int		iNbft,							
								double	dPhiFreqStep,
								double	dftFreqStep,

								double	dExTime,

								double	dLambdaX,
								int		NbSigTime,		/* Term Structure of g(t) */
								double	*SigTime,		
								double	*Sig,
								double	dTStar,			/* Tstar in years from today */						 
								double	dAlpha,
								double	dLambdaEps,
								double	dRho,
								double	PhiMean,

								double	MaxTime,
								int		SaveFile,

								/* Outputs */
								double	***Density,
								double	**DensitySpeq)
{
double	dCoef;
int		iNumPhiFreq, iNumftFreq, nbFtHalf, nbPhiHalf;
double	dftFreq, dPhifreq, dExpReal, dImag, PhiTrans, PhiTransStep, JumpPhi, JumpPhiTrans;
int		i, j, endi, nbp, nbT;
double	sig0;
double	ARe, AIm, DRe, DIm, DReTemp, DImTemp;
double	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im;
double	Coef1ReInit, Coef2ReInit, Coef2ImInit, Coef3ReInit, Coef3ImInit;
double	K1Re, K1Im, K2Re, K2Im, K3Re, K3Im, K4Re, K4Im;
double	dt, t1, t2, dt2;
double	expt, expt2, expt2init, coefexp;
double	CoefA;
double	*Coef1ReT	= NULL,
		*Coef2ReT	= NULL,
		*Coef2ImT	= NULL,
		*Coef3ReT	= NULL,
		*Coef3ImT	= NULL,
		*CoefAT		= NULL,
		*CoefexpT	= NULL;

long	*nbpT		= NULL;
Err		err = NULL;

	/* Normalisation such that rlft3 gives us the true ifft of TF(P(f_phi,f_ft))) */
	dCoef = 2.0 * dftFreqStep * dPhiFreqStep;	/*(iNbPhi*iNbft/2.0); */	

	dftFreqStep *= PI2;
	dPhiFreqStep *= PI2;
	
	nbFtHalf = iNbft / 2 + 1;
	nbPhiHalf = iNbPhi / 2;
	
	PhiTransStep = dPhiFreqStep * PhiMean;
	JumpPhi = iNbPhi * dPhiFreqStep;
	JumpPhiTrans = iNbPhi * PhiTransStep;

	/* Constant evaluation */

	/* First find the volatility index */	
	endi = Get_Index(dExTime, SigTime, NbSigTime);
	nbT = endi + 1;

	Coef1ReT = dvector(0, endi);
	Coef2ReT = dvector(0, endi);
	Coef2ImT = dvector(0, endi);
	Coef3ReT = dvector(0, endi);
	Coef3ImT = dvector(0, endi);
	CoefAT = dvector(0, endi);
	CoefexpT = dvector(0, endi);
	nbpT = lvector(0, endi);

	if (!Coef1ReT || !Coef2ReT || !Coef2ImT || !Coef3ReT || !Coef3ImT || !CoefAT || !CoefexpT || !nbpT)
	{
		err = "Memory allocation faillure in LGMSVCalculateTFDensity";
		goto FREE_RETURN;
	}

	expt2init = exp(dLambdaX * dExTime);
	Coef1ReInit = -0.5 * dAlpha * dAlpha;
	Coef2ReInit = dLambdaEps;
	expt = exp(-dLambdaX * dTStar);
	Coef2ImInit = -dAlpha * dRho * expt;
	Coef3ReInit = 0.5 * expt * expt;
	Coef3ImInit = -1.0 / expt2init / expt2init;	

	for (i=endi; i>=0; i--)
	{
		if (i>0)
		{
			t1 = SigTime[i-1];
		}
		else
		{
			/* First part */
			t1 = 0.0;
		}
		
		if (i==endi || 0==endi)
		{
			/* Last part */
			t2 = dExTime;
		}
		else
		{
			t2 = SigTime[i];
		}

		nbpT[i] = (int) ((t2 - t1) / MaxTime + 0.5);
		dt = (t2 - t1) / nbpT[i];
		dt2 = dt / 2.0;		
		sig0 = Sig[i];

		CoefexpT[i] = exp(-dLambdaX * dt2);
		Coef1ReT[i] = -Coef1ReInit * dt;
		Coef2ReT[i] = -Coef2ReInit * dt;
		Coef2ImT[i] = -Coef2ImInit * sig0 * dt;
		Coef3ReT[i] = -Coef3ReInit * sig0 * sig0 * dt;
		Coef3ImT[i] = -Coef3ImInit * sig0 * sig0 * dt;
		CoefAT[i] = Coef2ReInit * dt2;
	}

	dPhifreq = 0.0;
	PhiTrans = 0.0;

	for (iNumPhiFreq=1; iNumPhiFreq<=iNbPhi; iNumPhiFreq++)
	{					
		dftFreq = 0.0;		

		for (iNumftFreq=1; iNumftFreq<=nbFtHalf; iNumftFreq++)
		{
			/* Initialisation */	
			ARe = 0.0;
			AIm = -PhiTrans;
			DRe = 0.0;
			DIm = 0.0;
			expt2 = expt2init;
			
			for (i=endi; i>=0; i--)
			{					
				/* Solve the ODE between t1 and t2 */												
				coefexp = CoefexpT[i];
				Coef1Re = Coef1ReT[i];
				Coef2Re = Coef2ReT[i];
				Coef2Im = Coef2ImT[i] * dftFreq;
				Coef3Re = Coef3ReT[i] * dftFreq * dftFreq;
				Coef3Im = Coef3ImT[i] * dPhifreq;
				CoefA = CoefAT[i];
				nbp = nbpT[i];

				for (j=0; j<nbp; j++)
				{
					DReTemp = DRe;
					DImTemp = DIm;								
					
					LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
												expt2, DRe, DIm, &K1Re, &K1Im);
					expt2 *= coefexp;
					LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
												expt2, DRe + 0.5 * K1Re, DIm + 0.5 * K1Im,  &K2Re, &K2Im);
					LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
												expt2, DRe + 0.5 * K2Re, DIm + 0.5 * K2Im,  &K3Re, &K3Im);
					expt2 *= coefexp;
					LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
												expt2, DRe + K3Re, DIm + K3Im,  &K4Re, &K4Im);

					DRe += (K1Re + K4Re) / 6.0 + (K2Re + K3Re) / 3.0;
					DIm += (K1Im + K4Im) / 6.0 + (K2Im + K3Im) / 3.0;

					ARe += CoefA * (DRe + DReTemp);
					AIm += CoefA * (DIm + DImTemp);					
				}
			}

			if (iNumftFreq < nbFtHalf)
			{
				dExpReal = dCoef * exp(ARe + DRe);
				dImag = fmod(AIm + DIm, PI2);

				/* Real Part */
				Density[1][iNumPhiFreq][2*iNumftFreq-1] = dExpReal * cos(dImag);		
				/* Imaginary Part */
				Density[1][iNumPhiFreq][2*iNumftFreq] = dExpReal * sin(dImag);
			}
			else
			{
				DensitySpeq[1][2*iNumPhiFreq-1] = dCoef * exp(ARe + DRe);
			}

			dftFreq += dftFreqStep;
		}

		dPhifreq += dPhiFreqStep;
		PhiTrans += PhiTransStep;

		if (iNumPhiFreq == nbPhiHalf)
		{
			dPhifreq -= JumpPhi;
			PhiTrans -= JumpPhiTrans;
		}
	}

	/* Save the TF density */
	if (SaveFile)
	{
		LGMSVSaveTFDensity( iNbPhi,
							iNbft,					 
							Density);
	}

FREE_RETURN:

	if (Coef1ReT) free_dvector(Coef1ReT, 0, endi);
	if (Coef2ReT) free_dvector(Coef2ReT, 0, endi);
	if (Coef2ImT) free_dvector(Coef2ImT, 0, endi);
	if (Coef3ReT) free_dvector(Coef3ReT, 0, endi);
	if (Coef3ImT) free_dvector(Coef3ImT, 0, endi);
	if (CoefAT) free_dvector(CoefAT, 0, endi);
	if (CoefexpT) free_dvector(CoefexpT, 0, endi);
	if (nbpT) free_lvector(nbpT, 0, endi);

	return err;
}

/*	Calculation of the value of A and D */

void LGMSVSolveODELogTSRungeKutta(	double	dPhifreq,
									double	dFtFreq,
									double	T,
									
									double	dLambdaX,
									int		NbSigTime,		/* Term Structure of g(t) */
									double	*SigTime,		
									double	*Sig,
									double	dTStar,			/* Tstar in years from today */						 
									double	dAlpha,
									double	dLambdaEps,
									double	dRho,
									double	PhiMean,

									int		NbPoint,

									double	*LogTSRe,
									double	*LogTSIm)
{
int		i, j, endi, nbp;
double	sig0;
double	ARe, AIm, DRe, DIm, DReTemp, DImTemp;
double	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im;
double	Coef1ReInit, Coef2ReInit, Coef2ImInit, Coef3ReInit, Coef3ImInit;
double	K1Re, K1Im, K2Re, K2Im, K3Re, K3Im, K4Re, K4Im;
double	dt_ref, dt, t1, t2, dt2;
double	expt, expt2, coefexp;
double	CoefA;

	/* First find the volatility index */
	endi = Get_Index(T, SigTime, NbSigTime);

	/* Initialisation */	
	ARe = 0.0;
	AIm = 0.0;
	DRe = 0.0;
	DIm = 0.0;

	expt2 = exp(dLambdaX * T);

	Coef1ReInit = -0.5 * dAlpha * dAlpha;
	Coef2ReInit = dLambdaEps;	
	expt = exp(-dLambdaX * dTStar) * dFtFreq;
	Coef2ImInit = -dAlpha * dRho * expt;
	Coef3ReInit = 0.5 * expt * expt;
	Coef3ImInit = -dPhifreq / expt2 / expt2;

	dt_ref = T / NbPoint;
	
	for (i=endi; i>=0; i--)
	{
		if (i>0)
		{
			t1 = SigTime[i-1];
		}
		else
		{
			/* First part */
			t1 = 0.0;
		}
		
		if (i==endi || endi==0)
		{
			/* Last part */
			t2 = T;
		}
		else
		{
			t2 = SigTime[i];
		}		

		/* Solve the ODE between t1 and t2 */

		nbp = (int) ((t2 - t1) / dt_ref);
		dt = (t2 - t1) / nbp;
		dt2 = dt / 2.0;
		coefexp = exp(-dLambdaX * dt2);
		sig0 = Sig[i];

		Coef1Re = -Coef1ReInit * dt;
		Coef2Re = -Coef2ReInit * dt;
		Coef2Im = -Coef2ImInit * sig0 * dt;
		Coef3Re = -Coef3ReInit * sig0 * sig0 * dt;
		Coef3Im = -Coef3ImInit * sig0 * sig0 * dt;
		CoefA = Coef2ReInit * dt2;

		for (j=0; j<nbp; j++)
		{
			DReTemp = DRe;
			DImTemp = DIm;

			LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
										expt2, DRe, DIm, &K1Re, &K1Im);
			expt2 *= coefexp;
			LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
										expt2, DRe + 0.5 * K1Re, DIm + 0.5 * K1Im,  &K2Re, &K2Im);
			LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
										expt2, DRe + 0.5 * K2Re, DIm + 0.5 * K2Im,  &K3Re, &K3Im);
			expt2 *= coefexp;
			LGMSVSolveODELogTSEvalFunc(	Coef1Re, Coef2Re, Coef2Im, Coef3Re, Coef3Im,
										expt2, DRe + K3Re, DIm + K3Im,  &K4Re, &K4Im);
			DRe += (K1Re + K4Re) / 6.0 + (K2Re + K3Re) / 3.0;
			DIm += (K1Im + K4Im) / 6.0 + (K2Im + K3Im) / 3.0 ;

			ARe += CoefA * (DRe + DReTemp);
			AIm += CoefA * (DIm + DImTemp);
		}
	}

	*LogTSRe = ARe + DRe;
	*LogTSIm = AIm + DIm - dPhifreq * PhiMean;
}




/* -------------------------------------------------------------------------------------------------------------
	LGMSVClosedForm	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVClosedForm2(	
						/* Parameter of diffusion */
						/* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
						double	dLambdaX,
						int		iNbSigTime,		/* Term Structure of g(t) */
						double	*SigTime,		
						double	*Sig,
						double	dTStar,			/* Tstar in years from today */						 
						double	dAlpha,
						double	dLambdaEps,
						double	dRho,

						/* Product description */
						long	lExDate,		/* Exercice date of the swaption  */
						double	dExTime,	/* Exercice of the swaption in years from today */
						int		iNbCoupon,		/* Description of the cashflows */
						double	*CouponTime,
						long	*CouponDate,
						double	*Coupon,
						char	*cYieldCurve,	/* Yield Curve */

						/* Parameter of grids */
						int		iNbPhi,			/* Number of Phi : Should be a power of two */
						int		iNbft,			/* Number of ft : Should be a power of two */	
						double	iNbSigmaPhiGridLeft,
						double	iNbSigmaPhiGridRight,
						double	iNbSigmaftLeft,
						double	iNbSigmaftRight,

						double	iRatioPhi,
						double	iRatioFt,
						int		iPriorityFreqPhi,
						int		iPriorityFreqFt,
						double	iMaxTime,
						int		SaveFile,
						
						/* Outputs */
						double *Price)
{
/* Declaration of locals variables */
Err				err = NULL;
double			lambdaArray[10];

/* For moments calculation */
LGMSVSolFunc	FuncPhi_, *FuncPhi = &FuncPhi_;
LGMSVSolFunc	FuncPhi2_, *FuncPhi2 = &FuncPhi2_;
LGMSVSolFunc	FuncV2_, *FuncV2 = &FuncV2_;
LGMSVSolFunc	FuncPhiV_, *FuncPhiV = &FuncPhiV_;
double			ExpectPhi, ExpectPhi2, ExpevtV2, ExpectPhiV;

double			***Density			= NULL,
				**Densityspeq		= NULL,
				**PayOff			= NULL;

double			dPhitMean, dPhitStd, dPhiStep;
double			dftStd, dftStep;
double			dPhiFreqStep, dftFreqStep;
int				iIndexPhiMean, iIndexft0;

double			LimitPhi, Limitft;
int				i, endi;
double			Coef1ReInit, Coef2ReInit, Coef2ImInit, Coef3ReInit, Coef3ImInit;
double			t1, t2, dt, dt2, expt, expt2init, sig1, AlphaEq;

double			*Coef1ReT	= NULL,
				*Coef2ReT	= NULL,
				*Coef2ImT	= NULL,
				*Coef3ReT	= NULL,
				*Coef3ImT	= NULL,
				*CoefAT		= NULL,
				*CoefexpT	= NULL;

long			*nbpT		= NULL;

clock_t			time1, time2;

	time1 = clock();

	/* Initialisation and memory allocation */	
	PayOff = dmatrix(1, iNbPhi,
					 1, iNbft);

	Density =f3tensor(	1, 1,
						1, iNbPhi,
						1, iNbft);
	
	Densityspeq = dmatrix(	1, 1,
							1, 2 * iNbPhi);

	endi = Get_Index(dExTime, SigTime, iNbSigTime);

	Coef1ReT = dvector(0, endi);
	Coef2ReT = dvector(0, endi);
	Coef2ImT = dvector(0, endi);
	Coef3ReT = dvector(0, endi);
	Coef3ImT = dvector(0, endi);
	CoefAT = dvector(0, endi);
	CoefexpT = dvector(0, endi);
	nbpT = lvector(0, endi);

	/* Gestion of allocation errors */
	if (!Density || !Densityspeq ||
		!Coef1ReT || !Coef2ReT || !Coef2ImT || !Coef3ReT || !Coef3ImT || !CoefAT || !CoefexpT || !nbpT)
	{
		err = "Memory allocation error (1) in LGMSVClosedForm2";
		goto FREE_RETURN;
	}

	/* Precalculations */

	/* Constant independent on time */
	expt2init = exp(dLambdaX * dExTime);
	Coef1ReInit = -0.5 * dAlpha * dAlpha;
	Coef2ReInit = dLambdaEps;
	expt = exp(-dLambdaX * dTStar);
	Coef2ImInit = -dAlpha * dRho * expt;
	Coef3ReInit = 0.5 * expt * expt;
	Coef3ImInit = -1.0 / expt2init / expt2init;

	Limitft = log(iRatioFt);
	LimitPhi = log(iRatioPhi);

	LGMSVMomentInit2(dLambdaX, dAlpha, dLambdaEps, dRho, FuncPhi, FuncPhi2, FuncV2, FuncPhiV);
	ExpectPhi = 0.0;
	ExpectPhi2 = 0.0;
	ExpevtV2 = 1.0;
	ExpectPhiV = 0.0;
	
	t1 = 0.0;
	for (i=0; i<=endi; i++)
	{
		/* Precalculation on the option i*/

		/* First the time discretisation */			
		if (i < endi)
		{
			t2 = SigTime[i];
		}
		else
		{
			t2 = dExTime;
		}

		nbpT[i] = (int) ((t2 - t1) / iMaxTime + 0.5);
		dt = (t2 - t1) / nbpT[i];
		dt2 = dt / 2.0;

		

		AlphaEq = dAlpha / 2.0;
		if (dLambdaEps > 1.0E-16)
		{
			AlphaEq *= sqrt((1.0 - exp(-dLambdaEps * t2)) / (dLambdaEps * t2));
		}		

		/* Calculate constant values */
		sig1 = Sig[i];
		CoefexpT[i] = exp(-dLambdaX * dt2);
		Coef1ReT[i] = -Coef1ReInit * dt;
		Coef2ReT[i] = -Coef2ReInit * dt;
		Coef2ImT[i] = -Coef2ImInit * sig1 * dt;
		Coef3ReT[i] = -Coef3ReInit * sig1 * sig1 * dt;
		Coef3ImT[i] = -Coef3ImInit * sig1 * sig1 * dt;
		CoefAT[i] = Coef2ReInit * dt2;

		LGMSVMomentCalculation2(sig1 * sig1, t2 - t1, ExpectPhi, ExpectPhi2, ExpevtV2, ExpectPhiV,
			FuncPhi, FuncPhi2, FuncV2, FuncPhiV, lambdaArray, &ExpectPhi, &ExpectPhi2, &ExpevtV2, &ExpectPhiV);

		t1 = t2;
	}

	dPhitMean = ExpectPhi;
	dPhitStd = sqrt(ExpectPhi2 - ExpectPhi * ExpectPhi);
	dftStd = sqrt(ExpectPhi) * exp(-dLambdaX * (dTStar - dExTime));

	/* Find Frequence */		
	LGMSVFindFreq(	iNbPhi, iNbft, AlphaEq, dRho, endi, expt2init,
					Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT, CoefexpT, nbpT,
					dPhitMean, dPhitStd, dftStd, iNbSigmaPhiGridLeft, iNbSigmaPhiGridRight, iNbSigmaftLeft, iNbSigmaftRight,
					LimitPhi, Limitft, iPriorityFreqPhi, iPriorityFreqFt,
					&dPhiFreqStep, &dPhiStep, &iIndexPhiMean, &dftFreqStep, &dftStep, &iIndexft0);

	/* First pricing */
	LGMSVOptionPricePrecalc(iNbPhi, iNbft, dLambdaX, dTStar, endi, expt2init,
							Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT, CoefexpT, nbpT,
							dPhitMean, SaveFile, Density, Densityspeq, PayOff,
							dExTime, iNbCoupon, CouponTime, Coupon, 
							dPhiFreqStep, dPhiStep, iIndexPhiMean, dftFreqStep, dftStep, iIndexft0, Price);
		
	
	/* Save the Payoff and Density */
	if (SaveFile)
	{
		LGMSVSaveDensity(iNbPhi,
						 iNbft,
						 iIndexPhiMean,						 
						 iIndexft0,
						 dPhiStep,
						 dftStep,
						 dPhitMean,
						 dPhitStd,
						 Density);

		LGMSVSavePayOff(iNbPhi,
						iNbft,
						iIndexPhiMean,
						iIndexft0,
						dPhiStep,
						dftStep,
						dPhitMean,
						PayOff);
	}
	

	time2 = clock();
	
FREE_RETURN:

	/* free memory */
	if (Density)
		free_f3tensor(Density, 1, 1,1, iNbPhi,1, iNbft);

	if (Densityspeq)
		free_dmatrix(Densityspeq, 1, 1,1, iNbPhi);

	if (PayOff)
		free_dmatrix(PayOff, 1, iNbPhi, 1,iNbft);

	if (Coef1ReT) free_dvector(Coef1ReT, 0, endi);
	if (Coef2ReT) free_dvector(Coef2ReT, 0, endi);
	if (Coef2ImT) free_dvector(Coef2ImT, 0, endi);
	if (Coef3ReT) free_dvector(Coef3ReT, 0, endi);
	if (Coef3ImT) free_dvector(Coef3ImT, 0, endi);
	if (CoefAT) free_dvector(CoefAT, 0, endi);
	if (CoefexpT) free_dvector(CoefexpT, 0, endi);
	if (nbpT) free_lvector(nbpT, 0, endi);
}

Err LGMSVOption(
				char			*und_name,						/*	Name of the underlying */
				char			*yc_name,						/*	Name of the yield curve */
				char			*ref_rate_name,					/*	Name of the reference rate */
				char			*swaption_freq,					/*	Frequency and basis of underlying swaptions */
				char			*swaption_basis,

				long			lExDate,
				long			lStartDate,
				long			lEndDate,
				double			dStrike,

				double			iMaxTime,
				double			iRatioPhi,
				double			iRatioFt,
				int				iPriorityFreqPhi,
				int				iPriorityFreqFt,
				int				iSaveFile,
				long			iNbPhi,
				long			iNbft,

				double			iNbSigmaPhiGridLeft,
				double			iNbSigmaPhiGridRight,
				double			iNbSigmaXGridLeft,
				double			iNbSigmaXGridRight,

				/* Output */
				double			*pSwaptionPrice)
{
	int				i,   ncpn;
	SrtCompounding	ifreq;
	SrtBasisCode	ibasis;
	long			cpn_date[MAX_CPN];
	double			cpn_time[MAX_CPN],
					cpn_cvg[MAX_CPN];
				
	long			theo_date, act_date, temp_date;
	long			today; 

	double			lExTime;
	double			swp_rte, swp_cash, level;

	SrtCurvePtr		yc_ptr;
	Err				err				= NULL;

	int				dNbSig;
	double			*dSigTime		= NULL,
					*dSig			= NULL,
					*dAlphaTS		= NULL,
					*dRhoTS			= NULL,
					*dLambdaEpsTS	= NULL;

	double			dTau;
	double			dAlpha;
	double			dRho;
	double			dLambdaEps;
	double			dTStar;
	int				one2F;

	yc_ptr = lookup_curve (yc_name);
	if (!yc_ptr)
	{
		err = "Yield Curve not found";
		goto FREE_RETURN;
	}

	today = get_today_from_curve (yc_ptr);	

	/*	1.)	Setup the bond schedule and its coupons */
	
	/*	Coupons */

	err = interp_compounding (swaption_freq, &ifreq);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = interp_basis (swaption_basis, &ibasis);
	if (err)
	{
		goto FREE_RETURN;
	}

	theo_date = lEndDate;
	act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
	ncpn = 1;

	while (act_date >= lStartDate)
	{
		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
		act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		ncpn++;
	}
	ncpn--;

	if (ncpn < 2)
	{
		err = "Not enough coupons";
		goto FREE_RETURN;		
	}

	theo_date = lEndDate;
	act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
	i = ncpn - 1;

	while (i >= 0)
	{
		cpn_time[i] = (act_date - today) * YEARS_IN_DAY; 
		cpn_date[i] = act_date;

		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

		temp_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		cpn_cvg[i] = coverage (temp_date, act_date, ibasis);
		act_date = temp_date;

		i--;
	}
	
	/* Calculate the Basis Spread */
	err = swp_f_ForwardRate(
							cpn_date[0],
							cpn_date[ncpn-1],
							swaption_freq,
							swaption_basis,
							yc_name,
							ref_rate_name,
							&swp_rte);
	
	level = 0.0;
	for (i=1; i<ncpn; i++)
	{
		level += cpn_cvg[i] * swp_f_df(today, cpn_date[i], yc_name);
	}
	swp_cash = (swp_f_df(today, cpn_date[0], yc_name) - swp_f_df(today, cpn_date[ncpn-1], yc_name)) / level;

	dStrike -= swp_rte - swp_cash;

	/* +S*LVL = (if S=SwapCash Rate) Df(Start)-Df(end) */	
	/* -K*Lvl */

	cpn_cvg[0] = 1.0 * swp_f_df(today, cpn_date[0], yc_name);	

	for (i=1; i<ncpn; i++)
	{
		cpn_cvg[i] *= -dStrike * swp_f_df(today, cpn_date[i], yc_name);
	}

	cpn_cvg[ncpn-1] -= 1.0 * swp_f_df(today, cpn_date[ncpn-1], yc_name);
	lExTime = (lExDate - today) * YEARS_IN_DAY;

	/* Get the TS */		

	err = Get_LGMSV_TermStructure(	und_name,
									&dSigTime,
									&dSig,
									&dAlphaTS,
									&dRhoTS,
									&dLambdaEpsTS,
									&dTStar,
									&dNbSig,
									&dTau,
									&one2F,
									NULL,
									NULL,
									NULL);

	if (err)
	{
		goto FREE_RETURN;
	}

	/* for now no TS */

	dAlpha = dAlphaTS[0];
	dRho = dRhoTS[0],
	dLambdaEps = dLambdaEpsTS[0];

	for (i=1; i<dNbSig; i++)
	{
		if ((fabs(dAlphaTS[i] - dAlpha) > 1.0E-08) || (fabs(dRhoTS[i] - dRho) > 1.0E-08) || (fabs(dLambdaEpsTS[i] - dLambdaEps) > 1.0E-08))
		{
			err = "no alpha, rho, lameps TS allowed for now";
			goto FREE_RETURN;
		}
	}
	
	LGMSVClosedForm2(
					1.0 / dTau,
					dNbSig,		/* Term Structure of g(t) */
					dSigTime,		
					dSig,
					dTStar,			/* Tstar in years from today */
					dAlpha,		/* Alpha of V = Eps^2 */
					dLambdaEps,	/* LambdaEps of V = Eps^2 */
					dRho,
					lExDate,							/* Exercice date of the swaption  */
					lExTime,		/* Exercice of the swaption in years from today */					 
					ncpn,		/* Description of the cashflows */
					cpn_time,
					cpn_date,
					cpn_cvg,
					yc_name,	/* Yield Curve */					 
					iNbPhi,			/* Number of Phi : Should be a power of two */
					iNbft,			/* Number of ft : Should be a power of two */	
					iNbSigmaPhiGridLeft,
					iNbSigmaPhiGridRight,
					iNbSigmaXGridLeft,
					iNbSigmaXGridRight,
					iRatioPhi,
					iRatioFt,
					iPriorityFreqPhi,
					iPriorityFreqFt,
					iMaxTime,
					iSaveFile,
					pSwaptionPrice);

FREE_RETURN:

	if (dSigTime) free (dSigTime);
	if (dSig) free (dSig);
	if (dAlphaTS) free (dAlphaTS);
	if (dRhoTS) free (dRhoTS);
	if (dLambdaEpsTS) free (dLambdaEpsTS);

	return err;
}


/* -------------------------------------------------------------------------------------------------------------
	LGMSVClosedForm	

  -------------------------------------------------------------------------------------------------------------- */
void LGMSVCalib2(	
						/* Parameter of diffusion */
						/* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
						double	dLambdaX,
						double	dTStar,
						double	dAlpha,
						double	dLambdaEps,
						double	dRho,

						/* Product description */
						char	*cYieldCurve,
						double	*dExTime,
						int		nbEx,
						int		*iNbCoupon,
						double	**CouponTime,
						double	**Coupon,
						double	*TargetPrice,
						double	*TargetVol,

						/*	Additional informations to estimate Vega */
						double	*Fwd,
						double	*Strike,
						double	*Level,
						
						/* Parameter of grids */
						int		iNbPhi,			/* Number of Phi : Should be a power of two */
						int		iNbft,			/* Number of ft : Should be a power of two */	
						double	iNbSigmaPhiGridLeft,
						double	iNbSigmaPhiGridRight,
						double	iNbSigmaftLeft,
						double	iNbSigmaftRight,

						double	iRatioPhi,
						double	iRatioFt,
						int		iPriorityFreqPhi,
						int		iPriorityFreqFt,
						double	iMaxTime,
						
						/* Newton parameters */
						double	Precision,
						int		NbIterMax,

						/* Outputs */
						double	*sigma)
{
/* Declaration of locals variables */
Err				err = NULL;
double			lambdaArray[10];

/* For moments calculation */
LGMSVSolFunc	FuncPhi_, *FuncPhi = &FuncPhi_;
LGMSVSolFunc	FuncPhi2_, *FuncPhi2 = &FuncPhi2_;
LGMSVSolFunc	FuncV2_, *FuncV2 = &FuncV2_;
LGMSVSolFunc	FuncPhiV_, *FuncPhiV = &FuncPhiV_;
double			ExpectPhi, ExpectPhi2, ExpevtV2, ExpectPhiV;
double			ExpectPhiTemp, ExpectPhi2Temp, ExpevtV2Temp, ExpectPhiVTemp;
double			CoefStd;

double			***Density			= NULL,
				**Densityspeq		= NULL,
				**PayOff			= NULL;

double			*Coef1ReT	= NULL,
				*Coef2ReT	= NULL,
				*Coef2ImT	= NULL,
				*Coef3ReT	= NULL,
				*Coef3ImT	= NULL,
				*CoefAT		= NULL,
				*CoefexpT	= NULL,
				**res_iter	= NULL;

long			*nbpT		= NULL;

double			Coef1ReInit, Coef2ReInit, Coef2ImInit, Coef3ReInit, Coef3ImInit;
double			dPhitMean, dPhitStd, dPhiStep;
double			dftStd, dftStep;
double			dPhiFreqStep, dftFreqStep;
double			expt, t1, t2, dt, dt2, expt2init;
double			UpdateRatio;
int				iIndexPhiMean, iIndexft0;

double			LimitPhi, Limitft;
double			AlphaEq;
double			var, newvar, vega, imp_vol, shift_price, beta_vol;
double			sig1, price1, sig2, price2, factor;
int				i, j, k, l;
clock_t			time1, time2;

	time1 = clock();

	/* Initialisation and memory allocation */
	PayOff = dmatrix(1, iNbPhi,
					 1, iNbft);

	Density =f3tensor(	1, 1,
						1, iNbPhi,
						1, iNbft);
	
	Densityspeq = dmatrix(	1, 1,
							1, 2 * iNbPhi);

	Coef1ReT = dvector(0, nbEx-1);
	Coef2ReT = dvector(0, nbEx-1);
	Coef2ImT = dvector(0, nbEx-1);
	Coef3ReT = dvector(0, nbEx-1);
	Coef3ImT = dvector(0, nbEx-1);
	CoefAT = dvector(0, nbEx-1);
	CoefexpT = dvector(0, nbEx-1);
	nbpT = lvector(0, nbEx-1);

	res_iter = dmatrix(0, NbIterMax, 0, 1);	

	/* Gestion of allocation errors */
	if (!Density || !Densityspeq ||
		!Coef1ReT || !Coef2ReT || !Coef2ImT || !Coef3ReT || !Coef3ImT || !CoefAT || !CoefexpT || !nbpT || !res_iter)
	{
		err = "Memory allocation error (1) in LGMSVClosedForm2";
		goto FREE_RETURN;
	}

	/* Constant independent on time */
	Coef1ReInit = -0.5 * dAlpha * dAlpha;
	Coef2ReInit = dLambdaEps;
	expt = exp(-dLambdaX * dTStar);
	Coef2ImInit = -dAlpha * dRho * expt;
	Coef3ReInit = 0.5 * expt * expt;

	LimitPhi = log(iRatioPhi);
	Limitft = log(iRatioFt);

	LGMSVMomentInit2(dLambdaX, dAlpha, dLambdaEps, dRho, FuncPhi, FuncPhi2, FuncV2, FuncPhiV);
	ExpectPhi = 0.0;
	ExpectPhi2 = 0.0;
	ExpevtV2 = 1.0;
	ExpectPhiV = 0.0;

	t1 = 0.0;
	var = 0.0;	

	for (i=0; i<nbEx; i++)
	{
		/* Precalculation on the option i*/

		/* First the time discretisation */			
		t2 = dExTime[i];
		nbpT[i] = (int) ((t2 - t1) / iMaxTime + 0.5);
		dt = (t2 - t1) / nbpT[i];
		dt2 = dt / 2.0;		

		k = 1;

		/* Constant calculation */
		expt2init = exp(dLambdaX * t2);
		Coef3ImInit = -1.0 / expt2init / expt2init;
		CoefStd = expt * expt2init;
		
		AlphaEq = dAlpha / 2.0;

		if (dLambdaEps > 1.0E-16)
		{			
			AlphaEq *= sqrt((1.0 - exp(-dLambdaEps * t2)) / (dLambdaEps * t2));
		}		

		/* Update constant values already calculated */		
		UpdateRatio = exp(-2.0 * dLambdaX * (t2 - t1));
		for (j=0; j<i; j++)
		{
			Coef3ImT[j] *= UpdateRatio;
		}

		/* First guess */		
		sig1 = sigma[i];		
		newvar = var + sig1 * sig1 * (t2 - t1);

		/* Moment evaluation */		
		LGMSVMomentCalculation2(sig1 * sig1, t2 - t1, ExpectPhi, ExpectPhi2, ExpevtV2, ExpectPhiV,
			FuncPhi, FuncPhi2, FuncV2, FuncPhiV, lambdaArray, &ExpectPhiTemp, &ExpectPhi2Temp, &ExpevtV2Temp, &ExpectPhiVTemp);

		dPhitMean = ExpectPhiTemp;
		dPhitStd = sqrt(ExpectPhi2Temp - ExpectPhiTemp * ExpectPhiTemp);
		dftStd = sqrt(ExpectPhiTemp) * CoefStd;

		/* Calculate constant values */
		CoefexpT[i] = exp(-dLambdaX * dt2);
		Coef1ReT[i] = -Coef1ReInit * dt;
		Coef2ReT[i] = -Coef2ReInit * dt;
		Coef2ImT[i] = -Coef2ImInit * sig1 * dt;
		Coef3ReT[i] = -Coef3ReInit * sig1 * sig1 * dt;
		Coef3ImT[i] = -Coef3ImInit * sig1 * sig1 * dt;
		CoefAT[i] = Coef2ReInit * dt2;

		/* Find Frequence */		
		LGMSVFindFreq(iNbPhi, iNbft, AlphaEq, dRho, i, expt2init,
								Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT, CoefexpT, nbpT,
								dPhitMean, dPhitStd, dftStd, iNbSigmaPhiGridLeft, iNbSigmaPhiGridRight, iNbSigmaftLeft, iNbSigmaftRight,
								LimitPhi, Limitft, iPriorityFreqPhi, iPriorityFreqFt,
								&dPhiFreqStep, &dPhiStep, &iIndexPhiMean, &dftFreqStep, &dftStep, &iIndexft0);

		/* First pricing */
		LGMSVOptionPricePrecalc(iNbPhi, iNbft, dLambdaX, dTStar, i, expt2init,
								Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT, CoefexpT, nbpT,
								dPhitMean, 0, Density, Densityspeq, PayOff,
								dExTime[i], iNbCoupon[i], CouponTime[i], Coupon[i], 
								dPhiFreqStep, dPhiStep, iIndexPhiMean, dftFreqStep, dftStep, iIndexft0, &price1);

		if (price1 < 0.0 || price1 > 10000)
		{
			smessage("Numerical problem at option %d, try to decrease MaxTime or number of Phi/ft", i+1);	
			for (j=i; j<nbEx; j++)
			{
				sigma[j] = 0.0;
			}
			goto FREE_RETURN;
		}
		
		if (fabs(TargetPrice[i] - price1) > Precision)
		{
			/* We do a Newton */
			res_iter[0][0] = sig1;
			res_iter[0][1] = price1;

			/* New guess */

			/* Try to estimate vega using SABR */

			err = srt_f_optimpvol(	price1,
									Fwd[i],
									Strike[i],
									t2,
									Level[i],
									SRT_CALL,
									SRT_NORMAL,
									&imp_vol);

			if (err)
			{				
				goto FREE_RETURN;
			}

			vol_conv(imp_vol,
					 SABR_STR_NORM,
					 &beta_vol,
					 SABR_ATM_BETA,
					 Fwd[i],
					 Strike[i],
					 t2,
					 AlphaEq,
					 0.0,
					 dRho);

			if (TargetPrice[i] > price1)
			{
				factor = 1;
			}
			else
			{
				factor = -1;
			}

			err = op_sabr_pricing(Fwd[i], Strike[i], t2, Level[i], SRT_CALL, AlphaEq, 0.0, dRho, beta_vol * (1.0 + factor * LGM_VEGA_SHIFT), SABR_ATM_BETA, &shift_price);
			
			if (err)
			{
				goto FREE_RETURN;
			}
									
			vega = (shift_price - price1) / (factor * LGM_VEGA_SHIFT);
			
			sig2 = sqrt((newvar * pow(1.0 + (TargetPrice[i] - price1) / vega, 2) - var) / (t2 - t1));
				
			/*
			sig2 = sqrt((newvar * TargetPrice[i] * TargetPrice[i] / price1 / price1 - var) / (t2 - t1));
			*/						

			while (fabs(TargetPrice[i] - price1) > Precision && k < NbIterMax)
			{
				k++;

				/* New caluclation pour sig2 */
				sigma[i] = sig2;

				/* Moment evaluation */
				LGMSVMomentCalculation2(sig2 * sig2, t2 - t1, ExpectPhi, ExpectPhi2, ExpevtV2, ExpectPhiV,
				FuncPhi, FuncPhi2, FuncV2, FuncPhiV, lambdaArray, &ExpectPhiTemp, &ExpectPhi2Temp, &ExpevtV2Temp, &ExpectPhiVTemp);

				dPhitMean = ExpectPhiTemp;
				dPhitStd = sqrt(ExpectPhi2Temp - ExpectPhiTemp * ExpectPhiTemp);
				dftStd = sqrt(ExpectPhiTemp) * CoefStd;

				/* Calculate constant values */
				UpdateRatio = sig2 / sig1;
				Coef2ImT[i] *= UpdateRatio;
				Coef3ReT[i] *= UpdateRatio * UpdateRatio;
				Coef3ImT[i] *= UpdateRatio * UpdateRatio;

				/* Find Frequence */		
				LGMSVFindFreq(iNbPhi, iNbft, AlphaEq, dRho, i, expt2init,
										Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT, CoefexpT, nbpT,
										dPhitMean, dPhitStd, dftStd, iNbSigmaPhiGridLeft, iNbSigmaPhiGridRight, iNbSigmaftLeft, iNbSigmaftRight,
										LimitPhi, Limitft, iPriorityFreqPhi, iPriorityFreqFt,
										&dPhiFreqStep, &dPhiStep, &iIndexPhiMean, &dftFreqStep, &dftStep, &iIndexft0);

				/* First pricing */
				LGMSVOptionPricePrecalc(iNbPhi, iNbft, dLambdaX, dTStar, i, expt2init,
										Coef1ReT, Coef2ReT, Coef2ImT, Coef3ReT, Coef3ImT, CoefAT, CoefexpT, nbpT,
										dPhitMean, 0, Density, Densityspeq, PayOff,
										dExTime[i], iNbCoupon[i], CouponTime[i], Coupon[i], 
										dPhiFreqStep, dPhiStep, iIndexPhiMean, dftFreqStep, dftStep, iIndexft0, &price2);
				
				l = 0;
				while (l < k-1 && res_iter[l][0] < sig2)
				{
					l++;
				}

				if (l < k-1)
				{
					for (j=k-2; j>=l; j--)
					{
						res_iter[j+1][0] = res_iter[j][0];
						res_iter[j+1][1] = res_iter[j][1];
					}
				}

				res_iter[l][0] = sig2;
				res_iter[l][1] = price2;

				if (fabs(sig2 - sig1) > 1.0E-10)
				{
					vega = (price2 - price1) / (sig2 - sig1);
				}
				else
				{
					vega = -10;
				}

				if (vega < 0.0 || price2 < 0.0 || price2 > 10000)
				{
					smessage("Numerical problem at option %d, try to increase MaxTime or number of Phi/ft", i+1);	
					for (j=i; j<nbEx; j++)
					{
						sigma[j] = 0.0;
					}
					goto FREE_RETURN;
				}

				sig1 = sig2;
				price1 = price2;

				sig2 = solve_for_next_coef(	res_iter,
											k,
											TargetPrice[i],
											2);
			}

			LGMSVMomentCalculation2(sig2 * sig2, t2 - t1, ExpectPhi, ExpectPhi2, ExpevtV2, ExpectPhiV,
				FuncPhi, FuncPhi2, FuncV2, FuncPhiV, lambdaArray, &ExpectPhi, &ExpectPhi2, &ExpevtV2, &ExpectPhiV);

			/* Update constant values */
			UpdateRatio = sig2 / sig1;
			Coef2ImT[i] *= UpdateRatio;
			Coef3ReT[i] *= UpdateRatio * UpdateRatio;
			Coef3ImT[i] *= UpdateRatio * UpdateRatio;
		}
		else
		{
			sig2 = sig1;
			ExpectPhi = ExpectPhiTemp;
			ExpectPhi2 = ExpectPhi2Temp;
			ExpevtV2 = ExpevtV2Temp;
			ExpectPhiV = ExpectPhiVTemp;
		}

		sigma[i] = sig2;
		var += sig2 * sig2 * (t2 - t1);		

		t1 = t2;

		if (fabs(TargetPrice[i] - price1) < Precision)
		{
			smessage("Option %d: success at iteration %d", i+1, k);
		}
		else
		{
			smessage("Option %d: may have failed", i+1);
		}

	}

	time2 = clock();

	smessage("Calibration time :%.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);


FREE_RETURN:

	/* free memory */
	if (Density)
		free_f3tensor(Density, 1, 1,1, iNbPhi,1, iNbft);

	if (Densityspeq)
		free_dmatrix(Densityspeq, 1, 1,1, iNbPhi);

	if (PayOff)
		free_dmatrix(PayOff, 1, iNbPhi, 1,iNbft);

	if (Coef1ReT) free_dvector(Coef1ReT, 0, nbEx-1);
	if (Coef2ReT) free_dvector(Coef2ReT, 0, nbEx-1);
	if (Coef2ImT) free_dvector(Coef2ImT, 0, nbEx-1);
	if (Coef3ReT) free_dvector(Coef3ReT, 0, nbEx-1);
	if (Coef3ImT) free_dvector(Coef3ImT, 0, nbEx-1);
	if (CoefAT) free_dvector(CoefAT, 0, nbEx-1);
	if (CoefexpT) free_dvector(CoefexpT, 0, nbEx-1);
	if (nbpT) free_lvector(nbpT, 0, nbEx-1);

	if (res_iter) free_dmatrix(res_iter, 0, NbIterMax, 0, 1);	
}

static Err get_end_date(
long				ex_date,
long				struct_end_date,
char				*tenor,
int					theo_act,
long				*end_date)
{
	long			start_date;
	Err				err;

	start_date = add_unit (ex_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	strupper (tenor);
	strip_white_space (tenor);
	if (*tenor == 'D')
	{
		*end_date = struct_end_date;
	}
	else
	{
		err = add_tenor (start_date, tenor, NO_BUSDAY_CONVENTION, end_date);
		if (err)
		{
			return err;
		}
	}

	if (theo_act)
	{
		*end_date = bus_date_method (*end_date, MODIFIED_SUCCEEDING);
	}

	return NULL;
}

Err cpd_calibSV_diagonal_2(
	/*	Market */
	char			*yc_name,						/*	Name of the yield curve */
	char			*vol_curve_name,				/*	Name of the market vol curve */
	char			*ref_rate_name,					/*	Name of the reference rate */
	Err				(*get_cash_vol)(				/*	Function to get cash vol from the market */
						char	*vol_curve_name,	
						double	start_date, 
						double	end_date,
						double	cash_strike,
						int		zero,
						char	*ref_rate_name,
						double	*vol,
						double	*power),
	char			*instr_freq,					/*	Frequency and basis of instruments */
	char			*instr_basis,
						/*	If ex_date is NULL,
						exercise dates will be generated 2bd before start */
	/*	Structure */
	int				num_ex_dates,					/*	Exercise dates, 
															all supposed to be on or after today */
	long			*ex_date,						/*	Supposed to be sorted */
	char			**end_tenor,					/*	Tenors of the underlying instruments
															or "DIAG" */
	long			end_date,						/*	End date for diagonal */
	double			*strike,						/*	Strikes 
														0: ATM */
	/*	Model */
	double			dLambdaX,
	double			dAlpha,
	double			dLambdaEps,
	double			dRho,
	double			dTStar,

	/* Parameter of grids */
	int				iNbPhi,			/* Number of Phi : Should be a power of two */
	int				iNbft,			/* Number of ft : Should be a power of two */	
	double			iNbSigmaPhiGridLeft,
	double			iNbSigmaPhiGridRight,
	double			iNbSigmaftLeft,
	double			iNbSigmaftRight,
					
	double			iRatioPhi,
	double			iRatioFt,
	int				iPriorityFreqPhi,
	int				iPriorityFreqFt,
	double			iMaxTime,
	
	/* Newton parameters */
	double			Precision,
	int				NbIterMax,
	/*	Output */
	int				*num_sig,						/*	Answer */
	double			**sig_time,
	double			**sig,
	/*	Parameters */
	CPD_DIAG_CALIB_PARAM param,
	/*	Calibration instrument data */
	CPD_CALIB_INST_DATA	inst_data)					/*	NULL = don't save calibration instrument data */
{
	int				i, j, k, l, nex, ncpn;
	SrtCompounding	ifreq;
	SrtBasisCode	ibasis;
	double			ex_time[MAX_CPN],
					ex_lfwd[MAX_CPN],
					ex_llvl[MAX_CPN],
					ex_lstrike[MAX_CPN],
					ex_lvol[MAX_CPN],
					ex_lprice[MAX_CPN],
					atm_price[MAX_CPN],
					ex_zeta[MAX_CPN];
	int				ex_cpn[MAX_CPN],
					ex_endcpn[MAX_CPN];
	long			cpn_date[MAX_CPN];
	double			cpn_time[MAX_CPN],
					cpn_cvg[MAX_CPN],
					cpn_df[MAX_CPN];
	long			tmplng1[MAX_CPN],
					tmplng2[MAX_CPN];
	long			*theo_end_dates,
					*act_end_dates;
	long			theo_date, act_date, temp_date, temp_date2;
	long			today; 
	double			lvl, dfi, dff;
	double			power;
	double			swp_rte, spr;
	double			std;
	SrtCurvePtr		yc_ptr;
	double			TStar_date;
	double			AlphaEq, atm_vol;
	int				*iNbCoupon;
	double			**CouponTime;
	double			**Coupon;	
	Err				err				= NULL;

	theo_end_dates = &(tmplng1[0]);
	act_end_dates = &(tmplng2[0]);

	*sig_time = NULL;
	*sig = NULL;

	yc_ptr = lookup_curve (yc_name);
	if (!yc_ptr)
	{
		err = "Yield Curve not found";
		goto FREE_RETURN;
	}
	today = get_today_from_curve (yc_ptr);
	

	/*	1.)	Setup the bond schedule and its coupons */
	
	/*	Coupons */

	err = interp_compounding (instr_freq, &ifreq);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = interp_basis (instr_basis, &ibasis);
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Find the end date as the longest total maturity */
	theo_date = end_date;
	act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
	for (i=0; i<num_ex_dates; i++)
	{
		err = get_end_date (ex_date[i], end_date, end_tenor[i], 0, &(theo_end_dates[i]));
		if (err)
		{
			goto FREE_RETURN;
		}
		act_end_dates[i] = bus_date_method (theo_end_dates[i], MODIFIED_SUCCEEDING);
	}
	for (i=0; i<num_ex_dates; i++)
	{
		if (theo_end_dates[i] > theo_date || act_end_dates[i] > act_date)
		{
			theo_date = theo_end_dates[i];
			act_date = act_end_dates[i];
		}
	}
	ncpn = 1;
	temp_date = theo_date;
	temp_date2 = act_date;

	while (act_date > today)
	{
		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
		act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		ncpn++;
	}
	ncpn--;

	if (ncpn < 2)
	{
		err = "Not enough coupons in cpd_calib_diagonal";
		goto FREE_RETURN;		
	}

	theo_date = temp_date;
	act_date = temp_date2;
	i = ncpn - 1;

	while (i >= 0)
	{
		cpn_time[i] = (act_date - today) * YEARS_IN_DAY; 
		cpn_date[i] = act_date;

		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

		temp_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		cpn_cvg[i] = coverage (temp_date, act_date, ibasis);
		cpn_df[i] = swp_f_df (today, act_date, yc_name);
		act_date = temp_date;

		i--;
	}
	cpn_cvg[0] = 0.0;

	/*	Exercise */

	/*	Remove past dates */
	while (ex_date[0] <= today)
	{
		ex_date++;
		end_tenor++;
		theo_end_dates++;
		act_end_dates++;
		strike++;
		num_ex_dates--;
		if (num_ex_dates == 0)
		{
			"All exercise dates are past in cpd_calib_diagonal";
			goto FREE_RETURN;
		}
	}

	/*	Remove redundant dates */
	j = ncpn - 1;
	l = ncpn + 1;
	for (i=num_ex_dates-1; i>=0; i--)
	{
		while (j > 0 && cpn_date[j] > ex_date[i])
		{
			j--;
		}
		if (cpn_date[j] < ex_date[i])
		{
			j++;
		}
		
		if (j >= ncpn - 1 || j == l)
		{
			for (k=i-1; k>=0; k--)
			{
				ex_date[k+1] = ex_date[k];
				strcpy (end_tenor[k+1], end_tenor[k]);
				theo_end_dates[k+1] = theo_end_dates[k];
				act_end_dates[k+1] = act_end_dates[k];
				strike[k+1] = strike[k];
			}

			ex_date++;
			end_tenor++;
			theo_end_dates++;
			act_end_dates++;
			strike ++;
			num_ex_dates--;			
			if (num_ex_dates < 1)
			{
				err = "All exercise dates are past in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
		else
		{
			l = j;
		}
	}

	/*	Remove close dates */
	j = num_ex_dates - 1;
	for (i=num_ex_dates-2; i>=0; i--)
	{
		if ((ex_date[j] - ex_date[i]) * YEARS_IN_DAY < param->min_time - ONE_MONTH)
		{
			for (k=i-1; k>=0; k--)
			{
				ex_date[k+1] = ex_date[k];
				strcpy (end_tenor[k+1], end_tenor[k]);
				theo_end_dates[k+1] = theo_end_dates[k];
				act_end_dates[k+1] = act_end_dates[k];
				strike[k+1] = strike[k];
			}

			ex_date++;
			end_tenor++;
			theo_end_dates++;
			act_end_dates++;
			strike ++;
			num_ex_dates--;	
			j--;
			if (num_ex_dates < 1)
			{
				err = "All exercise dates are past in cpd_calib_diagonal";
				goto FREE_RETURN;
			}
		}
		else
		{
			j = i;
		}
	}

	/*	Remove last? */
	if (param->skip_last && num_ex_dates > 1)
	{
		num_ex_dates--;
	}

	nex = num_ex_dates;
	j = 0;		
	for (i=0; i<nex; i++)
	{
		while (cpn_date[j] < ex_date[i])
		{
			j++;
		}

		ex_cpn[i] = j;
		ex_time[i] = (ex_date[i] - today) * YEARS_IN_DAY;

		k = j;
		while (cpn_date[k] < act_end_dates[i])
		{
			k++;
		}
		if (k > 0 && cpn_date[k] - act_end_dates[i] > act_end_dates[i] - cpn_date[k-1])
		{
			k--;
		}

		if (k <= j)
		{
			k = j + 1;
		}
		ex_endcpn[i] = k;
		
		if (j >= ncpn || k >= ncpn)
		{
			err = "Coupon date bug in cpd_calib_diagonal";
			goto FREE_RETURN;
		}
	}

	/*	Underlyings */

	/*	Long */

	for (i=0; i<nex; i++)
	{
		j = ex_cpn[i];
		l = ex_endcpn[i];

		lvl = 0.0;
		for (k=j+1; k<=l; k++)
		{
			lvl += cpn_cvg[k] * cpn_df[k];
		}
		dfi = swp_f_df (today, cpn_date[j], yc_name);
		dff = swp_f_df (today, cpn_date[l], yc_name);

		ex_llvl[i] = lvl;
		ex_lfwd[i] = (dfi - dff) / lvl;

		/*	ATM std */
		err = get_cash_vol (
			vol_curve_name, 
			add_unit (ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
			theo_end_dates[i], 
			ex_lfwd[i], 
			0, 
			ref_rate_name, 
			&std, 
			&power);
		if (err)
		{
			goto FREE_RETURN;
		}
		std += (param->shift_type == 1? std * param->vol_shift: param->vol_shift);
		if (power > 0.5)
		{
			power = srt_f_optblksch(
				ex_lfwd[i], 
				ex_lfwd[i], 
				std, 
				ex_time[i], 
				1.0, 
				SRT_CALL, 
				PREMIUM);
			err = srt_f_optimpvol(
				power,
				ex_lfwd[i],
				ex_lfwd[i],
				ex_time[i],
				1.0,
				SRT_CALL,
				SRT_NORMAL,
				&std);
		}
		std *= sqrt (ex_time[i]);
		
		/*	Strike */
		if (param->strike_type == 0 || strike[i] < 1.0e-04)
		{
			ex_lstrike[i] = ex_lfwd[i];
		}
		else
		if (param->strike_type == 1)
		{
			ex_lstrike[i] = strike[i];
		}
		else
		if (param->strike_type == 2)
		{
			if (err = swp_f_ForwardRate(
				cpn_date[j],
				theo_end_dates[i],
				instr_freq,
				instr_basis,
				yc_name,
				ref_rate_name,
				&swp_rte))
			{
				goto FREE_RETURN;
			}

			spr = swp_rte - ex_lfwd[i];
			
			ex_lstrike[i] = strike[i] - spr;
		}
		else
		if (param->strike_type == 3)
		{
			ex_lstrike[i] = ex_lfwd[i] + strike[i] * std;
		}

		/*	Apply max std */
		if (ex_lstrike[i] > ex_lfwd[i] + param->max_std * std)
		{
			ex_lstrike[i] = ex_lfwd[i] + param->max_std * std;
		}
		else if (ex_lstrike[i] < ex_lfwd[i] - param->max_std * std)
		{
			ex_lstrike[i] = ex_lfwd[i] - param->max_std * std;
		}

		/*	Make sure strikes are positive (actually more than 1bp)
				otherwise use ATM	*/
		if (ex_lstrike[i] < 1.0e-04)
		{
			ex_lstrike[i] = ex_lfwd[i];
		}

		err = get_cash_vol (
			vol_curve_name, 
			add_unit (ex_date[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING),
			theo_end_dates[i], 
			ex_lstrike[i], 
			0, 
			ref_rate_name, 
			&(ex_lvol[i]), 
			&power);
		if (err)
		{
			goto FREE_RETURN;
		}
		ex_lvol[i] += (param->shift_type == 1? ex_lvol[i] * param->vol_shift: param->vol_shift);

		if (power > 0.5)
		{
			ex_lprice[i] = srt_f_optblksch(
				ex_lfwd[i],
				ex_lstrike[i],
				ex_lvol[i],
				ex_time[i],
				ex_llvl[i],
				SRT_CALL,
				PREMIUM);
		}
		else
		{
			ex_lprice[i] = srt_f_optblknrm(
				ex_lfwd[i],
				ex_lstrike[i],
				ex_lvol[i],
				ex_time[i],
				ex_llvl[i],
				SRT_CALL,
				PREMIUM);
		}
	}
	
	CouponTime = (double **) calloc(nex, sizeof(double *));
	Coupon = (double **) calloc(nex, sizeof(double *));	
	iNbCoupon = (int *) calloc(nex, sizeof(int));

	if (!CouponTime || !Coupon || !iNbCoupon)
	{
		err = "Memory allocation failure";
		goto FREE_RETURN;
	}

	TStar_date = today + dTStar * DAYS_IN_YEAR;

	/* DF(Today,Tstar) 
	dDFTodayTstar = swp_f_df(today, TStar_date, yc_name);
	*/

	for (i=0; i<nex; i++)
	{
		CouponTime[i] = &(cpn_time[ex_cpn[i]]);
		
		iNbCoupon[i] = ex_endcpn[i] - ex_cpn[i] + 1;
		Coupon[i] = (double *) calloc (iNbCoupon[i], sizeof (double));
		Coupon[i][0] = 1.0;
		
		for (j=1; j<iNbCoupon[i]; j++)
		{
			Coupon[i][j] = -ex_lstrike[i] * cpn_cvg[j+ex_cpn[i]];
		}

		Coupon[i][iNbCoupon[i]-1] -= 1.0;
				
		/*
		dDFFTexTstar = dDFTodayTstar / swp_f_df(lExDate, today + DAYS_IN_YEAR * ex_time[i], yc_name);
		*/

		for (j=0; j<iNbCoupon[i]; j++)
		{
			Coupon[i][j] *= cpn_df[j+ex_cpn[i]];
		}
	}

	
	*num_sig = nex;
	*sig_time = (double*) calloc (nex, sizeof (double));
	*sig = (double*) calloc (nex, sizeof (double));
	
	if (!sig_time || !sig)
	{
		err = "Allocation error (3) in cpd_calib_diagonal";
		goto FREE_RETURN;
	}

	for (i=0; i<nex; i++)
	{
		(*sig_time)[i] = ex_time[i];
	}

	/* 3 ) Computation of the first guess */

	/* First find the equivalent ATM volatilities */	

	if (power > 0.5)
	{
		for (i=0; i<nex; i++)
		{
			AlphaEq = dAlpha;
			if (dLambdaEps > 1.0E-16)
			{
				AlphaEq = dAlpha * sqrt((1.0 - exp(-2.0 * dLambdaEps * ex_time[i])) / (2.0 * dLambdaEps * ex_time[i]));
			}

			vol_conv(ex_lvol[i],
					 SABR_STR_LOG,
					 &atm_vol,
					 SABR_ATM_LOG,
					 ex_lfwd[i],
					 ex_lstrike[i],
					 ex_time[i],
					 AlphaEq,
					 0.0,
					 dRho);

			atm_price[i] = srt_f_optblksch(	ex_lfwd[i],
											ex_lfwd[i],
											atm_vol,
											ex_time[i],
											ex_llvl[i],
											SRT_CALL,
											PREMIUM);
		}
	}
	else
	{
		for (i=0; i<nex; i++)
		{
			AlphaEq = dAlpha;
			if (dLambdaEps > 1.0E-16)
			{
				AlphaEq = dAlpha * sqrt((1.0 - exp(-2.0 * dLambdaEps * ex_time[i])) / (2.0 * dLambdaEps * ex_time[i]));
			}

			vol_conv(ex_lvol[i],
					 SABR_STR_NORM,
					 &atm_vol,
					 SABR_ATM_NORM,
					 ex_lfwd[i],
					 ex_lstrike[i],
					 ex_time[i],
					 AlphaEq,
					 0.001,
					 dRho);

			atm_price[i] = srt_f_optblknrm(	ex_lfwd[i],
											ex_lfwd[i],
											atm_vol,
											ex_time[i],
											ex_llvl[i],
											SRT_CALL,
											PREMIUM);
		}
	}

	err = lgmprcapgivenlambda_2(ncpn,
								cpn_time,
								cpn_df,
								cpn_cvg,
								nex,
								ex_time,
								ex_cpn,
								ex_endcpn,
								ex_lfwd,
								atm_price,
								ex_zeta,
								dLambdaX,
								1,
								0.0,
								0.0,
								0.0,
								0);	
	if (err)
	{
		goto FREE_RETURN;
	}

	(*sig)[0] = sqrt (ex_zeta[0] * 2 * (dLambdaX) / (exp (2 * (dLambdaX) * ex_time[0]) - 1.0));

	for (i=1; i<nex; i++)
	{
		(*sig_time)[i] = ex_time[i];
		if (ex_zeta[i] > ex_zeta[i-1])
		{
			(*sig)[i] = sqrt ((ex_zeta[i] - ex_zeta[i-1]) * 2 * (dLambdaX)
								/ (exp (2 * (dLambdaX) * ex_time[i]) - exp ( 2 * (dLambdaX) * ex_time[i-1])));
		}
		else
		{
			smessage ("Diagonal calibration failed at exercise year %.2f - Calibration stopped", ex_time[i]);
			for (j=i; j<nex; j++)
			{
				(*sig)[j] = (*sig)[i-1];
			}
			i = nex;
		}
	}

	/* 4 ) calibration */

	dAlpha *= 2.0;
	dLambdaEps *= 2.0;

	LGMSVCalib2(							
				dLambdaX,
				dTStar,
				dAlpha,
				dLambdaEps,
				dRho,
				yc_name,
				ex_time,
				nex,
				iNbCoupon,
				CouponTime,
				Coupon,
				ex_lprice,
				ex_lvol,
				ex_lfwd,
				ex_lstrike,
				ex_llvl,
				iNbPhi,
				iNbft,
				iNbSigmaPhiGridLeft,
				iNbSigmaPhiGridRight,
				iNbSigmaftLeft,
				iNbSigmaftRight,
				iRatioPhi,
				iRatioFt,
				iPriorityFreqPhi,
				iPriorityFreqFt,
				iMaxTime,						
				Precision,
				NbIterMax,
				*sig);

	/*	5.)	Save instrument data if required */
	if (inst_data)
	{
		inst_data->num_inst = nex;
		inst_data->start_dates = (long*) calloc (nex, sizeof (long));
		inst_data->end_dates = (long*) calloc (nex, sizeof (long));
		inst_data->short_strikes = (double*) calloc (nex, sizeof (double));
		inst_data->long_strikes = (double*) calloc (nex, sizeof (double));

		if (!inst_data->start_dates 
			|| !inst_data->end_dates
			|| !inst_data->short_strikes
			|| !inst_data->long_strikes)
		{
			err = "Allocation error (4) in cpd_calib_diagonal";
			goto FREE_RETURN;
		}

		for (i=0; i<nex; i++)
		{
			inst_data->start_dates[i] = cpn_date[ex_cpn[i]];
			inst_data->end_dates[i] = cpn_date[ex_endcpn[i]];
			inst_data->short_strikes[i] = 0.0;
			inst_data->long_strikes[i] = ex_lstrike[i];
		}
	}

FREE_RETURN:

	if (err)
	{
		if (*sig_time) free (*sig_time);
		*sig_time = NULL;

		if (*sig) free (*sig);
		*sig = NULL;

		if (inst_data)
		{
			cpd_free_calib_inst_data (inst_data);
		}
	}

	if (CouponTime) free (CouponTime);
	if (Coupon)
	{
		for (i=0; i<nex; i++)
		{
			if (Coupon[i]) free (Coupon[i]);
		}
		free (Coupon);
	}
	if (iNbCoupon) free (iNbCoupon);

	return err;
}