


#include        "utallhdr.h"
#include        <NUM_H_ALLHDR.H>

/* ========================================================================== */

#ifndef		BGMInterpolVol_H
#define		BGMInterpolVol_H


Err srt_BGMInterpVols(double *OldMats, 
					  int nummats,
					  double *GivenVols,
					  SrtDiffusionType InputVol,
					  SrtDiffusionType OutputVol,
					  char *szRefRateCode,
					  char *szYieldCurveName,
					  double **NewVols,
					  long *MaxNumPeriod);


Err srt_f_BGMGetATMSwaptionVolFromFraVols(double *TenorMats,
										  int nTenorMats,
										  double *FraVols,
										  char *szRefRateCode,
										  char *szYieldCurveName,
										  char *VolAssumption,
										  SrtDiffusionType	StatioVolLnNorm,
										  char *CorrAssumption,
										  double **correl,
										  long	ncorr,
										  long  swaptionmat,
										  long  swaptionund,
										  char *swapFreq,
										  char *swapBasis,
										  double **swaptionvol);

Err srt_f_BGMGetATMSwaptionVolFromFraVols1Parameter(double *TenorMats,
										  int nTenorMats,
										  double *FraVols,
										  char *szRefRateCode,
										  char *szYieldCurveName,
										  char *VolAssumption,
										  SrtDiffusionType InputVol,
										  SrtDiffusionType OutputVol,
										  char *CorrAssumption,
										  double correl,
										  double VolMat,
										  double  swaptionmat,
										  double  swaptionund,
										  char *swapFreq,
										  char *swapBasis,
										  double **swaptionvol);



Err srt_f_StripCapletVols(double *CapMaturities,
						  int nCapMaturities,
						  double *CapStrikes,
						  double *CapVols,
						  double *CapletVols, /*this one already contains the first caplets provided by the user, we here are completing it*/
						  int ncapletvolsknown,
						  char *szRefRateCode,
						  char *szYieldCurveName,
						  double *Alpha,
						  double *Beta,
						  double *Rho);


Err srt_f_GenCapletVol(  double Maturity,
						 double HumpMaturity,
						 double HumpNormalVol,
						 double FirstMat,
						 double FirstSpread,
						 double SecondMat,
						 double SecondSpread,
						 double BackEndParam,
						 double CutoffMat,
						 double VarDecRatio,
						 double *normalvol);

#endif
