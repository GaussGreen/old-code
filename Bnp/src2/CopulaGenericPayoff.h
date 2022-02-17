#ifndef __COPULAGENERICPAYOFF_H
#define __COPULAGENERICPAYOFF_H


#include "srt_h_types.h"


/* Generic Copula structure */
/* ------------------------ */
typedef enum {GAUSSIAN, SV} CopulaType;

typedef struct
{
	CopulaType			iCopulaType;	/* Type of copula */
	
	long				iNbUnderlying;	/* Number of underlying */
	long				lNumPaths;		/* Number of paths */	
	long				lNbPoints;		/* Nb of points for approx of marginal cumulative */
	SrtMCSamType		iMCType;		/* Type of the used MC method */

	/* Correlation */
	int					iIsCorrMatrix;
	double				**dCorrMatrix;

	/* Other Params*/
	void				*OtherParams;

} Generic_Copula, *GENERIC_COPULA;

typedef struct
{
	
	double				dMaturity;
	double				*dFwds;
	double				*dVols;
	double				*dAlpha;
	double				*dBeta;
	double				*dRho;
	double				**dCorrMatrix;
	SrtDiffusionType	eTypeInput;

} Specific_SV_Copula, *SPECIFIC_SV_COPULA;

/* Specific Copula structure */
/* ------------------------- */


/* Generic Model structure */
/* ----------------------- */
typedef enum {BMM, SL} ModelType;

typedef struct
{
	ModelType	iModelType;

	/* Specific Model Parameter */
	void		*ModelParameter;
}
Generic_Model, *GENERIC_MODEL;


/* Specific Model structure */
/* ------------------------ */
typedef struct
{
	double	dForward1;
	double	dSigmaBeta1;
	double	dForward2;
	double	dSigmaBeta2;
	double	dBeta;
	double	dPi;

	/* For SABR Calibration */
	double	dNStdCalib;

	/* For Marginal calculation */
	double	dNStd;
	double	dMaxError;

} BMM_Model, *BMM_MODEL;


/* Payoff structures */
/* ----------------- */

typedef enum PayOffCallPutDigType_
{
	PAYOFF_CALL,
	PAYOFF_PUT,
	PAYOFF_IV,
	PAYOFF_STRADDLE,
	PAYOFF_DIGUP,
	PAYOFF_DIGDOWN,
	PAYOFF_DIGUP_CS,	/* Digital Call with Call Spread */
	PAYOFF_DIGDOWN_CS	/* Digital Put	with Call Spread */
} PayOffCallPutDigType;

/* PayOff_ProductofCalls 
	p underlying of value X[j] (j=0 to p-1)
	n number of product
	Payoff = product (i=0 to n-1) of Max(0 ; Sum (j=1->p) a[i,j] X[j]) - Strike[i]) if  iCallPut[i] = Call 
									 Max(0 ; Strike[i] - Sum (j=1->p) a[i,j] X[j]) ) if iCallPut[i] = Put */
typedef struct
{
	int		iNbProduct;			/* Nb of product of max */
	int		iNbUnderlying;		/* Nb of undrlying */
	
	double	*dStrike;			/* Strike  */
	PayOffCallPutDigType *iCallPut;	/* CallPut */
	double	**a;				/* coef for all underlying */
	double	*b;					/* coef on X/Y */
	double	call_spread;		/* Call Spread width used for  	PAYOFF_DIGUP_CS and PAYOFF_DIGDOWN_CS*/



} PRODUCT_Call_Payoff, *PRODUCT_CALL_PAYOFF;

/* PayOff_X_div_Y 

  max(max(X - Kx, 0) / Y - Kxdivy,0)
  1/Y can be replicate by put as sum Ni*(Ki-Y)+

 */
typedef struct
{
	int		iNbUnderlying;		/* Nb of undrlying */
	
	double	*dStrike;			/* Strike  */
	PayOffCallPutDigType *iCallPut;	/* CallPut */

	int		isYfloored;
	double	floorY;
	int		isYcapped;
	double	capY;

	int		iUseReplication;		/* Use Replication */
	int		iNbReplicatedStrike;	/* Nb of replicated strike */
	double	*dReplicatedStrike;		/* replicated strike */
	double	*dNotional;				/* Notional */

} XDIVY_Payoff;


typedef struct
{
	double	*fwd;
	double	*volsqrT;

	int		iNbUnderlying;		/* Nb of underlying */
	double	**dSqCorrMatrix;
	double	CoefDensity;
	
	void	*PayOffParam;
	double	*UnderlyingValue;
	Err (*PayOff) (/* Underlying value */
												   int		iNbUnderlying,
												   double	*UnderlyingValue,

												   /* For Payoff description */
												   int		iNbProduct,
												   void		*PayOffParam,
												   
												   /* OutPuts */
												   double	*dPV);
} IntGaussParam;


/* ----------------------------------------------------------------------------------------------------------*/
/*											Functions														 */
/* ----------------------------------------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------------------------------------
	initGenericModel
		Allocate inside 
   ----------------------------------------------------------------------------------------------------------- */
Err initGenericModel(/* Input */
					 int			iNbUnderlying,
					 ModelType		iModelType,

					 /* For BMM type */
					 double			dNStdBMMCalib,

					 /* OutPut */
					 Generic_Model	**GenericModel,
					 void			**SpecificModel);

Err interp_call_put_dig(const char *constStr, PayOffCallPutDigType *val);

/* ---------------------------------------------------------------------------
	
	PayOff_ProductofCalls
	
	p underlying of value X[j] (j=0 to p-1)
	n number of product
	Payoff = product (i=0 to n-1) of Max(0 ; Sum (j=1->p) a[i,j] X[j]) - Strike[i]) if  iCallPut[i] = Call 
									 Max(0 ; Strike[i] - Sum (j=1->p) a[i,j] X[j]) ) if iCallPut[i] = Put  

  ---------------------------------------------------------------------------- */
Err PayOff_ProductofCalls(/* Underlying value */
						  int		iNbUnderlying,
						  double	*UnderlyingValue,

						  /* For Payoff description */
						  int		iNbProduct,
						  void		*PayOffParam,
						  
						  /* OutPuts */
						  double	*dPV);

Err PayOff_ProductofCalls1(/* Underlying value */
							int		iNbUnderlying,
							double	*UnderlyingValue,

							/* For Payoff description */
							int		iNbProduct,
							void	*PayOffParam,
							
							/* OutPuts */
							double	*dPV);

Err PayOff_XDivY(/* Underlying value */
							int		iNbUnderlying,
							double	*UnderlyingValue,

							/* For Payoff description */
							int		iNbProduct,
							void	*PayOffParam,
							
							/* OutPuts */
							double	*dPV);

Err PayOff_MinofCalls(/* Underlying value */
							int		iNbUnderlying,
							double	*UnderlyingValue,

							/* For Payoff description */
							int		iNbProduct,
							void	*PayOffParam,
							
							/* OutPuts */
						  double	*dPV);
/* ---------------------------------------------------------------------------
	
	GaussianCopulaGetSamples
	
	From the marginales and a gaussian correlation get all the samples
	dSamples[iNumPath][iNumUnd]

  ---------------------------------------------------------------------------- */
Err	GaussianCopulaGetSamples(/* Copula Parameter */
							  Generic_Copula	*CopulaParam,

							  /* Marginales Distributions */
							  long		*lNbPoints,		/* Array of Nb points in the cumulative */
							  double	**dX,			/* dX[iNumUnderlying][iNumPoints]  */
							  double	**dCumulative,	/* dCumulative[iNumUnderlying][iNumPoints]  */

							  /* Ouput */
							  double	**dSamples		/* dSamples[iNumPath][iNumUnd]  */);

Err	SVCopulaGetSamples(/* Copula Parameter */
						  Generic_Copula	*CopulaParam,

						  /* Marginales Distributions */
						  long		*lNbPoints,		/* Array of Nb points in the cumulative */
						  double	**dX,			/* dX[iNumUnderlying][iNumPoints]  */
						  double	**dCumulative,	/* dCumulative[iNumUnderlying][iNumPoints]  */

						  /* Ouput */
						  double	**dSamples		/* dSamples[iNumPath][iNumUnd]  */);

/* ---------------------------------------------------------------------------
	
	Copula_GenericPayOff_Price

  ---------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------
	
	Copula_GenericPayOff_Price
	dPremium should be allocated

  ---------------------------------------------------------------------------- */
Err Copula_GenericPayOff_Price (	/* Copula Parameter */
									Generic_Copula			*CopulaParam,
									
									/* PayOff parameters */
									int						iNbProduct,
									void					*PayOffParam,
						  								  
									/* For generating Marginales Distributions */
									double					dMaturity,
									double					*dForward,
									Generic_Model			*UnderlyingModelParam,

									/* Copula function */
									Err	(*CopulaGetSamples) (/* Copula Parameter */
															Generic_Copula	*CopulaParam,

															/* Marginales Distributions */
															long		*lNbPoints,		/* Array of Nb points in the cumulative */
															double	**dX,			/* dX[iNumUnderlying][iNumPoints]  */
															double	**dCumulative,	/* dCumulative[iNumUnderlying][iNumPoints]  */

															/* Ouput */
															double	**dSamples		/* dSamples[iNumPath][iNumUnd]  */),

														  
									/* Payoff Function */
									Err (*PayOff) (/* Underlying value */
												   int		iNbUnderlying,
												   double	*UnderlyingValue,

												   /* For Payoff description */
												   int		iNbProduct,
												   void		*PayOffParam,
												   
												   /* OutPuts */
												   double	*dPV),

									/* Results */
									double					*dPremium);

/* ---------------------------------------------------------------------------
	
	Copula_SABRCalib_GenericPayOff_Price

  ---------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
	
	Copula_SABRCalib_GenericPayOff_Price

  ---------------------------------------------------------------------------- */

Err Copula_SABRCalib_GenericPayOff_Price(
										/* Copula Parameter */
										Generic_Copula			*CopulaParam,
										 
										/* PayOff parameters */
										int						iNbProduct,
										void					*PayOffParam,

										/* For generating Marginales Distributions */
										double					dMaturity,
										double					*dForward,
										Generic_Model			*UnderlyingModelParam,

										/* Sabr Parameters */
										double					*dSigma,
										double					*dAlpha,
										double					*dBeta,
										double					*dRho,
										
										SrtDiffusionType		eTypeInput,

										/* Copula function */
										Err	(*CopulaGetSamples) (/* Copula Parameter */
																Generic_Copula	*CopulaParam,

																/* Marginales Distributions */
																long	*lNbPoints,		/* Array of Nb points in the cumulative */
																double	**dX,			/* dX[iNumUnderlying][iNumPoints]  */
																double	**dCumulative,	/* dCumulative[iNumUnderlying][iNumPoints]  */

																/* Ouput */
																double	**dSamples		/* dSamples[iNumPath][iNumUnd]  */),

															  
										/* Payoff Function */
										Err (*PayOff) (/* Underlying value */
												   int		iNbUnderlying,
												   double	*UnderlyingValue,

												   /* For Payoff description */
												   int		iNbProduct,
												   void		*PayOffParam,
												   
												   /* OutPuts */
												   double	*dPV),

										/* Results */
										double					*dPremium);


/* Multi trapezoidal Integration */
double	TrapIntegration(/* limite */
						double	a,
						double	b,
						

						int		iNbPoints,
						int		iNumIntegration,

						/* Function parameters */
						int		iNbX,
						double	*X,
						IntGaussParam	*Param,

						double (*f)(int	iNbX,
									double *X,
									IntGaussParam *Param));


double	PayOffGaussianDensity(int	iNbX,
							  double *X,
						      IntGaussParam *Param);


#endif
