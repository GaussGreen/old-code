
#ifndef __VOLBOND_PROD_STRUCT_H
#define __VOLBOND_PROD_STRUCT_H


//---------------------------------------------------------------------------
//----------------------------	mylgm2f	-------------------------------------
//---------------------------------------------------------------------------

typedef struct
{
	/* Df parameters at evt_date */
	long					NDates;
	double                  *lam1;
	double                  *lam2;
//	double                  *sigmaSquare;
//	double                  *initDf;
	double                  *Df;
//	double                  rho;
	double                  *AdjDf;
	
	/* distribution of X1 and X2 */
	double                  esp1;
	double                  esp2;
	double                  var1;
	double                  var2;
	double                  covar;

} my_lgm2F, *MYLGM2F;


Err Fill_MY_LGM2F(long NDates, MYLGM2F mylgm2f);

Err Free_MY_LGM2F(MYLGM2F mylgm2f);

Err MY_LGM2F_Preliminary_Computation_For_Df(long today, long eventDate, 
											long *Dates, MYLGM2F mylgm2f, 
											TermStruct* l, char *yc);

void MY_LGM2F_Compute_Df_From_Factors_Value(double fact1, 
											double fact2, 
											MYLGM2F mylgm2f);

//void MY_LGM2F_Compute_Df_From_Factors_Value(long today, long eventDate, 
//										   double fact1, double fact2, 
//										   MYLGM2F mylgm2f, TermStruct* l, 
//										   char *yc);

Err MY_LGM2F_Compute_Mean_Variance_Covariance(long today, long eventDate, 
											  long probaDate, MYLGM2F mylgm2f, 
											  TermStruct* l, char *yc);

//---------------------------------------------------------------------------
//----------------------	swap_for_lgm2f	---------------------------------
//---------------------------------------------------------------------------
typedef struct
{

	long SwapStartDate;
	long SwapEndDate;

//		float leg
	long float_NDates;
	long *float_pay_dates;
	long *float_start_dates;
	long *float_end_dates;

	double *float_cvgs;
	double *float_spreads;

	my_lgm2F float_pay_lgm2f;
	my_lgm2F float_start_lgm2f;
	my_lgm2F float_end_lgm2f;

//		fix leg
	long fix_NDates;
	long *fix_pay_dates;

	double *fix_cvgs;

	my_lgm2F fix_pay_lgm2f;

//		Swap Rate
	double swapRate;

} swap_for_lgm2f, *SWAPFORLGM2F;


Err Fill_SWAP_FOR_LGM2F(long today, long SwapStartDate, long SwapEndDate, 
						char *SwapFreq, char *SwapBasis, char *RefRate,
						SWAPFORLGM2F swapForLgm2f);

Err SWAP_FOR_LGM2F_Preliminary_Computation(long today, long eventDate, 
											swap_for_lgm2f *swapForLgm2f,
											TermStruct* l, char *yc);

double SWAP_FOR_LGM2F_SwapRate_From_Factors_Value(double fact1, double fact2, swap_for_lgm2f *swapForLgm2f);


//double SWAP_FOR_LGM2F_SwapRate_From_Factors_Value(long today, long eventDate, 
//											double fact1, double fact2, 
//											swap_for_lgm2f *swapForLgm2f,
//											TermStruct* l, char *yc);

Err Free_SWAP_FOR_LGM2F(SWAPFORLGM2F swapForLgm2f);



//---------------------------------------------------------------------------
//---------------------		VolBond Struct		-----------------------------
//---------------------------------------------------------------------------

typedef struct
{
	double notional;

// First Swap rate parameters
	char *first_SwapFreq;
	char *first_SwapBasis;
	char *first_SwapRefRate;
	long first_fixingDate;	// fixing date of the first swap rate
	long first_startDate;	// start date of the first swap rate
	long first_endDate;	// end date of the first swap rate
	swap_for_lgm2f firstSwap;

// Second Swap rate parameters
	char *second_SwapFreq;
	char *second_SwapBasis;
	char *second_SwapRefRate;
	long second_fixingDate;	// fixing date of the second swap rate
	long second_startDate;	// start date of the second swap rate
	long second_endDate;	// end date of the second swap rate
	swap_for_lgm2f secondSwap;

// Coupon's payment date
	long payment_Date;

// Payoff Parameters
	double alphaC;
	double alphaP;
	double betaC;
	double betaP;
	double strikeC;
	double strikeP;

} volBondCoupon, *VOLBONDCOUPON;


typedef struct
{
	long NbOfCoupons;
	double notional;
	volBondCoupon *VolBondCpn;

} volBondStruct, *VOLBONDSTRUCT;

Err Fill_VolBondCoupon(long today, 
					   char *first_SwapFreq, char *first_SwapBasis, 
					   char *first_SwapRefRate,
					   char *second_SwapFreq, char *second_SwapBasis, 
					   char *second_SwapRefRate,
					   double notional,
					   long first_fixingDate, long first_startDate, 
					   long first_endDate,
					   long second_fixingDate, long second_startDate, 
					   long second_endDate,
					   long payment_Date, double alphaC, double alphaP,
					   double betaC, double betaP, double strikeC, double strikeP,
					   volBondCoupon *volBondCpn);

Err PayOff_VolBondCoupon(long today, double fact1, double fact2, 
							volBondCoupon *volBondCpn, 
							double CMSAdj, 
							double sabrAlpha,
							double sabrBeta,
							double sabrBetaVol,
							double sabrRho,
							double *Payoff,
							TermStruct* l, char *yc);

Err Fill_VolBondStruct(long today, long NbCoupons, double *notional,
					   char **first_SwapFreq, char **first_SwapBasis, 
					   char **first_SwapRefRate,
					   char **second_SwapFreq, char **second_SwapBasis, 
					   char **second_SwapRefRate,
					   long *first_fixingDate, long *first_startDate, 
					   long *first_endDate,
					   long *second_fixingDate, long *second_startDate, 
					   long *second_endDate,
					   long *payment_Date, 
					   double *alphaC, double *alphaP, double *betaC, double *betaP, 
					   double *strikeC, double *strikeP,
					   volBondStruct *volBond);

Err Free_VolBondCoupon(volBondCoupon *volBondCpn);


Err Free_VolBondStruct(volBondStruct *volBond);


#endif
