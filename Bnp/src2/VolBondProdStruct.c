#include "srt_h_all.h"
#include "opfnctns.h"
#include "math.h"
#include "VolBondProdStruct.h"


//-----------------------------------------------------------------------------
//----------------------	Functions for MY_LGM2F	---------------------------
//-----------------------------------------------------------------------------
Err Init_MY_LGM2F(MYLGM2F mylgm2f)
{
	Err err = NULL;

	mylgm2f->lam1 = NULL;
	mylgm2f->lam2 = NULL;
	mylgm2f->Df = NULL;
	mylgm2f->AdjDf = NULL;

	return err;
}


Err Fill_MY_LGM2F(long NDates, MYLGM2F mylgm2f)
{
	Err err = NULL;

	mylgm2f->NDates = NDates;

	mylgm2f->lam1 = NULL;
	mylgm2f->lam1 = (double*) calloc(NDates, sizeof(double));
	if (!mylgm2f->lam1)
	{
		err = "Allocation error in Fill_MY_LGM2F";
		goto FREE_RETURN;
	}
	mylgm2f->lam2 = NULL;
	mylgm2f->lam2 = (double*) calloc(NDates, sizeof(double));
	if (!mylgm2f->lam2)
	{
		err = "Allocation error in Fill_MY_LGM2F";
		goto FREE_RETURN;
	}
//	mylgm2f->sigmaSquare = (double*) calloc(NDates, sizeof(double));
//	if (!mylgm2f->sigmaSquare)
//	{
//		err = "Allocation error in Fill_MY_LGM2F";
//		goto FREE_RETURN;
//	}
//	mylgm2f->initDf = (double*) calloc(NDates, sizeof(double));
//	if (!mylgm2f->initDf)
//	{
//		err = "Allocation error in Fill_MY_LGM2F";
//		goto FREE_RETURN;
//	}
	mylgm2f->Df = NULL;
	mylgm2f->Df = (double*) calloc(NDates, sizeof(double));
	if (!mylgm2f->Df)
	{
		err = "Allocation error in Fill_MY_LGM2F";
		goto FREE_RETURN;
	}
	mylgm2f->AdjDf = NULL;
	mylgm2f->AdjDf = (double*) calloc(NDates, sizeof(double));
	if (!mylgm2f->AdjDf)
	{
		err = "Allocation error in Fill_MY_LGM2F";
		goto FREE_RETURN;
	}

FREE_RETURN:
	if (err)
	{
		err = Free_MY_LGM2F(mylgm2f);
	}

	return err;
}


Err Free_MY_LGM2F(MYLGM2F mylgm2f)
{
	if (mylgm2f->lam1)
	{
		free (mylgm2f->lam1);
		mylgm2f->lam1 = NULL;
	}
	if (mylgm2f->lam2)
	{
		free (mylgm2f->lam2);
		mylgm2f->lam2 = NULL;
	}
//	if (mylgm2f->sigmaSquare)
//	{
//		free (mylgm2f->sigmaSquare);
//		mylgm2f->sigmaSquare = NULL;
//	}
//	if (mylgm2f->initDf)
//	{
//		free (mylgm2f->initDf);
//		mylgm2f->initDf = NULL;
//	}
	if (mylgm2f->Df)
	{
		free (mylgm2f->Df);
		mylgm2f->Df = NULL;
	}
	if (mylgm2f->AdjDf)
	{
		free (mylgm2f->AdjDf);
		mylgm2f->AdjDf = NULL;
	}

	return NULL;
}


Err MY_LGM2F_Preliminary_Computation_For_Df(long today, long eventDate, 
											long *Dates, MYLGM2F mylgm2f, 
											TermStruct* l, char *yc)
{
	Err err = NULL;
	SrtTFTSVec FVec, PsiVec1, PsiVec2;
	SrtTFTSMat GVec;
	double Phi1, Phi2, Phi12;
	long i;
	double sigmaSquare;
	double timeEventDate, timeMaturity;

	timeEventDate = (eventDate - today) / 365.00;

	err = get_2f_F_funcs(timeEventDate, l, &FVec);
	if(err)
	{
		return err;
	}

	err = get_2f_G_funcs(timeEventDate, l, &GVec);
	if(err)
	{
		return err;
	}

	err = get_2f_Psi_funcs(timeEventDate, l, &PsiVec1);
	if(err)
	{
		return err;
	}

	Phi1 = FVec[0] * FVec[0] * GVec[0][0];
	Phi2 = FVec[1] * FVec[1] * GVec[1][1];
	Phi12 = FVec[0] * FVec[1] * GVec[0][1];

	for(i=0;i<mylgm2f->NDates;++i)
	{
		timeMaturity = (Dates[i] - today) / 365.00;
		err = get_2f_Psi_funcs(timeMaturity, l, &PsiVec2);
		if(err)
		{
			return err;
		}

		mylgm2f->lam1[i] = (PsiVec2[0] - PsiVec1[0])/FVec[0];
		mylgm2f->lam2[i] = (PsiVec2[1] - PsiVec1[1])/FVec[1];

		sigmaSquare = mylgm2f->lam1[i] * mylgm2f->lam1[i] * Phi1
					+ mylgm2f->lam2[i] * mylgm2f->lam2[i] * Phi2
					+ 2 * mylgm2f->lam1[i] * mylgm2f->lam2[i] * Phi12;

		mylgm2f->AdjDf[i] = (swp_f_df(today, Dates[i], yc)/swp_f_df(today, eventDate, yc))
						* exp(-0.5 * sigmaSquare);

	}

	return err;

}


Err MY_LGM2F_Compute_Mean_Variance_Covariance(long today, long eventDate, 
											  long probaDate, MYLGM2F mylgm2f, 
											  TermStruct* l, char *yc)
{
	SrtTFTSVec FVec, PsiVec1, PsiVec2;
	SrtTFTSMat GVec;
	double lam1, lam2;
	double Phi1, Phi2, Phi12;
	double timeEventDate, timeProbaDate;
	Err err = NULL;

	timeEventDate = (eventDate - today) / 365.00;
	timeProbaDate = (probaDate - today) / 365.00;

	if(err)
	{
		return err;
	}
	err = get_2f_F_funcs(timeEventDate, l, &FVec);
	if(err)
	{
		return err;
	}
	err = get_2f_G_funcs(timeEventDate, l, &GVec);
	if(err)
	{
		return err;
	}
	err = get_2f_Psi_funcs(timeEventDate, l, &PsiVec1);
	if(err)
	{
		return err;
	}

	Phi1 = FVec[0] * FVec[0] * GVec[0][0];
	Phi2 = FVec[1]*FVec[1]* GVec[1][1];
	Phi12 = FVec[0]*FVec[1]*GVec[0][1];

	err = get_2f_Psi_funcs(timeProbaDate, l, &PsiVec2);
	if(err)
	{
		return err;
	}

	lam1 = (PsiVec2[0] - PsiVec1[0])/FVec[0];
	lam2 = (PsiVec2[1] - PsiVec1[1])/FVec[1];

	mylgm2f->esp1 = - lam1 * Phi1 - lam2 * Phi12;
	mylgm2f->esp2 = - lam2 * Phi2 - lam1 * Phi12;

	mylgm2f->var1 = Phi1;
	mylgm2f->var2 = Phi2;

	mylgm2f->covar = Phi12;

	return err;

}


void MY_LGM2F_Compute_Df_From_Factors_Value(double fact1, double fact2, 
											MYLGM2F mylgm2f)//(long today, long eventDate, double fact1, double fact2, MYLGM2F mylgm2f, TermStruct* l, char *yc)
{
	long i;

	for(i=0;i<mylgm2f->NDates;++i)
	{
/*		mylgm2f->Df[i] = (mylgm2f->initDf[i] / swp_f_df(today, eventDate, yc))
						* exp(-0.5 * mylgm2f->sigmaSquare[i])
						* exp(-mylgm2f->lam1[i] * fact1)
						* exp(-mylgm2f->lam2[i] * fact2);
*/
		mylgm2f->Df[i] = mylgm2f->AdjDf[i] * exp(-mylgm2f->lam1[i] * fact1) 
										   * exp(-mylgm2f->lam2[i] * fact2);
	}

}


//-----------------------------------------------------------------------------
//------------------	Functions for SWAP_FOR_LGM2F	-----------------------
//-----------------------------------------------------------------------------

Err Fill_SWAP_FOR_LGM2F(long today, long SwapStartDate, long SwapEndDate, 
						char *SwapFreq, char *SwapBasis, char *RefRate,
						SWAPFORLGM2F swapForLgm2f)
{
	Err err = NULL;
	SwapDP swapdp;
	long float_nb_dates, float_nb_pay_dates;
	long *float_fixing_dates=NULL, *float_pay_dates=NULL;
	long fix_nb_dates, fix_nb_pay_dates;
	long *fix_start_dates=NULL, *fix_end_dates=NULL, *fix_pay_dates=NULL;
	long i;

	swapForLgm2f->SwapStartDate = SwapStartDate;
	swapForLgm2f->SwapEndDate = SwapEndDate;

	swapForLgm2f->float_pay_dates = NULL;
	swapForLgm2f->float_start_dates = NULL;
	swapForLgm2f->float_end_dates = NULL;

	swapForLgm2f->float_cvgs = NULL;
	swapForLgm2f->float_spreads = NULL;

	swapForLgm2f->fix_pay_dates = NULL;
	swapForLgm2f->fix_cvgs = NULL;

	Init_MY_LGM2F(&(swapForLgm2f->fix_pay_lgm2f));
	Init_MY_LGM2F(&(swapForLgm2f->float_start_lgm2f));
	Init_MY_LGM2F(&(swapForLgm2f->float_end_lgm2f));
	Init_MY_LGM2F(&(swapForLgm2f->float_pay_lgm2f));
	
	err = swp_f_initSwapDP(SwapStartDate, SwapEndDate, SwapFreq, SwapBasis, &swapdp);
	if(err)
	{
//		err = "Error in Fill_SWAP_FOR_LGM2F : Cannot initialize SwapDP";
		goto FREE_RETURN;
	}
	err = swp_f_make_FloatLegDatesCoveragesAndSpreads(&swapdp, today, RefRate, 
												&float_pay_dates,
												&float_nb_pay_dates, 
												&float_fixing_dates, 
												&(swapForLgm2f->float_start_dates),
												&(swapForLgm2f->float_end_dates),
												&(swapForLgm2f->float_cvgs), 
												&(swapForLgm2f->float_spreads), 
												&float_nb_dates);
	if(err)
	{
//		err = "Error in Fill_SWAP_FOR_LGM2F : Cannot make float leg";
		goto FREE_RETURN;
	}

	swapForLgm2f->float_pay_dates = (long *) calloc(float_nb_pay_dates-1, sizeof(long));
	if(!(swapForLgm2f->float_pay_dates))
	{
		err = "Error in Fill_SWAP_FOR_LGM2F : Allocation failed";
		goto FREE_RETURN;
	}

	swapForLgm2f->float_NDates = float_nb_dates;


	err = swp_f_make_FixedLegDatesAndCoverages(&swapdp, today, 
											&fix_pay_dates, 
											&fix_nb_pay_dates, 
											&fix_start_dates, 
											&fix_end_dates, 
											&(swapForLgm2f->fix_cvgs), 
											&fix_nb_dates);
	if(err)
	{
//		err = "Error in Fill_SWAP_FOR_LGM2F : Cannot make fixed leg";
		goto FREE_RETURN;
	}

	swapForLgm2f->fix_pay_dates = (long *) calloc(fix_nb_pay_dates-1, sizeof(long));
	if(!(swapForLgm2f->fix_pay_dates))
	{
		err = "Error in Fill_SWAP_FOR_LGM2F : Allocation failed";
		goto FREE_RETURN;
	}

	swapForLgm2f->fix_NDates = fix_nb_dates;

	for(i=1;i<max(fix_nb_pay_dates,float_nb_pay_dates);++i)
	{
		if(i<fix_nb_pay_dates)
		{
			swapForLgm2f->fix_pay_dates[i-1] = fix_pay_dates[i];
		}

		if(i<float_nb_pay_dates)
		{
			swapForLgm2f->float_pay_dates[i-1] = float_pay_dates[i];
		}
	}

	err = Fill_MY_LGM2F(fix_nb_dates, &(swapForLgm2f->fix_pay_lgm2f));
	if(err)
	{
//		err = "Error in Fill_SWAP_FOR_LGM2F : Cannot fill fix_pay MY_LGM2";
		goto FREE_RETURN;
	}

	err = Fill_MY_LGM2F(float_nb_dates, &(swapForLgm2f->float_start_lgm2f));
	if(err)
	{
//		err = "Error in Fill_SWAP_FOR_LGM2F : Cannot fill float start MY_LGM2";
		goto FREE_RETURN;
	}

	err = Fill_MY_LGM2F(float_nb_dates, &(swapForLgm2f->float_end_lgm2f));
	if(err)
	{
//		err = "Error in Fill_SWAP_FOR_LGM2F : Cannot fill float end MY_LGM2";
		goto FREE_RETURN;
	}

	err = Fill_MY_LGM2F(float_nb_dates, &(swapForLgm2f->float_pay_lgm2f));
	if(err)
	{
//		err = "Error in Fill_SWAP_FOR_LGM2F : Cannot fill float pay MY_LGM2";
		goto FREE_RETURN;
	}


FREE_RETURN:
	if(err)
	{
		Free_SWAP_FOR_LGM2F(swapForLgm2f);
	}

	if (float_fixing_dates)
	{
		free (float_fixing_dates);
		float_fixing_dates = NULL;
	}
	if (float_pay_dates)
	{
		free (float_pay_dates);
		float_pay_dates = NULL;
	}
	if (fix_start_dates)
	{
		free (fix_start_dates);
		fix_start_dates = NULL;
	}
	if (fix_end_dates)
	{
		free (fix_end_dates);
		fix_end_dates = NULL;
	}
	if (fix_pay_dates)
	{
		free (fix_pay_dates);
		fix_pay_dates = NULL;
	}

	return err;
}



Err Free_SWAP_FOR_LGM2F(SWAPFORLGM2F swapForLgm2f)
{
	Err err = NULL;
	
	if (swapForLgm2f->float_pay_dates)
	{
		free (swapForLgm2f->float_pay_dates);
		swapForLgm2f->float_pay_dates = NULL;
	}
	if (swapForLgm2f->float_start_dates)
	{
		free (swapForLgm2f->float_start_dates);
		swapForLgm2f->float_start_dates = NULL;
	}
	if (swapForLgm2f->float_end_dates)
	{
		free (swapForLgm2f->float_end_dates);
		swapForLgm2f->float_end_dates = NULL;
	}
	if (swapForLgm2f->float_cvgs)
	{
		free (swapForLgm2f->float_cvgs);
		swapForLgm2f->float_cvgs = NULL;
	}
	if (swapForLgm2f->float_spreads)
	{
		free (swapForLgm2f->float_spreads);
		swapForLgm2f->float_spreads = NULL;
	}

	Free_MY_LGM2F(&(swapForLgm2f->float_start_lgm2f));
	Free_MY_LGM2F(&(swapForLgm2f->float_end_lgm2f));
	Free_MY_LGM2F(&(swapForLgm2f->float_pay_lgm2f));

	if (swapForLgm2f->fix_pay_dates)
	{
		free (swapForLgm2f->fix_pay_dates);
		swapForLgm2f->fix_pay_dates = NULL;
	}
	if (swapForLgm2f->fix_cvgs)
	{
		free (swapForLgm2f->fix_cvgs);
		swapForLgm2f->fix_cvgs = NULL;
	}

	Free_MY_LGM2F(&(swapForLgm2f->fix_pay_lgm2f));
	
	return err;
}


Err SWAP_FOR_LGM2F_Preliminary_Computation(long today, long eventDate, 
											swap_for_lgm2f *swapForLgm2f,
											TermStruct* l, char *yc)
{
	Err err = NULL;

	err = MY_LGM2F_Preliminary_Computation_For_Df(today, eventDate, 
											swapForLgm2f->fix_pay_dates, 
											&(swapForLgm2f->fix_pay_lgm2f), 
											l, yc);
	if(err)
	{
		err = "Error in SWAP_FOR_LGM2F_Preliminary_Computation : fix pay Df";
		return err;
	}

	err = MY_LGM2F_Preliminary_Computation_For_Df(today, eventDate, 
											swapForLgm2f->float_start_dates, 
											&(swapForLgm2f->float_start_lgm2f), 
											l, yc);
	if(err)
	{
		err = "Error in SWAP_FOR_LGM2F_Preliminary_Computation : float start Df";
		return err;
	}

	err = MY_LGM2F_Preliminary_Computation_For_Df(today, eventDate, 
											swapForLgm2f->float_end_dates, 
											&(swapForLgm2f->float_end_lgm2f),
											l, yc);
	if(err)
	{
		err = "Error in SWAP_FOR_LGM2F_Preliminary_Computation : float end Df";
		return err;
	}

	err = MY_LGM2F_Preliminary_Computation_For_Df(today, eventDate, 
											swapForLgm2f->float_pay_dates, 
											&(swapForLgm2f->float_pay_lgm2f),
											l, yc);
	if(err)
	{
		err = "Error in SWAP_FOR_LGM2F_Preliminary_Computation : float pay Df";
		return err;
	}

	return err;
}


double SWAP_FOR_LGM2F_SwapRate_From_Factors_Value(double fact1, double fact2, swap_for_lgm2f *swapForLgm2f)//(long today, long eventDate, double fact1, double fact2, swap_for_lgm2f *swapForLgm2f, TermStruct* l, char *yc)
{
	long i;
	double FloatPV;
	double LVL;

	MY_LGM2F_Compute_Df_From_Factors_Value(fact1, fact2, 
											&(swapForLgm2f->fix_pay_lgm2f));//(today, eventDate, fact1, fact2, &(swapForLgm2f->fix_pay_lgm2f), l, yc);

	MY_LGM2F_Compute_Df_From_Factors_Value(fact1, fact2, 
											&(swapForLgm2f->float_start_lgm2f));//(today, eventDate, fact1, fact2, &(swapForLgm2f->float_start_lgm2f), l, yc);

	MY_LGM2F_Compute_Df_From_Factors_Value(fact1, fact2, 
											&(swapForLgm2f->float_end_lgm2f));//(today, eventDate, fact1, fact2, &(swapForLgm2f->float_end_lgm2f), l, yc);

	MY_LGM2F_Compute_Df_From_Factors_Value(fact1, fact2, 
											&(swapForLgm2f->float_pay_lgm2f));//(today, eventDate, fact1, fact2, &(swapForLgm2f->float_pay_lgm2f), l, yc);

	LVL=0;
	FloatPV=0;
	for(i=0;i<max(swapForLgm2f->fix_NDates, swapForLgm2f->float_NDates);++i)
	{
		if(i<swapForLgm2f->fix_NDates)
		{
			LVL += (swapForLgm2f->fix_pay_lgm2f).Df[i] * (swapForLgm2f->fix_cvgs)[i];
		}

		if(i<swapForLgm2f->float_NDates)
		{
			FloatPV += (swapForLgm2f->float_pay_lgm2f).Df[i]*(swapForLgm2f->float_cvgs)[i]*swapForLgm2f->float_spreads[i];
			FloatPV += (swapForLgm2f->float_pay_lgm2f).Df[i]*((swapForLgm2f->float_start_lgm2f).Df[i]/(swapForLgm2f->float_end_lgm2f).Df[i] - 1);
		}
	}

	return FloatPV / LVL;
}


//-----------------------------------------------------------------------------
//----------------	Functions for VolBondCoupon	and VolBondStruct	-----------
//-----------------------------------------------------------------------------

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
					   volBondCoupon *volBondCpn)
{
	Err err;

	volBondCpn->notional = notional;

	volBondCpn->first_SwapBasis = first_SwapBasis;
	volBondCpn->first_SwapFreq = first_SwapFreq;
	volBondCpn->first_SwapRefRate = first_SwapRefRate;

	volBondCpn->first_fixingDate = first_fixingDate;
	volBondCpn->first_startDate = first_startDate;
	volBondCpn->first_endDate = first_endDate;

	volBondCpn->second_SwapBasis = second_SwapBasis;
	volBondCpn->second_SwapFreq = second_SwapFreq;
	volBondCpn->second_SwapRefRate = second_SwapRefRate;

	volBondCpn->second_fixingDate = second_fixingDate;
	volBondCpn->second_startDate = second_startDate;
	volBondCpn->second_endDate = second_endDate;

	volBondCpn->payment_Date = payment_Date;

	volBondCpn->alphaC = alphaC;
	volBondCpn->alphaP = alphaP;
	volBondCpn->betaC = betaC;
	volBondCpn->betaP = betaP;
	volBondCpn->strikeC = strikeC;
	volBondCpn->strikeP = strikeP;

	err = Fill_SWAP_FOR_LGM2F(today, first_startDate, first_endDate, 
						first_SwapFreq, first_SwapBasis, first_SwapRefRate,
						&(volBondCpn->firstSwap));
	if(err)
	{
//		err = "Error in Fill_VolBondCoupon : Cannot fill first swap";
		goto FREE_RETURN;
	}

	err = Fill_SWAP_FOR_LGM2F(today, second_startDate, second_endDate, 
						second_SwapFreq, second_SwapBasis, second_SwapRefRate,
						&(volBondCpn->secondSwap));
	if(err)
	{
//		err = "Error in Fill_VolBondCoupon : Cannot fill second swap";
		goto FREE_RETURN;
	}


FREE_RETURN:
	if(err)
	{
		Free_VolBondCoupon(volBondCpn);
	}

	return err;
}


Err Free_VolBondCoupon(volBondCoupon *volBondCpn)
{
	Err err = NULL;

	err = Free_SWAP_FOR_LGM2F(&(volBondCpn->firstSwap));
	err = Free_SWAP_FOR_LGM2F(&(volBondCpn->secondSwap));
	
	return err;
}


Err Free_VolBondStruct(volBondStruct *volBond)
{
	Err err = NULL;
	long i;

	if(volBond->VolBondCpn)
	{
		for(i=0;i<volBond->NbOfCoupons;++i)
		{
			err = Free_VolBondCoupon(&(volBond->VolBondCpn[i]));
		}

		free(volBond->VolBondCpn);
		volBond->VolBondCpn=NULL;
	}

	return err;
}


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
					   volBondStruct *volBond)
{
	long i;
	Err err;

	volBond->NbOfCoupons = NbCoupons;

	volBond->VolBondCpn = NULL;
	volBond->VolBondCpn = (volBondCoupon*) calloc(NbCoupons, sizeof(volBondCoupon));
	if (!volBond->VolBondCpn)
	{
		err = "Error in Fill_VolBondStruct : volBondCpn allocation failed ";
		goto FREE_RETURN;
	}

	for(i=0;i<NbCoupons;++i)
	{
		err = Fill_VolBondCoupon(today, first_SwapFreq[i], first_SwapBasis[i], first_SwapRefRate[i], 
					   second_SwapFreq[i], second_SwapBasis[i], second_SwapRefRate[i], 
					   notional[i],
					   first_fixingDate[i], first_startDate[i], 
					   first_endDate[i],
					   second_fixingDate[i], second_startDate[i], 
					   second_endDate[i],
					   payment_Date[i], 
					   alphaC[i], alphaP[i], betaC[i], betaP[i],
					   strikeC[i], strikeP[i],
					   &(volBond->VolBondCpn[i]));
		if (err)
		{
			goto FREE_RETURN;
		}

	}

FREE_RETURN:
	if (err)
	{
		Free_VolBondStruct(volBond);
	}

	return err;
}


Err PayOff_VolBondCoupon(long today, double fact1, double fact2, 
							volBondCoupon *volBondCpn, 
							double CMSAdj, 
							double sabrAlpha,
							double sabrBeta,
							double sabrBetaVol,
							double sabrRho,
							double *Payoff,
							TermStruct* l, char *yc)
{
	Err err=NULL;
	double firstswap;
	double secondswap;

	firstswap = max(0.00000001, SWAP_FOR_LGM2F_SwapRate_From_Factors_Value(fact1, fact2, &(volBondCpn->firstSwap)));

	secondswap = max(0.00000001, SWAP_FOR_LGM2F_SwapRate_From_Factors_Value(fact1, fact2, &(volBondCpn->secondSwap)));

	*Payoff = 	volBondCpn->alphaC * srt_f_optblkschbetastochquick(secondswap * CMSAdj, 
					(volBondCpn->betaC / volBondCpn->alphaC) * firstswap + (volBondCpn->strikeC / volBondCpn->alphaC), 
					(volBondCpn->second_fixingDate - volBondCpn->first_fixingDate)/365.0, 
					sabrBetaVol, 
					sabrAlpha, 
					sabrBeta, 
					sabrRho, 
					1.0, 
					SRT_BETAVOL, SRT_LOGNORMAL, SRT_CALL, PREMIUM)
				
				+volBondCpn->alphaP * srt_f_optblkschbetastochquick(secondswap * CMSAdj, 
					(volBondCpn->betaP / volBondCpn->alphaP) * firstswap + (volBondCpn->strikeP / volBondCpn->alphaP), 
					(volBondCpn->second_fixingDate - volBondCpn->first_fixingDate)/365.0, 
					sabrBetaVol, 
					sabrAlpha, 
					sabrBeta, 
					sabrRho, 
					1.0, 
					SRT_BETAVOL, SRT_LOGNORMAL, SRT_PUT, PREMIUM);

	return err;

}


