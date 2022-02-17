#include "ARMKernel\glob\firsttoinc.h"

/*********************************************************************************/
/*! \class  CreditManager CreditManager_Pricing_1F.cpp "CreditManager_FFT_1F.cpp"
 *  \author:	L. JACQUEL
 *	\version:	1.0
 *	\date:		July 2005
 *	\brief:		Implementation of the 1 Factor Copula Model - Pricing Legs Part:
				whatever Models behind
/***********************************************************************************/


#include "ICMKernel\cair\credit_manager.h"

#include "ICMKernel\util\icm_integrator.h"
#include "ICMKernel\glob\icm_maths.h"

//#include <nags.h>		//	s15abc


void CreditManager::PriceBasket_1F(DoubleVector& Outputs)
{
	double	NPV		=	0.0;
	double	NPVPrem	=	0.0;
	double	NPVDef	=	0.0;
	double	NPVPremATM	=	0.0;
	double	PremATM	=	0.0;

	double	TheCash	=	0.0;
	double	TheAccrued	=	0.0;

	DoubleVector	DefaultLegOutputs;
	DoubleVector	PremiumLegOutputs;

	Outputs.clear();
	DefaultLegOutputs.clear();
	PremiumLegOutputs.clear();

	if (its_PricingLegsType != CBN_PREMIUMLEGONLY)
	{
		if (its_CreditObservationType == CO_LOSSES)
			PriceDefaultLeg_1F_Losses(DefaultLegOutputs);
		else	// NB DEF
			PriceDefaultLeg_1F_NbDef(DefaultLegOutputs);

		// NPV Default Leg
		NPVDef	=	DefaultLegOutputs[0];
	}

	if (its_PricingLegsType != CBN_DEFAULTLEGONLY)
	{
		if (its_CreditObservationType == CO_LOSSES)
			PricePremiumLeg_1F_Losses(PremiumLegOutputs);
		else
			PricePremiumLeg_1F_NbDef(PremiumLegOutputs);

		// NPV Premium Leg
		NPVPrem	=	PremiumLegOutputs[0];		
		// NPV ATM Premium Leg
		NPVPremATM	=	PremiumLegOutputs[1];
		// without Notio
		PremATM		=	PremiumLegOutputs[2];

		// Cash
		TheCash	=	PremiumLegOutputs[3];
		// Accrued
		TheAccrued	=	PremiumLegOutputs[4];
	}
			
	// -------------------------------------------------------
	NPV	=	PremNPVFlag * NPVPrem + DefNPVFlag * NPVDef;
	
	Outputs.resize(5);
	
	Outputs[0]	=	NPV;
	Outputs[1]	=	NPVDef;
	Outputs[2]	=	NPVPrem;
	Outputs[3]	=	NPVPremATM;
	Outputs[4]	=	PremATM;
}


// -----------------------------------------------------------------
// 1F Losses Computation: DEFAULT LEG - LOSSES
// -----------------------------------------------------------------

void CreditManager::PriceDefaultLeg_1F_Losses(DoubleVector& Outputs)
{
	double	NPVDef;
	double	LMin, LMax;
	double	TMin, TMax;

	double	ToPay;
	double	TheDF;

	double	T1, T2, TCenter;
	double	LossT1, LossT2;

	RelativeDate	NextDate;

	// Split in two parts
	// Default Leg and Premium Leg Pricing
	NPVDef = 0.0;

	int	i;

	LMin = its_DL_LossMin * its_BasketNotional;
	LMax = its_DL_LossMax * its_BasketNotional;

	TMin	=	its_DL_CreditWindowLow;
	TMax	=	its_DL_CreditWindowUp;

	i=0;
	
	ToPay = 0.0;
	
	DoubleVector	TheDateToFound;
	TheDateToFound.push_back(0.0);

	DoubleVector::iterator	iterToFound;
//	DoubleVector::iterator	iterFound;
	iterToFound	=	TheDateToFound.begin();

	CreditLossComputation	TheCreditLossComputation;

	TheCreditLossComputation	=	(its_CreditModelType == CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F) ? CLC_EXPECTATION_INTERPOLATION : CLC_EXPECTATION;

	switch (its_DL_PaymentType)
	{
		case CEP_ATDEFAULTDATE:

			T1 = TMin; 

			while (T1 < 0.0)
				T1 +=	(RelativeDate)its_Time_Step_Prorata_1F;

			// Compute Loss at T1
			ComputeLossesBeforeT_1F(T1, TheCreditLossComputation, LMin, LMax, LossT1);
			
			if (T1)
			{
				TCenter = 0.5*T1;
				TheDF = its_DF[TCenter];
				NPVDef += LossT1 * TheDF;
			}

			// then according to the Time Step go through the interval
			while (T1 < TMax)
			{
				T2 = T1 + (RelativeDate)its_Time_Step_Prorata_1F; 
				if (T2 > TMax) T2 = TMax;
				TCenter = 0.5*(T1+T2);

				TheDF = its_DF[TCenter];
				
				// Next Point
				ComputeLossesBeforeT_1F(T2, TheCreditLossComputation, LMin, LMax, LossT2);

				NPVDef += (LossT2-LossT1) * TheDF;

				T1 = T2;
				LossT1 = LossT2;
			}

			
			break;

		case CEP_ATNEXTCREDITDATE:
		case CEP_ATNEXTPAYMENTDATE:
	
			T1 = TMin;

			while (T1 < TMax)
			{
				TheDateToFound[0]	=	T1+1;
				NextDate	=	T1+1;

/*				if (its_DL_PaymentType == CEP_ATNEXTCREDITDATE)
					iterFound	=	find_if(its_PL_CreditWindowUps.begin(), its_PL_CreditWindowUps.end(), GreaterThan(T1+1));
				if (its_DL_PaymentType == CEP_ATNEXTPAYMENTDATE)	
					iterFound	=	upper_bound(its_PL_CreditWindowUps.begin(), its_PL_CreditWindowUps.end(), iterToFound);
*/			
//				NextDate	= *iterFound;

				T2 = NextDate;
				if (T2 > TMax)
					T2 = TMax;

				// Flag 0 means, Computes Expectation and Not Proba
				ComputeLossesBeforeT_1F(T1, TheCreditLossComputation, LMin, LMax, LossT1);
				ComputeLossesBeforeT_1F(T2, TheCreditLossComputation, LMin, LMax, LossT2);

				TheDF	=	its_DF[(int)NextDate];

				ToPay += (LossT2-LossT1) * TheDF;

				T1 = T2;
			}

			NPVDef = ToPay;

			break;

		case CEP_ATFIXEDDATE:

			// Flag 0 means, Computes Expectation and Not Proba
			ComputeLossesBeforeT_1F(TMin, TheCreditLossComputation, LMin, LMax, LossT1);
			ComputeLossesBeforeT_1F(TMax, TheCreditLossComputation, LMin, LMax, LossT2);

			TheDF	=	TheDFAtAFixedDate;
			NPVDef	=	(LossT2-LossT1) * TheDF;

			break;

		default:
			break;
	}
		
	Outputs.clear();
	Outputs.push_back(NPVDef);

}


// -----------------------------------------------------------------
// 1F Losses Computation: DEFAULT LEG - NB DEF
// -----------------------------------------------------------------

void CreditManager::PriceDefaultLeg_1F_NbDef(DoubleVector& Outputs)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"PriceDefaultLeg_1F_NbDef not implemented yet!");
}


// -----------------------------------------------------------------
// 1F Losses Computation: PREMIUM LEG - LOSSES
// -----------------------------------------------------------------

void CreditManager::PricePremiumLeg_1F_Losses(DoubleVector& Outputs)
{
	int	i;
	double	NPVPrem, ProbT1, ProbT2, Int, T1, T2, ProbaPremiumLeg1Fi;
	double	TCenter;

	double	TheDF;
	double	TheValue, NPVPremATM, PremATM;
	double	LMin, LMax;

	NPVPrem		=	0.0;
	NPVPremATM	=	0.0;
	PremATM		=	0.0;

	for (i=0;i<its_PL_NbFlows;i++)
	{
		// Flows are paid if at the Credit Fixing Date, the event has occurred
		// using 'static' variables
		its_PL_StartDate	=	its_PL_StartDates[i];
		its_PL_EndDate		=	its_PL_EndDates[i];

		// Ratio
		its_PL_Ratio			=	its_PL_Ratios[i];
			
		// LOSS OBSERVATIONS
		its_PL_Notio			=	its_PL_Notios[i];
		its_PL_Spread			=	its_PL_Spreads[i];
		
		// CREDIT OBSERVATIONS
		its_PL_CreditWindowLow	=	its_PL_CreditWindowLows[i];
		its_PL_CreditWindowUp	=	its_PL_CreditWindowUps[i];

		its_PL_PaymentDate		=	its_PL_PaymentDates[i];

		// Credit Basket Data
		its_CurrentIndex		=	i;

		// LOSS OBSERVATIONS
		its_PL_PaymentDate		=	its_PL_PaymentDates[i];
		
		its_PL_CreditFlag		=	its_PL_CreditFlags[i];
		its_PL_LossMin			=	its_PL_LossMins[i];
		its_PL_LossMax			=	its_PL_LossMaxs[i];
		its_PL_CreditSpreadCap	=	its_PL_CreditSpreadCaps[i];
		its_PL_Redemption		=	its_PL_Redemptions[i];

		// I need notional amounts
		LMin	=	its_PL_LossMin * its_BasketNotional;
		LMax	=	its_PL_LossMax * its_BasketNotional;

		// Discount Factor
		TheDF	=	its_PL_PaymentDates_DF[i];

		 // 1Factor Pricing: For the moment, please, keep NbMin=0 for Premium Leg (Type=0), or USE MC Method. Think if it makes sense to have NbMin > 0"

		if (its_PL_EndDate <= 0.0)
			continue;

		// Non prorata part of the Premium Leg
		// TO BE REVIEWED IN ORDER TO TAKE INTO ACCOUNT 'its_PL_CreditWindowLow'
		ProbaPremiumLeg_Credit_i_1F(its_PL_CreditWindowUp, its_PL_CreditFlag, 0, 0, LMin, LMax, ProbaPremiumLeg1Fi);

		// to be optimized
		TheValue	=	its_PL_Ratio * TheDF * ProbaPremiumLeg1Fi;		
		PremATM		+=	TheValue;

		TheValue	*=	its_PL_Notio;
		NPVPremATM	+=	TheValue;
		
		TheValue	*=	its_PL_Spread;
		NPVPrem		+=	TheValue;

		// Prorata part of the Premium Leg
		if (its_CreditPremiumLegAccrued == CAP_PRORATA)
		{
			Int	= 0.0;							
			T1  = its_PL_StartDate; 

			ProbaPremiumLeg_Credit_i_1F(T1, its_PL_CreditFlag, 0.0, 0.0, its_PL_LossMin, its_PL_LossMax, ProbT1);
			
			while (T1 < its_PL_EndDate)
			{

				T2 = T1 + its_Time_Step_Prorata_1F; 
				if (T2 > its_PL_EndDate) T2 = its_PL_EndDate; 

				TCenter = 0.5*(T1+T2);

				// flag to be changed?
				if (its_PL_PaymentType == CEP_ATDEFAULTDATE)
					TheDF	=	its_DF[(int)TCenter];

				ProbaPremiumLeg_Credit_i_1F(T2, its_PL_CreditFlag, 0.0, 0.0, its_PL_LossMin, its_PL_LossMax, ProbT2);

				Int += (TCenter - its_PL_StartDate)/(its_PL_EndDate-its_PL_StartDate) * TheDF * (ProbT2-ProbT1);

				T1 = T2;
				ProbT1 = ProbT2;

			}

			TheValue	=	its_PL_Ratio * Int;		
			PremATM		-=	TheValue;

			TheValue	*=	its_PL_Notio;
			NPVPremATM	-=	TheValue;
			
			TheValue	*=	its_PL_Spread;
			NPVPrem		-=	TheValue;
		}
	}

	Outputs.clear();

	// NPV Premium Leg
	Outputs.push_back(NPVPrem);
	// NPV ATM Premium Leg
	Outputs.push_back(NPVPremATM);
	// without Notio
	Outputs.push_back(PremATM);
}


// -----------------------------------------------------------------
// 1F Losses Computation: PREMIUM LEG - NB DEF
// -----------------------------------------------------------------

void CreditManager::PricePremiumLeg_1F_NbDef(DoubleVector& Outputs)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"PricePremiumLeg_1F_NbDef not implemented yet!");
}



// -----------------------------------------------------------------
// 1F Losses Computation: T must be transformed in year fractions

//		PayoffFLag = 0
//							Expectation of [LossMin < Losses < LossMax] at time T
//		PayoffFlag = 1  
//							Proba that [LossMin < Losses < LossMax] at time T

// -----------------------------------------------------------------
void CreditManager::ComputeLossesBeforeT_1F(double T, CreditLossComputation PayOffType, double LossMin, double LossMax, double& Result)
{
	if (T <= 0.)
	{
		Result = 0.0;
		return;
	}

	map<Loss_Key, double>::iterator	iter;
	Loss_Key	current_loss_key(T, LossMin, LossMax);

	// already computed?
	iter	=	its_Loss_Distrib_Map.find(current_loss_key);

	if (iter != its_Loss_Distrib_Map.end())
	{
		Result	=	iter->second;
		return;
	}

	switch (its_CreditModelType)
	{
	case CMT_ANALYTIC_RECURSIVE_1F:
	case CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F:

		ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_1F(MATHTIME(T), PayOffType, LossMin, LossMax, Result);
		break;

	case CMT_ANALYTIC_LARGE_PORTFOLIO_1F:
		
		ComputeLossesBeforeT_LARGE_PORTFOLIO_1F(MATHTIME(T), PayOffType, LossMin, LossMax, Result);
		break;

	case CMT_ANALYTIC_LHP:
		
		ComputeLossesBeforeT_LHP(MATHTIME(T), PayOffType, LossMin, LossMax, Result);
		break;

	case CMT_ANALYTIC_LHP_JPM:
		
		ComputeLossesBeforeT_LHP_JPM(MATHTIME(T), PayOffType, LossMin, LossMax, Result);
		break;

	case CMT_ANALYTIC_LHP_PLUS:
		
		ComputeLossesBeforeT_LHP_PLUS(MATHTIME(T), PayOffType, LossMin, LossMax, Result);
		break;

	case CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F:

		ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F(MATHTIME(T), PayOffType, LossMin, LossMax, Result);
		break;

	}

	its_Loss_Distrib_Map.insert(make_pair(current_loss_key, Result));
}



// --------------------------------------------------------------------------------------
// Analytical calculation of the Call on Loss : E[(L-K)1{L>K} ]
// --------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------
// RECURSIVE ALGORITHM
// --------------------------------------------------------------------------------------

double CreditManager::LossCall(int NbCredits, double Strike, double MaxLoss, double ExpectedLoss)
{
	double	Result;


	if (Strike <= 0.0)
		return	ExpectedLoss - Strike;
	else if (MaxLoss < Strike)
		return 0.0;
	else if (NbCredits == 1)
		return  (ExpectedLoss - Strike * its_Current_DefProb[0]);
	else
	{
		double	CurrentLoss	=	its_Input_Losses[NbCredits - 1];
		double	p_def	=	its_Current_DefProb[NbCredits - 1];

		// Last Credit has defaulted with proba p_def
		Result = p_def * LossCall(NbCredits - 1, Strike - CurrentLoss, MaxLoss - CurrentLoss, ExpectedLoss - p_def * CurrentLoss);
		
		// Last Credit has NOT defaulted
		Result += (1.0 - p_def) * LossCall(NbCredits - 1, Strike, MaxLoss - CurrentLoss, ExpectedLoss - p_def * CurrentLoss);
		
		return Result;
	}
}


// --------------------------------------------------------------------------------------
// Analytical calculation of the Digit on Loss : E[1{L>K} ]
// --------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------
// RECURSIVE ALGORITHM
// --------------------------------------------------------------------------------------

double CreditManager::LossDigit(int NbCredits, double Strike, double MaxLoss)
{
	double	Result;

	if (Strike <= 0.0)
		return 1.0;
	else if (MaxLoss < Strike)
		return 0.0;
	else if (NbCredits == 1)
		return its_Current_DefProb[0];
	else
	{
		double	CurrentLoss	=	its_Input_Losses[NbCredits - 1];
		double	p_def	=	its_Current_DefProb[NbCredits - 1];

		// Last Credit has defaulted with proba p_def
		Result = p_def * LossDigit(NbCredits - 1, Strike - CurrentLoss, MaxLoss - CurrentLoss);
		
		// Last Credit has NOT defaulted
		Result += (1.0 - p_def) * LossDigit(NbCredits - 1, Strike, MaxLoss - CurrentLoss);
		
		return Result;
	}
}


void CreditManager::ComputeLossesBeforeT_FFT_1F(double T, CreditLossComputation PayOffType, double LossMin, double LossMax, double& Result)
{
	ICMTHROW(ERR_INVALID_DATA,"ComputeLossesBeforeT_FFT_1F not implemented yet!");
}


// --------------------------------------------------------------------------------------
// Calculate Probability of Coupon Payment during period i of Premium Leg
// --------------------------------------------------------------------------------------

void CreditManager::ProbaPremiumLeg_Credit_i_1F(double T, int Typei, int NbMax, int NbMin, double LMin, double LMax, double& Result)
{
	double	ProbT, LossT;
	CreditLossComputation	TheCreditLossComputation;

	if (Typei == 1)
	{
		// Condition on Losses
		TheCreditLossComputation	=	(its_CreditModelType == CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F) ? CLC_PROBABILITY_INTERPOLATION : CLC_PROBABILITY;

		ComputeLossesBeforeT_1F(T, TheCreditLossComputation, LMin, LMax, ProbT);
		Result = ProbT;

	}
	else if (Typei == 2)
	{
		// Decreasing Losses		
		TheCreditLossComputation	=	(its_CreditModelType == CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F) ? CLC_EXPECTATION_INTERPOLATION : CLC_EXPECTATION;

		ComputeLossesBeforeT_1F(T, TheCreditLossComputation, LMin, LMax, LossT);

		Result = (1.0 - LossT/(LMax-LMin));
	}
	else if (Typei == 0)
	{
		// Guaranteed
		Result = 1.0;;
	}

}