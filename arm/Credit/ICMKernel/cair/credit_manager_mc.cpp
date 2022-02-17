#include "ARMKernel\glob\firsttoinc.h"


#include "ICMKernel\cair\credit_manager.h"
#include "ICMKernel\glob\icm_maths.h"
#include "ICMKernel\cair\credit_manager.h"


// ----------------------------------------------------------------------------
// Generation of default times with Betas
// ----------------------------------------------------------------------------

void CreditManager::PriceBasketForThisSimulation(DoubleVector& Outputs)
{
	double	NPV		=	0.0;
	double	NPVPrem	=	0.0;
	double	NPVDef	=	0.0;
	double	NPVPremATM	=	0.0;
	double	PremATM	=	0.0;

	DoubleVector	DefaultLegOutputs;
	DoubleVector	PremiumLegOutputs;

	Outputs.clear();
	DefaultLegOutputs.clear();
	PremiumLegOutputs.clear();

	if (its_PricingLegsType != CBN_PREMIUMLEGONLY)
	{
		if (its_CreditObservationType == CO_LOSSES)
			PriceDefaultLegForThisSimulation_Losses(DefaultLegOutputs);
		else	// NB DEF
			PriceDefaultLegForThisSimulation_NbDef(DefaultLegOutputs);

		// NPV Default Leg
		NPVDef	=	DefaultLegOutputs[0];
	}

	if (its_PricingLegsType != CBN_DEFAULTLEGONLY)
	{
		if (IsCDOStandard())
			ComputeLossesAmountBeforeCreditWindowsDates();
		else
			ComputeLossesAmountBeforeCreditWindowsDates_CDO_SQUARE();

		if (its_CreditObservationType == CO_LOSSES)
		{
			PricePremiumLegForThisSimulation_Losses(PremiumLegOutputs);
		}
		else
			PricePremiumLegForThisSimulation_NbDef(PremiumLegOutputs);

		// NPV Premium Leg
		NPVPrem	=	PremiumLegOutputs[0];		
		// NPV ATM Premium Leg
		NPVPremATM	=	PremiumLegOutputs[1];
		// without Notio
		PremATM		=	PremiumLegOutputs[2];
	}

	// -------------------------------------------------------
	// VARIANCE REDUCTION
	// -------------------------------------------------------

	// only available for LOSSES
	double	TheValue;
	double	The_LossInPct;		// in Pct of Basket Notional
	
	if (its_MC_Variance_Reduction != CMCVR_NONE)
	{
		if (its_CreditObservationType == CO_LOSSES)
		{
			TheValue	=	1.0;					
			
			switch (its_MC_Variance_Reduction)
			{					
				case CMCVR_IS_PURE_FACTORS:
					TheValue	=	IS_Individual_FactorShift;
					break;

				case CMCVR_IS_FACTORS:
					TheValue	=	IS_Individual_FactorShift;

				case CMCVR_IS_IDIOSYNCRATIC:
					
					if (IsCDOStandard())
						GetHowMuchLossBeforeADate(MATHTIME(its_MaxCreditDate), The_LossInPct);
					else // if (IsCDOSquare())
						GetHowMuchLossBeforeADate_CDO_SQUARE(MATHTIME(its_MaxCreditDate), The_LossInPct);

	//				GetHowMuchLossBeforeADate(its_IS_Loss_Maturity, The_LossInPct);
					The_LossInPct	/=	TheBasketNotional;
					TheValue		*=	exp(- its_IS_Theta * The_LossInPct);		// scale because it is in %
					TheValue		*=	IS_Common_exp_twists;

	//				NPVPrem		=	TheValue;
					break;

				default:

					break;
			}

			NPVPrem		*=	TheValue;
			NPVPremATM	*=	TheValue;
			PremATM		*=	TheValue;
			NPVDef		*=	TheValue;

		}
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


// ----------------------------------------------------------------------------
// Default Leg Pricing for a Simulation Losses
// ----------------------------------------------------------------------------

void CreditManager::PriceDefaultLegForThisSimulation_Losses(DoubleVector& Outputs)
{
	// ---------------------------------------------------------------------------------
//	if (IsCDOSquare())
//	{
//		PriceDefaultLegForThisSimulation_Losses_CDO_SQUARE(Outputs);
//		return;
//	}
	// ---------------------------------------------------------------------------------

	double	NPVDef;
	double	LMin, LMax;
	double	TMin, TMax;
	double	tau_value;

	double	ToPay;
	double	Li, TotalLoss;
	double	TheDF;

	RelativeDate	NextDate, tau_relative_date;

	int		the_id;

	// Split in two parts
	// Default Leg and Premium Leg Pricing
	NPVDef = 0.0;

	int	i;

	if (IsCDOStandard())
	{
		LMin = its_DL_LossMin * its_BasketNotional;
		LMax = its_DL_LossMax * its_BasketNotional;
	}
	else
	{
		LMin = its_DL_LossMin * itsCDOSquareNotional;
		LMax = its_DL_LossMax * itsCDOSquareNotional;
	}

	TMin	=	its_DL_CreditWindowLow;
	TMax	=	its_DL_CreditWindowUp;

	i=0;

	// I want to find all items between Tmin and Tmax
	set<Tau_Item>::iterator iter;
	Tau_Item	LowTime(0,MATHTIME(TMin+T_EPS),0.0);		// T_EPS, just to avoid a default matching exactly the date TMin
	Tau_Item	UpTime(0,MATHTIME(TMax),0.0);

	for (iter = itsSortedDefaultTimes.lower_bound(LowTime);
				iter != itsSortedDefaultTimes.lower_bound(UpTime);
				++iter)
	{
		// id of Non-Sorted Credit
		the_id	=	iter->id;
		tau_value	=	iter->tau;
	
		ToPay = 0.0;
		
		// ------------------------------------------------------------
		// OBSERVATION LOSSES
		// get Loss Amount Before date T = tau(i) excluded

		if (IsCDOStandard())
		{
			TotalLoss	=	iter->cumul_loss;
			// current Loss
			Li = its_LossesAmount[the_id];		// its_losses * its_notionals
		}
		else
		{
			// if IsCDOSquare()
			// some computation must have been carried out
			// ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE
			TotalLoss	=	iter->prev_cumul_loss_square;

			Li	=	iter->equiv_current_loss_square;
		}
		// ------------------------------------------------------------

		// Check if it lies in the interval
		if (TotalLoss < LMin)
		{
			if (TotalLoss + Li < LMin)
				ToPay = 0.0;
			else if (TotalLoss + Li < LMax)
				ToPay = TotalLoss + Li - LMin;
			else
				ToPay = LMax - LMin;
		}
		else if (TotalLoss < LMax)
		{
			if (TotalLoss + Li < LMax)
				ToPay = Li;
			else
				ToPay = LMax - TotalLoss;
		}
		else
			ToPay = 0.0;

		// Now, get the Discount Factor
		if (ToPay)
		{
			tau_relative_date	=	DAYSTIME(tau_value);

			switch (its_DL_PaymentType)
			{
				case CEP_ATDEFAULTDATE:
					// tau_value is not an integer --> approximation in order to avoid this computation?
					// what approximation, low or up bound?
					TheDF = its_DF[tau_relative_date];
					break;

				case CEP_ATNEXTCREDITDATE:
			
					Get_DFAtNextCreditDate(tau_relative_date, NextDate, TheDF);
					break;

				case CEP_ATNEXTPAYMENTDATE:

					Get_DFAtNextPaymentDate(tau_relative_date, NextDate, TheDF);
					break;

				case CEP_ATFIXEDDATE:

					TheDF	=	TheDFAtAFixedDate;
					break;

				default:
					// throw
					TheDF	=	-1.0;
			}
		
			NPVDef += ToPay * TheDF;
		}	// end ToPay
	}	// end loop Get_NbCredits()

	//Long Short CDO
	NPVDef = MAX(NPVDef,0.);
	Outputs.clear();
	Outputs.push_back(NPVDef);
}



// ----------------------------------------------------------------------------
// Default Leg Pricing for a Simulation NbDef
// ----------------------------------------------------------------------------

void CreditManager::PriceDefaultLegForThisSimulation_NbDef(DoubleVector& Outputs)
{
	// ---------------------------------------------------------------------------------
//	if (IsCDOSquare())
//	{
//		PriceDefaultLegForThisSimulation_NbDef_CDO_SQUARE(Outputs);
//		return;
//	}
	// ---------------------------------------------------------------------------------

	double	NPVDef;
	double	TMin, TMax;
	double	tau_value;

	double	ToPay;
	double	Li, TotalNbDef;
	double	TheDF;

	RelativeDate	NextDate, tau_relative_date;

	int		the_id;

	// Split in two parts
	// Default Leg and Premium Leg Pricing
	NPVDef = 0.0;

	TMin	=	its_DL_CreditWindowLow;
	TMax	=	its_DL_CreditWindowUp;

	// I want to find all items between Tmin and Tmax
	set<Tau_Item>::iterator iter;
	Tau_Item	LowTime(0,MATHTIME(TMin),0.0);
	Tau_Item	UpTime(0,MATHTIME(TMax),0.0);

	for (iter = itsSortedDefaultTimes.lower_bound(LowTime);
				iter != itsSortedDefaultTimes.lower_bound(UpTime);
				++iter)
	{
		// id of Non-Sorted Credit
		the_id	=	iter->id;
		tau_value	=	iter->tau;
	
		ToPay = 0.0;
		
		// OBSERVATION NB DEF
		// get Nb of Defaults before date T = tau(i) excluded
		if (IsCDOStandard())
			TotalNbDef	=	iter->cumul_nbdef + 1;
		else // (IsCDOSquare())
			// some computation must have been carried out
			// ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE
			TotalNbDef	=	iter->cumul_nbdef_square + 1;

		// current Loss
		Li = its_LossesAmount[the_id];		// its_losses * its_notionals

		// Check if it lies in the interval
		if ((TotalNbDef >= its_DL_NbDefMin) && (TotalNbDef <= its_DL_NbDefMax))
		{
			// Think over Fixed Recovery
			ToPay	=	Li;

			// Now, get the Discount Factor
			tau_relative_date	=	DAYSTIME(tau_value);

			switch (its_DL_PaymentType)
			{
				case CEP_ATDEFAULTDATE:
					// tau_value is not an integer --> approximation in order to avoid this computation?
					// what approximation, low or up bound?
					TheDF = its_DF[tau_relative_date];
					break;

				case CEP_ATNEXTCREDITDATE:
			
					Get_DFAtNextCreditDate(tau_relative_date, NextDate, TheDF);
					break;

				case CEP_ATNEXTPAYMENTDATE:

					Get_DFAtNextPaymentDate(tau_relative_date, NextDate, TheDF);
					break;

				case CEP_ATFIXEDDATE:

					TheDF	=	TheDFAtAFixedDate;
					break;

				default:
					// throw
					TheDF	=	-1.0;
			}
		
			NPVDef += ToPay * TheDF;
		}	// end ToPay
	}	// end loop Get_NbCredits()

	//Long Short CDO
	NPVDef = MAX(NPVDef,0.);

	Outputs.clear();
	Outputs.push_back(NPVDef);
}


// ----------------------------------------------------------------------------
// Premium Leg Pricing for a Simulation LOSSES
// ----------------------------------------------------------------------------

void CreditManager::PricePremiumLegForThisSimulation_Losses(DoubleVector& Outputs)
{
	// ---------------------------------------------------------------------------------
//	if (IsCDOSquare())
//	{
//		PricePremiumLegForThisSimulation_Losses_CDO_SQUARE(Outputs);
//		return;
//	}

	// ---------------------------------------------------------------------------------

	int	i;

	double	NPVPrem;
	double	NPVPremATMMargin;
	double	PremATMMargin;
	double	NPVFlow;
	double	NPVFlowForATMMargin;
	double	NPVFlow0, NPVFlow1;
	double	NPVFlowForATMMargin0, NPVFlowForATMMargin1;
	double	FlowForATMMargin;
	double	FlowForATMMargin0, FlowForATMMargin1;
	double	TheDF;
	double	Local_ZEPS = 1.e-6;

	DoubleVector	FlowOutputs;
	FlowOutputs.clear();

	int		NbDefCreditStart, NbDefCreditEnd;
	RelativeDate	Tau_idef, Last_idef;

	set<Tau_Item>::iterator		iter_start, iter_end;

	// -------------------------
	// Default Times are Sorted
	// -------------------------
	NPVPrem = 0.0;
	NPVPremATMMargin	=	0.0;
	PremATMMargin	=	0.0;

	// -------------------------
	// specific Credit TARN
	// -------------------------
	its_Cumulative_Coupon	=	0.0;
	its_Toggle_Redemption	=	false;
	its_Final_Redemption	=	false;

	int j;

	j=0;

	for (i=0;i<its_PL_NbFlows;i++)
	{
		its_iFlow	=	i;
		// Flows are paid if at the Credit Fixing Date, the event has occurred
		// using 'static' variables
		its_PL_StartDate	=	its_PL_StartDates[i];
		its_PL_EndDate		=	its_PL_EndDates[i];

		if (its_PL_EndDate <= 0.0)
			continue;

		its_PL_CreditWindowLow	=	its_PL_CreditWindowLows[i];
		its_PL_CreditWindowUp	=	its_PL_CreditWindowUps[i];

		its_PL_PaymentDate		=	its_PL_PaymentDates[i];

		// Ration
		its_PL_Ratio			=	its_PL_Ratios[i];
			
		// Credit Basket Data
		its_CurrentIndex		=	i;

		// LOSS OBSERVATIONS
		its_PL_Notio			=	its_PL_Notios[i];
		its_PL_Spread			=	its_PL_Spreads[i];
		its_PL_PaymentDate		=	its_PL_PaymentDates[i];
		
		its_PL_CreditFlag		=	its_PL_CreditFlags[i];
		its_PL_LossMin			=	its_PL_LossMins[i];
		its_PL_LossMax			=	its_PL_LossMaxs[i];
		its_PL_CreditSpreadCap	=	its_PL_CreditSpreadCaps[i];
		its_PL_Redemption			=	its_PL_Redemptions[i];

		// Discount Factor
		TheDF	=	its_PL_PaymentDates_DF[i];
		
		// PRORATA or Not

		if ((its_CreditPremiumLegAccrued == CAP_NON_PRORATA) || 
				((its_CreditPremiumLegAccrued == CAP_PRORATA) && (its_PL_CreditFlag != CPLT_OUTSTANDING)))
		{
			Price_A_Flow_ForPremiumLegForThisSimulation_Losses(FlowOutputs);
			
			// NPV Flow
			NPVFlow	=	FlowOutputs[0];
			// NPV Flow for ATM Margin
			NPVFlowForATMMargin	=	FlowOutputs[1];
			// Flow for ATM Margin
			FlowForATMMargin	=	FlowOutputs[2];

			NPVPrem				+=	TheDF * NPVFlow;
			NPVPremATMMargin	+=	TheDF * NPVFlowForATMMargin;
			PremATMMargin		+=	TheDF * FlowForATMMargin;

		}
		else	// CAP_PRORATA --> means that I consider [Begin ; End] interval
		{		// only matters with OUTSTANDING FLAG

			// --------------------------------------------------------------------------------------------------
			// I want to know if some defaults occured between [its_PL_CreditWindowLow ; its_PL_CreditWindowUp]

			double	mt_PL_StartDate	=	MATHTIME(its_PL_StartDate);
			double	mt_PL_EndDate	=	MATHTIME(its_PL_EndDate);
			double	mt_Last_idef;

			if (IsCDOStandard())
			{
				// Number of Defaults before Start Date
				GetHowManyDefaultsBeforeADate(mt_PL_StartDate, NbDefCreditStart, iter_start);
				
				// Number of Defaults before End Date
				GetHowManyDefaultsBeforeADate(mt_PL_EndDate, NbDefCreditEnd, iter_end);
			}
			else // if (IsCDOSquare())
			{
				// Number of Defaults before Start Date
				GetHowManyDefaultsBeforeADate_CDO_SQUARE(mt_PL_StartDate, NbDefCreditStart, iter_start);
				
				// Number of Defaults before End Date
				GetHowManyDefaultsBeforeADate_CDO_SQUARE(mt_PL_EndDate, NbDefCreditEnd, iter_end);
			}
			// --------------------------------------------------------------------------------------------------
			
			// So, Flow Evaluation at End Date
			Price_A_Flow_ForPremiumLegForThisSimulation_Losses(mt_PL_EndDate, FlowOutputs);

			NPVFlow1				= FlowOutputs[0];;
			NPVFlowForATMMargin1	= FlowOutputs[1];
			FlowForATMMargin1		= FlowOutputs[2];

			Price_A_Flow_ForPremiumLegForThisSimulation_Losses(FlowOutputs);
			
			// NPV Flow
			NPVFlow	=	FlowOutputs[0];
			// NPV Flow for ATM Margin
			NPVFlowForATMMargin	=	FlowOutputs[1];
			// Flow for ATM Margin
			FlowForATMMargin	=	FlowOutputs[2];

			if (iter_start != iter_end)
			{			
				// DF at PaymentDate(i) is already computed
				// while we are before EndDate add the Prorata Flow!

				Last_idef	=	its_PL_EndDate;
				mt_Last_idef	=	MATHTIME(its_PL_EndDate);

				// let the last sub-interval
				while (iter_end != iter_start)
				{
					--iter_end;

					Tau_idef = iter_end->tau;
					
					if ((Tau_idef > mt_PL_StartDate) && (Tau_idef < mt_PL_EndDate))
					{
						Price_A_Flow_ForPremiumLegForThisSimulation_Losses(Tau_idef + Local_ZEPS, FlowOutputs);
						
						// NPV Flow
						NPVFlow0				=	FlowOutputs[0];
						// NPV Flow for ATM Margin
						NPVFlowForATMMargin0	=	FlowOutputs[1];
						// Flow for ATM Margin
						FlowForATMMargin0	=	FlowOutputs[2];

						NPVPrem				+= (NPVFlow0 - NPVFlow1) * (mt_Last_idef - Tau_idef) / (its_PL_EndDate - its_PL_StartDate) * TheDF;
						NPVPremATMMargin	+= (NPVFlowForATMMargin0 - NPVFlowForATMMargin1)*(mt_Last_idef - Tau_idef)/(its_PL_EndDate - its_PL_StartDate) * TheDF;
						PremATMMargin		+= (FlowForATMMargin0 - FlowForATMMargin1)*(mt_Last_idef - Tau_idef)/(its_PL_EndDate - its_PL_StartDate) * TheDF;
						
						NPVFlow1				= NPVFlow0;
						NPVFlowForATMMargin1	= NPVFlowForATMMargin0;
						FlowForATMMargin1		= FlowForATMMargin0;

						// Next DF
						if (its_PL_PaymentType == CEP_ATDEFAULTDATE)
							TheDF	=	its_DF[(int)(Tau_idef*365)];	// Tau_idef must be a RelativeDate

						mt_Last_idef	=	Tau_idef;
					}

				}

				if (mt_Last_idef != mt_PL_StartDate)	// within DB_TOL ?
				{
					// final sub-interval: [its_PL_StartDate ; tau]
					Price_A_Flow_ForPremiumLegForThisSimulation_Losses(mt_PL_StartDate, FlowOutputs);

					// NPV Flow
					NPVFlow0				=	FlowOutputs[0];
					// NPV Flow for ATM Margin
					NPVFlowForATMMargin0	=	FlowOutputs[1];
					// Flow for ATM Margin
					FlowForATMMargin0		=	FlowOutputs[2];

					NPVPrem				+= (NPVFlow0 - NPVFlow1) * (mt_Last_idef - mt_PL_StartDate) / (mt_PL_EndDate - mt_PL_StartDate) * TheDF;
					NPVPremATMMargin	+= (NPVFlowForATMMargin0 - NPVFlowForATMMargin1)*(mt_Last_idef - mt_PL_StartDate)/(mt_PL_EndDate - mt_PL_StartDate) * TheDF;
					PremATMMargin		+= (FlowForATMMargin0 - FlowForATMMargin1)*(mt_Last_idef - mt_PL_StartDate)/(mt_PL_EndDate - mt_PL_StartDate) * TheDF;
					
					NPVFlow1				= NPVFlow0;
					NPVFlowForATMMargin1	= NPVFlowForATMMargin0;
					FlowForATMMargin1		= FlowForATMMargin0;
				}
			}

			NPVPrem				+= TheDF * NPVFlow1;
			NPVPremATMMargin	+= TheDF * NPVFlowForATMMargin1;
			PremATMMargin		+= TheDF * FlowForATMMargin1;
		}

		// early exit
		if (its_Final_Redemption)
			i	=	its_PL_NbFlows;

//		tmpPrem[its_iFlow]	+=	NPVPrem;
	}

	Outputs.clear();
	Outputs.push_back(NPVPrem);
	Outputs.push_back(NPVPremATMMargin);
	Outputs.push_back(PremATMMargin);

}


// ----------------------------------------------------------------------------
// A Flow of Premium Leg Pricing for a Simulation
// ----------------------------------------------------------------------------

void CreditManager::Price_A_Flow_ForPremiumLegForThisSimulation_Losses(DoubleVector& Outputs)
{
	double	NPVFlow;
	double	NPVFlowForATMMargin;
	double	CurrentLoss, CurrentLossLow, CurrentLossUp;
	double	Coupon;

	// ---------------------------------------------------
	// I need to know the total amount of loss before the date CreditWindowUp
	if (its_PL_CreditFlag != CPLT_GUARANTEED)
	{
		if (IsCDOStandard())
		{
			GetHowMuchLossBeforeACreditObservation(0, CurrentLossLow);		// the date is its_PL_CreditWindowLow
			if (CurrentLossLow != 0.0)
				ICMLOG("Loss Low different from 0.0 as expected!" << CurrentLossLow << " - for simulation id: " << its_SimulId << ". First Default Time is: " << (itsSortedDefaultTimes.begin())->tau);

			GetHowMuchLossBeforeACreditObservation(1, CurrentLossUp);		// the date is its_PL_CreditWindowUp
		}
		else // if (IsCDOSquare())
		{
			GetHowMuchLossBeforeACreditObservation_CDO_SQUARE(0, CurrentLossLow);		// the date is its_PL_CreditWindowLow
			if (CurrentLossLow != 0.0)
				ICMLOG("Loss Low different from 0.0 as expected!" << CurrentLossLow << " - for simulation id: " << its_SimulId << ". First Default Time is: " << (itsSortedDefaultTimes.begin())->tau);

			GetHowMuchLossBeforeACreditObservation_CDO_SQUARE(1, CurrentLossUp);		// the date is its_PL_CreditWindowUp
		}

	// ----------------------------------------------------------------------------
//	if ((its_fOut = fopen("c:\\test\\mc_Price_A_Flow_ForPremiumLegForThisSimulation_Losses.txt", "a+")) == NULL) return;
//	fprintf(its_fOut, "simul:\t\t%u\t\tFlow:\t\t%u\t\tLossLow:\t\t%.lf\t\tLossUp:\t\t%.lf ----------------- \n", its_SimulId, its_iFlow, CurrentLossLow, CurrentLossUp);
	// ----------------------------------------------------------------------------
//	fclose(its_fOut);
	
		// resize
		CurrentLossLow	/= TheBasketNotional;
		CurrentLossUp	/= TheBasketNotional;

		CurrentLoss	=	CurrentLossUp - CurrentLossLow;
	}
	// ---------------------------------------------------
	
	// ---------------------------------------------------
	// SWITCH ACCORDING TO TYPES SELECTION
	NPVFlow = 0.0;
	NPVFlowForATMMargin	=	0.0;
	// ---------------------------------------------------

	switch (its_PL_CreditFlag)
	{
		case CPLT_GUARANTEED:
				NPVFlowForATMMargin	=	its_PL_Notio * its_PL_Ratio;
				NPVFlow = NPVFlowForATMMargin * its_PL_Spread;

			break;

		case CPLT_NOTIONAL:

			if ((CurrentLoss >= its_PL_LossMin) && (CurrentLoss <= its_PL_LossMax))
			{
				NPVFlowForATMMargin	=	its_PL_Notio * its_PL_Ratio;
				NPVFlow = NPVFlowForATMMargin * its_PL_Spread;
			}
			
			break;

		case CPLT_OUTSTANDING:
/*			
			NPVFlowForATMMargin	=	its_PL_Notio * its_PL_Ratio;
			NPVFlowForATMMargin	*=	(FMAX(its_PL_LossMax - CurrentLoss, 0.0) - FMAX(its_PL_LossMin - CurrentLoss, 0.0))/(its_PL_LossMax-its_PL_LossMin);
			NPVFlow = NPVFlowForATMMargin * its_PL_Spread;
			
			break;

		case CPLT_CUMULCOUPON:
*/			
			NPVFlowForATMMargin	=	its_PL_Notio * its_PL_Ratio;
			NPVFlow = NPVFlowForATMMargin;
			Coupon	=	(FMAX(its_PL_LossMax - CurrentLoss, 0.0) - FMAX(its_PL_LossMin - CurrentLoss, 0.0))/(its_PL_LossMax-its_PL_LossMin);
			NPVFlowForATMMargin	*=	Coupon;
			its_Last_FlowForATMMargin	=	Coupon;

			Coupon	*=	its_PL_Spread;
			
			its_Cumulative_Coupon	+=	Coupon;
			its_Last_Coupon			=	Coupon;

			NPVFlow *= Coupon;
			
			break;

		case CPLT_CUMULREDEMPTION:
			
			NPVFlowForATMMargin	=	its_PL_Notio * its_PL_Ratio;
			if (!its_Toggle_Redemption)
			{
				// first passage
				its_Toggle_Redemption	=	true;
				Coupon	=	(FMAX(its_PL_LossMax - CurrentLoss, 0.0) - FMAX(its_PL_LossMin - CurrentLoss, 0.0))/(its_PL_LossMax-its_PL_LossMin);
				NPVFlowForATMMargin	*=	Coupon;
				NPVFlow = NPVFlowForATMMargin	*	its_PL_Spread;
				
				its_Last_FlowForATMMargin	=	Coupon;
				its_Last_Coupon	=	Coupon	*	its_PL_Spread;
			}
			else
			{
				Coupon	=	its_Last_Coupon;
				NPVFlow = NPVFlowForATMMargin	*	Coupon;
				NPVFlowForATMMargin	*=	its_Last_FlowForATMMargin;
			}

			its_Cumulative_Coupon	+=	its_Last_Coupon;
			
			if ((its_Cumulative_Coupon >= its_PL_CreditSpreadCap) || (its_CurrentIndex == its_PL_NbFlows-1))
			{
				// final redemption
				NPVFlow	+=	its_PL_Redemption	*	its_PL_Notio;
				NPVFlowForATMMargin	+=	NPVFlow;
				its_Final_Redemption	=	true;
			}
			
			break;
	}

	// ---------------------------------------------------
	Outputs.clear();
	Outputs.push_back(NPVFlow);
	Outputs.push_back(NPVFlowForATMMargin);		// keep Notio for spread if amortizing structure
	Outputs.push_back(NPVFlowForATMMargin	/ its_PL_Notio);		// for RPV01
}


// ----------------------------------------------------------------------------
// A Flow of Premium Leg Pricing for a Simulation
// ----------------------------------------------------------------------------

void CreditManager::Price_A_Flow_ForPremiumLegForThisSimulation_Losses(RelativeDate TheDate, DoubleVector& Outputs)
{
	double	NPVFlow;
	double	NPVFlowForATMMargin;
	double	FlowForATMMargin;
	double	CurrentLoss;

	// ---------------------------------------------------
	// I need to know the total amount of loss before the date CreditWindowUp
	if (its_PL_CreditFlag != CPLT_GUARANTEED)
	{
		if (IsCDOStandard())
			GetHowMuchLossBeforeADate(TheDate, CurrentLoss);
		else // if (IsCDOSquare())
			GetHowMuchLossBeforeADate_CDO_SQUARE(TheDate, CurrentLoss);

		// resize
		CurrentLoss	/= TheBasketNotional;
	}
	// ---------------------------------------------------
	
	// ---------------------------------------------------
	// SWITCH ACCORDING TO TYPES SELECTION
	FlowForATMMargin	=	its_PL_Ratio;
	NPVFlowForATMMargin	=	its_PL_Notio * its_PL_Ratio;
	NPVFlow = 0.0;
	// ---------------------------------------------------

	switch (its_PL_CreditFlag)
	{
		case CPLT_GUARANTEED:
				NPVFlow = NPVFlowForATMMargin * its_PL_Spread;

			break;
		case CPLT_NOTIONAL:

			if ((CurrentLoss >= its_PL_LossMin) && (CurrentLoss <= its_PL_LossMax))
				NPVFlow = NPVFlowForATMMargin * its_PL_Spread;
			
			break;
		case CPLT_OUTSTANDING:
			
			NPVFlowForATMMargin	*=	(FMAX(its_PL_LossMax - CurrentLoss, 0.0) - FMAX(its_PL_LossMin - CurrentLoss, 0.0))/(its_PL_LossMax-its_PL_LossMin);
			NPVFlow = NPVFlowForATMMargin * its_PL_Spread;
			
			break;
	}

	// ---------------------------------------------------
	Outputs.clear();
	Outputs.push_back(NPVFlow);
	Outputs.push_back(NPVFlowForATMMargin);		// keep Notio for spread if amortizing structure
	Outputs.push_back(NPVFlowForATMMargin	/ its_PL_Notio);		// for RPV01

}


// ----------------------------------------------------------------------------
// Premium Leg Pricing for a Simulation NBDEF
// ----------------------------------------------------------------------------

void CreditManager::PricePremiumLegForThisSimulation_NbDef(DoubleVector& Outputs)
{
	// ---------------------------------------------------------------------------------
//	if (IsCDOSquare())
//	{
//		PricePremiumLegForThisSimulation_NbDef_CDO_SQUARE(Outputs);
//		return;
//	}
	// ---------------------------------------------------------------------------------
	int	i;

	double	NPVPrem;
	double	NPVPremATMMargin;
	double	PremATMMargin;
	double	NPVFlow;
	double	NPVFlowForATMMargin;
	double	FlowForATMMargin;
	double	TheDF;
	double	Local_ZEPS = 1.e-6;

	DoubleVector	FlowOutputs;
	FlowOutputs.clear();

	set<Tau_Item>::iterator		iter_start, iter_end;

	// -------------------------
	// Default Times are Sorted
	// -------------------------
	NPVPrem = 0.0;
	NPVPremATMMargin	=	0.0;
	PremATMMargin	=	0.0;
		
	for (i=0;i<its_PL_NbFlows;i++)
	{
		// Flows are paid if at the Credit Fixing Date, the event has occurred
		// using 'static' variables
		its_PL_StartDate	=	its_PL_StartDates[i];
		its_PL_EndDate		=	its_PL_EndDates[i];

		its_PL_CreditWindowLow	=	its_PL_CreditWindowLows[i];
		its_PL_CreditWindowUp	=	its_PL_CreditWindowUps[i];

		its_PL_PaymentDate		=	its_PL_PaymentDates[i];

		// Ration
		its_PL_Ratio			=	its_PL_Ratios[i];
			
		// Credit Basket Data
		its_CurrentIndex		=	i;

		// LOSS OBSERVATIONS
		its_PL_Notio			=	its_PL_Notios[i];
		its_PL_Spread			=	its_PL_Spreads[i];
		its_PL_PaymentDate		=	its_PL_PaymentDates[i];
		its_PL_NbDefMin			=	its_PL_NbDefMins[i];
		its_PL_NbDefMax			=	its_PL_NbDefMaxs[i];
		
		its_PL_CreditFlag		=	its_PL_CreditFlags[i];		// useless?
		its_PL_LossMin			=	its_PL_LossMins[i];
		its_PL_LossMax			=	its_PL_LossMaxs[i];
		its_PL_CreditSpreadCap	=	its_PL_CreditSpreadCaps[i];

		// Discount Factor
		TheDF	=	its_PL_PaymentDates_DF[i];
		
		Price_A_Flow_ForPremiumLegForThisSimulation_NbDef(FlowOutputs);
		
		// NPV Flow
		NPVFlow	=	FlowOutputs[0];
		// NPV Flow for ATM Margin
		NPVFlowForATMMargin	=	FlowOutputs[1];
		// Flow for ATM Margin
		FlowForATMMargin	=	FlowOutputs[2];

		NPVPrem				+=	TheDF * NPVFlow;
		NPVPremATMMargin	+=	TheDF * NPVFlowForATMMargin;
		PremATMMargin		+=	TheDF * FlowForATMMargin;

	}

	Outputs.push_back(NPVPrem);
	Outputs.push_back(NPVPremATMMargin);
	Outputs.push_back(PremATMMargin);
}


// ----------------------------------------------------------------------------
// A Flow of Premium Leg Pricing for a Simulation NBDEF
// ----------------------------------------------------------------------------

void CreditManager::Price_A_Flow_ForPremiumLegForThisSimulation_NbDef(DoubleVector& Outputs)
{
	int		NbDefLow;
	int		NbDefUp;
	int		NbDef, NbPrevDef;

	double	spread, TheRatio;

	double	NPVFlow;
	double	NPVFlowForATMMargin;
	double	FlowForATMMargin;

	RelativeDate	tau_i;

	set<Tau_Item>::iterator	iter_low;
	set<Tau_Item>::iterator	iter_up;
	set<Tau_Item>::iterator	iter_prev;

	NPVFlow = 0.0;
	NPVFlowForATMMargin	=	0.0;
	FlowForATMMargin	=	0.0;

	// ---------------------------------------------------
	if (IsCDOStandard())
	{
		GetHowManyDefaultsBeforeACreditObservation(0, NbDefLow, iter_low);
		GetHowManyDefaultsBeforeACreditObservation(1, NbDefUp, iter_up);
	}
	else // if (IsCDOSquare())
	{
		GetHowManyDefaultsBeforeACreditObservation_CDO_SQUARE(0, NbDefLow, iter_low);
		GetHowManyDefaultsBeforeACreditObservation_CDO_SQUARE(1, NbDefUp, iter_up);
	}

	NbDef	=	NbDefUp - NbDefLow;
	spread = 0.0;
	// ---------------------------------------------------

	// test in range or not
	if ((NbDef >= its_PL_NbDefMin) && (NbDef <= its_PL_NbDefMax))
	{
		spread	=	its_PL_Spread;
		TheRatio	=	1.0;
	}
	else
	{
		// I have to pay if the Number of Defaults is not in the interval at date 
		//	Observation Date T(i), but was at Observation Date T(i-1)
		if (its_CreditPremiumLegAccrued == CAP_PRORATA)
		{
			// Correct Schedule!
			if (IsCDOStandard())
				GetHowManyDefaultsBeforeADate(its_PL_StartDate, NbPrevDef, iter_prev);
			else // if (IsCDOSquare())
				GetHowManyDefaultsBeforeADate_CDO_SQUARE(its_PL_StartDate, NbPrevDef, iter_prev);

			// get the previous number of defaults
			if ((NbPrevDef >= its_PL_NbDefMin) && (NbPrevDef <= its_PL_NbDefMax))
			{
				spread = its_PL_Spread;
				// I pay: [tau(1) - T(i-1)] * m(i)

				// I get the last relevant Default Time
//				--iter_up;	// test on PrevDef
				
				tau_i	=	iter_prev->tau;
				// Ratio
				TheRatio	=	(tau_i - its_PL_StartDate) / (its_PL_EndDate - its_PL_StartDate);				
			}
		}
	}

	if (spread)
	{
		FlowForATMMargin	=	its_PL_Ratio;
		NPVFlowForATMMargin	=	its_PL_Notio * its_PL_Ratio;
		
		NPVFlow = NPVFlowForATMMargin;
		NPVFlow *=	spread;
		
		if ((its_CreditPremiumLegAccrued == CAP_PRORATA) && (TheRatio != 1.0))
		{
			// Adjustment
			NPVFlow	*=	TheRatio;
			NPVFlowForATMMargin	*=	TheRatio;
			FlowForATMMargin	*=	TheRatio;
		}
	}

	Outputs.clear();
	Outputs.push_back(NPVFlow);
	Outputs.push_back(NPVFlowForATMMargin);
	Outputs.push_back(FlowForATMMargin);

}

// ----------------------------------------------------------------------------
// ComputeLossesAmountBeforeCreditWindowsDates
// ----------------------------------------------------------------------------

void CreditManager::ComputeLossesAmountBeforeCreditWindowsDates()
{
	// to be optimized
	// for each Credit Date, find the cumulative loss
	set<Tau_Item>::iterator iter;
	set<Tau_Item>::iterator iter_DefTime;
	set<Tau_Item>::iterator iter_prev_end;
	set<Tau_Item>::iterator iter_DefTimeBis;

//	set<Tau_Item>::iterator iter_DefTime_Lower;
//	set<Tau_Item>::iterator iter_DefTime_Upper;
//	pair<set<Tau_Item>::iterator, set<Tau_Item>::iterator>	pair_DefTime;

	double	CreditDate, tau_date, PrevCreditDate;
	double	TheCumulativeLoss	=	0.0;
	int		TheCumulativeNbDef	=	0;

	double	tmpDate;

	iter	= its_SortedCreditObservationDatesAndLosses.begin();
	while (iter != its_SortedCreditObservationDatesAndLosses.end())
	{
		// current Credit Date
		CreditDate	=	iter->tau;

		// RESET!!!
		iter->loss_at_that_time	=	0.0;
/*		
		iter_DefTime_Lower	=	pair_DefTime.first;
		iter_DefTime_Upper	=	pair_DefTime.second;

		// test
		if (iter_DefTime_Upper
		!=	itsSortedDefaultTimes.end())
		{
			tau_date	=	iter_DefTime_Upper->tau;
			
			if ((iter_DefTimeBis	=	itsSortedDefaultTimes.find(*iter)) != itsSortedDefaultTimes.end()) // tau_date == CreditDate
			{
				iter->loss_at_that_time	=	iter_DefTimeBis->loss_at_that_time;
			}
*/
		iter_DefTime	=	itsSortedDefaultTimes.upper_bound(*iter);
		
		// tes
		if (iter_DefTime !=	itsSortedDefaultTimes.end())
		{
			tau_date	=	iter_DefTime->tau;

			if ((iter_DefTimeBis	=	itsSortedDefaultTimes.find(*iter)) != itsSortedDefaultTimes.end()) // tau_date == CreditDate
			{
				tmpDate	=	iter_DefTimeBis->tau;
				iter->loss_at_that_time	=	iter_DefTimeBis->loss_at_that_time;

				// just ugly, in order to gain some code lines!
//				iter_DefTime	=	iter_DefTimeBis;
			}
			
			// Cumulative Loss + Current Loss
			TheCumulativeLoss	=	iter_DefTime->cumul_loss;
			TheCumulativeNbDef	=	iter_DefTime->cumul_nbdef;

			// Set the Cumulative Loss
			iter->cumul_loss	=	TheCumulativeLoss;
			iter->cumul_nbdef	=	TheCumulativeNbDef;	// just posterior

			++iter;

		}
		else
		{
			// I am already out of the range, so that for all posterior Credit Dates (the current included)
			// the Loss Amount, is given by the last computed CumulativeLoss
			if (itsSortedDefaultTimes.size() > 0)
			{
				iter_prev_end	=	itsSortedDefaultTimes.end();
				--iter_prev_end;
				PrevCreditDate		=	iter_prev_end->tau;
				TheCumulativeLoss	= iter_prev_end->cumul_loss + its_LossesAmount[iter_prev_end->id];
				TheCumulativeNbDef	= iter_prev_end->cumul_nbdef + 1;
			}

			while (iter != its_SortedCreditObservationDatesAndLosses.end())
			{
				iter->cumul_loss	=	TheCumulativeLoss;
				iter->cumul_nbdef	=	TheCumulativeNbDef;		// to be checked!
				++iter;
			}
		}
	}

	// ---------------------------------------------------
	if (its_Display_LossesAmount_Before_CreditDates_Flag)
	{
		//	set<Tau_Item>::iterator iter;
		if (its_temp_index <= 0)		// either Price or First hedge Credit
		{
			// ----------------------------------------------------------------------------
			if ((its_fOut = fopen("c:\\test\\mc_ComputeLossesAmountBeforeCreditWindowsDates.txt", "a+")) == NULL) return;
			fprintf(its_fOut, " ----------------- ComputeLossesAmountBeforeCreditWindowsDates:\t\t%u ----------------- \n", its_SimulId);
			// ----------------------------------------------------------------------------
			for (iter = its_SortedCreditObservationDatesAndLosses.begin();
						iter != its_SortedCreditObservationDatesAndLosses.end();
						++iter)
			{
				fprintf(its_fOut, "Id:\t%u\t\tCredit Date:\t%.2lf\t\tCumul loss:\t%.2lf\t\tCumul NbDef:\t%u\n", iter->id, iter->tau, iter->cumul_loss, iter->cumul_nbdef);
			}
			fclose(its_fOut);
		}
	}
}


void CreditManager::GetHowMuchLossBeforeACreditObservation(int LowOrUpFlag, double& TheLoss)
{
	// hypothesis: TheDate must be in 'its_SortedCreditObservationDatesAndLosses'
	// uses the current id: its_CurrentIndex
	// otherwise:	implementation of a search in its_Default

	set<Tau_Item>::iterator	iter;

	// LowOrUpFlag = 0 --> Low
	if (LowOrUpFlag == 0)
	{
		iter	=	its_CreditObervationDatesSortedLowsIds[its_CurrentIndex];
		TheLoss	=	iter->cumul_loss	-	iter->loss_at_that_time;
	}
	else
	{
		iter	=	its_CreditObervationDatesSortedUpsIds[its_CurrentIndex];
		TheLoss	=	iter->cumul_loss;
	}
}


void CreditManager::GetHowMuchLossBeforeADate(RelativeDate TheDate, double& TheLoss)
{
	if (itsSortedDefaultTimes.size() == 0)
	{
		TheLoss	=	0.0;
		return;
	}

	set<Tau_Item>::iterator	iter	=	itsSortedDefaultTimes.begin();
	set<Tau_Item>::iterator	iter_prev_end	=	itsSortedDefaultTimes.end();	
	--iter_prev_end;

	// the Item to be searched
	Tau_Item	SearchItem(0, (double) TheDate, 0.0);

	// iter points at the first item which is above or equal TheDate
	iter	=	itsSortedDefaultTimes.lower_bound(SearchItem);
	
	if (iter == itsSortedDefaultTimes.end())
	{		
		TheLoss	=	iter_prev_end->cumul_loss;
		TheLoss	+=	its_LossesAmount[iter_prev_end->id];
	}
	else
		TheLoss	=	iter->cumul_loss;

}


void CreditManager::GetHowManyDefaultsBeforeACreditObservation(int LowOrUpFlag, int& TheNbDef, set<Tau_Item>::iterator&	iter)
{
	// hypothesis: TheDate must be in 'its_SortedCreditObservationDatesAndLosses'
	// uses the current id: its_CurrentIndex
	// otherwise:	implementation of a search in its_Default

	// LowOrUpFlag = 0 --> Low
	if (LowOrUpFlag == 0)
		iter	=	its_CreditObervationDatesSortedLowsIds[its_CurrentIndex];
	else
		iter	=	its_CreditObervationDatesSortedUpsIds[its_CurrentIndex];

	TheNbDef	=	iter->cumul_nbdef;
}


void CreditManager::GetHowManyDefaultsBeforeADate(RelativeDate TheDate, int& TheNbDef, set<Tau_Item>::iterator&	iter)
{
	iter	=	itsSortedDefaultTimes.begin();

	if (itsSortedDefaultTimes.size() == 0)
	{
		TheNbDef	=	0;
		return;
	}

	set<Tau_Item>::iterator	iter_prev_end	=	itsSortedDefaultTimes.end();
	--iter_prev_end;

	// the Item to be searched
	Tau_Item	SearchItem(0, (double) TheDate, 0.0);

	// iter points at the first item which is above or equal TheDate
	iter	=	itsSortedDefaultTimes.lower_bound(SearchItem);
	
	if (iter == itsSortedDefaultTimes.end())
	{
		TheNbDef	=	iter_prev_end->cumul_nbdef;
		TheNbDef++;
	}
	else
		TheNbDef	=	iter->cumul_nbdef;

}


void CreditManager::GetHowManyDefaultsBeforeADate(RelativeDate TheDate, double& TheNbDef)
{
	set<Tau_Item>::iterator	iter;
	
	iter	=	itsSortedDefaultTimes.begin();

	if (itsSortedDefaultTimes.size() == 0)
	{
		TheNbDef	=	0;
		return;
	}

	set<Tau_Item>::iterator	iter_prev_end	=	itsSortedDefaultTimes.end();
	--iter_prev_end;

	// the Item to be searched
	Tau_Item	SearchItem(0, (double) TheDate, 0.0);

	// iter points at the first item which is above or equal TheDate
	iter	=	itsSortedDefaultTimes.lower_bound(SearchItem);
	
	if (iter == itsSortedDefaultTimes.end())
	{
		TheNbDef	=	iter_prev_end->cumul_nbdef;
		TheNbDef++;
	}
	else
		TheNbDef	=	iter->cumul_nbdef;

}



// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
//	MONTE-CARLO DELTA DERIVATIVES COMPUTATION
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void	CreditManager::HedgesFastSpreadParallel_Monte_Carlo(double& NPV)
{
	if (!its_HedgesRunning) return;

	int	i;
	DoubleVector	Outputs;
	
	ICM_DefaultCurve*	DefCurve	=	NULL;
	ICM_DefaultCurve*	ShiftedDefCurve	=	NULL;

	DoubleVector	TheLossesAmountShifted;

	double	PremLeg;
	double	DefLeg;

	double	delta, shift;
	DoubleVector	valueDefault;
	DoubleVector	valueNonDefault;
	double	defaultBarrierDerivative;
	double	modifiedDefaultTimesBeforeShift;
	double	modifiedDefaultTimesAfterShift;
	DoubleVector	valueModDefault;
	DoubleVector	valueModDefaultWithShift;
	double	defaultBarrier;

	double	delta_PremLeg, delta_DefLeg;

	if (its_Bump_Choice != CHB_FAST_SPREAD_PARALLEL)
		ICMTHROW(ERR_INVALID_DATA,"Only Spread Sensitivities is allowed for FAST MC HEDGES!");

	double	TimeHorizon	=	MATHTIME(its_MaxDate);

	int	scenario_id;
	scenario_id	=	0;

	// First passage
	if (its_CurrentHedgesIndex	==	-1)
	{
		fill(its_AllNPVS.begin(), its_AllNPVS.end(), 0.0);
		fill(its_AllDefLegs.begin(), its_AllDefLegs.end(), 0.0);
		fill(its_AllPremLegs.begin(), its_AllPremLegs.end(), 0.0);

		// -----------------------------------------------------
		// FAST HEDGE DATA ALLOCATIONS

		its_CommonFactors.resize(1);		// 1 Factor model
		its_BarrierDerivatives_FastHedge.resize(Get_NbCredits());

		// WHAT VALUE???
		shift	=	1e-3;
		// -----------------------------------------------------

		// -------------------------------------------------------
		// Compute the Directional Shifts
		ComputeDirectionalShift(scenario_id);
		// -------------------------------------------------------

		// -------------------------------------------------------
		// Reset Outputs
		// -------------------------------------------------------
		ResetOutputs();	// NPV to 0.0

		// -------------------------------------------------------
		// Schedule Checkings
		// -------------------------------------------------------
//		ScheduleCheckings();	// already done in Central NPV

		// -------------------------------------------------------
		// Compute Useful DF Data
		// -------------------------------------------------------
//		ComputeAllDF();			// already done in Central NPV

		// -------------------------------------------------------
		// Compute Useful Data for the Basket
		// -------------------------------------------------------
//		ComputeBasketNotional(); // already done in Central NPV

		// -------------------------------------------------------
		// Prepare for Pricing
		// -------------------------------------------------------
//		SortAllCreditObservationDates(); // already done in Central NPV

		// -------------------------------------------------------
		// All Credit Curves have already been generated
		// -------------------------------------------------------
	//	GenerateAllDefaultCurves();

		// -------------------------------------------------------
		// Initialization of the Random Generator
		// -------------------------------------------------------
		SetRandomGenerator();

		// -------------------------------------------------------
		// Compute Barriers for Default Times
		// -------------------------------------------------------
//		ComputeMaxCreditDate(); // already done in Central NPV
		
		ComputeBarriers(its_MaxCreditDate);
		ComputeBarriersShifted(its_MaxCreditDate, scenario_id);

		// -------------------------------------------------------
		// According to Correlation
		// -------------------------------------------------------
		switch (its_CorrelationType)
		{
			case CT_FLAT:
			case CT_BETA:				
			case CT_FACTOR_LOADING_2:
				break;

			case CT_MATRIX:
				// Cholesky computation
				CholeskyComputation();
				break;

		}

		// -------------------------------------------------------
		// The Pricing Loop  
		// -------------------------------------------------------
		for (its_SimulId=0;its_SimulId<its_NbSimul;its_SimulId++)
		{
			// -------------------------------------------------------
			// According to Correlation
			// -------------------------------------------------------
			switch (its_CorrelationType)
			{
				case CT_FLAT:
					// take care of correl = 1
				case CT_BETA:					
				case CT_FACTOR_LOADING_2:
					// FAST HEDGES = true
					GenerateDefaultTimes_Betas_ForHedge(scenario_id, true);
					break;

				case CT_MATRIX:
					// Cholevsky
					break;
			}

			// -------------------------------------------------------
			// ComputeCumulativeLossesAndNbDefaults
			// with initial Curves (non shifted)
			// -------------------------------------------------------
			if (IsCDOStandard())
				ComputeCumulativeLossesAndNbDefaults();
			else // if (IsCDOSquare())
				ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE();

			Keep_SortedDefaultTimes();

			// -------------------------------------------------------
			// Here is the difference for the Hedge
			// -------------------------------------------------------

//			Display_DefaultTimes();		// and shifted

			for (i=0;i<Get_NbCredits();i++)
			{
				its_temp_index	=	i;

				// -------------------------------------------------------
				// Get Value with Default of Current Credit at MATURITY
				// -------------------------------------------------------

				// -------------------------------------------------------
				// UpdateCumulativeLossesAndNbDefaults
				// with initial Curves and the shifted one
				// Price Basket For This Simul
				// -------------------------------------------------------
				SetShiftedDefaultTime(i, TimeHorizon);
				UpdateCumulativeLossesAndNbDefaults(i);
				PriceBasketForThisSimulation(valueDefault);				
				
				// TO BE IMPROVED: restore Sorted Default Times
				Restore_SortedDefaultTimes();
				// -------------------------------------------------------

				// -------------------------------------------------------
				// Get Value with NON Default of Current Credit at MATURITY
				// -------------------------------------------------------

				// -------------------------------------------------------
				SetShiftedDefaultTime(i, TimeHorizon + 1.0);
				UpdateCumulativeLossesAndNbDefaults(i);
				PriceBasketForThisSimulation(valueNonDefault);				

				// TO BE IMPROVED: restore Sorted Default Times
				Restore_SortedDefaultTimes();
				// -------------------------------------------------------

				// -------------------------------------------------------
				// Get Derivative of Default Barrier for the Current Credit
				// -------------------------------------------------------
				defaultBarrierDerivative	=	ComputeBarrierDerivative(i, scenario_id, TimeHorizon);

				// -------------------------------------------------------
				// Get Modified Default Times
				// -------------------------------------------------------
				modifiedDefaultTimesBeforeShift	=	ComputeModifiedDefaultTime(i, 0.0);
				modifiedDefaultTimesAfterShift	=	ComputeModifiedDefaultTime(i, shift);

				// -------------------------------------------------------
				// Get Values with modified Default Times with and without Shift
				// -------------------------------------------------------

				// -------------------------------------------------------
				SetShiftedDefaultTime(i, modifiedDefaultTimesBeforeShift);
				UpdateCumulativeLossesAndNbDefaults(i);
				PriceBasketForThisSimulation(valueModDefault);				

				// TO BE IMPROVED: restore Sorted Default Times
				Restore_SortedDefaultTimes();
				// -------------------------------------------------------

				// -------------------------------------------------------
				SetShiftedDefaultTime(i, modifiedDefaultTimesAfterShift);
				UpdateCumulativeLossesAndNbDefaults(i);
				PriceBasketForThisSimulation(valueModDefaultWithShift);				

				// TO BE IMPROVED: restore Sorted Default Times
				Restore_SortedDefaultTimes();
				// -------------------------------------------------------

				// Get Default Barrier
				defaultBarrier	=	0.0;

				// -------------------------------------------------------
				// Calculate delta
				delta	=	valueDefault[0] - valueNonDefault[0] * defaultBarrierDerivative;
				delta	+=	(valueModDefaultWithShift[0] - valueModDefault[0]) * defaultBarrier / shift;

				delta_DefLeg	=	valueDefault[1] - valueNonDefault[1] * defaultBarrierDerivative;
				delta_DefLeg	+=	(valueModDefaultWithShift[1] - valueModDefault[1]) * defaultBarrier / shift;

				delta_PremLeg	=	valueDefault[2] - valueNonDefault[2] * defaultBarrierDerivative;
				delta_PremLeg	+=	(valueModDefaultWithShift[2] - valueModDefault[2]) * defaultBarrier / shift;
				// -------------------------------------------------------

				its_AllNPVS[i]		+=	delta;
				its_AllDefLegs[i]	+=	delta_DefLeg;
				its_AllPremLegs[i]	+=	delta_PremLeg;

				// restore Sorted Default Times
				Restore_SortedDefaultTimes();

			}

		}

		NPV	=	0.0;

	}
	else
	{
		// simply get the right value
		NPV	=	its_AllNPVS[its_CurrentHedgesIndex];
		NPV	/=	its_NbSimul;

		PremLeg	=	its_AllPremLegs[its_CurrentHedgesIndex];
		PremLeg	/=	its_NbSimul;

		DefLeg	=	its_AllDefLegs[its_CurrentHedgesIndex];
		DefLeg	/=	its_NbSimul;
	}

}




// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
//	Compute Directional Shift in Default Intensities
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
void	CreditManager::ComputeDirectionalShift(int col_id)
{
	int	i, k;
	int	TheDim;

	ARM_Vector	TheLambdas;
	ARM_Vector	TheShiftedLambdas;

	ARM_Vector	DirectionalShift;

	// -----------------------------------------------------------
	// FAST HEDGES
	// -----------------------------------------------------------

	for (i=0; i<Get_NbCredits(); i++)
	{
		TheLambdas			=	*(its_ArrayDefaultCrv[i])->GetLambdas();
		TheShiftedLambdas	=	*(*its_MatrixShiftedDefaultCrv)(i, col_id)->GetLambdas();

		TheDim	=	TheLambdas.GetSize();
		
		DirectionalShift.Resize(TheDim);

		for (k=0; k<TheDim; k++)
			DirectionalShift.Elt(k)	=	TheShiftedLambdas.Elt(k) - TheLambdas.Elt(k);

		(its_ArrayDefaultCrv[i])->SetDirectionalLambdasShift((ARM_Vector*)DirectionalShift.Clone());
	}
}
 

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
//	Compute Barrier Derivative for a given Credit Num at time T
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
double	CreditManager::ComputeBarrierDerivative(int Num, int col_id, double T)
{
	double	TheBarrier;

	TheBarrier	=	its_BarrierDerivatives_FastHedge[Num];
	return	TheBarrier * (1.0 - exp((its_ArrayDefaultCrv[Num])->SurvivalProba(its_MaxDate) / ((*its_MatrixShiftedDefaultCrv)(Num, col_id))->SurvivalProba(its_MaxDate)));

}


double	CreditManager::ComputeModifiedDefaultTime(int Num, double Shift)
{
	double	uniform;
	double	aux;

	double	shifted_Zi;		// Zi = ri
	double	bari;
	
	double	shifted_tau;

	bari	=	itsBarriers_Standard[Num];

	DoubleVector	TheBetas;
	DoubleVector	TheComplementaryBetas;
	DoubleVector	TheCommonFactors;
	
	double			TheCommonFactor;	// 1F !!!
	double			TheComplementaryBeta;
	double			TheBeta;

	TheBetas				=	GetBetas();
	TheComplementaryBetas	=	GetComplementaryBetas();	
	TheCommonFactors		=	GetCommonFactors();

	TheBeta					=	TheBetas[Num];
	TheComplementaryBeta	=	TheComplementaryBetas[Num];
	TheCommonFactor			=	TheCommonFactors[0];		// 1 Factor

	double	pp;
	double	ConditionalCumulativeDefaultDistribution;

	pp	=	(bari - TheCommonFactor * TheBeta);

	if (CHECK_EQUAL(TheComplementaryBeta, 0.0))
		if (pp > 0.0)
			ConditionalCumulativeDefaultDistribution	=	1.0;
		else
			ConditionalCumulativeDefaultDistribution	=	0.0;
	else
		ConditionalCumulativeDefaultDistribution	=	 NAG_cumul_normal(pp / TheComplementaryBeta);

	// uniform
	uniform	=	NAG_random_continuous_uniform_ab(0.0, ConditionalCumulativeDefaultDistribution); 

	// to Normal
	aux	=	NAG_deviates_normal( uniform);

	shifted_Zi	=	TheCommonFactor * TheBeta + TheComplementaryBeta * aux;

//	switch (its_CopulaType)
//	case CCT_GAUSSIAN:

	aux = NAG_cumul_normal(shifted_Zi);

	shifted_tau	=	(its_ArrayDefaultCrv[Num])->DefProbInverse(aux);

	// I could update the Shifted default Times vector right now
//	itsSortedDefaultTimes.insert(Tau_Item(Num,t));

	return	shifted_tau;
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
//	IMPOSE DEFAULT TIME
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

void	CreditManager::SetShiftedDefaultTime(int Num, double DefaultTimeT)
{
	its_DefaultTimesShifted[Num]	=	DefaultTimeT;
}