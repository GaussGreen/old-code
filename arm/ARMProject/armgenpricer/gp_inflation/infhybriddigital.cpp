/*!
 *
 * Copyright (c) CDC IXIS CIB February 2006 Paris
 *	\file infcapfloorrielyield.cpp
 *
 *  \brief inflation capFloor on riel yield object
 *	\author  Y. Khlif
 *	\version 1.0
 *	\date February 2006
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpinflation/infhybriddigital.h"

/// gpinflation
#include "gpinflation/infcurv.h"
#include "gpinflation/infbsmodel.h"
#include "gpinflation/infmultibsmodel.h"
#include "gpinflation/infleg.h"
#include "gpinflation/infidx.h"
#include "gpinflation/assetinfo.h"
#include "gpinflation/infcapfloorrielyield.h"


/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/stringconvert.h"

///gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/spreadoption_lognormal_interface.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"
#include "gpclosedforms/gaussian_integrals.h"

/// kernel
#include <inst/swaption.h>
#include <inst/cmsleg.h>
#include <inst/fixleg.h>
#include <glob/dates.h>
#include <pricer/ipricer.h>


CC_USING_NS(std,string)


CC_BEGIN_NAMESPACE( ARM )

const double RESET_LAG_INF_IR = 7;
const double STEP_FOR_INTEGRAL = 120;
const double CALL_SPREAD_SHIFT_IR  =	0.000001;
const double CALL_SPREAD_SHIFT_INF =	0.00005;

///////////////////////////////////////////////////
///	Class  : ARM_InfHybridDigital
///	Routine: Constructor
///	Returns: 
///	Action : Constructor from a swap object
////////////////////////////////////////////////////
ARM_InfHybridDigital::ARM_InfHybridDigital(ARM_SwapLeg* payLeg, 
										   ARM_SwapLeg* digitLeg, 
										   int payOffType, 
										   ARM_ReferenceValue *barrierProfile,
										   int CF, 
										   int RecOrPay)
:
	/// cannot used the swaption constructor with swap as it
	/// imposes one fixed leg at least!
	ARM_Swaption(),
	itsPayLeg( NULL ),
	itsDigitLeg( NULL ),
	itsPayoffType(payOffType),
	itsCapOrFloor(CF), 
	itsPricingStrikesForView(NULL),
	itsInflationVolsForView(NULL),
	itsDigitLegVolsForView(NULL),
	itsCorrelationsForView(NULL)
{
	if( (!(dynamic_cast<ARM_InfLeg*>(payLeg))) && (!(dynamic_cast<ARM_InfLeg*>(digitLeg))))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": For GP_INF_Digital, one of the two legs should be an inflation leg!");

	ARM_SwapLeg* payLegClone = (ARM_SwapLeg*)(payLeg->Clone());
	ARM_SwapLeg* digitLegClone = (ARM_SwapLeg*)(digitLeg->Clone());

	// we use the attribute of ARM_SwapLeg to differ the payLeg from the Digitleg in the swap 
	//and not use for the PV pricing
	//so the payLeg is usually setted to K_PAY 
	// and the  digitLeg is usually setted to K_RCV
	payLegClone->SetInitRcvOrPay(K_PAY); 
	digitLegClone->SetInitRcvOrPay(K_RCV); 

	ARM_Swap* swap = new ARM_Swap(payLegClone, digitLegClone);
    ARM_Swap::Copy( (ARM_Object*) swap);
	SetName( ARM_INFHYBRIDDIGITAL );
	SetOptionType(RecOrPay);
	double barrier           = -1000000.0;
	SetStrikes(barrierProfile);


	/// fix and inflation
	if( (GetFirstLeg()->GetRcvOrPay())==K_PAY )
	{
		itsPayLeg	= GetFirstLeg();
		itsDigitLeg	= GetSecondLeg();		
	}
	else
	{
		itsPayLeg	= GetSecondLeg() ;
		itsDigitLeg	= GetFirstLeg();
	}

	if(itsPayLeg->GetPaymentFreq()!= itsDigitLeg->GetPaymentFreq())
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": the two legs should have the same pay Freqency and the same reset freqency!");

	if( strcmp( itsDigitLeg->GetCurrencyUnit()->GetCcyName(), itsPayLeg->GetCurrencyUnit()->GetCcyName() ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"inflation Digital only allowed for legs with same interest rate currency!");

	delete payLegClone;
	delete digitLegClone;
	delete swap;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfHybridDigital
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_InfHybridDigital::ARM_InfHybridDigital( const ARM_InfHybridDigital& rhs )
:	ARM_Swaption( rhs ),
	itsPricingStrikesForView(rhs.itsPricingStrikesForView),
	itsInflationVolsForView(rhs.itsInflationVolsForView),
	itsDigitLegVolsForView(rhs.itsDigitLegVolsForView),
	itsCorrelationsForView(rhs.itsCorrelationsForView)
{
	itsPayLeg = rhs.itsPayLeg ?  (ARM_SwapLeg*) rhs.itsPayLeg->Clone() : NULL;
	itsDigitLeg = rhs.itsDigitLeg ?  (ARM_SwapLeg*) rhs.itsDigitLeg->Clone() : NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfHybridDigital
///	Routine: Assignment operator
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_InfHybridDigital& ARM_InfHybridDigital::operator=( const ARM_InfHybridDigital& rhs )
{
	if( this != &rhs )
	{
		ARM_Swaption::operator=( rhs );
		itsPayLeg	= rhs.itsPayLeg;
		itsDigitLeg	= rhs.itsDigitLeg;
		itsPricingStrikesForView	= rhs.itsPricingStrikesForView;
		itsInflationVolsForView		= rhs.itsInflationVolsForView;
		itsDigitLegVolsForView		= rhs.itsDigitLegVolsForView;
		itsCorrelationsForView		= rhs.itsCorrelationsForView;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfHybridDigital
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_InfHybridDigital::~ARM_InfHybridDigital()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_InfHybridDigital
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_InfHybridDigital::toString() const
{
	return string();
}


////////////////////////////////////////////////////
///	Class  : ARM_InfHybridDigital
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_InfHybridDigital::Clone()
{
	return new ARM_InfHybridDigital(*this);
}
///////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: ComputePrice
///	Returns: ARM_Object*
///	Action : Computes the inflation swaption price
////////////////////////////////////////////////////
double ARM_InfHybridDigital::ComputePrice(int)
{
	ARM_Model* model		= GetModel();
	ARM_Date modelAsOfDate	= model->GetStartDate();
	double price			= 0.0;
	if( model->CanPriceInflation() >= PRICE_FWDNOPTION )
	{
		//////////////////////////////Model Information
		ARM_InfBSModel* pricingMod = dynamic_cast<ARM_InfBSModel*>(model);
		if( !pricingMod )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"model is of wrong type.. is supposed to price option on CPI but failed to do so because it does  not inherit from InfOptionModel.. Please advise");
	
		int dayCount			= pricingMod->GetInfFwdCurv()->GetMonthlyInterpType();
		ARM_Date lastKnownDate	= pricingMod->GetVolatility()->GetLastKnownDate();
		ARM_Date modelAsOfDate	= pricingMod->GetStartDate();

		ARM_InfBSModelPtr firstPricingMod = ARM_InfBSModelPtr((ARM_InfBSModel*) pricingMod->Clone());
		ARM_InfMultiBSModel* multiPricingMod = dynamic_cast<ARM_InfMultiBSModel*>(pricingMod);
		///FIN Model Information
		
		/////////////////////////////option Information
		int callPut	= GetCapOrFloor();
		int lonOrShort = GetOptionType();
		ARM_ReferenceValue* barrierProfile = GetStrikes();
		ARM_ReferenceValue notional(*itsPayLeg->GetAmount());
		int payoffType= GetPayoffType();
		///FIN option Information	

		if( dynamic_cast<ARM_InfLeg*>(itsPayLeg) )
		{
			ARM_InfLeg* infPayLeg = dynamic_cast< ARM_InfLeg* >(itsPayLeg->Clone());				
			int swapType		= infPayLeg->GetSwapType();
			if( swapType != K_YEARTOYEAR_LEG)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + ": only year on year leg are supported to price INF cap!");

			ARM_InfIdx* infPayidx	= (ARM_InfIdx*) infPayLeg->GetIRIndex();
			if( !infPayidx )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						ARM_USERNAME + ": null inflation index!");
			
			string indexPayName	= infPayidx->GetIndexName();

			if( multiPricingMod )
				firstPricingMod = multiPricingMod->GetCorrespondingInflationModel( infPayidx );
		}

		switch( payoffType )
		{			/// 1) inflation vs fix
		case K_SIMPLE:
			{
				//First Case pay Fix
				if( GetFixedLeg() )
				{				
					ARM_SwapLeg* fixLeg		= (ARM_SwapLeg*) (itsPayLeg->Clone());
				
					double optLet, CPIForward, timeToStart, tenor, discountfactor;
					int nbCapLet = fixLeg->GetPaymentDates()->GetSize();


					//////////////DigitalLeg Information
					if(!(dynamic_cast<ARM_InfLeg*>(itsDigitLeg)) )
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							ARM_USERNAME + ": if payLeg is fixed the digit leg have to be an inflation Leg!");

					ARM_InfLeg* infdigitLeg = dynamic_cast< ARM_InfLeg* >(itsDigitLeg->Clone());				
					int swapType		= infdigitLeg->GetSwapType();
					if( swapType != K_YEARTOYEAR_LEG)
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							ARM_USERNAME + ": only year on year leg are supported to price INF cap!");

					ARM_InfIdx* infDigitIdx	= (ARM_InfIdx*) infdigitLeg->GetIRIndex();
					if( !infDigitIdx )
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
								ARM_USERNAME + ": null inflation index!");
					
					string indexDigitName	= infDigitIdx->GetIndexName();

					if( multiPricingMod )
						firstPricingMod = multiPricingMod->GetCorrespondingInflationModel( infDigitIdx );
					///FIN DigitalLeg Information


					int infRoP	= infdigitLeg->GetRcvOrPay();
					int fixRoP	= fixLeg->GetRcvOrPay();
										
					double leverage = infdigitLeg->GetMultiple();
					double constant = infdigitLeg->GetConstant();
					double spread = leverage+constant;

					//For View
					itsPricingStrikesForView = ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
					itsInflationVolsForView = ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));									

					ARM_Date numDate, denomDate;
					int k;				
					for (k = 0; k<nbCapLet; k++) 
					{
						numDate			= ARM_Date( infdigitLeg->GetNumResetDates()->Elt(k) );
						denomDate		= ARM_Date( infdigitLeg->GetDenomResetDates()->Elt(k) );

						double interestTerm	= infdigitLeg->GetInterestTerms()->Elt(k);
						
						///////////the infdigitLegInfo by flow
						
						///GetFwdRates() return the (leverage*forward + spread)*interestterm*infRoP
						///where forward = CPI(i)/CPI(i-1) -1
						///here we compute the CPIForward=leverage*(forward+1)
						CPIForward		= infdigitLeg->GetFwdRates()->Elt(k) / CC_NS( ARM_Constants, rateBase )/(interestTerm*infRoP) - constant;  
						discountfactor	= infdigitLeg->GetDiscountFactors()->Elt(k);
						double payTime	= infdigitLeg->GetPaymentDates()->Elt(k);
					
						///////////the fixLegInfo by flow
						double fixlegRate = fixLeg->GetFwdRates()->Elt(k)/ CC_NS( ARM_Constants, rateBase );
						double fixInterestTerm	= fixLeg->GetInterestTerms()->Elt(k);
						
						
						////compute the strike for the caplet
						double barrier = barrierProfile->CptReferenceValue(payTime);
																				
						/// to account for fixing difference!
						double	renormalisationFactor = 1.0;
						timeToStart		= firstPricingMod->GetModelTimeWPublishLag( denomDate, infDigitIdx );
						if( timeToStart < 0.0 )
						{
							/// in this particular case we need to account for a fixing different from the one
							/// of the zero coupon option
							/// the renormalisation is YtYFixingCPI/CurrentCPI
							/// or the denomCPIRates / Reference CPI of the Curve
							renormalisationFactor = infdigitLeg->GetDenomCPIRates()->Elt(k) / firstPricingMod->GetCPIIndexValue();
						}

						double timeToExpiry		= firstPricingMod->GetModelTimeWPublishLag( numDate, infDigitIdx);
						tenor = timeToExpiry	- timeToStart;
						double remainingTenor	= tenor;
						if( timeToStart < 0.0 )
							remainingTenor = remainingTenor+timeToStart <0? 0 : remainingTenor+timeToStart;

									
						////compute the volatilit of the CPIFwd
						double TjVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, denomDate ); 
						ARM_Date denomDatewPublishLag = firstPricingMod->GetModelDateWPublishLag( denomDate, infDigitIdx );
						double TjForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, denomDatewPublishLag );

						double TiVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, numDate ); 
						double TiFromAsOf		= CountYearsWithoutException( dayCount, modelAsOfDate, numDate ); 
						ARM_Date numDatewPublishLag = firstPricingMod->GetModelDateWPublishLag( numDate, infDigitIdx );
						double TiForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, numDatewPublishLag );
						tenor					= TiVolLookup-TjVolLookup;


						double pricingStrike = barrier/CC_NS( ARM_Constants, rateBase ) - constant;
						double amount = notional.CptReferenceValue(payTime);

						if( TiForMaturity < 0.0 || tenor < 0.0 || ((tenor+TjForMaturity) < 0.0))
						{
							tenor = -1;
							optLet = (callPut*(CPIForward-pricingStrike) > 0.0) ? 1.0 : 0.0;
											
							price+= discountfactor*lonOrShort*optLet*fixInterestTerm*amount*fixlegRate;
							
						}
						else
						{						
							if(TjForMaturity <= 0.0)
							{
								tenor		   += TjForMaturity;
								TjForMaturity	= StringMaturityToYearTerm( "1d" );
								TjVolLookup		= StringMaturityToYearTerm( "1d" );
							}
							if( TjVolLookup <= 0.0 )							
								TjVolLookup		= StringMaturityToYearTerm( "1d" );
							
							if( TiForMaturity < K_NEW_DOUBLE_TOL)
								TiForMaturity = StringMaturityToYearTerm( "1d" );							

							double totalvolatility_plus		= K_NEW_DOUBLE_TOL;
							double totalvolatility_moins	= K_NEW_DOUBLE_TOL;

							if(leverage>1e-15)
							{
								double lookUpStrike_plus = pow( renormalisationFactor * (pricingStrike+CALL_SPREAD_SHIFT_INF)/leverage, 1.0/tenor) -1.0;
								double vol			= firstPricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike_plus , tenor )/CC_NS( ARM_Constants, volBase );
								totalvolatility_plus	= vol * sqrt( tenor );
								double lookUpStrike_moins = pow( renormalisationFactor * (pricingStrike-CALL_SPREAD_SHIFT_INF)/leverage, 1.0/tenor) -1.0;
								vol			= firstPricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike_moins , tenor )/CC_NS( ARM_Constants, volBase );
								totalvolatility_moins	= vol * sqrt( tenor );
							}				
							////////////////////////////End of compute the volatilit of the CPIFwd

							//BS Pricing
							if( fabs(pricingStrike) < K_NEW_DOUBLE_TOL)
								pricingStrike = K_NEW_DOUBLE_TOL;
							optLet = DigitalBlackSholesSmooth_Formula( CPIForward, totalvolatility_plus, totalvolatility_moins, discountfactor, pricingStrike, callPut, CALL_SPREAD_SHIFT_INF);
							price		   += lonOrShort*optLet*fixInterestTerm*amount*fixlegRate;

							//for View
							(*itsPricingStrikesForView)[k] = pricingStrike-leverage;
							(*itsInflationVolsForView)[k] = totalvolatility_plus;
						}
					}
					delete fixLeg;
					delete infdigitLeg;

				} // End Case Pay Fix

				else if( dynamic_cast<ARM_InfLeg*>(itsPayLeg) ) //Second Case Pay INF
				{
					if( GetFixedLeg())
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							"digital leg can not be fixed ");
					
				
					//////////////PayLeg Information
					ARM_InfLeg* infLeg = dynamic_cast< ARM_InfLeg* >(itsPayLeg->Clone());				
					int swapType		= infLeg->GetSwapType();
					if( swapType != K_YEARTOYEAR_LEG)
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							ARM_USERNAME + ": only year on year leg are supported to price INF cap!");

					ARM_InfIdx* infIdx	= (ARM_InfIdx*) infLeg->GetIRIndex();
					if( !infIdx )
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
								ARM_USERNAME + ": null inflation index!");
					
					string indexName	= infIdx->GetIndexName();

					if( multiPricingMod )
						firstPricingMod = multiPricingMod->GetCorrespondingInflationModel( infIdx );
					///FIN PayLeg Information

					/////////////Data For infLeg
					int infRoP	= infLeg->GetRcvOrPay();				
					double leverage = infLeg->GetMultiple();
					double constant = infLeg->GetConstant();
					double spread = leverage+constant;

					int nbCapLet = infLeg->GetPaymentDates()->GetSize();
					double CPIForward, timeToStart,tenor, discountfactor;
					ARM_Date numDate, denomDate;
					int k;
					
					//For View
					itsPricingStrikesForView	= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
					itsInflationVolsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
					itsDigitLegVolsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
					itsCorrelationsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));	
					
					if((dynamic_cast<ARM_InfLeg*>(itsDigitLeg)) ) //INFPay and INFDigit
					{
						ARM_InfLeg* infLeg_digit = dynamic_cast< ARM_InfLeg* >(itsDigitLeg->Clone());				
						int swapType_digit		= infLeg_digit->GetSwapType();
						if( swapType_digit != K_YEARTOYEAR_LEG)
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
								ARM_USERNAME + ": only year on year leg are supported to price INF cap!");

						ARM_InfIdx* infIdx_digit	= (ARM_InfIdx*) infLeg_digit->GetIRIndex();
						if( !infIdx_digit )
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
									ARM_USERNAME + ": null inflation index!");
						
						string indexName_digit	= infIdx_digit->GetIndexName();

						if(indexName_digit!= indexName)
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							 ARM_USERNAME + ": the case of INFPayLeg And INFDigitLeg must have the same inflation index!");

									
						double leverage_digit = infLeg_digit->GetMultiple();
						double constant_digit = infLeg_digit->GetConstant();
						double spread_digit = leverage_digit+constant_digit;
						if(fabs(leverage_digit)< K_NEW_DOUBLE_TOL)
						{
							for (k = 0; k<nbCapLet; k++) 
							{
								numDate			= ARM_Date( infLeg->GetNumResetDates()->Elt(k) );
								denomDate		= ARM_Date( infLeg->GetDenomResetDates()->Elt(k) );
								double interestTerm	= infLeg->GetInterestTerms()->Elt(k);
								
								///////////the infdigitLegInfo by flow
								
								///GetFwdRates() return the (leverage*forward + spread)*interestterm 
								///where forward = CPI(i)/CPI(i-1) -1
								///here we compute the CPIForward=leverage*(forward+1)
								double payFlow		= infRoP*infLeg->GetFwdRates()->Elt(k) / CC_NS( ARM_Constants, rateBase );  
								discountfactor	= infLeg->GetDiscountFactors()->Elt(k);
								double payTime	= infLeg->GetPaymentDates()->Elt(k);
								double amount = notional.CptReferenceValue(payTime);
								double barrier = barrierProfile->CptReferenceValue(payTime)/CC_NS( ARM_Constants, rateBase );

								double digit = (callPut*(spread_digit-barrier)>0 ) ? 1.0 : 0.0;
								double oplet = lonOrShort*discountfactor*amount*payFlow*digit;
								price += oplet;											
							}
						}
						else
						{
							
							ARM_ReferenceValue* capStrike = (ARM_ReferenceValue*) barrierProfile->Clone();
							int sizeRefValue = barrierProfile->GetDiscreteValues()->GetSize();
							ARM_ReferenceValue* digitalCoupon = (ARM_ReferenceValue*) barrierProfile->Clone();
							for (k = 0; k<sizeRefValue; k++) 
							{
								double barrier = barrierProfile->GetDiscreteValues()->Elt(k)/CC_NS( ARM_Constants, rateBase );
								double capStrikeValue = infRoP*(leverage*(barrier - spread_digit)/leverage_digit+spread)*CC_NS( ARM_Constants, rateBase );
								capStrike->GetDiscreteValues()->Elt(k) = capStrikeValue;
								digitalCoupon->GetDiscreteValues()->Elt(k) = infRoP*capStrikeValue;
							}
							int capCallPut_RespectLeverage = (leverage > 0.0 ) ? callPut : - callPut;
							int capCallPut=capCallPut_RespectLeverage*infRoP;
							ARM_Date startDate =infLeg->GetStartDateNA();
							ARM_Date endDate=infLeg->GetEndDateNA();
							
							double fixedRate = 0.0; 
							int freqPay = infLeg->GetPaymentFreq();
							ARM_FixLeg fixLegForCap(startDate,
													   endDate, 
													   fixedRate, 
													   K_RCV, 
													   freqPay);
							ARM_Swap swapForCap(infLeg,&fixLegForCap); 
						
							ARM_InfCapFloorRielYield infCap(&swapForCap, capCallPut, capStrike);

							int infDayCount   =infLeg->GetDayCount();
							int decompFreq = infLeg->GetDecompFreq();
						    int payTiming  = infIdx->GetPayTiming();
							int intRule    = infIdx->GetIntRule();
														   
							ARM_FixLeg fixedLegForDigital(startDate,
														   endDate,
														   digitalCoupon,
														   K_RCV, 
														   freqPay,
														   infDayCount,
														   decompFreq,
														   payTiming ,
														   intRule);
							
							ARM_InfHybridDigital digital(&fixedLegForDigital, 
										   infLeg_digit, 
										   K_SIMPLE, 
										   barrierProfile,
										   callPut, 
										   lonOrShort);

							
							ARM_IFPricer pricerCap(&infCap, &(*firstPricingMod));
							double capPrice = pricerCap.Price();

							ARM_IFPricer pricerDigit(&digital, &(*firstPricingMod));
							double digitPrice = pricerDigit.Price();

							price = lonOrShort*capCallPut_RespectLeverage*capPrice+digitPrice;

							delete digitalCoupon;
							delete capStrike;
						}
						
						
					}//End INFPay and INFDigit
				
					else //INFPay and  is LiborDigit
					{
						ARM_SwapLeg* irLeg	= (ARM_SwapLeg*) (itsDigitLeg->Clone());
						ARM_IRIndex* irIdx	= irLeg->GetIRIndex();
						int irRorP = irLeg->GetRcvOrPay();
						ARM_ReferenceValue* spreadRefValue = irLeg->GetSpreads();
						if(spreadRefValue == NULL)
							spreadRefValue = new ARM_ReferenceValue(irLeg->GetSpread());

						for (k = 0; k<nbCapLet; k++) 
						{
							numDate			= ARM_Date( infLeg->GetNumResetDates()->Elt(k) );
							denomDate		= ARM_Date( infLeg->GetDenomResetDates()->Elt(k) );
							double interestTerm	= infLeg->GetInterestTerms()->Elt(k);
							
							///////////the infdigitLegInfo by flow
							
							///GetFwdRates() return the (leverage*forward + spread)*interestterm 
							///where forward = CPI(i)/CPI(i-1) -1
							///here we compute the CPIForward=leverage*(forward+1)
							CPIForward		= infLeg->GetFwdRates()->Elt(k) / CC_NS( ARM_Constants, rateBase )/(interestTerm*infRoP) - constant;  
							discountfactor	= infLeg->GetDiscountFactors()->Elt(k);
							double payTime	= infLeg->GetPaymentDates()->Elt(k);				
															
							/// to account for fixing difference!
							double	renormalisationFactor = 1.0;
							timeToStart		= pricingMod->GetModelTimeWPublishLag( denomDate, infIdx );
							
							double timeToExpiry		= pricingMod->GetModelTimeWPublishLag( numDate, infIdx);
							tenor = timeToExpiry	- timeToStart;
													

							////////////////////////////////////compute the volatilit of the CPIFwd
							double TjVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, denomDate ); 
							ARM_Date denomDatewPublishLag = pricingMod->GetModelDateWPublishLag( denomDate, infIdx );
							double TjForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, denomDatewPublishLag );

							double TiVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, numDate ); 
							double TiFromAsOf		= CountYearsWithoutException( dayCount, modelAsOfDate, numDate ); 
							ARM_Date numDatewPublishLag = pricingMod->GetModelDateWPublishLag( numDate, infIdx );
							double TiForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, numDatewPublishLag );
							tenor					= TiVolLookup-TjVolLookup;

							double irResetDate	= irLeg->GetResetDates()->Elt(k);
							double irResetTerm = (irResetDate-modelAsOfDate.GetJulian())/K_YEAR_LEN;

							if( TiForMaturity < 0.0 || tenor < 0.0 || ((tenor+TjForMaturity) < 0.0) || (irResetTerm<0.0))
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
										ARM_USERNAME + ": reset dates are before the AsOfDate can't priced the LiborLeg we don't have a reset manager for Libor Rates !" );


							else
							{
								if( TjForMaturity <= 0.0 )
								{
									tenor		   += TjForMaturity;
									TjForMaturity	= StringMaturityToYearTerm( "1d" );
									TjVolLookup		= StringMaturityToYearTerm( "1d" );
								}
								if( TjVolLookup <= 0.0 )
								{							
									TjVolLookup		= StringMaturityToYearTerm( "1d" );
								}

								double infResetTerm = (numDatewPublishLag.GetJulian()-modelAsOfDate.GetJulian())/K_YEAR_LEN;

								if( TiForMaturity < K_NEW_DOUBLE_TOL)
								{
										TiForMaturity = StringMaturityToYearTerm( "1d" );
										infResetTerm = 1.0/365.0;
								}

									///////////the irLegInfo by flow
								double irlegRate = irLeg->GetFwdRates()->Elt(k)/CC_NS( ARM_Constants, rateBase );
								
								double irInterestTerm	= irLeg->GetInterestTerms()->Elt(k);
								
												
								////compute the strike for the caplet
								double barrier = barrierProfile->CptReferenceValue(payTime);
								double spread_IR = spreadRefValue->CptReferenceValue(irResetDate)/CC_NS( ARM_Constants, rateBase );
								
								double pricingStrike = barrier/CC_NS( ARM_Constants, rateBase ) - spread_IR;

								double resetLag=fabs(irResetTerm-infResetTerm)*K_YEAR_LEN;
								
								double totalVol= K_NEW_DOUBLE_TOL;										
								if(leverage>1e-15)
								{
									double lookUpStrike = 0.0; //Liquid strike for market vol
									double vol			= pricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike , tenor )/CC_NS( ARM_Constants, rateBase );
									totalVol	= vol * sqrt( tenor );
								}
								////////////////////////////End of compute the volatilit of the CPIFwd
								double irVol = 0.0;
								double correl = 0.0;
								
								if(resetLag>RESET_LAG_INF_IR)
									throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
										ARM_USERNAME + ": resetlag between Libor and INF not allowed to exceed 7 days !" );
								
								else
								{
									double irTenor = irIdx->GetYearTerm();

									if ( irLeg->GetName() == ARM_CMSLEG )
										irTenor = ((ARM_CMSLeg *) irLeg)->GetSwapYearTerm();
									double irmoyeness = (pricingStrike-irlegRate)*CC_NS( ARM_Constants, rateBase ); //ATM Vol for the pay underlying
									irVol = pricingMod->GetIRModel()->GetVolatility()->ComputeVolatility( irResetTerm, irmoyeness,irTenor )/CC_NS( ARM_Constants, rateBase );


									double irmoyeness_plus = ((pricingStrike+CALL_SPREAD_SHIFT_IR)-irlegRate)*CC_NS( ARM_Constants, rateBase ); //ATM Vol for the pay underlying
									double irVol_plus = pricingMod->GetIRModel()->GetVolatility()->ComputeVolatility( irResetTerm, irmoyeness_plus,irTenor )/CC_NS( ARM_Constants, rateBase );
									
									double irmoyeness_moins = ((pricingStrike-CALL_SPREAD_SHIFT_IR)-irlegRate)*CC_NS( ARM_Constants, rateBase ); //ATM Vol for the pay underlying
									double irVol_moins = pricingMod->GetIRModel()->GetVolatility()->ComputeVolatility( irResetTerm, irmoyeness_moins,irTenor )/CC_NS( ARM_Constants, rateBase );

									correl = pricingMod->GetInfIRCorrel( TiForMaturity,irTenor, infIdx, irIdx, "INF_YoY/IR_FWD" );

									if( fabs(pricingStrike) < K_NEW_DOUBLE_TOL)
										pricingStrike = K_NEW_DOUBLE_TOL;													
									
									double totalirVol = irVol*sqrt(irResetTerm);
									double totalirVol_plus	=  irVol_plus*sqrt(irResetTerm);
									double totalirVol_moins = irVol_moins*sqrt(irResetTerm);
									double optLet;								
									double optLetDigitPayConstant;
									double optLetDigitPayFloat;
									
									optLetDigitPayConstant = constant*DigitalBlackSholesSmooth_Formula( irlegRate, totalirVol_plus, totalirVol_moins, discountfactor, pricingStrike, callPut, CALL_SPREAD_SHIFT_IR);
									
									double irlegRate_ADJ = irlegRate*exp(correl*totalirVol*totalVol);
									optLetDigitPayFloat = CPIForward*DigitalBlackSholesSmooth_Formula( irlegRate_ADJ, totalirVol_plus, totalirVol_moins, discountfactor, pricingStrike, callPut, CALL_SPREAD_SHIFT_IR);
									optLet = optLetDigitPayConstant+optLetDigitPayFloat;
									double amount = notional.CptReferenceValue(payTime);				
									price		   += lonOrShort*optLet*interestTerm*amount;
								}

								//for View
								(*itsPricingStrikesForView)[k] = pricingStrike;
								(*itsInflationVolsForView)[k] = totalVol;
								(*itsDigitLegVolsForView)[k] = irVol;
								(*itsCorrelationsForView)[k] = correl;
							}
						}						
						delete irLeg;			
					}
					delete infLeg;
				}  //End Second Case Pay INF	

				else //Third case Pay IR index
				{
					ARM_SwapLeg* irLeg		= (ARM_SwapLeg*) (itsPayLeg->Clone());
					ARM_IRIndex* irIdx	= irLeg->GetIRIndex();
					int irRorP = irLeg->GetRcvOrPay();
					ARM_ReferenceValue* spreadRefValue = irLeg->GetSpreads();
					if(spreadRefValue == NULL)
					{
						spreadRefValue = new ARM_ReferenceValue(irLeg->GetSpread());
					}

					int nbCapLet = irLeg->GetPaymentDates()->GetSize();
					double CPIForward, timeToStart, tenor, discountfactor;
				
					//////////////DigitalLeg Information
					if(!(dynamic_cast<ARM_InfLeg*>(itsDigitLeg)) )
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							ARM_USERNAME + ": if payLeg is Libor the digit leg have to be an inflation Leg!");

					ARM_InfLeg* infdigitLeg = dynamic_cast< ARM_InfLeg* >(itsDigitLeg->Clone());				
					int swapType		= infdigitLeg->GetSwapType();
					if( swapType != K_YEARTOYEAR_LEG)
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							ARM_USERNAME + ": only year on year leg are supported to price INF cap!");

					ARM_InfIdx* infDigitIdx	= (ARM_InfIdx*) infdigitLeg->GetIRIndex();
					if( !infDigitIdx )
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
								ARM_USERNAME + ": null inflation index!");
					
					string indexDigitName	= infDigitIdx->GetIndexName();

					if( multiPricingMod )
						firstPricingMod = multiPricingMod->GetCorrespondingInflationModel( infDigitIdx );
					///FIN DigitalLeg Information

					/////////////Data For infdigitLeg
					int infRoP	= infdigitLeg->GetRcvOrPay();				
					double leverage = infdigitLeg->GetMultiple();
					double constant = infdigitLeg->GetConstant();
					double spread = leverage+constant;

					ARM_Date numDate, denomDate;
					int k;
					
					//For View
					itsPricingStrikesForView	= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
					itsInflationVolsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
					itsDigitLegVolsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
					itsCorrelationsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));	
								
					for (k = 0; k<nbCapLet; k++) 
					{
						numDate			= ARM_Date( infdigitLeg->GetNumResetDates()->Elt(k) );
						denomDate		= ARM_Date( infdigitLeg->GetDenomResetDates()->Elt(k) );

						double interestTerm	= infdigitLeg->GetInterestTerms()->Elt(k);
						
						///////////the infdigitLegInfo by flow
						
						///GetFwdRates() return the (leverage*forward + spread)*interestterm 
						///where forward = CPI(i)/CPI(i-1) -1
						///here we compute the CPIForward=leverage*(forward+1)
						CPIForward		= infdigitLeg->GetFwdRates()->Elt(k) / CC_NS( ARM_Constants, rateBase )/(interestTerm*infRoP) - constant;  
						discountfactor	= infdigitLeg->GetDiscountFactors()->Elt(k);
						double payTime	= infdigitLeg->GetPaymentDates()->Elt(k);				
														
						/// to account for fixing difference!
						double	renormalisationFactor = 1.0;
						timeToStart		= pricingMod->GetModelTimeWPublishLag( denomDate, infDigitIdx );
						if( timeToStart < 0.0 )
						{
							renormalisationFactor = infdigitLeg->GetDenomCPIRates()->Elt(k) / pricingMod->GetCPIIndexValue();
						}

						double timeToExpiry		= pricingMod->GetModelTimeWPublishLag( numDate, infDigitIdx);
						tenor = timeToExpiry	- timeToStart;
						double remainingTenor	= tenor;
						if( timeToStart < 0.0 )
							remainingTenor = remainingTenor+timeToStart <0? 0 : remainingTenor+timeToStart;

						/// to handle weird case of CPIForward negative!					
						

						////////////////////////////////////compute the volatilit of the CPIFwd
						double TjVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, denomDate ); 
						ARM_Date denomDatewPublishLag = pricingMod->GetModelDateWPublishLag( denomDate, infDigitIdx );
						double TjForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, denomDatewPublishLag );

						double TiVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, numDate ); 
						double TiFromAsOf		= CountYearsWithoutException( dayCount, modelAsOfDate, numDate ); 
						ARM_Date numDatewPublishLag = pricingMod->GetModelDateWPublishLag( numDate, infDigitIdx );
						double TiForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, numDatewPublishLag );
						tenor					= TiVolLookup-TjVolLookup;

						double irInterestTerm	= irLeg->GetInterestTerms()->Elt(k);
						double irResetDate	= irLeg->GetResetDates()->Elt(k);
						double irResetTerm = (irResetDate-modelAsOfDate.GetJulian())/K_YEAR_LEN;

						if( TiForMaturity < 0.0 || tenor < 0.0 || ((tenor+TjForMaturity) < 0.0) || (irResetTerm < 0.0))
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
										ARM_USERNAME + ": reset dates are before the AsOfDate can't priced the LiborLeg we don't have a reset manager for Libor Rates !" );

						else
						{
							if( TjForMaturity <= 0.0 )
							{
								tenor		   += TjForMaturity;
								TjForMaturity	= StringMaturityToYearTerm( "1d" );
								TjVolLookup		= StringMaturityToYearTerm( "1d" );
							}
							if( TjVolLookup <= 0.0 )
							{							
								TjVolLookup		= StringMaturityToYearTerm( "1d" );
							}

							double infResetTerm = (numDatewPublishLag.GetJulian()-modelAsOfDate.GetJulian())/K_YEAR_LEN;

							if( TiForMaturity < K_NEW_DOUBLE_TOL)
							{
									TiForMaturity = StringMaturityToYearTerm( "1d" );
									infResetTerm = 1.0/365.0;
							}

								///////////the irLegInfo by flow 
							double irlegRate = irLeg->GetFwdRates()->Elt(k)/CC_NS( ARM_Constants, rateBase );
																	
							////compute the strike for the caplet
							double barrier = barrierProfile->CptReferenceValue(payTime);
							double spread_IR = spreadRefValue->CptReferenceValue(irResetDate)/CC_NS( ARM_Constants, rateBase );
							
							double pricingStrike = barrier/CC_NS( ARM_Constants, rateBase ) - constant;

							double resetLag=fabs(irResetTerm -infResetTerm)*K_YEAR_LEN;
							double totalVol = K_NEW_DOUBLE_TOL;
							double totalVol_plus = K_NEW_DOUBLE_TOL;
							double totalVol_moins = K_NEW_DOUBLE_TOL;
											
							if(leverage>1e-15)
							{
								double lookUpStrike = pow( renormalisationFactor * pricingStrike/leverage, 1.0/tenor)-1.0;
								double vol			= pricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike , tenor )/CC_NS( ARM_Constants, rateBase );
								totalVol	= vol * sqrt( tenor );

								double lookUpStrike_plus = pow( renormalisationFactor * (pricingStrike+CALL_SPREAD_SHIFT_INF)/leverage, 1.0/tenor)-1.0;
								vol			= pricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike_plus , tenor )/CC_NS( ARM_Constants, rateBase );
								totalVol_plus	= vol * sqrt( tenor );

								double lookUpStrike_moins = pow( renormalisationFactor * (pricingStrike-CALL_SPREAD_SHIFT_INF)/leverage, 1.0/tenor)-1.0;
								vol			= pricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike_moins , tenor )/CC_NS( ARM_Constants, rateBase );
								totalVol_moins	= vol * sqrt( tenor );
							}
							////////////////////////////End of compute the volatilit of the CPIFwd
							
							double irTenor = irIdx->GetYearTerm();

							if ( irLeg->GetName() == ARM_CMSLEG )
								irTenor = ((ARM_CMSLeg *) irLeg)->GetSwapYearTerm();
							

							double irmoyeness = 0.0; //ATM Vol for the pay underlying
							double irVol = pricingMod->GetIRModel()->GetVolatility()->ComputeVolatility( irResetTerm, irmoyeness,irTenor )/CC_NS( ARM_Constants, rateBase );
							
							double correl = pricingMod->GetInfIRCorrel( TiForMaturity,irTenor, infDigitIdx, irIdx, "INF_YoY/IR_FWD" );
						
							if(resetLag>RESET_LAG_INF_IR)
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
									ARM_USERNAME + ": resetlag between Libor and INF not allowed to exceed 7 days !" );
							else
							{
								if( fabs(pricingStrike) < K_NEW_DOUBLE_TOL)
									pricingStrike = K_NEW_DOUBLE_TOL;													
								
								double totalirVol =  irVol*sqrt(irResetTerm);
								double optLet;								
								double optLetDigitPaySpread;
								double optLetDigitPayFloat;
								
								optLetDigitPaySpread = spread_IR*DigitalBlackSholesSmooth_Formula( CPIForward, totalVol_plus, totalVol_moins, discountfactor, pricingStrike, callPut, CALL_SPREAD_SHIFT_INF);
								
								double CPIForward_ADJ = CPIForward*exp(correl*totalirVol*totalVol);
								optLetDigitPayFloat = irlegRate*DigitalBlackSholesSmooth_Formula( CPIForward_ADJ, totalVol_plus, totalVol_moins, discountfactor, pricingStrike, callPut, CALL_SPREAD_SHIFT_INF);
								optLet = optLetDigitPaySpread+optLetDigitPayFloat;
								double amount = notional.CptReferenceValue(payTime);				
								price		   += lonOrShort*optLet*interestTerm*amount;
							}

							//for View
							(*itsPricingStrikesForView)[k] = pricingStrike-leverage;
							(*itsInflationVolsForView)[k] = totalVol;
							(*itsDigitLegVolsForView)[k] = irVol;
							(*itsCorrelationsForView)[k] = correl;
						}
					}	
					delete irLeg;
					delete infdigitLeg;
				}//End Third case Pay IR index							
			}
			break;
			
			/// 2) case of libor vs inflation
		case K_PRODUCT:
			{
			}			
			break;
			
		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": Digital Payoff are K_SIMPLE or K_PRODUCT !");
		}		
	
	}		
	else
	{ 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		 ARM_USERNAME + ": model is of wrong type.. cannot price option on CPI.. Please advise");
	}

	SetPrice(price);
	return (price);
}



////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: View
///	Returns: void
///	Action : creates a file in which it stores the object
///			 information
////////////////////////////////////////////////////
void ARM_InfHybridDigital::View(char* id, FILE* ficOut )
{
    FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( NULL == ficOut )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	/// use the method to string
	/// and just says the type and what is in it
	fprintf(fOut, " \n ===========================================================================================================" );
	fprintf(fOut, "	\n ==========================================> INFLATION CapFloor  <==========================================" );
	fprintf(fOut, " \n ===========================================================================================================\n\n\n\n" );

	int payoffType= GetPayoffType();

	
	
	char d1[20];
		
	/// for readability split the fprintf function

	fprintf(fOut, " \n =================================================================" );
    fprintf(fOut, " \n =====================> Quick Pricing Info  <=====================" );
	fprintf(fOut, "	\n ================================================================= \n\n" );
	switch( payoffType )
	{
		/// 1) inflation vs fix
		case K_SIMPLE:
			{
				fprintf(fOut, " \n ==============> INFLATION CapFloor : FIXED STRIKE Case  <============== \n" );

				if (itsPricingStrikesForView != ARM_VectorPtr(NULL)
					&& itsInflationVolsForView != ARM_VectorPtr(NULL))
				{
					if (itsDigitLegVolsForView != ARM_VectorPtr(NULL)
					&& itsCorrelationsForView != ARM_VectorPtr(NULL))
					{
						fprintf(fOut, "\n StartDates\t Strike\t INFVol\t SecondLegVol\t Correlation \n");
						int nbCapLet = GetNumFlows();
						int  i;
						for (i = 0; i < nbCapLet; ++i )
						{					
							((ARM_Date) (*(itsPayLeg->GetFlowStartDates()))[i]).JulianToStrDate(d1);
													
							fprintf(fOut, " %s\t %.3lf\t %.3lf\t %.3lf\t\t %.2lf \n", 
								d1, (*itsPricingStrikesForView)[i]*100.0,(*itsInflationVolsForView)[i]*100.0, (*itsDigitLegVolsForView)[i]*100.0, (*itsCorrelationsForView)[i]*100.0 );						
						}
					}
					else
					{
						fprintf(fOut, "\n StartDates\t Strike\t TotalVol \n");
						int nbCapLet = GetNumFlows();
						int  i;
						for (i = 0; i < nbCapLet; ++i )
						{					
							((ARM_Date) (*(itsPayLeg->GetFlowStartDates()))[i]).JulianToStrDate(d1);
													
							fprintf(fOut, " %s\t %.3lf\t %.3lf \n", 
								d1, (*itsPricingStrikesForView)[i]*100.,(*itsInflationVolsForView)[i]*100.0 );						
						}
					}
				}
				else
				{
					fprintf(fOut, "\n -----> Not Priced Yet \n" );
				}
			}
		break;
				
			/// 2) case of libor vs inflation
		case K_PRODUCT:
			{
				ARM_SwapLeg* irLeg		= (ARM_SwapLeg*) (itsDigitLeg->Clone());

				if( dynamic_cast<ARM_InfLeg*>(irLeg) )
				{
					fprintf(fOut, " \n ==============> INFLATION CapFloor : INFLATION SPREAD Case  \n\n" );						
				}
				else
				{
					fprintf(fOut, " \n ==============> INFLATION CapFloor : INFLATION HYBRIDE Case \n\n" );
				}
				if (itsPricingStrikesForView != ARM_VectorPtr(NULL)
					&& itsInflationVolsForView != ARM_VectorPtr(NULL)
					&& itsDigitLegVolsForView != ARM_VectorPtr(NULL)
					&& itsCorrelationsForView != ARM_VectorPtr(NULL))
				{
					fprintf(fOut, "\n StartDates\t Strike\t INFVol\t SecondLegVol\t Correlation \n");
					int nbCapLet = GetNumFlows();
					int  i;
					for (i = 0; i < nbCapLet; ++i )
					{					
						((ARM_Date) (*(itsPayLeg->GetFlowStartDates()))[i]).JulianToStrDate(d1);
												
						fprintf(fOut, " %s\t %.3lf\t %.3lf\t %.3lf\t\t %.2lf \n", 
							d1, (*itsPricingStrikesForView)[i]*100.0,(*itsInflationVolsForView)[i]*100.0, (*itsDigitLegVolsForView)[i]*100.0, (*itsCorrelationsForView)[i]*100.0 );						
					}
				}
				else
				{
					fprintf(fOut, "\n -----> Not Priced Yet \n" );
				}
			}
		break;
			
		default:
			{
				fprintf(fOut, "\n ==============> No More Information \n\n" );
			}
	}
	
	
	fprintf(fOut, " \n\n\n\n ================================================================== " );
    fprintf(fOut, " \n =====================> Underlying Swap Legs <===================== " );
	fprintf(fOut, "	\n ================================================================== \n\n" );

	ARM_Swap::View(id,fOut);
	
	/// to allow to have nested view
    if ( ficOut == NULL )
       fclose(fOut);
}



////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: CptCashFlowValues
///	Returns: void
///	Action : computes the various cash flow value for the different leg!
////////////////////////////////////////////////////
void ARM_InfHybridDigital::CptCashFlowValues()
{}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/