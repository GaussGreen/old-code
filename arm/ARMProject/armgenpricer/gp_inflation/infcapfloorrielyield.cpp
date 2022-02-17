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

#include "gpinflation/infcapfloorrielyield.h"

/// gpinflation
#include "gpinflation/infcurv.h"
#include "gpinflation/infbsmodel.h"
#include "gpinflation/infmultibsmodel.h"
#include "gpinflation/infleg.h"
#include "gpinflation/infidx.h"
#include "gpinflation/assetinfo.h"


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
#include <glob/dates.h>


CC_USING_NS(std,string)
#define CC_MIN(a, b)  ((a) < (b) ? (a) : (b))

CC_BEGIN_NAMESPACE( ARM )

const double RESET_LAG_INF_IR = 5;
const double STEP_FOR_INTEGRAL = 120;

///////////////////////////////////////////////////
///	Class  : ARM_InfCapFloorRielYield
///	Routine: Constructor
///	Returns: 
///	Action : Constructor from a swap object
////////////////////////////////////////////////////
ARM_InfCapFloorRielYield::ARM_InfCapFloorRielYield( 
	ARM_Swap* swap, 
	int CoF, 
	ARM_ReferenceValue *strikeProfile)
:
	/// cannot used the swaption constructor with swap as it
	/// imposes one fixed leg at least!
	ARM_Swaption(),
	itsAssetsInfo( NULL ),
	itsFirstInfLeg( NULL ),
	itsSecondLeg( NULL ),
	itsPricingStrikesForView(NULL),
	itsInflationVolsForView(NULL),
	itsSecondLegVolsForView(NULL),
	itsCorrelationsForView(NULL)

{
	/// priority rule for rcv or Pay!
	/// ugly but hey the only way to derive from swaption without modifying the already
	/// too long code of swaption!
    ARM_Swap::Copy( (ARM_Object*) swap);
	SetName( ARM_INFCAPFLOORRIELYIELD );

	/// set the underlying swaphttp://www.guadeloupe-fr.com/at.1/navigPhotoThequPGunePhoto/PG001025.gif/
	/*double strike = strikeProfile->CptReferenceValue(0.0);
    SetExpiryDate( (ARM_Date) ARM_DEFAULT_DATE);
	SetStrike( strike );*/
	SetOptionType( CoF );
	double strike           = -1000000.0;
	SetStrikes(strikeProfile);

	if( (!(dynamic_cast<ARM_InfLeg*>(GetFirstLeg()))) && (!(dynamic_cast<ARM_InfLeg*>(GetSecondLeg()))))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": For GP_INF_CapFloor, one of the two legs should be an inflation leg!");


	/// fix and inflation
	if( GetFixedLeg() )
	{
		itsFirstInfLeg	= GetFloatLeg();
		itsSecondLeg	= GetFixedLeg();

		if(itsFirstInfLeg->GetPaymentFreq()!= itsSecondLeg->GetPaymentFreq())
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": the two legs should have the same pay Freqency and the same reset freqency!");
		ARM_SingleAssetInfo fixedLedInfo( itsSecondLeg->GetCurrencyUnit()->GetCcyName(), itsSecondLeg->GetRcvOrPay(), K_FIXED_LEG );
		ARM_InfIdx* infidx		= (ARM_InfIdx*) itsFirstInfLeg->GetIRIndex();
		ARM_SingleAssetInfo infLedInfo(  infidx->GetIndexName(), itsFirstInfLeg->GetRcvOrPay(), K_GENERICINF_LEG );
		itsAssetsInfo = new ARM_TwoAssetsInfo( infLedInfo, fixedLedInfo, strike );
	}
	/// libor and inflation
	else
	{
		if( dynamic_cast<ARM_InfLeg*>(GetFirstLeg()) )
		{
			itsFirstInfLeg	= GetFirstLeg();
			itsSecondLeg	= GetSecondLeg();
		}
		else
		{
			itsSecondLeg	= GetFirstLeg();
			itsFirstInfLeg	= GetSecondLeg();
		}
		
		if( !dynamic_cast<ARM_AverageSwapLeg*> 	(itsSecondLeg ) ){
			if((itsFirstInfLeg->GetIRIndex()->GetResetFrequency()!= itsSecondLeg->GetIRIndex()->GetResetFrequency()) || (itsFirstInfLeg->GetPaymentFreq()!= itsSecondLeg->GetPaymentFreq()))
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + ": the two legs should have the same pay Freqency and the same reset freqency!");
		}

		ARM_InfIdx* infidx		= (ARM_InfIdx*) itsFirstInfLeg->GetIRIndex();
		ARM_SingleAssetInfo infLedInfo(   infidx->GetIndexName(), itsFirstInfLeg->GetRcvOrPay(),	K_GENERICINF_LEG );
		ARM_SingleAssetInfo floatLedInfo( itsSecondLeg->GetCurrencyUnit()->GetCcyName(),   itsSecondLeg->GetRcvOrPay(),		K_FLOATING_LEG  );
		itsAssetsInfo = new ARM_TwoAssetsInfo( infLedInfo, floatLedInfo, strike );
	}

	/// set the currency according to the fixed leg
	if( strcmp( itsSecondLeg->GetCurrencyUnit()->GetCcyName(), itsFirstInfLeg->GetCurrencyUnit()->GetCcyName() ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"inflation CapFloor only allowed for legs with same interest rate currency!");
}

////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloorRielYield
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_InfCapFloorRielYield::ARM_InfCapFloorRielYield( const ARM_InfCapFloorRielYield& rhs )
:	ARM_Swaption( rhs ),
	itsAssetsInfo( rhs.itsAssetsInfo? (ARM_TwoAssetsInfo*) rhs.itsAssetsInfo->Clone() : NULL ),
	itsPricingStrikesForView(rhs.itsPricingStrikesForView),
	itsInflationVolsForView(rhs.itsInflationVolsForView),
	itsSecondLegVolsForView(rhs.itsSecondLegVolsForView),
	itsCorrelationsForView(rhs.itsCorrelationsForView)
{
	itsFirstInfLeg = rhs.itsFirstInfLeg ?  (ARM_SwapLeg*) rhs.itsFirstInfLeg->Clone() : NULL;
	itsSecondLeg = rhs.itsSecondLeg ?  (ARM_SwapLeg*) rhs.itsSecondLeg->Clone() : NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloorRielYield
///	Routine: Assignment operator
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_InfCapFloorRielYield& ARM_InfCapFloorRielYield::operator=( const ARM_InfCapFloorRielYield& rhs )
{
	if( this != &rhs )
	{
		delete itsAssetsInfo;
		ARM_Swaption::operator=( rhs );
		itsAssetsInfo	= rhs.itsAssetsInfo == NULL ? NULL : new ARM_TwoAssetsInfo( *rhs.itsAssetsInfo );
		itsFirstInfLeg	= rhs.itsFirstInfLeg;
		itsSecondLeg	= rhs.itsSecondLeg;
		itsPricingStrikesForView	= rhs.itsPricingStrikesForView;
		itsInflationVolsForView		= rhs.itsInflationVolsForView;
		itsSecondLegVolsForView		= rhs.itsSecondLegVolsForView;
		itsCorrelationsForView		= rhs.itsCorrelationsForView;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloorRielYield
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_InfCapFloorRielYield::~ARM_InfCapFloorRielYield()
{
	delete itsAssetsInfo;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloorRielYield
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_InfCapFloorRielYield::toString() const
{
	return string();
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloorRielYield
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_InfCapFloorRielYield::Clone()
{
	return new ARM_InfCapFloorRielYield(*this);
}
///////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: ComputePrice
///	Returns: ARM_Object*
///	Action : Computes the inflation swaption price
////////////////////////////////////////////////////
double ARM_InfCapFloorRielYield::ComputePrice(int)
{
	ARM_Model*	model			= GetModel();
	ARM_Date	modelAsOfDate	= model->GetStartDate();
	double		price			= 0.0;
	double		base			= CC_NS( ARM_Constants, rateBase );

	if( model->CanPriceInflation() >= PRICE_FWDNOPTION ){

		ARM_InfLeg*		infLeg = dynamic_cast< ARM_InfLeg* >(itsFirstInfLeg->Clone());
		if	(!infLeg )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": try to price an inflation cap/floor on a leg that is not an inflation leg... Please advise");

		int swapType		=	infLeg->GetSwapType();
		if( swapType		!=	K_YEARTOYEAR_LEG)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": only year on year leg are supported to price INF cap!");

		ARM_InfIdx* infidx	= (ARM_InfIdx*) infLeg->GetIRIndex();
		if( !infidx )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": null inflation index!");

		ARM_InfBSModel* pricingMod = dynamic_cast<ARM_InfBSModel*>(model);
		if( !pricingMod )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"model is of wrong type.. is supposed to price option on CPI but failed to do so because it does  not inherit from InfOptionModel.. Please advise");

		
		int						dayCount			= pricingMod -> GetInfFwdCurv() -> GetMonthlyInterpType();
		ARM_Date				lastKnownDate		= pricingMod -> GetVolatility() -> GetLastKnownDate();
		ARM_Date				modelAsOfDate		= pricingMod -> GetStartDate();
		
		int						callPut				= GetOptionType();
		string					indexName			= infidx -> GetIndexName();
		ARM_ReferenceValue*		strikeProfile		= GetStrikes();
		int						secondLegType		= itsAssetsInfo -> GetOptionType();

		ARM_ReferenceValue		notional			= *itsFirstInfLeg -> GetAmount() ;

		switch( secondLegType )
		{
			/// 1) inflation vs fix
		case K_FIXED_LEG:
			{
				ARM_InfBSModelPtr firstPricingMod = ARM_InfBSModelPtr((ARM_InfBSModel*) pricingMod->Clone());
				ARM_InfMultiBSModel* multiPricingMod = dynamic_cast<ARM_InfMultiBSModel*>(pricingMod);
				if( multiPricingMod )
					firstPricingMod = multiPricingMod->GetCorrespondingInflationModel( infidx );
				
				ARM_SwapLeg* fixLeg		= (ARM_SwapLeg*) (itsSecondLeg->Clone());
				
				double optLet, CPIForward, timeToStart, tenor, discountfactor;
				int nbCapLet = infLeg->GetPaymentDates()->GetSize(); //GetNumFlows();

				int infRoP	= infLeg->GetRcvOrPay();
				int fixRorP = fixLeg->GetRcvOrPay();			
				
				double leverage = infLeg->GetMultiple();
				double constant = infLeg->GetConstant();
				double spread = leverage+constant;

				//For View
				itsPricingStrikesForView = ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
				itsInflationVolsForView = ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));									

				ARM_Date numDate, denomDate;
				int k;				
				for (k = 0; k<nbCapLet; k++) 
				{
					numDate			= ARM_Date( infLeg->GetNumResetDates()->Elt(k) );
					denomDate		= ARM_Date( infLeg->GetDenomResetDates()->Elt(k) );

					double interestTerm	= infLeg->GetInterestTerms()->Elt(k);
					
					///////////the infLegInfo by flow
					
					///GetFwdRates() return the (leverage*forward + spread)*interestterm 
					///where forward = CPI(i)/CPI(i-1) -1
					///here we compute the CPIForward=leverage*(forward+1)
					CPIForward		= infLeg->GetFwdRates()->Elt(k) / CC_NS( ARM_Constants, rateBase )/(interestTerm*infRoP) - constant;  
					discountfactor	= infLeg->GetDiscountFactors()->Elt(k);
					double payTime	= infLeg->GetPaymentDates()->Elt(k);
				
					///////////the fixLegInfo by flow
					double fixlegRate = fixLeg->GetFwdRates()->Elt(k);
					double fixInterestTerm	= fixLeg->GetInterestTerms()->Elt(k);
					
					
					////compute the strike for the caplet
					double strike = strikeProfile->CptReferenceValue(payTime);
					strike = strike/infRoP-fixlegRate*(fixRorP*fixInterestTerm)/(infRoP*interestTerm);
					
				
													
					/// to account for fixing difference!
					double	renormalisationFactor = 1.0;
					timeToStart		= firstPricingMod->GetModelTimeWPublishLag( denomDate, infidx );
					if( timeToStart < 0.0 )
					{
						/// in this particular case we need to account for a fixing different from the one
						/// of the zero coupon option
						/// the renormalisation is YtYFixingCPI/CurrentCPI
						/// or the denomCPIRates / Reference CPI of the Curve
						renormalisationFactor = infLeg->GetDenomCPIRates()->Elt(k) / firstPricingMod->GetCPIIndexValue();
					}

					double timeToExpiry		= firstPricingMod->GetModelTimeWPublishLag( numDate, infidx);
					tenor = timeToExpiry	- timeToStart;
					double remainingTenor	= tenor;
					if( timeToStart < 0.0 )
						remainingTenor = remainingTenor+timeToStart <0? 0 : remainingTenor+timeToStart;

					/// to handle weird case of CPIForward negative!
					if( CPIForward < 0 )
					{
						strike		= -strike;
						callPut		= -callPut;
						CPIForward	= -CPIForward;
					}
					

					////compute the volatilit of the CPIFwd
					double TjVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, denomDate ); 
					ARM_Date denomDatewPublishLag = firstPricingMod->GetModelDateWPublishLag( denomDate, infidx );
					double TjForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, denomDatewPublishLag );

					double TiVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, numDate ); 
					double TiFromAsOf		= CountYearsWithoutException( dayCount, modelAsOfDate, numDate ); 
					ARM_Date numDatewPublishLag = firstPricingMod->GetModelDateWPublishLag( numDate, infidx );
					double TiForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, numDatewPublishLag );
					tenor					= TiVolLookup-TjVolLookup;

					double pricingStrike = strike/CC_NS( ARM_Constants, rateBase ) - constant;
					if( fabs(pricingStrike) < K_NEW_DOUBLE_TOL)
							pricingStrike = K_NEW_DOUBLE_TOL;

					double amount = notional.CptReferenceValue(payTime);
					double totalVol;
					
					if( (TiForMaturity < 0.0) || (tenor < 0.0) || ((tenor+TjForMaturity) < 0.0))
					{
						tenor = -1;
						totalVol = K_NEW_DOUBLE_TOL;
						optLet = BlackSholes_Formula( CPIForward, totalVol, discountfactor, pricingStrike, callPut*infRoP);

						optLet += optLet*interestTerm*amount;
						price+= optLet;
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
						
						if(leverage>1e-15)
						{
							double lookUpStrike = pow( renormalisationFactor * pricingStrike/leverage, 1.0/tenor) -1.0;
							double vol			= firstPricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike , tenor )/CC_NS( ARM_Constants, volBase );
							totalVol	= vol * sqrt( tenor );
						}

						else
						{
							totalVol = K_NEW_DOUBLE_TOL;
						}
								
						////////////////////////////End of compute the volatilit of the CPIFwd

						//BS Pricing						
						optLet = BlackSholes_Formula( CPIForward, totalVol, discountfactor, pricingStrike, callPut*infRoP);										
						price		   += optLet*interestTerm*amount;

						//for View
						(*itsPricingStrikesForView)[k] = pricingStrike-leverage;
						(*itsInflationVolsForView)[k] = totalVol;
					}
				}
				delete fixLeg;
			}
			break;
			
			/// 2) case of libor vs inflation
		case K_FLOATING_LEG:
			{
				//commun Data 
				double optLet, discountfactor;
				int nbCapLet = infLeg->GetPaymentDates()->GetSize();//GetNumFlows();
				int k;
				
				//For View
				itsPricingStrikesForView	= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
				itsInflationVolsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
				itsSecondLegVolsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));
				itsCorrelationsForView		= ARM_VectorPtr(new ARM_GP_Vector(nbCapLet));			
				
				//Data For IrLeg
				ARM_SwapLeg* irLeg		= (ARM_SwapLeg*) (itsSecondLeg->Clone());

				if( dynamic_cast<ARM_InfLeg*>(irLeg) )
				{					
					ARM_InfLeg* secondInfLeg = dynamic_cast< ARM_InfLeg* >(irLeg->Clone());
					ARM_InfIdx* infidx1	= (ARM_InfIdx*) secondInfLeg->GetIRIndex();
										
					ARM_InfMultiBSModel* multiPricingMod = dynamic_cast<ARM_InfMultiBSModel*>(pricingMod);
					if( !multiPricingMod )
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							"model is of wrong type.. is supposed to price option on CPI spreads but failed to do so because it does  not inherit from ARM_InfMultiBSModel.. Please advise");
					
					int infRoP	= infLeg->GetRcvOrPay();
					int infRoP1	= secondInfLeg->GetRcvOrPay();

					if( infRoP == infRoP1)
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							"The two inflations legs have to be 'one Rec and the other Pay'.. Please advise");

					ARM_InfBSModelPtr firstPricingMod = multiPricingMod->GetCorrespondingInflationModel( infidx );
					ARM_InfBSModelPtr secondPricingMod = multiPricingMod->GetCorrespondingInflationModel( infidx1 );

					//get model for each leg 
					dayCount		= firstPricingMod->GetInfFwdCurv()->GetMonthlyInterpType();
					lastKnownDate	= firstPricingMod->GetVolatility()->GetLastKnownDate();
					modelAsOfDate	= firstPricingMod->GetStartDate();

					int dayCount1			= secondPricingMod->GetInfFwdCurv()->GetMonthlyInterpType();
					ARM_Date lastKnownDate1	= secondPricingMod->GetVolatility()->GetLastKnownDate();
					ARM_Date modelAsOfDate1	= secondPricingMod->GetStartDate();

					
					
					/////////////Data For first InfLeg
					double CPIForward, timeToStart, tenor;									
					double leverage = infLeg->GetMultiple();
					double constant = infLeg->GetConstant();
					double spread = leverage+constant;
					ARM_Date numDate, denomDate;
					

					/////////////Data For second InfLeg
					double CPIForward1, timeToStart1, tenor1;									
					double leverage1 = secondInfLeg->GetMultiple();
					double constant1 = secondInfLeg->GetConstant();
					double spread1 = leverage1+constant1;
					ARM_Date numDate1, denomDate1;
								
					for (k = 0; k<nbCapLet; k++) 
					{
						////////////////////////////////////////////////
						///////////the FIRST infLegInfo by flow
						////////////////////////////////////////////////
						numDate			= ARM_Date( infLeg->GetNumResetDates()->Elt(k) );
						denomDate		= ARM_Date( infLeg->GetDenomResetDates()->Elt(k) );

						double interestTerm	= infLeg->GetInterestTerms()->Elt(k);						
						///here we compute the CPIForward=leverage*(forward+1)
						CPIForward		= infLeg->GetFwdRates()->Elt(k) / CC_NS( ARM_Constants, rateBase )/(interestTerm*infRoP) - constant;  
						discountfactor	= infLeg->GetDiscountFactors()->Elt(k);
						double payTime	= infLeg->GetPaymentDates()->Elt(k);				
						double	renormalisationFactor = 1.0;
						timeToStart		= firstPricingMod->GetModelTimeWPublishLag( denomDate, infidx );
						if( timeToStart < 0.0 )
						{
							renormalisationFactor = infLeg->GetDenomCPIRates()->Elt(k) / firstPricingMod->GetCPIIndexValue();
						}

						double timeToExpiry		= firstPricingMod->GetModelTimeWPublishLag( numDate, infidx);
						tenor = timeToExpiry	- timeToStart;
						double remainingTenor	= tenor;
						if( timeToStart < 0.0 )
							remainingTenor = remainingTenor+timeToStart <0? 0 : remainingTenor+timeToStart;

						
						double TjVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, denomDate ); 
						ARM_Date denomDatewPublishLag = firstPricingMod->GetModelDateWPublishLag( denomDate, infidx );
						double TjForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, denomDatewPublishLag );

						double TiVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, numDate ); 
						double TiFromAsOf		= CountYearsWithoutException( dayCount, modelAsOfDate, numDate ); 
						ARM_Date numDatewPublishLag = firstPricingMod->GetModelDateWPublishLag( numDate, infidx );
						double TiForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, numDatewPublishLag );
						tenor					= TiVolLookup-TjVolLookup;

						////////////////////////////////////////////////
						///////////the SECOND infLegInfo by flow
						////////////////////////////////////////////////
						numDate1		= ARM_Date( secondInfLeg->GetNumResetDates()->Elt(k) );
						denomDate1		= ARM_Date( secondInfLeg->GetDenomResetDates()->Elt(k) );

						double interestTerm1	= secondInfLeg->GetInterestTerms()->Elt(k);						
						CPIForward1		= secondInfLeg->GetFwdRates()->Elt(k) / CC_NS( ARM_Constants, rateBase )/(interestTerm1*infRoP1) - constant1;  
						double payTime1	= secondInfLeg->GetPaymentDates()->Elt(k);				
														
						double	renormalisationFactor1 = 1.0;
						timeToStart1		= secondPricingMod->GetModelTimeWPublishLag( denomDate1, infidx1 );
						if( timeToStart < 0.0 )
						{
							renormalisationFactor1 = secondInfLeg->GetDenomCPIRates()->Elt(k) / secondPricingMod->GetCPIIndexValue();
						}

						double timeToExpiry1		= secondPricingMod->GetModelTimeWPublishLag( numDate1, infidx1);
						tenor1 = timeToExpiry1	- timeToStart1;
						double remainingTenor1	= tenor1;
						if( timeToStart1 < 0.0 )
							remainingTenor1 = remainingTenor1+timeToStart1 <0? 0 : remainingTenor1+timeToStart1;

						double TjVolLookup1		= CountYearsWithoutException( dayCount1, lastKnownDate1, denomDate1 ); 
						ARM_Date denomDatewPublishLag1 = secondPricingMod->GetModelDateWPublishLag( denomDate1, infidx1 );
						double TjForMaturity1	= CountYearsWithoutException( dayCount1, modelAsOfDate, denomDatewPublishLag1 );

						double TiVolLookup1		= CountYearsWithoutException( dayCount1, lastKnownDate1, numDate1 ); 
						double TiFromAsOf1		= CountYearsWithoutException( dayCount1, modelAsOfDate, numDate1 ); 
						ARM_Date numDatewPublishLag1 = secondPricingMod->GetModelDateWPublishLag( numDate1, infidx1 );
						double TiForMaturity1	= CountYearsWithoutException( dayCount1, modelAsOfDate1, numDatewPublishLag1 );
						tenor1					= TiVolLookup1-TjVolLookup1;
						///////////////////////////////////////////////

						double strike = strikeProfile->CptReferenceValue(payTime);
						double pricingStrike = (strike/infRoP)/CC_NS( ARM_Constants, rateBase )-constant-constant1*(infRoP1*interestTerm1)/(infRoP*interestTerm);
						if( fabs(pricingStrike) < K_NEW_DOUBLE_TOL)
								pricingStrike = K_NEW_DOUBLE_TOL;
						
						double totalVol;
						double totalVol1;
						
						double correl = 1.0;

						if( TiForMaturity < 0.0 || tenor < 0.0 || ((tenor+TjForMaturity) < 0.0) || TiForMaturity1 < 0.0 || tenor1 < 0.0 || ((tenor1+TjForMaturity1) < 0.0))
						{
							tenor = -1;
							tenor1 = -1;
							totalVol1 = K_NEW_DOUBLE_TOL;
							totalVol = K_NEW_DOUBLE_TOL;

							double optionMaturity = 1.0;//fixtif is integrated in the total vol of each underlying
								
							double optionType = ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION;						
							/// Uses the closed formula library on the spread option valuation		
							int n			= STEP_FOR_INTEGRAL;													
							
							optLet = Export_LogNormal_SpreadOption(CPIForward, CPIForward1*(infRoP1*interestTerm1)/(-infRoP*interestTerm), totalVol, totalVol1, 
																			correl, pricingStrike, optionMaturity, infRoP*callPut, optionType, n);

							double amount = notional.CptReferenceValue(payTime);				
							price		   += optLet*interestTerm*amount*discountfactor;						
						}
						else
						{
							if( TjForMaturity <= 0.0)
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
							
							
							if( TjForMaturity1 <= 0.0 )
							{
								tenor1		   += TjForMaturity1;
								TjForMaturity1	= StringMaturityToYearTerm( "1d" );
								TjVolLookup1		= StringMaturityToYearTerm( "1d" );
							}
							if( TjVolLookup1 <= 0.0 )
							{
								TjVolLookup1		= StringMaturityToYearTerm( "1d" );
							}

							double infResetTerm1 = (numDatewPublishLag1.GetJulian()-modelAsOfDate.GetJulian())/K_YEAR_LEN;
							if( TiForMaturity1 < K_NEW_DOUBLE_TOL)
							{
									TiForMaturity1 = StringMaturityToYearTerm( "1d" );
									infResetTerm1 = 1.0/365.0;
							}

							
							
							double  strikeForInf = strike/infRoP-(CPIForward1+constant1)*(infRoP1*interestTerm1)/(infRoP*interestTerm)*CC_NS( ARM_Constants, rateBase );
							double pricingStrikeForINF = strikeForInf/CC_NS( ARM_Constants, rateBase ) - constant;

							double  strikeForInf1 = strike*interestTerm/(infRoP1*interestTerm1)-(CPIForward+constant)*(infRoP*interestTerm)/(infRoP1*interestTerm1)*CC_NS( ARM_Constants, rateBase );
							double pricingStrikeForINF1 = strikeForInf1/CC_NS( ARM_Constants, rateBase ) - constant1;
						
						

							double resetLag=(infResetTerm1 -infResetTerm)*K_YEAR_LEN;							
							
							if( CPIForward < 0 )
							{
								strikeForInf		= -strikeForInf;
								callPut		= -callPut;
								CPIForward	= -CPIForward;
							}
							if(leverage>1e-15)
							{
								double lookUpStrike = pow( renormalisationFactor * pricingStrikeForINF/leverage, 1.0/tenor)-1.0;
								double vol			= firstPricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike , tenor )/CC_NS( ARM_Constants, rateBase );
								totalVol	= vol * sqrt( tenor );
							}
							else
								totalVol = K_NEW_DOUBLE_TOL;							

							if( CPIForward1 < 0 )
							{
								strikeForInf1		= -strikeForInf1;
								CPIForward1	= -CPIForward1;
							}
							if(leverage1>1e-15)
							{
								double lookUpStrike1 = pow( renormalisationFactor1 * pricingStrikeForINF1/leverage1, 1.0/tenor1)-1.0;
								double vol1			= secondPricingMod->GetVolatility()->ComputeVolatility( TjVolLookup1, lookUpStrike1 , tenor1 )/CC_NS( ARM_Constants, rateBase );
								totalVol1	= vol1 * sqrt( tenor1 );
							}
							else
								totalVol1 = K_NEW_DOUBLE_TOL;
							
							double yoyTenor = 1.0;					
							correl = multiPricingMod->GetInfInfCorrel( TiForMaturity,yoyTenor, infidx, infidx1, "INF1_YoY/INF2_YoY" );
						
							//First Pricing Case INF and IR reset at approximatively at the same date (7days error tolerance) 
							// pricing as a spread option 
							if(fabs(resetLag)< RESET_LAG_INF_IR)
							{
								double optionMaturity = 1.0;//fixtif is integrated in the total vol of each underlying
								
								double optionType = ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION;						
								/// Uses the closed formula library on the spread option valuation		
								int n			= STEP_FOR_INTEGRAL;														
								
								optLet = Export_LogNormal_SpreadOption(CPIForward, CPIForward1*(infRoP1*interestTerm1)/(-infRoP*interestTerm), totalVol, totalVol1, 
																				correl, pricingStrike, optionMaturity, infRoP*callPut, optionType, n);

								double amount = notional.CptReferenceValue(payTime);				
								price		   += optLet*interestTerm*amount*discountfactor;
							}

							//Second Pricing Case INF reset before IR reset 
							// double integral pricing
							else if(resetLag>RESET_LAG_INF_IR) //Ti<(irResetDate-RESET_LAG_INF_IR)
							{							
								double totalVolCase2 = sqrt(infResetTerm/infResetTerm1)*totalVol1;						
								double totalvol_ForBS =  totalVol1*sqrt(infResetTerm1-infResetTerm )/sqrt(infResetTerm1);
							
								int step = STEP_FOR_INTEGRAL;
								
								double optLet;
								
								GaussLegendre_Coefficients c(step,-1,1) ; 
								double max = 7.0;
								double min = -7.0;

								double delta = (max-min)/2;
								double mid = (max+min)/2;
								
								double Sum = 0.0;
								for(int i=0;i<step;i++)
								{
									double x = mid+delta*c.get_point(i);
									double w_i = delta*c.get_weight(i);

									double fdwINF = CPIForward1*exp(-0.5*totalVolCase2*totalVolCase2+totalVolCase2*x);
									for (int j=0; j<step; j++)
									{								
										double y= mid+delta*c.get_point(j);
										double orth_y_x = (correl*x+sqrt(1-correl*correl)*y);
										double w_j = delta*c.get_weight(j);								
										double fwdIR = CPIForward*exp(-0.5*totalVol*totalVol+totalVol*orth_y_x);
										double a_k = interestTerm/(infRoP1*interestTerm1)*strike/CC_NS( ARM_Constants, rateBase )-constant1-(infRoP*interestTerm)/(infRoP1*interestTerm1)*constant;
										double strike_y = a_k-(infRoP*interestTerm)/(infRoP1*interestTerm1)*fwdIR;
										double f_xy = BlackSholes_Formula(fdwINF, totalvol_ForBS, 1.0, strike_y, callPut*infRoP);
										Sum+=w_j*f_xy*exp(-0.5*y*y)*w_i*exp(-0.5*x*x);
									}
								}
								
								optLet =1.0/(2.0*PI)*Sum;
							

								double amount = notional.CptReferenceValue(payTime);				
								price		   += discountfactor*optLet*interestTerm1*amount;
							}

							//third Pricing Case  reset IR before INF reset 
							// double integral pricing
							else //resetLag< (-RESET_LAG_INF_IR) irResetDate<(Ti-RESET_LAG_INF_IR)
							{									
								double totalVolCase2 = sqrt(infResetTerm1/infResetTerm)*totalVol;						
								double totalvol_ForBS =  totalVol*sqrt(infResetTerm-infResetTerm1)/sqrt(infResetTerm);
								
								int step = STEP_FOR_INTEGRAL;
								
								double optLet;
								
								GaussLegendre_Coefficients c(step,-1,1) ; 
								double max = 7.0;
								double min = -7.0;

								double delta = (max-min)/2;
								double mid = (max+min)/2;
								
								double Sum = 0.0;
								for(int i=0;i<step;i++)
								{
									double x = mid+delta*c.get_point(i);
									double w_i = delta*c.get_weight(i);

									double fdwINF = CPIForward*exp(-0.5*totalVolCase2*totalVolCase2+totalVolCase2*x);
									for (int j=0; j<step; j++)
									{								
										double y= mid+delta*c.get_point(j);
										double orth_y_x = (correl*x+sqrt(1-correl*correl)*y);
										double w_j = delta*c.get_weight(j);								
										double fwdIR = CPIForward1*exp(-0.5*totalVol*totalVol+totalVol*orth_y_x);
										double a_k = 1/infRoP*strike/CC_NS( ARM_Constants, rateBase )-constant-(infRoP1*interestTerm1)/(infRoP*interestTerm)*constant1;
										double strike_y = a_k-(infRoP1*interestTerm1)/(infRoP*interestTerm)*fwdIR;
										double f_xy = BlackSholes_Formula(fdwINF, totalvol_ForBS, 1.0, strike_y, callPut*infRoP);
										Sum+=w_j*f_xy*exp(-0.5*y*y)*w_i*exp(-0.5*x*x);
									}
								}
								
								optLet =1.0/(2.0*PI)*Sum;
							

								double amount = notional.CptReferenceValue(payTime);				
								price		   += discountfactor*optLet*interestTerm*amount;
							}

							//for View
							(*itsPricingStrikesForView)[k] = strike/CC_NS( ARM_Constants, rateBase )-infRoP*(constant+leverage)-(constant1+leverage1)*infRoP1*interestTerm1/interestTerm;
							(*itsInflationVolsForView)[k] = totalVol;
							(*itsSecondLegVolsForView)[k] = totalVol1;
							(*itsCorrelationsForView)[k] = correl;
						}
					}
					delete secondInfLeg;
				}
				else
				{
	
					//	IRS DATA	//
					ARM_IRIndex*			irIdx			= irLeg -> GetIRIndex();
					int						irRorP			= irLeg -> GetRcvOrPay();
					double					irTenor			= irIdx -> GetYearTerm();
					double					irlegRate;					
					double					irInterestTerm;
			
					double					irResetTerm;
					double					spread;


					ARM_ReferenceValue*		spreadRefValue	= irLeg -> GetSpreads();
					if( !spreadRefValue )	spreadRefValue	= new ARM_ReferenceValue(irLeg->GetSpread());

					
					//	INF DATA	//

					ARM_Date	numDate;
					ARM_Date	denomDate;

					double		payTime;
					double		interestTerm;
					double		CPIForward; 
					double		timeToStart;
				
					double		renormalisationFactor;
					double		timeToExpiry;
					double		remainingTenor;
					double		tenor;

					int			infRoP		= infLeg -> GetRcvOrPay();				
					double		leverage	= infLeg -> GetMultiple();
					double		constant	= infLeg -> GetConstant();
					double		coleverage	= infLeg -> GetCoMultiple();
					double		correl		= K_NEW_DOUBLE_TOL;

					double		strike;
					double		amount;
					
					//	INF VOL		//

					double		TjVolLookup;	
					ARM_Date	denomDatewPublishLag;	
					double		TjForMaturity;			

					double		TiVolLookup;	
					ARM_Date	numDatewPublishLag;
					double		TiForMaturity;	
					double		irVol;	
					double		irLevVol;	
					double		strikeForInf;
					double		pricingStrikeForINF;				
					double		resetLag;						
					double		totalVol;
					double		infResetTerm;

					ARM_GP_VectorPtr	InfIntTerm	( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	InfFwdRate	( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	IrFwdRate	( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	InfSpread	( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	IrResTerm	( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	IrIntTerm	( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	IrVolTerm	( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	IrLevVolTerm( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	PriceStrike ( new ARM_GP_Vector() );
					ARM_GP_VectorPtr	DiscFact	( new ARM_GP_Vector() );

					if ( dynamic_cast <ARM_AverageSwapLeg* > ( irLeg) ) {
						ComputeDataFromLivretA(	pricingMod,
												irLeg,
												infLeg,   
												InfIntTerm,
												InfFwdRate,
												IrFwdRate,
												InfSpread,
												IrResTerm,
												IrIntTerm,
												IrVolTerm,
												IrLevVolTerm,
												PriceStrike,
												DiscFact);
					}
					else{
						ComputeDataFromSwapLeg(	pricingMod,	
												irLeg,
												infLeg,
												InfIntTerm,
												InfFwdRate,
												IrFwdRate,
												InfSpread,
												IrResTerm,
												IrIntTerm,
												IrVolTerm,
												IrLevVolTerm,
												PriceStrike,
												DiscFact);
					}

	

					for (k = 0; k<nbCapLet; k++) {

						numDate				= ARM_Date( infLeg -> GetNumResetDates()	-> Elt(k) );
						denomDate			= ARM_Date( infLeg -> GetDenomResetDates()	-> Elt(k) );
						interestTerm		= InfIntTerm -> Elt(k);
						CPIForward			= InfFwdRate -> Elt(k);
						discountfactor		= DiscFact	 -> Elt(k);
						payTime				= infLeg	 -> GetPaymentDates()	-> Elt(k);				
														
						/// to account for fixing difference!
						renormalisationFactor		=	1.0;
						timeToStart					=	pricingMod -> GetModelTimeWPublishLag( denomDate, infidx );
						if( timeToStart < 0.0 )
							renormalisationFactor	=	infLeg -> GetDenomCPIRates() -> Elt(k)/	pricingMod -> GetCPIIndexValue();
			
						timeToExpiry				= pricingMod->GetModelTimeWPublishLag( numDate, infidx);
						tenor						= timeToExpiry	- timeToStart;
						remainingTenor				= tenor;
						if( timeToStart < 0.0 )
							remainingTenor			= remainingTenor+timeToStart <0? 0 : remainingTenor + timeToStart;

						/// to handle weird case of CPIForward negative!					
						

						////////////////////////////////////compute the volatility of the CPIFwd
						TjVolLookup				=	CountYearsWithoutException( dayCount,	lastKnownDate,	denomDate ); 
						denomDatewPublishLag	=	pricingMod->GetModelDateWPublishLag(	denomDate,		infidx );
						TjForMaturity			=	CountYearsWithoutException( dayCount,	modelAsOfDate,	denomDatewPublishLag );

						TiVolLookup				=	CountYearsWithoutException( dayCount,	lastKnownDate,	numDate ); 
						numDatewPublishLag		=	pricingMod->GetModelDateWPublishLag(	numDate,		infidx );
						TiForMaturity			=	CountYearsWithoutException( dayCount,	modelAsOfDate,	numDatewPublishLag );
				
						tenor					=	TiVolLookup-TjVolLookup;
						strike					=	strikeProfile->CptReferenceValue(payTime);
						amount					=	notional.CptReferenceValue(payTime);



						double irTenor			= irIdx->GetYearTerm();  //uniquement pour correl et fwd

						if ( irLeg->GetName() == ARM_CMSLEG )
							irTenor = ((ARM_CMSLeg *) irLeg)->GetSwapYearTerm();

						irlegRate				=	IrFwdRate -> Elt(k);					
						irInterestTerm			=	IrIntTerm -> Elt(k);
						irResetTerm				=	IrResTerm -> Elt(k);
						spread					=	InfSpread -> Elt(k);
						
						irVol					=	IrVolTerm -> Elt(k);
						irLevVol				=	IrLevVolTerm -> Elt(k);	
						
						base					=	CC_NS( ARM_Constants, rateBase );
						strikeForInf			=	PriceStrike -> Elt(k);
						pricingStrikeForINF		=	strikeForInf/base - constant;	

						if( TiForMaturity < 0.0 || tenor < 0.0 || ((tenor+TjForMaturity) < 0.0))
						{
							tenor			= -1;
							infResetTerm	= StringMaturityToYearTerm( "1d" );
							resetLag		=(irResetTerm -infResetTerm)*K_YEAR_LEN;
							totalVol		= K_NEW_DOUBLE_TOL; 
							
						}
						else
						{
							if( TjForMaturity <= 0.0 )
							{
								tenor				+=	TjForMaturity;
								TjForMaturity		=	StringMaturityToYearTerm( "1d" );
								TjVolLookup			=	StringMaturityToYearTerm( "1d" );
							}
							if( TjVolLookup <= 0.0 )
								TjVolLookup			=	StringMaturityToYearTerm( "1d" );
							
							infResetTerm			= ( numDatewPublishLag.GetJulian() - modelAsOfDate.GetJulian())/K_YEAR_LEN;
							resetLag				= ( irResetTerm -infResetTerm)*K_YEAR_LEN;

							if( TiForMaturity < K_NEW_DOUBLE_TOL)
							{
									TiForMaturity	=	StringMaturityToYearTerm( "1d" );
									infResetTerm	=	1.0/365.0;
							}											
								
							if( CPIForward < 0 )
							{
								strikeForInf		=	-strikeForInf;
								pricingStrikeForINF	=	strikeForInf/base - constant;
								callPut				=	-callPut;
								CPIForward			=	-CPIForward;
							}
							if(leverage>1e-15)
							{
								double lookUpStrike =	pow( renormalisationFactor * pricingStrikeForINF/leverage, 1.0/tenor)-1.0;
								double vol			=	pricingMod->GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike , tenor )/base;
								totalVol			=	vol * sqrt( tenor );
							}

							else
								totalVol = K_NEW_DOUBLE_TOL;
							
							////////////////////////////End of compute the volatilit of the CPIFwd							
							
							
							correl = pricingMod->GetInfIRCorrel( TiForMaturity,irTenor, infidx, irIdx, "INF_YoY/IR_FWD" );
						}									


						if(resetLag>K_NEW_DOUBLE_TOL) //Ti<(irResetDate-RESET_LAG_INF_IR)
						{		
								
						
							double totalirVol =  irLevVol*sqrt(irResetTerm -infResetTerm);
												
							int step = STEP_FOR_INTEGRAL;
						
							double optLet;
													
							GaussLegendre_Coefficients c(step,-1,1) ; 
							double max = 7.0;
							double min = -7.0;

							double delta = (max-min)/2;
							double mid = (max+min)/2;
							double Sum = 0.0;
							for(int i=0;i<step;i++)
							{
								double x = mid+delta*c.get_point(i);
								double w_i = delta*c.get_weight(i);

								double fwdIR = irlegRate/CC_NS( ARM_Constants, rateBase )*exp(-0.5*infResetTerm*irVol*irVol+irVol*sqrt(infResetTerm)*x) ;
								for (int j=0; j<step; j++)
								{								
									double y= mid+delta*c.get_point(j);
									double w_j = delta*c.get_weight(j);
									double fdwINF = CPIForward*exp(-0.5*totalVol*totalVol+totalVol*(correl*x+sqrt(1-correl*correl)*y));
									double a_k = (interestTerm/(irRorP*irInterestTerm)*strike-spread)/base-constant*(infRoP*interestTerm)/(irRorP*irInterestTerm);
									double strike_y = a_k+(-infRoP*interestTerm)/(irRorP*irInterestTerm)*fdwINF;

// on renormalise le strike et on retablie le prix du cap en multipliant par le coleverage:
									
									double f_xy = coleverage*BlackSholes_Formula( fwdIR, totalirVol, 1.0, strike_y/coleverage, callPut*irRorP);
									Sum+=w_j*f_xy*exp(-0.5*y*y)*w_i*exp(-0.5*x*x);
								}
							}
							
							optLet	= 1.0/(2.0*PI)*Sum;
							
							price	+= discountfactor*optLet*irInterestTerm*amount;

						}

						//third Pricing Case  reset IR before INF reset 
						// double integral pricing
						else //resetLag< (-RESET_LAG_INF_IR) irResetDate<(Ti-RESET_LAG_INF_IR)
						{
							double totalVolCase2 = sqrt(irResetTerm/infResetTerm)*totalVol;						
							double totalvol_ForBS =  totalVol*sqrt(infResetTerm-irResetTerm )/sqrt(infResetTerm);
							
							int step = STEP_FOR_INTEGRAL;
							
							double optLet;
							
							GaussLegendre_Coefficients c(step,-1,1) ; 
							double max = 7.0;
							double min = -7.0;

							double delta = (max-min)/2;
							double mid = (max+min)/2;
							
							double Sum = 0.0;
							for(int i=0;i<step;i++)
							{
								double x = mid+delta*c.get_point(i);
								double w_i = delta*c.get_weight(i);

								double fdwINF = CPIForward*exp(-0.5*totalVolCase2*totalVolCase2+totalVolCase2*x);
								for (int j=0; j<step; j++)
								{								
									double y= mid+delta*c.get_point(j);
									double orth_y_x = (correl*x+sqrt(1-correl*correl)*y);
									double w_j = delta*c.get_weight(j);								
									double fwdIR = irlegRate/base *exp(-0.5*irResetTerm*irVol*irVol+irVol*sqrt(irResetTerm)*orth_y_x);
									
									// ligne ajoutée pour prendre en compte le leverage du fwdIR daand le strike du call sur Inf index

									fwdIR *= coleverage;
									double a_k = (1/infRoP*strike-spread*(irRorP*irInterestTerm)/(infRoP*interestTerm))/base-constant;
									double strike_y = a_k+(irRorP*irInterestTerm)/(-infRoP*interestTerm)*fwdIR;
									double f_xy = BlackSholes_Formula(fdwINF, totalvol_ForBS, 1.0, strike_y, callPut*infRoP);
									Sum+=w_j*f_xy*exp(-0.5*y*y)*w_i*exp(-0.5*x*x);
								}
							}
							
							optLet =1.0/(2.0*PI)*Sum;
										
							price		   += discountfactor*optLet*interestTerm*amount;
						}

						//for View
						(*itsPricingStrikesForView)[k] = (strike-spread*irRorP*irInterestTerm/interestTerm)/CC_NS( ARM_Constants, rateBase )-infRoP*(constant+leverage);
						(*itsInflationVolsForView)[k] = totalVol;
						(*itsSecondLegVolsForView)[k] = irVol;
						(*itsCorrelationsForView)[k] = correl;
					}					
				}
				delete irLeg;
			}
			break;
			
		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": capfloor can be inflation vs fixed, inflation vs float!");
		}

		delete infLeg;
	
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
void ARM_InfCapFloorRielYield::View(char* id, FILE* ficOut )
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
	/*if( itsAssetsInfo )
	{
		fprintf( fOut, " Quick View\n\n");
		fprintf( fOut, "%s\n", itsAssetsInfo->toString().c_str() );
	}*/

	int secondLegType = itsAssetsInfo->GetOptionType();

	
	
	char d1[20];
		
	/// for readability split the fprintf function

	fprintf(fOut, " \n =================================================================" );
    fprintf(fOut, " \n =====================> Quick Pricing Info  <=====================" );
	fprintf(fOut, "	\n ================================================================= \n\n" );
	switch( secondLegType )
	{
		/// 1) inflation vs fix
		case K_FIXED_LEG:
			{
				fprintf(fOut, " \n ==============> INFLATION CapFloor : FIXED STRIKE Case  <============== \n" );

				if (itsPricingStrikesForView != ARM_VectorPtr(NULL)
					&& itsInflationVolsForView != ARM_VectorPtr(NULL))
				{
					fprintf(fOut, "\n StartDates\t Strike\t TotalVol \n");
					int nbCapLet = itsFirstInfLeg->GetPaymentDates()->GetSize();//GetNumFlows();
					int  i;
					for (i = 0; i < nbCapLet; ++i )
					{					
						((ARM_Date) (*(itsFirstInfLeg->GetFlowStartDates()))[i]).JulianToStrDate(d1);
												
						fprintf(fOut, " %s\t %.3lf\t %.3lf \n", 
							d1, (*itsPricingStrikesForView)[i]*100.,(*itsInflationVolsForView)[i]*100.0 );						
					}
				}
				else
				{
					fprintf(fOut, "\n -----> Not Priced Yet \n" );
				}
			}
		break;
				
			/// 2) case of libor vs inflation
		case K_FLOATING_LEG:
			{
				ARM_SwapLeg* irLeg		= (ARM_SwapLeg*) (itsSecondLeg->Clone());

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
					&& itsSecondLegVolsForView != ARM_VectorPtr(NULL)
					&& itsCorrelationsForView != ARM_VectorPtr(NULL))
				{
					fprintf(fOut, "\n StartDates\t Strike\t INFVol\t SecondLegVol\t Correlation \n");
					int nbCapLet = irLeg->GetPaymentDates()->size();;
					int  i;
					for (i = 0; i < nbCapLet; ++i )
					{					
						((ARM_Date) (*(itsFirstInfLeg->GetFlowStartDates()))[i]).JulianToStrDate(d1);
												
						fprintf(fOut, " %s\t %.3lf\t %.3lf\t %.3lf\t\t %.2lf \n", 
							d1, (*itsPricingStrikesForView)[i]*100.0,(*itsInflationVolsForView)[i]*100.0, (*itsSecondLegVolsForView)[i]*100.0, (*itsCorrelationsForView)[i]*100.0 );						
					}
				}
				else
				{
					fprintf(fOut, "\n -----> Not Priced Yet \n" );
				}

				delete irLeg;
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
void ARM_InfCapFloorRielYield::CptCashFlowValues()
{}

void ARM_InfCapFloorRielYield::ComputeDataFromSwapLeg(	ARM_InfBSModel*			pricingMod,
												   		ARM_SwapLeg*			irLeg,
														ARM_InfLeg*				infLeg,											
														ARM_GP_VectorPtr		InfIntTerm,
														ARM_GP_VectorPtr		InfFwdRate,
														ARM_GP_VectorPtr		IrFwdRate,
														ARM_GP_VectorPtr		InfSpread,
														ARM_GP_VectorPtr		IrResTerm,
														ARM_GP_VectorPtr		IrIntTerm,
														ARM_GP_VectorPtr		IrVolTerm,
														ARM_GP_VectorPtr		IrLevVolTerm,
														ARM_GP_VectorPtr		PriceStrike,
														ARM_GP_VectorPtr		DiscFact	){

	double		tmp;
	double		tmpTenor;
	double		tmpStrike;
	ARM_Date	tmpDate;

	ARM_Date	InitDate	=	GetModel()->GetStartDate();	
	double		base		=	CC_NS( ARM_Constants, rateBase );
	double		constant	=	infLeg	-> GetConstant();
	double		coleverage	=	infLeg	-> GetCoMultiple();
	int			nbCaplet	=	infLeg	-> GetPaymentDates() -> GetSize();
	int			infRorP		=	infLeg  -> GetRcvOrPay();
	int			irRorP		=	irLeg   -> GetRcvOrPay();

	InfIntTerm	->	resize(nbCaplet);
	InfFwdRate	->	resize(nbCaplet);
	InfSpread	->	resize(nbCaplet);

	IrFwdRate	->	resize(nbCaplet);
	IrResTerm	->	resize(nbCaplet);
	IrIntTerm	->	resize(nbCaplet);
	IrVolTerm	->	resize(nbCaplet);
	IrLevVolTerm->	resize(nbCaplet);
	PriceStrike	->	resize(nbCaplet);
	DiscFact	->	resize(nbCaplet);


	ARM_ReferenceValue*			spreadRefValue	= irLeg -> GetSpreads();
	if( !spreadRefValue )		spreadRefValue	= new ARM_ReferenceValue(irLeg->GetSpread());

	for ( int k=0; k< nbCaplet; k++){

		InfIntTerm -> Elt(k)	=	infLeg -> GetInterestTerms() -> Elt(k);

	///////////the infLegInfo by flow
						
///GetFwdRates() return the (leverage*forward + spread)*interestterm 
///where forward = CPI(i)/CPI(i-1) -1
///here we compute the CPIForward=leverage*(forward+1)

		tmp						=	infLeg ->	GetFwdRates() -> Elt(k) / base;
		tmp						/=	infRorP *	InfIntTerm	  -> Elt(k);
		tmp						-=	constant;  
		InfFwdRate -> Elt(k)	=	tmp;
	
		IrFwdRate  -> Elt(k)	=	irLeg -> GetFwdRates() -> Elt(k);
				
		tmpDate					=	irLeg -> GetResetDates() -> Elt(k);
		InfSpread  -> Elt(k)	=	spreadRefValue -> CptReferenceValue( tmpDate	);
		IrResTerm  -> Elt(k)	=	( tmpDate.GetJulian() -  InitDate.GetJulian() )/K_YEAR_LEN;

		IrIntTerm -> Elt(k)		=	irLeg -> GetInterestTerms() -> Elt(k);

		tmpDate					=	infLeg -> GetPaymentDates()	-> Elt(k);				
		tmpStrike				=	GetStrikes() -> CptReferenceValue( tmpDate );

		tmpTenor				=	irLeg -> GetIRIndex() -> GetYearTerm();
		if(irLeg->GetName()		==	ARM_CMSLEG )
				tmpTenor		=	((ARM_CMSLeg *) irLeg) -> GetSwapYearTerm();

		tmp						=	tmpStrike; 
		tmp						-=	infRorP * base * ( InfFwdRate -> Elt(k)  + constant );
		tmp						*=	irRorP * ( InfIntTerm -> Elt(k) / IrIntTerm -> Elt(k) );
		tmp						-=	IrFwdRate -> Elt(k) + InfSpread -> Elt(k);		
		IrVolTerm -> Elt(k)		=	pricingMod->GetIRModel()->GetVolatility()->ComputeVolatility( IrResTerm -> Elt(k), tmp, tmpTenor )/base;


		tmp						=	tmpStrike; 
		tmp						-=	infRorP * base * ( InfFwdRate -> Elt(k)  + constant );
		tmp						*=	irRorP * ( InfIntTerm -> Elt(k) / IrIntTerm -> Elt(k) );
		tmp						-=	InfSpread -> Elt(k);
		tmp						/=	coleverage;
		tmp						-=	IrFwdRate -> Elt(k);		
		IrLevVolTerm -> Elt(k)	=	pricingMod->GetIRModel()->GetVolatility()->ComputeVolatility( IrResTerm -> Elt(k), tmp, tmpTenor )/base;


		tmp						=	tmpStrike *infRorP;
		tmp						-=	irRorP*infRorP*( coleverage *IrFwdRate -> Elt(k) + InfSpread -> Elt(k) )*( IrIntTerm -> Elt(k) )/( InfIntTerm -> Elt(k) );
		PriceStrike->Elt(k)		=	tmp;

		DiscFact -> Elt(k)		=	infLeg -> GetDiscountFactors()-> Elt(k);
	}
}






void ARM_InfCapFloorRielYield::ComputeDataFromLivretA(	ARM_InfBSModel*			pricingMod,
												   		ARM_SwapLeg*			irLeg,
														ARM_InfLeg*				infLeg,											   
														ARM_GP_VectorPtr		InfIntTerm,
														ARM_GP_VectorPtr		InfFwdRate,
														ARM_GP_VectorPtr		IrFwdRate,
														ARM_GP_VectorPtr		InfSpread,
														ARM_GP_VectorPtr		IrResTerm,
														ARM_GP_VectorPtr		IrIntTerm,
														ARM_GP_VectorPtr		IrVolTerm,
														ARM_GP_VectorPtr		IrLevVolTerm,
														ARM_GP_VectorPtr		PriceStrike,
														ARM_GP_VectorPtr		DiscFact	){

	const double TOL= 1e-10;

	double		tmp;
	double		tmp1;
	double		tmp2;
	double		tmpTenor;
	double		tmpStrike;
	ARM_Date	tmpDate;
	ARM_Date	InitDate;

	ARM_Date	AsOfDate	=	GetModel()	-> GetStartDate();

	int			nbResDate	=	irLeg		-> GetResetDates()	-> GetSize();
	int			nbCaplet	=	infLeg		-> GetPaymentDates()-> GetSize();
	double		constant	=	infLeg		-> GetConstant();
	double		coleverage	=	infLeg		-> GetCoMultiple();
	int			infRorP		=	infLeg		-> GetRcvOrPay();
	int			irRorP		=	irLeg		-> GetRcvOrPay();

	double		base		=	CC_NS( ARM_Constants, rateBase );

	vector< vector< double >	>	V_IrIntTerm; 
	vector< vector< double >	>	V_IrFwdRate; 
	vector< vector< ARM_Date >	>	V_IrResDate; 

	V_IrIntTerm.clear();
	V_IrFwdRate.clear();
	V_IrResDate.clear();

	vector< double	>				V_tmpVol;
	vector< double	>				V_tmpLevVol;
	vector< double	>				V_tmpInt;
	vector< double	>				V_tmpFwd;
	vector< ARM_Date>				V_tmpRes;

	V_tmpInt.clear();
	V_tmpFwd.clear();
	V_tmpRes.clear();

	double							tmpExpFwd;
	double							tmpVarFwd;		
	double							tmpFwdRate;
	double							tmpIntTerm;
	double							tmpIrResDate;
	double							tmpIrVolTerm;
	double							tmpIrLevVolTerm;
	double							tmpIrSpread;
	double							tmpCouru;

	ARM_Date						tmpResDate;
	ARM_Date						listDate;
	ARM_Date						tmpStartDate;

	int								tmpMonth;
	int								tmpYear;	
	int								tmpNb;

	InitDate		=	AsOfDate;
	tmpIntTerm		=	0;


	for (int i = 0; i<nbCaplet; i++) {

		tmpDate		=	infLeg->GetFlowStartDates()->Elt(i);
		tmpMonth	=	tmpDate.GetMonth();
		tmpYear		=	tmpDate.GetYear();
	
		if	( tmpMonth <7 )	tmpStartDate = ARM_Date( 1, 12,	tmpYear-1	); 
		else				tmpStartDate = ARM_Date( 1, 6,	tmpYear		); 

		for (int j = 0; j<nbResDate; j++) {

			listDate		=	irLeg -> GetResetDates()	->Elt(j);

			if(			tmpStartDate.GetJulian()	<=  listDate.GetJulian()	){
				if(		tmpStartDate.GetMonth()		==	listDate.GetMonth() 
													&&	tmpStartDate.GetYear()		==	listDate.GetYear()		){
					if( InitDate != listDate){

						tmpResDate	= listDate;		
						tmpFwdRate	= irLeg->GetFwdRates()->Elt(j);
						tmpIntTerm	= (tmpResDate.GetJulian() -AsOfDate.GetJulian())/K_YEAR_LEN;

						V_tmpRes.push_back(tmpResDate);
						V_tmpFwd.push_back(tmpFwdRate);
						V_tmpInt.push_back(tmpIntTerm);

						InitDate = listDate;
						}
					}
					else
						break;
				}
			}

			V_IrIntTerm.push_back(	V_tmpInt	);	
			V_IrFwdRate.push_back(	V_tmpFwd	);
			V_IrResDate.push_back(	V_tmpRes	);

			tmpIntTerm	 = 0;
			V_tmpInt.clear();
			V_tmpFwd.clear();
			V_tmpRes.clear();
		}

	InfIntTerm	->	resize(nbCaplet);
	InfFwdRate	->	resize(nbCaplet);
	InfSpread	->	resize(nbCaplet);

	IrFwdRate	->	resize(nbCaplet);
	IrResTerm	->	resize(nbCaplet);
	IrIntTerm	->	resize(nbCaplet);
	IrVolTerm	->	resize(nbCaplet);
	IrLevVolTerm->	resize(nbCaplet);
	PriceStrike	->	resize(nbCaplet);
	DiscFact	->	resize(nbCaplet);

	ARM_ReferenceValue*			spreadRefValue	= irLeg -> GetSpreads();
	if( !spreadRefValue )		spreadRefValue	= new ARM_ReferenceValue(irLeg->GetSpread());

	for ( int k=0; k< nbCaplet; k++){

		tmpCouru				=	0.0;
		
		tmp						=	infLeg -> GetInterestTerms() -> Elt(k);
		InfIntTerm -> Elt(k)	=	tmp;
		IrIntTerm -> Elt(k)		=	tmp;

		tmp						=	infLeg ->	GetFwdRates() -> Elt(k) / base;
		tmp						/=	infRorP *	InfIntTerm	  -> Elt(k);
		tmp						-=	constant;  
		InfFwdRate -> Elt(k)	=	tmp;

		tmpDate					=	infLeg -> GetPaymentDates()	-> Elt(k);				
		tmpStrike				=	GetStrikes() -> CptReferenceValue( tmpDate );
		
		tmpTenor				=	irLeg -> GetIRIndex() -> GetYearTerm();
		if(irLeg->GetName()		== ARM_CMSLEG )
				tmpTenor		= ((ARM_CMSLeg *) irLeg) -> GetSwapYearTerm();
		tmpNb					=	V_IrIntTerm[k].size();
	
		V_tmpVol.clear();
		V_tmpLevVol.clear();

		for (	i=0; i<tmpNb;	i++){

			tmpIrResDate		=	(V_IrResDate[k][i].GetJulian() - AsOfDate.GetJulian() )/K_YEAR_LEN;
			tmpIrSpread			=	spreadRefValue->CptReferenceValue( tmpIrResDate );

			if( tmpIrResDate >0 ){
				tmp					=  tmpStrike; 
				tmp					-= infRorP * irRorP * base * ( InfFwdRate -> Elt(k)  + constant );
				tmp					-= tmpIrSpread;
				tmp1				=- V_IrFwdRate[k][i] + tmp;		
				tmpIrVolTerm		=  pricingMod->GetIRModel()->GetVolatility()->ComputeVolatility( tmpIrResDate, tmp1 ,tmpTenor )/base;

				tmp2				=- V_IrFwdRate[k][i] + tmp/coleverage;	
				tmpIrLevVolTerm		=	pricingMod->GetIRModel()->GetVolatility()->ComputeVolatility( tmpIrResDate, tmp2 ,tmpTenor )/base;
			
			}
			else{
				tmpCouru = dynamic_cast<ARM_AverageSwapLeg*> ( irLeg)->GetCouru();
				if (tmpCouru==0.0)
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
					ARM_USERNAME + ": Since The LivretA Index starts before the AsofDate: we should input a non zero couru!");
				
				tmpIrVolTerm		= 0.0;
				tmpIrLevVolTerm		= 0.0;
			}
			V_tmpVol.push_back( tmpIrVolTerm );
			V_tmpLevVol.push_back( tmpIrLevVolTerm );
		}
		

	tmpExpFwd				=	tmpCouru;
	tmpVarFwd				=	0;
	for (	i=0; i<tmpNb;	i++){
		for (int	j=0; j<tmpNb;	j++){
			tmp				=	V_IrFwdRate[k][i]*V_IrFwdRate[k][j];
			tmpVarFwd		+=	tmp*exp( V_tmpVol[i]*V_tmpVol[j]*CC_MIN(V_IrIntTerm[k][i],V_IrIntTerm[k][j]) );
		}
		tmpExpFwd += V_IrFwdRate[k][i];
	}

	IrFwdRate -> Elt(k)	=	tmpExpFwd/tmpNb;
	if( V_IrIntTerm[k][tmpNb-1]<0.0 )
		IrResTerm -> Elt(k) = 0.0;
	else
		IrResTerm -> Elt(k)	=	V_IrIntTerm[k][tmpNb-1];

	tmpDate				=	V_IrResDate[k][tmpNb-1].GetJulian();
	InfSpread -> Elt(k)	=	spreadRefValue->CptReferenceValue( tmpDate );

	if ( tmpVarFwd == 0.0 )
		IrVolTerm -> Elt(k)	= 0.0;

	else if ( tmpVarFwd/((tmpExpFwd -tmpCouru )*(tmpExpFwd-tmpCouru)) <1.0-TOL )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		ARM_USERNAME + ": Problem concerning the determination on the vol of livret A");
	
	else if (IrResTerm -> Elt(k) !=0 ) {
		tmp					= 	log( tmpVarFwd/( (tmpExpFwd-tmpCouru) * (tmpExpFwd -tmpCouru)) )/ IrResTerm -> Elt(k) ;
		tmp					=	sqrt( fabs(tmp) );
		IrVolTerm -> Elt(k)	=	tmp;
	}
	else
		IrVolTerm -> Elt(k)	=	0.0;
	
	tmpExpFwd				=	tmpCouru;
	tmpVarFwd				=	0;

	for (	i=0; i<tmpNb;	i++){
		for (int	j=0; j<tmpNb;	j++){
			tmp				=	V_IrFwdRate[k][i]*V_IrFwdRate[k][j];
			tmpVarFwd		+=	tmp*exp( V_tmpLevVol[i]*V_tmpLevVol[j]*CC_MIN(V_IrIntTerm[k][i],V_IrIntTerm[k][j]) );
		}
		tmpExpFwd += V_IrFwdRate[k][i];
	}
	if ( tmpVarFwd == 0.0 )
		IrLevVolTerm -> Elt(k)	= 0.0;
	
	else if ( tmpVarFwd/(( tmpExpFwd -tmpCouru )*(tmpExpFwd-tmpCouru))<1.0-TOL  )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		ARM_USERNAME + ": Problem concerning the determination on the vol of livret A");

	else if (IrResTerm -> Elt(k) !=0 ) {
		tmp					= 	log( tmpVarFwd/( (tmpExpFwd -tmpCouru )*(tmpExpFwd-tmpCouru) ) )/ IrResTerm -> Elt(k) ;
		tmp					=	sqrt( fabs(tmp) );
		IrLevVolTerm->Elt(k)=	tmp;
	}
	else
		IrVolTerm -> Elt(k)	=	0.0;


	tmp					=	tmpStrike *infRorP;
	tmp					-=	irRorP*infRorP*( coleverage*IrFwdRate -> Elt(k) + InfSpread -> Elt(k) )*( IrIntTerm -> Elt(k) )/( InfIntTerm -> Elt(k) );
	PriceStrike->Elt(k) =	tmp;

	DiscFact -> Elt(k)	= infLeg -> GetDiscountFactors()-> Elt(k);

	}

}


double ARM_InfCallSpreadYield::ComputePrice(int){ 
		
		double price;
		double initPrice = ARM_InfCapFloorRielYield::ComputePrice();

		ARM_ReferenceValue*		tmpstrike			= GetStrikes();

		*tmpstrike += (double) EPSILON;
		SetStrikes( tmpstrike);
		
		price = -( ARM_InfCapFloorRielYield::ComputePrice()-initPrice )/EPSILON;

		*tmpstrike +=  - (double) EPSILON;
		SetStrikes( tmpstrike );

		return price;
}








































CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/



