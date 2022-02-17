/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *	\file infswaption.cpp
 *
 *  \brief inflation swpation object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpinflation/infswaption.h"

/// gpinflation
#include "gpinflation/infcurv.h"
#include "gpinflation/infbsmodel.h"
#include "gpinflation/infleg.h"
#include "gpinflation/infidx.h"
#include "gpinflation/assetinfo.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"

///gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"

/// kernel
#include <inst/swaption.h>
#include <glob/dates.h>


CC_USING_NS(std,string)


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : InfSwopt_StoreVolInfo
///	Routine: class to store all the data used
///				in the pricing step
////////////////////////////////////////////////////
class InfSwopt_StoreVolInfo : public StoreInfoObj
{
public:
	InfSwopt_StoreVolInfo( ARM_InfSwaption* infSwaption )
	:	itsInfSwaption( infSwaption )
	{}

	virtual void Store( double* data )
	{
		double fwd1				= data[0];
		double vol1				= data[1];
		double maturity1		= data[2];
		double fwd2				= data[3];
		double vol2				= data[4];
		double maturity2		= data[5];
		double correlation		= data[6];
		double annuity			= data[7];
		double optionMaturity	= data[8];
		double pricingStrike	= data[9];

		itsInfSwaption->GetAssetsInfo()->SetFwdAndVolAsset1( fwd1, vol1 );
		itsInfSwaption->GetAssetsInfo()->SetFwdAndVolAsset2( fwd2, vol2 );
		itsInfSwaption->GetAssetsInfo()->SetCorrelation( correlation );
	}

private:
	ARM_InfSwaption* itsInfSwaption;
};


////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: Constructor
///	Returns: 
///	Action : Constructor from a swap object
////////////////////////////////////////////////////
ARM_InfSwaption::ARM_InfSwaption( 
	ARM_Swap* swap, 
	int rcvPay, 
	double strike, 
	const ARM_Date& expiryDate )
:
	/// cannot used the swaption constructor with swap as it
	/// imposes one fixed leg at least!
	ARM_Swaption(),
	itsAssetsInfo( NULL ),
	itsFirstInfLeg( NULL ),
	itsSecondLeg( NULL )
{
	/// priority rule for rcv or Pay!
	size_t nbInfLeg = ARM_InfSwaption::NbOfInflationLeg( swap );
	if( !nbInfLeg )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": found no inflation leg in the swap!");

	/// ugly but hey the only way to derive from swaption without modifying the already
	/// too long code of swaption!
    ARM_Swap::Copy( (ARM_Object*) swap);
	SetName( ARM_INFSWAPTION );

	/// set the underlying swap
    SetExpiryDate( (ARM_Date&) expiryDate);
	SetStrike( strike );
	SetOptionType( rcvPay );

	switch( nbInfLeg )
	{
	case 1:
		{
			/// fix and inflation
			if( GetFixedLeg() )
			{
				itsFirstInfLeg	= GetFloatLeg();
				itsSecondLeg	= GetFixedLeg();

				itsSecondLeg->SetInitRcvOrPay(rcvPay);
				itsSecondLeg->SetFixedRate( strike );

				/// the SetInitRcvOrPay does compute the cash flow, hence used instead of SetRcvOrPay!
				itsFirstInfLeg->SetInitRcvOrPay(-rcvPay);
				ARM_SingleAssetInfo fixedLedInfo( itsSecondLeg->GetCurrencyUnit()->GetCcyName(), rcvPay, K_FIXED_LEG );
				ARM_InfIdx* infidx		= (ARM_InfIdx*) itsFirstInfLeg->GetIRIndex();
				ARM_SingleAssetInfo infLedInfo(  infidx->GetIndexName(), -rcvPay, K_GENERICINF_LEG );
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
				
				if( itsFirstInfLeg->GetRcvOrPay() * itsSecondLeg->GetRcvOrPay() > 0 )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						ARM_USERNAME + ": one leg should be receive and one pay!");

				/// the SetInitRcvOrPay does compute the cash flow, hence used instead of SetRcvOrPay!
				itsFirstInfLeg->SetInitRcvOrPay( rcvPay );
				itsSecondLeg->SetInitRcvOrPay( -rcvPay );

				ARM_InfIdx* infidx		= (ARM_InfIdx*) itsFirstInfLeg->GetIRIndex();
				ARM_SingleAssetInfo infLedInfo(   infidx->GetIndexName(), itsFirstInfLeg->GetRcvOrPay(),	K_GENERICINF_LEG );
				ARM_SingleAssetInfo floatLedInfo( itsSecondLeg->GetCurrencyUnit()->GetCcyName(),   itsSecondLeg->GetRcvOrPay(),		K_FLOATING_LEG  );
				itsAssetsInfo = new ARM_TwoAssetsInfo( infLedInfo, floatLedInfo, strike );
			}
		}
		break;
		
	case 2:
		/// inflation inflation
		{
			itsFirstInfLeg	= GetFirstLeg();
			itsSecondLeg	= GetSecondLeg();

			itsFirstInfLeg->SetInitRcvOrPay( rcvPay );
			itsSecondLeg->SetInitRcvOrPay( -rcvPay );


			ARM_InfIdx* infidx1	= (ARM_InfIdx*) itsFirstInfLeg->GetIRIndex();
			ARM_InfIdx* infidx2	= (ARM_InfIdx*) itsSecondLeg->GetIRIndex();

			ARM_SingleAssetInfo infLedInfo1(   infidx1->GetIndexName(),	itsFirstInfLeg->GetRcvOrPay(),	K_GENERICINF_LEG );
			ARM_SingleAssetInfo infLedInfo2(   infidx2->GetIndexName(),	itsSecondLeg->GetRcvOrPay(),	K_GENERICINF_LEG );
			itsAssetsInfo = new ARM_TwoAssetsInfo( infLedInfo1, infLedInfo2, strike );

			ARM_InfLeg* secondInfLeg = dynamic_cast<ARM_InfLeg*>( itsSecondLeg );
			if( !secondInfLeg )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Both legs are supposed to be inflation ");

			if( secondInfLeg->GetNumPublicationDates()->Elt(0) < expiryDate.GetJulian() )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Expiry must be before first cash flow is paid. ");
		}
		break;
	}

	ARM_InfLeg* firstInfLeg = dynamic_cast<ARM_InfLeg*>( itsFirstInfLeg );

	if( firstInfLeg->GetNumPublicationDates() && firstInfLeg->GetNumPublicationDates()->Elt(0) < expiryDate.GetJulian() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Expiry must be before first cash flow is paid. ");


	/// set the currency according to the fixed leg
	if( strcmp( itsSecondLeg->GetCurrencyUnit()->GetCcyName(), itsFirstInfLeg->GetCurrencyUnit()->GetCcyName() ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"inflation swaption only allowed for legs with same interest rate currency!");
}


	
////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: NbOfInflationLeg
///	Returns: 
///	Action : for a given swap computes the nb of inflation leg!
////////////////////////////////////////////////////
size_t ARM_InfSwaption::NbOfInflationLeg( ARM_Swap* swap )
{
	size_t NbOfInflationLeg = 0;
	ARM_SwapLeg* leg = swap->GetFirstLeg();
	if( dynamic_cast<ARM_InfLeg*>(leg) )
		++NbOfInflationLeg;
	leg = swap->GetSecondLeg();
	if( dynamic_cast<ARM_InfLeg*>(leg) )
		++NbOfInflationLeg;
	return NbOfInflationLeg;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_InfSwaption::ARM_InfSwaption( const ARM_InfSwaption& rhs )
:	ARM_Swaption( rhs ),
	itsAssetsInfo( rhs.itsAssetsInfo? (ARM_InfSwaptionContext*) rhs.itsAssetsInfo->Clone() : NULL )
	
{
	itsFirstInfLeg = rhs.itsFirstInfLeg ?  (ARM_SwapLeg*) rhs.itsFirstInfLeg->Clone() : NULL;
	itsSecondLeg = rhs.itsSecondLeg ?  (ARM_SwapLeg*) rhs.itsSecondLeg->Clone() : NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: Assignment operator
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_InfSwaption& ARM_InfSwaption::operator=( const ARM_InfSwaption& rhs )
{
	if( this != &rhs )
	{
		delete itsAssetsInfo;
		ARM_Swaption::operator=( rhs );
		itsAssetsInfo	= rhs.itsAssetsInfo == NULL ? NULL : new ARM_TwoAssetsInfo( *rhs.itsAssetsInfo );
		itsFirstInfLeg	= rhs.itsFirstInfLeg;
		itsSecondLeg	= rhs.itsSecondLeg;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_InfSwaption::~ARM_InfSwaption()
{
	delete itsAssetsInfo;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_InfSwaption::toString() const
{
	return string();
}


////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_InfSwaption::Clone()
{
	return new ARM_InfSwaption(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfSwaption
///	Routine: ComputePrice
///	Returns: ARM_Object*
///	Action : Computes the inflation swaption price
////////////////////////////////////////////////////
double ARM_InfSwaption::ComputePrice(int)
{
	ARM_Model* model		= GetModel();

	/// do not confuse model->GetStartDate that is equal to the as of date
	/// where as GetStartDate() gives the start date of the swap leg
	ARM_Date modelAsOfDate	= model->GetStartDate();
	double price			= 0.0;
	bool success;

	if( model->CanPriceInflation() >= PRICE_FWDNOPTION )
	{
		InfOptionModel* pricingMod = dynamic_cast<InfOptionModel*>(model);
		if( !pricingMod )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"model is of wrong type.. is supposed to price option on CPI but failed to do so because it does  not inherit from InfOptionModel.. Please advise");

		/// discriminate the various cases!
		int secondLegType = itsAssetsInfo->GetOptionType();
		switch( secondLegType )
		{
			/// 1) inflation vs fix
		case K_FIXED_LEG:
			{
				ARM_InfIdx* infidx		= (ARM_InfIdx*) itsFirstInfLeg->GetIRIndex();
				ARM_InfBSModel* PartialMod = dynamic_cast<ARM_InfBSModel*>(pricingMod);
				
				double asOf			= PartialMod->GetZeroCurve()->GetAsOfDate().GetJulian();
				int dayCount		= PartialMod->GetInfFwdCurv()->GetMonthlyInterpType();

				ARM_InfLeg* firstInfLeg = dynamic_cast< ARM_InfLeg* >(itsFirstInfLeg->Clone());
				ARM_SwapLeg* fixLeg		= (ARM_SwapLeg*) (itsSecondLeg->Clone());

				int flowSize = firstInfLeg->GetNumResetDates()->size();
				ARM_Date lastNumDate = ARM_Date((*firstInfLeg->GetNumResetDates())[flowSize-1]);
				ARM_Date firstDenomDate = ARM_Date((*firstInfLeg->GetDenomResetDates())[0]);

				double numFromAsof	= CountYearsWithoutException( dayCount, asOf, lastNumDate );
				double denomFromAsof	= CountYearsWithoutException( dayCount, asOf, firstDenomDate );


				double multiplier = firstInfLeg->GetMultiple();
				double constant = firstInfLeg->GetConstant();
				double spread = multiplier+constant;

				double newConstant	= constant-spread; //as if spread is zero

				firstInfLeg->SetConstant( newConstant );
				firstInfLeg->PropagateModel( PartialMod );
				firstInfLeg->SetModel( PartialMod );
				double infWithOutSpread = firstInfLeg->ComputePrice();				

				fixLeg->SetFixedRate(1.0);
				fixLeg->SetModel( PartialMod);
				double fixLegPVBP		= fixLeg->ComputePrice()* fixLeg->GetRcvOrPay();

				ARM_Swap* swapSpread_Null		= new ARM_Swap( firstInfLeg, fixLeg );
				swapSpread_Null->SetModel( PartialMod );
				double targetPrice		= 0.0;
				double CPISwapForward_WithoutSpread	= swapSpread_Null->PriceToRate( modelAsOfDate, targetPrice )/ CC_NS( ARM_Constants, rateBase );

				double CPISwapForward = ARM_Swap::PriceToRate( modelAsOfDate, targetPrice ) / CC_NS( ARM_Constants, rateBase );

				double swapSpread = CPISwapForward-CPISwapForward_WithoutSpread;

				//Compute FloatO1
				ARM_InfLeg* infLegSpread1 = dynamic_cast< ARM_InfLeg* >(firstInfLeg->Clone());
				double spread1 = 0.01;
				double newConstant_spread1 = newConstant + spread1; 
				infLegSpread1->SetConstant( newConstant_spread1 );
				infLegSpread1->PropagateModel( PartialMod );
				infLegSpread1->SetModel( PartialMod );
				double infWithSpread1 = infLegSpread1->ComputePrice();
				double floatLegPVBP = infWithSpread1-infWithOutSpread;
				double coeffSpread = multiplier*fabs(floatLegPVBP/fixLegPVBP);
				
				double strike = GetStrike()/ CC_NS( ARM_Constants, rateBase )-swapSpread;
				double fwdPricing = CPISwapForward_WithoutSpread+coeffSpread;

				double pricingStrike = strike+ coeffSpread;
				double strikeForVol  = strike;
			
				/// a pay fixed swaption is a call on the inflation
				int callPut				= -GetOptionType();
				double swapTenor		= (GetEndDate().GetJulian()-GetStartDate().GetJulian())/K_YEAR_LEN;//numFromAsof-denomFromAsof;
				double swapExpiry		= (((ARM_SwapLeg*) itsSecondLeg)->GetResetDates()->Elt(0)-modelAsOfDate.GetJulian())/K_YEAR_LEN; //(((ARM_InfLeg*) itsFirstInfLeg)->GetFlowStartDates()->Elt(0)-modelAsOfDate.GetJulian())/K_YEAR_LEN;
				double optionMaturity	= (GetExpiryDate()-modelAsOfDate)/K_YEAR_LEN;
			
				InfSwopt_StoreVolInfo StoringFunc( this );
				/// the other information concerning the pricing is delegating to the storing func

				/// set the various element to the assets info
				itsAssetsInfo->SetOptionMaturity( optionMaturity );
				itsAssetsInfo->SetTenorAsset1( swapTenor );
				itsAssetsInfo->SetExpiryAsset1( swapExpiry );
				itsAssetsInfo->SetAnnuity(fixLegPVBP);
				itsAssetsInfo->SetPricingStrike( strike * CC_NS( ARM_Constants, rateBase ));
				itsAssetsInfo->SetFwdAsset1(CPISwapForward_WithoutSpread);				
				
				int notType = itsFirstInfLeg->GetAmount()->GetValueType();				
				if (notType==0)
				{
					
					itsFirstInfLeg->SetModel(PartialMod);

					ARM_ReferenceValue notional(*itsFirstInfLeg->GetAmount());
					notional.SetCalcFrequency(-1000);				

					ARM_SwapLeg* secondLeg				= (ARM_SwapLeg*) itsSecondLeg->Clone();
					secondLeg->GetIRIndex()->SetPayFrequency( itsFirstInfLeg->GetPaymentFreq());
					secondLeg->GetIRIndex()->SetResetFrequency(itsFirstInfLeg->GetPaymentFreq());
					secondLeg->CptCashFlowDates();

					secondLeg->SetAmount(&notional);					
					secondLeg->SetModelVariable(NULL);
					secondLeg->SetModel(PartialMod);

					ARM_GP_Vector floatStartTimes	= To_ARM_GP_Vector(itsFirstInfLeg->GetFlowStartDates());
					ARM_GP_Vector fixPayTimes		= To_ARM_GP_Vector(itsFirstInfLeg->GetPaymentDates());
					ARM_GP_Vector floatEndTimes		= To_ARM_GP_Vector(itsFirstInfLeg->GetFlowEndDates());
					ARM_GP_Vector fixPayPeriods		= To_ARM_GP_Vector(secondLeg->GetInterestTerms());
					ARM_GP_Vector* dfTerms		    = dynamic_cast< ARM_InfLeg* >(itsFirstInfLeg)->GetDiscountFactors();
				
					//ARM_GP_Vector fwdRates			= To_ARM_GP_Vector(itsFirstInfLeg->GetFwdRates());	
					
					//For the case Of VAr Notional we compute it here with firstInfLeg to have fwdRates without spread
					ARM_GP_Vector fwdRates			= To_ARM_GP_Vector(firstInfLeg->GetFwdRates());
					fwdRates/=100.0;

					ARM_GP_Vector CommunNotional(floatStartTimes.size(), 0.);
					for (int i=0; i<CommunNotional.size(); i++)
					{
						double paydate_i = fixPayTimes[i];
						CommunNotional[i] = itsFirstInfLeg->GetAmount()->CptReferenceValue(paydate_i);
					}
										
					int infsize = CommunNotional.size();
					ARM_GP_Vector infcoeffs(infsize);
					ARM_GP_Vector inftenors(infsize);
					ARM_GP_Vector infvols(infsize);
					double floatStartTime = floatStartTimes[0];

					double volAsset1 = 0.0;
					if(multiplier>1e-15)
					{
						fwdRates/=multiplier;
						strikeForVol/=multiplier;
						double basketInfVol  = PartialMod->ComputeInfVolWithVarNotional(
											CommunNotional,	
											floatStartTime,
											floatEndTimes,
											fixPayTimes,
											dfTerms,
											fixPayPeriods,
											fwdRates,
											strikeForVol,
											infidx,
											itsAssetsInfo,
											infcoeffs,
											inftenors,
											infvols);
						
						ARM_Swap infSwap( firstInfLeg, secondLeg );
						infSwap.SetModel( PartialMod );
						double CPISwapForwardForVol= fabs(infSwap.PriceToRate( GetStartDate(), targetPrice )/CC_NS( ARM_Constants, rateBase ))/multiplier;
						/////////////////

						volAsset1 = basketInfVol/(CPISwapForwardForVol+1.0);
						
						int callput		= (fwdPricing>pricingStrike) ? -1 : 1;
						double normalPrice = VanillaOption_N(CPISwapForwardForVol+1.0, 
														basketInfVol,
														strikeForVol+1.0, 
														optionMaturity,  
														callput );


						success = true;
						volAsset1 = VanillaImpliedVol_BS(CPISwapForwardForVol+1.0, 
															strikeForVol+1.0, 
															optionMaturity, 
															normalPrice, 
															callput,
															&volAsset1,
															&success);

						
					}
					else
					{
						volAsset1 = K_NEW_DOUBLE_TOL;
					}

					itsAssetsInfo->SetVolAsset1(volAsset1);					
					
					delete secondLeg;					
				}

				itsAssetsInfo->SetStrikeForVol1(strikeForVol* CC_NS( ARM_Constants, rateBase ));
				
				/// computes the price
				price	= pricingMod->SingleAssetOptionPrice( fwdPricing, pricingStrike, callPut, fixLegPVBP,  infidx, itsAssetsInfo, StoringFunc );
				delete firstInfLeg;
				delete fixLeg;
				delete swapSpread_Null;
				delete infLegSpread1;

			}
			break;
			
			/// 2) case of libor vs inflation
		case K_FLOATING_LEG:
			{
				/// set the currency according to the fixed leg
				if( strcmp( itsSecondLeg->GetCurrencyUnit()->GetCcyName(), itsFirstInfLeg->GetCurrencyUnit()->GetCcyName() ) )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						"fix vs inflation only allowed for inflation leg with same currency as the fixed leg!");

				/// Call or Put?
				int callPut				= GetOptionType();

				/// create the corresponding fixed legs
				ARM_InfIdx*  infidx		= (ARM_InfIdx*) itsFirstInfLeg->GetIRIndex();
				ARM_IRIndex* irIdx		= itsSecondLeg->GetIRIndex();
				double fixedRate		= 1.0;
				ARM_Currency* infLegCcy	= itsFirstInfLeg->GetCurrencyUnit();
				ARM_Currency* libLegCcy = itsSecondLeg->GetCurrencyUnit();

				int infFreq = itsFirstInfLeg->GetPaymentFreq();
				int secondLegFreq = itsSecondLeg->GetPaymentFreq();

				int fixlegFreq = (infFreq <secondLegFreq) ? infFreq : secondLegFreq;
				
				/// create the fixed leg of the inflation rate
				ARM_SwapLeg* InfFixedLeg= new ARM_SwapLeg( 
					itsFirstInfLeg->GetStartDateNA(), itsFirstInfLeg->GetEndDateNA(), fixedRate, -itsFirstInfLeg->GetRcvOrPay(),
					fixlegFreq, infLegCcy->GetFixedDayCount(), K_COMP_PROP, K_ARREARS,
					K_ADJUSTED, K_SHORTSTART, infLegCcy, infLegCcy->GetCcyName(), K_NX_NONE, "NULL", K_YES );
				ARM_InfLeg* firstInfLeg2 = dynamic_cast< ARM_InfLeg* >(itsFirstInfLeg->Clone());
				double multiplier = firstInfLeg2->GetMultiple();
				
				ARM_SwapLeg* floatMatchingFixedLeg = new ARM_SwapLeg( 
					itsFirstInfLeg->GetStartDateNA(), itsFirstInfLeg->GetEndDateNA(), fixedRate, itsFirstInfLeg->GetRcvOrPay(),
					fixlegFreq, infLegCcy->GetFixedDayCount(), K_COMP_PROP, K_ARREARS,
					K_ADJUSTED, K_SHORTSTART, infLegCcy, infLegCcy->GetCcyName(), K_NX_NONE, "NULL", K_YES );
							
				/// create the fixed leg of the libor rate
				ARM_SwapLeg* LibFixedLeg= new ARM_SwapLeg
					(itsSecondLeg->GetStartDateNA(), itsSecondLeg->GetEndDateNA(), fixedRate, itsFirstInfLeg->GetRcvOrPay(),
					fixlegFreq, libLegCcy->GetFixedDayCount(), K_COMP_PROP, K_ARREARS,
					K_ADJUSTED, K_SHORTSTART, libLegCcy, libLegCcy->GetCcyName(),K_NX_NONE, "NULL", K_YES );
							
				/// Recover the inflation spread part
				
				
				if( firstInfLeg2->GetSwapType() != K_YEARTOYEAR_LEG && firstInfLeg2->GetSwapType() != K_OATTYPE_LEG)
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						ARM_USERNAME + ": only year on year or an OAT leg supported!");

				double spreadInf = 0., spreadInfEq = 0., InfFixLegPVBP = 0., floatSpreadPV_INF	= 0.0, infWithSpread = 0.0, infWithOutSpread=0.0 ;
				
				/// Recover the inflation spread part in the case of the year on year swaption
				if( firstInfLeg2->GetSwapType() == K_YEARTOYEAR_LEG )
				{
					spreadInf	= firstInfLeg2->GetConstant() + firstInfLeg2->GetMultiple();
					spreadInfEq = spreadInf * CC_NS( ARM_Constants, rateBase )*itsFirstInfLeg->GetRcvOrPay();			

					double newConstant	= firstInfLeg2->GetConstant()-spreadInf;
					firstInfLeg2->PropagateModel( model );
					firstInfLeg2->SetModel( model );
					infWithSpread = firstInfLeg2->ComputePrice();

					firstInfLeg2->SetConstant( newConstant );
					firstInfLeg2->PropagateModel( model );
					firstInfLeg2->SetModel( model );
					infWithOutSpread = firstInfLeg2->ComputePrice();

					floatSpreadPV_INF = infWithSpread-infWithOutSpread;

					InfFixedLeg->SetModel( model);
					InfFixLegPVBP = InfFixedLeg->ComputePrice() / InfFixedLeg->GetFixedRate()*InfFixedLeg->GetRcvOrPay();		
										
				}
				
				
				/// Recover the Libor spread part
				double spreadLibor = itsSecondLeg->GetSpread();
				itsSecondLeg->SetModel( model );
				double priceFloatLegWithSpread		= itsSecondLeg->ComputePrice();
				ARM_SwapLeg* floatLeg				= (ARM_SwapLeg*) itsSecondLeg->Clone();
				floatLeg->SetSpread(0.0);
				floatLeg->SetModel( model );
				double priceFloatLegWithOutSpread	= floatLeg->ComputePrice();

				double floatSpreadPV_IR = priceFloatLegWithSpread-priceFloatLegWithOutSpread;
		

				ARM_SwapLeg* floatLeg2				= (ARM_SwapLeg*) itsSecondLeg->Clone();
				floatLeg2->SetSpread(0.0);
				floatLeg2->SetModel( model );
				double priceFloatLegWithOutSpread2	= floatLeg2->ComputePrice();
				ARM_ReferenceValue notionalfloat(*floatLeg2->GetAmount());				
				//
				notionalfloat.SetCalcFrequency(-1000);
				LibFixedLeg->SetAmount(&notionalfloat);

				LibFixedLeg->SetModelVariable(NULL);
///YK
				LibFixedLeg->SetModel( model);
				double LibFixLegPVBP = LibFixedLeg->ComputePrice() / LibFixedLeg->GetFixedRate()*LibFixedLeg->GetRcvOrPay();

				floatMatchingFixedLeg->SetAmount(&notionalfloat);
				floatMatchingFixedLeg->SetModelVariable(NULL);
				floatMatchingFixedLeg->SetModel(model);
				

				/// Recover the fwd nominal swap rate
				double targetPrice		= 0.0;
				ARM_Swap* floatSwap		= new ARM_Swap( floatLeg2, floatMatchingFixedLeg );
				floatSwap->SetModel( model );
				double floatSwapRate	= fabs(floatSwap->PriceToRate( GetStartDate(), targetPrice )/ CC_NS( ARM_Constants, rateBase ));

			
				InfFixedLeg->SetModel( model );
				ARM_ReferenceValue notional(*itsFirstInfLeg->GetAmount());
				notional.SetCalcFrequency(-1000);
				InfFixedLeg->SetAmount(&notional);
				
				InfFixedLeg->SetModelVariable(NULL);
				InfFixedLeg->SetModel(model);

				
				ARM_Swap* infSwap		= new ARM_Swap( firstInfLeg2, InfFixedLeg );
				infSwap->SetModel( model );
				double CPISwapForward	= infSwap->PriceToRate( GetStartDate(), targetPrice )/CC_NS( ARM_Constants, rateBase );

				/// to compute the equivalent strike, we regroup the inflation and the Libor 
				/// spread part over over the fixLegPVBP!

				/// Renormalize the inflation spread to the frequence and the basis of the Libor one
				spreadLibor	= (priceFloatLegWithSpread - priceFloatLegWithOutSpread) / LibFixLegPVBP;
				
				double strike = GetStrike() + (floatSpreadPV_IR*itsSecondLeg->GetRcvOrPay()- floatSpreadPV_INF*itsFirstInfLeg->GetRcvOrPay())/LibFixLegPVBP;

				strike = strike/CC_NS( ARM_Constants, rateBase );
				
		
			
				ARM_InfBSModel* PartialMod = dynamic_cast<ARM_InfBSModel*>(pricingMod);
				////////////////////////////////
				ARM_InfLeg* infLegspread1 = dynamic_cast< ARM_InfLeg* >(firstInfLeg2->Clone());
				double spread1 = 0.01;
				double newConstant_spread1 = infLegspread1->GetConstant() + spread1; 
				infLegspread1->SetConstant( newConstant_spread1 );
				infLegspread1->PropagateModel( PartialMod );
				infLegspread1->SetModel( PartialMod );
				double infWithSpread1 = infLegspread1->ComputePrice();
				double floatLegPVBP = infWithSpread1-infWithOutSpread;
				double coeffSpread = multiplier*fabs(floatLegPVBP/LibFixLegPVBP);
				double fwdPricing = CPISwapForward+coeffSpread;
				double pricingStrike = strike+coeffSpread;
				double strikeForVolInf  = (floatSwapRate+strike);
				double strikeForVolIR  = CPISwapForward-strike;
				itsAssetsInfo->SetStrikeForVol1(strikeForVolInf);
				itsAssetsInfo->SetStrikeForVol2(strikeForVolIR);				
				/////////////////////////////////////

				double optionMaturity	= (GetExpiryDate()-modelAsOfDate)/K_YEAR_LEN;
				double swapTenor		= (GetEndDate().GetJulian()-GetStartDate().GetJulian())/K_YEAR_LEN;
				double tenorAsset1		= (itsFirstInfLeg->GetEndDate().GetJulian()-itsFirstInfLeg->GetStartDate().GetJulian())/K_YEAR_LEN;
				double tenorAsset2		= (itsSecondLeg->GetEndDate().GetJulian()-itsSecondLeg->GetStartDate().GetJulian())/K_YEAR_LEN;
				double swapExpiry		= (((ARM_InfLeg*) itsFirstInfLeg)->GetFlowStartDates()->Elt(0)-modelAsOfDate.GetJulian())/K_YEAR_LEN;
				InfSwopt_StoreVolInfo StoringFunc( this );

				/// set the various element to the assets info
				itsAssetsInfo->SetOptionMaturity( optionMaturity );
				itsAssetsInfo->SetAnnuity(LibFixLegPVBP);
				itsAssetsInfo->SetPricingStrike( strike * CC_NS( ARM_Constants, rateBase ));
				itsAssetsInfo->SetTenorAsset1( tenorAsset1 );
				itsAssetsInfo->SetExpiryAsset1( swapExpiry);
				itsAssetsInfo->SetExpiryAsset2( optionMaturity );
				itsAssetsInfo->SetTenorAsset2( tenorAsset2 );
				itsAssetsInfo->SetFirstAssetType( firstInfLeg2->GetSwapType() );
				itsAssetsInfo->SetFirstAssetCoupon( firstInfLeg2->GetMultiple() );
				itsAssetsInfo->SetFwdAsset1(CPISwapForward*CC_NS( ARM_Constants, rateBase ));
				itsAssetsInfo->SetFwdAsset2(floatSwapRate*CC_NS( ARM_Constants, rateBase ));
				
			
				
				//if (itsFirstInfLeg->GetAmortization())
				int notType = itsFirstInfLeg->GetAmount()->GetValueType();
				if (notType==0)
				{
					if( !(PartialMod->GetCorrelMatrixValidated()))
					{ 
						PartialMod->ValidateCorrelManagerForVarNotionalSwaption(infidx,irIdx);
					}

					if(!(PartialMod->GetisGoodCorrelMatrix()))
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						"Correltion Matrix for the univers INF/IR is not of good type eigen values are negatives");
					
					floatMatchingFixedLeg->SetAmount(&notionalfloat);
					floatMatchingFixedLeg->SetModelVariable(NULL);
					floatMatchingFixedLeg->SetModel(model);
				
				
					ARM_GP_Vector infFloatStartTimes	= To_ARM_GP_Vector(itsFirstInfLeg->GetFlowStartDates());
					ARM_GP_Vector infFixPayTimes		= To_ARM_GP_Vector(itsFirstInfLeg->GetPaymentDates());
					ARM_GP_Vector infFloatEndTimes		= To_ARM_GP_Vector(itsFirstInfLeg->GetFlowEndDates());
					ARM_GP_Vector* infDfTerms			= firstInfLeg2->GetDiscountFactors();
					ARM_GP_Vector infFwdRates			= To_ARM_GP_Vector(firstInfLeg2->GetFwdRates());
					
					infFwdRates/= 100.0;
					delete InfFixedLeg;
					InfFixedLeg= new ARM_SwapLeg( 
					itsFirstInfLeg->GetStartDateNA(), itsFirstInfLeg->GetEndDateNA(), fixedRate, -itsFirstInfLeg->GetRcvOrPay(),
					firstInfLeg2->GetPaymentFreq(), infLegCcy->GetFixedDayCount(), K_COMP_PROP, K_ARREARS,
					K_ADJUSTED, K_SHORTSTART, infLegCcy, infLegCcy->GetCcyName(), K_NX_NONE, "NULL", K_YES );

					ARM_GP_Vector infFixPayPeriods		= To_ARM_GP_Vector(InfFixedLeg->GetInterestTerms());

					////////////////////
					InfFixedLeg->SetAmount(&notional);
					InfFixedLeg->SetModelVariable(NULL);
					InfFixedLeg->SetModel(model);

					delete infSwap;
					infSwap		= new ARM_Swap( firstInfLeg2, InfFixedLeg );
					infSwap->SetModel( model );
					double CPISwapForwardForVol= 1.0;
					/////////////////

					ARM_GP_Vector infNotional(infFloatStartTimes.size(), 0.);
					for (int i=0; i<infNotional.size(); i++)
					{
						double paydate_i = infFixPayTimes[i];
						infNotional[i] = itsFirstInfLeg->GetAmount()->CptReferenceValue(paydate_i);
					}
					
					int infsize = infNotional.size();
					ARM_GP_Vector infcoeffs(infsize);
					ARM_GP_Vector inftenors(infsize);
					ARM_GP_Vector infvols(infsize);
					double infFloatStartTime = infFloatStartTimes[0];

					double basketInfVol = K_NEW_DOUBLE_TOL;

					if(multiplier>1e-15)
					{
						infFwdRates/=multiplier;
						strikeForVolInf  /= multiplier;			

						basketInfVol = PartialMod->ComputeInfVolWithVarNotional(
											infNotional,	
											infFloatStartTime,
											infFloatEndTimes,
											infFixPayTimes,
											infDfTerms,
											infFixPayPeriods,
											infFwdRates,
											strikeForVolInf,
											infidx,
											itsAssetsInfo,
											infcoeffs,
											inftenors,
											infvols);

						CPISwapForwardForVol= infSwap->PriceToRate( GetStartDate(), targetPrice )/CC_NS( ARM_Constants, rateBase );
					}
					
					itsAssetsInfo->SetStrikeForVol1(strikeForVolInf);
									
				
					ARM_GP_Vector irFloatStartTimes	= To_ARM_GP_Vector(floatLeg->GetFlowStartDates());
					ARM_GP_Vector irFixPayTimes		= To_ARM_GP_Vector(floatLeg->GetPaymentDates());
					ARM_GP_Vector irFloatEndTimes	= To_ARM_GP_Vector(floatLeg->GetFlowEndDates());

					delete floatMatchingFixedLeg;
					floatMatchingFixedLeg = new ARM_SwapLeg( 
					itsFirstInfLeg->GetStartDateNA(), itsFirstInfLeg->GetEndDateNA(), fixedRate, itsFirstInfLeg->GetRcvOrPay(),
					secondLegFreq, infLegCcy->GetFixedDayCount(), K_COMP_PROP, K_ARREARS,
					K_ADJUSTED, K_SHORTSTART, infLegCcy, infLegCcy->GetCcyName(), K_NX_NONE, "NULL", K_YES );
				
					ARM_GP_Vector irFixPayPeriods	= To_ARM_GP_Vector(floatMatchingFixedLeg->GetInterestTerms());

					//////////////////////////
					floatMatchingFixedLeg->SetAmount(&notionalfloat);
					floatMatchingFixedLeg->SetModelVariable(NULL);
					floatMatchingFixedLeg->SetModel(model);
				

					/// Recover the fwd nominal swap rate
					delete floatSwap;
					floatSwap		= new ARM_Swap( floatLeg2, floatMatchingFixedLeg );
					floatSwap->SetModel( model );
					double floatSwapRateForVol	= fabs(floatSwap->PriceToRate( GetStartDate(), targetPrice )/ CC_NS( ARM_Constants, rateBase ));			

					///////////////////////////
									
					double startDate  = itsSecondLeg->GetStartDateNA().DMYToJulian();
					double endDate	  = itsSecondLeg->GetEndDateNA().DMYToJulian();

					ARM_GP_Vector irNotional(irFloatStartTimes.size(), 0.);
					for (int j=0; j<irNotional.size(); j++)
					{
						double paydate_j = irFixPayTimes[j];
						irNotional[j] = itsSecondLeg->GetAmount()->CptReferenceValue(paydate_j);
					}

					int irsize = irNotional.size();
					ARM_GP_Vector ircoeffs(irsize);
					ARM_GP_Vector irtenors(irsize);
					ARM_GP_Vector irvols(irsize);
					
					double strikeIr = -(floatSwapRate-CPISwapForward+strike);
					double basketIrVol = PartialMod->ComputeIRVolWithVarNotional(
										irNotional, 
										startDate,
										endDate,
										irFloatStartTimes,
										irFloatEndTimes,
										irFixPayTimes,
										irFixPayPeriods,
										strikeIr,
										irIdx,
										itsAssetsInfo,
										ircoeffs,
										irtenors,
										irvols);
					 // Stores Ir Vol


					int callput		= (strikeIr<0.0) ? -1 : 1;
					double strikeOption = floatSwapRateForVol+strikeIr;
					if (strikeOption < 1e-3) strikeOption = 1e-3;

					double normalPrice = VanillaOption_N(floatSwapRateForVol, 
													basketIrVol,
													strikeOption, 
													optionMaturity,  
													callput );


					double volAsset2		= basketIrVol/floatSwapRateForVol;

					success = true;
					volAsset2 = VanillaImpliedVol_BS(floatSwapRateForVol, 
														strikeOption, 
														optionMaturity, 
														normalPrice, 
														callput,
														&volAsset2,
														&success);

					double volAsset1 = basketInfVol/(CPISwapForwardForVol+1.0);
						
					int callputInf		= (CPISwapForwardForVol>strikeForVolInf) ? -1 : 1;
					double normalPriceInf = VanillaOption_N(CPISwapForwardForVol+1.0, 
													basketInfVol,
													strikeForVolInf+1.0, 
													optionMaturity,  
													callputInf );


					success = true;
					volAsset1 = VanillaImpliedVol_BS(CPISwapForwardForVol+1.0, 
														strikeForVolInf+1.0, 
														optionMaturity, 
														normalPriceInf, 
														callputInf,
														&volAsset1,
														&success);
										
					itsAssetsInfo->SetVolAsset1(volAsset1);
					itsAssetsInfo->SetVolAsset2(volAsset2);


					double covar = PartialMod->ComputeInfIrCovarWithVarNotional(infcoeffs,inftenors,infvols,
						ircoeffs,irtenors,irvols,infidx,irIdx);
					
					double correl = 0.0; 

					if( (basketIrVol>1e-15) && (basketInfVol>1e-15) )
						correl = covar/(basketIrVol*basketInfVol);

					itsAssetsInfo->SetCorrelation(correl);					
				}

				/*price	= pricingMod->TwoAssetsOptionPrice( CPISwapForward, floatSwapRate,
					strike, callPut, LibFixLegPVBP, infidx, irIdx, itsAssetsInfo, StoringFunc );*/

				price	= pricingMod->TwoAssetsOptionPrice( fwdPricing, floatSwapRate,
					pricingStrike, callPut, LibFixLegPVBP, infidx, irIdx, itsAssetsInfo, StoringFunc );


				
				/// restore spread

				if( price < 0.0 )
					price = 0.0;
				
				/// the other information concerning the pricing is delegating to the storing func
				delete InfFixedLeg;
				delete LibFixedLeg;
				delete floatMatchingFixedLeg;
				delete infSwap;
				delete floatSwap;
				delete floatLeg;
				delete floatLeg2;
				delete firstInfLeg2;
				delete infLegspread1;
			}
			break;
			
			/// 3) case of inflation vs inflation
		case K_GENERICINF_LEG:
			{
				double fixedRate					= 1.0;
				ARM_Currency* firstInfLegCcy		= itsFirstInfLeg->GetCurrencyUnit();
				ARM_SwapLeg* firstInfFixedLeg		= new ARM_SwapLeg( 
					itsFirstInfLeg->GetStartDateNA(), itsFirstInfLeg->GetEndDateNA(), fixedRate, -itsFirstInfLeg->GetRcvOrPay(),
					firstInfLegCcy->GetFixedPayFreq(), firstInfLegCcy->GetFixedDayCount(), K_COMP_PROP, K_ARREARS,
					K_ADJUSTED, K_SHORTSTART, firstInfLegCcy, firstInfLegCcy->GetCcyName(), 
					K_NX_NONE, "NULL", K_YES );
				
				ARM_SwapLeg* secondInfMirrorFixedLeg = new ARM_SwapLeg( *firstInfFixedLeg );
				secondInfMirrorFixedLeg->SetRcvOrPay( -firstInfFixedLeg->GetRcvOrPay() );
				
				ARM_Swap* firstSwap		= new ARM_Swap( itsFirstInfLeg, firstInfFixedLeg );
				firstSwap->SetModel( model );
				double targetPrice		= 0.0;
				double CPISwapForward1	= firstSwap->PriceToRate( GetStartDate(), targetPrice )/ CC_NS( ARM_Constants, rateBase );

				ARM_Swap* secondSwap	= new ARM_Swap( itsSecondLeg, secondInfMirrorFixedLeg  );
				secondSwap->SetModel( model );
				double CPISwapForward2	= secondSwap->PriceToRate( GetStartDate(), targetPrice )/ CC_NS( ARM_Constants, rateBase );
				int callPut				= GetOptionType();

				double optionMaturity	= (GetExpiryDate()-modelAsOfDate)/K_YEAR_LEN;
				double swapTenor		= (GetEndDate().GetJulian()-GetStartDate().GetJulian())/K_YEAR_LEN;
				double tenorAsset1		= (itsFirstInfLeg->GetEndDate().GetJulian()-itsFirstInfLeg->GetStartDate().GetJulian())/K_YEAR_LEN;
				double tenorAsset2		= (itsSecondLeg->GetEndDate().GetJulian()-itsSecondLeg->GetStartDate().GetJulian())/K_YEAR_LEN;
				double swapExpiry1		= (((ARM_InfLeg*) itsFirstInfLeg)->GetFlowStartDates()->Elt(0)-modelAsOfDate.GetJulian())/K_YEAR_LEN;
				double swapExpiry2		= (((ARM_InfLeg*) itsSecondLeg)->GetFlowStartDates()->Elt(0)-modelAsOfDate.GetJulian())/K_YEAR_LEN;
				firstInfFixedLeg->SetModel( model );
				double fixLegPVBP		= firstInfFixedLeg->ComputePrice() / firstInfFixedLeg->GetFixedRate() * firstInfFixedLeg->GetRcvOrPay();
				double strike			= GetStrike() / CC_NS( ARM_Constants, rateBase );
				InfSwopt_StoreVolInfo StoringFunc( this );

				/// set the various element to the assets info
				itsAssetsInfo->SetOptionMaturity( optionMaturity );
				itsAssetsInfo->SetAnnuity(fixLegPVBP);

				itsAssetsInfo->SetPricingStrike( GetStrike() );
				itsAssetsInfo->SetTenorAsset1( tenorAsset1 );
				itsAssetsInfo->SetTenorAsset2( tenorAsset2 );
				itsAssetsInfo->SetExpiryAsset1( swapExpiry1);
				itsAssetsInfo->SetExpiryAsset2( swapExpiry2);

				ARM_InfIdx* infidx		= (ARM_InfIdx*) itsFirstInfLeg->GetIRIndex();
				price	= pricingMod->TwoAssetsOptionPrice( CPISwapForward1, CPISwapForward2,
					strike, callPut, fixLegPVBP, infidx, itsSecondLeg->GetIRIndex(), itsAssetsInfo, StoringFunc );
				if( price < 0.0 )
					price = 0.0;

				delete firstInfFixedLeg;
				delete secondInfMirrorFixedLeg;
				delete firstSwap;
				delete secondSwap;
			}
			break;
			/// other cases!
		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": swaption can be inflation vs fixed, inflation vs float, inflation vs inflation!");
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
void ARM_InfSwaption::View(char* id, FILE* ficOut )
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
    fprintf(fOut, " =======> INFLATION SWAPTION  <====== \n\n" );

	if( itsAssetsInfo )
	{
		fprintf( fOut, " Quick View\n\n");
		fprintf( fOut, "%s\n", itsAssetsInfo->toString().c_str() );
	}

    fprintf(fOut, "\n\n =======> UNDERLYING INFLATION SWAP <====== \n\n" );

	ARM_Swaption::View(id,fOut);
	
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
void ARM_InfSwaption::CptCashFlowValues()
{}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/