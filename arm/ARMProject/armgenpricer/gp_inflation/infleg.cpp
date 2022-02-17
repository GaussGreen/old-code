/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file infleg.cpp
 *
 *  \brief implements the ARM_InfLeg class, a class to deal with SwapLegs
 *		Inherit from ARM_SwapLeg
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

#include "gpinflation/infleg.h"

/// gpinflation
#include "gpinflation/infidx.h"			/// for inflation Index
#include "gpinflation/infdata.h"		/// for default data
#include "gpinflation/infcurvmodel.h"	/// to get access to the inflation curve model

/// gpbase
#include "gpbase/datestrip.h"			/// for ARM_DateStrip
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gpvector.h"			/// for ARM_GP_Vector
#include "gpbase/gplinalgconvert.h"		/// for conversion

/// kernel
#include <glob/paramview.h>				/// to access the paramview method

CC_BEGIN_NAMESPACE( ARM )


/// function to convert date into gap if it is not really a gap but rather a 
///	date! uses a calendar
int ARM_InfLeg::GetGapFromGapOrDate( int rollDate, const char* calendar, const ARM_Date& refDate )
{
	/// test whether it is really a julian date!
	const int JULIANDATEMIN = 2000000.0;
	if(rollDate<JULIANDATEMIN )
		return rollDate;
	ARM_Date targetGapDate( (double)rollDate );	

	/// computes naive gap
	int gap=targetGapDate-const_cast<ARM_Date&>(refDate);
	/// if gap>0... gap is potentially 
	int sign,period;
	if(gap>0)
	{
		sign	= 1;
		period  = K_DAILY;
	}
	else
	{
		sign	=-1;
		period  =-K_DAILY;
		gap		=-gap;
	}


	ARM_Date gapDate = refDate;
	gapDate.AddPeriodMult(period,gap,const_cast<char*>(calendar));

	while((targetGapDate-gapDate)*sign<0)
	{
		--gap;
		gapDate=refDate;
		gapDate.AddPeriodMult(period,gap,const_cast<char*>(calendar));
	}
	return gap*sign;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: Constructor
///	Returns: 
///	Action : General constructor with 2 dateStrips
///				ARM_DateStrip gives a very condensed form of the constructor
///				this should be often preferred
///
///			Because of the non constcorrectness of ARM
///			we are forced to const_cast... not for pleasure
///			but for legacy reasons
////////////////////////////////////////////////////

ARM_InfLeg::ARM_InfLeg(
        const string& indexName, 
		int rcvOrPay,
		int interpType,
		double multiple,
        double constant,
		int finalNotionalType,
		const ARM_DateStrip* numDateStrip,
		const ARM_DateStrip* denomDateStrip,
		const ARM_Currency*	discountCcy )
:
	ARM_SwapLeg(),
	/// member data initialisation
	itsSwapType( K_GENERICINF_LEG ),
	itsInterpType( interpType	),
	itsMultiple( multiple		),
	itsConstant( constant		),
	itsNumResetDates( NULL		),
	itsDenomResetDates( NULL	),
	itsNumCPIRates( NULL		),
	itsDenomCPIRates( NULL		),
	itsDiscountFactors( NULL	),
	itsNumPublicationDates( NULL),
	itsFirstReset( GETDEFAULTVALUE )
{
	///  validate the datestrips
	/// checks that num and denom datstrip have the same payment dates
	if( *(numDateStrip->GetPaymentDates()) != 
		*(denomDateStrip->GetPaymentDates()) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		   "denom strip dates and reset strip dates do not have the same payment dates");

	/// creates a dummy index
	ARM_InfIdx* tmpInfIdx = new ARM_InfIdx( indexName );
	Set( tmpInfIdx, rcvOrPay, finalNotionalType, discountCcy, numDateStrip, denomDateStrip );
	delete tmpInfIdx;

	/// Initialisation
	SetName( ARM_INFLEG );
}



////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: Set
///	Returns: 
///	Action : Set function. Initialises 
///
///			- the reset numerator dates
///			- the reset denominator dates
///			- the underlying swap
///
////////////////////////////////////////////////////
void ARM_InfLeg::Set(
        const ARM_InfIdx* infIdx, 
		int rcvOrPay,
		int finalNotionalType,
		const ARM_Currency*	discountCcy,
		const ARM_DateStrip* numDateStrip,
		const ARM_DateStrip* denomDateStrip )
{
	/// validation: check that we have same
	/// nb of num and denom resetDates
	if( numDateStrip->GetResetDates()->size() !=
		denomDateStrip->GetResetDates()->size() )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			   "denom strip dates and reset strip dates are not of the same size!");

	if( finalNotionalType != K_NX_NONE	 && 
		finalNotionalType != K_NX_INFEND && 
		finalNotionalType != K_NX_END	 &&
 		finalNotionalType != K_NX_ASINFEND )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		   "invalid notional type");

	SetInitRcvOrPay( rcvOrPay );

	/// the default currency is either the one of the inflation index
	/// or if not provided the default currency
	if( discountCcy == NULL )
	{
		if( GetIRIndex() )
			SetCurrencyUnit( (ARM_Currency*) GetIRIndex()->GetCurrencyUnit()->Clone() );
		else
			SetCurrencyUnit( ARM_DEFAULT_CURRENCY );
	}
	else
		SetCurrencyUnit( const_cast<ARM_Currency*>(discountCcy)  );

	setIndexName( infIdx->GetIndexName() );
	setPayCurrencyName( GetCurrencyUnit()->GetCcyName() );
	
	/// copies the reset dates from the denominator strip dates
	itsDenomResetDates = new ARM_GP_Vector( *(denomDateStrip->GetResetDates() ) );

	/// if the numerator strip date is not equal to zero
	/// can initialise properly the underlying swap
	itsNumResetDates = new ARM_GP_Vector( *(numDateStrip->GetResetDates() ) );

	/// Sets the underlying swap
	/// Computes the various data to initialise the underlying swap leg
	int size = itsNumResetDates->size();

	ARM_GP_Vector* tmpForwards			= new ARM_GP_Vector( size, 0.0 );
	ARM_GP_Vector* tmpCashFlowValues	= new ARM_GP_Vector( size, 0.0 );
	ARM_GP_Vector* tmpNotional			= new ARM_GP_Vector( size, 100.0 );

	if( itsSwapType == K_YEARTOYEAR_LEG )
	{
		for (size_t j=0; j<size; j++)
		{
			//tmpNotional->Elt(j) = GetAmount()->GetDiscreteValues()->Elt(j);
		}
	}
	/// similarly for the spread and fixRate arguments
	double spread	= 0.0;
	double fixRate	= K_FLOAT_RATE;

	/// to get exactly the same interest days and base
	ARM_GP_Vector* InterestDays = numDateStrip->GetInterestTerms();
	double basis = FromDayCountToBasis(infIdx->GetDayCount());
	for( int i=0; i<InterestDays->size(); ++i)
		(*InterestDays)[i] *= basis;

	CC_NS(std,auto_ptr)<ARM_ReferenceValue> pAmount(new ARM_ReferenceValue(
		To_pARM_Vector(numDateStrip->GetPaymentDates()),
		To_pARM_Vector(tmpNotional)));
	delete tmpNotional;
	
	/// and finally use the public function Set (Summit style) of SwapLeg
	/// no need to clone the underlying vector since
	/// set clones for us
	ARM_Vector* tmpForwards2 = To_pARM_Vector(tmpForwards);
	delete tmpForwards;

	ARM_Vector* FlowStartDates		= To_pARM_Vector(numDateStrip->GetFlowStartDates());
	ARM_Vector* FlowEndDates		= To_pARM_Vector(numDateStrip->GetFlowEndDates());
	ARM_Vector* FlowPaymentDates	= To_pARM_Vector(numDateStrip->GetPaymentDates());
	ARM_Vector* FlowResetDates		= To_pARM_Vector(numDateStrip->GetResetDates());
	ARM_Vector* InterDays			= To_pARM_Vector(InterestDays);

	ARM_SwapLeg::Set(
		FlowStartDates,
		FlowEndDates,
		FlowPaymentDates,
		FlowResetDates,
		InterDays,
		tmpForwards2,
		pAmount.get(),
		const_cast<ARM_InfIdx*>( infIdx ),
		rcvOrPay,
		spread,
		fixRate,
		GetCurrencyUnit() );
	
	if( FlowStartDates	)	{ delete FlowStartDates;	FlowStartDates	=NULL; }
	if( FlowEndDates	)	{ delete FlowEndDates;		FlowEndDates	=NULL; }
	if( FlowPaymentDates)	{ delete FlowPaymentDates;	FlowPaymentDates=NULL; }
	if( FlowResetDates	)	{ delete FlowResetDates;	FlowResetDates	=NULL; }
	if( InterDays		)	{ delete InterDays;			InterDays		=NULL; }
	/// set the notional and we do not use the swaplet set with notional
	/// because of the implied computation!
	SetNxFlag( finalNotionalType );

	itsNumCPIRates		= new ARM_GP_Vector( size, 0.0 );
	itsDenomCPIRates	= new ARM_GP_Vector( size, 0.0 );
	itsDiscountFactors	= new ARM_GP_Vector( size, 0.0 );

	/// no daycount for zero coupon leg
	if( itsSwapType == K_ZEROCOUPON_LEG )
	{
		/// Interest Terms should be always 1
		/// because it clones it use an object!
		ARM_Vector tmpInterestTerms(1, 1);
		SetInterestTerms( &tmpInterestTerms );
	}

	SetFwdRates( tmpForwards2 );
	delete tmpForwards2;

	ARM_Vector* tmpCashFlowValues2 = To_pARM_Vector(tmpCashFlowValues);
	delete tmpCashFlowValues;
	SetCashFlowValues(tmpCashFlowValues2	);	/// do not clone it so do not delete it!

	if( itsSwapType == K_YEARTOYEAR_LEG )
		UpdateLegAmortization();
}



////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: ARM_InfLeg (constructor)
///	Returns: 
///	Action : Version to specify many of the inputs
///				without using a dateStrip object
////////////////////////////////////////////////////

ARM_InfLeg::ARM_InfLeg(
	const ARM_Date& startDate,
	const ARM_Date& endDate,
	const string& indexName,
	int swapType,
	int rcvOrPay,
	int interpType,
	double multiple, 
	double Comultiple, 
    double constant, 
	int resetFreq,
	int dayCount,
	const char* resetCalendar,
	int resetRule,
	int intRule,
	int stubRule,
	int resetNumGaporRollDate,
	int resetDenomGaporRollDate,
	int payFreq,
	int payGap,
	const char* payCalendar,
	int adjFirstDate,
	int finalNotionalType,
	double firstReset,
	const ARM_Currency*	discountCcy )
:
	/// member data initialisation
	itsFirstReset( firstReset ),
	itsInterpType( interpType ),
	itsSwapType( swapType	),
	itsNumResetDates( NULL	),
	itsDenomResetDates( NULL),
	itsNumCPIRates( NULL	),
	itsDenomCPIRates( NULL	),
	itsNumPublicationDates( NULL ),
	itsResetNumGaporRollDate( resetNumGaporRollDate ),
	itsResetDenomGaporRollDate( resetDenomGaporRollDate )

{
	const int JULIANDATEMIN = 2000000.0;
	/// creates a dummy index
	ARM_InfIdx* tmpInfIdx = new ARM_InfIdx( indexName );

	char* ourPayCalendar;
	if( strcmp( payCalendar, GETDEFAULTVALUESTR ) == 0 )
		ourPayCalendar = tmpInfIdx->GetCurrencyUnit()->GetCcyName();
	else
		ourPayCalendar = const_cast<char*>(payCalendar);

	switch( itsSwapType )
	{
	case K_YEARTOYEAR_LEG:
		{
			if( resetFreq == K_DEF_FREQ )
				payFreq = resetFreq = K_ANNUAL;
			break;
		}
	case K_ZEROCOUPON_LEG:
		{
			if( resetFreq == K_DEF_FREQ )
				payFreq = resetFreq = K_ZEROCOUPON;
			break;
		}
	case K_OATTYPE_LEG:
		{
			if( resetFreq == K_DEF_FREQ )
				payFreq = resetFreq = K_ANNUAL;
			break;
			
		}
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"unsupported swap type");
		}
	}

	/// gets the default for payFreq
	if( payFreq == K_DEF_FREQ )
		payFreq = resetFreq;

	/// some validation
	/// to be exception safe validation has to be
	/// the first thing
	if( swapType == K_ZEROCOUPON_LEG && 
		( resetFreq != K_ZEROCOUPON ||payFreq != K_ZEROCOUPON ) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			   "Zero Coupon Leg with not ZERO Coupon freq in reset and pay");


	/// because we do not have any date Strip
	ARM_DateStrip* tmpNumDateStrip;
	
	/// the logic of datestrip is to use roll date if and only
	///	if the frequency is not annual and we are given a roll date
	ARM_Date numRefStartDate;
	if(resetFreq!=K_ZEROCOUPON )
	{
		numRefStartDate=startDate;
		numRefStartDate.AddPeriod(resetFreq,const_cast<char*>(resetCalendar));
	}
	else
		numRefStartDate=endDate;
	int numGap = ARM_InfLeg::GetGapFromGapOrDate(resetNumGaporRollDate,resetCalendar,numRefStartDate);

	tmpNumDateStrip= new ARM_DateStrip(
		startDate, endDate, resetFreq,
		dayCount, resetCalendar,
		resetRule, intRule, stubRule,
		numGap, payFreq, payGap, ourPayCalendar ,
		K_ARREARS, K_ARREARS, adjFirstDate );

	if(		resetFreq!=K_ANNUAL && resetFreq!=K_ZEROCOUPON
		&&  resetNumGaporRollDate>JULIANDATEMIN)
	{
		int period=numGap<0? -K_DAILY:K_DAILY;
		ARM_Date numStartDate=startDate;
		numStartDate.AddPeriodMult(period,abs(numGap),const_cast<char*>(resetCalendar));
		ARM_Date numEndDate = endDate;
		numEndDate.AddPeriodMult(period,abs(numGap),const_cast<char*>(resetCalendar));

		/// this datestrip is only to compute the reset date!
		/// a little bit brute force... but that is the way it is!
		ARM_DateStrip resetDateStrip(
			numStartDate, numEndDate, resetFreq,
			dayCount, resetCalendar,
			resetRule, intRule, stubRule,
			0.0, payFreq, payGap, ourPayCalendar ,
			K_ARREARS, K_ARREARS, adjFirstDate );

		/// change the reset dates!
		tmpNumDateStrip->SetResetDates( resetDateStrip.GetResetDates() );
	};

	/// Sets NumPublicationDates (for calibration purpose)
	ARM_GP_Vector* numReset = tmpNumDateStrip->GetResetDates();
	size_t i,CFSize = numReset->size();
	itsNumPublicationDates = new ARM_GP_Vector( CFSize );
	for(i=0; i<CFSize; i++)
		(*itsNumPublicationDates)[i] = (*numReset)[i]-numGap;

	ARM_DateStrip* tmpDenomDateStrip;
	int denomGap = ARM_InfLeg::GetGapFromGapOrDate(resetDenomGaporRollDate,resetCalendar,startDate);

	switch( itsSwapType )
	{
	case K_YEARTOYEAR_LEG:
		{
			itsMultiple = multiple;
			itsCoMultiple = Comultiple;
			itsConstant = -multiple + constant;
			tmpDenomDateStrip = new ARM_DateStrip(
				startDate, endDate, resetFreq,
				dayCount, resetCalendar,
				resetRule, intRule, stubRule,
				denomGap, payFreq, payGap, ourPayCalendar ,
				K_ADVANCE, K_ARREARS, adjFirstDate );
			
			if(		resetFreq!=K_ANNUAL && resetFreq!=K_ZEROCOUPON
				&&	resetDenomGaporRollDate>JULIANDATEMIN)
			{
				ARM_Date denomStartDate( (double)resetDenomGaporRollDate);

				ARM_Date denomEndDate=endDate;
				int period=denomGap<0? -K_DAILY:K_DAILY;
				denomEndDate.AddPeriodMult(period,abs(denomGap),const_cast<char*>(resetCalendar));

				ARM_DateStrip resetDateStrip(
					denomStartDate, denomEndDate, resetFreq,
					dayCount, resetCalendar,
					resetRule, intRule, stubRule,
					0.0, payFreq, payGap, ourPayCalendar,
					K_ADVANCE, K_ARREARS, adjFirstDate );
				
				/// change the reset dates!
				tmpDenomDateStrip->SetResetDates( resetDateStrip.GetResetDates() );
			}
	
			break;
		}
	case K_ZEROCOUPON_LEG:
		{
			itsMultiple = multiple;
			itsCoMultiple = Comultiple;
			itsConstant = -multiple + constant;
			
			tmpDenomDateStrip = new ARM_DateStrip(
				startDate, endDate, resetFreq,
				dayCount, resetCalendar,
				resetRule, intRule, stubRule,
				denomGap, payFreq, payGap, ourPayCalendar ,
				K_ADVANCE, K_ARREARS, adjFirstDate );

			if(		resetFreq!=K_ANNUAL && resetFreq!=K_ZEROCOUPON
				&&	resetDenomGaporRollDate>JULIANDATEMIN)
			{
				ARM_Date denomStartDate( (double)resetDenomGaporRollDate);

				ARM_Date denomEndDate=endDate;
				int period=denomGap<0? -K_DAILY:K_DAILY;
				denomEndDate.AddPeriodMult(period,abs(denomGap),const_cast<char*>(resetCalendar));

				ARM_DateStrip resetDateStrip(
					denomStartDate, denomEndDate, resetFreq,
					dayCount, resetCalendar,
					resetRule, intRule, stubRule,
					0.0, payFreq, payGap, ourPayCalendar ,
					K_ADVANCE, K_ARREARS, adjFirstDate );

				/// change the reset dates!
				tmpDenomDateStrip->SetResetDates( resetDateStrip.GetResetDates() );
			}

			break;

		}
	case K_OATTYPE_LEG:
		{
			// multiple is given in 100 basis hence we need to divide
			itsMultiple = multiple/CC_NS( ARM_Constants, rateBase );
			itsCoMultiple = Comultiple;
			itsConstant = constant;

			tmpDenomDateStrip = new ARM_DateStrip(
				startDate, endDate, resetFreq,
				dayCount, resetCalendar,
				resetRule, intRule, stubRule,
				denomGap, payFreq, payGap, ourPayCalendar ,
				K_ADVANCE, K_ARREARS, adjFirstDate );
			
			double T0 = tmpDenomDateStrip->GetResetDates()->Elt(0);
			int size = tmpDenomDateStrip->GetResetDates()->size();
			
			ARM_GP_Vector* tmpVec = new ARM_GP_Vector( size, T0 );
			tmpDenomDateStrip->SetResetDates( tmpVec );
			
			break;
			
		}
	case K_GENERICINF_LEG: default:
		{
			itsMultiple = multiple;
			itsCoMultiple = Comultiple;
			itsConstant = constant;
			tmpDenomDateStrip = new ARM_DateStrip(
				startDate, endDate, resetFreq,
				dayCount, resetCalendar,
				resetRule, intRule, stubRule,
				denomGap, payFreq, payGap, ourPayCalendar,
				K_ADVANCE, K_ARREARS, adjFirstDate );

			if(		resetFreq!=K_ANNUAL && resetFreq!=K_ZEROCOUPON
				&&	resetDenomGaporRollDate>JULIANDATEMIN)
			{
				ARM_Date denomStartDate( (double)resetDenomGaporRollDate);

				ARM_Date denomEndDate=endDate;
				int period=denomGap<0? -K_DAILY:K_DAILY;
				denomEndDate.AddPeriodMult(period,abs(denomGap),const_cast<char*>(resetCalendar));

				ARM_DateStrip resetDateStrip(
					denomStartDate, denomEndDate, resetFreq,
					dayCount, resetCalendar,
					resetRule, intRule, stubRule,
					0.0, payFreq, payGap, ourPayCalendar,
					K_ADVANCE, K_ARREARS, adjFirstDate );

				/// change the reset dates!
				tmpDenomDateStrip->SetResetDates( resetDateStrip.GetResetDates() );
			}

			break;
		}

	}
	tmpInfIdx->SetPayFrequency(payFreq);
	tmpInfIdx->SetResetFrequency(resetFreq);
	tmpInfIdx->SetIntRule(intRule);
	tmpInfIdx->SetDayCount(dayCount);
	SetDayCount(dayCount);

	Set( tmpInfIdx , rcvOrPay, finalNotionalType, discountCcy, tmpNumDateStrip, tmpDenomDateStrip  );


	delete tmpInfIdx;
	delete tmpNumDateStrip;
	delete tmpDenomDateStrip;

	SetName( ARM_INFLEG );
};


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : CopyNoCleanUp to factorise code
///				WARNING: this copies only the part of the 
///				Derived object and not the base object
////////////////////////////////////////////////////

void ARM_InfLeg::CopyNoCleanUp( const ARM_InfLeg& infLeg )
{
	itsSwapType			= infLeg.itsSwapType;
	itsInterpType		= infLeg.itsInterpType;
	itsMultiple			= infLeg.itsMultiple;
	itsCoMultiple		= infLeg.itsCoMultiple;
	itsConstant			= infLeg.itsConstant;
	itsFirstReset		= infLeg.itsFirstReset;
	
	itsPayCurrencyName	= infLeg.itsPayCurrencyName;
	itsIndexName		= infLeg.itsIndexName;
	itsResetNumGaporRollDate =	infLeg.itsResetNumGaporRollDate;
	itsResetDenomGaporRollDate = infLeg.itsResetDenomGaporRollDate;

	itsNumResetDates	= new ARM_GP_Vector( *infLeg.itsNumResetDates );
	itsDenomResetDates	= new ARM_GP_Vector( *infLeg.itsDenomResetDates );
	itsNumCPIRates		= new ARM_GP_Vector( *infLeg.itsNumCPIRates );
	itsDenomCPIRates	= new ARM_GP_Vector( *infLeg.itsDenomCPIRates );
	itsDiscountFactors	= new ARM_GP_Vector( *infLeg.itsDiscountFactors );
	if( infLeg.itsNumPublicationDates )
		itsNumPublicationDates = new ARM_GP_Vector( *infLeg.itsNumPublicationDates );
	else
		itsNumPublicationDates = NULL;

}


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: Copy constructor
///	Returns: 
///	Action : Copy constructor
////////////////////////////////////////////////////

ARM_InfLeg::ARM_InfLeg( const ARM_InfLeg& infLeg )
	/// Base class copy construction
:	ARM_SwapLeg( infLeg )
{
	CopyNoCleanUp( infLeg );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: operator =
///	Returns: 
///	Action : operator =
////////////////////////////////////////////////////

ARM_InfLeg& ARM_InfLeg::operator = ( const ARM_InfLeg& infLeg )
{
	//check for self assignment
	if( this != &infLeg )
	{
		/// delete used pointors
		CleanUp();

		ARM_SwapLeg::operator=( infLeg );
		/// assign to all data members
		CopyNoCleanUp( infLeg );	
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: CleanUp
///	Returns: 
///	Action : CleanUp function
////////////////////////////////////////////////////

void ARM_InfLeg::CleanUp()
{
	delete itsNumResetDates;
	itsNumResetDates = NULL;
	delete itsDenomResetDates;
	itsDenomResetDates = NULL;
	delete itsNumCPIRates;
	itsNumCPIRates = NULL;
	delete itsDenomCPIRates;
	itsDenomCPIRates = NULL;
	delete itsDiscountFactors;
	itsDiscountFactors= NULL;
	delete itsNumPublicationDates;
	itsNumPublicationDates = NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_InfLeg::~ARM_InfLeg()
{
	CleanUp();
}



////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: Clone
///	Returns: 
///	Action : Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfLeg::Clone(void)
{
	// this is faster ... hence prefered to the old ARM solution
	return new ARM_InfLeg(*this);
}




///////////////////////////////////////////////////////
///			PRICING AND VIEW PART
//////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: CptExpectedFwdRates
///	Returns: 
///	Action : Compute the expected fwd Rates
////////////////////////////////////////////////////
void ARM_InfLeg::CptExpectedFwdRates()
{
	/// Get the model
	ARM_Model* model = GetModel();

	/// we test if the functionality of the model at least allow
	/// to price fwd CPI
	if( model->CanPriceInflation() >= PRICE_FWD_CPI )
	{

		InfFwdModel* pricingMod = dynamic_cast<InfFwdModel*>(model);
		if( !pricingMod )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				"model is of wrong type.. is supposed to price CPI fwd but failed to do so because it does  not inherit from InfFwdModel");

		int nFlows = itsNumResetDates->size();
		
		/// by construction 
		/// itsNumCPIRates and itsDenomCPIRates
		/// are either null or initialised with new
		delete itsNumCPIRates;
		delete itsDenomCPIRates;

		ARM_GP_Vector tmpFwds( nFlows, 0.0);
		itsNumCPIRates		= new ARM_GP_Vector( nFlows, 0.0);
		itsDenomCPIRates	= new ARM_GP_Vector( nFlows, 0.0);
		ARM_GP_Vector CPIRatio;
		ARM_InfIdx* infIdx	= (ARM_InfIdx*) GetIRIndex();
		if( !infIdx )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				ARM_USERNAME + ": infIdx is null");

		
		/// the first reset can be overwritten
		/// but the function FwdCPIRatio which returns an ARM_GP_Vector
		/// of dimension 3 can take this overwritten reset!
		int i;
		
		for( i=0; i<nFlows; ++i )
		{
			/// FwdCPIRatio returns in the following order
			/// CPIRatio[0] is the ratio of CPI
			/// CPIRatio[1] is the numerator CPI
			/// CPIRatio[2] is the denominator CPI

			/// the convention is that if the last argument is different from GETDEFAULTVALUE 
			/// then the reset is used!
			CPIRatio = pricingMod->FwdCPIRatio( 
				itsNumResetDates->Elt(i), 
				itsDenomResetDates->Elt(i),
				(*GetPaymentDates())[i],
				itsMultiple, 
				itsConstant, 
				itsInterpType, 
				itsFirstReset, 
				infIdx );
			tmpFwds[i]				= GetRcvOrPay() * (CPIRatio[0] * CC_NS( ARM_Constants, rateBase ) + GetSpread())* GetInterestTerms()->Elt(i);
			itsNumCPIRates->Elt(i)	= CPIRatio[1];
			itsDenomCPIRates->Elt(i)= CPIRatio[2];
		}

		switch( GetNxFlag() )
		{
		case K_NX_END:
			tmpFwds[nFlows-1] += CC_NS( ARM_Constants, rateBase );
			break;
		case K_NX_INFEND:
			{
				CPIRatio = pricingMod->FwdCPIRatio( 
					itsNumResetDates->Elt(nFlows-1), 
					itsDenomResetDates->Elt(0),
					(*GetPaymentDates())[nFlows-1],
					1.0, 
					0.0, 
					itsInterpType, 
					itsFirstReset, 
					(ARM_InfIdx*) GetIRIndex() );
				tmpFwds[nFlows-1] += GetRcvOrPay() * CPIRatio[0] * CC_NS( ARM_Constants, rateBase );
				break;
			}
		case K_NX_ASINFEND:
			{
				CPIRatio = pricingMod->FwdCPIRatio( 
					itsNumResetDates->Elt(nFlows-1), 
					itsDenomResetDates->Elt(0),
					(*GetPaymentDates())[nFlows-1],
					1.0, 
					-1.0, 
					itsInterpType, 
					itsFirstReset, 
					(ARM_InfIdx*) GetIRIndex()  
				);
				tmpFwds[nFlows-1] += GetRcvOrPay() * CPIRatio[0] * CC_NS( ARM_Constants, rateBase );
				break;
			}
		}
		
		ARM_Vector tmpFwds2 = To_ARM_Vector(tmpFwds);
		SetFwdRates(&tmpFwds2);		
		/// properly clean the tmp variables
		/// since SetFwdRates clones the ARM_Vector*
	}
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		 ARM_USERNAME + ": model is of wrong type.. cannot price fwd CPI");
}


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: CptCashFlowValues
///	Returns: 
///	Action : Compute the expected cash flows
////////////////////////////////////////////////////
void ARM_InfLeg::CptCashFlowValues()
{
	/// Test the amotistization
	int amortType = this->GetAmount()->GetValueType();

    if ( ( amortType == 0 )
         && ( this->GetIRIndex()->GetPayFrequency() ) 
         < this->GetAmount()->GetCalcFrequency() )
    {
       this->CptCashFlowValuesIfAmortization();

       return;
    }

	double payJulDate, amount;

	int nFlows = itsNumResetDates->size();
	ARM_GP_Vector* tmpFlowValues = new ARM_GP_Vector( nFlows, 0.0);
    ARM_Vector* cfTerms			 = GetYearTerms();
	ARM_Model* model			 = GetModel();

	int i;
	for( i=0; i<nFlows; ++i )
	{
		/// ignore past cash flows
		if( cfTerms->Elt(i) < 0.0 )
		{
			itsDiscountFactors->Elt(i)	= 0.0;
			tmpFlowValues->Elt(i)		= 0.0;
		}
		else
		{
			// test zero curve to avoid crashing!
			if( model->GetZeroCurve()->GetDiscountFactors() )
			{
				itsDiscountFactors->Elt(i)	= model->ZeroPrice( 0.0, cfTerms->Elt(i), GetDiscountingYC() );
				tmpFlowValues->Elt(i)		= GetFwdRates()->Elt(i);
				if( itsSwapType == K_YEARTOYEAR_LEG )
				{
					payJulDate				=  GetPaymentDates()->Elt(i);
					amount					=  GetAmount()->CptReferenceValue(payJulDate)/100;
					tmpFlowValues->Elt(i)	*= amount;
				}
			}
			else
			{
				itsDiscountFactors->Elt(i)	= 0.0;
				tmpFlowValues->Elt(i)		= 0.0;
			}
		}
	}

    SetDecompRates(GetFwdRates());		/// Clone it so need to clone it
	ARM_Vector* tmpFlowValues2 = To_pARM_Vector(tmpFlowValues);
	delete tmpFlowValues;
    SetCashFlowValues(tmpFlowValues2);	/// do not clone it
}



////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: PropagateModel
///	Returns: 
///	Action : the propagate model simply computes the forward cashflow
////////////////////////////////////////////////////
void ARM_InfLeg::PropagateModel(ARM_Model* model )
{
    CptExpectedFwdRates();
    CptCashFlowValues();
}



////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: ComputePrice
///	Returns: 
///	Action : Compute the Price
////////////////////////////////////////////////////

double ARM_InfLeg::ComputePrice(int mode)
{
	if( itsSwapType == K_YEARTOYEAR_LEG )
		UpdateLegAmortization();
	
	double price = 0.0;
	size_t nFlows = itsNumResetDates->size();

	for( size_t i=0.0; i<nFlows; ++i )
		price += itsDiscountFactors->Elt(i) * GetCashFlowValues()->Elt(i);
	SetPrice(price);
	
	return(price);
}

////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: View
///	Returns: 
///	Action : View method
////////////////////////////////////////////////////

void ARM_InfLeg::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;
	
	/// since the view method is common
	/// for all inflationswap
	/// get the type of the swap
    fprintf( fOut, "\n\n =======> INFLATION SWAP LEG  <====== \nInflation Leg Type : %s \n", ARM_ParamView::GetMappingName( S_LEGTYPE, itsSwapType ) );
	
	if( GetIRIndex() )
		GetIRIndex()->View(id, fOut);
	
	// receive or pay
    fprintf(fOut, " Receive or Pay    : %s\n\n", ARM_ParamView::GetMappingName( S_RECEIVE_PAY, GetRcvOrPay() ) );
	
	int nFlows = (GetFlowStartDates())->GetSize();
	
	char d1[20];
	char d2[20];
	char d3[20];
	char d4[20];
	char d5[20];
	
	int  i;
	
	/// for readability split the fprintf function
	fprintf(fOut, "\n Start Dates\t End Dates\t Interest Days\t Terms");
	fprintf(fOut, "\t Num Dates\t Num Resets\t Denom Dates\t Denom Resets");
	fprintf(fOut, "\t Payment Dates\t DF \t\tCF Forwards\tNominal\n ");
	
	for (i = 0; i < nFlows; ++i )
	{
		((ARM_Date) (*GetFlowStartDates())[i]).JulianToStrDate(d1);
		((ARM_Date) (*GetFlowEndDates())[i]).JulianToStrDate(d2);
		((ARM_Date) (*GetNumResetDates())[i]).JulianToStrDate(d3);
		((ARM_Date) (*itsDenomResetDates)[i]).JulianToStrDate(d4);
		((ARM_Date) (*GetPaymentDates())[i]).JulianToStrDate(d5);
		
		fprintf(fOut, " %s\t %s\t %3.0lf\t\t %.4lf", 
			d1, d2, (*GetInterestDays())[i],  (*GetInterestTerms())[i] );
		fprintf(fOut, "\t %s\t %.5lf\t %s\t %.5lf", 
			d3, itsNumCPIRates->Elt(i), d4, itsDenomCPIRates->Elt(i) );
		fprintf(fOut, "\t %s\t %.7lf \t%.7lf \t%.2lf\n", 
			d5, itsDiscountFactors->Elt(i), GetCashFlowValues()->Elt(i),
			GetAmount()->CptReferenceValue((ARM_Date) GetPaymentDates()->Elt(i)) );

	}
	
	/// to allow to have nested view
	if ( ficOut == NULL )
		fclose(fOut);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: GetDiscountFactors
///	Returns: 
///	Action : Accessor to discount factors
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_InfLeg::GetDiscountFactors() const
{
	return itsDiscountFactors;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: GetNumResetDates
///	Returns: 
///	Action : Accessor to num reset dates
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_InfLeg::GetNumResetDates() const
{
	return itsNumResetDates;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: GetDenomResetDates
///	Returns: 
///	Action : Accessor to DenomResetDates
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_InfLeg::GetDenomResetDates() const
{
	return itsDenomResetDates;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: GetSwapType
///	Returns: 
///	Action : Accessor to SwapType
////////////////////////////////////////////////////
int ARM_InfLeg::GetSwapType() const
{
	return itsSwapType;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: GetNumCPIRates
///	Returns: 
///	Action : Accessor to NumCPIRates
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_InfLeg::GetNumCPIRates() const
{
	return itsNumCPIRates;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfLeg
///	Routine: GetDenomCPIRates
///	Returns: 
///	Action : Accessor to DenomCPIRates
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_InfLeg::GetDenomCPIRates() const
{
	return itsDenomCPIRates;
}



ARM_AverageSwapLeg::ARM_AverageSwapLeg(	ARM_Date&		startDate, 
										ARM_Date&		endDate, 
										ARM_IRIndex*	irIndex,
										int				rcvOrPay, 
										double			spread, 
										int				stub, 
										int				decompFreq,
										ARM_Currency*	ccyName, 
										int				dayCount,
										char*			payCalName,
										int				decompPricingFlag,
										int				nxChange,
										char*			refDate,
										int				adjStartDate,
										double			couru)
{
    if ( endDate <= startDate )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Start date must be > end date");
    }

    Init();

	itsCouru = couru;

    SetName(ARM_SWAPLEG);

    Set(	startDate,
			endDate, 
			irIndex, 
			rcvOrPay, 
			spread, 
			stub, 
			decompFreq, 
			ccyName, 
			dayCount,
			payCalName, 
			decompPricingFlag, 
			nxChange,
			refDate, 
			adjStartDate);
}

void ARM_AverageSwapLeg::Set(	ARM_Date& startDate, 
								ARM_Date& endDate, 
								ARM_IRIndex* irIndex, 
								int rcvOrPay, 
								double spread, 
								int stub, 
								int decompFreq,
								ARM_Currency* ccy, 
								int dayCount,
								char* payCalName,
								int decompPricingFlag,
								int nxChange,
								char* refDate,
								int adjStartDate)
{
    ARM_IRIndex* iri = NULL;

    itsAdjStartDateFlag = adjStartDate;

    //---> Reset and Pay calendars

    if (refDate)
    {
       if (itsRefDate)
       {
          delete itsRefDate;

          itsRefDate = NULL;
       }

       if (strcmp(refDate, "NULL"))
       {
          itsRefDate = new char [strlen(refDate)+1];

          strcpy(itsRefDate, refDate);
       }
    }
    else
    {
       if (itsRefDate)
       {
          delete itsRefDate;

          itsRefDate = NULL;
       }
    }
    char* tmpPayCal;


    itsNotionalExchangeFlag = nxChange;

    if ( payCalName == NULL )
    {
	   tmpPayCal = ccy->GetPayCalName(ccy->GetVanillaIndexType());
		
       strcpy(itsPayCalendar,	tmpPayCal);
       strcpy(itsResetCalendar, tmpPayCal);

       delete tmpPayCal;
    }
    else
    {
       strcpy(itsPayCalendar, payCalName);
    }


    itsDecompPricingFlag = decompPricingFlag;

    // End Calendars


    itsStartDateNA = startDate;

    itsEndDateNA  = endDate;


    itsReceiveOrPay = rcvOrPay;

    itsDiscountingYC = K_DOMESTIC;

    if (irIndex)
    {
       SetIRIndex(irIndex);
 
       itsSpread = spread;
       
       if ( irIndex->GetIndexType() == IDXFIXED )
       {
          itsLegType = K_FIXED_LEG;

          itsIndexType = IDXFIXED;

          itsFixedRate = spread;

          itsSpread    = 0.0;
       }
       else
       { 
          itsFixedRate = K_FLOAT_RATE;
          itsLegType   = K_FLOATING_LEG;
       }

       itsDecompFrequency = irIndex->GetDecompFreq();
    }
    else
    {
       iri = new ARM_IRIndex();

       SetIRIndex(iri);

       delete iri;

       itsSpread = 0.0;

       itsFixedRate = spread;
       itsLegType   = K_FIXED_LEG;

       itsDecompFrequency = decompFreq;
    }

		ARM_Date	tmp;
		int			m		= startDate.GetMonth();
		int			y		= startDate.GetYear();
		char*		tmpCcy	= ccy->GetCcyName();
		
		if ( m <7 ) 
			tmp = ARM_Date(1, 12, y-1); 
		else
			tmp = ARM_Date (1, 6, y); 

		if ( 	!tmp.IsBusinessDay(	 tmpCcy )  )
					tmp = tmp.NextBusinessDay( tmpCcy ) ;
			
		int resetGap = -CountBusinessDays(tmp,	startDate, tmpCcy);

 if ( resetGap <= 0 )
       itsIRIndex->SetResetGap(resetGap);

    itsPaymentFreq = irIndex->GetPayFrequency();

    itsIndexType = irIndex->GetIndexType();

    if ( irIndex->GetIndexStyle() == IN_ARREARS )
    {
       itsIndexStyle = IN_ARREARS;
    }
     

    int fwdRule = itsIRIndex->GetFwdRule();

    char ccyName[4];

    if (ccy)
       strcpy(ccyName, ccy->GetCcyName());
    else if (ARM_DEFAULT_CURRENCY)
    {
       strcpy(ccyName, ARM_DEFAULT_CURRENCY->GetCcyName());
    }
    else
    {
       strcpy(ccyName, ARM_DEFAULT_COUNTRY);
    }

    itsStartDate = startDate;

    if (adjStartDate)
    {
       if ( GetIRIndex()->GetIntRule() == K_ADJUSTED )
          itsStartDate.GoodBusinessDay(K_FOLLOWING, ccyName);
    }

    SetEndDate(endDate);

    ARM_Date YCmaturity = itsIRIndex->SettlementToMaturity(endDate);

    SetMaturity(YCmaturity);
    
    SetYCSensitivityDate(YCmaturity);

    if (ccy)
    {
       SetCurrencyUnit(ccy);

       itsDiscountingYC = DiscountWithDomYC();
    }
        
    InitQuanto();

        
    if ( dayCount == -1 )
       SetDayCount(itsIRIndex->GetDayCount());
    else
       SetDayCount(dayCount);

    itsStubMeth = stub;

    SetNumFlows(0);
    
    // compute cash flows Start, End Reset and Payment dates
        
    CptCashFlowDates();

    UpdateLegAmortization();

	
	ARM_Date indDate= ARM_Date( GetResetDates()->Elt(0) );
		if ( GetResetDates()->size() >1 ){
			for ( int i=0; i< GetResetDates()->size();i++){
				tmp = ARM_Date( GetResetDates()->Elt(i) );
				if ( tmp.GetMonth()== 6 || tmp.GetMonth()== 12) 
					indDate= ARM_Date( GetResetDates()->Elt(i) );
				else
					GetResetDates()->Elt(i)=indDate.GetJulian();

			}
		}    
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/


