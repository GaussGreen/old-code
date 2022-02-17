/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file infcapfloor.cpp
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \brief object for inflation cap and floor
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */


#include "gpinflation/infcapfloor.h"

/// gpinflation
#include "gpinflation/infdata.h"
#include "gpinflation/infbsmodel.h"

/// gpbase
#include "gpbase/gpvector.h"
/// To conversion
#include "gpbase/gplinalgconvert.h"		/// for conversion

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_InfCapFloor::ARM_InfCapFloor( 
	const ARM_Date& startDate,
	const ARM_Date& endDate,
	const string& indexName,
	int capOrFloor,
	double strike,
	double leverage,
	double spread,
	int swapType,
	int rcvOrPay,
	int interpType,	
	int resetFreq,
	int dayCount,
	const char* resetCalendar,
	int fwdRule,
	int intRule,
	int stubRule,
	int resetGap,
	int payFreq,
	int payGap,
	const char* payCalendar,
	int adjFirstDate,
	double firstReset,
	const ARM_Currency*	discountCcy )
{
	double CoMultiple = 1.0;
	ARM_InfLeg* tmpInfleg = new ARM_InfLeg( startDate, endDate, indexName,
		swapType, rcvOrPay, interpType, leverage, CoMultiple, spread, resetFreq, dayCount,  resetCalendar,
		fwdRule, intRule, stubRule, resetGap, resetGap, payFreq, payGap, payCalendar, adjFirstDate,
		K_NX_NONE, firstReset, discountCcy );

	ARM_CapFloor::Set( tmpInfleg, capOrFloor, strike );

	Init();
	SetName( ARM_INFCAPFLOOR);
	delete tmpInfleg;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: Constructor
///	Returns: 
///	Action : Constructor from a swapLeg
////////////////////////////////////////////////////
ARM_InfCapFloor::ARM_InfCapFloor( ARM_SwapLeg* swapLeg, 
	int capFloor, 
	double strike )
: ARM_CapFloor( swapLeg, capFloor, strike )
{
	Init();
	SetName( ARM_INFCAPFLOOR);
}

////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: Constructor
///	Returns: 
///	Action : Constructor from a swapLeg
///				with a strike with a ref value
////////////////////////////////////////////////////
ARM_InfCapFloor::ARM_InfCapFloor( ARM_SwapLeg* swapLeg, int capFloor, ARM_ReferenceValue *strike)
: ARM_CapFloor( swapLeg, capFloor, strike )
{
	Init();
	SetName( ARM_INFCAPFLOOR);
}




////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: Init
///	Returns: 
///	Action : Init method
////////////////////////////////////////////////////
void ARM_InfCapFloor::Init(void)
{
	int nbCapLet = GetNumFlows();
	itsVolInfo   = new ARM_Matrix( nbCapLet, 5, 0.0 );
	itsCFValues	 = new ARM_GP_Vector( nbCapLet, 0.0 );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_InfCapFloor::ARM_InfCapFloor(const ARM_InfCapFloor& rhs)
: ARM_CapFloor( rhs )	/// should add the member of the class
{
	Init();
	CopyNoCleanUp( rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: Operator = 
///	Returns: 
///	Action : Operator = 
////////////////////////////////////////////////////
ARM_InfCapFloor& ARM_InfCapFloor::operator = (const ARM_InfCapFloor &rhs )
{
	if( this !=	 &rhs )
	{
		CleanUp();
		ARM_CapFloor::operator = ( rhs );
		CopyNoCleanUp( rhs );
	}

	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_InfCapFloor::~ARM_InfCapFloor()
{
	CleanUp();
}




////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: CleanUp
///	Returns: 
///	Action : CleanUp
////////////////////////////////////////////////////

void ARM_InfCapFloor::CleanUp()
{
	delete itsVolInfo;
	itsVolInfo = NULL;
	delete itsCFValues;
	itsCFValues = NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : CopyNoCleanUp
////////////////////////////////////////////////////

void ARM_InfCapFloor::CopyNoCleanUp( const ARM_InfCapFloor& rhs )
{
	itsVolInfo = (ARM_Matrix* ) rhs.itsVolInfo->Clone();
	itsCFValues	 = ( ARM_GP_Vector*) rhs.itsCFValues->Clone();
}




////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: Clone
///	Returns: 
///	Action : Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfCapFloor::Clone()
{
	return new ARM_InfCapFloor( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: View
///	Returns: 
///	Action : View method for inf cap and floor
////////////////////////////////////////////////////
void ARM_InfCapFloor::View(char* id, FILE* ficOut)
{
	ARM_SwapLeg* tmpInfLeg = GetSwapLeg();
	ARM_InfLeg* infLeg = dynamic_cast<ARM_InfLeg*>(tmpInfLeg);
	if( !infLeg )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"inflation cap/floor use a leg that is not an inflation leg... Please advise");

	
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
	fprintf(fOut, "\n\n =======> INFLATION CAP AND FLOOR  <====== \n" );
	
	fprintf(fOut, "\n\n Corresponding leg\n" );
	
	GetSwapLeg()->View(id, fOut );

	/// Pricing information
	fprintf(fOut, "\n\n =======> PRICING INFORMATION  <====== \n\n" );

	fprintf(fOut, "Denom Dates\t Num Dates\t Tenor\t Raw Strike\t Lookup Strike\t Forward\t Vol Level\t ForwardOption\n" );

	int nbCapLet = GetNumFlows();
	int i;
	char d1[20];
	char d2[20];
	double numDate, denomDate;

	for (i = 0; i < nbCapLet; ++i )
	{
		numDate			= infLeg->GetNumResetDates()->Elt(i);
		denomDate		= infLeg->GetDenomResetDates()->Elt(i);
		((ARM_Date) denomDate).JulianToStrDate(d1);
		((ARM_Date) numDate).JulianToStrDate(d2);
		fprintf(fOut, "%s\t %s\t %.3lf\t %.3lf\t\t %.6lf\t %.6lf\t %.6lf\t %.6lf\n", d1, d2, 
			itsVolInfo->Elt(i,4),
			itsVolInfo->Elt(i,2), 
			itsVolInfo->Elt(i,1), 
			infLeg->GetCashFlowValues()->Elt(i),
			itsVolInfo->Elt(i,0),
			itsCFValues->Elt(i) 
		);
	}
	
	/// to allow to have nested view
	if ( ficOut == NULL )
		fclose(fOut);
}




////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: PropagateModel
///	Returns: 
///	Action : the propagate model simply computes the forward cashflow
////////////////////////////////////////////////////
void ARM_InfCapFloor::PropagateModel(ARM_Model* model )
{
	GetSwapLeg()->SetModel( model );

	CptCashFlowValues();
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: CptCashFlowValues
///	Returns: 
///	Action : Computation of the single forward cashflow
////////////////////////////////////////////////////
void ARM_InfCapFloor::CptCashFlowValues()
{
	/// nothing right now
	
	ARM_Vector* tmpCashFlows  = To_pARM_Vector(itsCFValues);
	SetCashFlowValues(tmpCashFlows);
	
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: StoreVol
///	Returns: 
///	Action : Store Vol information
////////////////////////////////////////////////////
void ARM_InfCapFloor::StoreVol( double vol, double volLookupStrike, double strike, double timeToStart, double tenor, int k )
{
	itsVolInfo->Elt(k,0 ) = vol;
	itsVolInfo->Elt(k,1 ) = volLookupStrike;
	itsVolInfo->Elt(k,2 ) = strike;
	itsVolInfo->Elt(k,3 ) = timeToStart;
	itsVolInfo->Elt(k,4 ) = tenor;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: GetConstantWithNoSpread
///	Returns: 
///	Action : Get the spread ignoring spread effect!
////////////////////////////////////////////////////
double ARM_InfCapFloor::GetConstantWithNoSpread(int swapType) const
{
	switch(swapType)
	{
	case K_ZEROCOUPON_LEG: case K_YEARTOYEAR_LEG:
		{
			ARM_InfLeg* infLeg = (ARM_InfLeg*) (GetSwapLeg());
			double constant = infLeg->GetConstant();
			return constant;
		}
		
	case K_OATTYPE_LEG:
		return -0.0;
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"unknown type of option... can either be a zeroCoupon or a Year to year.. Please advise");
	}
}	



////////////////////////////////////////////////////
///	Class  : ARM_InfCapFloor
///	Routine: ComputePrice
///	Returns: 
///	Action : Computes cap/floor price at itsSetllement
////////////////////////////////////////////////////
double ARM_InfCapFloor::ComputePrice(int)
{
	ARM_Model* model = GetModel();

	if( model->CanPriceInflation() >= PRICE_FWDNOPTION )
	{
		InfOptionModel* pricingMod = dynamic_cast<InfOptionModel*>(model);
		if( !pricingMod )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": model is of wrong type.. is supposed to price option on CPI but failed to do so because it does  not inherit from InfOptionModel.. Please advise");
		
		ARM_SwapLeg* tmpInfLeg = GetSwapLeg();
		ARM_InfLeg* infLeg = dynamic_cast<ARM_InfLeg*>(tmpInfLeg);
		if( !infLeg )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": try to price an inflation cap/floor on a leg that is not an inflation leg... Please advise");

		/// this code
		double price = 0.0, optLet, CPIForward, timeToStart, tenor, discountfactor;
		int nbCapLet = GetNumFlows();

		int callPut			= GetCapFloorType();
		ARM_ReferenceValue* strikeProfile = (ARM_ReferenceValue*) GetStrikes()->Clone();
		double strike;
		if(strikeProfile==NULL)
		{
			strike = GetStrike();
			strikeProfile = new ARM_ReferenceValue(strike);
		}
				
		int swapType		= infLeg->GetSwapType();
		double constant		= GetConstantWithNoSpread(swapType);
		ARM_InfIdx* infidx	= (ARM_InfIdx*) infLeg->GetIRIndex();
		if( !infidx )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": null inflation index!");

		int infRoP	= infLeg->GetRcvOrPay();

		ARM_Date numDate, denomDate;
		int k;
		for (k = 0; k<nbCapLet; k++) 
		{
			numDate			= ARM_Date( infLeg->GetNumResetDates()->Elt(k) );
			denomDate		= ARM_Date( infLeg->GetDenomResetDates()->Elt(k) );

			double interestTerm	= infLeg->GetInterestTerms()->Elt(k);
			/// we need to add 1 to the forward since it is really the ratio of CPI
			/// however the Constant of the infleg is not anymore 1 if it is a leg with a spread, hence the notation!
			CPIForward		= infLeg->GetFwdRates()->Elt(k) / CC_NS( ARM_Constants, rateBase )/(infRoP*interestTerm) - constant;
			discountfactor	= infLeg->GetDiscountFactors()->Elt(k);
			double payTime	= infLeg->GetPaymentDates()->Elt(k);

			/// to account for fixing difference!
			double	renormalisationFactor = 1.0;
			timeToStart		= pricingMod->GetModelTimeWPublishLag( denomDate, (ARM_InfIdx*) infLeg->GetIRIndex() );
			if( timeToStart < 0.0 )
			{
				/// in this particular case we need to account for a fixing different from the one
				/// of the zero coupon option
				/// the renormalisation is YtYFixingCPI/CurrentCPI
				/// or the denomCPIRates / Reference CPI of the Curve
				renormalisationFactor = infLeg->GetDenomCPIRates()->Elt(k) / pricingMod->GetCPIIndexValue();
			}

			double timeToExpiry		= pricingMod->GetModelTimeWPublishLag( numDate, (ARM_InfIdx*) infLeg->GetIRIndex() );
			tenor = timeToExpiry	- timeToStart;
			double remainingTenor	= tenor;
			if( timeToStart < 0.0 )
				remainingTenor = remainingTenor+timeToStart <0? 0 : remainingTenor+timeToStart;

			strike = strikeProfile->CptReferenceValue(payTime);
		
			StoreVolInfo StoringFunc( strike, timeToStart, remainingTenor, k, this );

			/// to handle weird case of CPIForward negative!
			if( CPIForward < 0 )
			{
				strike		= -strike;
				callPut		= -GetCapFloorType();
				CPIForward	= -CPIForward;
			}

			ARM_InfCapFloorContext context( numDate, denomDate, swapType, renormalisationFactor );
			optLet				= pricingMod->SingleAssetOptionPrice( CPIForward, infRoP*strike, (infRoP*callPut), discountfactor, infidx, &context, StoringFunc );
			itsCFValues->Elt(k)	= optLet/discountfactor*interestTerm;
			price		   += optLet*interestTerm;
		}

		delete strikeProfile;
		SetPrice(price);
		return (price);
	}

	else
	{ 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"model is of wrong type.. cannot price option on CPI.. Please advise");
	}
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

