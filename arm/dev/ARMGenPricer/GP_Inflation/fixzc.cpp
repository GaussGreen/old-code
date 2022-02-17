/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file fixzc.cpp
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \brief implements a fixed zero coupon leg
 *		Inherit from ARM_SwapLeg
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

#include "gpinflation/fixzc.h"
#include "gpinflation/infdata.h"

/// gpbse
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"

#include <glob/paramview.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_FixZC
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_FixZC::ARM_FixZC( 
		const ARM_Date& startDate,
		const ARM_Date& endDate,
		double fixRate,
		int rcvOrPay,
		int dayCount,
		int intRule,
		int stubRule,
		int payGap,
		const char* payCalendar,
		const ARM_Currency*	discountCcy )
:
	ARM_SwapLeg( 
		const_cast<ARM_Date&>(startDate),
		const_cast<ARM_Date&>(endDate),
		fixRate,
		rcvOrPay,
		K_ZEROCOUPON,
		dayCount,
		K_COMP_PROP,
        K_ARREARS,
        intRule,
        stubRule,
        const_cast<ARM_Currency*>(discountCcy),
		const_cast<char*>(payCalendar)  )
{
	/// rewrittes properly the dates
	ARM_Date tmpPaymentDate	= endDate;
	tmpPaymentDate.GapBusinessDay( payGap, const_cast<char*>(payCalendar) );
	ARM_Vector* tmpPaymentDates  = new ARM_Vector( 1, tmpPaymentDate.GetJulian() );
	SetPaymentDates( tmpPaymentDates );

	ARM_Vector* tmpResetDates = new ARM_Vector( 1, endDate.GetJulian() );
	SetResetDates( tmpResetDates );
	delete tmpResetDates;
	
	itsDiscountFactors = new ARM_GP_Vector( 1.0, GETDEFAULTVALUE  );

	/// interest Terms
	ARM_Vector* tmpInterestTermValues = new ARM_Vector( 1.0, 1.0 );
	SetInterestTermValues( tmpInterestTermValues );
	delete tmpInterestTermValues;

	ARM_Vector* tmpCashFlowValues = new ARM_Vector( 1.0, 0.0 );
	SetCashFlowValues(tmpCashFlowValues	);	/// do not clone it so do not delete it!

	SetName( ARM_FIXZC );
}



////////////////////////////////////////////////////
///	Class  : ARM_FixZC
///	Routine: Copy constructor, assignement operator
///	Returns: 
///	Action : Copy constructor, assignement operator
////////////////////////////////////////////////////
ARM_FixZC::ARM_FixZC( const ARM_FixZC& rhs )
: ARM_SwapLeg( rhs ), itsDiscountFactors( new ARM_GP_Vector( *rhs.itsDiscountFactors ) )
{}


ARM_FixZC& ARM_FixZC::operator=( const ARM_FixZC& rhs)
{
	if( this !=	 &rhs )
	{
		delete itsDiscountFactors;
		ARM_SwapLeg::operator = ( rhs );
		itsDiscountFactors = rhs.itsDiscountFactors? (ARM_GP_Vector*) rhs.itsDiscountFactors->Clone() : NULL;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_FixZC
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_FixZC::~ARM_FixZC()
{
	delete itsDiscountFactors;
	itsDiscountFactors = NULL;
}


	
	
////////////////////////////////////////////////////
///	Class  : ARM_FixZC
///	Routine: Clone
///	Returns: 
///	Action : Clone method
////////////////////////////////////////////////////
ARM_Object* ARM_FixZC::Clone()
{
	return new ARM_FixZC( * this );
}


////////////////////////////////////////////////////
///	Class  : ARM_FixZC
///	Routine: ComputePrice
///	Returns: 
///	Action : call the model to compute the discounted cashflows
////////////////////////////////////////////////////
double ARM_FixZC::ComputePrice(int)
{
	double price = GetCashFlowValues()->Elt(0) * itsDiscountFactors->Elt(0);
    SetPrice(price);
    return(price);
}

////////////////////////////////////////////////////
///	Class  : ARM_FixZC
///	Routine: PropagateModel
///	Returns: 
///	Action : the propagate model simply computes the forward cashflow
////////////////////////////////////////////////////
void ARM_FixZC::PropagateModel(ARM_Model* )
{
    CptCashFlowValues();
}


////////////////////////////////////////////////////
///	Class  : ARM_FixZC
///	Routine: CptCashFlowValues
///	Returns: 
///	Action : Computation of the single forward cashflow
////////////////////////////////////////////////////

void ARM_FixZC::CptCashFlowValues()
{
	int numFlows = GetFlowStartDates()->GetSize();
	double interestTerm = 0.0;

	for (int i = 0; i < numFlows; i++)
	{
		interestTerm += GetInterestTerms()->Elt(i);
	}

	ARM_GP_Vector* tmpFlowValues	= new ARM_GP_Vector( 1, 0.0);
	tmpFlowValues->Elt(0) = GetRcvOrPay() * ( pow( 
		(1.0+ GetFixedRate() / CC_NS( ARM_Constants, rateBase ) ), interestTerm ) - 1.0 ) * CC_NS( ARM_Constants, rateBase );
	
    ARM_Model* model = GetModel(); 

    if ( model == NULL )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Yield Curve Model must be set before <ComputePrice> method");

	itsDiscountFactors->Elt(0)	= ( model->ZeroPrice( 0.0, GetYearTerms()->Elt(0), GetDiscountingYC() ) );
    
	ARM_Vector* tmpFlowValues2 = To_pARM_Vector(tmpFlowValues);
	delete tmpFlowValues;
	SetCashFlowValues(tmpFlowValues2); /// do not clone it
}



////////////////////////////////////////////////////
///	Class  : ARM_FixZC
///	Routine: View
///	Returns: 
///	Action : View method
////////////////////////////////////////////////////

void ARM_FixZC::View(char* id, FILE* ficOut)
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
    fprintf(fOut, "\n\n =======> FIX ZERO COUPON LEG  <====== \n" );

	/// currency
    fprintf(fOut, "\n Currency\t: %s", GetCurrencyUnit()->GetCcyName() );
	/// fix rate
	fprintf(fOut, "\n Fix Rate\t: %.2f", GetFixedRate() );
	// receive or pay
    fprintf(fOut, "\n Rcv or Pay\t: %s\n\n", ARM_ParamView::GetMappingName( S_RECEIVE_PAY, GetRcvOrPay() ) );
	
	/// for readability split the fprintf function
	fprintf(fOut, "\n Start Dates\t End Dates\t Interest Days\t Terms");
	fprintf(fOut, "\t Payment Dates\t DF \t\tFwd CF\n");

	char d1[20];
	char d2[20];
	char d3[20];

	int numFlows = GetFlowStartDates()->GetSize();

	((ARM_Date) (*GetFlowStartDates())[0]).JulianToStrDate(d1);
	((ARM_Date) (*GetFlowEndDates())[numFlows-1]).JulianToStrDate(d2);
	((ARM_Date) (*GetPaymentDates())[0]).JulianToStrDate(d3);

	double interestTerm = 0.0;
	double interestDays = 0.0;

	for (int i = 0; i < numFlows; i++)
	{
		interestTerm += GetInterestTerms()->Elt(i);
		interestDays += GetInterestDays()->Elt(i);
	}
	
	fprintf(fOut, " %s\t %s\t %3.0lf\t\t %.2lf", 
		d1, d2, interestDays,  interestTerm);
	fprintf(fOut, "\t %s\t %.7lf \t%.7lf\n", 
		d3, itsDiscountFactors->Elt(0),  GetCashFlowValues()->Elt(0) );
	/// to allow to have nested view
    if ( ficOut == NULL )
       fclose(fOut);
}

CC_END_NAMESPACE()


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
