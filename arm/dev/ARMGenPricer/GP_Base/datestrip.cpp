/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file datestrip.cpp
 *	\brief Date strip is an object that can generate a strip of dates
 *		Handles only dates and compute 
 *			-the reset dates, 
 *			-the payment dates
 *			-the interest terms
 *
 *		This object is a very simple object that can be used in order
 *		to compute the date ... Use rather datestrip than swap leg
 *		if you are only interested in dates generations!
 *
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

#include "gpbase/datestrip.h"
#include "gpbase/gpvector.h"
#include "gpbase/argconvdefault.h"
#include "gpbase/globalconstant.h"

#include "gpbase/utilityport.h"  /// for CC_Min

#include <iomanip> /// for setprecision & co.
CC_USING_NS(std,setw);
CC_USING_NS(std,fixed);
CC_USING_NS(std,setprecision);

/// kernel
#include <glob/expt.h>
#include <ccy/currency.h>
#include <glob/paramview.h>


CC_BEGIN_NAMESPACE( ARM )

const int NB_BUSINESS_DAYS_PER_YEAR = 261;

struct EqualWithSevenDays : public CC_NS( std, binary_function )<double, double,bool>
{
    
	bool operator()( const double& value1, const double& value2) const
    {
        /// equality is only based on the type!
        return (fabs(value1 - value2) <= ARM_GlobalConstant::ARM_SEVENDAYS_LAG);
    }
};

/////////////////////////////////////////////////////////////////
///	static variable to say that it is a blank data!
/////////////////////////////////////////////////////////////////
const double ARM_DateStrip::DateStripCombiner_BlankData = -1.;
	
/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_DateStrip::Clone() const
{
	///faster than the standard ARM code
	return new ARM_DateStrip( *this );	
};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: Copy constructor
/////////////////////////////////////////////////////////////////
ARM_DateStrip::ARM_DateStrip( const ARM_DateStrip& rhs )
: ARM_RootObject( rhs ),
	itsStartDate( rhs.itsStartDate),
	itsEndDate( rhs.itsEndDate),
	itsResetFreq( rhs.itsResetFreq),
	itsDayCount	( rhs.itsDayCount),
	itsFwdRule	( rhs.itsFwdRule),
	itsIntRule	( rhs.itsIntRule),
	itsStubRule	( rhs.itsStubRule),
	itsResetGap	( rhs.itsResetGap),
	itsPayFreq	( rhs.itsPayFreq),
	itsPayGap	( rhs.itsPayGap),
	itsResetTiming( rhs.itsResetTiming),
	itsPayTiming( rhs.itsPayTiming),
	itsAdjFirstdate( rhs.itsAdjFirstdate),
	itsFirstDateFwdRule( rhs.itsFirstDateFwdRule),
	itsRefDate( rhs.itsRefDate),
    itsStdSpotDays  ( rhs.itsStdSpotDays),
    itsIndexFreq    ( rhs.itsIndexFreq),
    itsAccruedOrFull ( rhs.itsAccruedOrFull),
	itsFlowStartDates	( CreateClone(rhs.itsFlowStartDates) ),
    itsFlowEndDates		( CreateClone(rhs.itsFlowEndDates) ),
	itsResetDates		( CreateClone(rhs.itsResetDates) ),
	itsPaymentDates		( CreateClone(rhs.itsPaymentDates) ),	
    itsInterestDays		( CreateClone(rhs.itsInterestDays) ),		
    itsInterestTerms	( CreateClone(rhs.itsInterestTerms) ),	
	itsFwdRateStartDates( CreateClone(rhs.itsFwdRateStartDates) ),
	itsFwdRateEndDates	( CreateClone(rhs.itsFwdRateEndDates) ) ,
	itsIsBuiltFromOutside(rhs.itsIsBuiltFromOutside)
{
	itsPayCalendar		= new char[ strlen( rhs.itsPayCalendar ) + 1 ];
	strcpy( itsPayCalendar, rhs.itsPayCalendar );
	itsResetCalendar	= new char[ strlen( rhs.itsResetCalendar )+ 1 ];
	strcpy( itsResetCalendar, rhs.itsResetCalendar );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: CheckSizeAndSorted
///	Returns: void
///	Action : checks size and that it is sorted!
/////////////////////////////////////////////////////////////////
void ARM_DateStrip::CheckSizeAndSorted( size_t defaultSize, ARM_GP_Vector* vec, const char* objName )
{
	CheckSize( defaultSize, vec, objName );
	/// test that it is sorted!
	if( *vec != vec->sort() )
	{
		char msg[255];
		sprintf( msg, "vector %s is not sorted, please advise!", objName );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: CheckSize
///	Returns: void
///	Action : checks size is equal to default size
/////////////////////////////////////////////////////////////////
void ARM_DateStrip::CheckSize( size_t defaultSize, ARM_GP_Vector* vec, const char* objName )
{
	if( defaultSize != vec->size() )
	{
		char msg[255];
		sprintf( msg, "expected size of %s to be %o but found %o, please advise!",
				objName, defaultSize, vec->size() );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: ARM_DateSrip
///	Returns: void
///	Action : constructor with the already given ARM_GP_Vector
///				clones everything for more safety!
/////////////////////////////////////////////////////////////////

ARM_DateStrip::ARM_DateStrip( 
	ARM_GP_Vector* FlowStartDates,			/// Flow start dates and Flow 
										/// end dates are used to compute the
	ARM_GP_Vector* FlowEndDates,			/// period of interest 
	ARM_GP_Vector* FwdRateStartDates,		/// fwd start dates
	ARM_GP_Vector* FwdRateEndDates,		/// fwd end dates
	ARM_GP_Vector* ResetDates,				/// resetDates
	ARM_GP_Vector* PaymentDates,			/// paymentDates
	ARM_GP_Vector* InterestDays,			/// numbers of days between 2 periods
	ARM_GP_Vector* InterestTerms )			/// interest term... conversion of InterestDays 
:
	itsPayCalendar( NULL ),
	itsResetCalendar( NULL ),
	itsFlowStartDates( NULL ),
    itsFlowEndDates( NULL  ),
	itsFwdRateStartDates( NULL ),
    itsFwdRateEndDates( NULL ),
	itsResetDates( NULL  ),
	itsPaymentDates( NULL  ),
    itsInterestDays( NULL  ),
	itsInterestTerms( NULL ),
	itsIsBuiltFromOutside( true )		/// says that it is built from input and not from outside
{
	/// the validation is to check that data are in increasing order with same size!
	size_t size = FlowStartDates->size();
	CheckSizeAndSorted( size, FlowEndDates,		"FlowEndDates"		);
	if (FwdRateStartDates)
		CheckSizeAndSorted( size, FwdRateStartDates,"FwdRateStartDates"	);
	if (FwdRateEndDates)
		CheckSizeAndSorted( size, FwdRateEndDates,	"FwdRateEndDates"	);
	CheckSizeAndSorted( size, ResetDates,		"ResetDates"		);
	CheckSizeAndSorted( size, PaymentDates,		"ResetDates"		);
	CheckSize(			size, InterestDays,		"InterestDays"		);
	CheckSize(			size, InterestTerms,	"InterestTerms"		);

    /// For compatibility
	itsResetCalendar= GetCalendar( GETDEFAULTVALUESTR );
	itsPayCalendar	= GetCalendar( GETDEFAULTVALUESTR );

	/// passed all validation test	
	itsFlowStartDates	= (ARM_GP_Vector*) FlowStartDates->Clone();
    itsFlowEndDates		= (ARM_GP_Vector*) FlowEndDates->Clone();
	if (FwdRateStartDates)
		itsFwdRateStartDates= (ARM_GP_Vector*) FwdRateStartDates->Clone();
	if (FwdRateEndDates)
		itsFwdRateEndDates	= (ARM_GP_Vector*) FwdRateEndDates->Clone();
	itsResetDates		= (ARM_GP_Vector*) ResetDates->Clone();
	itsPaymentDates		= (ARM_GP_Vector*) PaymentDates->Clone();
    itsInterestDays		= (ARM_GP_Vector*) InterestDays->Clone();
	itsInterestTerms	= (ARM_GP_Vector*) InterestTerms->Clone();
	Init();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: ARM_DateSrip
///	Returns: void
///	Action : default constructor
/////////////////////////////////////////////////////////////////

ARM_DateStrip::ARM_DateStrip() 
:
	itsPayCalendar( NULL ),
	itsResetCalendar( NULL ),
	itsFlowStartDates( NULL ),
    itsFlowEndDates( NULL  ),
	itsFwdRateStartDates( NULL ),
    itsFwdRateEndDates( NULL ),
	itsResetDates( NULL  ),
	itsPaymentDates( NULL  ),
    itsInterestDays( NULL  ),
	itsInterestTerms( NULL ),
	itsIsBuiltFromOutside( true )		/// says that it is built from input and not from outside
{}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: ~ARM_DateStrip
///	Action : Destructor
/////////////////////////////////////////////////////////////////
ARM_DateStrip::~ARM_DateStrip()
{
	delete itsPayCalendar;
	itsPayCalendar = NULL;
	delete itsResetCalendar;
	itsResetCalendar = NULL;
	
	delete itsFlowStartDates;
	itsFlowStartDates = NULL;
    delete itsFlowEndDates;
	itsFlowEndDates = NULL;

	delete itsFwdRateStartDates;
	itsFwdRateStartDates = NULL;
    delete itsFwdRateEndDates;
	itsFwdRateEndDates = NULL;

	delete itsResetDates;
	itsResetDates = NULL;
	delete itsPaymentDates;
	itsPaymentDates = NULL;

    delete itsInterestDays;
	itsInterestDays = NULL;
    delete itsInterestTerms;
	itsInterestTerms = NULL;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: GenerateDatesFromTimingAndGap
///	Returns: void
///	Action : depending on the timing ARREARS and ADVANCE
///				and the gap, generate a strip of dates based
///				on itsflowStartDates or itsFlowEndDates
///				with the gap
/////////////////////////////////////////////////////////////////

ARM_GP_Vector* ARM_DateStrip::GenerateDatesFromTimingAndGap(
	int timing,
	int gap,
	const char* calendar,
	ARM_GP_Vector* FlowStartDates,
	ARM_GP_Vector* FlowEndDates ) const
{
	ARM_GP_Vector* flowDates;

    switch( timing ) 
    {
        case K_ADVANCE:
			flowDates= FlowStartDates;
	        break;

        case K_ARREARS:
			flowDates= FlowEndDates;
	        break;

		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			   "reset Timing invalid");
	}

	int nFlows = itsFlowStartDates->size();
	ARM_GP_Vector* result = new ARM_GP_Vector( nFlows );
	int i;

    for (i = 0; i < nFlows; ++i) 
    {
        ARM_Date tempDate = ARM_Date(flowDates->Elt(i));
        tempDate.GapBusinessDay(gap, const_cast<char*>(calendar) ); /// UGLY but forced to const cast
        result->Elt(i) = tempDate.GetJulian();
    }

	return result;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: GetCalendar
///	Returns: char*
///	Action : find the proper calendar
/////////////////////////////////////////////////////////////////

char* ARM_DateStrip::GetCalendar( const char* calendar )
{
	char*	defaultCalendar  = ARM_DEFAULT_CURRENCY->GetCcyName();

	char* tempCalendar = new char[ strlen( calendar ) + 1 ];
	if( !strcmp( calendar, "NULL" ) ||
		!strcmp( calendar, "DEFAULT" ) )
		strcpy( tempCalendar, defaultCalendar );
	else
		strcpy( tempCalendar, calendar );

	return tempCalendar;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: ARM_DateStrip
///	Returns: constructor
///	Action : initialises correctly all the ARM_GP_Vector member objects
///				on purpose we do not allow to have array of ARM_DateStrip
/////////////////////////////////////////////////////////////////

ARM_DateStrip::ARM_DateStrip(
		const ARM_Date& startDate,
		const ARM_Date& endDate,
		int resetFreq,						/// frequency of the strip
		int dayCount,						/// dayCount method used for computing the accrued
		const char* resetCalendar,			/// calendar used for reset
		int fwdRule,						/// whether fwds are with adjusted dates
		int intRule,						/// whether interest are K_ADJUSTED
		int stubRule,						/// ability to have K_SHORTSTART etc
		int resetGap,						/// reset gap
		int payFreq,						/// frequency of the strip
		int payGap,							/// pay gap
		const char* payCalendar,			/// calendar used for payment
		int resetTiming,					/// whether reset are in arrears or in advance
		int payTiming,						/// whether payment are in arrears or in advance
		int adjFirstdate,					/// adjust the first date to business day
		const char* refDateChar,			//// reference Date
        int stdSpotDays,                    /// standard spot days
        int indexFreq ,                     /// term of the forward
        int accruedOrFull,                  /// Accrud mode or full mode to generat FwdDates 
		int firstDateFwdRule                    
)
:
	itsStartDate( startDate ),			/// start date 
	itsEndDate( endDate ),				/// end date
	itsResetFreq( resetFreq ),			/// frequency of the strip
	itsDayCount( dayCount ),			/// dayCount itsmethod used for computing the accrued
	itsResetCalendar( NULL ),			/// calendar used for reset
	itsFwdRule(	fwdRule ),				/// whether fwds are with adjusted dates
	itsIntRule(	intRule ),				/// whether interest are K_ADJUSTED
	itsStubRule( stubRule ),			// ability to have K_SHORTSTARt itsetc

	itsResetGap( resetGap == GETDEFAULTVALUE? -ARM_DEFAULT_CURRENCY->GetSpotDays(): resetGap ),
										/// reset gap
	itsPayFreq(	payFreq  == K_DEF_FREQ?			resetFreq: payFreq ),
										/// frequency of the strip
	itsPayGap( payGap   == GETDEFAULTVALUE? 0: payGap ),				
										/// pay gap
	itsPayCalendar( NULL ),				/// calendar used for payment
	itsResetTiming( resetTiming ),		/// whether reset are in arrears or in advance
	itsPayTiming( payTiming ),			/// whether payment itsare in arrears or in advance
	itsAdjFirstdate( adjFirstdate ),	/// adjust itsthe first itsdate to business day
	itsFirstDateFwdRule( firstDateFwdRule == 9999 ? itsFwdRule : firstDateFwdRule ),

	/// initialise to NULL to avoid deleting incorrectly pointor
	itsFlowStartDates( NULL ),
    itsFlowEndDates( NULL ),
    itsFwdRateStartDates( NULL ),
    itsFwdRateEndDates( NULL ),
	itsResetDates( NULL ),
	itsPaymentDates( NULL ),
    itsInterestDays( NULL ),
    itsInterestTerms( NULL ),
	itsIsBuiltFromOutside( false ),		/// says that it is built from input and not from outside

	itsStdSpotDays( stdSpotDays == GETDEFAULTVALUE ? ARM_DEFAULT_CURRENCY->GetSpotDays(): stdSpotDays ),
										/// standard spot days

	itsIndexFreq( indexFreq == GETDEFAULTVALUE ? resetFreq : indexFreq ), /// term of the forward
    itsAccruedOrFull(accruedOrFull)

{
    /// common part (mainly for view)

    /// get the resetCalendar and payCalendar
	/// if "NULL" takes the one given by the currency of the DEFAULT_CURRENCY
	/// of the country
	itsResetCalendar= GetCalendar( resetCalendar );
	itsPayCalendar	= GetCalendar( payCalendar );

    /// in the special case of startDate >= endDate
    /// returns an empty datestrip
    if( startDate >= endDate )
    {
	    itsFlowStartDates   = new ARM_GP_Vector;
        itsFlowEndDates     = new ARM_GP_Vector;
        itsFwdRateStartDates= new ARM_GP_Vector;
        itsFwdRateEndDates  = new ARM_GP_Vector;
	    itsResetDates       = new ARM_GP_Vector;
	    itsPaymentDates     = new ARM_GP_Vector;
        itsInterestDays     = new ARM_GP_Vector;
        itsInterestTerms    = new ARM_GP_Vector;
    }
    else
    {
        /// reference date if not provided is equal to startDate
	    if( strcmp( refDateChar, GETDEFAULTVALUESTR ) == 0 )
		    itsRefDate = startDate;
	    else
		    itsRefDate = ARM_Date( const_cast<char*>(refDateChar) );

	    /// validate is first to be exception safe:
	    /// to avoid memory leak
	    Validate();

	    //// nested class for exception safety
	    struct Hold
	    {
		    Hold( char* a1, char* a2 )
			    :	arg1( a1 ), arg2(a2 ), released( false ){}
		    ~Hold( ){
			    if( !released ){
				    delete[] arg1;
				    delete[] arg2;
			    }
		    }
		    void release() { released = true; }
	    private:
		    bool released;
		    char* arg1;
		    char* arg2;
	    };

	    Hold KeepData( itsResetCalendar, itsPayCalendar );
	    CptCashFlowDates();
	    KeepData.release();
    }

    Init();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: Validate
///	Returns: void
///	Action : Validation for the DateStrip... throw exception if necessary
/////////////////////////////////////////////////////////////////

void ARM_DateStrip::Validate()
{
	/// validation
	/// empty for the time being
	if( itsStartDate > itsEndDate )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		   string( "start date " ) + itsStartDate.toString() + string( " after end date " )
		   + itsEndDate.toString() );

	/// first test the case of zero coupon
	/// zero coupon allowed if 
	///		payFreq = K_ZEROCOUPON  
	///		or if resetFreq= K_ZEROCOUPON  and pay = K_ZEROCOUPON  
	if( itsResetFreq == K_ZEROCOUPON  && itsPayFreq != K_ZEROCOUPON )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		   "Reset Frequency is zero coupon while pay is not! only both zero coupon reset and pay frequency allowed!" );

	/// test different pay and reset frequency
	if( itsPayFreq == K_ZEROCOUPON )
	{
		if( itsPayTiming != K_ARREARS )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "payment zero coupon only allowed for payment timing in arrears!" );
	}
	else
	{
		// Non sense for daily case
		if (itsPayFreq != K_DAILY)
		{
			if( itsPayFreq > itsResetFreq )
				if( ( itsPayFreq % itsResetFreq ) != 0 )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					   "Pay Frequency is not a multiple of reset frequency, the date strip cannot handle this, use two date strips instead!" );
		}

		// Non sense for daily case
		if (itsResetFreq != K_DAILY)
		{
			if( itsResetFreq > itsPayFreq )
				if( ( itsResetFreq  % itsPayFreq) != 0 )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					   "Reset Frequency is not a multiple of pay frequency, the date strip cannot handle this, use two date strips instead!" );
		}
	}
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: ARM_GP_Vector* 
///	Returns: void
///	Action : Function to adjust for the first date
/////////////////////////////////////////////////////////////////

ARM_GP_Vector* ARM_DateStrip::AdjForRefDate( ARM_GP_Vector* FlowDates, char* calendar )
{
    ARM_Date AdjStartDate = itsStartDate;
	ARM_Date adjRefDate   = ARM_Date(itsRefDate);

	int theRule = ((itsAdjFirstdate == K_ADJUSTED) && (itsIntRule == K_ADJUSTED)) ? itsFirstDateFwdRule : K_UNADJUSTED;
    AdjStartDate.GoodBusinessDay( theRule, calendar );
	adjRefDate.GoodBusinessDay(theRule, itsPayCalendar);

	 // Adjust the Reference date
	FlowDates->Elt(0) = adjRefDate.GetJulian();
	if( AdjStartDate.GetJulian() != FlowDates->Elt(0) ){
		FlowDates->push_back(AdjStartDate.GetJulian() );
		FlowDates->sort();
		CC_NS( std, unique )(FlowDates->begin(), FlowDates->end(),EqualWithSevenDays());
	}

	return FlowDates;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: CptCashFlowDates
///	Returns: void
///	Action : compute Dates
/////////////////////////////////////////////////////////////////

void ARM_DateStrip::CptCashFlowDates()
{
	/// first compute crudely itsFlowStartDates
	/// UGLY but forced to const cast
	int usedPayFreq = itsPayFreq == K_ZEROCOUPON? itsResetFreq : itsPayFreq;

	ARM_GP_Vector* tmpFlowStartDates = CptStripOfStartDates( 
		const_cast<ARM_Date&>(itsRefDate), 
		const_cast<ARM_Date&>(itsEndDate),
		usedPayFreq, itsFwdRule, itsStubRule, itsIntRule,
		itsPayCalendar, itsAdjFirstdate );

	/// adjusting for the first date using the payCalendar
	tmpFlowStartDates = AdjForRefDate( tmpFlowStartDates, itsPayCalendar );

	/// generate dates according to the itsResetCalendar
	/// temporary variable

	ARM_GP_Vector* tmpResetFlowStartDates;

	/// tests whether the reset calendar is inflation
	/// which allows to reset the same day!

	/// Inflation calendar is the only case where calendar
	/// used to compute the reset date is not the pay calendar
	/// for any other cases, the tmpResetDates computed are using 
	///	the pay Calendar!
	if( strcmp( ARM_UPPERCASE(itsResetCalendar), "INF" ) == 0 )
	{
		tmpResetFlowStartDates = CptStripOfStartDates( 
				const_cast<ARM_Date&>(itsRefDate), 
				const_cast<ARM_Date&>(itsEndDate),
				itsResetFreq, itsFwdRule, itsStubRule, itsIntRule,
				itsResetCalendar, itsAdjFirstdate );

		tmpResetFlowStartDates = AdjForRefDate( tmpResetFlowStartDates, itsResetCalendar );
	}
	else
	{
		// here we use the payCalendar but with the reset frequency
		if(itsResetFreq != itsPayFreq)
		{
			tmpResetFlowStartDates = CptStripOfStartDates( 
					const_cast<ARM_Date&>(itsRefDate), 
					const_cast<ARM_Date&>(itsEndDate),
					itsResetFreq, itsFwdRule, itsStubRule, itsIntRule,
					itsPayCalendar, itsAdjFirstdate );

			tmpResetFlowStartDates = AdjForRefDate( tmpResetFlowStartDates, itsResetCalendar );
		}
		else
			tmpResetFlowStartDates = (ARM_GP_Vector*) tmpFlowStartDates->Clone();
	}

	/// initialise all vectors
	int nFlows			= MAX(tmpFlowStartDates->size(), tmpResetFlowStartDates->size());

	itsFlowStartDates	= new ARM_GP_Vector( nFlows );
	itsFlowEndDates		= new ARM_GP_Vector( nFlows );
	itsInterestDays		= new ARM_GP_Vector( nFlows );
	itsInterestTerms	= new ARM_GP_Vector( nFlows );


	/// temporary variable
	ARM_GP_Vector* resetFlowStartDates = new ARM_GP_Vector( nFlows );
	ARM_GP_Vector* resetFlowEndDates	= new ARM_GP_Vector( nFlows );

    /// Fill cash flow period start and end dates, reset and payment dates
	int i, j = 0;
	int flowStartDatesSize = tmpFlowStartDates->size();

	if (flowStartDatesSize == 1)
	{
		if(itsResetFreq <= itsPayFreq)
		{
			for (i=0; i<nFlows; ++i) 
			{
				itsFlowStartDates->Elt(i) = tmpFlowStartDates->Elt(0);
				itsFlowEndDates->Elt(i) = ARM_Date(itsEndDate).GoodBusinessDay( itsFwdRule*itsIntRule, itsPayCalendar).GetJulian();
				
				resetFlowStartDates->Elt(i) = tmpResetFlowStartDates->Elt(0);
				resetFlowEndDates->Elt(i) = ARM_Date(itsEndDate).GoodBusinessDay( itsFwdRule*itsIntRule, itsResetCalendar).GetJulian();
			}
		}
		else
		{
			/// New case not supported before
			for (j=0; j<nFlows; j++)
			{
				// flow start and end dates are equal to start and end dates of the period
				itsFlowStartDates->Elt(j) = tmpFlowStartDates->Elt(0);
				itsFlowEndDates->Elt(j) = ARM_Date(itsEndDate).GoodBusinessDay( itsFwdRule*itsIntRule, itsPayCalendar).GetJulian();

				// reset start and reset end dates (used to deduce reset dates)
				resetFlowStartDates->Elt(j) = tmpResetFlowStartDates->Elt(j);
				if (j != nFlows-1)
					resetFlowEndDates->Elt(j) = tmpResetFlowStartDates->Elt(j+1);
				else
					resetFlowEndDates->Elt(j) = ARM_Date(itsEndDate).GoodBusinessDay( itsFwdRule*itsIntRule, itsResetCalendar).GetJulian();

				itsInterestTerms->Elt(j)= CountYears(itsDayCount, itsFlowStartDates->Elt(j), itsFlowEndDates->Elt(j));
				itsInterestDays->Elt(j)	= DaysBetweenDates(itsDayCount, itsFlowStartDates->Elt(j), itsFlowEndDates->Elt(j));
			}
		}
	}
	else
	{
		if (itsResetFreq >= itsPayFreq)
		{
			//start and end dates
			for (i=0; i<flowStartDatesSize-2; ++i) 
			{
				//in the test, adjust resetFlowStart with PayCalendar to match dates in case of INFLATION !
				while (ARM_Date(tmpResetFlowStartDates->Elt(j)).GoodBusinessDay( itsFwdRule*itsIntRule, itsPayCalendar).GetJulian() < tmpFlowStartDates->Elt(i+1))
				{
					itsFlowStartDates->Elt(j) = tmpFlowStartDates->Elt(i);
					itsFlowEndDates->Elt(j) = tmpFlowStartDates->Elt(i+1);
					j++;
				}
				// specific treatment for WEEKLY reset:
				// for next period, check that first reset is not after start !!
				if ((itsResetFreq == K_WEEKLY) && (itsPayFreq != K_ZEROCOUPON))
				{
					if (ARM_Date(tmpResetFlowStartDates->Elt(j)).GapBusinessDay( itsResetGap, itsResetCalendar, itsFwdRule*itsIntRule).GetJulian() >= tmpFlowStartDates->Elt(i+1))
					{
						itsFlowStartDates->Elt(j-1) = tmpFlowStartDates->Elt(i+1);
						itsFlowEndDates->Elt(j-1) = tmpFlowStartDates->Elt(i+2);
					}
				}
			}
			for(; j<nFlows; ++j)
			{
				//in the test, adjust resetFlowStart with PayCalendar to match dates in case of INFLATION !
				if ( (ARM_Date(tmpResetFlowStartDates->Elt(j)).GoodBusinessDay( itsFwdRule*itsIntRule, itsPayCalendar).GetJulian() < tmpFlowStartDates->Elt(flowStartDatesSize-1))
					&& //WEEKLY case:
					! ( (itsResetFreq == K_WEEKLY) && (itsPayFreq != K_ZEROCOUPON) 
						&& ARM_Date(tmpResetFlowStartDates->Elt(j)).GapBusinessDay( itsResetGap, itsResetCalendar, itsFwdRule*itsIntRule).GetJulian() >= tmpFlowStartDates->Elt(i+1)) )
				{
					itsFlowStartDates->Elt(j) = tmpFlowStartDates->Elt(i);
					itsFlowEndDates->Elt(j) = tmpFlowStartDates->Elt(i+1);
				}
				else
				{
					itsFlowStartDates->Elt(j) = tmpFlowStartDates->Elt(i+1);
					itsFlowEndDates->Elt(j) = ARM_Date(itsEndDate).GoodBusinessDay( itsFwdRule*itsIntRule, itsPayCalendar).GetJulian();
				}
			}
			// reset start and reset end dates (used to deduce reset dates)
			for (j=0; j<nFlows; j++)
			{
				resetFlowStartDates->Elt(j) = tmpResetFlowStartDates->Elt(j);
				if (j != nFlows-1)
					resetFlowEndDates->Elt(j) = tmpResetFlowStartDates->Elt(j+1);
				else
					resetFlowEndDates->Elt(j) = ARM_Date(itsEndDate).GoodBusinessDay( itsFwdRule*itsIntRule, itsResetCalendar).GetJulian();

				itsInterestTerms->Elt(j)= CountYears(itsDayCount, itsFlowStartDates->Elt(j), itsFlowEndDates->Elt(j));
				itsInterestDays->Elt(j)	= DaysBetweenDates(itsDayCount, itsFlowStartDates->Elt(j), itsFlowEndDates->Elt(j));
			}
		}
		else
		{
			//start and end dates
			for (i=0; i<nFlows; ++i) 
			{
				itsFlowStartDates->Elt(i) = tmpFlowStartDates->Elt(i);
				if (i != nFlows-1)
					itsFlowEndDates->Elt(i) = tmpFlowStartDates->Elt(i+1);
				else
					itsFlowEndDates->Elt(i) = ARM_Date(itsEndDate).GoodBusinessDay( itsFwdRule*itsIntRule, itsPayCalendar).GetJulian();
			}
			// reset start and reset end dates (used to deduce reset dates)
			for (i=0; i<nFlows-1 && j<tmpResetFlowStartDates->size()-1; j++) 
			{
				//in the test, adjust reset with PayCalendar to match dates in case of INFLATION !
				while (tmpFlowStartDates->Elt(i) < ARM_Date(tmpResetFlowStartDates->Elt(j+1)).GoodBusinessDay( itsFwdRule*itsIntRule, itsPayCalendar).GetJulian())
				{
					resetFlowStartDates->Elt(i) = tmpResetFlowStartDates->Elt(j);
					resetFlowEndDates->Elt(i) = tmpResetFlowStartDates->Elt(j+1);
					i++;
				}
			}
			for (; i<nFlows; i++) 
			{
				resetFlowStartDates->Elt(i) = tmpResetFlowStartDates->Elt(tmpResetFlowStartDates->size()-1);
				resetFlowEndDates->Elt(i) = ARM_Date(itsEndDate).GoodBusinessDay( itsFwdRule*itsIntRule, itsResetCalendar).GetJulian();
			}
		}
	}

	for (j=0; j<nFlows; j++)
	{
		itsInterestTerms->Elt(j)= CountYears(itsDayCount, itsFlowStartDates->Elt(j), itsFlowEndDates->Elt(j));
		itsInterestDays->Elt(j)	= DaysBetweenDates(itsDayCount, itsFlowStartDates->Elt(j), itsFlowEndDates->Elt(j));
	}

	/// be aware that we do not need 
	/// a new ARM_GP_Vector( nFlows ) for itsResetDates
	/// and since it is done in the GenerateDatesFromTimingAndGap
	itsResetDates	= GenerateDatesFromTimingAndGap( itsResetTiming, itsResetGap, itsResetCalendar, resetFlowStartDates, resetFlowEndDates );
	if( itsPayFreq == K_ZEROCOUPON )
	{
		ARM_Date finalPayDate( (*itsFlowEndDates)[ itsFlowEndDates->size()-1 ] );
		finalPayDate.GapBusinessDay(itsPayGap, itsPayCalendar );
		itsPaymentDates		= new ARM_GP_Vector( itsFlowEndDates->size(), finalPayDate.GetJulian() );
	}
	else
	{
		itsPaymentDates	= GenerateDatesFromTimingAndGap( itsPayTiming, itsPayGap, itsPayCalendar, itsFlowStartDates, itsFlowEndDates );
	}

	/// generate fwdRateStart and end Dates
	/// in the particular case of zero coupon
	/// fwd start and end are the same as the flow dates
	if( itsResetFreq == K_ZEROCOUPON )
	{
		itsFwdRateStartDates = (ARM_GP_Vector*) itsFlowStartDates->Clone();
		itsFwdRateEndDates	 = (ARM_GP_Vector*) itsFlowEndDates->Clone();
	}
	else 
	/// otherwise the logic comes from the reset timing 
	/// and the fwdRateEndDates are computed adding the resetFreq
	{
		itsFwdRateStartDates = (ARM_GP_Vector*) itsResetDates->Clone();
		itsFwdRateEndDates	 = new ARM_GP_Vector(itsResetDates->size());

		if( itsResetTiming == K_ADVANCE )
		{
			if(itsAccruedOrFull == K_ACCRUED)
            {
				for(i=0; i<itsFwdRateStartDates->size(); ++i )
				{
					ARM_Date tmpStartDate(itsFwdRateStartDates->Elt(i));
					tmpStartDate.GapBusinessDay(itsStdSpotDays,itsResetCalendar);
					itsFwdRateStartDates->Elt(i) = tmpStartDate.GetJulian();
					itsFwdRateEndDates->Elt(i) = itsFwdRateStartDates->Elt(i)+itsFlowEndDates->Elt(i) 
						- itsFlowStartDates->Elt(i);
				}
            }
            else
            {
				for(i=0; i<itsFwdRateStartDates->size(); ++i )
				{
					ARM_Date tmpStartDate(itsFwdRateStartDates->Elt(i));
					tmpStartDate.GapBusinessDay(itsStdSpotDays,itsResetCalendar);
					itsFwdRateStartDates->Elt(i) = tmpStartDate.GetJulian();
			
					ARM_Date tmpEndDate(itsFwdRateStartDates->Elt(i));
					tmpEndDate.AddPeriod(itsIndexFreq,itsResetCalendar);
					tmpEndDate.GoodBusinessDay(itsFwdRule*itsIntRule,itsResetCalendar);
					itsFwdRateEndDates->Elt(i) = tmpEndDate.GetJulian();
				}
            }
		}
		else /// in arrears case
		{
			for(i=0; i<itsFwdRateStartDates->size(); ++i )
			{
				ARM_Date tmpStartDate(itsFwdRateStartDates->Elt(i));
				tmpStartDate.GapBusinessDay(itsStdSpotDays,itsResetCalendar);
				itsFwdRateStartDates->Elt(i) = tmpStartDate.GetJulian();

				ARM_Date tmpEndDate(itsFwdRateStartDates->Elt(i));
				tmpEndDate.AddPeriod(itsIndexFreq,itsResetCalendar);
				tmpEndDate.GoodBusinessDay(itsFwdRule*itsIntRule,itsResetCalendar);
				itsFwdRateEndDates->Elt(i) = tmpEndDate.GetJulian();
			}
		}
	}

	delete resetFlowStartDates;
	delete resetFlowEndDates;
	delete tmpFlowStartDates;
	delete tmpResetFlowStartDates;
}
 


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: CptCashFlowDates
///	Returns: void
///	Action : Init function. We do not initialise pointor to NULL
///				because this is unnecesary, although some people may not
///				know!
/////////////////////////////////////////////////////////////////

void ARM_DateStrip::Init()
{
	CC_ARM_SETNAME(ARM_DATESTRIP);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: GetMemberData
///	Returns: ARM_GP_Vector*
///	Action : method to access data
/////////////////////////////////////////////////////////////////

ARM_GP_Vector* ARM_DateStrip::GetMemberData( int type )
{
	switch( type )
    {
        case K_START_DATES:
			return itsFlowStartDates;

        case K_END_DATES:
			return itsFlowEndDates;

        case K_RESET_DATES:
			return itsResetDates;

        case K_PAY_DATES:
			return itsPaymentDates;

		case K_INT_DAYS:
			return itsInterestDays;

		case K_INT_TERMS:
			return itsInterestTerms;

		case K_FWD_START_DATES:
			return itsFwdRateStartDates;
		
		case K_FWD_END_DATES:
			return itsFwdRateEndDates;

        default :
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			   "unknown type");
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: CptStritOfStartDates
///	Returns: ARM_GP_Vector*
///	Action : computes the strip of Start Dates
/////////////////////////////////////////////////////////////////

ARM_GP_Vector* ARM_DateStrip::CptStripOfStartDates(
	ARM_Date& StartDate, 
	ARM_Date& EndDate, 
	int frequency,
	int fwdRule, 
	int TypeStub, 
	int intRule, 
	char* ccy,
	int adjFirstdate)
{
	/// part to handle correctly the case of zero coupon
	ARM_Date AdjStartDate;
	
	int theRule = intRule == 0? intRule : fwdRule;
	
	if ((adjFirstdate) && ( theRule != 0 ))
	{
		AdjStartDate = StartDate;
		AdjStartDate.GoodBusinessDay(theRule, ccy);
	}
	else
		AdjStartDate = StartDate;
	
	if ( frequency == K_ZEROCOUPON )
		return new ARM_GP_Vector(1, StartDate.GetJulian());
	
	/// other cases
	int fin, nbFlow = 0;
	int IsStub = 0;
	ARM_Date flowDate, AdjFlowDate, AdjEndDate;
	ARM_Date tmpDate;
	int i, SizeMax;
	
	AdjEndDate = EndDate;
	AdjEndDate.GoodBusinessDay(theRule, ccy);
	
	// Case of zero coupon : put in result the right start date
	SizeMax = int(((EndDate - StartDate)/365+1)*frequency);
	ARM_GP_Vector result(SizeMax);
	
	// if the stub must be at the begining, we decrease <period> months
	// from EndDate until StartDate
	int step = 2;
	
	if ( TypeStub == K_SHORTSTART || TypeStub == K_LONGSTART )
	{
		flowDate = EndDate;
		
		flowDate.AddPeriod(-1*frequency, ccy);
		
		AdjFlowDate = flowDate;
		
		AdjFlowDate.GoodBusinessDay(theRule, ccy);
		
		while ( AdjStartDate.GetJulian() < AdjFlowDate.GetJulian() )
		{
			result[nbFlow] = flowDate.GetJulian();
			nbFlow++;
			
			if (( frequency != K_DAILY ) && ( frequency != K_WEEKLY ) )
			{
				flowDate = EndDate;
				flowDate.AddPeriodMult(-1*frequency, step, ccy);
			}
			else
				flowDate.AddPeriod(-1*frequency, ccy);
			
			AdjFlowDate = flowDate;
			AdjFlowDate.GoodBusinessDay(theRule, ccy);
			step++;
		}

		if (AdjFlowDate.GetJulian() == AdjStartDate.GetJulian())
		{
			IsStub = 0;
			result[nbFlow] = flowDate.GetJulian();
			nbFlow++;
		}
		else 
		{
			IsStub = 1; 
			result[nbFlow] = StartDate.GetJulian();
			nbFlow++;
		}
	}
	// if the stub must be at the end, we increase <period> months
	// from StartDate until EndDate
	else
	{
		flowDate		= StartDate;
		result[nbFlow]	= flowDate.GetJulian();
		flowDate.AddPeriod(frequency, ccy);
		AdjFlowDate		= flowDate;
		AdjFlowDate.GoodBusinessDay(theRule, ccy);
		nbFlow++;
		
		while ( AdjFlowDate.GetJulian() < AdjEndDate.GetJulian())
		{
			result[nbFlow] = flowDate.GetJulian();
			nbFlow++;
			flowDate = StartDate;
			flowDate.AddPeriodMult(frequency, step, ccy);
			AdjFlowDate = flowDate;
			AdjFlowDate.GoodBusinessDay(theRule, ccy);
			step++;
		}
		
		IsStub = AdjFlowDate.GetJulian() == AdjEndDate.GetJulian() ? 0 : 1;
	}
	
	fin = nbFlow;
	
	// if TypeStub is LONG (start or end) we merge the last buckets
	if (IsStub && TypeStub == K_LONGSTART && nbFlow >= 2)
	{
		nbFlow--;
		result[nbFlow-1] = result[nbFlow];
	}
	
	if (IsStub && TypeStub == K_LONGEND && nbFlow >= 2)
		nbFlow--;
	
	///  Adjusting the dates
	ARM_GP_Vector resultAdj(nbFlow);
	
	if ( TypeStub == K_SHORTSTART || TypeStub == K_LONGSTART )
	{
		for (i = 0; i < nbFlow; i++)
		{
			tmpDate = ARM_Date(result[i]);
			resultAdj[nbFlow-(i+1)] = tmpDate.GoodBusinessDay(theRule, ccy).GetJulian();
		}
	}
	else if (TypeStub == K_SHORTEND || TypeStub == K_LONGEND)
	{
		for (i = 0; i < nbFlow; i++)
		{
			tmpDate = ARM_Date(result[i]);
			resultAdj[i] = tmpDate.GoodBusinessDay(theRule, ccy).GetJulian();
		}
	}

	//  Delete doublons
	int size = 0;
	double prevDate = 0.0;
	
	for (i = 0; i < nbFlow; i++)
	{
		if (prevDate != resultAdj[i])
		{
			prevDate = resultAdj[size] = resultAdj[i];
			size++;
		}
	}
// FIXMEFRED: mig.vc8 (22/05/2007 15:51:14): explicit cast
	return new ARM_GP_Vector(size, &(*resultAdj.begin()));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: All the set method
///	Returns: void
///	Action : Set the objects! with cloning them!
/// set directly the ARM_GP_Vector pointor with cloning it
/// it has been a design decision to clone dates .. 
/// use the noclone version if required
/////////////////////////////////////////////////////////////////

void ARM_DateStrip::SetFlowStartDates( ARM_GP_Vector* FlowStartDates )
{
	delete itsFlowStartDates; 
	itsFlowStartDates = (ARM_GP_Vector*) FlowStartDates->Clone(); 
}

void ARM_DateStrip::SetFlowEndDates( ARM_GP_Vector* FlowEndDates )
{
	delete itsFlowEndDates; 
	itsFlowEndDates = (ARM_GP_Vector*) FlowEndDates->Clone();
}

void ARM_DateStrip::SetFwdRateStartDates( ARM_GP_Vector* FwdRateStartDates )
{
	delete itsFwdRateStartDates;
	itsFwdRateStartDates = (ARM_GP_Vector*) FwdRateStartDates->Clone();
}

void ARM_DateStrip::SetFwdRateEndDates( ARM_GP_Vector* FwdRateEndDates)
{
	delete itsFwdRateEndDates;
	itsFwdRateEndDates = (ARM_GP_Vector*) FwdRateEndDates->Clone();
}

void ARM_DateStrip::SetResetDates( ARM_GP_Vector* ResetDates )
{
	delete itsResetDates;
	itsResetDates = (ARM_GP_Vector*) ResetDates->Clone();
}

void ARM_DateStrip::SetPaymentDates( ARM_GP_Vector* PaymentDates )
{
	delete itsPaymentDates;
	itsPaymentDates = (ARM_GP_Vector*) PaymentDates->Clone();
}

void ARM_DateStrip::SetInterestDays( ARM_GP_Vector* InterestDays )
{
	delete itsInterestDays; 
	itsInterestDays  = (ARM_GP_Vector*) InterestDays->Clone(); 
}

void ARM_DateStrip::SetInterestTerms( ARM_GP_Vector* InterestTerms ) 
{
	delete itsInterestTerms;
	itsInterestTerms = (ARM_GP_Vector*) InterestTerms->Clone(); 
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: All the set method
///	Returns: void
///	Action : Set the objects! without cloning them!
/////////////////////////////////////////////////////////////////

void ARM_DateStrip::SetFlowStartDatesNoClone( ARM_GP_Vector* FlowStartDates ) 
{
	delete itsFlowStartDates;
	itsFlowStartDates =  FlowStartDates; 
}

void ARM_DateStrip::SetFlowEndDatesNoClone( ARM_GP_Vector* FlowEndDates ) 
{
	delete itsFlowEndDates;
	itsFlowEndDates =  FlowEndDates; 
}

void ARM_DateStrip::SetFwdRateStartDatesNoClone( ARM_GP_Vector* FwdRateStartDates ) 
{
	delete itsFwdRateStartDates;
	itsFwdRateStartDates = FwdRateStartDates; 
}

void ARM_DateStrip::SetFwdRateEndDatesNoClone( ARM_GP_Vector* FwdRateEndDates) 
{
	delete itsFwdRateEndDates; 
	itsFwdRateEndDates = FwdRateEndDates; 
}

void ARM_DateStrip::SetResetDatesNoClone( ARM_GP_Vector* ResetDates ) 
{
	delete itsResetDates; 
	itsResetDates =  ResetDates; 
}

void ARM_DateStrip::SetPaymentDatesNoClone( ARM_GP_Vector* PaymentDates ) 
{
	delete itsPaymentDates; 
	itsPaymentDates =  PaymentDates;
}

void ARM_DateStrip::SetInterestDaysNoClone( ARM_GP_Vector* InterestDays ) 
{
	delete itsInterestDays;
	itsInterestDays  =  InterestDays; 
}

void ARM_DateStrip::SetInterestTermsNoClone( ARM_GP_Vector* InterestTerms ) 
{
	delete itsInterestTerms;
	itsInterestTerms =  InterestTerms; 
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: InsertDate
///	Returns: void
///	Action : insert a new date in the datestrip
/////////////////////////////////////////////////////////////////
void ARM_DateStrip::InsertDate(
		int idx,
		double flowStartDate,
		double flowEndDate,
		double fwdRateStartDate,
		double fwdRateEndDate,
		double resetDate,
		double paymentDate,
		double interestDays,
		double interestTerm)
{
	itsFlowStartDates->insert(itsFlowStartDates->begin()+idx,flowStartDate);
	itsFlowEndDates->insert(itsFlowEndDates->begin()+idx,flowEndDate);
	itsFwdRateStartDates->insert(itsFwdRateStartDates->begin()+idx,fwdRateStartDate);
	itsFwdRateEndDates->insert(itsFwdRateEndDates->begin()+idx,fwdRateEndDate);
	itsResetDates->insert(itsResetDates->begin()+idx,resetDate);
	itsPaymentDates->insert(itsPaymentDates->begin()+idx,paymentDate);
	itsInterestDays->insert(itsInterestDays->begin()+idx,interestDays);
	itsInterestTerms->insert(itsInterestTerms->begin()+idx,interestTerm);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: InsertDate
///	Returns: void
///	Action : insert a new date in the datestrip
/////////////////////////////////////////////////////////////////
void ARM_DateStrip::fill(const ARM_DateStrip& datestrip)
{
	itsFlowStartDates->fill(*datestrip.GetFlowStartDates());
	itsFlowEndDates->fill(*datestrip.GetFlowEndDates());
	itsFwdRateStartDates->fill(*datestrip.GetFwdRateStartDates());
	itsFwdRateEndDates->fill(*datestrip.GetFwdRateEndDates());
	itsResetDates->fill(*datestrip.GetResetDates());
	itsPaymentDates->fill(*datestrip.GetPaymentDates());
	itsInterestDays->fill(*datestrip.GetInterestDays());
	itsInterestTerms->fill(*datestrip.GetInterestTerms());
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: ResizeAndBuilt
///	Returns: void
///	Action : Rezise and rebuilt a new object
/////////////////////////////////////////////////////////////////
// FIXMEFRED: mig.vc8 (22/05/2007 16:01:06): a lot of missing explicit casts
void ARM_DateStrip::ResizeAndBuilt( size_t begin, size_t end)
{
	ARM_GP_Vector FlowStartDates(&(*itsFlowStartDates->begin())+begin, &(*itsFlowStartDates->begin())+end);
	delete itsFlowStartDates;
	itsFlowStartDates = (ARM_GP_Vector*)FlowStartDates.Clone();

	ARM_GP_Vector FlowEndDates(&(*itsFlowEndDates->begin())+begin, &(*itsFlowEndDates->begin())+end);
	delete itsFlowEndDates;
	itsFlowEndDates = (ARM_GP_Vector*)FlowEndDates.Clone();

	ARM_GP_Vector FwdRateStartDates(&(*itsFwdRateStartDates->begin())+begin, &(*itsFwdRateStartDates->begin())+end);
	delete itsFwdRateStartDates;
	itsFwdRateStartDates = (ARM_GP_Vector*)FwdRateStartDates.Clone();

	ARM_GP_Vector FwdRateEndDates(&(*itsFwdRateEndDates->begin())+begin, &(*itsFwdRateEndDates->begin())+end);
	delete itsFwdRateEndDates;
	itsFwdRateEndDates = (ARM_GP_Vector*)FwdRateEndDates.Clone();

	ARM_GP_Vector ResetDates(&(*itsResetDates->begin())+begin, &(*itsResetDates->begin())+end);
	delete itsResetDates;
	itsResetDates = (ARM_GP_Vector*)ResetDates.Clone();

	ARM_GP_Vector PaymentDates(&(*itsPaymentDates->begin())+begin, &(*itsPaymentDates->begin())+end);
	delete itsPaymentDates;
	itsPaymentDates = (ARM_GP_Vector*)PaymentDates.Clone();

	ARM_GP_Vector InterestDays(&(*itsInterestDays->begin())+begin, &(*itsInterestDays->begin())+end);
	delete itsInterestDays;
	itsInterestDays = (ARM_GP_Vector*)InterestDays.Clone();

	ARM_GP_Vector InterestTerms(&(*itsInterestTerms->begin())+begin, &(*itsInterestTerms->begin())+end);
	delete itsInterestTerms;
	itsInterestTerms = (ARM_GP_Vector*)InterestTerms.Clone();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_DateStrip
///	Routine: toString
///	Returns: string
///	Action : toString method for the date strip
/////////////////////////////////////////////////////////////////

string ARM_DateStrip::toString( const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << "\n Date Strip Object:\n";

	size_t i,nFlows=100000;
    nFlows = (itsFlowStartDates     ? CC_Min(nFlows,itsFlowStartDates->size()) : nFlows);
    nFlows = (itsFlowEndDates       ? CC_Min(nFlows,itsFlowEndDates->size()) : nFlows);
    nFlows = (itsResetDates         ? CC_Min(nFlows,itsResetDates->size()) : nFlows);
    nFlows = (itsPaymentDates       ? CC_Min(nFlows,itsPaymentDates->size()) : nFlows);
    nFlows = (itsFwdRateStartDates  ? CC_Min(nFlows,itsFwdRateStartDates->size()) : nFlows);
    nFlows = (itsFwdRateEndDates    ? CC_Min(nFlows,itsFwdRateEndDates->size()) : nFlows);
    nFlows = (itsInterestDays       ? CC_Min(nFlows,itsInterestDays->size()) : nFlows);
    nFlows = (itsInterestTerms      ? CC_Min(nFlows,itsInterestTerms->size()) : nFlows);
    if(nFlows == 100000)
        nFlows=0;

    double notional = 0.0, cfValue = 0.0;
	
	if( itsIsBuiltFromOutside )
		os << "\n\nCreated with given vector.. no conventions\n\n\n";
	else
	{
		 os << "\n\nDetails\n\n";
		 os << "Start Date		: " <<itsStartDate.toString() << "\n";
		 os << "End Date		: " << itsEndDate.toString() << "\n";
		 os << "Ref Date		: " << itsRefDate.toString() << "\n";
		 os << "Reset Freq		: " << ARM_ArgConvReverse_LgNameFrequency.GetString( itsResetFreq ) << "\n";
		 os << "Day Count		: " << ARM_ArgConvReverse_LgNameDayCount.GetString( itsDayCount ) << "\n";
		 os << "Reset Calendar		: " << itsResetCalendar << "\n";
		 os << "Fwd Rule		: "	<< ARM_ArgConvReverse_FwdRules.GetString(itsFwdRule) << "\n";
		 os << "Int Rule		: "	<< ARM_ArgConvReverse_InterestRules.GetString(itsIntRule)  << "\n";
		 os << "Stub Rule		: "	<< ARM_ArgConvReverse_StubRules.GetString( itsStubRule) << "\n";
		 os << "Reset Gap		: " << itsResetGap << "days\n";
		 os << "Pay Freq		: "		<< ARM_ArgConvReverse_LgNameFrequency.GetString( itsPayFreq ) << "\n";
		 os << "Pay Gap  		: " << itsPayGap << "days\n";
		 os << "PayCalendar		: " << itsPayCalendar << "\n";
		 os << "Reset Timing		: " << ARM_ArgConvReverse_LgTimingMod.GetString( itsResetTiming) << "\n";
		 os << "Pay Timing		: " << ARM_ArgConvReverse_LgTimingMod.GetString( itsPayTiming) << "\n";
		
		os << "Adjust first date	:";
		if( itsAdjFirstdate )
			os << " true (not relevant if IntRule = UNADJUSTED)\n";
		else
			os << " false\n";
	}
	
    os << "\n\nStart Dates\t End Dates\t Fixing Dates\t Payment Dates\t FwdStart Dates\t FwdEnd Dates\t Interest Days\t Interest Terms\n";

    bool isFlowStartDates,isFlowEndDates,isResetDates,isPaymentDates,isFwdRateStartDates,isFwdRateEndDates;
    for (i = 0; i < nFlows; ++i )
    {

        isFlowStartDates=(itsFlowStartDates        && (*itsFlowStartDates)[i] != ARM_DateStrip::DateStripCombiner_BlankData);
        isFlowEndDates=(itsFlowEndDates            && (*itsFlowEndDates)[i] != ARM_DateStrip::DateStripCombiner_BlankData);
        isResetDates=(itsResetDates                && (*itsResetDates)[i] != ARM_DateStrip::DateStripCombiner_BlankData);
        isPaymentDates=(itsPaymentDates            && (*itsPaymentDates)[i] != ARM_DateStrip::DateStripCombiner_BlankData);
        isFwdRateStartDates=(itsFwdRateStartDates  && (*itsFwdRateStartDates)[i] != ARM_DateStrip::DateStripCombiner_BlankData);
        isFwdRateEndDates=(itsFwdRateEndDates      && (*itsFwdRateEndDates)[i] != ARM_DateStrip::DateStripCombiner_BlankData);

        if( isFlowStartDates || isFlowEndDates || isResetDates || isPaymentDates || isFwdRateStartDates ||isFwdRateEndDates)
		{
            os  << (isFlowStartDates        ? ARM_Date((*itsFlowStartDates)[i]).toString()      : "          ") << "\t " 
				<< (isFlowEndDates          ? ARM_Date((*itsFlowEndDates)[i]).toString()        : "          ") << "\t " 
				<< (isResetDates            ? ARM_Date((*itsResetDates)[i]).toString()          : "          ") << "\t " 
				<< (isPaymentDates          ? ARM_Date((*itsPaymentDates)[i]).toString()        : "          ") << "\t " 
				<< (isFwdRateStartDates     ? ARM_Date((*itsFwdRateStartDates)[i]).toString()   : "          ") << "\t " 
				<< (isFwdRateEndDates       ? ARM_Date((*itsFwdRateEndDates)[i]).toString()     : "          ") << "\t "  
				<< fixed << setw(15) << setprecision(10)
				<< (itsInterestDays         ? (*itsInterestDays)[i]                             : 0.0) << "\t " 
				<< (itsInterestTerms        ? (*itsInterestTerms)[i]                            : 0.0) << "\n";
		}
	}
	return os.str();
}


CC_END_NAMESPACE()

/*----------------------------------------------------------------------*/
/*---- End of file ----*/
