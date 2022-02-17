
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/vanillacorridor.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/typedef.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/pricingmodel.h"
#include "gpbase/port.h"
#include "gpinfra/typedef.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/dealdescription.h"
#include "gpbase/env.h"
#include "gpbase/gpvector.h"
#include "gpbase/checkarg.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/countedptr.h"
#include "gpbase/autocleaner.h"

/// gpnummethods
#include "gpnummethods/treemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )

const int TREE_NBSTEPS_PER_YEAR	=30;
const double STD_DEV_RATIO		=5.0;
const double MIN_STD_DEV		=0.001;


const double ARM_VanillaSpreadOptionArg::PayDateToleranceInDays = 3.0;

const string ARM_VanillaSpreadOptionArg::SpreadOptionColNamesTable [] =
{
    "EventDate",	
   	"PayDate",
    "RateLongStartDate",
	"RateLongEndDate",
	"RateShortStartDate",
	"RateShortEndDate",
	"RateLong",
    "RateShort",	
	"IT",	   
    "Notional",
	"Discount",	 
	"CoeffLong",
    "CoeffShort",
	"Strike",	
    "SpreadOption"
};

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSpreadOptionArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaSpreadOptionArg::CopyNoCleanUp(const ARM_VanillaSpreadOptionArg& rhs)
{
    itsStartTime        =   rhs.itsStartTime;
    itsEndTime          =   rhs.itsEndTime;

#if defined(__GP_STRICT_VALIDATION)
	ThrowErrorOnNullObject( "ResetTimes",						rhs.itsResetTimes );
	ThrowErrorOnNullObject( "PayTimes",							rhs.itsPayTimes );
	ThrowErrorOnNullObject( "PayPeriods",						rhs.itsPayPeriods );
	ThrowErrorOnNullObject( "Notional",							rhs.itsNotional );
	ThrowErrorOnNullObject( "CoeffLong",						rhs.itsCoeffLong );
	ThrowErrorOnNullObject( "CoeffShort",						rhs.itsCoeffShort );
	ThrowErrorOnNullObject( "Strikes",							rhs.itsStrikes );
	ThrowErrorOnNullObject( "swapLongFloatStartTime",	        rhs.itsSwapLongFloatStartTime );
	ThrowErrorOnNullObject( "swapLongFloatEndTime",				rhs.itsSwapLongFloatEndTime );
	/*ThrowErrorOnEmptyVector( "swapLongFixPayTimes",				rhs.itsSwapLongFixPayTimes );
	ThrowErrorOnEmptyVector( "swapLongFixPayPeriods",	        rhs.itsSwapLongFixPayPeriods );*/

	ThrowErrorOnNullObject( "swapShortFloatStartTime",	        rhs.itsSwapShortFloatStartTime );
	ThrowErrorOnNullObject( "swapShortFloatEndTime",	        rhs.itsSwapShortFloatEndTime );
	/*ThrowErrorOnEmptyVector( "swapShortFixPayTimes",			rhs.itsSwapShortFixPayTimes );
	ThrowErrorOnEmptyVector( "swapShortFixPayPeriods",	        rhs.itsSwapShortFixPayPeriods );*/

//	ThrowErrorOnNullObject( "swapPayFloatStartTime",	        rhs.itsSwapPayFloatStartTime );
//	ThrowErrorOnNullObject( "swapPayFloatEndTime",				rhs.itsSwapPayFloatEndTime );
	/*ThrowErrorOnEmptyVector( "swapPayFixPayTimes",				rhs.itsSwapPayFixPayTimes );
	ThrowErrorOnEmptyVector( "swapPayFixPayPeriods",	        rhs.itsSwapPayFixPayPeriods );*/

//	ThrowErrorOnNullObject( "periodIndex",						rhs.itsPeriodIndex );
//	ThrowErrorOnNullObject( "payIndexResetTime",				rhs.itsPayIndexResetTimes );
#endif

	/// for fast access no test of NULL pointor in release!	
	itsResetTimes						=  static_cast<ARM_GP_Vector*>( rhs.itsResetTimes->Clone());
	itsPayTimes							=  static_cast<ARM_GP_Vector*>( rhs.itsPayTimes->Clone());
	itsPayPeriods						=  static_cast<ARM_GP_Vector*>( rhs.itsPayPeriods->Clone());
	itsNotional							=  static_cast<ARM_GP_Vector*>( rhs.itsNotional->Clone());
	itsCoeffLong						=  static_cast<ARM_GP_Vector*>( rhs.itsCoeffLong->Clone());
	itsCoeffShort						=  static_cast<ARM_GP_Vector*>( rhs.itsCoeffShort->Clone());
	itsStrikes							=  static_cast<ARM_GP_Vector*>( rhs.itsStrikes->Clone());

	itsSwapLongFloatStartTime			=  static_cast<ARM_GP_Vector*>( rhs.itsSwapLongFloatStartTime->Clone());
	itsSwapLongFloatEndTime				=  static_cast<ARM_GP_Vector*>( rhs.itsSwapLongFloatEndTime->Clone());
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsSwapLongFixPayTimes,itsSwapLongFixPayTimes);
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsSwapLongFixPayPeriods,itsSwapLongFixPayPeriods);

	itsSwapShortFloatStartTime			=  static_cast<ARM_GP_Vector*>( rhs.itsSwapShortFloatStartTime->Clone());
	itsSwapShortFloatEndTime			=  static_cast<ARM_GP_Vector*>( rhs.itsSwapShortFloatEndTime->Clone());
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsSwapShortFixPayTimes,itsSwapShortFixPayTimes);
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsSwapShortFixPayPeriods,itsSwapShortFixPayPeriods);

	itsSwapPayFloatStartTime			=  rhs.itsSwapPayFloatStartTime ? static_cast<ARM_GP_Vector*>( rhs.itsSwapPayFloatStartTime->Clone()) : NULL;
	itsSwapPayFloatEndTime				=  rhs.itsSwapPayFloatEndTime ? static_cast<ARM_GP_Vector*>( rhs.itsSwapPayFloatEndTime->Clone()) : NULL;
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsSwapPayFixPayTimes,itsSwapPayFixPayTimes);
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsSwapPayFixPayPeriods,itsSwapPayFixPayPeriods);

	itsIsDigital = rhs.itsIsDigital;
	itsSpread1 = rhs.itsSpread1;
	itsSpread2 = rhs.itsSpread2;
	itsFixValues 						= rhs.itsFixValues ? static_cast<ARM_GP_Vector*>( rhs.itsFixValues->Clone()) : NULL;
	itsPayIndexLeverages				= rhs.itsPayIndexLeverages ? static_cast<ARM_GP_Vector*>( rhs.itsPayIndexLeverages->Clone()) : NULL;

	itsShortIndexIsLibor	= rhs.itsShortIndexIsLibor;
	itsLongIndexIsLibor		= rhs.itsLongIndexIsLibor;
	itsPayIndexType			= rhs.itsPayIndexType;

	itsPeriodIndex			= rhs.itsPeriodIndex ? static_cast<ARM_IntVector*>(rhs.itsPeriodIndex->Clone()) : NULL;
	itsPayIndexResetTimes	= rhs.itsPayIndexResetTimes ? static_cast<ARM_GP_Vector*>(rhs.itsPayIndexResetTimes->Clone()) : NULL;
}


ARM_VanillaSpreadOptionArg::ARM_VanillaSpreadOptionArg(const ARM_VanillaSpreadOptionArg& arg)
:	ARM_VanillaArg(arg),
	itsResetTimes(NULL),
	itsPayTimes(NULL),
	itsPayPeriods(NULL),
	itsNotional(NULL),
	itsCoeffLong(NULL),
	itsCoeffShort(NULL),
	itsStrikes(NULL),
	itsSwapLongFloatStartTime(NULL),
	itsSwapLongFloatEndTime(NULL),
	itsSwapShortFloatStartTime(NULL),
	itsSwapShortFloatEndTime(NULL),
	itsSwapPayFloatStartTime(NULL),
	itsSwapPayFloatEndTime(NULL),
	itsSwapLongFixPayTimes(0),
	itsSwapLongFixPayPeriods(0),
	itsSwapShortFixPayTimes(0),
	itsSwapShortFixPayPeriods(0),
	itsSwapPayFixPayTimes(0),
	itsSwapPayFixPayPeriods(0),
	itsIsDigital(false),
	itsSpread1(0.0),
	itsSpread2(0.0),
	itsFixValues(NULL),
	itsPayIndexLeverages(NULL),
	itsShortIndexIsLibor (false),
	itsLongIndexIsLibor (false),
	itsPayIndexType (K_FIXED)
{
    CopyNoCleanUp(arg);
}

ARM_VanillaSpreadOptionArg& ARM_VanillaSpreadOptionArg::operator=(const ARM_VanillaSpreadOptionArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
	 	CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}


void ARM_VanillaSpreadOptionArg::CleanUp()
{
	delete itsResetTimes;
	delete itsPayTimes;
	delete itsPayPeriods;
	delete itsNotional;
	delete itsCoeffLong;
	delete itsCoeffShort;
	delete itsStrikes;
	delete itsSwapLongFloatStartTime;
	delete itsSwapLongFloatEndTime;
	DeletePointorVector<ARM_GP_Vector>(itsSwapLongFixPayTimes);
	DeletePointorVector<ARM_GP_Vector>(itsSwapLongFixPayPeriods);
	delete itsSwapShortFloatStartTime;
	delete itsSwapShortFloatEndTime;
	DeletePointorVector<ARM_GP_Vector>(itsSwapShortFixPayTimes);
	DeletePointorVector<ARM_GP_Vector>(itsSwapShortFixPayPeriods);
	delete itsSwapPayFloatStartTime;
	delete itsSwapPayFloatEndTime;
	DeletePointorVector<ARM_GP_Vector>(itsSwapPayFixPayTimes);
	DeletePointorVector<ARM_GP_Vector>(itsSwapPayFixPayPeriods);
	delete itsFixValues;
	delete itsPayIndexLeverages;
	delete itsPeriodIndex;
	delete itsPayIndexResetTimes;
}


ARM_VanillaSpreadOptionArg::~ARM_VanillaSpreadOptionArg()
{
	CleanUp();
}

ARM_Object* ARM_VanillaSpreadOptionArg::Clone() const
{
	return new ARM_VanillaSpreadOptionArg(*this);
}

ARM_GenSecurityPtr ARM_VanillaSpreadOptionArg::VanillaSpreadOptionToGenSec() const
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ARM_RowInfo ARM_VanillaSpreadOptionArg::NumArgColumnNames() const
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*size_t colNamesSize = sizeof(SpreadOptionColNamesTable)/sizeof(SpreadOptionColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = SpreadOptionColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ARM_RowInfo ARM_VanillaSpreadOptionArg::NumArgMiddleRows( size_t eventIdx, const ARM_GP_VectorPtr& eventDates) const
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	size_t descSize = sizeof(SpreadOptionColNamesTable)/sizeof(SpreadOptionColNamesTable[0]);

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

    string zeroValue("0");

	rowDescVec[SpreadOption] = zeroValue;
    rowTypeVec[SpreadOption] = ARM_DOUBLE;

	string curveName = GetCurveName();
	
	// EventDate 
    double eventDate=eventDates->Elt(eventIdx);
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;

	// PayDate 
    double payDate=(*(itsPayDates))[eventIdx];
    CC_Ostringstream payDateDesc;
    payDateDesc << CC_NS(std,fixed) << payDate;
    rowDescVec[PayDate] = payDateDesc.str();
    rowTypeVec[PayDate] = ARM_DATE_TYPE;

	// Rate1StartDate 
    double rate1StartDate=(*(itsRate1StartDates))[eventIdx];
    CC_Ostringstream rate1StartDateDesc;
    rate1StartDateDesc << CC_NS(std,fixed) << rate1StartDate;
    rowDescVec[Rate1StartDate] = rate1StartDateDesc.str();
    rowTypeVec[Rate1StartDate] = ARM_DATE_TYPE;

	// Rate1EndDate 
    double rate1EndDate=(*(itsRate1EndDates))[eventIdx];
    CC_Ostringstream rate1EndDateDesc;
    rate1EndDateDesc << CC_NS(std,fixed) << rate1EndDate;
    rowDescVec[Rate1EndDate] = rate1EndDateDesc.str();
    rowTypeVec[Rate1EndDate] = ARM_DATE_TYPE;

	// Rate2StartDate 
    double rate2StartDate=(*(itsRate2StartDates))[eventIdx];
    CC_Ostringstream rate2StartDateDesc;
    rate2StartDateDesc << CC_NS(std,fixed) << rate2StartDate;
    rowDescVec[Rate2StartDate] = rate2StartDateDesc.str();
    rowTypeVec[Rate2StartDate] = ARM_DATE_TYPE;
	
	// Rate2EndDate 
    double rate2EndDate=(*(itsRate2EndDates))[eventIdx];
    CC_Ostringstream rate2EndDateDesc;
    rate2EndDateDesc << CC_NS(std,fixed) << rate2EndDate;
    rowDescVec[Rate2EndDate] = rate2EndDateDesc.str();
    rowTypeVec[Rate2EndDate] = ARM_DATE_TYPE;

    //Rate1
	CC_Ostringstream rate1Desc;
	rate1Desc << "SwapRate(" << curveName << "," << SpreadOptionColNamesTable[Rate1StartDate] << "[i],";
	rate1Desc << SpreadOptionColNamesTable[Rate1EndDate] << "[i])";
	rowDescVec[Rate1] = rate1Desc.str();
	rowTypeVec[Rate1] = ARM_STRING;

	//Rate2
	CC_Ostringstream rate2Desc;
	rate2Desc << "SwapRate(" << curveName << "," << SpreadOptionColNamesTable[Rate2StartDate] << "[i],";
	rate2Desc << SpreadOptionColNamesTable[Rate2EndDate] << "[i])";
	rowDescVec[Rate2] = rate2Desc.str();
	rowTypeVec[Rate2] = ARM_STRING;

	//IT
	CC_Ostringstream itDesc;
	itDesc << CC_NS(std,fixed) << (*itsPayPeriods)[eventIdx];
	rowDescVec[IT] = itDesc.str();
	rowTypeVec[IT] = ARM_DOUBLE;

   
	//Notional
	CC_Ostringstream notionalDesc;
	notionalDesc << CC_NS(std,fixed) << (*itsNotional)[eventIdx];
	rowDescVec[Notional] = notionalDesc.str();
	rowTypeVec[Notional] = ARM_DOUBLE;

    //Discount
	CC_Ostringstream discountDesc;
	discountDesc << "DF(" << curveName << "," << SpreadOptionColNamesTable[PayDate] << "[i])*";
	discountDesc << SpreadOptionColNamesTable[IT] << "[i]*";
	discountDesc << SpreadOptionColNamesTable[Notional] << "[i]";
	rowDescVec[Discount] = discountDesc.str();
	rowTypeVec[Discount] = ARM_STRING;

	
	//Coeff1
	CC_Ostringstream coeff1Desc;
	coeff1Desc << CC_NS(std,fixed) << (*itsCoeff1)[eventIdx];
	rowDescVec[Coeff1] = coeff1Desc.str();
	rowTypeVec[Coeff1] = ARM_DOUBLE;

	//Coeff2
	CC_Ostringstream coeff2Desc;
	coeff2Desc << CC_NS(std,fixed) << (*itsCoeff2)[eventIdx];
	rowDescVec[Coeff2] = coeff2Desc.str();
	rowTypeVec[Coeff2] = ARM_DOUBLE;

    //Strike
	CC_Ostringstream strikeDesc;
	strikeDesc << CC_NS(std,fixed) << (*itsStrikes)[eventIdx];
	rowDescVec[Strike] = strikeDesc.str();
	rowTypeVec[Strike] = ARM_DOUBLE;

    //SpreadOption
	CC_Ostringstream spreadOptionDesc;
	int callPut = GetCallPut();
	spreadOptionDesc << "Max(" ;
	spreadOptionDesc << CC_NS(std,fixed) << callPut <<"*(";
	spreadOptionDesc << SpreadOptionColNamesTable[Coeff2] << "[i]" << "*" << SpreadOptionColNamesTable[Rate2] << "[i]-";
	spreadOptionDesc << SpreadOptionColNamesTable[Coeff1] << "[i]"<< "*" << SpreadOptionColNamesTable[Rate1] << "[i])-";
	spreadOptionDesc << SpreadOptionColNamesTable[Strike]<< "[i]";
	spreadOptionDesc << ",0)" ;
	spreadOptionDesc << "*" << SpreadOptionColNamesTable[Discount] << "[i]";
	
	rowDescVec[SpreadOption] = spreadOptionDesc.str();
	rowTypeVec[SpreadOption] = ARM_STRING;

	return ARM_RowInfo(rowDescVec,rowTypeVec);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//ARM_GP_VectorPtr ARM_VanillaSpreadOptionArg::DatesStructure() const
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>( itsResetDates->Clone() ) );*/

	size_t descSize = sizeof(SpreadOptionColNamesTable);
	size_t ColumnSize = sizeof(SpreadOptionColNamesTable)/sizeof(SpreadOptionColNamesTable[0]);
	size_t RowSize = sizeof(SpreadOptionColNamesTable[0]);
	vector<string> text(descSize);
	vector<ARM_GP_VALUE_TYPE> format(descSize, ARM_MISSING_TYPE );
	ARM_DealDescriptionPtr pDealDescription( new ARM_DealDescription( text, format, RowSize,ColumnSize) );		
	return ARM_GenSecurityPtr( new ARM_GenSecurity( pDealDescription, GetCurveName())  );
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: AnalyticPrice
///	Returns: 
///	Action : analytic price of spreadoption with a model checking 
////////////////////////////////////////////////////
double ARM_VanillaSpreadOptionArg::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
	double price;
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    
    if(IRModel)
	{
		if (itsIsDigital)
		{
			ARM_VectorPtr PriceDown,PriceUp;
			ARM_GP_Vector strikesDown(*itsStrikes);
			strikesDown += itsSpread1;
			ARM_GP_Vector strikesUp(*itsStrikes);
			strikesUp += itsSpread2;
			PriceDown = IRModel->VanillaSpreadOptionScalar(
							GetCurveName(),
							GetEvalTime(),
							GetCallPut(),
							itsStartTime,
							itsEndTime,
							*itsResetTimes,
							*itsPayTimes,
							*itsPayPeriods,
							*itsNotional,
							*itsCoeffLong,
							*itsCoeffShort,
							strikesDown,
							*itsSwapLongFloatStartTime,
							*itsSwapLongFloatEndTime,
							itsSwapLongFixPayTimes,
							itsSwapLongFixPayPeriods,
							*itsSwapShortFloatStartTime,
							*itsSwapShortFloatEndTime,
							itsSwapShortFixPayTimes,
							itsSwapShortFixPayPeriods,
							dumStates);
			PriceUp = IRModel->VanillaSpreadOptionScalar(
							GetCurveName(),
							GetEvalTime(),
							GetCallPut(),
							itsStartTime,
							itsEndTime,
							*itsResetTimes,
							*itsPayTimes,
							*itsPayPeriods,
							*itsNotional,
							*itsCoeffLong,
							*itsCoeffShort,
							strikesUp,
							*itsSwapLongFloatStartTime,
							*itsSwapLongFloatEndTime,
							itsSwapLongFixPayTimes,
							itsSwapLongFixPayPeriods,
							*itsSwapShortFloatStartTime,
							*itsSwapShortFloatEndTime,
							itsSwapShortFixPayTimes,
							itsSwapShortFixPayPeriods,
							dumStates);
			price = ((*PriceDown)[0]-(*PriceUp)[0])/(itsSpread2-itsSpread1);

		}
		else
		{
			ARM_VectorPtr Price;
			Price = IRModel->VanillaSpreadOptionScalar(
							GetCurveName(),
							GetEvalTime(),
							GetCallPut(),
							itsStartTime,
							itsEndTime,
							*itsResetTimes,
							*itsPayTimes,
							*itsPayPeriods,
							*itsNotional,
							*itsCoeffLong,
							*itsCoeffShort,
							*itsStrikes,
							*itsSwapLongFloatStartTime,
							*itsSwapLongFloatEndTime,
							itsSwapLongFixPayTimes,
							itsSwapLongFixPayPeriods,
							*itsSwapShortFloatStartTime,
							*itsSwapShortFloatEndTime,
							itsSwapShortFixPayTimes,
							itsSwapShortFixPayPeriods,
							dumStates);
			price = (*Price)[0];
		}
	}
	else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price SpreadOption, please advise");
  
    return price;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSpreadOptionArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaSpreadOptionArg::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}
////////////////////////////////////////////////////
///	Struct : ARM_VanillaSpreadOptionArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaSpreadOptionArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaSpreadOptionArg"; 
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSpreadOptionArg
///	Routine: ComputeIndexSchedulesAndAdjustDates
///	Returns: string
///	Action : adjusts start & end dates and computes index fix pay times & periods
////////////////////////////////////////////////////
void ARM_VanillaSpreadOptionArg::ComputeIndexSchedulesAndAdjustDates (ARM_Currency* ccy, double asOfDate) 
{
	///----------------------------------------------------------------------
	/// NOTE : if indexes are LIBOR, we assume that (one-dim) schedule
	///		   have already been generated and that dates are already adjusted
	///	--> generate schedules only in CMS case

	char fixCalendar[100];
	ccy->CalcFixPayCal(fixCalendar);
	int fixFreq		= ccy->GetFixedPayFreq();
	int fixDayCount = ccy->GetFixedDayCount();

	size_t size = GetSwapLongFloatStartTime()->size();

	if (!itsLongIndexIsLibor)
	{
		// free existing vectors
		DeletePointorVector<ARM_GP_Vector>(itsSwapLongFixPayTimes);
		DeletePointorVector<ARM_GP_Vector>(itsSwapLongFixPayPeriods);
		itsSwapLongFixPayTimes.resize(size);
		itsSwapLongFixPayPeriods.resize(size);

		for (size_t i(0); i<size; i++)
		{
			// get (unadjusted) start/ end dates
			double longStartDate  = asOfDate + (*itsSwapLongFloatStartTime)[i];
			double longEndDate    = asOfDate + (*itsSwapLongFloatEndTime)[i];

			// generate date strip
			ARM_DateStrip longDateStrip(longStartDate, longEndDate, fixFreq, fixDayCount, fixCalendar,
										K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
										fixCalendar );


			// push
			itsSwapLongFixPayTimes[i] = new ARM_GP_Vector(*longDateStrip.GetPaymentDates());
			itsSwapLongFixPayPeriods[i] = new ARM_GP_Vector(*longDateStrip.GetInterestTerms());

			// overwrite start and end dates with adjusted dates
			(*itsSwapLongFloatStartTime)[i] = (*longDateStrip.GetFlowStartDates())[0] - asOfDate;
			(*itsSwapLongFloatEndTime)[i]   = (*longDateStrip.GetFlowEndDates())[longDateStrip.GetFlowEndDates()->size()-1] - asOfDate;

			// convert dates into times
			size_t schedSize = itsSwapLongFixPayTimes[i]->size();
			for (size_t j=0; j<schedSize; j++)
				(*itsSwapLongFixPayTimes[i])[j] -= asOfDate ;
		}
	}

	if (!itsShortIndexIsLibor)
	{
		// free existing vectors
		DeletePointorVector<ARM_GP_Vector>(itsSwapShortFixPayTimes);
		DeletePointorVector<ARM_GP_Vector>(itsSwapShortFixPayPeriods);
		itsSwapShortFixPayTimes.resize(size);
		itsSwapShortFixPayPeriods.resize(size);

		for (size_t i(0); i<size; i++)
		{
			// get (unadjusted) start/ end dates
			double shortStartDate  = asOfDate + (*itsSwapShortFloatStartTime)[i];
			double shortEndDate    = asOfDate + (*itsSwapShortFloatEndTime)[i];

			// generate date strip
			ARM_DateStrip shortDateStrip(shortStartDate, shortEndDate, fixFreq, fixDayCount, fixCalendar,
										K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
										fixCalendar );


			// push
			itsSwapShortFixPayTimes[i] = new ARM_GP_Vector(*shortDateStrip.GetPaymentDates());
			itsSwapShortFixPayPeriods[i] = new ARM_GP_Vector(*shortDateStrip.GetInterestTerms());

			// overwrite start and end dates with adjusted dates
			(*itsSwapShortFloatStartTime)[i] = (*shortDateStrip.GetFlowStartDates())[0] - asOfDate;
			(*itsSwapShortFloatEndTime)[i]   = (*shortDateStrip.GetFlowEndDates())[shortDateStrip.GetFlowEndDates()->size()-1] - asOfDate;

			// convert dates into times
			size_t schedSize = itsSwapShortFixPayTimes[i]->size();
			for (size_t j=0; j<schedSize; j++)
				(*itsSwapShortFixPayTimes[i])[j] -= asOfDate ;
		}
	}

	if (itsPayIndexType == K_CMS)
	{
		/// size = nb of payment periods (not nb of observation dates)
		size = itsSwapPayFloatStartTime->size();

		// free existing vectors
		DeletePointorVector<ARM_GP_Vector>(itsSwapPayFixPayTimes);
		DeletePointorVector<ARM_GP_Vector>(itsSwapPayFixPayPeriods);
		itsSwapPayFixPayTimes.resize(size);
		itsSwapPayFixPayPeriods.resize(size);

		for (size_t i(0); i<size; i++)
		{
			// get (unadjusted) start/ end dates
			double payStartDate  = asOfDate + (*itsSwapPayFloatStartTime)[i];
			double payEndDate    = asOfDate + (*itsSwapPayFloatEndTime)[i];

			// generate date strip
			ARM_DateStrip payDateStrip(payStartDate, payEndDate, fixFreq, fixDayCount, fixCalendar,
										K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
										fixCalendar );


			// push
			itsSwapPayFixPayTimes[i]   = new ARM_GP_Vector(*payDateStrip.GetPaymentDates());
			itsSwapPayFixPayPeriods[i] = new ARM_GP_Vector(*payDateStrip.GetInterestTerms());

			// overwrite start and end dates with adjusted dates
			(*itsSwapPayFloatStartTime)[i] = (*payDateStrip.GetFlowStartDates())[0] - asOfDate;
			(*itsSwapPayFloatEndTime)[i]   = (*payDateStrip.GetFlowEndDates())[payDateStrip.GetFlowEndDates()->size()-1] - asOfDate;

			// convert dates into times
			size_t schedSize = itsSwapPayFixPayTimes[i]->size();
			for (size_t j=0; j<schedSize; j++)
				(*itsSwapPayFixPayTimes[i])[j] -= asOfDate ;
		}
	}

}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


