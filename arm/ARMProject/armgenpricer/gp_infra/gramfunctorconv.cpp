/*!
 *
 *	\file gramfunctorconv.cpp
 *
 *  \brief gramfunctorconv is the file to put all convention based and default argument
 *	computation used for the evaluation of gramnode
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

#include "gpinfra/gramfunctorconv.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gramfunctiontable.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/pricingmodel.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/stringconvert.h"

/// kernel
//#include <inst/irindex.h>
#include <ccy/currency.h>

/// for easy debugging of shared nodes!
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
	#include "gpbase/pair.h"
#endif


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Routine: GetAndCheckEvalTime
///	Returns: void
///	Action : checks that the eval date is not in the past
///				throw an exception if it is the case!
////////////////////////////////////////////////////

double GetAndCheckEvalTime( double evalDate, ARM_PricingModel* mod, const string& funcName )
{
	double evalTime	= mod->GetTimeFromDate( evalDate );
	if( evalTime < 0 )
	{
		CC_Ostringstream os;
		os << "Trying to use a " << funcName << " in the past: eval date " 
			<< ARM_Date( evalDate ).toString() 
			<< " while model as of date " << mod->GetAsOfDate().toString() << " "
			<< ARM_USERNAME  << ": please advise";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	return evalTime;
}


////////////////////////////////////////////////////
///	Routine: ComputePayDate
///	Returns: ARM_Date
///	Action : Compute a pay date
////////////////////////////////////////////////////
ARM_Date ComputeResetOrPayDate(const ARM_Date& advDate,const ARM_Date& arrDate,
                        int timing,char* calendar,
                        const ARM_GramFctorArg& argGapOrDate,
                        ARM_PricingModel* mod, const string& curveName,
                        const string& funcName)
{
    ARM_Date resetOrPayDate(timing == K_ADVANCE ? advDate : arrDate);
    if( GFAT_DOUBLE_TYPE == argGapOrDate.GetType() )
    {
        double gap = argGapOrDate.GetDouble();
	    if( DEFAULT_CCY_DOUBLE == gap  )
	    {
            /// Default value for resetGap !
		    ARM_Currency* ccy	= mod->GetCurrency( curveName );
		    gap        = -ccy->GetSpotDays();
	    }
		resetOrPayDate.GapBusinessDay(gap,calendar);
    }
    else if(GFAT_DATE_TYPE == argGapOrDate.GetType())
        resetOrPayDate=argGapOrDate.GetDate();
    else
    {
		CC_Ostringstream os;
		os << "In " << funcName << " can't compute a reset or payment date using start & end dates " 
			<< ARM_USERNAME  << ": please advise";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

    return resetOrPayDate;
}


////////////////////////////////////////////////////
///	Routine: ComputeSettlementDate
///	Returns: ARM_Date
///	Action : Compute a settlement date
////////////////////////////////////////////////////
ARM_Date ComputeSettlementDate( 
	const ARM_Date& startDate, 
	char* calendar, 
	const ARM_GramFctorArg& argGapOrDate,
	int defaultGap,
    const string& funcName )
{
    if( GFAT_DOUBLE_TYPE == argGapOrDate.GetType() )
    {
		ARM_Date settlementDate( startDate );
        double gap = argGapOrDate.GetDouble()== DEFAULT_CCY_DOUBLE? defaultGap : argGapOrDate.GetDouble();
		settlementDate.GapBusinessDay(gap,calendar);
		return settlementDate;
    }
    else if( GFAT_DATE_TYPE == argGapOrDate.GetType())
    {
		return argGapOrDate.GetDate();
	}
    else
    {
		CC_Ostringstream os;
		os << "In " << funcName << " can't compute a settlement date " 
			<< ARM_USERNAME  << ": please advise";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }
}

/// function to set the settlement date
void SetSettlementDate( 
	vector< ARM_ExpNodePtr >& nodes,
	const ARM_Date& startDate, 
	const string& calendar,
	int defaultGap,
	ARM_GramFctorArgVector& args, /// changed potentially
	size_t i,
    const string& funcName )
{
    if( GFAT_DOUBLE_TYPE == args[i].GetType() )
    {
		ARM_Date settlementDate( startDate );
        double gap = args[i].GetDouble()== DEFAULT_CCY_DOUBLE? defaultGap : args[i].GetDouble();
		settlementDate.GapBusinessDay( gap, const_cast<char*>(calendar.c_str()) );
        nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeDate( settlementDate ) );
		args[i].SetDate( settlementDate );
    }
    else if( GFAT_DATE_TYPE == args[i].GetType())
    {
		return;
	}
    else
    {
		CC_Ostringstream os;
		os << "In " << funcName << " can't compute a settlement date " 
			<< ARM_USERNAME  << ": please advise";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }
}


double SetSettlementDate( 
						const double& julian, 
						const string& calendar,
						int defaultGap,
						double value
						)
{
	ARM_Date settlementDate( julian );
    double gap = ((value== DEFAULT_CCY_DOUBLE)? defaultGap : value);
	settlementDate.GapBusinessDay( gap, const_cast<char*>(calendar.c_str()) );
    return settlementDate.GetJulian();
}

////////////////////////////////////////////////////
///	Routine: ComputeIndexTerm
///	Returns: string
///	Action : Compute the index term
////////////////////////////////////////////////////

string ComputeIndexTerm( const ARM_GramFctorArg& argTerm, const ARM_GramFctorArg& argStart,
                         const ARM_GramFctorArg& argEnd,const string& funcName)
{
    string termString = argTerm.GetString();

	if(termString == "" )
	{
        if( GFAT_STRING_TYPE == argEnd.GetType() )
            //// use this maturity as the index term
            termString = argEnd.GetString();
        else if(GFAT_DATE_TYPE == argEnd.GetType())
        {
            /// compute a maturity from start & end
            int termInMonths=(int)floor((argEnd.GetDate().GetJulian()-argStart.GetDate().GetJulian())/30.44 + 0.5);
            if(termInMonths<1)
                termInMonths=1;
            CC_Ostringstream termInMonthsStr;
            termInMonthsStr << termInMonths << "M";
            termString=termInMonthsStr.str();
        }
        else
        {
		    CC_Ostringstream os;
		    os << "In " << funcName << " can't compute a index term using start & end dates " 
			    << ARM_USERNAME  << ": please advise";
		    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
        }
	}
	return termString;
}


////////////////////////////////////////////////////
///	Routine: SetIndexTerm
///	Returns: void
///	Action : Sets in the appropriate node the index term
////////////////////////////////////////////////////

void SetIndexTerm( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i,
                   size_t startIdx, size_t endIdx, ARM_PricingModel* mod, const string& curveName,
                   const string& funcName)
{
	nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeString( ComputeIndexTerm(args[i],args[startIdx],args[endIdx],funcName) ) );
}


////////////////////////////////////////////////////
///	Routine: GetResetDaysGap
///	Returns: double (the resetdaysGap)
///	Action : according to the type of the arg, either
///				gets the default reset day gap from the currency
///				or takes the double value
////////////////////////////////////////////////////

double GetResetDaysGap( const ARM_GramFctorArg& arg, ARM_PricingModel* mod, const string& curveName)
{
	double resetDaysGap = arg.GetDouble();
	if( DEFAULT_CCY_DOUBLE == resetDaysGap  )
	{
		ARM_Currency* ccy	= mod->GetCurrency( curveName );
		resetDaysGap        = -ccy->GetSpotDays();
	}
	return resetDaysGap;
}



////////////////////////////////////////////////////
///	Routine: SetResetDaysGap
///	Returns: void
///	Action : Sets in the appropriate node the reset Default
///			Modifies the parsed node if appropriate
////////////////////////////////////////////////////

void SetResetDaysGap( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& curveName )
{
    if( GFAT_DOUBLE_TYPE == args[i].GetType() )
        nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeDateOrDouble( GetResetDaysGap(args[i],mod,curveName)  ) );
}


////////////////////////////////////////////////////
///	Routine: SetEndDateResetCal, SetEndDatePayCal
///	Returns: void
///	Action : Sets in the appropriate node the end date
///			Modifies the parsed node if appropriate
////////////////////////////////////////////////////
void SetEndDateResetCal( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& curveName )
{
	/// Modification of nodes if appropriate
	/*if( args[i].GetType() == GFAT_STRING_TYPE )
	{
		ARM_Currency* ccy			= mod->GetCurrency( curveName );
		char* ccyName				= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( ccyName );
		char* resetCalendar		    = ccy->GetResetCalName(defaultIndex);

        nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeDate( ComputeEndDate(args[i-1].GetDate(), args[i], resetCalendar ) ) );

        delete resetCalendar;
	};*/
}

void SetEndDatePayCal( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& curveName )
{
	/// Modification of nodes if appropriate
	/*if( args[i].GetType() == GFAT_STRING_TYPE )
	{
		ARM_Currency* ccy			= mod->GetCurrency( curveName );
		char* ccyName				= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( ccyName );
		char* payCalendar		    = ccy->GetPayCalName(defaultIndex);

        nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeDate( ComputeEndDate(args[i-1].GetDate(), args[i], payCalendar ) ) );

		delete payCalendar;
	}*/;
}


////////////////////////////////////////////////////
///	Routine: ComputeEndDate
///	Returns: ARM_Date
///	Action : using the arg, the startdate and the calendar
///				computes the end date
////////////////////////////////////////////////////

ARM_Date ComputeEndDate( const ARM_Date& startDate, const ARM_GramFctorArg& arg,const char* calendar ) 
{
	/// the second argument can either be a date or a maturity
	ARM_Date endDate;

	if( GFAT_STRING_TYPE == arg.GetType() )
	{
		endDate = startDate;
		string tenorString	= arg.GetString();
		endDate.AddPeriod( tenorString, calendar );
		if (calendar)
			endDate.GoodBusinessDay( K_MOD_FOLLOWING, (char*) calendar );
	}
	else
		endDate = arg.GetDate();

	return endDate;
}


////////////////////////////////////////////////////
///	Routine: GetDayCountFtor related functions
////////////////////////////////////////////////////
long GetFixedDayCountFtor::operator() () const
{
    return itsCcy->GetFixedDayCount();
}

long GetFloatDayCountFtor::operator() () const
{
    return itsCcy->GetLiborIndexDayCount();
}

////////////////////////////////////////////////////
///	Routine: GetDayCount
///	Returns: long
///	Action : using the arg and the ccy returns either the default dayCount
///             or the dayCount corresponding to the string
////////////////////////////////////////////////////

int GetDayCount( const ARM_GramFctorArg& arg, ARM_PricingModel* mod, const string& curveName, const GetDayCountFtor& getDayCount)
{
    long dayCount;
    string dayCountString = arg.GetString();

	if(dayCountString == DEFAULT_CCY_CHAR )
        dayCount = getDayCount();
	else
		dayCount = ARM_ArgConv_DayCount.GetNumber( dayCountString );

    return dayCount;
}




////////////////////////////////////////////////////
///	Routine: SetDayCount
///	Returns: long
///	Action : set method for day count
////////////////////////////////////////////////////

bool SetDayCount( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& curveName, const GetDayCountFtor& getDayCount)
{
    bool isDefault;
    long dayCount;
    string dayCountString = args[i].GetString();

	if( (isDefault = (dayCountString == DEFAULT_CCY_CHAR)) )
	{
		ARM_Currency* ccy = mod->GetCurrency( curveName );

        dayCount = getDayCount();

		string dayCountString = ARM_ArgConvReverse_DayCount.GetString( dayCount );
		nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeString( dayCountString ) );
	}

    return isDefault;
}


////////////////////////////////////////////////////
///	Routine: GetFrequencyFtor related functions
////////////////////////////////////////////////////
double GetFixedFrequencyFtor::operator() () const
{
    return itsCcy->GetFixedPayFreq();
}

double GetFloatFrequencyFtor::operator() () const
{
    return itsCcy->GetLiborTerm();
}


////////////////////////////////////////////////////
///	Routine: GetTerm
///	Returns: long
///	Action : using the arg and the ccy returns either the default term
///             or the term corresponding to the string
////////////////////////////////////////////////////

double GetTerm( const ARM_GramFctorArg& arg, ARM_PricingModel* mod, const string& curveName, string functionName, const GetFrequencyFtor& getFrequency)
{
	double term;

	std::string periodString = arg.GetString();

	if (periodString == DEFAULT_CCY_CHAR)
		term = 1.0/getFrequency();
	else
    {
		term = StringMaturityToYearTerm(periodString);

	    if( term <= 0 )
	    {
		    CC_Ostringstream os;
		    os << functionName + " : invalid frequency: found " << periodString
			    << "corresponding to " << term << " "
			    << ARM_USERNAME << ": please advise!";
		    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	    }
    }

	return term;
}


////////////////////////////////////////////////////
///	Routine: SetTerm
///	Returns: long
///	Action : set method for term with the functor
////////////////////////////////////////////////////
bool SetTerm( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& curveName, const GetFrequencyFtor& getFrequency)
{
	return SetTerm(nodes,args,i,mod,curveName,1.0/getFrequency());
}


////////////////////////////////////////////////////
///	Routine: SetTerm
///	Returns: long
///	Action : set method for frequency 
////////////////////////////////////////////////////
bool SetTerm(vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, double term)
{
	bool isDefault;
	std::string periodString = args[i].GetString();

	if( (isDefault = (periodString == std::string( DEFAULT_CCY_CHAR ))) )
	{
		nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeString( YearTermToStringMaturity( term  ) ) );
	}

    return isDefault;
}


////////////////////////////////////////////////////
///	Routine: GetFrequency
///	Returns: long
///	Action : using the arg and the ccy returns either the default frequency
///             or the frequency corresponding to the string
////////////////////////////////////////////////////

int GetFrequency( const ARM_GramFctorArg& arg, ARM_PricingModel* mod, const string& curveName, string functionName, const GetFrequencyFtor& getFrequency)
{
	int freq;

	std::string periodString = arg.GetString();

	if (periodString == DEFAULT_CCY_CHAR)
		freq = getFrequency();
	else
    {
		///
		/// case #1 : input freq of type 1Y / 6M / 3M / 1M etc...
		///
		if (isdigit(periodString[0]))
		{
			double term = 1.0/StringMaturityToYearTerm(periodString);

			if( term <= 0 )
			{
				CC_Ostringstream os;
				os << functionName + " : invalid frequency: found " << periodString
					<< "corresponding to " << term << " "
					<< ARM_USERNAME << ": please advise!";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
			freq = (int)term;
		}
		///
		/// case #2 : input freq of type A, S, Q, B, M, W, D, Z
		///
		else
		{
			if(periodString == "A")
			{
				freq = K_ANNUAL;
			}
			else if(periodString == "S")
			{
				freq = K_SEMIANNUAL;
			}
			else if(periodString == "Q")
			{
				freq = K_QUARTERLY;
			}
			else if(periodString == "B")
			{
				freq = K_BIMONTHLY;
			}
			else if(periodString == "M")
			{
				freq = K_MONTHLY;
			}
			else if(periodString == "W")
			{
				freq = K_WEEKLY;
			}
			else if(periodString == "D")
			{
				freq = K_DAILY;
			}
			else if(periodString == "Z")
			{
				freq = K_ZEROCOUPON;
			}
			else
			{
				Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetFrequency : valid frequencies are 1Y, 6M, 3M, 1M, 1W, 1D or A, S, Q, B, M, W, D, Z" );
			}
		}
    }

	return freq;
}


////////////////////////////////////////////////////
///	Routine: SetFrequency
///	Returns: long
///	Action : set method for frequency 
////////////////////////////////////////////////////
bool SetFrequency( vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& curveName, const GetFrequencyFtor& getFrequency)
{
	return SetFrequency(nodes,args,i,mod,curveName,getFrequency());
}

////////////////////////////////////////////////////
///	Routine: SetFrequency
///	Returns: long
///	Action : set method for frequency 
////////////////////////////////////////////////////
bool SetFrequency(vector< ARM_ExpNodePtr >& nodes, const ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, int frequency)
{
	bool isDefault;
	std::string periodString = args[i].GetString();

	if( (isDefault = (periodString == std::string( DEFAULT_CCY_CHAR ))) )
	{
		nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeString( YearTermToStringMaturity( 1.0/frequency  ) ) );
	}

    return isDefault;
}

////////////////////////////////////////////////////
///	Routine: SetFrequency
///	Returns: long
///	Action : set method for frequency 
////////////////////////////////////////////////////
bool SetFrequencyToGramFctorArg(vector< ARM_ExpNodePtr >& nodes, ARM_GramFctorArgVector& args, size_t i, ARM_PricingModel* mod, const string& CurveName, const string& freqStr)
{
	bool isDefault;
	std::string periodString = args[i].GetString();

	if( (isDefault = (periodString == std::string( DEFAULT_CCY_CHAR ))) )
	{
		///nodes[i] = ARM_ExpNodePtr( new ARM_ExpNodeString( freqStr ) );
		((ARM_GramFctorArgVector&)args)[i].SetString(freqStr);
	}

    return isDefault;
}


////////////////////////////////////////////////////
///	Routine: GetCurve
///	Returns: long
///	Action : Get the curve in the node
////////////////////////////////////////////////////
ARM_GP_CurvePtr GetCurve(const ARM_GramFctorArg& arg)
{
	if (arg.GetType() == GFAT_CURVE_TYPE)
	{
		return arg.GetCurve();
	}
	else if (arg.GetType() == GFAT_DOUBLE_TYPE)
	{
		std::vector<double> abs(1,0.0);
		std::vector<double> ord(1,arg.GetDouble());
		return ARM_GP_CurvePtr(NULL);//new ARM_Curve(abs,ord,new ARM::ARM_StepUpRightOpenCstExtrapolDble));
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "This cannot be convert into a curve." );
	}
}

////////////////////////////////////////////////////
///	Routine: GetCurve
///	Returns: long
///	Action : Get the curve in the node
////////////////////////////////////////////////////
ARM_GP_VectorPtr GetProfileValues(const ARM_GramFctorArg& arg, size_t offset,
                               const std::vector<double>& times, const string& msg)
{
    ARM_GP_VectorPtr values;

    size_t i,nbOutValues = times.size();;

    if( arg.GetType() ==  GFAT_DOUBLE_TYPE )
		values	= ARM_GP_VectorPtr(new ARM_GP_Vector(nbOutValues,arg.GetDouble()));

	else if( arg.GetType() ==  GFAT_VECTOR_TYPE )
    {
        values = ARM_GP_VectorPtr(new ARM_GP_Vector(nbOutValues));

        ARM_VectorPtr inValues = arg.GetVector();
        size_t nbInValues = inValues->size();
        size_t nbMin = (nbInValues < nbOutValues+offset ? nbInValues : nbOutValues+offset);
        for(i=0;i+offset<nbMin;++i)
            (*values)[i] = (*inValues)[i+offset];
        for(;i<nbOutValues;++i)
            (*values)[i] = (*inValues)[nbInValues-1];
    }
    else if( arg.GetType() ==  GFAT_CURVE_TYPE )
	{
        ARM_GP_Vector outValues(nbOutValues);

		ARM_GP_CurvePtr inValuesCurve = arg.GetCurve();
        for(i=0; i<nbOutValues; ++i)
            outValues[i] = inValuesCurve->Interpolate(times[i]);

		values = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(outValues.Clone()));
	}
	else
    {
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + msg +" should be double or vector" );
    }

    return values;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


