//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MaturityFromRollDate.cpp
//
//   Description : Defines floating expiries used to define yield curve & 
//                 vol surface points e.g. 1M, 5Y
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_MATURITY_FROM_ROLL_DATE_CPP
#include "edginc/MaturityFromRollDate.hpp"
#include "edginc/Format.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Hashtable.hpp"

DRLIB_BEGIN_NAMESPACE


/** Override clone method to copy our extra data over */
IObject* MaturityFromRollDate::clone() const{
    int& count = getRefCount();
    if (count == 0){
        return new MaturityFromRollDate(*this);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}

// implemented for fast clone
MaturityFromRollDate::MaturityFromRollDate(const MaturityFromRollDate& matPeriod):
    Expiry(TYPE), maturityPeriod(matPeriod.maturityPeriod),period(matPeriod.period), 
    rollDay(matPeriod.rollDay), rollMonths(matPeriod.rollMonths) {}


/** Returns true if given expiry matches this */
bool MaturityFromRollDate::equals(const Expiry* expiry) const{
    return maturityPeriod->equalTo(expiry);
}



/** this gets called after an object is constructed from a data dictionary.
    Not after an object has been copied (see override of clone method below) */
void MaturityFromRollDate::validatePop2Object()
{
    static const string method = "MaturityFromRollDate::validatePop2Object";
	try
	{
		int i;


		// need to set maturityPeriod if hasn't been done already
		// need to do this here since can't in private constructor because don't have period yet
		
		if(!maturityPeriod) maturityPeriod = MaturityPeriodSP(new MaturityPeriod(period));

		// check rollMonths
		for(i = 0;i< rollMonths.size();i++)
		{
			if(rollMonths[i] < 1 || rollMonths[i] > 12) throw ModelException("rollMonths out of bound [1,12]");
		}

		// check if dates are valid by constructing MonthDayYear for non-leap year = 2005
		// then convery to DateTime. This will fail if invalid date
		for(i = 0;i< rollMonths.size();i++)
		{
			DateTime::MonthDayYear dmy(rollDay,rollMonths[i],2005);
			DateTime test = dmy.toDateTime();
		}


	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
	
}

MaturityFromRollDate::MaturityFromRollDate(const string& period) : 
    Expiry(TYPE), maturityPeriod(0), period(period), rollDay(20)
	{
    
		// initialise rollMonths to credit IMM months
		IntArray months(4);
		months[0] = 3;
		months[1] = 6;
		months[2] = 9;
		months[3] = 12;

		const_cast<IntArray&>(rollMonths) = months;

    validatePop2Object();
}

MaturityFromRollDate::MaturityFromRollDate(int count, const string& interval) :
    Expiry(TYPE), period(Format::toString(count)+interval), rollDay(20)
{

		maturityPeriod = MaturityPeriodSP(new MaturityPeriod(count,interval));

		// initialise rollMonths to credit IMM months
		IntArray months(4);
		months[0] = 3;
		months[1] = 6;
		months[2] = 9;
		months[3] = 12;

		const_cast<IntArray&>(rollMonths) = months;

		validatePop2Object();	 
}


MaturityFromRollDate::~MaturityFromRollDate() {
    // empty
}

string MaturityFromRollDate::toString() const {
    return period;
}

DateTime MaturityFromRollDate::toDate(const DateTime& aDate) const {
    static const string method = "MaturityFromRollDate::toDate";
    try {
		DateTime temp = maturityPeriod->toDate(aDate);
		DateTime afterTemp = getNextRollDate(maturityPeriod->toDate(aDate));
        return afterTemp;
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed");
    } 
}


/** get first Roll date on or after date*/
DateTime MaturityFromRollDate::getNextRollDate(const DateTime & date) const
{
	return getNextRollDate(date, rollDay, rollMonths);
}

/** get first Roll date on or after date*/
DateTime MaturityFromRollDate::getNextRollDate(const DateTime & date, int _rollDay, const IntArray & _rollMonths) 
{
	 int day = date.toMDY().day;
	 int month = date.toMDY().month;
	 int year = date.toMDY().year;

	 //check if day > _rollDay. If it is move month by 1
	 if(day> _rollDay) month += 1;
	 if(month > 12)
	 {
		 month = 1;
		 year += 1;
	 }

	 // find first IMM month greater or equal to nextMonth
	int i =0;
	while(i < _rollMonths.size() &&  _rollMonths[i] < month) i++;

	/** first IMM month greater or equal to month */
	int nextIMMmonth;
	/** first IMM year greater or equal to year */
	int nextIMMyear = year;
	if(i < _rollMonths.size())
	{
		nextIMMmonth = _rollMonths[i];
	}
	else
	{
		nextIMMmonth = _rollMonths[0];
		nextIMMyear += 1;
	}


	DateTime::MonthDayYear out(_rollDay,nextIMMmonth,nextIMMyear);

	return out.toDateTime();

}

/* for reflection */
MaturityFromRollDate::MaturityFromRollDate():Expiry(TYPE), maturityPeriod(0), rollDay(20)
{
	// initialise rollMonths to credit IMM months
	IntArray months(4);
	months[0] = 3;
	months[1] = 6;
	months[2] = 9;
	months[3] = 12;

	const_cast<IntArray&>(rollMonths) = months;
}

class MaturityFromRollDateHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MaturityFromRollDate, clazz);
        SUPERCLASS(Expiry);
        EMPTY_SHELL_METHOD(defaultMaturityFromRollDate);
        clazz->enableCloneOptimisations(); // assignment operator clones
        FIELD(period, "period eg 1W");
		FIELD(rollDay, "roll day in month (default = 20)");
		FIELD_MAKE_OPTIONAL(rollDay);
		FIELD(rollMonths, "months in which have roll date (default = {3,6,9,12})");
		FIELD_MAKE_OPTIONAL(rollMonths);
    }
    
    static IObject* defaultMaturityFromRollDate(){
        return new MaturityFromRollDate();
    }
};

CClassConstSP const MaturityFromRollDate::TYPE = CClass::registerClassLoadMethod(
    "MaturityFromRollDate", typeid(MaturityFromRollDate), MaturityFromRollDateHelper::load);

DEFINE_TEMPLATE_TYPE(MaturityFromRollDateArray);




DRLIB_END_NAMESPACE
