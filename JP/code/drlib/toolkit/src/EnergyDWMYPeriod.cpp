//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyDWMYPeriod.cpp
//
//   Description : Defines and decompose nD(WMY) periods used in energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/EnergyDWMYPeriod.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

#define MAX_BUFFER 32
#define VALID_MATURITIES "DMYW"

EnergyDWMYPeriod::EnergyDWMYPeriod(const EnergyDWMYPeriod& thePeriod):
    CObject(TYPE), period(thePeriod.period), 
    interval(thePeriod.interval), count(thePeriod.count){}

void EnergyDWMYPeriod::validatePop2Object()
{
    static const string routine = "EnergyDWMYPeriod::validatePop2Object";
   
    // chop string (e.g. 6M) into count (i.e. 6) and interval (i.e. M)

    const char *inp = period.c_str();
    char  numberBuff[MAX_BUFFER];
    char *nump = numberBuff;            /* Pointer to number */

    /* Copy sign,if any, to numberBuff */
    if (*inp == '-' || *inp == '+')
    {
        *nump++ = *inp++;
    }

    /* Copy digits, if any, to numberBuff */
    while (isdigit(*inp)) { 
        *nump++ = *inp++;
    }

    *nump = '\0';                       /* Null terminate */

    if (inp != period.c_str())
	{
        // Found some digits 
        const_cast<int&>(count) = atoi(numberBuff);
    }
    else
	{
        const_cast<int&>(count) = 1;
    }
    // dangerous to build a string from a char which is 0
    char intervalAsChar = toupper(*inp);
    if (intervalAsChar == '\0' || *(inp+1) != '\0' ||
        !strchr(VALID_MATURITIES, intervalAsChar))
	{
        if (intervalAsChar != '\0')
		{
            const_cast<string&>(interval) = intervalAsChar;
        }
        throw ModelException(routine, "Invalid period: '"+interval+
                             "' for EnergyDWMYPeriod "+period+". "+
                                 "Period must be D, M,Y,or W");
    }

	if (intervalAsChar == 'W')
	{
		intervalAsChar = 'D';
		count *=7;
	}
	else if (intervalAsChar == 'Y')
	{
		intervalAsChar = 'M';
		count *=12;
	}
		
    const_cast<string&>(interval) = intervalAsChar;
   
}

EnergyDWMYPeriod::EnergyDWMYPeriod(const string& period) : 
    CObject(TYPE), period(period), count(0)
{
    // populate count and interval
    validatePop2Object();
}

EnergyDWMYPeriod& EnergyDWMYPeriod::operator=(const EnergyDWMYPeriod& thePeriod)
{
    if( this != &thePeriod )
	{
		period = thePeriod.period; 
        interval = thePeriod.interval;
		count= thePeriod.count;
	}
	return *this;
}

/****
IObject* EnergyDWMYPeriod::clone() const {
    static const string routine = "EnergyDWMYPeriod::clone";
    try {
        IObject* Copy = CObject::clone(); // call parent
        EnergyDWMYPeriod* copy = dynamic_cast<EnergyDWMYPeriod*>(Copy);
        if(!copy) {
                throw ModelException("Clone method failed");
        }
		copy->period = period;
        copy->interval = interval;
        copy->count = count;
       
        return copy;
    } catch(exception& e) {
            throw ModelException(e, routine);
    }
}
**/

EnergyDWMYPeriod::~EnergyDWMYPeriod()
{
}

string EnergyDWMYPeriod::getPeriod() const
{
    return period;
}

int EnergyDWMYPeriod::getCount() const
{
    return count;
}
string EnergyDWMYPeriod::getInterval() const
{
    return interval;
}

// Overwrite the period to that between contracts. Hidden from constructing interface.
void EnergyDWMYPeriod::setToContractPeriod(int theContractType)
{
	if (theContractType == 0 )
        interval = "I"; // for EnergyUnderlyer::INDEX
	else
		interval = "O"; // for EnergyUnderlyer::OPTION
}

EnergyDWMYPeriod::EnergyDWMYPeriod():CObject(TYPE), count(0){}

/* for reflection */
class EnergyDWMYPeriodHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyDWMYPeriod, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEnergyDWMYPeriod);
        FIELD(period, "period eg 1W");
		FIELD(interval, "interval");
		FIELD_MAKE_TRANSIENT(interval);
		FIELD(count, "count");
		FIELD_MAKE_TRANSIENT(count);
    }
    
    static IObject* defaultEnergyDWMYPeriod(){
        return new EnergyDWMYPeriod();
    }
};

CClassConstSP const EnergyDWMYPeriod::TYPE = CClass::registerClassLoadMethod(
    "EnergyDWMYPeriod", typeid(EnergyDWMYPeriod), EnergyDWMYPeriodHelper::load);


DRLIB_END_NAMESPACE
