//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyTermPeriodBenchmark.cpp
//
//   Description : Defines a deal period with a benchmark for energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/EnergyTermPeriodBenchmark.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE


void EnergyTermPeriodBenchmark::validatePop2Object()
{
    static const string routine = "EnergyTermPeriodBenchmark::validatePop2Object"; 
	
	// First check if it's PROMPT from FutureSwapSchedule interface
	if ( benchmark == "PROMPT")
	{
		// Pass responsibility to outside for determining the labels
		prompt = true;
		return;
	}

	// Per logic in DRConvertBenchmark...

	char *dash_ptr, dash;
	int dashPosition;
	char marker;
	int month, year;
	PeriodFormat periodFormat;
	const char *in = benchmark.c_str();
	int len = benchmark.size();
	char * out = new char[len+1];
	int i;

	for (i=0; i<len; ++i)
	{
		char c = in[i];
		c = toupper(c);

		if (isgraph(c))
			out[i] = c;
	}
	out[i] = '\n'; 

	if (!(dash_ptr = (/*const_cast*/char*)strchr(out, '-')))  // if no - character return
		throw ModelException(routine,"no - character.");

	dashPosition = dash_ptr - out;

	if (dashPosition == 1)  // A-year format
	{
		char period;
		if (sscanf(out, "%c%c%d%c", &period, &dash, &year, &marker) != 4)
		{
			throw ModelException(routine,"Must be 4 chars for A-year format");
		}
		else if (period != 'A' || dash != '-' ||
			     marker != '\n' || !adjustYear(&year) )
		{
			throw ModelException(routine,"Some wrong char in A-year format");
		}
		
		month = 12;  // month set to end of year (december)
		periodFormat = ANNUAL;
	}
	else if (dashPosition == 2)  // (Q|S)1-year format or IN-year (injection) format
	{
		char label[3];

		if (sscanf(out, "%c%c%c%d%c", &label[0], &label[1], &dash, &year, &marker) != 5)
		{
			throw ModelException(routine,"Must be 5 chars for (Q|S)1-year format or IN-year (injection) format");
		}
		else if (dash != '-' || marker != '\n' || !adjustYear(&year) )
		{
			throw ModelException(routine,"Some wrong char in (Q|S)1-year format or IN-year (injection) format");
		}

		if (label[0] == 'Q')
		{
			label[2] = '\0';  // must terminate with null for atoi
			month = atoi(&label[1]) * 3;
			if (!month)
				throw ModelException(routine,"wrong number for Q in (Q|S)1-year format or IN-year (injection) format");
			
			periodFormat = QUARTER;
		}
		else if (label[0] == 'S')
		{
			label[2] = '\0';
			month = atoi(&label[1]) * 6;
			if (!month)
				throw ModelException(routine,"wrong number for S in (Q|S)1-year format or IN-year (injection) format");
			periodFormat = SEMIANNUAL;
		}
		else if (label[0] == 'I' && label[1] == 'N')  // injection season
		{
			month = 10;
			periodFormat = INJECTION;
		}
		else
			throw ModelException(routine,"Something wrong in (Q|S)1-year format or IN-year (injection) format");
		

		if (month < 0 || month > 12)
			throw ModelException(routine,"Some wrong number in (Q|S)1-year format or IN-year (injection) format");
		
	}
	else if (dashPosition == 3)  // OUT-0x or Cal-02
	{
		char label[3];
		if (sscanf(out, "%c%c%c%c%d%c", &label[0], &label[1], &label[2], &dash, &year, &marker) != 6)
		{
			throw ModelException(routine,"Must be 6 chars for OUT-0x or Cal-02 format");
		}
		else if (dash != '-' || marker != '\n' || !adjustYear(&year))
		{
			throw ModelException(routine,"Some wrong char in OUT-0x or Cal-02 format");
		}

		if (label[0] == 'O' && label[1] == 'U' && label[2] == 'T')
		{
			// withdrawl period always ends in March of the following year
			month = 3;
			year++;
			periodFormat = WITHDRAWL;
		}
		else if (label[0] == 'C' && label[1] == 'A' && label[2] == 'L')
		{
			month = 12;
			periodFormat = ANNUAL;
		}
		else
			throw ModelException(routine,"Something wrong in OUT-0x or Cal-02 format");
		
	}
	else
		throw ModelException(routine,"dash '-' in a wrong position");  // dash in invalid position;

	DateTime outDate = DateTime::MonthDayYear(1,month,year).toDateTime();

	// populate endDate
	endDate = outDate.returnEndOfMonth(false);

	// try to populate startDate...
	int numMonths;

	// number of months to subtract to get starting month of periods
	switch (periodFormat)
	{
		case MONTH:
			numMonths = 0;  // same month
			break;

		case QUARTER:
			numMonths = 2;
			break;

		case SEMIANNUAL:
			numMonths = 5;
			break;

		case ANNUAL:
			numMonths = 11;
			break;

		case INJECTION:
			numMonths = 6;
			break;

		case WITHDRAWL:
			numMonths = 4;
			break;
	}

    if(month<=numMonths)
	{
	    month += 12;
	    year -= 1;
	}
	month -= numMonths;

	startDate = DateTime::MonthDayYear(1,month,year).toDateTime();

	delete [] out;

}

EnergyTermPeriodBenchmark::EnergyTermPeriodBenchmark() : EnergyTermPeriod(TYPE)
{
}

EnergyTermPeriodBenchmark::~EnergyTermPeriodBenchmark()
{
}

bool EnergyTermPeriodBenchmark::adjustYear(int* year)
{

	if (*year>=0 && *year<70)
	{
		*year += 2000;
		return true;
	}
	else if (*year>=70 && *year<100)
	{
		*year += 1900;
		return true;	
	}
	else if (*year>=1900 && *year<2300)
	{
		return true;
	}
	
	return false;
}


/* for reflection */
class EnergyTermPeriodBenchmarkHelper
{

public:

   /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyTermPeriodBenchmark, clazz);
        SUPERCLASS(EnergyTermPeriod);
        EMPTY_SHELL_METHOD(defaultEnergyTermPeriodBenchmark);
        FIELD(benchmark, "Term Benchmark");
		
    }
    
    static IObject* defaultEnergyTermPeriodBenchmark(){
        return new EnergyTermPeriodBenchmark();
    }
};

CClassConstSP const EnergyTermPeriodBenchmark::TYPE = CClass::registerClassLoadMethod(
    "EnergyTermPeriodBenchmark", typeid(EnergyTermPeriodBenchmark), EnergyTermPeriodBenchmarkHelper::load);


DRLIB_END_NAMESPACE
