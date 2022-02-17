//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesBootstrapperTimeOnlyTimeOnly.cpp
//
//   Description : CDOQuotesBootstrapperTimeOnlyTimeOnly defines an 'iterator' over CDO quotes. 
//					It differs from CDOQuotesBootstrapperTimeOnly in that it is for boostrapping 
//					in the time dimension only (not strike).
//
//   Author      : Matthias Arnsdorf
//
//   Date        : September 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CDOQuotesBootstrapperTimeOnly.hpp"
#include "edginc/CDO.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/TimePoint2D.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE

/** Virtual destructor */
CDOQuotesBootstrapperTimeOnly::~CDOQuotesBootstrapperTimeOnly() {}

/** Constructor */
CDOQuotesBootstrapperTimeOnly::CDOQuotesBootstrapperTimeOnly(CDOQuotesConstSP cdoQuotes) :
CDOQuotesBootstrapper(cdoQuotes, false), // do not ignore 100pc point
bootstrapMaturities(0),
currentMatIdx(0),
instruments(0)
{
	static const string method = "CDOQuotesBootstrapperTimeOnly::CDOQuotesBootstrapperTimeOnly";

	// construct array of unique maturities
	bootstrapMaturities = DateTimeArraySP(new DateTimeArray(*CDOQuotesBootstrapper::maturities));
	DateTime::removeDuplicates(*bootstrapMaturities, true);

	// construct instruments collection
	int numDates = bootstrapMaturities->size();
	instruments.resize(numDates);
	int i = 0;
	for(; i < numDates; i++)
	{
		instruments[i] = CInstrumentArraySP(new CInstrumentArray());
	}
	
	int d = 0;
	DateTime date;
	// loop through all instruments
	for (CDOQuotesBootstrapper::init(); !CDOQuotesBootstrapper::end(); CDOQuotesBootstrapper::next()) 
	{
#ifdef DEBUG
		if(d >= (int) CDOQuotesBootstrapper::maturities->size())
		{
			throw ModelException("Index d is out of range", method);
		}
#endif
		date = (*CDOQuotesBootstrapper::maturities)[d];
		// generate CDO instruments and add to instruments
		instruments[date.find(*bootstrapMaturities)]->push_back(
			CDOQuotesBootstrapper::buildCurrentInstrument()
			);

		d++;
	}


}

/** Constructor (internal - used by reflection) */
CDOQuotesBootstrapperTimeOnly::CDOQuotesBootstrapperTimeOnly() :
CDOQuotesBootstrapper(CDOQuotesConstSP(0), false),
bootstrapMaturities(0),
currentMatIdx(0),
instruments(0)
{}

/**
* Method called before first step of the loop
* [Implements IBootstrapper]
* */
void CDOQuotesBootstrapperTimeOnly::init() 
{
	currentMatIdx = 0;
}

/**
* Method called after each step of the loop (ie go to the next quote)
* [Implements IBootstrapper]
* */
void CDOQuotesBootstrapperTimeOnly::next() {
	currentMatIdx++;
}

/**
* Method called to test the end of the loop (returns true when there is no more quote)
* [Implements IBootstrapper]
* */
bool CDOQuotesBootstrapperTimeOnly::end() const 
{
	return (currentMatIdx >= bootstrapMaturities->size());
}

/**
* Returns "state" corresponding to current step of the loop
* this is just the current maturity 
* [Implements IBootstrapper]
* */
IObjectSP CDOQuotesBootstrapperTimeOnly::getCurrentState() const 
{
	return DateTimeSP(new DateTime((*bootstrapMaturities)[currentMatIdx]));
}

/** Returns the instruments corresponding to the current state */
CInstrumentArraySP CDOQuotesBootstrapperTimeOnly::buildCurrentInstruments() {
	static string method = "CDOQuotesBootstrapperTimeOnly::buildCurrentInstruments";
	try {
		if (end()) {
			// no more quotes
			throw ModelException(method,
				"Internal error: no more quotes to build instrument.");
		} else 
		{
			if(instruments.size() == 0)
			{
				throw ModelException("instruments array has 0 size");
			}
			if(currentMatIdx >= (int)instruments.size())
			{
				throw ModelException("currentMatIdx out of bounds");
			}
			if(!instruments[currentMatIdx])
			{
				throw ModelException("instruments are not initialised");
			}
			return instruments[currentMatIdx];
		}
	} catch (exception& e){
		throw ModelException(e, method);
	}
}

/** Invoked when Class is 'loaded' */
void CDOQuotesBootstrapperTimeOnly::load(CClassSP& clazz) {
	clazz->setPrivate(); // don't make visible to EAS/spreadsheet
	REGISTER(CDOQuotesBootstrapperTimeOnly, clazz);
	SUPERCLASS(CDOQuotesBootstrapper);
	EMPTY_SHELL_METHOD(defaultConstructor);
	FIELD_NO_DESC(bootstrapMaturities);
	FIELD_NO_DESC(currentMatIdx);
}

/** Default constructor */
IObject* CDOQuotesBootstrapperTimeOnly::defaultConstructor() {
	return new CDOQuotesBootstrapperTimeOnly();
}

/** TYPE (for reflection) */
CClassConstSP const CDOQuotesBootstrapperTimeOnly::TYPE =
CClass::registerClassLoadMethod(
								"CDOQuotesBootstrapperTimeOnly",
								typeid(CDOQuotesBootstrapperTimeOnly),
								CDOQuotesBootstrapperTimeOnly::load);

DRLIB_END_NAMESPACE

