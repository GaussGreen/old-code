//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : ITimelineSpec.cpp
//
//   Description : 
//
//   Date        : Dec 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/ITimelineSpec.hpp"
#include "edginc/IHasMaturityDate.hpp"

DRLIB_BEGIN_NAMESPACE

/*
IObject* ITimelineSpec::defaultITimelineSpec()
{
	return NULL;
};

*/

static void ITimelineSpecLoad(CClassSP& clazz)
{
    REGISTER_INTERFACE(ITimelineSpec, clazz);
	EXTENDS(IObject);
	clazz->setPublic(); // make visible to EAS/spreadsheet
};

CClassConstSP const ITimelineSpec::TYPE = 
    CClass::registerInterfaceLoadMethod("ITimelineSpec", 
                                        typeid(ITimelineSpec), 
                                        ITimelineSpecLoad);

bool ITimelineSpecload()
{
	return ITimelineSpec::TYPE != NULL;
};


CClassConstSP const TimelineMaturity::TYPE = CClass::registerClassLoadMethod(
    "TimelineMaturity", typeid(TimelineMaturity), TimelineMaturity::load);


IObject* TimelineMaturity::defaultTimelineMaturity()
{
	return new TimelineMaturity();
};

void TimelineMaturity::validatePop2Object()
{}; 

void TimelineMaturity::load(CClassSP& clazz)
{
	clazz->setPublic();

	REGISTER(TimelineMaturity, clazz);
	
	SUPERCLASS(CObject);

	IMPLEMENTS(ITimelineSpec);

	EMPTY_SHELL_METHOD(defaultTimelineMaturity);
};

DateTimeArrayConstSP TimelineMaturity::timeline(const IInstrument &instrument) const
{
	DateTimeArraySP temp = DateTimeArraySP(new DateTimeArray());

	const IHasMaturityDate *gcdo = dynamic_cast<const IHasMaturityDate *> (&instrument);

	if (!gcdo)
		throw ModelException("Not a generalisedCDO");

	temp->push_back(instrument.getValueDate());

	temp->push_back(gcdo->maturityDate());

	return temp;
	
};

///////////////////////////////////////////////


CClassConstSP const TimelineDates::TYPE = CClass::registerClassLoadMethod(
    "TimelineDates", typeid(TimelineDates), TimelineDates::load);


IObject* TimelineDates::defaultTimelineDates()
{
	return new TimelineDates();
};

void TimelineDates::validatePop2Object()
{}; 

void TimelineDates::load(CClassSP& clazz)
{
	clazz->setPublic();

	REGISTER(TimelineDates, clazz);
	
	SUPERCLASS(CObject);

	FIELD(dates, "List of dates");

	IMPLEMENTS(ITimelineSpec);

	EMPTY_SHELL_METHOD(defaultTimelineDates);
};

TimelineDates::TimelineDates(
	DateTimeArray dts,
	CClassConstSP clazz)
	: CObject(clazz)
{
	dates = DateTimeArraySP(new DateTimeArray(dts));
};

DateTimeArrayConstSP TimelineDates::timeline(const IInstrument &instrument) const
{
	return dates;
};

DRLIB_END_NAMESPACE

