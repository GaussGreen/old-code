//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : ITimelineSpec.hpp
//
//   Description : 
//
//   Date        : Dec 2006
//
//----------------------------------------------------------------------------

#ifndef ITIMELINESPEC_HPP
#define ITIMELINESPEC_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

class MCARLO_DLL ITimelineSpec	: public virtual IObject
{

public:

	static CClassConstSP const TYPE;
	
	virtual DateTimeArrayConstSP timeline(const IInstrument &instrument) const = 0;
//	{
//		throw ModelException("timeline function undefined for interface");
//	}

	//static void load(CClassSP& clazz);

//	ITimelineSpec() {};

//	~ITimelineSpec() {};

	// static IObject* defaultITimelineSpec();

};

DECLARE(ITimelineSpec);

class MCARLO_DLL TimelineMaturity : public CObject,
						 public virtual ITimelineSpec
{

public:

	static CClassConstSP const TYPE;

	virtual DateTimeArrayConstSP timeline(const IInstrument &instrument) const;

	static void load(CClassSP& clazz);

	static IObject* defaultTimelineMaturity();

	virtual void validatePop2Object(); 

	TimelineMaturity(CClassConstSP clazz = TYPE) : CObject(clazz) {};

	~TimelineMaturity(){};
};

DECLARE(TimelineMaturity);

class MCARLO_DLL TimelineDates : public CObject,
						 public virtual ITimelineSpec
{

	DateTimeArraySP dates;

public:

	static CClassConstSP const TYPE;

	virtual DateTimeArrayConstSP timeline(const IInstrument &instrument) const;

	static void load(CClassSP& clazz);

	static IObject* defaultTimelineDates();

	virtual void validatePop2Object(); 

	TimelineDates(CClassConstSP clazz = TYPE) : CObject(clazz) {};

	TimelineDates(DateTimeArray dts,
				  CClassConstSP clazz = TYPE) ;

	~TimelineDates(){};
};

DECLARE(TimelineDates);


DRLIB_END_NAMESPACE

#endif
