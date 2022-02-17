//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditDebug.hpp
//
//   Description : Debug info for credit
//
//   Author      : Jay Blumenstein
//
//   Date        : 18 Sep 2002
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_DEBUG_HPP
#define CREDIT_DEBUG_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** simulated path values */
class CREDIT_DLL CreditDebug: public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
	
    static IObject* defaultCreditDebug(){
        return new CreditDebug();
    }

	CreditDebug();
	
    smartPtr< DateTimeArrayArray > pathDates; //path dates
    DoubleArrayArrayArraySP pathSpots; // spots on each path date
    DoubleArrayArrayArraySP pathValues; // values on each path date
    DoubleArrayArraySP exposureWithMargin; // [path][exp date] = exposure after margin is considered.
	double calcTime;
};

typedef smartConstPtr<CreditDebug> CreditDebugConstSP;
typedef smartPtr<CreditDebug> CreditDebugSP;

typedef array<CreditDebugSP, CreditDebug> CreditDebugArray;

DRLIB_END_NAMESPACE

#endif
