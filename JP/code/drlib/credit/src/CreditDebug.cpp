//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditDebug.cpp
//
//   Description : Debug info for credit
//
//   Author      : Jay Blumenstein
//
//   Date        : 18 Sep 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditDebug.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const CreditDebug::TYPE = CClass::registerClassLoadMethod(
							"CreditDebug", typeid(CreditDebug), load);

DEFINE_TEMPLATE_TYPE(CreditDebugArray);

CreditDebug::CreditDebug() : CObject(TYPE), 
                    pathDates(smartPtr< DateTimeArrayArray >(0)), 
                    pathSpots(DoubleArrayArrayArraySP(   )), 
                    pathValues(DoubleArrayArrayArraySP(   )), 
                    exposureWithMargin(DoubleArrayArraySP(   )),
                    calcTime(0.0) {}


void CreditDebug::load(CClassSP& clazz)
{
        REGISTER(CreditDebug, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCreditDebug);
        clazz->setPublic(); // make visible to EAS/spreadsheet
		FIELD(pathDates, "path dates");
		FIELD(pathSpots, "spots on each path date");
		FIELD(pathValues, "values on each path date");
        FIELD(exposureWithMargin, "exposure for each path and exposure date, after simulating margin acct");
		FIELD(calcTime, "time for computation of calcPeakAndAverage.");
}

DRLIB_END_NAMESPACE
