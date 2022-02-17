//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarginAcct.cpp
//
//   Description : Parameters for a Margin Collateral Account
//
//   Author      : Jay Blumenstein
//
//   Date        : 13 Sep 2002
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/MarginAcct.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////
//
//  class MarginAcct
//
//  contains: All Credit inputs
//  methods: compute grids
//////////////////////////////////////////////////////////////////////

MarginAcct :: MarginAcct() : CObject(TYPE)
{
	marginCallHol = HolidaySP( Holiday::noHolidays() );
}

CClassConstSP const MarginAcct::TYPE = CClass::registerClassLoadMethod(
							"MarginAcct", typeid(MarginAcct), load);

void MarginAcct::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(MarginAcct, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultMarginAcct);
    FIELD(upfrontMargin, " var. margin in account already");
    FIELD(initialMargin, " init margin acct. Never repaid");
    FIELD(thresholdMarginCpty, " unsecured margin counter party");
    FIELD(thresholdMarginJPM, " unsecured margin JPM is never asked for");
    FIELD(minMarginCall, " (applied after threshold)");
    FIELD(marginCallDelay, " days btn margin call and receipt");
    FIELD(rehypothecate, " if counter party can use margin (only TRUE allowed currently)");
	FIELD(marginCallHol, "optional field: used to find the last MTM date for an exposure date.");
	FIELD_MAKE_OPTIONAL(marginCallHol);    
       
}

void MarginAcct::validatePop2Object(){
    const string method = "MarginAcct::validatePop2Object";

	if( !rehypothecate )
		throw ModelException(method, "Currently only TRUE allowed for rehypothecate");
}

DRLIB_END_NAMESPACE
