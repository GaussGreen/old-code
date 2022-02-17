/**
 * @file StdDevForLossMapping.hpp
	Author   : Juan Carlos Porras      Dec 06
	Definitions for std dev tweak to loss mapping
 */

#ifndef DRLIB_StdDevForLossMapping_H
#define DRLIB_StdDevForLossMapping_H

#include "edginc/Void.hpp"


DRLIB_BEGIN_NAMESPACE

 
struct RISKMGR_DLL StdDevForLossMapping : CObject {

    static CClassConstSP const TYPE;

	StdDevForLossMapping(); ~StdDevForLossMapping(); 

    typedef Void Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
