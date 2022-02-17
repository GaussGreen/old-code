/**
 * @file VolAtmPointwiseConstConvx.hpp
 */

#ifndef QLIB_VolAtmPointwiseConstConvx_H
#define QLIB_VolAtmPointwiseConstConvx_H

#include "edginc/ExpiryWindow.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

struct RISKMGR_DLL VolAtmPointwiseConstConvx: CObject {
    static CClassConstSP const TYPE;
    VolAtmPointwiseConstConvx(); ~VolAtmPointwiseConstConvx();

    typedef ExpiryWindow Qualifier;

    enum { discrete = 0 };
};

DRLIB_END_NAMESPACE

#endif
