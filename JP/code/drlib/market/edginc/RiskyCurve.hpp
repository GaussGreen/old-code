//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RiskyCurve.hpp
//
//   Description : Interface for risky curves
//
//   Author      : Andrew J Swain
//
//   Date        : 17 November 2003
//
//
//----------------------------------------------------------------------------

#ifndef _IRISKYCURVE_HPP
#define _IRISKYCURVE_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for risky curves */
class MARKET_DLL IRiskyCurve{
public:
    static CClassConstSP const TYPE;

#if 0
    /** this function calculates the discount factor based on the assumption
        that the on which the recovery is based is provided externally. 
        This allows to use different methodologies 
        (PV, face value + accrued etc.) to be included easily */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           double          cashFlow,
                           double          recoveryNotional) const = 0;

    /** this function calculates the discount factor based on the assumption
        that the on which the recovery is based is provided externally. 
        This allows to use different methodologies (PV, face value + accrued etc.)
        to be included easily  - this function will use the externally given
        recovery rather than the underlying risky curve's recovery */
    virtual double riskyPV(const DateTime& lodate,
                           const DateTime& hidate,
                           const double    cashFlow,
                           double          recoveryNotional,
                           bool            useAssetRecovery,
                           double          assetRecovery) const;
#endif

    virtual ~IRiskyCurve();
protected:
    IRiskyCurve();
};


DRLIB_END_NAMESPACE
#endif
