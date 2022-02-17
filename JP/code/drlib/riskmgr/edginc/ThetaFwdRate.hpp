//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ThetaFwdRate.hpp
//
//   Description : Theta Forward Rate - preserve forward zero rates (df's)
//
//   Author      : Ning Shen
//
//   Date        : 5 Nov 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_THETA_FWD_RATE_HPP
#define EDR_THETA_FWD_RATE_HPP
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Theta fwd rate shift */
class RISKMGR_DLL ThetaFwdRate: public Theta{
public:
    static CClassConstSP const TYPE;
    static const string NAME;

    /** Public constructor */
    ThetaFwdRate(int offset, HolidaySP hols);

private:

    ThetaFwdRate();
    ThetaFwdRate(const ThetaFwdRate &rhs);
    ThetaFwdRate& operator=(const ThetaFwdRate& rhs);

    static void load(CClassSP& clazz);
    static IObject* defaultThetaFwdRate();
};

typedef smartConstPtr<ThetaFwdRate> ThetaFwdRateConstSP;
typedef smartPtr<ThetaFwdRate> ThetaFwdRateSP;


DRLIB_END_NAMESPACE

#endif
