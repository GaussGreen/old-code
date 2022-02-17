//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Theta.hpp
//
//   Description : Theta Forward Spot shift
//
//   Author      : Mark A Robson
//
//   Date        : 26 June 2001
//
//
//----------------------------------------------------------------------------


#ifndef EDR_THETA_FWD_SPOT_HPP
#define EDR_THETA_FWD_SPOT_HPP
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Theta shift */
class RISKMGR_DLL ThetaFwdSpot: public Theta{
public:
    static CClassConstSP const TYPE;
    static const string NAME;

    /** Returns true  */
    virtual bool useAssetFwds() const;

    /** Public constructor */
    ThetaFwdSpot(int offset, HolidaySP hols);

private:
    // no additional fields
    /** for reflection */
    ThetaFwdSpot();
    ThetaFwdSpot(const ThetaFwdSpot &rhs);
    ThetaFwdSpot& operator=(const ThetaFwdSpot& rhs);
    friend class ThetaHelper;
    static void load(CClassSP& clazz);
    static IObject* defaultThetaFwdSpot();
};

typedef smartConstPtr<ThetaFwdSpot> ThetaFwdSpotConstSP;
typedef smartPtr<ThetaFwdSpot> ThetaFwdSpotSP;

DRLIB_END_NAMESPACE

#endif
