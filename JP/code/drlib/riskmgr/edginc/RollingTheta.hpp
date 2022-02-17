//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RollingTheta.hpp
//
//   Description : Rolling Theta shift
//
//   Author      : André Segger
//
//   Date        : 05 October 2001
//
//
//----------------------------------------------------------------------------


#ifndef ROLLING_THETA_HPP
#define ROLLING_THETA_HPP
#include "edginc/Interval.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/SensControl.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Theta shift */
class RISKMGR_DLL RollingTheta: public Sensitivity,
                    virtual public Additive {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static int    DEFAULT_SHIFT;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns true */
    virtual bool discreteShift() const;

    virtual void validatePop2Object();

    /** identifies the name used for storing associated results in the output*/
    virtual const string& getSensOutputName() const;

    /** identifies the packet in which the results are stored. Theta
        results are stored in the instrument packet */
    virtual const string& getPacketName() const;

    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*  tweakGroup,
                           CResults*    results);

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

protected:
    RollingTheta(const CClassConstSP& clazz,
                 const string&        outputName);
private:
    IntervalArraySP  thetaInterval;
    double           deltaShift;
    int              offset;
    bool             useAssetFwds; // optional
    HolidayWrapper   hols;

    /** for reflection */
    RollingTheta();
    RollingTheta(const RollingTheta &rhs);
    RollingTheta& operator=(const RollingTheta& rhs);
    friend class RollingThetaHelper;
};

typedef smartConstPtr<RollingTheta>             RollingThetaConstSP;
typedef smartPtr<RollingTheta>                  RollingThetaSP;


DRLIB_END_NAMESPACE

#endif
