//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : DecretionCurve.hpp
//
//   Description : A decretion curve base for ABS
//
//   Author      : Keith Jia
//
//   Date        : 03 January 2006
//
//
//----------------------------------------------------------------------------

#ifndef DECRETIONCURVE_HPP
#define DECRETIONCURVE_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/IDecretionCurve.hpp"

DRLIB_BEGIN_NAMESPACE

/** No object of DecrectionCurve can be created.  User should create objects
 ** of derived classes
 */

class MARKET_DLL DecretionCurve : public MarketObject,
                       virtual public IDecretionCurve 
{
public:
    static CClassConstSP const TYPE;
    virtual ~DecretionCurve();

    //------------------------------------------
    // MarketObject methods
    //------------------------------------------
    /** Returns a name */
    string getName() const;

protected:
    DecretionCurve(CClassConstSP clazz = TYPE);
    string name;

private:
    static void load(CClassSP& clazz);
};

typedef smartConstPtr<DecretionCurve> DecretionCurveConstSP;
typedef smartPtr<DecretionCurve>      DecretionCurveSP;
typedef MarketWrapper<DecretionCurve> DecretionCurveWrapper;

DRLIB_END_NAMESPACE
#endif
