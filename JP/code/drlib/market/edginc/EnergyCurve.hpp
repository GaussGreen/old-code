//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyCurve.hpp
//
//   Description : Base for all Energy Curves. Based on drcc.h and
//                 drcc.cpp in FXLIB.
//
//   Author      : Sean Chen
//
//   Date        : April 18, 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergyCurve_
#define _EnergyCurve_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/EnergyUnderlyer.hpp"

#include <string>
using namespace std;


DRLIB_BEGIN_NAMESPACE

class EnergyCurve;
typedef smartPtr<EnergyCurve> EnergyCurveSP;
typedef smartConstPtr<EnergyCurve> EnergyCurveConstSP;

class MARKET_DLL EnergyCurve : public MarketObject
{

public:

    static CClassConstSP const TYPE;


    virtual ~EnergyCurve();

    virtual DateTime getToday() const;
    virtual string getName() const;
    
protected:
    
    EnergyCurve(const CClassConstSP& clazz);

    static void load(CClassSP& clazz);


    DateTime             today; // $unregistered
    string                 name; // $unregistered
private:
    EnergyCurve(const EnergyCurve& v);
};


typedef MarketWrapper<EnergyCurve> EnergyCurveWrapper;

DRLIB_END_NAMESPACE

#endif
