//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyFuture.hpp
//
//   Description : Base for all Energy Curves. Based on drcommodityfuture.h and
//                 drcommodityfuture.cpp in FXLIB.
//
//   Author      : Sean Chen
//
//   Date        : April 18, 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergyFuture_
#define _EnergyFuture_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/ClosedFormEnergy.hpp"

#include "edginc/Instrument.hpp"
#include "edginc/Addin.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/SensControl.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/FirstSensDate.hpp"


#include <string>
using namespace std;


DRLIB_BEGIN_NAMESPACE

class EnergyFuture;
typedef smartPtr<EnergyFuture> EnergyFutureSP;
typedef smartConstPtr<EnergyFuture> EnergyFutureConstSP;
class DateTime;

class PRODUCTS_DLL EnergyFuture :    public CInstrument, 
                        public virtual ClosedFormEnergy::IIntoProduct,
                        public LastSensDate,
                        virtual public MonteCarlo::IIntoProduct,
			            public FirstSensDate
{

public:

    static CClassConstSP const TYPE;

    friend class EnergyFutureHelper;
    friend class EnergyFutureMC;
    virtual ~EnergyFuture();
    
//    EnergyFuture& operator=(const EnergyFuture& v);
    
    string getName() const;
    string discountYieldCurveName() const;
    
    double getPrice() const;
    DateTime getValueDate() const;
    double getNumContracts() const;
    string getBenchmark() const;
    
    virtual DateTime getValuedDate(const SensControl* sensControl) const;

    virtual DateTime endDate(const Sensitivity* sensControl) const;
    virtual ExpirySP getExpiry(const SensControl* sensControl) const;
    virtual DateTime beginDate(const SensControl* sensControl) const;

    EnergyFuturesCurveConstSP getEnergyFuturesCurve() const;
    EnergyUnderlyerConstSP getEnergyUnderlyer() const;
    
    virtual void GetMarket(const IModel* theMarket, const CMarketDataSP theSP);
    virtual void Validate(){}
    void validatePop2Object();
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;
    ClosedFormEnergy::IProduct* createProduct(ClosedFormEnergy* model) const;

private:
    
    EnergyFuture();
    
    static void load(CClassSP& clazz);

    string                name;
    double                price;
    DateTime             maturityDate;
    string                 benchmark;
    double                numContracts;
    string                 description;
    DateTime             valueDate;
    EnergyUnderlyerWrapper        energyUnderlyerWrapper;
    EnergyFuturesCurveWrapper    curveWrapper;
};


typedef MarketWrapper<EnergyFuture> EnergyFutureWrapper;

DRLIB_END_NAMESPACE

#endif

