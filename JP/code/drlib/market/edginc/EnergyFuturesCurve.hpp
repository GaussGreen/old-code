//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyFuturesCurve.hpp
//
//   Description : Energy Fututes Curve. Based on drstdcc.h and
//                 drstdcc.cpp in FXLIB.
//
//   Author      : Sean Chen
//
//   Date        : April 18, 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergyFuturesCurve_
#define _EnergyFuturesCurve_

#include "edginc/EnergyCurve.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/LinearInterpolator.hpp"

#include "edginc/MarketFactor.hpp"

#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/EnergyFuturesCurveParallel.hpp"
#include "edginc/PointwiseTweakableWith.hpp"
#include "edginc/EnergyFuturesCurvePointwise.hpp"
//#include "edginc/EnergyImpliedVolSurface.hpp"
#include "edginc/EnergyInstVolBase.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyFuturesCurve :  public EnergyCurve, 
                            public virtual IGetMarket,
                            virtual public ITweakableWithRespectTo<EnergyFuturesCurveParallel>,
                            virtual public ITweakableWithRespectTo<EnergyFuturesCurvePointwise>,
                            public virtual IMarketFactor
{
public:
    
    friend class EnergyFuturesCurveHelper;
    friend class GetEnergyFuturesCurveAddin;

    static CClassConstSP const TYPE;

    ~EnergyFuturesCurve();

    virtual void validatePop2Object();
    void buildLinearInterpolant();

    double interpolate(const DateTime& theContractDate) const;
    double fixing(const DateTime& theContractDate) const;

    double interpolate(const string& theContractLabel) const;
    double fixing(const string& theContractLabel) const;

    EnergyUnderlyerConstSP getEnergyUnderlyer() const;
	EnergyInstVolBaseConstSP getEnergyInstVolBase() const;

	 // Parallel tweak (sensitivity) support
    virtual string sensName(const EnergyFuturesCurveParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<EnergyFuturesCurveParallel>&);

    // Pointwise tweak (sensitivity) support
    virtual string sensName(const EnergyFuturesCurvePointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<EnergyFuturesCurvePointwise>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const EnergyFuturesCurvePointwise*) const;

    string getName() const;

    void getMarket(const IModel* model, const MarketData* market);

	//int getSize() const;

    // string getExpiryLabel(int i) const;
    const DateTimeArray& getFutureMaturityDates() const;
    double getFwdRate(int i) const;

	double getFwdRateForLabel(const string& label) const;

private:
    EnergyFuturesCurve();

 
    //string    name;
    EnergyUnderlyerWrapper    energyUnderlyerWrapper;
	EnergyInstVolBaseWrapper energyInstVolBaseWrapper; //optional, only use by diffusion engine
    mutable DateTimeArray futureMaturityDates;
    CStringArray expiryLabels;
    CDoubleArray rates ;
	CDoubleArray weights ;
    CDoubleArray expirysInDouble;
    Interpolator::InterpolantSP interpolantSP;

    // for tweaking ....
    ExpiryArraySP     expiries;
    //CDoubleArray    origRates;
	//CDoubleArray    shockedRates;
    //bool shocked;
    //Interpolator::InterpolantSP interpolantOrigSP;
    //bool tweaked;
    //double savedRate;
    //int savedI;
};

typedef MarketWrapper<EnergyFuturesCurve> EnergyFuturesCurveWrapper;
DECLARE(EnergyFuturesCurve);


DRLIB_END_NAMESPACE

#endif
