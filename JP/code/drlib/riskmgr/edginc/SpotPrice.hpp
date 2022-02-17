//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotPrice.hpp
//
//   Description : SpotPrice sensitivity
//
//   Author      : Stephen Hope
//
//   Date        : 15 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_SPOT_PRICE_H
#define EDG_SPOT_PRICE_H
#include "edginc/Sensitivity.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE
class SensControlPerName;
class Results;

/** SpotPrice Sensitivity. Effectively does a DELTA shift then reports
    the initial spot price as the sensitivity. This functionality is 
    needed by analytics */
class RISKMGR_DLL SpotPrice: public Sensitivity {
public:
    friend class SpotPriceHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
   
    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns false */
    virtual bool discreteShift() const;

    /** Combines spot prices between results packets (ie does a merge) */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;
    
    /** identifies the name used storing associated results in the output */
    const string& getSensOutputName() const;

    SpotPrice();
protected:
    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*      tweakGroup,
                           Results*         results);
    
private:
    SpotPrice(const SpotPrice &rhs);
    SpotPrice& operator=(const SpotPrice& rhs);

    void subCalculate(const OutputNameArray&        names,
                      IObjectSP                     objToTweak,
                      SensControlPerName&           delta,
                      CBoolArray&                   spotPriceCalcFlag,
                      Results*                      results);
 
};

typedef SpotPrice CSpotPrice;
typedef smartConstPtr<SpotPrice> SpotPriceConstSP;
typedef smartPtr<SpotPrice> SpotPriceSP;
typedef SpotPriceConstSP CSpotPriceConstSP;
typedef SpotPriceSP CSpotPriceSP;


DRLIB_END_NAMESPACE

#endif

