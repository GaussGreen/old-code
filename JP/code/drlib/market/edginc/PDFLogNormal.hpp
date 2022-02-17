//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFLogNormal.cpp
//
//   Description : generate model pdf for log-normal
//
//   Author      : Ning Shen
//
//   Date        : 15 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef PDFLogNormal_HPP
#define PDFLogNormal_HPP

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/Model.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolProcessedBS.hpp"

DRLIB_BEGIN_NAMESPACE

/** log-normal model pdf is left here, as it's only used for testing */
class MARKET_DLL PDFLogNormal : public PDFCalculator
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz){
        REGISTER(PDFLogNormal, clazz);
        SUPERCLASS(PDFCalculator);
    }

    virtual ~PDFLogNormal(){};

    // just for compilation now
    PDFLogNormal(const DateTime& valueDate,     // value date
                 const CAsset* asset,
                 const CVolProcessed* vol) 
                 : PDFCalculator(TYPE),
                 valueDate(valueDate), asset(asset){
        volBS = dynamic_cast<const CVolProcessedBS*>(vol);
        if (!volBS)
            throw ModelException("PDFLogNormal", "failed to construct the object, BS volrequired");
    }

	/** Calculate the cumulative probability at each strike.
        Each prob is estimated as a "small" call spread (i.e., by tweaking the strike)
        Returns false if if the difference between 2 consecutive probs is negative.
    */
    virtual void probabilities(const DoubleArray& strikes,
                               const DateTime&    maturity,
                               DoubleArray&       probs) const;

    virtual void probabilities(const CLatticeDouble&   strikes,
                               const DateTimeArray&    maturities,
                               CLatticeDouble&         probs) const
    {
        // do nothing for now
    }

    /** Calculate the "local" density at each strike.
        Each density is estimated as a "small" butterfly (using same strikes 
        as above together with the center strike).
        Returns false if any of the densities is negative
    */
    virtual void localDensity(const DoubleArray& strikes,
                              const DateTime&    maturity,
                              DoubleArray&       density) const
    {
        // do nothing for now
    }

    /** Calculate the "integrated" density at each strike.
        Each density is computed by taking the difference between 2 consecutive
        probabilities divided by the difference between the 2 consecutive strikes.
        false should be returned if any of the so-computed densities is negative;
    */        
    virtual void integratedDensity(const DoubleArray& strikes,
                                   const DateTime&    maturity,
                                   DoubleArray&       density) const
    {
        // do nothing for now
    }

private:
	const DateTime				valueDate; // $unregistered
	const CAsset*               asset; // $unregistered
    const CVolProcessedBS*           volBS; // $unregistered
};

DRLIB_END_NAMESPACE

#endif
