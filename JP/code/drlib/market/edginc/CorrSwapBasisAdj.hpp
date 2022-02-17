//----------------------------------------------------------------------------
//
//   Group       : QR Equities London
//
//   Filename    : CorrSwapBasisAdj.hpp
//
//   Description : 
//----------------------------------------------------------------------------
#ifndef QLIB_CorrSwapBasisAdj_HPP
#define QLIB_CorrSwapBasisAdj_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/Class.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/CorrSwapBasisAdjAbsolute.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CorrSwapBasisAdj: public MarketObject,
                                   virtual public ITweakableWithRespectTo<CorrSwapBasisAdjAbsolute> {
public:
    static CClassConstSP const TYPE;
	virtual ~CorrSwapBasisAdj();

	friend class CorrSwapBasisAdjHelper;

    /** Constructor */
    CorrSwapBasisAdj();    
    
    /** Constructor for zero object */
    CorrSwapBasisAdj(bool isZeroObj);    

    /** Validation */
    virtual void validatePop2Object();

    /** Override default implementation in order to set up interpolator */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Get squeeze for certain expiry */ 
	virtual double getSqueezeForExpiry(double time) const;
    
    /** Returns the CorrSwapBasisAdj's name, which corrsponds to the region, for simplicity */
    virtual string getName() const;

    /** Tweak (sensitivity) support */
    virtual string sensName(const CorrSwapBasisAdjAbsolute* name) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CorrSwapBasisAdjAbsolute>&);
    
private:
    string              region;					// mandatory
    ExpiryArraySP		expiries;				// mandatory
	DoubleArraySP		squeezes;				// mandatory

    // transient fields
	LinearInterpolantNonVirtualSP corrSwapBasisInterp;
    DoubleArraySP basisTradYears;
    bool isZeroObj;
};

typedef smartPtr<CorrSwapBasisAdj> CorrSwapBasisAdjSP;
typedef smartConstPtr<CorrSwapBasisAdj> CorrSwapBasisAdjConstSP;

DRLIB_END_NAMESPACE
#endif // CorrSwapBasisAdj_HPP
