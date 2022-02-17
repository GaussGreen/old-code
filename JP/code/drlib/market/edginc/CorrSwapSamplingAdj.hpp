//----------------------------------------------------------------------------
//
//   Group       : QR Equities London
//
//   Filename    : CorrSwapSamplingAdj.hpp
//
//   Description : 
//----------------------------------------------------------------------------
#ifndef QLIB_CorrSwapSamplingAdj_HPP
#define QLIB_CorrSwapSamplingAdj_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/Class.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/CorrSwapSamplingAdjAbsolute.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CorrSwapSamplingAdj:  public MarketObject,
                            virtual public ITweakableWithRespectTo<CorrSwapSamplingAdjAbsolute> {
public:
    static CClassConstSP const TYPE;
	virtual ~CorrSwapSamplingAdj();

	friend class CorrSwapSamplingAdjHelper;
    
    /** Constructor */
    CorrSwapSamplingAdj();    
    
    /** Constructor for zero object */
    CorrSwapSamplingAdj(bool isZeroObj);    

    /** Validation */
    virtual void validatePop2Object();

    /** Returns the squeeze */ 
    virtual double getSqueeze() const;	
    
    /** Returns the CorrSwapSamplingAdj's name */
    virtual string getName() const;

    /** Tweak (sensitivity) support */
    virtual string sensName(const CorrSwapSamplingAdjAbsolute*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CorrSwapSamplingAdjAbsolute>&);

    /** Initialises this piece of market data - 
        records the pair of names idenitfying the correlation */
    virtual void initialise(MarketData* market);
    
private: 
    string              name;                   // optional
    string              region1;				// mandatory
	string				region2;				// mandatory
	double				squeeze;

    /** transient fields */
    bool isZeroObj;
};

typedef smartPtr<CorrSwapSamplingAdj> CorrSwapSamplingAdjSP;
typedef smartConstPtr<CorrSwapSamplingAdj> CorrSwapSamplingAdjConstSP;

DRLIB_END_NAMESPACE
#endif // CorrSwapSamplingAdj_HPP
