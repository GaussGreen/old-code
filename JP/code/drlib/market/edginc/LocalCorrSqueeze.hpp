//----------------------------------------------------------------------------
//
//   Group       : QR Equities London
//
//   Filename    : LocalCorrSqueeze.hpp
//
//   Description : Holds market data for local correlation skew model
//
//   Author      : Eva Strasser
//
//----------------------------------------------------------------------------
#ifndef LOC_CORR_SQUEEZE_HPP
#define LOC_CORR_SQUEEZE_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/LocalCorrExpiryAndStrikewise.hpp"
#include "edginc/LocalCorrExpiry.hpp"
#include "edginc/LocalCorrVoid.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for market data for local correlation skew model */
class MARKET_DLL LocalCorrSqueeze: 
    public MarketObject,
    virtual public ITweakableWithRespectTo<LocalCorrExpiryAndStrikewise>,
    virtual public ITweakableWithRespectTo<LocalCorrExpiry>, 
    virtual public ITweakableWithRespectTo<LocalCorrVoid> {

public:
    static CClassConstSP const TYPE;

    static const string LOCAL_CORR_SQUEEZE_STEP;
    static const string LOCAL_CORR_SQUEEZE_LINEAR;

    ~LocalCorrSqueeze();

    virtual void validatePop2Object();
    
    /** dont plan to override these methods, thus not virtual */
    string getName() const;    
    double getMinSqueeze() const;
    double getMaxSqueeze() const;
    bool isZeroObject() const;

    double computeSqueeze(double marketFactor,
                          double tradingTime) const;
    
    virtual void getMarket(const IModel*        model, 
                           const MarketData*    market);
    
    /** for construction of dummy object */
    LocalCorrSqueeze(bool isZeroObj);
    
    static IObject* defaultLocalCorrSqueeze();

    /** Support for TWEAK sensitivity */
    virtual string sensName(const LocalCorrExpiryAndStrikewise*) const;
    virtual ExpiryAndStrikeArrayConstSP sensQualifiers(const LocalCorrExpiryAndStrikewise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<LocalCorrExpiryAndStrikewise>&);

    /** Support for POINTWISE sensitivity */
    virtual string sensName(const LocalCorrExpiry*) const;
    virtual ExpiryWindowArrayConstSP sensQualifiers(const LocalCorrExpiry*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<LocalCorrExpiry>&);

    /** Support for PARALLEL sensitivity */
    virtual string sensName(const LocalCorrVoid*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<LocalCorrVoid>&);

protected:    
    LocalCorrSqueeze();    
    
    static void load(CClassSP& clazz);

    /** registered fields */
    string              name;
    ExpiryArraySP       expiries;
    CDoubleArraySP      strikes;
    CDoubleMatrixSP     squeezes;    
    string              interpolationType;

    /** transient fields */
    double  minSqueeze; // min squeeze across whole matrix
    double  maxSqueeze; // max squeeze across whole matrix
    bool    isZeroObj;  // boolean to indicate whether or not dummy object
        
    DoubleArraySP   expiriesInTradingTime; // expiry converted into trading time using dummy metric
    DoubleArraySP   strikeIndices; // helper for sensitivities
    /** asd */
    int                 interpolationTypeInt; // to be converted into ENUMS
            
    /** only used in case interpolation type is LINEAR */
    InterpolantArray linearInterpSqueeze;
};
typedef smartPtr<LocalCorrSqueeze> LocalCorrSqueezeSP;
typedef array<LocalCorrSqueezeSP, LocalCorrSqueeze> LocalCorrSqueezeArray;

#ifndef LOC_CORR_SQUEEZE_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<LocalCorrSqueezeSP _COMMA_ LocalCorrSqueeze>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<LocalCorrSqueezeSP _COMMA_ LocalCorrSqueeze>);
#endif

DRLIB_END_NAMESPACE
#endif // LOC_CORR_SQUEEZE_HPP
