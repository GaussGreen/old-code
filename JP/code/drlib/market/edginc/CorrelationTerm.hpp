//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrelationTerm.hpp
//
//   Description : Holds implied correlation parameters for term structure
//
//   Author      : Oliver Brockhaus
//
//   Date        : 05 Dec 2002
//
//
//----------------------------------------------------------------------------
#ifndef CORRELATIONTERM_HPP
#define CORRELATIONTERM_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/Class.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/ShortTermSqueezeTweak.hpp"
#include "edginc/LongTermSqueezeTweak.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/Correlation.hpp"

DRLIB_BEGIN_NAMESPACE

class CorrelationTerm;

typedef smartPtr<CorrelationTerm> CorrelationTermSP;
typedef smartConstPtr<CorrelationTerm> CorrelationTermConstSP;
#ifndef QLIB_CORRELATIONTERM_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CorrelationTerm>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CorrelationTerm>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CorrelationTerm>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CorrelationTerm>);
#endif
// support for arrays of correlationTerms (note array of smart pointers)
typedef array<CorrelationTermSP, CorrelationTerm> CorrelationTermArray;
#ifndef QLIB_CORRELATIONTERM_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<CorrelationTermSP _COMMA_ CorrelationTerm>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CorrelationTermSP _COMMA_ CorrelationTerm>);
#endif
// support for smart pointer to array
typedef smartConstPtr<CorrelationTermArray> CorrelationTermArrayConstSP;
typedef smartPtr<CorrelationTermArray> CorrelationTermArraySP;
#ifndef QLIB_CORRELATIONTERM_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CorrelationTermArray>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CorrelationTermArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CorrelationTermArray>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CorrelationTermArray>);
#endif

// place holder for corr term data
struct MARKET_DLL CorrTermData {
    /** Constructor */
    CorrTermData(ExpiryArraySP      benchmarkTermExpiries,
                 ExpiryArraySP      shortTermExpiries,
                 ExpiryArraySP      longTermExpiries,
                 CDoubleMatrixSP    correlationMatrix,
                 CDoubleMatrixSP    shortTermSqueeze,
                 CDoubleMatrixSP    longTermSqueeze);

    ExpiryArraySP   benchmarkTermExpiries;
    ExpiryArraySP   shortTermExpiries;
    ExpiryArraySP   longTermExpiries;
    CDoubleMatrixSP correlationMatrix;
    CDoubleMatrixSP shortTermSqueeze;
    CDoubleMatrixSP longTermSqueeze;
};
typedef refCountPtr<CorrTermData> CorrTermDataSP;



class MARKET_DLL CorrelationTerm:  public MarketObject,
                        public virtual ShortTermSqueezeTweak::IShift,
                        public virtual LongTermSqueezeTweak::IShift
{
public:
    static CClassConstSP const TYPE;
    friend class CorrelationTermHelper;

    static const string SHORT_TERM; 
    static const string LONG_TERM; 
    static const int TERM_TIME;
    static const double EIGEN_VALUE_FLOOR;
    static const double MAX_SQ_ERROR;
   
    /** Used to indicate the type of calculation should be performed when
        an array of dates is passed to CorrelationTermMatrix;
        Currently hard-coded in MultiAsset to "forward" */
    typedef enum _TCalcType{
        fromFirst = 0, // calculate correlation from first date
        forward,       // calculate correlation between successive dates
        toLast         // calculate correlation to last date
    } TCalcType;

    virtual ~CorrelationTerm();

    /** Validation */
    virtual void validatePop2Object();

    virtual void getMarket(const IModel* model, const MarketData* market);

    /** returns squeeze & expiry */
    virtual double      getCorrShortTermSqueeze() const;
    virtual ExpirySP    getShortTermExpiry() const;
    virtual double      getCorrLongTermSqueeze() const;
    virtual ExpirySP    getLongTermExpiry() const;
    
    /** returns the correlationTerm's name */
    virtual string getName() const;

    /** Initialises this piece of market data - 
        records the pair of names idenitfying the correlation */
    virtual void initialise(MarketData* market);

    // Sensitivity methods
    
    /** Returns the name - used to determine whether to tweak the object */
    virtual string sensName(ShortTermSqueezeTweak* shift)const;
    
    /** Shifts the object using given shift */    
    virtual bool sensShift(ShortTermSqueezeTweak* shift);
    virtual bool sensShift(ShortTermSqueezeTweak* shift, bool useShiftSizeSign);

    /** Returns the name - used to determine whether to tweak the object */
    virtual string sensName(LongTermSqueezeTweak* shift)const;
    
    /** Shifts the object using given shift */    
    virtual bool sensShift(LongTermSqueezeTweak* shift);
    virtual bool sensShift(LongTermSqueezeTweak* shift, bool useShiftSizeSign);

    /** returns the short resp long term squeeze as double */
    static CorrTermDataSP getCorrelationTermSqueezesAndExpiries(
        const DateTime&                 refDate,
        int                             numAssets,
        const CorrelationCommonArray&   corrObjects, 
        const CorrelationTermArray&     corrTermArray);

public:
    string              name;                   // optional
    string              category1;
    string              category2;
    
    double              corrShortTermSqueeze;
    ExpirySP            corrShortTermExpiry;    // optional
    
    double              corrLongTermSqueeze;
    ExpirySP            corrLongTermExpiry;     // optional

    // transient fields
    bool                isZeroObject;           // whether or not its a zero object
    
    CorrelationTerm();
    CorrelationTerm(bool isZeroObject);
    
    static DoubleArray spotCorrelation(
                double                      refCorrelation,
                ExpirySP                    benchmarkTermExpiry,
                double                      corrShortTermSqueeze,
                ExpirySP                    corrShortTermExpiry,
                double                      corrLongTermSqueeze,
                ExpirySP                    corrLongTermExpiry,
                const DateTime&             valueDate,
                const DateTimeArray&        toDates,
                const TimeMetricConstSP&    timeMetricOne, 
                const TimeMetricConstSP&    timeMetricTwo);

    static DoubleMatrixArraySP  CorrelationTermMatrix( 
                const DateTime&             valueDate, 
                const DateTimeArray&        toDates,
                const DoubleMatrix&         fwdVarAtDates,                
                const TimeMetricArray&      timeMetricArray,
                CorrTermDataSP              data,
                _TCalcType                  calcType,
                double                      eigenValueFloor,
                double                      maxSqError,
                const DoubleArray&          barriers,
                const BoolArray &           skipFwdCorrelation);
};

// support for wrapper class
typedef MarketWrapper<CorrelationTerm> CorrelationTermWrapper;
#ifndef QLIB_CORRELATIONTERM_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CorrelationTerm>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CorrelationTerm>);
#endif

/** specialisations of arrayObjectCast */
template <> class MARKET_DLL arrayObjectCast<CorrelationTermWrapper>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const CorrelationTermWrapper& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(CorrelationTermWrapper& value);

    /** Turns the IObjectSP into a DateTime */
    static CorrelationTermWrapper fromIObject(IObjectSP& value);

};
// arrays of wrappers (note array of structures)
typedef array<CorrelationTermWrapper, CorrelationTermWrapper> CorrelationTermWrapperArray;
#ifndef QLIB_CORRELATIONTERM_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<CorrelationTermWrapper _COMMA_ CorrelationTermWrapper>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CorrelationTermWrapper _COMMA_ CorrelationTermWrapper>);
#endif


DRLIB_END_NAMESPACE
#endif // CORRELATIONTERM_HPP
