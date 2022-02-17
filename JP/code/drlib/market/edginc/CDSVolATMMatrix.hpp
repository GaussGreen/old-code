/**
 * @file CDSVolATMMatrix.hpp
 */

#ifndef QR_CDSVOLATMMATRIX_HPP
#define QR_CDSVOLATMMATRIX_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/ICDSVol.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Model.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/CRVolParallel.hpp"
#include "edginc/CRVolPointwise.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(CDSVolATMMatrix)
FORWARD_DECLARE(MarketData)
FORWARD_DECLARE(TimeMetric)
FORWARD_DECLARE(IVolProcessed)
FORWARD_DECLARE(CVolRequest)
FORWARD_DECLARE(ICDSParSpreads)
FORWARD_DECLARE(ICDS)
FORWARD_DECLARE(Expiry)
FORWARD_DECLARE(ExpiryPair)

/** Concrete class defining ATM volatility matrix for European CDS options. */
class MARKET_DLL CDSVolATMMatrix : public MarketObject,
                        virtual public ITweakableWithRespectTo<CRVolParallel>,
                        virtual public ITweakableWithRespectTo<CRVolPointwise>,
                        public virtual ICDSVol,
						public virtual Theta::Shift {
public:
    static CClassConstSP const TYPE;

    /**Volatilities represent the forward volatility of (1-F), where F is the
       forward upfront-fee payable by the protection buyer.*/
    static const string CDS_VOL_TYPE_CLEAN_PRICE;
    /**Volatilities represent the forward volatility of the CDS par spread.*/
    static const string CDS_VOL_TYPE_SPREAD;
    static const string CDS_VOL_INTERPOLATION_LINEAR;
    static const string CDS_VOL_INTERPOLATION_SPLINE;

    /** Returns name of vol */
    virtual string getName() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    virtual void validatePop2Object();

    /*=========================================================================
     * ICDSVol interface methods
     * - for now, just CDSVolRequestSimpleEuropean, which is suitable for
     *   European options with closed form forward distributions (e.g. BS 
     *   implied vol or multi-Q type models).
     *=======================================================================*/
    virtual IVolProcessed* getProcessedVol(const CVolRequest*    volRequest,
                                           const ICDSParSpreads* cds) const;

    /** Set the tweak name for CR vega tweaks. The tweak name will be
        identical with the name of the cube vol object that uses this vol
        object*/
    void setCRVegaTweakName(const string & tweakName) const;

    /**Find the metric time to exercise, the volatility, 
       and the indices into the expiry and ul maturity vectors. These are positive integers
       if there's an exact match. If not, the values are negative, and the entry should
       be inserted at -(r+1) to maintain the ordering.*/
    virtual void getInterpolatedVolatility(
        const DateTime& exDate, 
        const DateTime& ulMaturity,
        double* timeToExercise,
        double* volatility,
        int* expIndex,
        int* ulIndex1,
        int* ulIndex2) const;

    virtual ~CDSVolATMMatrix();

    const ExpiryArraySP getOptionExpiries() const;
    const ExpiryArraySP getUnderlyingExpiries() const;
    const TimeMetricSP getTimeMetric() const;
    const DoubleMatrix& getATMVolMatrix() const;
    const DateTime& getBaseDate() const;
    const bool areULExpiriesRelativeToBase() const;

    // Parallel credit vega tweak (sensitivity) support
    virtual string sensName(const CRVolParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CRVolParallel>&);
    
    // Pointwise credit vega tweak (sensitivity) support
    virtual string sensName(const CRVolPointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<CRVolPointwise>&);
    virtual ExpiryPairArrayConstSP sensQualifiers(const CRVolPointwise*) const;

    // Get Expiry window for credit vega tweak
    ExpiryArrayConstSP getExpiries() const;

	/*=========================================================================
     * Tweak interface methods
     *=======================================================================*/
    bool sensShift(Theta* theta);

protected:
    CDSVolATMMatrix(const CClassConstSP& clazz);
    CDSVolATMMatrix();
    
    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Name of the vol matrix*/
    string             name;    
    /**Defines how volatility time should be measured - trading days? simple? etc.
       Default should be simple time. */
    TimeMetricSP       metric;  
    /**Trade date/today*/
    DateTime           baseDate; 
    /**Is the vol a spread vol or a price vol? Must have value 
       "CDS_VOL_TYPE_SPREAD" or "CDS_VOL_TYPE_PRICE"*/
    string             volType;
    /**Sample CDS that is typical of the CDS whose vols are expressed by this matrix.*/
    ICDSSP             ulGenerator;
    /** IMPORTANT: atmInputVol[i][j] corresponds to the Black-Scholes implied volatility for the jth 
        option expiry, ith ul maturity. use negative numbers for missing values. */      
    DoubleMatrix       atmVol; 
    /**Option expiries in strictly increasing order. If these are relative expiries (e.g. "3M") as
       opposed to absolute dates, they will be interpreted as being relative to the baseDate.*/
    ExpiryArraySP      optionExpiries; 
    /**Underlying expiries in strictly increasing order. How these are interpreted is determined by
       ulExpiriesAreRelativeToOptionExpiries.*/
    ExpiryArraySP      ulExpiries;
    /**If false, the ulExpiries are relative to the optionExpiries; 
       if true, they are relative to the baseDate. Clearly, if the ulExpiries are absolute dates
       this makes no difference either way.*/
    bool               ulExpiriesRelativeToBase;
    /**Interpolation type for vols for underlying/option-exercise tenors which are not 
       specified exactly. Must have value "CDS_VOL_INTERPOLATION_LINEAR" 
       or "CDS_VOL_INTERPOLATION_SPLINE"*/
    string             interpolationType;

    //the name for CR vega tweak
    mutable string crTweakName;

    friend class CDSVolATMMatrixHelper;
};

DRLIB_END_NAMESPACE
#endif

