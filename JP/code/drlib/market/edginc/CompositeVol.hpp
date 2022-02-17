//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CompositeVol.hpp
//
//   Description : Organises the building of a composite vol from a correlated
//                 set of assets
//
//   Author      : Mark A Robson
//
//   Date        : 12 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_COMPOSITE_VOL_HPP
#define EDG_COMPOSITE_VOL_HPP

#include "edginc/Asset.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE


/** the composite vol does not exist as a separate object - the vol returned
    is an instance of another type of vol, currently a vol surface with only
    one strike (should be a vol curve). But the type of vol returned could
    depend on the component vols if desired.

    For clients who wish to create a composite vol there is only method,
    create, of interest. 

    Currenly, the CompositeVol only supports volatilites where the
    CProcessedVol is derived from CVolProcessedBS. The implementation works
    by passing an instance of a CompositeVol to each of the processed 
    volatilities. Each processed vol then supplies the relevant information
    to the the CompositeVol.
*/
class MARKET_DLL CompositeVol{
public:
    /** As above, but with correlation skew. Set correlationSkew = 0.0 to
        switch off correlation skew */
    static IVolProcessed* getProcessedVolCorrelationSkewAndTerm(
        const CVolRequest*        topVolRequest,    // top level request
        const CAsset*             topAsset,         // top level asset   
        bool                      fwdStartAsset,    // is asset fwd starting ?
        const DateTime&           assetStartDate,
        const string&             name,             // name for vol
        const TimeMetric*         metric,           // defines business time
        const CAssetWrapperArray& assets,           // array of assets 
        const CVolRequestLNArray& interps,          // array of vol interps 
        const CDoubleArray&       weights,          // weight for each asset
        const CDoubleMatrix&      correlations,     // correlations between assets
        const ExpiryArray&        corrBenchmarkExpiries,
        const double              correlationSkew,  // correlation skew of basket
        const double              correlationSkewPower,// correlation skew power
        const CDoubleMatrix&      corrShortTermSpread, // term structure of correlation
        const ExpiryArray&        corrShortTermExpiries,
        const CDoubleMatrix&      corrLongTermSpread,
        const ExpiryArray&        corrLongTermExpiries); // term structure of correlation

    /** add a list of dates to the CompositeVol. These dates determine
        when the composite vols should be calculated. Vols should
        use addVols when little calculation is needed for the vols */
    void addDates(
        const DateTime&       baseDate,    // base date for compDates
        const ExpiryArray*    compDates);  // dates for calculation of comp vol

    /* add a list of dates and the corresponding vol to the
       CompositeVol.  These dates determine when the composite vols
       should be calculated.  Vols should use addDates when lots of
       calculation is needed for the vols */
    void addVols(
        const DateTime&       baseDate,    // base date for compDates
        const ExpiryArray*    compDates,   // dates for calculation of comp vol
        const CDoubleArray&   vols);       // vols at compDates

    /** Constructs composite vol object which can be used to create multiple
        processed vols. Derived asset information is cached between calls */
    CompositeVol(const CAssetWrapperArray& assets);

    /** As above but takes CompositeVol object created by constructor */
    IVolProcessed* getProcessedVol(
        const CVolRequest*        topVolRequest,    // top level request
        const CAsset*             topAsset,         // top level asset   
        bool                      fwdStartAsset,    // is asset fwd starting ?
        const DateTime&           assetStartDate,
        const string&             name,             // name for vol
        const TimeMetric*         metric,           // defines business time
        const CVolRequestLNArray& interps,          // array of vol interps 
        const CDoubleArray&       weights,          // weight for each asset
        const CDoubleMatrix&      correlations,     // correlations between assets
        const ExpiryArray&        corrBenchmarkExpiries,
        const double              correlationSkew,  // correlation skew of basket
        const double              correlationSkewPower,// correlation skew power
        const CDoubleMatrix&      corrShortTermSpread, // term structure of correlation
        const ExpiryArray&        corrShortTermExpiries,
        const CDoubleMatrix&      corrLongTermSpread,
        const ExpiryArray&        corrLongTermExpiries); // term structure of correlation
    
    /** Returns array of expiries for last generated vol (will return
        an empty array if no vols yet generated) */
    const ExpiryArray& getBMExpiries();

    /** Return expiries as BenchmarkDates - handy for mix of fixed and rolling
        benchmarks */
    ExpiryArraySP fixedExpiries() const;

    /** Like getProcessedVol but returns the array of vols that would have 
        been fed into VolSurface constructor */
    DoubleArray  getVolsAtBMDates(
        const CVolRequest*        topVolRequest,    // top level request
        const CAsset*             topAsset,         // top level asset   
        bool                      fwdStartAsset,    // is asset fwd starting ?
        const DateTime&           assetStartDate,
        const string&             name,             // name for vol
        const TimeMetric*         metric,           // defines business time
        const CVolRequestLNArray& interps,          // array of vol interps 
        const CDoubleArray&       weights,          // weight for each asset
        const CDoubleMatrix&      correlations,     // correlations between assets
        const ExpiryArray&        corrBenchmarkExpiries,
        const double              correlationSkew,  // correlation skew of basket
        const double              correlationSkewPower,// correlation skew power
        const CDoubleMatrix&      corrShortTermSpread, // term structure of correlation
        const ExpiryArray&        corrShortTermExpiries,
        const CDoubleMatrix&      corrLongTermSpread,
        const ExpiryArray&        corrLongTermExpiries); // term structure of correlation

private:
    /** Gets component vols to send data to CompositeVol */
    void populate(
        const CVolRequestLNArray& interps);     /* (I) array of vol interps */
    
    /** function for reporting rubbish input data */
    static void CompVolTooLarge(const CAssetWrapper&  asset1, /* (I) */
                                const CAssetWrapper&  asset2, /* (I) */
                                double                vol1,   /* (I) */
                                double                vol2,   /* (I) */
                                double                corr,   /* (I) */
                                const DateTime&       date);  /* (I) */

    /** function for correlation skew */
    static void ProcessCorrelation( 
        const int                   iDate,                  //
        const double                yearFrac,               //
        const vector<CDoubleArray>& fwds,                   //
        const CDoubleArray&         weights,                //
        const CDoubleMatrix&        correlations,           // atm correlation matrix
        const double                correlationSkew,        //
        const double                correlationSkewPower,   //
        const double                refLevel,               // strike in percentage
        double*                     corrSqueeze,            /* (O) */
        double*                     limitSqueeze,           /* (O) */
        bool*                       useCorrelationSkew );   /* (O) */

    // static fields
    /** represents maximum size allowed of a necessary partial result
       in the calculation of the composite vol */
    static const double VOL_PRODUCT_TOO_LARGE;
    
    // fields
    DateTime              baseDate;
    DateTimeArray         compDates;     /* derived from compExpiries */
    ExpiryArray           compExpiries;  /* used as benchmark dates for 
                                            new vol */
    vector<CDoubleArray>  vols;          /* array of set of vols for each
                                            asset */
    int                   numAssets;     /* total number of assets */
    int                   currentAsset;  /* current asset */
    CAssetWrapperArray    assets;        // array of assets 
    vector<CDoubleArray>  fwds;          // array of set of fwd prices
};

DRLIB_END_NAMESPACE

#endif

