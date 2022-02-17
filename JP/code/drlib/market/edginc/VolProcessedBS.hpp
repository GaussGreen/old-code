//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CVolProcessedBS.hpp
//
//   Description : What a VolatilityBS can do
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOL_PROCESSED_BS_HPP
#define VOL_PROCESSED_BS_HPP
#include "edginc/VolProcessed.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE
class CompositeVol;
class CVolBase;
class CVolRequest;
class CAsset;

FORWARD_DECLARE(IPDFBoundaryProb);

/** Defines what a Black-Scholes processed volatility can do. Can be thought
    of as a Vol Curve */
class MARKET_DLL CVolProcessedBS: public CObject,
                       public virtual IVolProcessed{
public:
    static CClassConstSP const TYPE;
    static const double LEGACY_PDF_BOUNDARY_PROB;

    /** Used to indicate the type of calculation should be performed when
        an array of dates is passed to CalcVar/CalcVol */
    typedef enum _TCalcType{
        fromFirst = 0, // calculate vol/variance from first date
        forward,       // calculate vol/variance between successive dates
        toLast         // calculate vol/variance to last date
    } TCalcType;

    /** Calculates variance between 2 dates */
    virtual double CalcVar(const DateTime &date1,
                           const DateTime &date2) const = 0;
    
    /** Calculates variance between a series of dates. If the dateList
        has n dates in it, n-1 variances will be calculated. */
    virtual void CalcVar(const DateTimeArray& dateList,
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;

    /** Calculates variance beginning at dateFrom. If the dateList
        has n dates in it, n variances will be calculated. */
    virtual void CalcVar(const DateTime &dateFrom,
                         const DateTimeArray& datesTo,
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;


    /** Calculates volatility between 2 dates. If the dateList
        has n dates in it, n-1 vols will be calculated. */
    virtual double CalcVol(const DateTime& date1, 
                           const DateTime& date2) const = 0;

    /** Calculates vols between a series of dates */
    virtual void CalcVol(const DateTimeArray& dateList, 
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;

    /** Calculates vols between a start date and a series of dates */
    virtual void CalcVol(const DateTime &dateFrom,
                         const DateTimeArray& datesTo,
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;

    /** Populate a CompositeVol with data it requires (essentially
        benchmark dates and vols) */
    virtual void populateCompositeVol(CompositeVol* compositeVol) const = 0;

    /** utility method - interpolates vol then calculates vols from 
        dateFrom to each of datesTo */
    static void calcVol(const CVolBase*       vol,
                        const CVolRequest*    request,
                        const CAsset*         asset,
                        const DateTime&       dateFrom,
                        const DateTimeArray&  datesTo,
                        CDoubleArray&         vols);

    /** Generates an array of 'default' strikes given an asset. The
        strikes are chosen such that they are evenly distributed in
        probability space based upon the at-the-money 1 year vol. Strikes
        are also generated for probabilities "0+epsilon" and "1-epsilon". The
        maxNumBands is an indicative indicator for the number of strikes to
        generate. pdfBoundaryProb determines the range of strikes to return
        by limiting the probability distribution from above. */
    static DoubleArray defaultStrikes(int             maxNumBands,
                                      const DateTime& baseDate,
                                      const CAsset*   asset,
                                      const IPDFBoundaryProb* marketObj);

    /** As above but uses maxNumBands = 20 ie every 5% */
    static DoubleArray defaultStrikes(const DateTime& baseDate,
                                      const CAsset*   asset,
                                      const IPDFBoundaryProb* marketObj);

protected:
    CVolProcessedBS(const CClassConstSP& clazz);
private:
    CVolProcessedBS(const CVolProcessedBS &rhs);
    CVolProcessedBS& operator=(const CVolProcessedBS& rhs);
    
};

typedef smartConstPtr<CVolProcessedBS> CVolProcessedBSConstSP;
typedef smartPtr<CVolProcessedBS> CVolProcessedBSSP;
#ifndef QLIB_VOLPROCESSEDBS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CVolProcessedBS>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CVolProcessedBS>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CVolProcessedBS>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CVolProcessedBS>);
#endif

DRLIB_END_NAMESPACE
#endif
