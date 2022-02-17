//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedBSParam.hpp
//
//   Description : Processed BS parameterised vols. 
//                 Previously, was part of VolParam.hpp
//
//   Author      : Regis Guichard
//
//   Date        : 02 Mai 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_PROCESSED_BS_HPP
#define EDR_VOL_PROCESSED_BS_HPP

#include "edginc/VolParam.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/VolRequestLN.hpp"

DRLIB_BEGIN_NAMESPACE



/** Base class for Processed Parametrized Vol Surface.
    Default implementations are straight processing to a ProcessedSurface */
class MARKET_DLL CVolProcessedBSParam: public CVolProcessedBS
{
public:
    static CClassConstSP const TYPE;

    /** Calculates variance between 2 dates */
    virtual double CalcVar(const DateTime &date1,
                           const DateTime &date2) const;

    /** Calculates variance between a series of dates. If the dateList
        has n dates in it, n-1 variances will be calculated. */
    virtual void CalcVar(const DateTimeArray& dateList,
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;

    /** Calculates variance beginning at dateFrom. If the dateList
        has n dates in it, n variances will be calculated. */
    virtual void CalcVar(const DateTime&      dateFrom,
                         const DateTimeArray& datesTo,
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;

    /** Calculates volatility between 2 dates. If the dateList
        has n dates in it, n-1 vols will be calculated. */
    virtual double CalcVol(const DateTime& date1, 
                           const DateTime& date2) const;

    /** Calculates vols between a series of dates */
    virtual void CalcVol(const DateTimeArray& dateList, 
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;

    /** Calculates vols between a start date and a series of dates */
    virtual void CalcVol(const DateTime&      dateFrom,
                         const DateTimeArray& datesTo,
                         TCalcType            calcType, 
                         CDoubleArray&        vols) const;

    /** Populate a CompositeVol with data it requires (essentially
        benchmark dates and vols) */
    virtual void populateCompositeVol(CompositeVol* compositeVol) const;

    /** identifies the market data name of the volatility */
    virtual string getName() const;

    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const;

    /** retieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric() const;

    //// construct a vol processed BS for parameterised vol
    CVolProcessedBSParam(const CVolBase*         vol,
                         const CVolParamConstSP& volParam,
                         const CVolRequestLN*    volRequest, 
                         const CAsset*           asset);

    //// construct a vol processed BS for parameterised vol when struck
    CVolProcessedBSParam(const CVolBase*         vol,
                         const CVolParamConstSP& volParam,
                         const CVolRequestLN*    volRequest, 
                         const CAsset*           asset,
                         const FXAsset*          fxAsset,
                         const Correlation*      eqFXCorr);

    /** updates this VolProcessedBS with the given vol request and asset.
        Returns the updated VolProcessedBS - which may or may not be the
        same as 'this'. Use smart pointers to manage memory */
    CVolProcessedBSParam* update(const CVolRequestLN*    volRequest, 
                                 const CAsset*           asset);
protected:
    CVolProcessedBSParam(const CClassConstSP&    clazz, 
                         const CVolBase*         vol,
                         const CVolParamConstSP& volParam,
                         const CVolRequestLN*    volRequest, 
                         const CAsset*           asset,
                         const FXAsset*          fxAsset = 0,
                         const Correlation*      eqFXCorr = 0);

protected:

    void processToSurface() const;

    CVolBaseConstSP         myVol; // $unregistered
    CVolParamConstSP        myParamVol;
    CVolRequestLNConstSP    myVolRequestLN;
    AssetConstSP            myAsset;
    // the two below are potentially NULL
    FXAssetConstSP          myFXAsset;
    CorrelationConstSP      myEqFXCorr;
    mutable
    CVolProcessedBSConstSP  myProcessedSurface;
    mutable
    VolSurfaceConstSP       myVolSurf;
    mutable DoubleArraySP   strikes; // working area for collecting strikes $unregistered
private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<CVolProcessedBSParam> CVolProcessedBSParamSP;


DRLIB_END_NAMESPACE
#endif 
