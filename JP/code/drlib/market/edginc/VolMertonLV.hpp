//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolMertonLV.hpp
//
//   Description : VolMertonLV
//
//   Author      : Francois Lu
//
//   Date        : 2 Feb 2004
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VolMertonLV_HPP
#define EDR_VolMertonLV_HPP

#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/VolSpline.hpp"

DRLIB_BEGIN_NAMESPACE

typedef smartPtr<VolSpline> VolSplineSP;

class MARKET_DLL VolMertonLV: public CVolBase
{
private:

	friend class MertonLVCalib;
	friend class MertonLVOptimizer;
	friend class MertonLVCalibClosedFormProd;
	friend class VolMertonLVProcessed;

	void createMertonVolSurface(VolSurfaceSP& volSurface);

    static void acceptValueDateCollector(VolMertonLV* vol,
										CValueDateCollector* collector);

    // registered fields
    string name;    // name of the vol
    double ATMVol;
    double JumpRate;
    double JumpMean;
    double JumpWidth;
    DoubleMatrix    impliedMertonVols; // optional input

    // transient fields (won't appear in dd interface)
    DateTime      baseDate;
    TimeMetricSP  timeMetric;
	VolSurfaceSP  mertonSurface;
    VolSplineSP   mertonSurfaceSpline;

    VolMertonLV();

public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    static IObject* defaultCtor();
    static CClassConstSP const TYPE;

    // void validatePop2Object();

    /** Returns name of vol */
    virtual string getName() const;

    /** populate from market cache - default implementation provided */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    // version for quanto and struck
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

};
typedef smartPtr<VolMertonLV> VolMertonLVSP;

DRLIB_END_NAMESPACE

#endif
