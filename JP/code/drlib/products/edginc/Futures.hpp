//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Futures.hpp
//
//   Description : Futures price
//
//
//----------------------------------------------------------------------------

#ifndef EDG_FUTURES_HPP
#define EDG_FUTURES_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/VAsset.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CFutures : public CInstrument, 
                              public CClosedFormLN::IIntoProduct,
                              virtual public VIXFModel::IIntoProduct,
                              public LastSensDate, public Theta::IShift,
                              virtual public ISensitiveStrikes
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultFutures(){
        return new CFutures();
    }

    // override base implementation if required
    virtual void GetMarket(const IModel*, const CMarketDataSP);
    
    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** Vol Index Future */
    /** Implementation of VIXFModel::IntoProduct interface */
    virtual VIXFModel::IProduct* createProduct(const VIXFModel* model) const;

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual DateTime getValueDate() const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    /** false for VAsset, True otherwise */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which this instrument is sensitive */
    /** should be called only if asset = VAsset */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

private:
    // for product to access instrument data
    friend class CFuturesClosedFormProd;
    friend class VFutureProd;

    void requests(Control* control, CResults* results) const;
      
    CFutures();

    DateTime                valueDate;
	DateTime                matDate;

    CAssetWrapper           asset;
    string                  ccyTreatment;
    YieldCurveWrapper       discount;

    // just for pricing on maturity date, unregistered
    double spotAtMat; // $unregistered
};

DRLIB_END_NAMESPACE
#endif
