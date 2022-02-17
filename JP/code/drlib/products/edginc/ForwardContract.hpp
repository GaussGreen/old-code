//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ForwardContract.hpp
//
//   Description   forward contract
//
//
//----------------------------------------------------------------------------

#ifndef EDG_FORWARD_CONTRACT_HPP
#define EDG_FORWARD_CONTRACT_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/NumericalIntegrationLN.hpp"
#include "edginc/ForwardContractCreditSupport.hpp"
#include "edginc/Theta.hpp"
#include "edginc/VAsset.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CForwardContract : public CInstrument, 
                         public CClosedFormLN::IIntoProduct,
                         public NumericalIntegrationLN::IIntoProduct,
                         virtual public VIXFModel::IIntoProduct,
                         public Theta::IShift,
                         public LastSensDate,
                         public CreditSupport::Interface,
                         virtual public ISensitiveStrikes
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultForwardContract(){
        return new CForwardContract();
    }

    // override base implementation if required
    virtual void GetMarket(const IModel*, const CMarketDataSP);
    
    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** Implementation of NumericalIntegrationLN::IntoProduct interface */
    virtual NumericalIntegrationLN::IProduct* createProduct(NumericalIntegrationLN* model) const;

    /** Vol Index Forward */
    /** Implementation of VIXFModel::IntoProduct interface */
    virtual VIXFModel::IProduct* createProduct(const VIXFModel* model) const;

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    virtual DateTime getValueDate() const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    virtual CreditSupportSP createCreditSupport(CMarketDataSP market);

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
    void recordRequests(Control* control, Results* results) const;

    // for product to access instrument data
    friend class CForwardContractClosedFormProd;
    // for credit exposure support
    friend class ForwardContractCreditSupport;

    friend class ForwardContractNumerical;

    friend class VForwardProd;

    CForwardContract();

    DateTime                valueDate;
    DateTime                startDate;
    InstrumentSettlementSP  premiumSettle;
    InstrumentSettlementSP  instSettle;
    bool                    fwdStarting;
    bool                    oneContract;

    double                  notional;
    double                  initialSpot;
    double                  spotAtMaturity;
    ScheduleSP              exerciseSchedule;

    // old dividend re-invest (zeroing dividend) flag is back
    bool                    divReinvest;

    CAssetWrapper           asset;
    string                  ccyTreatment;
    YieldCurveWrapper       discount;
};

typedef smartConstPtr<CForwardContract> ForwardContractConstSP;
typedef smartPtr<CForwardContract> ForwardContractSP;

DRLIB_END_NAMESPACE
#endif
