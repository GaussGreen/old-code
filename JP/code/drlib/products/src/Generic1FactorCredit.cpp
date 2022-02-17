//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Generic1FactorCredit.cpp
//
//   Description : Base class for 1 factor generic credit instruments
//
//   Author      : Stephen Hope
//
//   Date        : 27 Sep 2002
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Generic1FactorCredit.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/ClosedFormCDSPS.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/ClosedFormBSImpliedSmile.hpp"
#include "edginc/ClosedFormMultiQSmile.hpp"

DRLIB_BEGIN_NAMESPACE

/* for reflection */
Generic1FactorCredit::Generic1FactorCredit(CClassConstSP clazz): 
    CInstrument(clazz), oneContract(false), notional(0.0)
{
    // empty
}

/** Do some asset specific validation */
void Generic1FactorCredit::validate() const{
    // currently empty - code previously here now done in getMarket call
}


/** Get the asset , par curve and discount market data */
void Generic1FactorCredit::GetMarket(const IModel*        model, 
                                     const CMarketDataSP  market)
{
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);

    if (!cdsParSpreads.isEmpty()) {
        ICDSParSpreads::getMarketData(model, 
                                      market.get(),
                                      discount.getName(),
                                      cdsParSpreads);
    }
    
    if (!!instSettle)
    {
        instSettle->getMarket(model, market.get());
    }
    if (!!premiumSettle) 
    {
        premiumSettle->getMarket(model, market.get());
    }

     // FirmAsset stuff ?????????????? ///
    // Dont need to get FirmAsset mkt data if we are only pricing via par curves
    if (!ClosedFormCDSPS::TYPE->isInstance(model) && !CClosedFormLN::TYPE->isInstance(model)
        && !ClosedFormBSImpliedSmile::TYPE->isInstance(model)
        && !ClosedFormMultiQSmile::TYPE->isInstance(model))
    {
        CAsset::getAssetMarketData(model, market.get(), "N", 
                                   discount, asset);
    }


    // call instrument specific getMarket routine if applicable
    if (IGetMarket::TYPE->isInstance(this)) {
        IGetMarket* imnt = dynamic_cast<IGetMarket*>(this);
        imnt->getMarket(model,market.get());
    }
}

DateTime Generic1FactorCredit::getValueDate()const
{
    return valueDate;
}

/** Get the sensitive strike for this volRequest.
    Will be called from the product that inherits 
    from this class not from the infrastructure **/
void Generic1FactorCredit::getSensStrikes(
    const OutputNameConstSP&         outputName,
    const CVolRequest*               volRequest,
    const SensitiveStrikeDescriptor& sensStrikeDesc,
    const DoubleArraySP&             sensitiveStrikes)
{
    asset->getSensitiveStrikes(volRequest, outputName, 
                               sensStrikeDesc, sensitiveStrikes);
}

//// roll through time 
bool Generic1FactorCredit::sensShift(Theta* theta)
{
    try 
    {
        valueDate = theta->rollDate(valueDate);
    }
    catch (exception& e) {
        throw ModelException(e, "Generic1FactorCredit::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool Generic1FactorCredit::avoidVegaMatrix(const IModel* model){
    // this does raise the question of whether we should have method(s)
    // on the model rather than the instrument
    if (MonteCarlo::TYPE->isInstance(model)){
        const MonteCarlo& mc = dynamic_cast<const MonteCarlo&>(*model);
        return !mc.vegaMatrixSupported(this);
    }
    return true;
}

/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP Generic1FactorCredit::getSensitiveStrikes(OutputNameConstSP outputName,
                                                  const IModel*      model){
    // this does raise the question of whether we should have method(s)
    // on the model rather than the instrument
    if (MonteCarlo::TYPE->isInstance(model)){
        const MonteCarlo& mc = dynamic_cast<const MonteCarlo&>(*model);
        return mc.getSensitiveStrikes(this, outputName);
    }
    throw ModelException("Generic1FactorCredit::getSensitiveStrikes",
                         "Sensitive Strikes not supported");
}


/** Returns the name of the instrument's discount currency. */
string Generic1FactorCredit::discountYieldCurveName() const {
    return discount.getName();
}


/* Returns the cdsParSpreads */
const ICDSParSpreads* Generic1FactorCredit::getICDSParSpreads() const {
    return const_cast<const ICDSParSpreads*>(cdsParSpreads.get());
}


class Generic1FactorCreditHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Generic1FactorCredit, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(ISensitiveStrikes);
        FIELD(valueDate,"valuation Date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(oneContract,"Is it one contract");
        // making oneContract optional is not ideal, however we want the
        // CreditDefaultSwap not to have it at all - this is the best 
        // compromise until/if we rewrite CreditDefaultSwap
        FIELD_MAKE_OPTIONAL(oneContract); // default false
        FIELD(notional,"Option notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(ccyTreatment,"Currency Treatment");
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD_MAKE_OPTIONAL(instSettle);
        FIELD(premiumSettle, "Premium settlement");
        FIELD_MAKE_OPTIONAL(premiumSettle);
        FIELD(asset,"Credit Underlying of option");
        FIELD_MAKE_OPTIONAL(asset);
        FIELD(cdsParSpreads, "Underlying of option");
        FIELD(discount,"Discount curve");
    }

};

CClassConstSP const Generic1FactorCredit::TYPE = 
CClass::registerClassLoadMethod(
    "Generic1FactorCredit", typeid(Generic1FactorCredit), 
    Generic1FactorCreditHelper::load);

DRLIB_END_NAMESPACE

