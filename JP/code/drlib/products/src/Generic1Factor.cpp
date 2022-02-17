//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Generic1Factor.cpp
//
//   Description : Base class for 1 factor generic instruments
//
//   Author      : Stephen Hope
//
//   Date        : 4 Sep 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/AssetUtil.hpp"

#include "edginc/IndexSpecEQ.hpp"

DRLIB_BEGIN_NAMESPACE

Generic1Factor::~Generic1Factor(){}

/* for reflection */
Generic1Factor::Generic1Factor(CClassConstSP clazz): 
    CInstrument(clazz), notional(0.0), initialSpot(0.0)
{
    // empty
}

/** Do some asset specific validation */
void Generic1Factor::validate()
{
    static const string method = "Generic1Factor::validate";
    try {
        AssetUtil::assetCrossValidate(asset.get(),
                                      fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Get the asset and discount market data */
void Generic1Factor::GetMarket(const IModel*          model, 
                               const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);
    CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                               discount, asset);

    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
    if( premiumSettle.get() )
        premiumSettle->getMarket(model, market.get());
}

DateTime Generic1Factor::getValueDate()const
{
    return valueDate;
}

/** Adds the DELAY_PRICE and FWD_AT_MAT if requested */
void Generic1Factor::addRequests(Control* control,
                                 Results* results,
                                 double fairValue,
                                 const DateTime& maturityDate)const
{
    // add the DELAY_PRICE if we have a premiumSettle
    InstrumentUtil::delayPriceHelper(control,
                                     results,
                                     fairValue,
                                     valueDate,
                                     discount.get(),
                                     asset.get(),
                                     premiumSettle.get());

    // FWD_AT_MAT
    InstrumentUtil::recordFwdAtMat(control,
                                   results,
                                   maturityDate,
                                   valueDate,
                                   asset.get());
}

/** Get the sensitive strike for this volRequest.
    Will be called from the product that inherits 
    from this class not from the infrastructure **/
void Generic1Factor::getSensStrikes(
    const OutputNameConstSP&         outputName,
    const CVolRequest*               volRequest,
    const SensitiveStrikeDescriptor& sensStrikeDesc,
    const DoubleArraySP&             sensitiveStrikes)
{
    asset->getSensitiveStrikes(volRequest, outputName, 
                               sensStrikeDesc, sensitiveStrikes);
}

//// roll through time 
bool Generic1Factor::sensShift(Theta* theta){
    // get the new date
    const DateTime& newDate = theta->rollDate(valueDate);
            
    /* If fwd start date falls between value date and theta date 
       have to set initial spot.*/
    if ( fwdStarting                           && 
         startDate.isGreaterOrEqual(valueDate) && 
         newDate.isGreaterOrEqual(startDate)   ) 
    {
        // returns the fwd price if shift is a theta fs
        initialSpot = asset->getThetaSpotOnDate(theta, startDate);
        
        // not fwd starting anymore 
        fwdStarting = false;
    }  
    
    // roll today 
    valueDate = newDate;
    
    return true; // continue to tweak components which implement Theta
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool Generic1Factor::avoidVegaMatrix(const IModel* model){
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
DoubleArraySP Generic1Factor::getSensitiveStrikes(OutputNameConstSP outputName,
                                                  const IModel*      model){
    // this does raise the question of whether we should have method(s)
    // on the model rather than the instrument
    if (MonteCarlo::TYPE->isInstance(model)){
        const MonteCarlo& mc = dynamic_cast<const MonteCarlo&>(*model);
        return mc.getSensitiveStrikes(this, outputName);
    }
    throw ModelException("Generic1Factor::getSensitiveStrikes",
                         "Sensitive Strikes not supported");
}


/** Returns the name of the instrument's discount currency. */
string Generic1Factor::discountYieldCurveName() const {
    return discount.getName();
}


class Generic1FactorHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Generic1Factor, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(ISensitiveStrikes);
        FIELD(valueDate,"valuation Date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(startDate,"Option start date");
        FIELD_MAKE_OPTIONAL(startDate);
        FIELD(fwdStarting,"Is it a fwd starting option");
        FIELD(oneContract,"Is it one contract");
        FIELD(notional,"Option notional");
        FIELD_MAKE_OPTIONAL(notional);
        FIELD(initialSpot,"Initial spot price");
        FIELD_MAKE_OPTIONAL(initialSpot);
        FIELD(ccyTreatment,"Currency Treatment");
        FIELD(instSettle, "Instrument settlement at maturity");
        FIELD(premiumSettle, "Premiumsettlement");
        FIELD_MAKE_OPTIONAL(premiumSettle);
        FIELD(asset,"Underlying of option");
        FIELD(discount,"Discount curve");
    }

};

CClassConstSP const Generic1Factor::TYPE = CClass::registerClassLoadMethod(
    "Generic1Factor", typeid(Generic1Factor), Generic1FactorHelper::load);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//          Credit support for Generic1Factor                           //     
//                                                                      //
//////////////////////////////////////////////////////////////////////////


CAssetSP Generic1FactorCreditSupport::getAsset() const
{
    return getInst()->asset.getSP();
}

string Generic1FactorCreditSupport::getInstCcyCode() const
{
    return getInst()->discount->getCcy();
}

Generic1FactorCreditSupport::Generic1FactorCreditSupport()
{}

DRLIB_END_NAMESPACE

