//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Futures.cpp
//
//   Description   value a stock given a quoted price
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Futures.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// helpers
void CFutures::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CFutures, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    IMPLEMENTS(VIXFModel::IIntoProduct);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(Theta::IShift);
    EMPTY_SHELL_METHOD(defaultFutures);
    FIELD(valueDate, "valuation Date");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(matDate, "Maturity date");
    FIELD(asset, "Underlying of option");
    FIELD(discount, "Discount curve");
    FIELD(ccyTreatment, "Currency Treatment");
}

CClassConstSP const CFutures::TYPE = CClass::registerClassLoadMethod(
    "FuturesContract", typeid(CFutures), load);
bool  CFuturesLoad() {
    return (CFutures::TYPE != 0);
}



void CFutures::Validate()
{
    static const string method = "Futures::Validate";
    try {
        // asset and discount curve could come from the market - ie. will
        // be NULL after pop2obj. Do not cross validate if either of them
        // is NULL.
        if (asset.get() && discount.get()) {
            AssetUtil::assetCrossValidate(asset.get(),
                                          false,
                                          valueDate,
                                          valueDate,
                                          discount,
                                          this);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// initiate GetMarket 
void CFutures::GetMarket(const IModel*       model, 
                         const CMarketDataSP market)
{
    static const string method = "CFutures::GetMarket";
    try {
        market->GetReferenceDate(valueDate);
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                                   discount, asset);

        
        const VIXFModel* vsModel = dynamic_cast<const VIXFModel*>(model);
        if (vsModel) { //VIX Futures
            discount.getData(vsModel->varSwapModel.get(), market);            
        }else{
            discount.getData(model, market);
        }

        // since the assets have changed, we have to validate the instrument
        validatePop2Object();

        // asset and discount curve must not be null 
        if (!asset) {
            throw ModelException(method, "Asset is NULL.");
        }

        if (!discount) {
            throw ModelException(method,  "Discount curve is NULL.");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// constructor
CFutures::CFutures(): CInstrument(TYPE), spotAtMat(0.0) {}

/** Rolls the value date and sets initial spot if rolling over start date */
bool CFutures::sensShift(Theta* shift)
{    
    DateTime aDate = shift->rollDate(valueDate);

    if ((aDate >= matDate && valueDate < matDate ) ||
        (valueDate == matDate && Maths::isZero(spotAtMat))) {
        spotAtMat = asset->getThetaSpotOnDate(shift, matDate);
    }

    // roll today 
    valueDate = aDate;
    
    return true;
};

DateTime CFutures::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime CFutures::endDate(const Sensitivity* sensControl) const {

    DateTime assetEnd = asset->settleDate(matDate); // this is for spot asset

    if (const LastSensDate* ptr = dynamic_cast<const LastSensDate*>(asset.get())){
        DateTime end = ptr->endDate(sensControl); // this is more general
        int daysDiff = end.daysDiff(valueDate);
        end = matDate.rollDate(daysDiff);
        assetEnd = assetEnd > end? assetEnd : end;
    }

    return assetEnd;
}

void CFutures::requests(Control* control, CResults* results) const {
    static const string method = "CFuturesClosedFormProd::requests";
    try {
        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(new DateTimeArray(1, matDate));
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }
        
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            if (valueDate.isGreaterOrEqual(matDate)) {
                CashFlowArraySP payoff(new CashFlowArray(1));
                (*payoff)[0].date   = matDate;
                (*payoff)[0].amount = spotAtMat;

                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        payoff.get()); 
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
/** false for VAsset, True otherwise */
bool CFutures::avoidVegaMatrix(const IModel* model){
    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());

    if (vAsset) {
        return false;
    }else{
        return true;
    }
}

/** returns all strikes on the vol surface to which this instrument is sensitive */
/** should be called only if asset = VAsset */
DoubleArraySP CFutures::getSensitiveStrikes(OutputNameConstSP outputName,
                                            const IModel*      model){

    const VAsset* vAsset = dynamic_cast<const VAsset*>(asset.get());

    if (vAsset) {
        const VIXFModel* m = dynamic_cast<const VIXFModel*>(model);
        if (!m) {
            throw ModelException("Model must be of type VIXFModel for Volatility Index Future.");
        }
        return vAsset->sensitiveStrikes(outputName, m);
    }else{ //shouldn't be called
        DoubleArraySP temp(new DoubleArray(asset->getSpot()));
        return temp;
    }
}

/** private class */
class CFuturesClosedFormProd: public CClosedFormLN::IProduct{
private:
    const CFutures*  instrFutures; // a reference

public:
    CFuturesClosedFormProd(const CFutures* instr): instrFutures(instr){}

    void price(CClosedFormLN*   model,
               Control*        control, 
               CResults*       results) const;

    //void requests(Control* control, CResults* results) const;
};

void CFuturesClosedFormProd::price(CClosedFormLN*  model,
                                   Control*        control, 
                                   CResults*       results) const
{
    static const string method = "CFuturesClosedFormProd::price";
    try {
        double fwdAtMat;  // forward at maturity

        if (instrFutures->valueDate == instrFutures->matDate) {
            fwdAtMat = instrFutures->spotAtMat;
        }
        else if (instrFutures->valueDate < instrFutures->matDate) {
            fwdAtMat = instrFutures->asset->fwdValue(instrFutures->matDate);
        }
        else {
            fwdAtMat = 0.0;
        }
 
        results->storePrice(fwdAtMat,instrFutures->discount->getCcy());

        if (control && control->isPricing()) {
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           instrFutures->matDate,
                                           instrFutures->valueDate,
                                           instrFutures->asset.get());

            instrFutures->requests(control, results);
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}
/*  
void CFuturesClosedFormProd::requests(Control* control, CResults* results) const {
    static const string method = "CFuturesClosedFormProd::requests";
    try {
        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(new DateTimeArray(1, instrFutures->matDate));
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }
        
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            if (instrFutures->valueDate.isGreaterOrEqual(instrFutures->matDate)) {
                CashFlowArraySP payoff(new CashFlowArray(1));
                (*payoff)[0].date   = instrFutures->matDate;
                (*payoff)[0].amount = instrFutures->spotAtMat;

                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        instrFutures->discount->getCcy(),
                                                        payoff.get()); 
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
*/
////////////////////////////Vol Index Future////////////////////////////
/** two models in order to call Var Swap BS and SV
for VIX Future
*/ 
class VFutureProd: virtual public VIXFModel::IProduct {
public:
    
    VFutureProd(const CFutures* inst): inst(inst){}
    
    /** This is the method responsible for calling vAsset->fwdValue
        and output additional components (BS, SV) 
    */
    void price(VIXFModel*       model,
               Control*      control, 
               CResults*     results) const;

private:
    const CFutures* inst; // a reference

};

/** This is the method responsible for calling vAsset->fwdValue
    and output additional components (BS, SV) 
*/
void VFutureProd::price(VIXFModel*       model,               
						  Control*          control, 
                          CResults*         results) const {
    static const string method = "VFutureProd::price";
    try {
        const VAsset* vAsset = dynamic_cast<const VAsset*>((inst->asset.get()));

        CControlSP ctrl(copy(control));
        vAsset->getAlgorithm()->setControl(ctrl);

        DateTime valueDate = inst->getValueDate();
        DateTime maturity =  inst->matDate;

        double totalValue = 100.0* vAsset->fwdValue(maturity);

        results->storePrice(totalValue, inst->discount->getCcy());

        /** vector "out" contains BS, SV components, only has values after calling vAsset->fwdValue*/
        /** record "out" into results based on the control */
        vAsset->recordExtraOutput(control, results);

        if (control && control->isPricing()) {
            inst->requests(control, results);

            // FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                results,
                maturity,
                valueDate,
                vAsset);
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
};
////////////////////////////////////////////////////////////////////////////////////////////

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* CFutures::createProduct(
    CClosedFormLN* model) const{
    return new CFuturesClosedFormProd(this);
}

/** Implementation of VIXFModel::IntoProduct interface */
VIXFModel::IProduct* CFutures::createProduct(const VIXFModel* model) const {
    return new VFutureProd(this);
}

/** Returns the name of the instrument's discount currency. */
string CFutures::discountYieldCurveName() const {
    return discount.getName();
}


DRLIB_END_NAMESPACE
