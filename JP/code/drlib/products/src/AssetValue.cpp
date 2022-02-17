//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetValue.cpp
//
//   Description : Asset instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 7 September 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AssetValue.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/VegaPointwise.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/AssetUtil.hpp"

DRLIB_BEGIN_NAMESPACE

void AssetValue::GetMarket(const IModel*       model, 
                           const CMarketDataSP market)
{
    static const string method("AssetValue::GetMarket");
    try {
        market->GetReferenceDate(valueDate);
        asset.getData(model, market);
        Validate();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void AssetValue::Validate() {
    static const string method("AssetValue::Validate");
    // asset could come from the market - ie. will be NULL after pop2obj. 
    // Do not validate type if it is NULL.
    if (asset.get()){
        // check that asset is a SimpleEquity
        if (!IAssetFairValue::TYPE->isInstance(asset.get())) {
            string m(asset->getClass()->getName() + " asset " + 
                     asset->getName() + " is not priceable");
            throw ModelException(method, m);
        }
    }
}

void AssetValue::validatePop2Object() {
    Validate();
}


bool AssetValue::sensShift(Theta* shift) {    
    // roll today 
    valueDate = shift->rollDate(valueDate);    
    return true;
};

DateTime AssetValue::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime AssetValue::endDate(const Sensitivity* sensControl) const 
{
    DateTime end = asset->settleDate(valueDate);
    return end;
}

/** private class */
class AssetValueClosedForm: public CClosedFormLN::IProduct{
private:
    const AssetValue* asset; // a reference

public:
    AssetValueClosedForm(const AssetValue* asset): asset(asset) {}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results) const;
};

void AssetValueClosedForm::price(CClosedFormLN* model,
                                 Control*       control,
                                 CResults*      results) const
{
    static const string method = "AssetValueClosedForm::price";
    try {
        const IObject* o = asset->asset.get();

        const IAssetFairValue* fv = dynamic_cast<const IAssetFairValue*>(o);

        if (!fv) {
            throw ModelException(method,
                                 asset->asset->getClass()->getName() + 
                                 " asset " + asset->asset->getName() + 
                                 " is not priceable");
        }


        double price = fv->fairValue();

        // store price in output request
        results->storePrice(price, AssetUtil::assetCcy(asset->asset.get()));

        // check whether VEGA sensitivities have been requested
        if (control && control->isPricing())  {
            // get array of all sensitivities
            SensitivityArrayConstSP allSens = control->getSens();

            // get array of all sensitivities
            for (int i = 0 ; i < allSens->size() ; ++i) {
                if (VegaParallel::TYPE->isInstance((*allSens)[i].get())      ||
                    VegaPointwise::TYPE->isInstance((*allSens)[i].get())     ||
                    VegaSkewParallel::TYPE->isInstance((*allSens)[i].get())  ||
                    VegaSkewPointwise::TYPE->isInstance((*allSens)[i].get()) ||
                    VegaMatrix::TYPE->isInstance((*allSens)[i].get())        ||
                    RootTimeVega::TYPE->isInstance((*allSens)[i].get())) {
                    results->storeNotApplicable((*allSens)[i].get());
                }
            }
        }

        // no output requests (forward at maturity, indicative vol or
        // delay price)
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* AssetValue::createProduct(
    CClosedFormLN* model) const {

    return new AssetValueClosedForm(this);
}


/** Returns the name of the instrument's discount currency. */
string AssetValue::discountYieldCurveName() const {
    return asset.getName();
}


// for reflection
AssetValue::AssetValue(): CInstrument(TYPE) {}

class AssetValueHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AssetValue, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::IShift);
        EMPTY_SHELL_METHOD(defaultAsset);
        FIELD(asset, "asset");
        FIELD(valueDate, "valuation Date");
        FIELD_MAKE_OPTIONAL(valueDate);
    }

    static IObject* defaultAsset(){
        return new AssetValue();
    }
};

CClassConstSP const AssetValue::TYPE = CClass::registerClassLoadMethod(
    "AssetValue", typeid(AssetValue), AssetValueHelper::load);
bool  AssetValueLoad() {
    return (AssetValue::TYPE != 0);
   }

   
DRLIB_END_NAMESPACE
