//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SyntheticAsset.cpp
//
//   Description : synthetic asset - prices just like its underlying
//
//   Author      : Andrew J Swain
//
//   Date        : 14 March 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/AssetUtil.hpp"

DRLIB_BEGIN_NAMESPACE

class SyntheticAsset: public Generic1Factor,
                      public CClosedFormLN::IIntoProduct,
                      public LastSensDate {
public:
    static CClassConstSP const TYPE; 

    virtual void validatePop2Object(){
        static const string method = "SyntheticAsset::validatePop2Object";
        // no such thing as fwd starting
        if (fwdStarting) {
            throw ModelException(method, "forward starting not allowed");
        }
    }

    virtual void Validate() {
        static const string method("SyntheticAsset::Validate");
        // check that asset is 'priceable'
        if (!IAssetFairValue::TYPE->isInstance(asset.get())) {
            throw ModelException(method, 
                                 asset->getClass()->getName() + " asset " + 
                                 asset->getName() + " is not priceable");
        }
        AssetUtil::assetCrossValidate(asset.get(),
                                      false,
                                      valueDate,
                                      valueDate,
                                      discount,
                                      this);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("synthetic asset");
        REGISTER(SyntheticAsset, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultSyntheticAsset);
        FIELD(maturity, "maturity");
        FIELD(finalLevel, "spot at maturity");
        FIELD_MAKE_OPTIONAL(finalLevel);
    }

    static IObject* defaultSyntheticAsset(){
        return new SyntheticAsset();
    }
    
private:
    friend class SyntheticAssetClosedForm;

    SyntheticAsset():Generic1Factor(TYPE) {}; 
    SyntheticAsset(const SyntheticAsset& rhs);
    SyntheticAsset& operator=(const SyntheticAsset& rhs);

    void requests(Control* control, CResults* results) const {
        static const string method = "CashFlowStream::requests";
        try {
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray dates(1, maturity);
                OutputRequestUtil::recordPaymentDates(control,results,&dates); 
            }
        
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request && valueDate.isGreaterOrEqual(maturity)) {
                CashFlowArray payments;
                payments.push_back(CashFlow(maturity, finalLevel));
                
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        &payments); 
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    void price(Control* control, CResults* results)const{
        static const string method = "SyntheticAsset::price";        
        try  {
            if (priceDeadInstrument(control, results)) {
                return; // all done;
            }
            const IObject* o = asset.get();
            const IAssetFairValue* fv = dynamic_cast<const IAssetFairValue*>(o);

            if (!fv) {
                throw ModelException(method,
                                     asset->getClass()->getName() + 
                                     " asset " + asset->getName() + 
                                     " is not priceable");
            }

            results->storePrice(fv->fairValue(), discount->getCcy());    
            requests(control, results);
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    
    DateTime endDate(const Sensitivity* sensControl) const {
        return asset->settleDate(maturity);
    }

    virtual bool priceDeadInstrument(CControl* control, CResults* results) const {       
        if (maturity.isGreater(valueDate)) {
            return false;
        }
        results->storePrice(0.0, discount->getCcy());
        requests(control, results);
        return true;   
    }   

    /** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

private:
    DateTime maturity;
    double   finalLevel;
};

CClassConstSP const SyntheticAsset::TYPE = CClass::registerClassLoadMethod(
    "SyntheticAsset", typeid(SyntheticAsset), SyntheticAsset::load);


/** private class */
class SyntheticAssetClosedForm: public CClosedFormLN::IProduct{
private:
    const SyntheticAsset* sa; // a reference

public:
    SyntheticAssetClosedForm(const SyntheticAsset* sa): sa(sa){}

    void price(CClosedFormLN* model,
               Control*    control, 
               CResults*   results)const{
        sa->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* SyntheticAsset::createProduct(CClosedFormLN* model) const
{
    return new SyntheticAssetClosedForm(this);
}

// for class loading 
bool SyntheticAssetLoad() {
    return (SyntheticAsset::TYPE != 0);
}

DRLIB_END_NAMESPACE
