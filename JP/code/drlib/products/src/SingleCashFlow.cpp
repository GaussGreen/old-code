// MUST always include this FIRST
#include "edginc/config.hpp"
// other include
#include "edginc/Instrument.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Addin.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

class SingleCashFlow: public CInstrument, 
                      public virtual ClosedForm::IIntoProduct,
                      public virtual Theta::Shift {

public:
    static CClassConstSP const TYPE;

    /** copy market data relevant to the instrument */
    virtual void GetMarket(const IModel*, const CMarketDataSP);

    virtual void Validate();

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual ClosedForm::IProduct* createProduct(ClosedForm* model) const;

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** Implementation of theta shift */
    virtual bool sensShift(Theta* shift);

    void price(Control* control, CResults* results) const;

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;
    
private:
    //// Handle any output requests that we support
    void requests(Control* control, CResults* results) const;

    // for reflection
    static void load(CClassSP& clazz);
    static IObject* defaultSingleCashFlow();

    SingleCashFlow();

    double amount;
    DateTime cfDate;
    DateTime valueDate;
    YieldCurveWrapper discount;
};

typedef smartPtr<SingleCashFlow> SingleCashFlowSP;


// some validation, here just for the example
void SingleCashFlow::Validate() {
    // do some validation (called just after market data has been collected)
}

// copy market data relevant to the instrument
void SingleCashFlow::GetMarket(const IModel* model, const CMarketDataSP market) {
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);
}


// main pricing method
void SingleCashFlow::price(Control* control, CResults* results) const {
    static const string method = "SingleCashFlow::price";
    try {
        double value;
        if (cfDate <= valueDate) {
            value = 0.0;
        } else {
            value = amount * discount->pv(cfDate);
        }
        results->storePrice(value, discount->getCcy());
        // then do any output requests
        if (control->isPricing()) {
            // but only on price run ie not during tweaks
            requests(control, results);
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

//// Handle any output requests that we support
void SingleCashFlow::requests(Control* control, CResults* results) const {
    static const string method = "SingleCashFlow::requests";
    try {
        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArray dates(1, cfDate);
            OutputRequestUtil::recordPaymentDates(control,results, &dates); 
        }
        
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            OutputRequestUtil::KnownCashFlows kcf;
            kcf.addKnownCashFlow(discount->getCcy(), cfDate, amount);
            kcf.recordKnownCashFlows(control, results);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** what's today ? */
DateTime SingleCashFlow::getValueDate() const {
    return valueDate;
}

/** private class */
class SingleCashFlowClosedForm: public ClosedForm::IProduct{
private:
    const SingleCashFlow* cf; // a reference

public:
    SingleCashFlowClosedForm(const SingleCashFlow* cf): cf(cf){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{
        cf->price(control, results);
    }
};
    
/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* SingleCashFlow::createProduct(
    ClosedForm* model) const{
    return new SingleCashFlowClosedForm(this);
}

bool SingleCashFlow::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);
    }
    catch (exception& e) {
        throw ModelException(e, "SingleCashFlow::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}


/** Returns the name of the instrument's discount currency */
string SingleCashFlow::discountYieldCurveName() const {
    return discount.getName();
}


/** Invoked when Class is 'loaded' */
void SingleCashFlow::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SingleCashFlow, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(ClosedForm::IIntoProduct);
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultSingleCashFlow);
    FIELD(amount, "cash flow amount");
    FIELD(cfDate, "cash flow date");
    FIELD(discount, "identifies discount curve");
    FIELD(valueDate, "valuation date");
    FIELD_MAKE_TRANSIENT(valueDate);
    Addin::registerConstructor(
        "SINGLE_CASH_FLOW",
        Addin::UTILITIES,
        "SingleCashFlow constructor",
        TYPE);
}

IObject* SingleCashFlow::defaultSingleCashFlow(){
    return new SingleCashFlow();
}

CClassConstSP const SingleCashFlow::TYPE = CClass::registerClassLoadMethod(
    "SingleCashFlow", typeid(SingleCashFlow), SingleCashFlow::load);

SingleCashFlow::SingleCashFlow() : CInstrument(TYPE) {}

/* to ensure class is linked in */
bool SingleCashFlowLoad() {
    return true;
}


DRLIB_END_NAMESPACE
