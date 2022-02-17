//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : SimpleCashFlowStream.cpp
//
//   Description : CashFlow instrument with no credit or tax support
//                 Taken from CashFlowStream.hpp
//
//   Author      : Gordon Stephens
//
//   Date        : 13 April 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SimpleCashFlowStream.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

void SimpleCashFlowStream::Validate() {
    static const string method = "SimpleCashFlowStream::Validate";
    try {
        for (int i = 1; i < cfl->getLength(); i++) {
            const DateTime& current  = (*cfl)[i].date;
            const DateTime& previous = (*cfl)[i-1].date;
            if (previous.isGreater(current)) {
                string m("cash flow date (" + current.toString() + ") " +
                         "must be greater than the previous cashflow date " + 
                         previous.toString());
                throw ModelException(method, m);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// copy market data relevant to the instrument
void SimpleCashFlowStream::GetMarket(const IModel* model, const CMarketDataSP market) {
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);
}


void SimpleCashFlowStream::price(Control* control, CResults* results) const {
    static const string method = "SimpleCashFlowStream::price";
    try {
        double value = 0.0;
        
        for (int i = 0; i < cfl->getLength(); i++) {
            DateTime cfDate = (*cfl)[i].date;
            if (cfDate.isGreater(valueDate)) {
                value += (*cfl)[i].amount * discount->pv(cfDate);
            }
        }
        
        results->storePrice(value, discount->getCcy());

        if (control && control->isPricing() ) {
            requests(control, results);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void SimpleCashFlowStream::requests(Control* control, CResults* results) const {
    static const string method = "SimpleCashFlowStream::requests";
    try {
        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(new DateTimeArray(cfl->getLength()));
            for (int i = 0; i < cfl->getLength(); i++) {
                (*dates)[i] = (*cfl)[i].date;
            }
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }
        
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfl.get()); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

SimpleCashFlowStream::SimpleCashFlowStream(CashFlowArray* cfl,
                                           const string&  discountName) : CInstrument(TYPE), discount(discountName), cfl(copy(cfl))
{
}

SimpleCashFlowStream::SimpleCashFlowStream(CClassConstSP clazz,
                                           CashFlowArray* cfl,
                                           const string&  discountName) : CInstrument(clazz), discount(discountName), cfl(copy(cfl))
{
}

/** private class */
class SimpleCashFlowStreamClosedForm: public ClosedForm::IProduct{
private:
    const SimpleCashFlowStream* cf; // a reference

public:
    SimpleCashFlowStreamClosedForm(const SimpleCashFlowStream* cf): cf(cf){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{
        cf->price(control, results);
    }
};
    
/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* SimpleCashFlowStream::createProduct(
    ClosedForm* model) const{
    return new SimpleCashFlowStreamClosedForm(this);
}

/** what's today ? */
DateTime SimpleCashFlowStream::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime SimpleCashFlowStream::endDate(const Sensitivity* sensControl) const {
    int last = cfl->getLength() -1;
    return ((*cfl)[last].date);
}

bool SimpleCashFlowStream::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);
    }
    catch (exception& e) {
        throw ModelException(e, "SimpleCashFlowStream::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}


/** Returns the name of the instrument's discount currency */
string SimpleCashFlowStream::discountYieldCurveName() const {
    return discount.getName();
}


// for reflection
SimpleCashFlowStream::SimpleCashFlowStream(): CInstrument(TYPE) {}
SimpleCashFlowStream::SimpleCashFlowStream(CClassConstSP clazz): CInstrument(clazz) {}

class SimpleCashFlowStreamHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SimpleCashFlowStream, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultSimpleCashFlowStream);
        FIELD(cfl, "cash flows");
        FIELD(discount, "identifies discount curve");
        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
    }

    static IObject* defaultSimpleCashFlowStream(){
        return new SimpleCashFlowStream();
    }
};

CClassConstSP const SimpleCashFlowStream::TYPE = CClass::registerClassLoadMethod(
    "SimpleCashFlowStream", typeid(SimpleCashFlowStream), SimpleCashFlowStreamHelper::load);
   

DRLIB_END_NAMESPACE
