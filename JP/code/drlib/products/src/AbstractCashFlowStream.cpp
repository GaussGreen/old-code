//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : AbstractCashFlowStream.cpp
//
//   Description : Price arbitrary cashflows
//
//   Author      : Gordon Stephens
//
//   Date        : 15 June 2005
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AbstractCashFlowStream.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/CashFlow.hpp"


DRLIB_BEGIN_NAMESPACE

/** private class */
class AbstractCashFlowStreamClosedForm: public ClosedForm::IProduct
{
private:
    const AbstractCashFlowStream* cf; // a reference

public:
    AbstractCashFlowStreamClosedForm(const AbstractCashFlowStream* cf): cf(cf){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{

        //supply just a closed form forward rate pricer
        ClosedFormForwardRatePricerSP fwdRateMdl =
            ClosedFormForwardRatePricerSP(
                new ClosedFormForwardRatePricer());

        cf->price(control, results, fwdRateMdl);
    }
};
    
void AbstractCashFlowStream::Validate()
{
    static const string method = "AbstractCashFlowStream::Validate";
    try
    {
        //nothing to do
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

// copy market data relevant to the instrument
void AbstractCashFlowStream::GetMarket(const IModel* model, const CMarketDataSP market)
{
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);
    //get for the cashflows
    for (int i = 0; i < cfl.getLength(); i++)
    {
        cfl[i]->getMarket(model,market.get());
    }
}


void AbstractCashFlowStream::price(Control* control, CResults* results, IForwardRatePricerSP model) const {
    static const string method = "AbstractCashFlowStream::price";
    try {
        double value = 0.0;
        
        for (int i = 0; i < cfl.getLength(); i++) {
            DateTime cfDate = cfl[i]->getPayDate();
            if (cfDate.isGreater(valueDate)) {
                value += cfl[i]->getAmount(model) * discount->pv(cfDate);
            }
        }
        
        results->storePrice(value, discount->getCcy());

        if (control && control->isPricing() ) {
            requests(control, results, model);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void AbstractCashFlowStream::requests(Control* control, CResults* results, IForwardRatePricerSP model) const {
    static const string method = "AbstractCashFlowStream::requests";
    try {
        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(new DateTimeArray(cfl.getLength()));
            for (int i = 0; i < cfl.getLength(); i++) {
                (*dates)[i] = cfl[i]->getPayDate();
            }
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }
        
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) 
        {
            CashFlowArraySP cfarray = CashFlowArraySP(new CashFlowArray());
            for(int i = 0; i < cfl.getLength(); i++)
            {
                CashFlow tmp(cfl[i]->getPayDate(), cfl[i]->getAmount(model));
                cfarray->push_back(tmp);
            }
            //throw ModelException(method, "KNOWN_CASHFLOWS not implemented");
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfarray.get()); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

AbstractCashFlowStream::AbstractCashFlowStream()
    : CInstrument(TYPE)
{
}

AbstractCashFlowStream::AbstractCashFlowStream(AbstractCashFlowArray cfl,
                                               const string&  discountName)
    : CInstrument(TYPE),
      cfl(cfl),
      discount(discountName)
{
}

/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* AbstractCashFlowStream::createProduct(ClosedForm* model) const
{
    return new AbstractCashFlowStreamClosedForm(this);
}

/** what's today ? */
DateTime AbstractCashFlowStream::getValueDate() const
{
    return valueDate;
}

/** when to stop tweaking */
DateTime AbstractCashFlowStream::endDate(const Sensitivity*) const
{
    DateTime maxDate; //initialise to 0
    for (int i = 0; i < cfl.getLength(); i++) {
        DateTime cfDate = cfl[i]->getPayDate();
        if (cfDate.isGreater(maxDate)) {
            maxDate = cfDate;
        }
    }
    return maxDate;
}

bool AbstractCashFlowStream::sensShift(Theta* shift)
{
    try
    {
        valueDate = shift->rollDate(valueDate);
    }
    catch (exception& e)
    {
        throw ModelException(e, "SimpleCashFlowStream::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}


/** Returns the name of the instrument's discount currency */
string AbstractCashFlowStream::discountYieldCurveName() const
{
    return discount.getName();
}


IObject* AbstractCashFlowStream::defaultAbstractCashFlowStream()
{
    return new AbstractCashFlowStream();
}

void AbstractCashFlowStream::load(CClassSP& clazz)
{
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AbstractCashFlowStream, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultAbstractCashFlowStream);
        FIELD(cfl, "cash flows");
        FIELD(discount, "identifies discount curve");
        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
    }

CClassConstSP const AbstractCashFlowStream::TYPE = CClass::registerClassLoadMethod(
    "AbstractCashFlowStream", typeid(AbstractCashFlowStream), load);


bool  AbstractCashFlowStreamLoad() {
    return (AbstractCashFlowStream::TYPE != 0);
       }

DRLIB_END_NAMESPACE
