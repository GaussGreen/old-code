//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AUD3YBondFuture.cpp
//
//   Description : AUD3YBondFuture instrument
//
//   Author      : Keiji Kitazawa
//
//   Date        : 27 Jun 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AUD3YBondFuture.hpp" // AJS - can do without header altogether
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

const double TICK_SHIFT = 0.0001; // 0.01% move

/* for reflection */
AUD3YBondFuture::AUD3YBondFuture(): GenericSimpleIR(TYPE)
{
    //empty
}

void AUD3YBondFuture::Validate()
{
    static const string method = "AUD3YBondFuture::Validate";

    try
    {
        if (valueDate > startDate)
        {
            throw ModelException(method, "valueDate (" + valueDate.toString() +
                                         ") is after startDate of bond (" + startDate.toString() +
                                         ")");
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

// copy market data relevant to the instrument
void AUD3YBondFuture::GetMarket(const IModel* model, const CMarketDataSP market){
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);

    PseudoBond->getMarket(model, market.get());
}

void AUD3YBondFuture::price(Control* control,
                        CResults* results) const
{
    static const string method = "AUD3YBondFuture::price";
    
    try
    {        
        double value = PseudoBond->presentValue(startDate, discount.getSP());
        value /= PseudoBond->getFaceValue();
        results->storePrice(value, discount->getCcy());
        
        if (control && control->isPricing())
                requests(control,results,value);

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** when to stop tweaking */
DateTime AUD3YBondFuture::endDate(const Sensitivity* sensControl) const 
{
    return PseudoBond->getMaturityDate();
}

void AUD3YBondFuture::requests(Control*        control, 
                               CResults*       results,
                               double          value) const {
    static const string method = "AUD3YBondFuture::requests";
    try {
        OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
        
        CashFlowArrayConstSP cfl = PseudoBond->getCashFlows(valueDate);

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

        // Store Tick Value
        request = control->requestsOutput(OutputRequest::TICK_VALUE);
        if (request) 
        {// Get implied yield from CTD Bond Price
            value *= PseudoBond->getFaceValue();    // Put back to contract value of Bond
            double x1 = 0.0;         //Initial value = 0.0%
            double x2 = 100.0;      //Initial value = 100%
            double f1 = PseudoBond->priceFromYield(x1, false, startDate) - value;
            double f2 = PseudoBond->priceFromYield(x2, false, startDate) - value;
            double dx = (x1-x2)/10.0; //factor to move
            for (int i=0;i<=200;i++)
            {
                if (f1*f2 >= 0.0) 
                {
                    if (fabs(f2) < 0.0001)
                        break;
                    else
                    {
                        x2 -= dx;       // push back to original step
                        dx /= 5.0;      // reduce the move size
                    }
                }
                x2 += dx; 
                f2 = PseudoBond->priceFromYield(x2, false, startDate) - value;
            }
            double impyield = x2;
            // Twkeak implied yield
            double tickYield = impyield + TICK_SHIFT;
            // price using our shifted yield
            double tickPrice = PseudoBond->priceFromYield(tickYield, false, startDate);
            double tick = tickPrice - value;
            results->storeRequestResult(request, tick);
        }    
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}



/** private class */
class AUD3YBondFutureClosedForm: public ClosedForm::IProduct{
private:
    const AUD3YBondFuture* lib; // a reference

public:
    AUD3YBondFutureClosedForm(const AUD3YBondFuture* lib): lib(lib){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{
        lib->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
ClosedForm::IProduct* AUD3YBondFuture::createProduct(
    ClosedForm* model) const{
    return new AUD3YBondFutureClosedForm(this);
}


class AUD3YBondFutureHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AUD3YBondFuture, clazz);
        SUPERCLASS(GenericSimpleIR);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultAUD3YBondFuture);
        FIELD(startDate, "Forward Starting Date of Bond");        
        FIELD(PseudoBond, "Underlying Bond, 3Y 6% semi-annual Bond");
    }

    static IObject* defaultAUD3YBondFuture(){
        return new AUD3YBondFuture();
    }
};

CClassConstSP const AUD3YBondFuture::TYPE = CClass::registerClassLoadMethod(
    "AUD3YBondFuture", typeid(AUD3YBondFuture), AUD3YBondFutureHelper::load);
bool  AUD3YBondFutureLoad() {
    return (AUD3YBondFuture::TYPE != 0);
}

   

DRLIB_END_NAMESPACE
