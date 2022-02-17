//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRFuture.cpp
//
//   Description : Interest rate future
//
//   Author      : Andrew J Swain
//
//   Date        : 11 January 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IRFuture.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/Results.hpp"
#include "edginc/DayCountConventionFactory.hpp"

DRLIB_BEGIN_NAMESPACE

// copy market data relevant to the instrument
void IRFuture::GetMarket(const IModel* model, const CMarketDataSP market) {
    market->GetReferenceDate(valueDate);
    // vol should come from YieldCurve
    discount.getData(model, market);
}


void IRFuture::Validate() {
    static const string method = "IRFuture::Validate";
    try {
        if (intStartDate.isGreater(intEndDate)) {
            throw ModelException(method,
                                 "start of borrowing period (" + 
                                 intStartDate.toString() + 
                                 ") is not before end of borrowing period (" +
                                 intEndDate.toString() + ")");
        }

        if (lastTradeDate.getDate() > intStartDate.getDate()) {
            throw ModelException(method,
                                 "last trading date (" + 
                                 lastTradeDate.toString() + 
                                 ") is after the start of the borrowing "
                                 "period (" + intStartDate.toString() + ")");
        }

        // this code only supports interest rate futures where the time between
        // the interest start and interest end dates is not greater than one 
        // year. i.e the forward rate is a simple interest rate rather than 
        // a swap rate 
        DateTime oneYear = MaturityPeriod::toDate(1, "Y", intStartDate);
        
        if (intEndDate.isGreater(oneYear)) {
            throw ModelException(method,
                                 "Algorithm cannot price futures on underlying"
                                 " rates of more than one year duration");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void IRFuture::price(CResults* results) const{
    static const string method = "IRFuture::price";
    try {
        double value = 0.0;
        DayCountConventionSP dayCount(DayCountConventionFactory::make(dcc));
        
        if (valueDate.isLess(lastTradeDate)) {
            double irvol = 0.0;

#if 0
            if (!vol) {
                // leave it as zero
            }
            else {
                MaturityPeriodSP tenor(MaturityPeriod::dateSubtract(intEndDate,
                                                                    intStartDate));

                CVolRequestSP volRequest(new SwapMaturityVolRequest(tenor.get()));

                CVolProcessedSP volCurve(vol->getProcessedVol(volRequest.get(), 0));

                CVolProcessedBS& volBS = dynamic_cast<CVolProcessedBS&>(*volCurve);

                irvol = volBS.CalcVol(valueDate, intStartDate);
            }
#endif
            double futRate = discount->future(intStartDate,
                                              intEndDate,
                                              dayCount.get(),
                                              CompoundBasis::SIMPLE,
                                              irvol);

            value = 1.0 - futRate;
        }       
        results->storePrice(value, discount->getCcy());
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** private class */
class IRFutureClosedForm: public CClosedFormLN::IProduct{
private:
    const IRFuture* irf; // a reference

public:
    IRFutureClosedForm(const IRFuture* irf): irf(irf){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results) const{
        irf->price(results);
    }
};
    
/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* IRFuture::createProduct(
    CClosedFormLN* model) const{
    return new IRFutureClosedForm(this);
}

/** what's today ? */
DateTime IRFuture::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime IRFuture::endDate(const Sensitivity* sensControl) const {
    return intEndDate;
}

bool IRFuture::sensShift(Theta* shift) {
    try {
        valueDate = shift->rollDate(valueDate);
    }
    catch (exception& e) {
        throw ModelException(e, "IRFuture::sensShift (theta)");
    }    
    return true; // our components have theta type sensitivity
}


/** Returns the name of the instrument's discount currency. */
string IRFuture::discountYieldCurveName() const {
    return discount.getName();
}


// for reflection
IRFuture::IRFuture(): CInstrument(TYPE) {}

class IRFutureHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IRFuture, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultIRFuture);
        FIELD(intStartDate, "startDate");
        FIELD(intEndDate, "endDate");
        FIELD(lastTradeDate, "lastTradeDate");
        FIELD(dcc, "day count convention");
        FIELD(vol, "vol");
        FIELD_MAKE_OPTIONAL(vol);
        FIELD(discount, "identifies discount curve");
        FIELD(valueDate, "valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
    }

    static IObject* defaultIRFuture(){
        return new IRFuture();
    }
};

CClassConstSP const IRFuture::TYPE = CClass::registerClassLoadMethod(
    "IRFuture", typeid(IRFuture), IRFutureHelper::load);
bool  IRFutureLoad() {
    return (IRFuture::TYPE != 0);
}


DRLIB_END_NAMESPACE
