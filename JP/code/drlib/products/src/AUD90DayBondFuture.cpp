//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AUD90DayBondFuture.cpp
//
//   Description : Interest rate future for AUD 90 days bond future.
//
//   Author      : Keiji Kitazaw
//
//   Date        : 20 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AUD90DayBondFuture.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

const double TICK_SHIFT = 0.0001; // 0.01% move

void AUD90DayBondFuture::Validate() {
    static const string method = "AUD90DayBondFuture::Validate";
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

void AUD90DayBondFuture::price(Control* control, CResults* results) const{
    static const string method = "AUD90DayBondFuture::price";
    try {
        double value = 0.0;
        double rate = 0.0;

        DayCountConventionSP dayCount(DayCountConventionFactory::make(dcc));
        
        if (valueDate <= lastTradeDate) {
            double irvol = 0.0;

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

            double futRate = discount->future(intStartDate,
                                              intEndDate,
                                              dayCount.get(),
                                              CompoundBasis::SIMPLE,
                                              irvol);

            rate = futRate;
            value = 1.0/(1.0+90.0/365.0*rate);   
        }    

        results->storePrice(value, discount->getCcy());

        
        if (control->isPricing()) 
        {// Store Tick Value
            OutputRequest* request =
                control->requestsOutput(OutputRequest::TICK_VALUE);
            if (request && valueDate <= lastTradeDate) 
            {
                double tickPrice = 1.0/(1.0+(rate+TICK_SHIFT)*90.0/365.0);
                double tick = tickPrice - value;
                results->storeRequestResult(request, tick);
            }
            // Add PAYMENT_DATES
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                //DateTime paymentDate = instSettle->settles(intEndDate, asset.get());
                DateTimeArray date(1, intEndDate);
                OutputRequestUtil::recordPaymentDates(control,results,&date); 
            }

            // Add KNOW_CASHFLOWS
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                // if there is known equity payoff, need to handle it.
                //DateTime paymentDate = instSettle->settles(intEndDate, asset.get());
                DateTimeArray date(1, intEndDate);
                // Pick up only known Fixed Leg.

                //CashFlowArraySP cfl (CashFlow(intEndDate, value));
                CashFlow cf (intEndDate, value);
                CashFlowArray cfl(1, cf);
                
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        &cfl);   
            }

        }

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** private class */
class AUD90DayBondFutureClosedForm: public virtual CClosedFormLN::IProduct, 
                                    public virtual ClosedFormIRLN::IProduct{
private:
    const AUD90DayBondFuture* irf; // a reference

public:
    AUD90DayBondFutureClosedForm(const AUD90DayBondFuture* irf): irf(irf){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results) const{
        irf->price(control, results);
    }

    void price(ClosedFormIRLN* model,
               Control*        control, 
               CResults*       results) const{
        irf->price(control, results);
    }
};
    
/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* AUD90DayBondFuture::createProduct(
    CClosedFormLN* model) const{
    return new AUD90DayBondFutureClosedForm(this);
}
ClosedFormIRLN::IProduct* AUD90DayBondFuture::createProduct(
    ClosedFormIRLN* model) const{
    return new AUD90DayBondFutureClosedForm(this);
}

/** when to stop tweaking */
DateTime AUD90DayBondFuture::endDate(const Sensitivity* sensControl) const {
    return intEndDate;
}

IRGridPointAbsArraySP AUD90DayBondFuture::getSensitiveIRVolPoints(
    OutputNameConstSP outputName,
    const IModel*      model) const {
    static const string method = "AUD90DayBondFuture::getSensitiveIRVolPoints";
    try {
        MaturityPeriodSP tenor(MaturityPeriod::dateSubtract(intEndDate,
                                                            intStartDate));
        
        SwapMaturityVolRequest volRequest(tenor.get());
        DateTimeArray dates(1, intStartDate);
        return VolProcessedBSIR::sensitiveIRVolPoints(vol.get(), 
                                                      discount.get(),
                                                      &volRequest, dates);
#if 0
        IRGridPointArraySP points(new IRGridPointArray(0));

        MaturityPeriodSP tenor(MaturityPeriod::dateSubtract(intEndDate,
                                                            intStartDate));

        CVolRequestSP    volRequest(new SwapMaturityVolRequest(tenor.get()));
        DateTimeArray    dates(1, intStartDate);
        
        vol->sensitiveIRVolPoints(volRequest.get(),
                                  outputName,
                                  dates,
                                  points);

        return points;
#endif
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// for reflection
AUD90DayBondFuture::AUD90DayBondFuture(): GenericSimpleIR(TYPE) {}

class AUD90DayBondFutureHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AUD90DayBondFuture, clazz);
        SUPERCLASS(GenericSimpleIR);;
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(ClosedFormIRLN::IIntoProduct);
        IMPLEMENTS(ISensitiveIRVolPoints);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultAUD90DayBondFuture);
        FIELD(intStartDate, "startDate");
        FIELD(intEndDate, "endDate");
        FIELD(lastTradeDate, "lastTradeDate");
        FIELD(dcc, "day count convention");
    }

    static IObject* defaultAUD90DayBondFuture(){
        return new AUD90DayBondFuture();
    }
};

CClassConstSP const AUD90DayBondFuture::TYPE = CClass::registerClassLoadMethod(
    "AUD90DayBondFuture", typeid(AUD90DayBondFuture), AUD90DayBondFutureHelper::load);

bool  AUD90DayBondFutureLoad() {
    return (AUD90DayBondFuture::TYPE != 0);
}

   

DRLIB_END_NAMESPACE
