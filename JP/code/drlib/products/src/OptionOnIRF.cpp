//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OptionOnIRF.cpp
//
//   Description : Option on interest rate future
//
//   Author      : Andrew J Swain
//
//   Date        : 2 December 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericSimpleIR.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Black.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/SensitiveIRVolPoints.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/ClosedFormIRLN.hpp"

DRLIB_BEGIN_NAMESPACE

class OptionOnIRF: public GenericSimpleIR,
                   public virtual CClosedFormLN::IIntoProduct,
                   public virtual ClosedFormIRLN::IIntoProduct,
                   public virtual ISensitiveIRVolPoints,
                   public virtual LastSensDate {
public:
    static CClassConstSP const TYPE; 

    virtual void validatePop2Object(){
        static const string method = "OptionOnIRF::validatePop2Object";
        // turn strike from price to rate dimension
        k  = (100.0 - strike)/100.0;

        // likewise fixing
        fix = (100.0 - fixing)/100.0;
    }

    virtual void Validate() {
        static const string method = "OptionOnIRF::Validate";
        try {
            if (!vol.get()) {
                throw ModelException(method, "no IR vol supplied");
            }

            if (futureExpiry.isLess(optionExpiry)) {
                throw ModelException(method,
                                     "future expiry (" + 
                                     futureExpiry.toString() + 
                                     ") is before option expiry (" + 
                                     optionExpiry.toString() + ")");
            }
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Option on interest rate future");
        REGISTER(OptionOnIRF, clazz);
        SUPERCLASS(GenericSimpleIR);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(ClosedFormIRLN::IIntoProduct);
        IMPLEMENTS(ISensitiveIRVolPoints);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultOptionOnIRF);
        FIELD(isCall, "isCall");
        FIELD(strike, "strike (e.g. 97 not 0.03)");
        FIELD(period, "period");
        FIELD(dcc, "dcc");
        FIELD(fixing, "fixing (e.g. 97 not 0.03)");
        FIELD_MAKE_OPTIONAL(fixing);
        FIELD(optionExpiry, "optionExpiry");
        FIELD(futureExpiry, "futureExpiry");
        FIELD(k, "k");
        FIELD(fix, "fix");
        FIELD_MAKE_TRANSIENT(k);
        FIELD_MAKE_TRANSIENT(fix);
    }

    static IObject* defaultOptionOnIRF(){
        return new OptionOnIRF();
    }
    
private:
    friend class OptionOnIRFClosedForm;

    OptionOnIRF():GenericSimpleIR(TYPE) {}; 
    OptionOnIRF(const OptionOnIRF& rhs);
    OptionOnIRF& operator=(const OptionOnIRF& rhs);

    void requests(Control* control, CResults* results) const {
        static const string method = "CashFlowStream::requests";
        try {
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray dates(1, optionExpiry);
                OutputRequestUtil::recordPaymentDates(control,results,&dates); 
            }
        
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                // if  historic fixing can be used
                if (valueDate.isGreaterOrEqual(optionExpiry)) {
                    double premium;
                    // inverted put/call logic
                    if (!isCall) {
                        premium = Maths::max(fix - k, 0.0);
                    }
                    else {
                        premium = Maths::max(k - fix, 0.0);
                    }
                    CashFlowArray payments;
                    payments.push_back(CashFlow(optionExpiry, premium));
                
                    OutputRequestUtil::recordKnownCashflows(control,
                                                            results,
                                                            discount->getCcy(),
                                                            &payments); 
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    void price(Control* control, CResults* results)const{
        static const string method = "OptionOnIRF::price";
        
        try  {
            if (priceDeadInstrument(control, results)) {
                return; // all done;
            }

            DayCountConventionSP daycount(DayCountConventionFactory::make(dcc));
            MaturityPeriod       futurePeriod(period);
            DateTime             futureEnd = futurePeriod.toDate(futureExpiry);

            CVolRequestSP volRequest(new SwapMaturityVolRequest(&futurePeriod));           
            CVolProcessedSP volCurve(vol->getProcessedVol(volRequest.get(), 0));            
            CVolProcessedBS& volBS = dynamic_cast<CVolProcessedBS&>(*volCurve);
            
            double variance = volBS.CalcVar(valueDate, optionExpiry);
            double vol      = volBS.CalcVol(valueDate, optionExpiry);
        
            double futRate = discount->future(futureExpiry,
                                              futureEnd,
                                              daycount.get(),
                                              CompoundBasis::SIMPLE,
                                              vol);

            double pv = discount->pv(optionExpiry);

            // quoted as options on future's PRICE - here we have an
            // option on future's RATE so need to invert call/put
            double value = Black::price(!isCall,
                                        futRate,
                                        k,
                                        pv,
                                        variance);

            results->storePrice(value, discount->getCcy());       

            if (control && control->isPricing() ) {               
                OutputRequest* request = control->requestsOutput(OutputRequest::IND_VOL);
                if (request) {
                    results->storeRequestResult(request, vol);
                }
                
                request = control->requestsOutput(OutputRequest::FWD_AT_MAT);
                if (request) {
                    CDoubleSP     fwd(CDouble::create(futRate));
                    OutputNameSP  name (new OutputName(discount->getCcy()));
                    results->storeRequestResult(request, fwd, name);
                }
                
                requests(control, results);
            }    
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    
    DateTime endDate(const Sensitivity* sensControl) const {
        MaturityPeriod  futurePeriod(period);
        return futurePeriod.toDate(futureExpiry);
    }

    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP outputName,
        const IModel*      model) const {
        static const string method = "OptionOnIRF::getSensitiveIRVolPoints";
        try {
            DateTimeArray      dates(1, optionExpiry);
            MaturityPeriod tenor(period);
            SwapMaturityVolRequest volRequest(&tenor);
            return VolProcessedBSIR::sensitiveIRVolPoints(vol.get(), 
                                                          discount.get(),
                                                          &volRequest, 
                                                          dates);
#if 0
            DateTimeArray      dates(1, optionExpiry);
            IRGridPointArraySP points(new IRGridPointArray(0));
            MaturityPeriodSP   tenor(new MaturityPeriod(period));
            CVolRequestSP      volRequest(new SwapMaturityVolRequest(tenor.get()));
            
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

    virtual bool priceDeadInstrument(CControl* control, CResults* results) const {       
        if (optionExpiry.isGreater(valueDate)) {
            return false;
        }
        else if (valueDate.isGreater(optionExpiry)) {
            results->storePrice(0.0, discount->getCcy());
            return true;   
        }            
        else {
            double premium;
            // inverted put/call logic
            if (!isCall) {
                premium = Maths::max(fix - k, 0.0);
            }
            else {
                premium = Maths::max(k - fix, 0.0);
            }            
            results->storePrice(premium, discount->getCcy());  

            if (control && control->isPricing()) {                            
                requests(control, results);
            }        
            return true;
        }
    }   

    virtual bool sensShift(Theta* shift) {
        try {
            DateTime newDate = shift->rollDate(valueDate);
            if ((newDate >= optionExpiry && valueDate < optionExpiry ) ||
                (valueDate == optionExpiry && Maths::isZero(fix))) {
                DayCountConventionSP daycount(DayCountConventionFactory::make(dcc));
                MaturityPeriod       futurePeriod(period);
                DateTime             futureEnd = futurePeriod.toDate(futureExpiry);
        
                // don't care about vol for tomorrow's rate
                fix = discount->fwd(futureExpiry,
                                    futureEnd,
                                    daycount.get(),
                                    CompoundBasis::SIMPLE);
            }  
            // roll the parent (updates value date etc)
            GenericSimpleIR::sensShift(shift);
        }
        catch (exception& e) {
            throw ModelException(e, "OptionOnIRF::sensShift (theta)");
        }    
        return true; // our components have theta type sensitivity
    }

/** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;
    ClosedFormIRLN::IProduct* createProduct(ClosedFormIRLN* model) const;

private:
    bool     isCall;
    double   strike;
    string   period;
    string   dcc;
    double   fixing;
    DateTime optionExpiry;
    DateTime futureExpiry;

    double   k;
    double   fix;
};

CClassConstSP const OptionOnIRF::TYPE = CClass::registerClassLoadMethod(
    "OptionOnIRF", typeid(OptionOnIRF), OptionOnIRF::load);


/** private class */
class OptionOnIRFClosedForm: public virtual CClosedFormLN::IProduct,
                             public virtual ClosedFormIRLN::IProduct {
private:
    const OptionOnIRF* oif; // a reference

public:
    OptionOnIRFClosedForm(const OptionOnIRF* oif): oif(oif){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results)const{
        oif->price(control, results);
    }

    void price(ClosedFormIRLN* model,
               Control*        control, 
               CResults*       results) const{
        oif->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* OptionOnIRF::createProduct(CClosedFormLN* model) const
{
    return new OptionOnIRFClosedForm(this);
}
ClosedFormIRLN::IProduct* OptionOnIRF::createProduct(ClosedFormIRLN* model) const
{
    return new OptionOnIRFClosedForm(this);
}

// for class loading 
bool OptionOnIRFLoad() {
    return (OptionOnIRF::TYPE != 0);
}

DRLIB_END_NAMESPACE
