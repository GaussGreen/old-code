//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Shrink.cpp
//
//   Description : Shrink  - figure this out later
//
//   Author      : Andrew J Swain
//
//   Date        : 20 November 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNAssetAndImnt.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Format.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/CashSettleDate.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

class Shrink: public GenericNAssetAndImnt,
              virtual public CClosedFormLN::IIntoProduct,
              virtual public ISensitiveStrikes,
              virtual public Theta::Shift,
              virtual public LastSensDate {
public:
    static CClassConstSP const TYPE; 

    virtual void Validate() {
        static const string method("Shrink::Validate");
        try {
            int i;
            int j;
            // let's not go there shall we
            if (instSettle->isMargin()) {
                throw ModelException(method, "margin settlement not supported");
            }
            if (instSettle->isPhysical()) {
                throw ModelException(method, "physical settlement not supported");
            }

            if (CashSettleDate::TYPE->isInstance(instSettle.get())) {
                throw ModelException(method, "cash settle date not supported");
            }

            if (feeDates.empty()) {
                throw ModelException(method, "no fee dates supplied");
            }

            // ensure dates for past values are increasing and distinct
            for (i = 1; i < pastNAV.size(); i++) {
                if (pastNAV[i].date <= pastNAV[i-1].date) {
                    throw ModelException(method, "sample dates must be in "
                                         "strictly increasing order (" + 
                                         pastNAV[i].date.toString() + 
                                         ") is on or before (" + 
                                         pastNAV[i-1].date.toString() +")");
                }
            }

            // ditto fee dates
            for (i = 1; i < feeDates.size(); i++) {
                if (feeDates[i] <= feeDates[i-1]) {
                    throw ModelException(method, "fee dates must be in "
                                         "strictly increasing order (" + 
                                         feeDates[i].toString() + 
                                         ") is on or before (" + 
                                         feeDates[i-1].toString() +")");
                }
            }

            if (pastNAV[pastNAV.size()-1].date > 
                feeDates[feeDates.size()-1]) {
                throw ModelException(method, "last sample date (" + 
                                     pastNAV[pastNAV.size()-1].date.toString() + 
                                     ") is after last fee date (" + 
                                     feeDates[feeDates.size()-1].toString() + ")");
            }

            samplesPerPeriod.resize(feeDates.size());

            // sample dates can't coincide with pay dates, and need at least
            // as many pay dates as samples
            if (feeDates.size() > pastNAV.size()) {
                throw ModelException(method, "more fee dates (" + 
                                     Format::toString(feeDates.size()) + 
                                     ") than sample dates (" + 
                                     Format::toString(pastNAV.size()));
            }

            for (i = 0, j = 0; i < feeDates.size(); i++) {
                int numSamples = 0;
                for (;j<pastNAV.size() && pastNAV[j].date <= feeDates[i];
                     j++) {
                    if (feeDates[i].equals(pastNAV[j].date)) {
                        throw ModelException(method,
                                             "fee on ("+feeDates[i].toString()+
                                             ") can't be a sample date too");
                    }
                    numSamples++;
                }
                samplesPerPeriod[i] = numSamples;
            }  

            // work out when each fee is paid
            payDates.resize(feeDates.size());
            for (i = 0; i < feeDates.size(); i++) {
                // can't be physical settlement, so can use null for asset
                payDates[i] = instSettle->settles(feeDates[i],0);
            }
        }    
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model) {
        ISensitiveStrikes* vegaMatrixImnt = 
            dynamic_cast<ISensitiveStrikes*>(imnt.get());

        if (vegaMatrixImnt) {
            return vegaMatrixImnt->avoidVegaMatrix(model);
        }
        else {
            return true;  // can't do it
        }
    }
        

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model) {
        if (avoidVegaMatrix(model)) {
            throw ModelException("Shrink::getSensitiveStrikes", 
                                 "VEGA_MATRIX is not valid for this instrument");
        }

        ISensitiveStrikes& strikes = dynamic_cast<ISensitiveStrikes&>(*imnt);

        return strikes.getSensitiveStrikes(outputName, model);
    }

    virtual bool sensShift(Theta* shift) {
        static const string method = "Shrink::sensShift";
        try {
            // calculate the new date
            DateTime newDate = shift->rollDate(valueDate);
            bool     needToPrice = true;
            bool     pastRollDate = false;
            double   value;
            int      i = 0;

            // set any past values to the value of the shrink itself
            while (!pastRollDate && i < pastNAV.size()) {
                if (pastNAV[i].date >  valueDate &&
                    pastNAV[i].date <= newDate) {
                    if (needToPrice) {
                        value = price();
                        needToPrice = false;
                    }
                    pastNAV[i].amount = value;
                }
                i++;
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }    
        return GenericNAssetAndImnt::sensShift(shift);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Shrink, clazz);
        SUPERCLASS(GenericNAssetAndImnt);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultShrink);
        FIELD(fee, "fee");
        FIELD(feeDates, "feeDates");
        FIELD(pastNAV, "pastNAV");
        FIELD(subtractUnderlyer, "subtractUnderlyer");
        FIELD_MAKE_OPTIONAL(subtractUnderlyer);
        FIELD_NO_DESC(samplesPerPeriod);
        FIELD_MAKE_TRANSIENT(samplesPerPeriod);
        FIELD_NO_DESC(payDates);
        FIELD_MAKE_TRANSIENT(payDates);
    }

    static IObject* defaultShrink(){
        return new Shrink();
    }
    
private:
    friend class ShrinkPricer;

    Shrink():GenericNAssetAndImnt(TYPE),fee(0.0),subtractUnderlyer(false) {}; 
    Shrink(const Shrink& rhs);
    Shrink& operator=(const Shrink& rhs);

    void recordOutputRequests(Control* control, 
                              Results* results) const {
        try {
            OutputRequest* request = NULL;
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                OutputRequestUtil::recordPaymentDates(control,results,&payDates); 
            }
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                CashFlowArray cfl;
                bool complete = false;
                int  prevIdx = 0;
                int  sampleIdx = 0;
                for (int i = 0; i < feeDates.size() && !complete; i++) {
                    // what's the last sample date for this fee
                    sampleIdx += samplesPerPeriod[i];
                    DateTime lastSample = pastNAV[sampleIdx-1].date;
                    if (lastSample <= valueDate) {
                        double pays = 0.0;
                        for (int j = 0; j < samplesPerPeriod[i]; j++) {
                            pays += pastNAV[prevIdx + j].amount;
                        }
                        pays *= fee/samplesPerPeriod[i];

                        CashFlow cf(payDates[i], pays);
                        cfl.push_back(cf);

                        prevIdx = sampleIdx;
                    }
                    else {
                        complete = true;
                    }
                }
                OutputRequestUtil::recordKnownCashflows(control,
                                                        results,
                                                        discount->getCcy(),
                                                        &cfl);   
            }

        }
        catch (exception) {
            // don't fail pricing for any of this
        }
    }    

    // the model interface
    void price(IModel* feeModel,Control* control,CResults* results) const {
        static const string method = "Shrink::price";
        try  {
            double value = price();

            results->storePrice(value, discount->getCcy());  

            if (control && control->isPricing()) {
                recordOutputRequests(control, results);
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // all we need to actually price this
    double price() const{
        static const string method = "Shrink::price";
        try  {
            double fvImnt = 0.0;
            double value = 0.0;
            double overallGrossUp = 1.0;
            double accrued = 0.0;
            int i;
            int j;

            // if there's a future payment, need value of underlying imnt
            if (payDates[payDates.size()-1] > valueDate) {
                CControlSP ctrl(Control::makeFromFlags("", 0.0));
                IModel*    model = imntModel.get();
                
                // 'Run' should route through control as appropriate
                ResultsSP baseResults(model->Run(imnt.get(), ctrl.get()));
                fvImnt = baseResults->retrievePrice();
            }                

            for (i = 0, j = 0; i < feeDates.size(); i++) {
                // for each sample in fee period
                double knownFee = 0.0;
                double sumGrowth = 0.0;
                for (; 
                     j < pastNAV.size() && pastNAV[j].date < payDates[i];
                     j++) {
                    if (payDates[i] <= valueDate) {
                        // do nothing, it's paid
                    }
                    else if (pastNAV[j].date <= valueDate) {
                        // known sample
                        knownFee+= pastNAV[j].amount;
                    }
                    else {
                        // sample and fee must be in the future
                        // imnt value in the future is simply the forward 
                        // value of the price, so just compute the sum of 
                        // the growth factors
                        sumGrowth += 1.0/discount->pv(pastNAV[j].date);
                    }
                }
                // work out any accrued fee
                knownFee /= samplesPerPeriod[i];
                knownFee *= discount->pv(payDates[i])*fee;
                accrued += knownFee;
                    
                // and the gross up for any future fees
                sumGrowth /= samplesPerPeriod[i];               
                double grossUp = 1.0+sumGrowth*discount->pv(payDates[i])*fee; 
                overallGrossUp *= grossUp;                
            }

            value = accrued + fvImnt*overallGrossUp;
            
            if (subtractUnderlyer) {
                value -= fvImnt;
            }
            return value;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    
    DateTime endDate(const Sensitivity* sensControl) const {
        LastSensDate* lsd = dynamic_cast<LastSensDate*>(imnt.get());

        if (lsd) {
            // always use the imnt's model to figure out the end date
            // rationale - the instrument knows its end date when it prices
            // and this avoid problems where the instrument end date only
            // works with its own model (i.e. GenericNFBase, I thank you)
            SensitivitySP shift(sensControl->spawn(imntModel.get()));

            return lsd->endDate(shift.get());
        }

        throw ModelException("Shrink::endDate",
                             "underlying imnt (" + 
                             imnt->getClass()->getName() + 
                             ") does not implement LastSensDate");
    }

    /** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

private:
    double        fee;
    DateTimeArray feeDates;
    CashFlowArray pastNAV;
    bool          subtractUnderlyer;

    // derived data
    IntArray      samplesPerPeriod;
    DateTimeArray payDates;

};

CClassConstSP const Shrink::TYPE = CClass::registerClassLoadMethod(
    "Shrink", typeid(Shrink), Shrink::load);


/** private class */
class ShrinkPricer: virtual public CClosedFormLN::IProduct{
private:
    const Shrink* opt; // a reference

public:
    ShrinkPricer(const Shrink* opt): opt(opt){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results) const{
        opt->price(model, control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* Shrink::createProduct(CClosedFormLN* model) const
{
    return new ShrinkPricer(this);
}


// for class loading 
bool ShrinkLoad() {
    return (Shrink::TYPE != 0);
}

DRLIB_END_NAMESPACE
