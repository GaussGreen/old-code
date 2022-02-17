//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Generic1FactorTester.cpp
//
//   Description : Public interface contains all current IMS supported
//                 data types + templates of templates
//
//   Author      : Andrew J Swain
//
//   Date        : 26 November 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/Results.hpp"
#include "edginc/Schedule.hpp"

DRLIB_BEGIN_NAMESPACE

class Generic1FactorTester: public Generic1Factor,
                            public CClosedFormLN::IIntoProduct {
public:
    static CClassConstSP const TYPE; 


    virtual void validatePop2Object(){
        static const string method = "Generic1FactorTester::validatePop2Object";
        

    }

    virtual void Validate(){
        static const string method = "Generic1FactorTester::Validate";

        if (american && 
            (toggles.size() != barrierLevels.size()   ||
             toggles.size() != monitorDates->size()   ||
             toggles.size() != annuityCoupons->size())) {
            throw ModelException("Generic1FactorTester::Validate",
                                 "lists must be same length for american");
        }
    }

/** Invoked when Class is 'loaded' */
static void load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Generic1FactorTester, clazz);
    SUPERCLASS(Generic1Factor);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultGeneric1FactorTester);
    FIELD(american, "american");
    FIELD(smoothing, "smoothing");
    FIELD(strike, "strike");
    FIELD(swapNotional, "swapNotional");
    FIELD(maturity, "maturity");
    FIELD(upDown, "upDown");
    FIELD(toggles, "toggles");
    FIELD(barrierLevels, "barrierLevels");
    FIELD(monitorDates, "monitorDates");
    FIELD(annuityCoupons, "annuityCoupons");
    FIELD(avgOut, "avgOut");
    FIELD(lowRangeSchedule, "lowRangeSchedule");
}

static IObject* defaultGeneric1FactorTester(){
    return new Generic1FactorTester();
}
    
private:
    friend class Generic1FactorTesterClosedForm;

    Generic1FactorTester():Generic1Factor(TYPE) {}; 
    Generic1FactorTester(const Generic1FactorTester& rhs);
    Generic1FactorTester& operator=(const Generic1FactorTester& rhs);

    void price(Control* control, CResults* results)const{
        static const string method = "Generic1FactorTester::price";
        
        try  {
            double value = 0.0;

            if (!american) {
                double fwd = asset->fwdValue(maturity);

                value = (fwd-strike)/swapNotional;
            }
            else {
                int i;
                double fwd;

                for (i = 0; i < monitorDates->size(); i++) {
                    fwd = asset->fwdValue((*monitorDates)[i]);
                    value += fwd * lowRangeSchedule->interpolate((*monitorDates)[i]);
                    value += toggles[i] * (*annuityCoupons)[i].amount;

                    if (upDown == "UP") {
                        value += barrierLevels[i];
                    }
                    else {
                        value -= barrierLevels[i];
                    }

                    value += smoothing;
                }

                value += avgOut->futureSampleSum(asset.get(), valueDate);

                value -= strike;
                value /= swapNotional;                      
            }
            results->storePrice(value, discount->getCcy());           
        }
        catch (exception& e) {
            throw ModelException(&e, method);
        }
    }
    

/** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

private:
    bool            american;
    int             smoothing;
    double          strike;
    double          swapNotional;
    DateTime        maturity;
    string          upDown;
    IntArray        toggles;
    DoubleArray     barrierLevels;
    DateTimeArraySP monitorDates;
    CashFlowArraySP annuityCoupons;
    SampleListSP    avgOut;
    ScheduleSP      lowRangeSchedule;
};


CClassConstSP const Generic1FactorTester::TYPE = CClass::registerClassLoadMethod(
    "Generic1FactorTester", typeid(Generic1FactorTester), Generic1FactorTester::load);


/** private class */
class Generic1FactorTesterClosedForm: public CClosedFormLN::IProduct{
private:
    const Generic1FactorTester* gft; // a reference

public:
    Generic1FactorTesterClosedForm(const Generic1FactorTester* gft): gft(gft){}

    void price(CClosedFormLN* model,
               Control*    control, 
               CResults*   results)const{
        gft->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* Generic1FactorTester::createProduct(CClosedFormLN* model) const
{
    return new Generic1FactorTesterClosedForm(this);
}

// for class loading 
bool Generic1FactorTesterLoad() {
    return (Generic1FactorTester::TYPE != 0);
}

DRLIB_END_NAMESPACE
