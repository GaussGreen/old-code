//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TimingRuleTester.cpp
//
//   Description : N-factor that tests EAS mapping of timing rules
//
//   Author      : Andrew J Swain
//
//   Date        : 26 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericNFactor.hpp"
#include "edginc/Results.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/MultiFactors.hpp"

DRLIB_BEGIN_NAMESPACE

class TimingRuleTester: public GenericNFactor, 
                        virtual public CClosedFormLN::IIntoProduct
{
private:
    /// fields ////////
    DateTimeArrayArray inDates;
    CashFlowArray      cfl; 
    StringArray        timeOfDate;
    StringArray        timeOfCFL;
    bool               isCall;

public:
    static CClassConstSP const TYPE;
    friend class GTRClosedForm;

    // validation
    void validatePop2Object(){
        int numAssets = assets->NbAssets();

        if (inDates.size() != numAssets) {
            throw ModelException("TimingRuleTester::validatePop2Object",
                                 "mismatch in dates & assets");
        }

        if (cfl.size() != numAssets) {
            throw ModelException("TimingRuleTester::validatePop2Object",
                                 "mismatch in cashflows & assets");
        }
    }  

    /** Implementation of ClosedFormLN::IntoProduct interface - the
        implementation of this is below */
    virtual CClosedFormLN::IProduct* createProduct(
        CClosedFormLN* model) const;

    void price(CResults* results) const {
        static const string method = "TimingRuleTester::price";
        try {
            double value = 0.0;

            for (int i = 0; i < assets->NbAssets(); i++) {
                string timeOfDay = DateTime::timeFormat(inDates[i][0].getTime());
                if (timeOfDay != timeOfDate[i]) {
                    throw ModelException(method,
                                         "for asset (" + assets->getName(i) + 
                                         ") expected " + timeOfDate[i] + 
                                         " time but got " + timeOfDay + 
                                         " for inDates");
                }

                timeOfDay = DateTime::timeFormat(cfl[i].date.getTime());
                if (timeOfDay != timeOfCFL[i]) {
                    throw ModelException(method,
                                         "for asset (" + assets->getName(i) + 
                                         ") expected " + timeOfCFL[i] + 
                                         " time but got " + timeOfDay + 
                                         " for cfl");
                }
            }

            value = isCall ? 1.0 : -1.0;
            value *= notional;
            results->storePrice(value, discount->getCcy());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }


private:
    TimingRuleTester(): GenericNFactor(TYPE), isCall(true) {} // for reflection
    TimingRuleTester(const TimingRuleTester& rhs); // not implemented
    TimingRuleTester& operator=(const TimingRuleTester& rhs); // not implemented

    static IObject* defaultTimingRuleTester(){
        return new TimingRuleTester();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TimingRuleTester, clazz);
        SUPERCLASS(GenericNFactor);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultTimingRuleTester);
        FIELD(inDates, "dates per asset");
        FIELD(cfl, "cashflows per asset");
        FIELD(timeOfDate, "timing of dates");
        FIELD(timeOfCFL, "timing of cashflows");
        FIELD(isCall, "is it a call option");
    }
};

/** private class */
class GTRClosedForm: public CClosedFormLN::IProduct{
private:
    const TimingRuleTester* gnf; // a reference

public:
    GTRClosedForm(const TimingRuleTester* gnf): gnf(gnf){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results) const{
        gnf->price(results);
    }
};


CClosedFormLN::IProduct* TimingRuleTester::createProduct(
    CClosedFormLN* model) const {
    return new GTRClosedForm(this);
}


CClassConstSP const TimingRuleTester::TYPE = CClass::registerClassLoadMethod(
    "TimingRuleTester", typeid(TimingRuleTester), TimingRuleTester::load);

// * for class loading (avoid having header file) */
bool TimingRuleTesterLoad() {
    return (TimingRuleTester::TYPE != 0);
}

DRLIB_END_NAMESPACE

