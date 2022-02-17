//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BadBoy.cpp
//
//   Description : Class that always leaks, has access errors
//                 and differences to prove that testcmp and
//                 purify are working
//
//   Author      : Andrew J Swain
//
//   Date        : 17 August 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BadBoy.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

void BadBoy::price(CResults* results) const{
    static const string method = "BadBoy::price";
    try {
        double x;

        // well the price will always change
        double value = clock();

        // allocate some memory and abuse it
        int size = 5;
        for (int i = 0; i < 10; i++){
            double* leaker = new double[size];

            for (int j = 0; j < size + 1; j++) {
                x = leaker[j];
            }

            value += x;
            leaker = 0;
        }
        results->storePrice(value, "XYZ");
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/** Returns the value date (aka today) the instrument is currently
    pricing for */
DateTime BadBoy::getValueDate() const{
    return DateTime(); // doesn't matter
}



/** Returns the name of the instrument's discount currency. */
string BadBoy::discountYieldCurveName() const {
    return "ABC"; // doesn't matter
}


class BadBoyClosedForm: public ClosedForm::IProduct{
private:
    const BadBoy* bb; // a reference

public:
    BadBoyClosedForm(const BadBoy* bb): bb(bb){}

    void price(ClosedForm* model,
               Control*    control,
               CResults*   results) const{
        bb->price(results);
    }
};

ClosedForm::IProduct* BadBoy::createProduct(
    ClosedForm* model) const{
    return new BadBoyClosedForm(this);
}

BadBoy::BadBoy(): CInstrument(TYPE){}

class BadBoyHelper{
public:
/** Invoked when this class is 'loaded' */
    static void load(CClassSP& clazz){
        //clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BadBoy, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(ClosedForm::IIntoProduct);
        EMPTY_SHELL_METHOD(defaultBadBoy);
    }
    static IObject* defaultBadBoy(){
        return new BadBoy();
    }
};

CClassConstSP const BadBoy::TYPE = CClass::registerClassLoadMethod(
    "BadBoy", typeid(BadBoy), BadBoyHelper::load);




DRLIB_END_NAMESPACE

