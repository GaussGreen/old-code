//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CashFlowStream.cpp
//
//   Description : CashFlow instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 16 February 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CashFlowStream.hpp"

DRLIB_BEGIN_NAMESPACE


/** Support for ITaxableInstBasic */
const DateTime CashFlowStream::getFinalPaymentDate() const {
    return cfl->back().date;
}

/** Support for ITaxableInstWithCoupons */
CashFlowArrayConstSP CashFlowStream::getCoupons() const {
    return cfl;
}


/** private class */
class CashFlowStreamClosedForm: public ClosedForm::IProduct{
private:
    const CashFlowStream* cf; // a reference

public:
    CashFlowStreamClosedForm(const CashFlowStream* cf): cf(cf){}

    void price(ClosedForm* model,
               Control*    control, 
               CResults*   results) const{
        cf->price(control, results);
    }
};
    

/** Constructor only needed for demo purposes */
CashFlowStream::CashFlowStream(CashFlowArray* cfl,
                               const string&  discountName):
    SimpleCashFlowStream(TYPE, cfl, discountName)
{
    // empty
}

// for reflection
CashFlowStream::CashFlowStream(): SimpleCashFlowStream(TYPE) {}

class CashFlowStreamHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CashFlowStream, clazz);
        SUPERCLASS(SimpleCashFlowStream);
        EMPTY_SHELL_METHOD(defaultCashFlowStream);
        //FIELD(cfl, "cash flows");
        //FIELD(discount, "identifies discount curve");
        //FIELD(valueDate, "valuation date");
        //FIELD_MAKE_OPTIONAL(valueDate);
    }

    static IObject* defaultCashFlowStream(){
        return new CashFlowStream();
    }
};

CClassConstSP const CashFlowStream::TYPE = CClass::registerClassLoadMethod(
    "CashFlowStream", typeid(CashFlowStream), CashFlowStreamHelper::load);

bool  CashFlowStreamLoad() {
    return (CashFlowStream::TYPE != 0);
   }

   

DRLIB_END_NAMESPACE
