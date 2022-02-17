//----------------------------------------------------------------------------
//
//   Group       : QR&D Interest Rates
//
//   Filename    : CallableTurbo.cpp
//
//   Description : instrument shell that makes use of components 
//                 with customized internal logic to fill optional data fields
//                 for each components. 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KOption.hpp"
#include "edginc/KFloatLeg.hpp"
#include "edginc/KTurbo.hpp"
#include "edginc/KSum.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

//*******************************  CallableTurbo *****************************/
// the shell instrument for CallableTurbo
class CallableTurbo : public CInstrument,
                      virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;

    /*********************** methods **********************/
    /* CInstrument:: */
    virtual void Validate(void) {}

    virtual DateTime getValueDate() const {
        return option->getValueDate();
    }
    virtual string discountYieldCurveName() const {
        return discount.getName();
    }
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const {
        return option->priceDeadInstrument(control, results);
    }
    virtual void GetMarket(const IModel *model, const CMarketDataSP market);
    virtual FDProductSP createProduct(FDModel * model) const;

private:
    CallableTurbo() : CInstrument(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new CallableTurbo(); }

    /****************** exported fields ************/
    KOptionSP           option;
    KFloatLegSP         floatLeg;
    KTurboSP            turboLeg;
    YieldCurveWrapper   discount;
};

FDProductSP CallableTurbo::createProduct(FDModel * model) const {
    return model->createProduct( option );
}

void CallableTurbo::GetMarket(const IModel *model, const CMarketDataSP market) {
    try {
        if (floatLeg->discount.getName().empty()) 
            floatLeg->discount = discount;
        floatLeg->outputName = "FUNDING_LEG";

        if (turboLeg->discount.getName().empty()) 
            turboLeg->discount = discount;
        turboLeg->outputName = "TURBO_LEG";

        if (turboLeg->legs.size()!=2)
            throw ModelException("The turbo component must have two legs");

        KSumSP sum(new KSum);
        sum->add(floatLeg, 1.);   
        sum->add(turboLeg, 1.);   
        sum->discount = discount;
        sum->outputName = "UNDERLYING";

        option->und = sum;

        if (!option->discount.getName().empty()) 
            throw ModelException("option.discount should not be provided. "
                                 "Use CallableTurbo.discount instead.");

        option->discount = discount;
        option->outputName = "Option";

        KTurboSP turbo = turboLeg;
        // cleanup fields so that iterator do not get lost
        floatLeg.reset(0);
        turboLeg.reset(0);

        option->GetMarket(model, market);
        discount.getData(model, market);

        {   // check that there is one domestic and one foreign turbo leg
            string ccy = discount->getCcy();
            string ccy1 = turbo->legs[0]->discount->getCcy();
            string ccy2 = turbo->legs[1]->discount->getCcy();
            if (ccy1==ccy && ccy2==ccy)
                throw ModelException("The two turbo legs are domestic, change one to foreign");
            if (ccy1!=ccy && ccy2!=ccy)
                throw ModelException("The two turbo legs are foreign, change one to domestic");
        }
    }
    catch (exception& e) {
        throw ModelException(e, "CallableTurbo::GetMarket");
    }
}

// export the members of the class through the library interface
void CallableTurbo::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CallableTurbo, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(option,"");
    FIELD(floatLeg,"");
    FIELD(turboLeg,"");
    FIELD(discount,"");
    Addin::registerConstructor(Addin::UTILITIES, CallableTurbo::TYPE);
}

// register the CallableTurbo in the library framework at startup
CClassConstSP const CallableTurbo::TYPE =  CClass::registerClassLoadMethod(
            "CallableTurbo", typeid(CallableTurbo), CallableTurbo::load);

// to ensure linker doesn't optimise out the class
bool CallableTurboLoad() {return CallableTurbo::TYPE != 0;}

DRLIB_END_NAMESPACE
