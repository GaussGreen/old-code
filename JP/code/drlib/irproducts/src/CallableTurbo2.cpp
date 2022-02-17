//----------------------------------------------------------------------------
//
//   Group       : QR&D Interest Rates
//
//   Filename    : CallableTurbo2.cpp
//
//   Description : instrument shell that makes use of components 
//                 with customized internal logic to fill optional data fields
//                 for each components. 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KOption.hpp"
#include "edginc/KFloatLeg.hpp"
#include "edginc/KRibTurbo.hpp"
#include "edginc/KSum.hpp"
#include "edginc/KCashflow.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

//*******************************  CallableTurbo2 *****************************/
// the shell instrument for CallableTurbo2
class CallableTurbo2 : public CInstrument,
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
    CallableTurbo2() : CInstrument(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new CallableTurbo2(); }

    /****************** exported fields ************/
    KOptionSP           option;
    KRibTurboSP         turboLeg;
    KRibTurboSP         fundingLeg;
    IProdCreatorSP      redeemer;
    YieldCurveWrapper   discount;
};

FDProductSP CallableTurbo2::createProduct(FDModel * model) const {
    return model->createProduct( option );
}

void CallableTurbo2::GetMarket(const IModel *model, const CMarketDataSP market) {
    try {
        if (fundingLeg->discount.getName().empty()) 
            fundingLeg->discount = discount;
		fundingLeg->outputName = "FUNDING_LEG";

        if (turboLeg->discount.getName().empty()) 
            turboLeg->discount = discount;
        turboLeg->outputName = "TURBO_LEG";

        if (turboLeg->legs->size()!=2)
            throw ModelException("turboLeg->legs->size()!=2");
        if (fundingLeg->legs->size()!=1)
            throw ModelException("fundLeg->legs->size()!=1");
        if (turboLeg->rib.get())
            throw ModelException("turboLeg->rib must be left empty");
        if (fundingLeg->rib.get())
            throw ModelException("fundLeg->rib must be left empty");

        KSumSP sum(new KSum);
        sum->add(fundingLeg, 1.);   
        sum->add(turboLeg, 1.);   
        sum->discount = discount;
        sum->outputName = "UNDERLYING";

        option->und = sum;

        if (!option->discount.getName().empty()) 
            throw ModelException("option.discount should not be provided. "
                                 "Use CallableTurbo2.discount instead.");
        option->discount = discount;
        option->outputName = "OPTION";

        if (redeemer.get()) {

            // at this stage only fixed cashflow ot redeemer option supported
            KCashflow* fixedRedeemer = dynamic_cast<KCashflow*>(redeemer.get());
            KOption* redeemerOption = dynamic_cast<KOption*>(redeemer.get());

            if (!(fixedRedeemer || redeemerOption))
                throw ModelException("Only redeemer of type KCashflow or KOption "
                                     "is currently supported - type supplied = " +
                                     redeemer->getClass()->getName());

            if (fixedRedeemer) {
                if (fixedRedeemer->discount.getName().empty())
                    throw ModelException("Redeemer component must supply a discount curve "
                                         "as it often pays in a different currency");

                if (fixedRedeemer->outputName.empty())
                    fixedRedeemer->outputName = "REDEEMER";
            }
            else {
                if (redeemerOption->discount.getName().empty())
                    throw ModelException("Redeemer component must supply a discount curve "
                                         "as it often pays in a different currency");

                if (redeemerOption->outputName.empty())
                    redeemerOption->outputName = "REDEEMER";

                // ??? check option type is actually REDEEMER
            }

            // simply add redeemer to the underlying sum component
            sum->add(redeemer, 1.0);
        }

        // cleanup fields so that iterator do not get lost
/*        fundingLeg.reset(0);
        turboLeg.reset(0);
        if (redeemer.get())
            redeemer.reset(0);*/

        option->GetMarket(model, market);
        discount.getData(model, market);

        // check component configurations and cross check components
        string ccy = discount->getCcy();
        string ccy1 = (*turboLeg->legs)[0]->discount->getCcy();
        string ccy2 = (*turboLeg->legs)[1]->discount->getCcy();
        if (ccy1==ccy && ccy2==ccy)
            throw ModelException("The two turbo legs are domestic, change one to foreign");
        if (ccy1!=ccy && ccy2!=ccy)
            throw ModelException("The two turbo legs are foreign, change one to domestic");
	}
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// export the members of the class through the library interface
void CallableTurbo2::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CallableTurbo2, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(option,"");
    FIELD(fundingLeg,"");
    FIELD(turboLeg,"");
    FIELD(redeemer,"");  FIELD_MAKE_OPTIONAL(redeemer);
    FIELD(discount,"");
    Addin::registerConstructor(Addin::UTILITIES, CallableTurbo2::TYPE);
}

// register the CallableTurbo2 in the library framework at startup
CClassConstSP const CallableTurbo2::TYPE =  CClass::registerClassLoadMethod(
            "CallableTurbo2", typeid(CallableTurbo2), CallableTurbo2::load);

// to ensure linker doesn't optimise out the class
bool CallableTurbo2Load() {return CallableTurbo2::TYPE != 0;}

DRLIB_END_NAMESPACE
