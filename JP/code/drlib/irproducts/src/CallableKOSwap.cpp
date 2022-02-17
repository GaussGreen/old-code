//----------------------------------------------------------------------------
//
//   Group       : QR&D Interest Rates
//
//   Filename    : CallableKOSwap.cpp
//
//   Description : instrument shell that makes use of components 
//                 with customized internal logic to fill optional data fields
//                 for each components. 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KOption.hpp"
#include "edginc/KKOSwap.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

// the shell instrument for CallableKOSwap
// It looks like we could factor out a new base class for "shell instruments"
// that are based on the KComponent architecture and move common functionality
// into the new base.
class CallableKOSwap : public CInstrument,
                    virtual public FDModel::IIntoProduct,
					virtual public LastSensDate
{
public:
    static CClassConstSP const TYPE;

    /* CInstrument:: */
    virtual void Validate(void) {}

    virtual DateTime getValueDate() const {
        return compRoot->getValueDate();
    }
    virtual string discountYieldCurveName() const {
        return discount.getName();
    }
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const {
        return compRoot->priceDeadInstrument(control, results);
    }
    virtual void GetMarket(const IModel *model, const CMarketDataSP market);

    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const{
        return model->createProduct( compRoot );
    }

	/* LastSensDate */
	virtual DateTime endDate(const Sensitivity* sensControl) const {
		return compRoot->endDate(sensControl);
	}

private:
    CallableKOSwap() : CInstrument(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new CallableKOSwap(); }

    /****************** exported fields ************/
    KOptionSP           option;
    KKOSwapSP           swapLeg;
    YieldCurveWrapper   discount;

    /****************** transiend fields ************/
    KComponentSP        compRoot;
};

// export the members of the class through the library interface
void CallableKOSwap::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CallableKOSwap, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(option, "option component");
    FIELD(swapLeg, "swapLeg component");
    FIELD(discount,"discount curve");
    FIELD(compRoot,""); FIELD_MAKE_TRANSIENT(compRoot);
    Addin::registerConstructor(Addin::UTILITIES, CallableKOSwap::TYPE);
}

// We link the CallableKOSwap components up in GetMarket
void CallableKOSwap::GetMarket(const IModel* model, const CMarketDataSP market) {
    try{

        if (swapLeg->outputName.empty())
            swapLeg->outputName = "SWAP_LEG_VALUE";
        if (swapLeg->discountYieldCurveName().empty()) {
            swapLeg->discount=discount;
        }

        option->und = swapLeg;

        compRoot = option;

        if (compRoot->discountYieldCurveName().empty())
            compRoot->discount = discount;

        // cleanup fields so that iterator do not get lost
        option.reset(0);
        swapLeg.reset(0);

        compRoot->GetMarket(model, market);
    }
    catch (exception& e) {
        throw ModelException(e, "CallableKOSwap::GetMarket");
    }
}

// register the CallableKOSwap in the library framework at startup
CClassConstSP const CallableKOSwap::TYPE =  CClass::registerClassLoadMethod(
            "CallableKOSwap", typeid(CallableKOSwap), CallableKOSwap::load);

// to ensure linker doesn't optimise out the class
bool CallableKOSwapLoad(void) {
    return (CallableKOSwap::TYPE != 0);
}

DRLIB_END_NAMESPACE
