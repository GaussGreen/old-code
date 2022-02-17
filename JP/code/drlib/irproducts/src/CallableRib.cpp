//----------------------------------------------------------------------------
//
//   Group       : QR&D Interest Rates
//
//   Filename    : CallableRib.cpp
//
//   Description : instrument shell that makes use of components 
//                 with customized internal logic to fill optional data fields
//                 for each components. 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KOption.hpp"
#include "edginc/KSum.hpp"
#include "edginc/KFloatLeg.hpp"
#include "edginc/KRibFloatLeg.hpp"
#include "edginc/KRibTurbo.hpp"
#include "edginc/KCashflow.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

// the shell instrument for CallableRib
// It looks like we could factor out a new base class for "shell instruments"
// that are based on the KComponent architecture and move common functionality
// into the new base.
class CallableRib : public CInstrument,
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
    CallableRib() : CInstrument(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new CallableRib(); }

    /****************** exported fields ************/
    KOptionSP           option;
    KRibFloatLegSP      floatLeg;
    KFloatLegSP         fixedLeg;
    KRibTurboSP         floatLeg2;
    KRibTurboSP         fixedLeg2;
    YieldCurveWrapper   discount;

    /****************** transiend fields ************/
    KComponentSP        compRoot;
};

// export the members of the class through the library interface
void CallableRib::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CallableRib, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(option, "option component");
    FIELD(floatLeg, "floatLeg component");   FIELD_MAKE_OPTIONAL(floatLeg);
    FIELD(fixedLeg, "fixedLeg component");   FIELD_MAKE_OPTIONAL(fixedLeg);
    FIELD(floatLeg2, "floatLeg component");  FIELD_MAKE_OPTIONAL(floatLeg2);
    FIELD(fixedLeg2, "fixedLeg component");  FIELD_MAKE_OPTIONAL(fixedLeg2);  
    FIELD(discount,"discount curve");
    FIELD(compRoot,""); FIELD_MAKE_TRANSIENT(compRoot);
    Addin::registerConstructor(Addin::UTILITIES, CallableRib::TYPE);
}

// We link the CallableRib components up in GetMarket
void CallableRib::GetMarket(const IModel* model, const CMarketDataSP market) {
    try{
        // TO DO: Make sure option exercise dates fall on Rib accrual dates or
        // the last Rib accrual end date.

        // We construct the internal KSum component which "sums" up the funding leg
        // and the floating (RIB) leg. 
        KSumSP sum(new KSum);
        sum->outputName="UNDERLYING";

        if (fixedLeg.get()) {
            if (fixedLeg2.get())
                throw ModelException("cannot use both fixedLeg and fixedLeg2 at the same time");
            sum->add(fixedLeg, 1.0);
            if (fixedLeg->outputName.empty())
                fixedLeg->outputName = "FIXED_LEG_VALUE";
            if (fixedLeg->discountYieldCurveName().empty())
                fixedLeg->discount=discount;
        }
        else {
            if (!fixedLeg2.get())
                throw ModelException("fixedLeg or fixedLeg2 must be provided");
             sum->add(fixedLeg2, 1.0);
            if (fixedLeg2->outputName.empty())
                fixedLeg2->outputName = "FIXED_LEG_VALUE";
            if (fixedLeg2->discountYieldCurveName().empty())
                fixedLeg2->discount=discount;
        }
    
        if (floatLeg.get()) {
            if (floatLeg2.get())
                throw ModelException("cannot use both floatLeg and floatLeg2 at the same time");

            if (floatLeg->outputName.empty())
                floatLeg->outputName = "FLOAT_LEG_VALUE";
            if (floatLeg->discountYieldCurveName().empty())
                floatLeg->discount=discount;
            sum->add(floatLeg,  1.0);
        }
        else {
            if (!floatLeg2.get())
                throw ModelException("floatLeg or floatLeg2 must be provided");

            if (floatLeg2->outputName.empty())
                floatLeg2->outputName = "FLOAT_LEG_VALUE";
            if (floatLeg2->discountYieldCurveName().empty())
                floatLeg2->discount=discount;
            sum->add(floatLeg2,  1.0);
        }

        option->und = sum;

        compRoot = option;

        if (compRoot->discountYieldCurveName().empty())
            compRoot->discount = discount;

        // cleanup fields so that iterator do not get lost
        option.reset(0);
        floatLeg.reset(0);
        floatLeg2.reset(0);
        fixedLeg.reset(0);
        fixedLeg2.reset(0);

        compRoot->GetMarket(model, market);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// register the CallableRib in the library framework at startup
CClassConstSP const CallableRib::TYPE =  CClass::registerClassLoadMethod(
            "CallableRib", typeid(CallableRib), CallableRib::load);

// to ensure linker doesn't optimise out the class
bool CallableRibLoad(void) {
    return (CallableRib::TYPE != 0);
}


// NEW version of the shell instrument for CallableRib
// It looks like we could factor out a new base class for "shell instruments"
// that are based on the KComponent architecture and move common functionality
// into the new base.
class CallableRib2 : public CInstrument,
                     virtual public FDModel::IIntoProduct,
                     virtual public LastSensDate {
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
    virtual FDProductSP createProduct(FDModel * model) const {
        return model->createProduct( compRoot );
    }

	/* LastSensDate */
	virtual DateTime endDate(const Sensitivity* sensControl) const {
		return compRoot->endDate(sensControl);
	}

private:
    CallableRib2() : CInstrument(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new CallableRib2(); }

    /****************** exported fields ************/
    KOptionSP           option;
    KRibTurboSP         complexLeg;
    KRibTurboSP         fundingLeg;
    IProdCreatorSP      redeemer;
    YieldCurveWrapper   discount;

    /****************** transiend fields ************/
    KComponentSP        compRoot;
};

// export the members of the class through the library interface
void CallableRib2::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CallableRib2, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    // exported fields
    FIELD(discount,"pricing discount curve for product");
    FIELD(option, "option component");
    FIELD(complexLeg, "RIB leg component");
    FIELD(fundingLeg, "funding leg component");
    FIELD(redeemer, "Redeemer component - either KCashflow or other component");
    FIELD_MAKE_OPTIONAL(redeemer);
    
    // transient fields
    FIELD(compRoot,""); 
    FIELD_MAKE_TRANSIENT(compRoot);

}

// Link the CallableRib components up in GetMarket
void CallableRib2::GetMarket(const IModel* model, const CMarketDataSP market) {
    try{

        KSumSP optionUnderlying(new KSum);
        optionUnderlying->outputName = "OPTION_UNDERLYING";

        // if redeemer is provided, we need 2 sum components
        // 1) complexRedeemer = complexLeg + redeemer
        // 2) OptionUnderlying = complexRedeemer + fundingLeg
        if (redeemer.get()) {

            // at this stage only fixed cashflow redeemer supported
            KCashflow* fixedRedeemer = dynamic_cast<KCashflow*>(redeemer.get());
            if (!fixedRedeemer)
                throw ModelException("Only redeemer of type KCashflow is currently "
                                     "supported - type supplied = " +
                                     redeemer->getClass()->getName());

            // sum 1
            KSumSP complexRedeemer(new KSum);
            complexRedeemer->outputName = "COMPLEX_REDEEMER";

            if (complexLeg->discountYieldCurveName().empty())
                complexLeg->discount = discount;
            complexRedeemer->add(complexLeg, 1.0);
            if (complexLeg->outputName.empty())
                complexLeg->outputName = "COMPLEX_LEG";
            complexRedeemer->add(redeemer, 1.0);
            if (fixedRedeemer->outputName.empty())
                fixedRedeemer->outputName = "REDEEMER";

            // sum 2
            if (optionUnderlying->discountYieldCurveName().empty())
                optionUnderlying->discount = discount;
            optionUnderlying->add(complexRedeemer, 1.0);
            optionUnderlying->add(fundingLeg, 1.0);   // ??? change this to +ve
            if (fundingLeg->outputName.empty())
                fundingLeg->outputName = "FUNDING_LEG";
        }
        else {
            // if redeemer not provided, just need 1 sum
            // 1) OptionUnderlying = complexLeg + fundingLeg
            if (optionUnderlying->discountYieldCurveName().empty())
                optionUnderlying->discount = discount;
                optionUnderlying->add(complexLeg, 1.0);
                optionUnderlying->add(fundingLeg, 1.0);   // ??? change this to +ve
                if (complexLeg->outputName.empty())
                    complexLeg->outputName = "COMPLEX_LEG";
                if (fundingLeg->outputName.empty())
                    fundingLeg->outputName = "FUNDING_LEG";
        }

        option->undList.resize(1);
        option->undList[0] = optionUnderlying;
        if (option->outputName.empty())
            option->outputName = "OPTION";
        
        compRoot = option;  // option is the top level component

        if (compRoot->discountYieldCurveName().empty())
            compRoot->discount = discount;

        // cleanup fields so that iterator do not get lost
        option.reset(0);
        complexLeg.reset(0);
        fundingLeg.reset(0);
        if (redeemer.get())
            redeemer.reset(0);

        compRoot->GetMarket(model, market);

        // Implement cross component checks
        // ??? toDo - need access to private members to do any sensible
        // checks, but mad to make all components friends of shell instruments
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// register the CallableRib2 in the library framework at startup
CClassConstSP const CallableRib2::TYPE =  CClass::registerClassLoadMethod(
            "CallableRib2", typeid(CallableRib2), CallableRib2::load);

// to ensure linker doesn't optimise out the class
bool CallableRib2Load(void) {
    return (CallableRib2::TYPE != 0);
}


DRLIB_END_NAMESPACE
