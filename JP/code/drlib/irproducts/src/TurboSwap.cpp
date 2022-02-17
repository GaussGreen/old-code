//----------------------------------------------------------------------------
//
//   Group       : QR&D Interest Rates
//
//   Filename    : TurboSwap.cpp
//
//   Description : instrument shell that makes use of components 
//                 with customized internal logic to fill optional data fields
//                 for each components. 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KFloatLeg.hpp"
#include "edginc/KRibTurbo.hpp"
#include "edginc/KSum.hpp"
#include "edginc/KCashflow.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

//*******************************  TurboSwap *****************************/
// the shell instrument for TurboSwap
class TurboSwap : public CInstrument,
                      virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;

    /*********************** methods **********************/
    /* CInstrument:: */
    virtual void Validate(void) {}

    virtual DateTime getValueDate() const {
        return topComponent->getValueDate();
    }
    virtual string discountYieldCurveName() const {
        return discount.getName();
    }
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const {
        return topComponent->priceDeadInstrument(control, results);
    }
    virtual void GetMarket(const IModel *model, const CMarketDataSP market);
    virtual FDProductSP createProduct(FDModel * model) const;

private:
    TurboSwap() : CInstrument(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new TurboSwap(); }

    /****************** exported fields ************/
    KRibTurboSP         fundLeg;
    KRibTurboSP         turboLeg;
    YieldCurveWrapper   discount;
    KCashflowSP         redeemer;
    /****************** transient fields ************/
    KComponentSP        topComponent;
};

FDProductSP TurboSwap::createProduct(FDModel * model) const {
    return model->createProduct( topComponent );
}

void TurboSwap::GetMarket(const IModel *model, const CMarketDataSP market) {
    try {
        if (fundLeg->discount.getName().empty()) fundLeg->discount = discount;
		fundLeg->outputName = "FUNDING_LEG";

        if (turboLeg->discount.getName().empty()) turboLeg->discount = discount;
        turboLeg->outputName = "TURBO_LEG";

        if (turboLeg->legs->size()!=2)
            throw ModelException("turboLeg->legs->size()!=2");
        if (fundLeg->legs->size()!=1)
            throw ModelException("fundLeg->legs->size()!=1");
        if (turboLeg->rib.get())
            throw ModelException("turboLeg->rib must be left empty");
        if (fundLeg->rib.get())
            throw ModelException("fundLeg->rib must be left empty");

        IProdCreatorSP funding;
        if (redeemer.get()) {
/*            if (redeemer->dates.size()!=1) 
                throw ModelException("redeemer->dates.size()!=1 (redeemer pays once)");*/

            KSumSP fundSum(new KSum);
            fundSum->add(fundLeg, 1.);   
            fundSum->add(redeemer, 1.);   
            fundSum->discount = discount;
            fundSum->outputName = "FUNDING_PLUS_REDEEMER";
            funding = fundSum;
        }
        else funding = fundLeg;

        KSumSP sum(new KSum);
        sum->add(funding, 1.);   
        sum->add(turboLeg, 1.);   
        sum->discount = discount;
        sum->outputName = "SWAP";

        topComponent = sum;

        if (!topComponent->discount.getName().empty()) throw ModelException(
            "topComponent.discount should not be provided. Use TurboSwap.discount instead.");
        topComponent->discount = discount;

        KRibTurboSP turbo = turboLeg;
        // cleanup fields so that iterator do not get lost
        fundLeg.reset(0);
        turboLeg.reset(0);

        topComponent->GetMarket(model, market);
        discount.getData(model, market);

        {   // check that there is one domestic and one foreign turbo leg
            string ccy = discount->getCcy();
            string ccy1 = (*turbo->legs)[0]->discount->getCcy();
            string ccy2 = (*turbo->legs)[1]->discount->getCcy();
            if (ccy1==ccy && ccy2==ccy)
                throw ModelException("The two turbo legs are domestic, change one to foreign");
            if (ccy1!=ccy && ccy2!=ccy)
                throw ModelException("The two turbo legs are foreign, change one to domestic");
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// export the members of the class through the library interface
void TurboSwap::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(TurboSwap, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(fundLeg,"");
    FIELD(turboLeg,"");
    FIELD(redeemer,"");  FIELD_MAKE_OPTIONAL(redeemer);
    FIELD(discount,"");
    FIELD(topComponent,""); FIELD_MAKE_TRANSIENT(topComponent);
    Addin::registerConstructor(Addin::UTILITIES, TurboSwap::TYPE);
}

// register the TurboSwap in the library framework at startup
CClassConstSP const TurboSwap::TYPE =  CClass::registerClassLoadMethod(
            "TurboSwap", typeid(TurboSwap), TurboSwap::load);

// to ensure linker doesn't optimise out the class
bool TurboSwapLoad() {return TurboSwap::TYPE != 0;}

DRLIB_END_NAMESPACE
