//----------------------------------------------------------------------------
//
//   Group       : QR&D Interest Rates
//
//   Filename    : Tarn.cpp
//
//   Description : instrument shell that makes use of components 
//                 with customized internal logic to fill optional data fields
//                 for each components. 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KAccumulatedCF.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

// the shell instrument for Tarn
// It looks like we could factor out a new base class for "shell instruments"
// that are based on the KComponent architecture and move common functionality
// into the new base.
class Tarn : public CInstrument,
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
    Tarn() : CInstrument(TYPE) {}
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new Tarn(); }

    /****************** exported fields ************/
    KAccumulatedCFSP    tarn;
    YieldCurveWrapper   discount;

    /****************** transiend fields ************/
    KComponentSP        compRoot;
};

// export the members of the class through the library interface
void Tarn::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(Tarn, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(tarn, "tarn component");
    FIELD(discount,"discount curve");
    FIELD(compRoot,""); FIELD_MAKE_TRANSIENT(compRoot);
    Addin::registerConstructor(Addin::UTILITIES, Tarn::TYPE);
}

// We link the Tarn components up in GetMarket
void Tarn::GetMarket(const IModel* model, const CMarketDataSP market) {
    try{

        if (tarn->outputName.empty())
            tarn->outputName = "TARN_LEG_VALUE";
        if (tarn->discountYieldCurveName().empty()) {
            tarn->discount=discount;
        }

        compRoot = tarn;

        if (compRoot->discountYieldCurveName().empty())
            compRoot->discount = discount;

        // cleanup fields so that iterator do not get lost
        tarn.reset(0);

        compRoot->GetMarket(model, market);
    }
    catch (exception& e) {
        throw ModelException(e, "Tarn::GetMarket");
    }
}

// register the Tarn in the library framework at startup
CClassConstSP const Tarn::TYPE =  CClass::registerClassLoadMethod(
            "Tarn", typeid(Tarn), Tarn::load);

// to ensure linker doesn't optimise out the class
bool TarnLoad(void) {
    return (Tarn::TYPE != 0);
}

DRLIB_END_NAMESPACE
