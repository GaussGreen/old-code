//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericNFactorTester.cpp
//
//   Description : N-factor that tests EAs interface
//
//   Author      : Andrew J Swain
//
//   Date        : 14 January 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericNFactor.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/ATMVolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

class GenericNFactorTester: public GenericNFactor, 
                            virtual public CClosedFormLN::IIntoProduct, 
                            virtual public IMCIntoProduct
{
protected:
    /// fields ////////
    CashFlowArray coupon;
    StringArray   perfType;
    DoubleArray   weight;       
    bool          isCall;
    double        overallStrike;
    StringArray   flags;

public:
    static CClassConstSP const TYPE;
    friend class GNFClosedForm;

    // validation
    void validatePop2Object(){
        // check that we've got as many weights as there are assets 
        AssetUtil::checkWeights(weight, assets->NbAssets());
    }
   
    /** Implementation of ClosedFormLN::IntoProduct interface - the
        implementation of this is below */
    virtual CClosedFormLN::IProduct* createProduct(
        CClosedFormLN* model) const;

    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

    void price(CResults* results) const {
        static const string method = "GenericNFactorTester::price";
        try {
            double value = 0.0;

            for (int i = 0; i < assets->NbAssets(); i++) {
                IMarketFactorConstSP factor(assets->getFactor(i));
                // to do: enhance this test
                if (!IGeneralAsset::TYPE->isInstance(factor)){
                    throw ModelException(method, "Only IGeneralAssets "
                                         "supported");
                }
                const IGeneralAsset* asset = DYNAMIC_CAST(IGeneralAsset,
                                                          factor.get());
                value += weight[i]*coupon[i].amount/asset->getSpot();
            }

            if (isCall) {
                value = Maths::max(value - overallStrike, 0.0);
            }
            else {
                value = Maths::max(overallStrike - value, 0.0);
            }
            value *= notional;
            results->storePrice(value, discount->getCcy());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }


private:
    friend class GNFMC;
    GenericNFactorTester(): GenericNFactor(TYPE), isCall(true) {} // for reflection
    GenericNFactorTester(const GenericNFactorTester& rhs); // not implemented
    GenericNFactorTester& operator=(const GenericNFactorTester& rhs); // not implemented

    static IObject* defaultGenericNFactorTester(){
        return new GenericNFactorTester();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(GenericNFactorTester, clazz);
        SUPERCLASS(GenericNFactor);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultGenericNFactorTester);
        FIELD(coupon,           "coupons");
        FIELD(perfType,         "perfType");
        FIELD(weight,           "weights");
        FIELD(isCall,           "is it a call option");
        FIELD_MAKE_OPTIONAL(isCall);
        FIELD(overallStrike,    "overallStrike");
        FIELD(flags,            "flags");
        FIELD_MAKE_OPTIONAL(flags);
    }
};

/** private class */
class GNFClosedForm: public CClosedFormLN::IProduct{
private:
    const GenericNFactorTester* gnf; // a reference

public:
    GNFClosedForm(const GenericNFactorTester* gnf): gnf(gnf){}

    void price(CClosedFormLN* model,
               Control*       control, 
               CResults*      results) const{
        gnf->price(results);
    }
};


CClosedFormLN::IProduct* GenericNFactorTester::createProduct(
    CClosedFormLN* model) const {
    return new GNFClosedForm(this);
}

class GNFMC: public IMCProduct, virtual public IMCProductLN {
private:
    const GenericNFactorTester* gnf; // a reference

public:
    GNFMC(const GenericNFactorTester* gnf, 
          const SimSeriesConstSP&     simSeries,
          const IRefLevelConstSP&     refLevel,
          const IPastValuesConstSP&   mcPastValues):
        IMCProduct(gnf->assets.get(),
                  gnf->valueDate,
                  gnf->discount.get(),
                  refLevel,
                  simSeries,
                  mcPastValues,
                  gnf->instSettle.get(),
                  gnf->coupon[0].date),
         gnf(gnf){}

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    void payoff(const IPathGenerator *pathGen, 
                IMCPrices &              prices) {
        double value = 0.0;

        for (int i = 0; i < gnf->assets->NbAssets(); i++) {
            value += gnf->weight[i]*gnf->coupon[i].amount/
                getMultiFactors()->assetGetSpot(i);
        }

        if (gnf->isCall) {
            value = Maths::max(value - gnf->overallStrike, 0.0);
        }
        else {
            value = Maths::max(gnf->overallStrike - value, 0.0);
        }
        value *= gnf->notional;

        prices.add(value);       
    }

    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray reqarr(1); // one interp level/path per asset here
        reqarr[0] = CVolRequestLNSP(new ATMVolRequest());
        return reqarr;
    }
};

IMCProduct* GenericNFactorTester::createProduct(const MonteCarlo* model) const {
    SimSeriesSP        simSeries(new SimSeries(assets->NbAssets()));
    DateTimeArray      dates(1);
    dates[0] = coupon[0].date;
    simSeries->addDates(dates);

    IRefLevelConstSP   refLevel(IRefLevel::Util::makeZero(valueDate));
    DoubleArray past(assets->NbAssets());
    for (int i = 0; i < past.size(); i++) {
        IMarketFactorConstSP factor(assets->getFactor(i));
        // to do: enhance this test
        if (!IGeneralAsset::TYPE->isInstance(factor)){
            throw ModelException("GenericNFactorTester::createProduct",
                                 "Only IGeneralAssets supported");
        }
        const IGeneralAsset* asset = DYNAMIC_CAST(IGeneralAsset,
                                                  factor.get());
        past[i] = asset->getSpot();
    }
    
    IPastValuesConstSP mcPastValues(IPastValues::Util::makeTrivial(valueDate, past));
    return new GNFMC(this, simSeries, refLevel, mcPastValues);
}

CClassConstSP const GenericNFactorTester::TYPE = CClass::registerClassLoadMethod(
    "GenericNFactorTester", typeid(GenericNFactorTester), GenericNFactorTester::load);

// * for class loading (avoid having header file) */
bool GenericNFactorTesterLoad() {
    return (GenericNFactorTester::TYPE != 0);
}

DRLIB_END_NAMESPACE

