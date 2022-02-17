//----------------------------------------------------------------------------
//
//   Group       : Energy Exotics Derivatives Research
//
//   Filename    : BasketAverage.cpp
//
//   Description : New Basket Average shell instrument 
//                 
//   Date        : November 2006
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Format.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/ITaxableInst.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/MCProductIQuicks.hpp"
#include "edginc/DateBuilder.hpp"
#include "edginc/ObservationBuilder.hpp"
#include "edginc/ContractsOffset.hpp"


DRLIB_BEGIN_NAMESPACE

/** BasketAverage product - option based upon performance of n assets
    sorted by their performance **/
class BasketAverage: public GenericNFBase,
                     public ITaxableInst::Basic{
protected:
    /// fields ////////
    string                          weightType;
    DoubleArrayArray                weights;
    bool                            isCall;
    double                          strike;
    IDateBuilderSP                  fixingDates;
    bool                            isRefAvgOut;
    ContractsOffsetArrayArray       offsets;

public:
    static CClassConstSP const TYPE;
    //friend class BasketAverageSVMC;

    // validation
    void validatePop2Object(){
        static const string routine("BasketAverage::validatePop2Object");
        GenericNFBase::validatePop2Object();
    
        // check either num of weights equals num of fixings or single weight that applied 
        // to all fixings for that asset
        for(int i=0;assets->NbAssets();++i) {
            if(weights[i].size()!=fixingDates->size() || weights[i].size()!=1) {
                throw ModelException(routine,"Incorrect number of weights for asset");
            }
        }
        // check num of assets equals num of arrays of weights
        if(assets->NbAssets() != weights.size()) {
            throw ModelException(routine, "Different number of assets ("+
                                          Format::toString(assets->NbAssets())+")"
                                          " to array of weights ("+
                                          Format::toString(weights.size())+")");
        }
        // check weight type either Percentage(P) or Unit(U)
        if (weightType != "P" || weightType != "U") {
            throw ModelException(routine, "BasketAverage: weightType must be U or P,"
                                 " but " + weightType + " given");
        }
        // check for negative strike
        if (Maths::isNegative(strike)){
            throw ModelException(routine, "strike ("+
                Format::toString(strike)+") is negative");
        }
        // check for empty fixings array
        if (fixingDates->size() == 0) {
            throw ModelException(routine, "No fixings supplied!");
        }
        if (isRefAvgOut && !isCall) {
            throw ModelException(routine, "Combination of Put with isRefAvgOut True is forbidden");
        }
        // check that num of assets matches num of offsets cols
        if (assets->NbAssets() != offsets.size()) {
            throw ModelException(routine, "Different number of assets ("+
                                          Format::toString(assets->NbAssets())+")"
                                          " to offsets ("+
                                          Format::toString(offsets.size())+")");
        }
        // check either num of offsets equals num of fixings or single offset is applied 
        //to all fixings for that asset
        for(int i=0;assets->NbAssets();++i) {
            if(offsets[i]->size()!=fixingDates->size() || offsets[i]->size()!=1) {
                throw ModelException(routine,"Incorrect number of offsets for asset");
            }
        }
    }

    virtual void GetMarket(const IModel* model, const CMarketDataSP market){
        static const string method = "BasketAverage::GetMarket";
        try {
            GenericNFBase::GetMarket(model, market);

            // get the market for the IDateBuilder
            fixingDates->getMarket(model, market.get());

        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    void Validate() {
        static const string method = "BasketAverage::Validate";
        try {
            GenericNFBase::Validate();

            // any extra validation to be put here

        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    const DateTimeArray samplingDates() const {
        return *(fixingDates->dates().get());
    }

    const DateTime getFinalPaymentDate() const {
        return instSettle->settles(fixingDates->end(), NULL);
    }

private:
    BasketAverage(): GenericNFBase(TYPE), isRefAvgOut(false) {}     // reflection
    BasketAverage(const BasketAverage& rhs);        // not implemented  
    BasketAverage& operator=(const BasketAverage& rhs);     // not implemented

    static IObject* defaultBasketAverage() {
        return new BasketAverage();
    }

    // Invoked then Class is 'loaded'
    static void load(CClassSP& clazz){
        REGISTER(BasketAverage, clazz);
        SUPERCLASS(GenericNFBase);
        EMPTY_SHELL_METHOD(defaultBasketAverage);
        FIELD(weightType, "Unit(U) or Percentage(P) basket");
        FIELD(weights, "Weights");
        FIELD(isCall, "is it a call option");
        FIELD(strike, "strike");
        FIELD(fixingDates, "fixingDates");
        FIELD(isRefAvgOut, "isRefAvgOut");
        FIELD_MAKE_OPTIONAL(isRefAvgOut);
        FIELD(offsets, "an array of offset arrays");
        clazz->setPublic();
    }
};

CClassConstSP const BasketAverage::TYPE = CClass::registerClassLoadMethod(
    "BasketAverage", typeid(BasketAverage), BasketAverage::load);

// * for class loading (avoid having header file) */
bool BasketAverageLoad() {
    return (BasketAverage::TYPE != 0);
}

DRLIB_END_NAMESPACE