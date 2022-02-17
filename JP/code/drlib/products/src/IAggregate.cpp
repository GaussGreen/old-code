//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IAggregate.cpp
//
//   Description : Captures collapse of dimension from array to single value in various ways
//
//   Date        : Nov 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_IAGGREGATE_CPP
#include "edginc/IAggregate.hpp"
#include "edginc/Algorithm.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

void IAggregateMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IAggregateMaker, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IAggregateMaker::TYPE = CClass::registerInterfaceLoadMethod(
    "IAggregateMaker", typeid(IAggregateMaker), load);

/*********************************************************************/
// Standard basket, but assumes components are supplied as perfs (for Pct)
class BasketAggregate : virtual public IAggregate {
public:
    BasketAggregate(BasketAggregateMaker* maker,
                    IDoubleArray*         components):
         maker(BasketAggregateMakerSP(copy(maker))), components(components) {
        if (maker->weights.size() != components->size()) {
            throw ModelException("BasketAggregate: Number of weights (" + Format::toString(maker->weights.size()) +
                ") must equal number of components (" + Format::toString(components->size()) + ")");
        }
    };

    virtual double aggregate() {
        double val = 0.0;
        for(int i=0; i<maker->weights.size(); i++) {
            val += maker->weights[i] * (*components)[i];
        }
        return val;
    };

private:
    BasketAggregateMakerSP maker;
    IDoubleArray*          components;
};

class BasketEqualAggregate : virtual public IAggregate {
public:
    BasketEqualAggregate(IDoubleArray* components):
         components(components) {
    };

    virtual double aggregate() {
        double val = 0.0;
        for(int i=0; i<components->size(); i++) {
            val += (*components)[i];
        }
        return val/components->size();
    };

private:
    IDoubleArray*          components;
};

class BasketEqualAggregateFiltered : virtual public IAggregate {
public:
    BasketEqualAggregateFiltered(IDoubleArray*  components,
                                 IAssetFilterSP filter):
        components(components),
        activeComponents(filter->getActiveAssets(false)) {}; // false -> not doing past ... not sure how this will play out

    virtual double aggregate() {
        double val = 0.0;
        for(int i=0; i<activeComponents.size(); i++) {
            int iComp = activeComponents[i];
            val += (*components)[iComp];
        }
        return val/activeComponents.size();
    };

private:
    IDoubleArray*   components;
    IAssetFilterSP  filter;
    const IntArray& activeComponents;
};


// The published face
BasketAggregateMaker::BasketAggregateMaker(const DoubleArray& weights):
    CObject(TYPE), weights(weights) {
    validatePop2Object();
}

void BasketAggregateMaker::validatePop2Object(){
// Dont make this check. In generalised world no need to be so fussy
//    if (weights.size()>0) {
//        Check::percWeights(weights, weights.size(), "");
//    }
}

IAggregate* BasketAggregateMaker::getAggregate(IDoubleArray* comps) {
    if (weights.size()>0) {
        return new BasketAggregate(this, comps);
    }
    return new BasketEqualAggregate(comps);
}

IAggregate* BasketAggregateMaker::getAggregate(IDoubleArray*   comps,
                                               IAssetFilterSP  filter) {
    // By identifying no filtering we avoid some restrictions
    // I know we could check for numToDrop=0 etc but that is invasion of 
    // the Filter's privacy...
    if (dynamic_cast<NoFilter*>(filter.get())) {
        return getAggregate(comps);
    }

    // we don't know how many components we'll have
    bool equalWeights = true;
    for (int i=1; i<weights.size(); i++) {
        if (!Maths::equals(weights[0], weights[i])){
            equalWeights = false;
            break;
        }
    }
    if (!equalWeights) {
        throw ModelException("BasketAggregateMaker::getAggregate",
                             "Only support Drop with equal weights");
    }
    return new BasketEqualAggregateFiltered(comps, filter);
}

class BasketAggregateMakerHelper{
public:
    static IObject* defaultBasketAggregateMaker(){
        return new BasketAggregateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BasketAggregateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAggregateMaker);
        EMPTY_SHELL_METHOD(defaultBasketAggregateMaker);
        FIELD(weights,          "weights");
        FIELD_MAKE_OPTIONAL(weights);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const BasketAggregateMaker::TYPE =
CClass::registerClassLoadMethod("BasketAggregateMaker", 
                                typeid(BasketAggregateMaker), BasketAggregateMakerHelper::load);

/*********************************************************************/
// Rainbowed version of above
class RainbowAggregate : virtual public IAggregate {
public:
    class OrderDoubles {
    public:
        int operator() (double d1, double d2) {
            return d1 > d2;
        }
    };

    RainbowAggregate(RainbowAggregateMaker* maker,
                     IDoubleArray*          components):
         maker(RainbowAggregateMakerSP(copy(maker))), components(components) {
        if (maker->weights.size() != components->size()) {
            throw ModelException("RainbowAggregate: Number of weights (" + Format::toString(maker->weights.size()) +
                ") must equal number of components (" + Format::toString(components->size()) + ")" );
        }
    };

    // XXX modifies components - may wish to have them safe and modify a local copy??
    virtual double aggregate() {
        Algorithm::shellSort(*components, components->begin(), components->end()-1);
        double val = 0.0;
        for(int i=0; i<maker->weights.size(); i++) {
            val += maker->weights[i] * (*components)[i];
        }
        return val;
    };

private:
    RainbowAggregateMakerSP maker;
    IDoubleArray*           components;
};

class RainbowAggregateFiltered : virtual public IAggregate {
public:
    RainbowAggregateFiltered(IDoubleArray*  components,
                             IAssetFilterSP filter,
                             bool           isBest):
        components(components),
        activeComponents(filter->getActiveAssets(false)),
        isBest(isBest) {}; // false -> not doing past ... not sure how this will play out

    virtual double aggregate() {
        if (activeComponents.size()<1) {
            throw ModelException("RainbowAggregateFiltered::aggregate",
                                 "Nothing to aggregate!");
        }
        double val = (*components)[activeComponents[0]];
        for(int i=1; i<activeComponents.size(); i++) {
            int iComp = activeComponents[i];
            if (isBest) {
                val = Maths::max(val, (*components)[iComp]);
            } else {
                val = Maths::min(val, (*components)[iComp]);
            }
        }
        return val;
    };

private:
    IDoubleArray*   components;
    IAssetFilterSP  filter;
    const IntArray& activeComponents;
    bool            isBest;
};


RainbowAggregateMaker::RainbowAggregateMaker(const DoubleArray& weights):
    CObject(TYPE), weights(weights) {
}

class RainbowAggregateMakerHelper {
public:
    static IObject* defaultRainbowAggregateMaker(){
        return new RainbowAggregateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RainbowAggregateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAggregateMaker);
        EMPTY_SHELL_METHOD(defaultRainbowAggregateMaker);
        FIELD(weights,          "weights");
        FIELD_MAKE_OPTIONAL(weights);     // for IMS only!
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

IAggregate* RainbowAggregateMaker::getAggregate(IDoubleArray* comps) {
    if (weights.size()<1) {
        throw ModelException("RainbowAggregateMaker::getAggregate",
                             "No weights supplied!");
    }
    return new RainbowAggregate(this, comps);
}


IAggregate* RainbowAggregateMaker::getAggregate(IDoubleArray*   comps,
                                                IAssetFilterSP  filter) {
    static const string method = "RainbowAggregateMaker::getAggregate";

    if (dynamic_cast<NoFilter*>(filter.get())) {
        return getAggregate(comps);
    }

    if (weights.size()<1) {
        throw ModelException(method,
                             "No weights supplied!");
    }
    // we don't know how many components we'll have - but can suppose at 
    // least one so offer best or worst. i.e. all weights 0 except one
    // which is best or worst position
    bool isBestOrWorst = true;
    int i;
    for(i=1; i<weights.size()-1; i++) {
        if (!Maths::isZero(weights[i])) {
            isBestOrWorst = false;
        }
    }
    if (!isBestOrWorst ||
        (!Maths::isZero(weights.front()) &&
         !Maths::isZero(weights.back()))) {
        throw ModelException(method,
                             "Only best or worst supported with drop");
    }
    // Maths::isZero(weights.back()) => isBest
    return new RainbowAggregateFiltered(comps, filter, 
                                        Maths::isZero(weights.back()));
}

CClassConstSP const RainbowAggregateMaker::TYPE =
CClass::registerClassLoadMethod("RainbowAggregateMaker", 
                                typeid(RainbowAggregateMaker), RainbowAggregateMakerHelper::load);

/*********************************************************************/
// Product (i.e. multiply) version
class ProductAggregate : virtual public IAggregate {
public:
    ProductAggregate(ProductAggregateMaker* maker,
                     IDoubleArray*          components):
         maker(ProductAggregateMakerSP(copy(maker))), components(components) {
        if (components->size()<1) {
            throw ModelException("ProductAggregate: Must have at least 1 component! Given : " + 
                                 Format::toString(components->size()));
        }
    };

    virtual double aggregate() {
        double val = 1.0;
        for(int i=0; i<components->size(); i++) {
            val *= (*components)[i];
        }
        return val;
    };

private:
    ProductAggregateMakerSP maker;
    IDoubleArray*           components;
};

class ProductAggregateFiltered : virtual public IAggregate {
public:
    ProductAggregateFiltered(IDoubleArray*  components,
                             IAssetFilterSP filter):
        components(components),
        activeComponents(filter->getActiveAssets(false)) {}; // false -> not doing past ... not sure how this will play out

    virtual double aggregate() {
        if (activeComponents.size()<1) {
            throw ModelException("ProductAggregateFiltered::aggregate",
                                 "Nothing to aggregate!");
        }
        double val = 1.0;
        for(int i=0; i<activeComponents.size(); i++) {
            int iComp = activeComponents[i];
            val *= (*components)[iComp];
        }
        return val;
    };

private:
    IDoubleArray*   components;
    IAssetFilterSP  filter;
    const IntArray& activeComponents;
};


ProductAggregateMaker::ProductAggregateMaker(): CObject(TYPE) {};

class ProductAggregateMakerHelper {
public:
    static IObject* defaultProductAggregateMaker(){
        return new ProductAggregateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ProductAggregateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAggregateMaker);
        EMPTY_SHELL_METHOD(defaultProductAggregateMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

IAggregate* ProductAggregateMaker::getAggregate(IDoubleArray* comps) {
    return new ProductAggregate(this, comps);
}

IAggregate* ProductAggregateMaker::getAggregate(IDoubleArray*   comps,
                                                IAssetFilterSP  filter) {
    if (dynamic_cast<NoFilter*>(filter.get())) {
        return getAggregate(comps);
    }
    return new ProductAggregateFiltered(comps, filter);
}

CClassConstSP const ProductAggregateMaker::TYPE =
CClass::registerClassLoadMethod("ProductAggregateMaker", 
                                typeid(ProductAggregateMaker), ProductAggregateMakerHelper::load);


/*********************************************************************/
// Reinvested version
class ReinvestAggregate : virtual public IAggregate {
public:
    ReinvestAggregate(ReinvestAggregateMaker* maker,
                     IDoubleArray*          components):
         maker(ReinvestAggregateMakerSP(copy(maker))), components(components) {
        if (maker->weights.size() != components->size()) {
            throw ModelException("ReinvestAggregate: Number of weights (" + Format::toString(maker->weights.size()) +
                ") must equal number of components (" + Format::toString(components->size()) + ")" );
        }
    };

    virtual double aggregate() {
        double val = 1.0;
        for(int i=0; i<maker->weights.size(); i++) {
            val *= 1.0 + maker->weights[i] * (*components)[i];
        }
        val -= 1.0;
        return val;
    };

private:
    ReinvestAggregateMakerSP maker;
    IDoubleArray*            components;
};

ReinvestAggregateMaker::ReinvestAggregateMaker(const DoubleArray& weights):
    CObject(TYPE), weights(weights) {
};

class ReinvestAggregateMakerHelper {
public:
    static IObject* defaultReinvestAggregateMaker(){
        return new ReinvestAggregateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ReinvestAggregateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAggregateMaker);
        EMPTY_SHELL_METHOD(defaultReinvestAggregateMaker);
        FIELD(weights,          "weights");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

IAggregate* ReinvestAggregateMaker::getAggregate(IDoubleArray* comps) {
    if (weights.size()<1) {
        throw ModelException("ReinvestAggregateMaker::getAggregate",
                             "No weights supplied!");
    }
    return new ReinvestAggregate(this, comps);
}

IAggregate* ReinvestAggregateMaker::getAggregate(IDoubleArray*   comps,
                                                 IAssetFilterSP  filter) {
    static const string method = "ReinvestAggregateMaker::getAggregate";
    if (!dynamic_cast<NoFilter*>(filter.get())) {
        throw ModelException("ReinvestAggregateMaker::getAggregate",
                             "Drop not supported");
    }
    return getAggregate(comps);
}

CClassConstSP const ReinvestAggregateMaker::TYPE =
CClass::registerClassLoadMethod("ReinvestAggregateMaker", 
                                typeid(ReinvestAggregateMaker), ReinvestAggregateMakerHelper::load);

/*********************************************************************/
// ReinvestRainbowed version
class ReinvestRainbowAggregate : virtual public IAggregate {
public:
    class OrderDoubles {
    public:
        int operator() (double d1, double d2) {
            return d1 > d2;
        }
    };

    ReinvestRainbowAggregate(ReinvestRainbowAggregateMaker* maker,
                     IDoubleArray*          components):
         maker(ReinvestRainbowAggregateMakerSP(copy(maker))), components(components) {
        if (maker->weights.size() != components->size()) {
            throw ModelException("ReinvestRainbowAggregate: Number of weights (" + Format::toString(maker->weights.size()) +
                ") must equal number of components (" + Format::toString(components->size()) + ")" );
        }
    };

    virtual double aggregate() {
        Algorithm::shellSort(*components, components->begin(), components->end()-1);
        double val = 1.0;
        for(int i=0; i<maker->weights.size(); i++) {
            val *= 1.0 + maker->weights[i] * (*components)[i];
        }
        val -= 1.0;
        return val;
    };

private:
    ReinvestRainbowAggregateMakerSP maker;
    IDoubleArray*                   components;
};


ReinvestRainbowAggregateMaker::ReinvestRainbowAggregateMaker(const DoubleArray& weights):
    CObject(TYPE), weights(weights) {
};

class ReinvestRainbowAggregateMakerHelper {
public:
    static IObject* defaultReinvestRainbowAggregateMaker(){
        return new ReinvestRainbowAggregateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ReinvestRainbowAggregateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAggregateMaker);
        EMPTY_SHELL_METHOD(defaultReinvestRainbowAggregateMaker);
        FIELD(weights,          "weights");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

IAggregate* ReinvestRainbowAggregateMaker::getAggregate(IDoubleArray* comps) {
    if (weights.size()<1) {
        throw ModelException("ReinvestRainbowAggregateMaker::getAggregate",
                             "No weights supplied!");
    }
    return new ReinvestRainbowAggregate(this, comps);
}

IAggregate* ReinvestRainbowAggregateMaker::getAggregate(IDoubleArray*   comps,
                                                        IAssetFilterSP  filter) {
    static const string method = "ReinvestRainbowAggregateMaker::getAggregate";

    if (!dynamic_cast<NoFilter*>(filter.get())) {
        throw ModelException("ReinvestRainbowAggregateMaker::getAggregate",
                             "Drop not supported");
    }
    return getAggregate(comps);
}


CClassConstSP const ReinvestRainbowAggregateMaker::TYPE =
CClass::registerClassLoadMethod("ReinvestRainbowAggregateMaker", 
                                typeid(ReinvestRainbowAggregateMaker), ReinvestRainbowAggregateMakerHelper::load);

/*********************************************************************/
// Conditional version - returns lowest/highest above/below some level
// if there are no such elements it returns the level itself
class ConditionalAggregate : virtual public IAggregate {
public:
    ConditionalAggregate(ConditionalAggregateMaker* maker,
                     IDoubleArray*          components):
         maker(ConditionalAggregateMakerSP(copy(maker))), components(components) {}

    virtual double aggregate() {
        Algorithm::shellSort(*components, components->begin(), components->end()-1);
        double val = maker->level;
        if (maker->isAbove) {
            for(int i=components->size() - 1; i>=0; i--) {
                if ((*components)[i] > maker->level) {
                    val = (*components)[i];
                    break;
                }
            }
        } else {
            for(int i=0; i<components->size(); i++) {
                if ((*components)[i] < maker->level) {
                    val = (*components)[i];
                    break;
                }
            }
        }
        return val;
    };

private:
    ConditionalAggregateMakerSP maker;
    IDoubleArray*                   components;
};

ConditionalAggregateMaker::ConditionalAggregateMaker():
    CObject(TYPE), level(0.0), isAbove(false) {
};

ConditionalAggregateMaker::ConditionalAggregateMaker(double level, bool isAbove):
    CObject(TYPE), level(level), isAbove(isAbove) {
};

class ConditionalAggregateMakerHelper {
public:
    static IObject* defaultConditionalAggregateMaker(){
        return new ConditionalAggregateMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ConditionalAggregateMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAggregateMaker);
        EMPTY_SHELL_METHOD(defaultConditionalAggregateMaker);
        FIELD(level,          "cut off level");
        FIELD_MAKE_OPTIONAL(level);
        FIELD(isAbove,          "Value above level?");
        FIELD_MAKE_OPTIONAL(isAbove);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

IAggregate* ConditionalAggregateMaker::getAggregate(IDoubleArray* comps) {
    return new ConditionalAggregate(this, comps);
}

IAggregate* ConditionalAggregateMaker::getAggregate(IDoubleArray*   comps,
                                                        IAssetFilterSP  filter) {
    static const string method = "ConditionalAggregateMaker::getAggregate";

    if (!dynamic_cast<NoFilter*>(filter.get())) {
        throw ModelException("ConditionalAggregateMaker::getAggregate",
                             "Drop not supported");
    }
    return getAggregate(comps);
}


CClassConstSP const ConditionalAggregateMaker::TYPE =
CClass::registerClassLoadMethod("ConditionalAggregateMaker", 
                                typeid(ConditionalAggregateMaker), ConditionalAggregateMakerHelper::load);

/*********************************************************************/

/* The ultimate wrapping of IDoubleArrayModifierMaker's mainly for use in Pyramid 
 */
#define AGGREGATE_TYPE_PLAINBASKET   "PlainBasket"
#define AGGREGATE_TYPE_RAINBOW       "RainbowBasket"
#define AGGREGATE_TYPE_PRODUCT       "ProductBasket"
#define AGGREGATE_TYPE_REINVEST      "ReinvestBasket"
#define AGGREGATE_TYPE_REINVEST_RBW  "ReinvestRainbowBasket"
#define AGGREGATE_TYPE_CONDITIONAL   "ConditionalBasket"

class AggregateMakerWrapper : public CObject,
                              virtual public IAggregateMaker {
public: // how can I have this protected or private?
    string                              aggregateMakerType;     // Plain or Rainbow
    BasketAggregateMakerSP              basketMaker;
    RainbowAggregateMakerSP             rainbowMaker;
    ProductAggregateMakerSP             productMaker;
    ReinvestAggregateMakerSP            reinvestMaker;
    ReinvestRainbowAggregateMakerSP     reinvestRainbowMaker;
    ConditionalAggregateMakerSP         conditionalMaker;

private:
    IAggregateMakerSP  realMaker;

public:
    static CClassConstSP const TYPE;

    virtual IAggregate* getAggregate(IDoubleArray* comps){
        return realMaker->getAggregate(comps);
    }

    virtual IAggregate* getAggregate(IDoubleArray*   comps,
                                     IAssetFilterSP  filter) {
        return realMaker->getAggregate(comps, filter);
    }

    // validation
    void validatePop2Object(){
        static const string routine = "AggregateMakerWrapper::validatePop2Object";

        if (aggregateMakerType.empty()){
            throw ModelException(routine, "Blank Aggregate Maker specified!");
        }
        if (aggregateMakerType==AGGREGATE_TYPE_PLAINBASKET) {
            if (basketMaker.get()) {
                realMaker = basketMaker;
            } else {
                throw ModelException(routine, "Expected plain basketMaker but none supplied!");
            }
        } else if (aggregateMakerType==AGGREGATE_TYPE_RAINBOW) {
            if (rainbowMaker.get()) {
                realMaker = rainbowMaker;
            } else {
                throw ModelException(routine, "Expected rainbowMaker but none supplied!");
            }
        } else if (aggregateMakerType==AGGREGATE_TYPE_PRODUCT) {
            if (productMaker.get()) {
                realMaker = productMaker;
            } else {
                throw ModelException(routine, "Expected productMaker but none supplied!");
            }
        } else if (aggregateMakerType==AGGREGATE_TYPE_REINVEST) {
            if (reinvestMaker.get()) {
                realMaker = reinvestMaker;
            } else {
                throw ModelException(routine, "Expected reinvestMaker but none supplied!");
            }
        } else if (aggregateMakerType==AGGREGATE_TYPE_REINVEST_RBW) {
            if (reinvestRainbowMaker.get()) {
                realMaker = reinvestRainbowMaker;
            } else {
                throw ModelException(routine, "Expected reinvestRainbowMaker but none supplied!");
            }
        } else if (aggregateMakerType==AGGREGATE_TYPE_CONDITIONAL   ) {
            if (conditionalMaker.get()) {
                realMaker = conditionalMaker;
            } else {
                throw ModelException(routine, "Expected conditionalMaker but none supplied!");
            }
        } else {
            throw ModelException(routine, "Unrecognised Aggregate Maker " + aggregateMakerType + 
                                 ". Expected " + AGGREGATE_TYPE_PLAINBASKET + " or " + AGGREGATE_TYPE_RAINBOW +
                                 " or " AGGREGATE_TYPE_PRODUCT + " or " + AGGREGATE_TYPE_REINVEST +
                                 " or " AGGREGATE_TYPE_REINVEST_RBW + " or " + AGGREGATE_TYPE_CONDITIONAL);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AggregateMakerWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAggregateMaker);
        EMPTY_SHELL_METHOD(defaultAggregateMakerWrapper);
        FIELD(aggregateMakerType, "PlainBasket or RainbowBasket");
        FIELD(basketMaker,  "basketMaker");
        FIELD_MAKE_OPTIONAL(basketMaker);
        FIELD(rainbowMaker,  "rainbowMaker");
        FIELD_MAKE_OPTIONAL(rainbowMaker);
        FIELD(productMaker,  "productMaker");
        FIELD_MAKE_OPTIONAL(productMaker);
        FIELD(reinvestMaker,  "reinvestMaker");
        FIELD_MAKE_OPTIONAL(reinvestMaker);
        FIELD(reinvestRainbowMaker,  "reinvestRainbowMaker");
        FIELD_MAKE_OPTIONAL(reinvestRainbowMaker);
        FIELD(conditionalMaker,  "conditionalMaker");
        FIELD_MAKE_OPTIONAL(conditionalMaker);
        FIELD(realMaker, "realMaker");
        FIELD_MAKE_TRANSIENT(realMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    AggregateMakerWrapper(): CObject(TYPE){}

    static IObject* defaultAggregateMakerWrapper(){
        return new AggregateMakerWrapper();
    }
};

typedef smartPtr<AggregateMakerWrapper> AggregateMakerWrapperSP;

CClassConstSP const AggregateMakerWrapper::TYPE =
CClass::registerClassLoadMethod("AggregateMakerWrapper", 
                                typeid(AggregateMakerWrapper), load);

DRLIB_END_NAMESPACE

    

