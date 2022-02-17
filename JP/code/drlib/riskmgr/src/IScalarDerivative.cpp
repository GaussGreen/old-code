/**
 * @file IScalarDerivative.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/Maths.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/AtomicHypothesis.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/HypotheticalQuantity.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"

DRLIB_BEGIN_NAMESPACE

IScalarDerivative::IScalarDerivative() {}
IScalarDerivative::~IScalarDerivative() {}

// 
// **********
//  oneSided
// **********
// 

class OneSidedFirstDerivative: public CObject,
                               public virtual IScalarDerivative {

    static void load(CClassSP& clazz) {
        REGISTER(OneSidedFirstDerivative, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScalarDerivative);
        EMPTY_SHELL_METHOD(DefaultConstructor<OneSidedFirstDerivative>::iObject);
    }

public:

    struct RQ: public RiskQuantity {

        static CClassConstSP const TYPE;

        static void load(CClassSP& clazz) {
            REGISTER(RQ, clazz);
            SUPERCLASS(RiskQuantity);
            EMPTY_SHELL_METHOD(DefaultConstructor<RQ>::iObject);
        }

        RQ(HypotheticalQuantityArrayConstSP parameters):
            RiskQuantity(TYPE, parameters)
        {}

        RQ(): RiskQuantity(TYPE) {}

        double value(const DoubleArray& vals,
                     const DoubleArray& dists) const {
            TRACE_METHOD;
            if (Maths::isZero(dists[1])) {
                TRACE("Oops, distance is zero");
                throw ModelException("Zero divisor in one-sided first deriv");
            }
            TRACE_SHOW((vals[1] - vals[0]) / dists[1]);
            return (vals[1] - vals[0]) / dists[1];
        }
    };

    static CClassConstSP const TYPE;

    OneSidedFirstDerivative(): CObject(TYPE) {}

    RiskQuantitySP riskQuantity(IResultsFunctionConstSP function,
                                IRiskAxisConstSP axis,
                                double coeff) const {
        TRACE_METHOD;

        HypotheticalQuantityArraySP ps(new HypotheticalQuantityArray(2));

        (*ps)[0] = HypotheticalQuantity::SP(axis->hypothesis(0), function);
        (*ps)[1] = HypotheticalQuantity::SP(axis->hypothesis(coeff), function);

        return RiskQuantitySP(new RQ(ps));
    }

    int order() const { return 1; }
};

CClassConstSP const OneSidedFirstDerivative::TYPE = CClass::registerClassLoadMethod(
    "OneSidedFirstDerivative", typeid(OneSidedFirstDerivative), load);

typedef OneSidedFirstDerivative::RQ OneSidedFirstDerivative_RQ;
CClassConstSP const OneSidedFirstDerivative_RQ::TYPE = CClass::registerClassLoadMethod(
    "OneSidedFirstDerivative::RQ", typeid(OneSidedFirstDerivative_RQ), OneSidedFirstDerivative_RQ::load);

IScalarDerivativeConstSP IScalarDerivative::oneSided() {
    static IScalarDerivativeConstSP it(new OneSidedFirstDerivative());
    return it;
}

// 
// **********
//  twoSided
// **********
// 

class TwoSidedFirstDerivative: public CObject,
                               public virtual IScalarDerivative {

    static void load(CClassSP& clazz) {
        REGISTER(TwoSidedFirstDerivative, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScalarDerivative);
        EMPTY_SHELL_METHOD(DefaultConstructor<TwoSidedFirstDerivative>::iObject);
    }

public:

    struct RQ: public RiskQuantity {

        static CClassConstSP const TYPE;

        static void load(CClassSP& clazz) {
            REGISTER(RQ, clazz);
            SUPERCLASS(RiskQuantity);
            EMPTY_SHELL_METHOD(DefaultConstructor<RQ>::iObject);
        }

        RQ(HypotheticalQuantityArrayConstSP parameters):
            RiskQuantity(TYPE, parameters)
        {}

        RQ(): RiskQuantity(TYPE) {}

        double value(const DoubleArray& vals,
                     const DoubleArray& dists) const {
            TRACE_METHOD;
            if (Maths::isZero(dists[1] - dists[0])) {
                TRACE("Oops, distance is zero");
                throw ModelException("Zero divisor in two-sided first deriv");
            }
            TRACE_SHOW((vals[1] - vals[0]) / (dists[1] - dists[0]));
            return (vals[1] - vals[0]) / (dists[1] - dists[0]);
        }
    };

    static CClassConstSP const TYPE;

    TwoSidedFirstDerivative(): CObject(TYPE) {}

    RiskQuantitySP riskQuantity(IResultsFunctionConstSP function,
                                IRiskAxisConstSP axis,
                                double coeff) const {
        TRACE_METHOD;

        HypotheticalQuantityArraySP ps(new HypotheticalQuantityArray(2));

        (*ps)[0] = HypotheticalQuantity::SP(axis->hypothesis(-coeff), function);
        (*ps)[1] = HypotheticalQuantity::SP(axis->hypothesis(coeff), function);

        return RiskQuantitySP(new RQ(ps));
    }

    int order() const { return 1; }
};

CClassConstSP const TwoSidedFirstDerivative::TYPE = CClass::registerClassLoadMethod(
    "TwoSidedFirstDerivative", typeid(TwoSidedFirstDerivative), load);

typedef TwoSidedFirstDerivative::RQ TwoSidedFirstDerivative_RQ;
CClassConstSP const TwoSidedFirstDerivative_RQ::TYPE = CClass::registerClassLoadMethod(
    "TwoSidedFirstDerivative::RQ", typeid(TwoSidedFirstDerivative_RQ), TwoSidedFirstDerivative_RQ::load);

IScalarDerivativeConstSP IScalarDerivative::twoSided() {
    static IScalarDerivativeConstSP it(new TwoSidedFirstDerivative());
    return it;
}

// 
// ********
//  second
// ********
// 

class SecondDerivative: public CObject,
                        public virtual IScalarDerivative {

    static void load(CClassSP& clazz) {
        REGISTER(SecondDerivative, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScalarDerivative);
        EMPTY_SHELL_METHOD(DefaultConstructor<SecondDerivative>::iObject);
    }

public:

    struct RQ: public RiskQuantity {

        static CClassConstSP const TYPE;

        static void load(CClassSP& clazz) {
            REGISTER(RQ, clazz);
            SUPERCLASS(RiskQuantity);
            EMPTY_SHELL_METHOD(DefaultConstructor<RQ>::iObject);
        }

        RQ(HypotheticalQuantityArrayConstSP parameters):
            RiskQuantity(TYPE, parameters)
        {}

        RQ(): RiskQuantity(TYPE) {}

        double value(const DoubleArray& vals,
                     const DoubleArray& dists) const {
            TRACE_METHOD;
            if (Maths::isZero(dists[2])) {
                TRACE("Oops, distance is zero");
                throw ModelException("Zero divisor in second deriv");
            }
            TRACE_SHOW((vals[1] + vals[2] - 2 * vals[0]) / (dists[/*2*/1] * dists[1]));
            return (vals[1] + vals[2] - 2 * vals[0]) / (dists[/*2*/1] * dists[1]);
        }
    };

    static CClassConstSP const TYPE;

    SecondDerivative(): CObject(TYPE) {}

    RiskQuantitySP riskQuantity(IResultsFunctionConstSP function,
                                IRiskAxisConstSP axis,
                                double coeff) const {
        TRACE_METHOD;

        HypotheticalQuantityArraySP ps(new HypotheticalQuantityArray(3));

		IHypothesisSP hypothesisPlus(copy(axis->hypothesis(coeff).get()));
		hypothesisPlus->setApproxOrder(2);
		IHypothesisSP hypothesisMinus(copy(axis->hypothesis(-coeff).get()));
		hypothesisPlus->setApproxOrder(2);

        (*ps)[0] = HypotheticalQuantity::SP(axis->hypothesis(0), function);
        (*ps)[1] = HypotheticalQuantity::SP(hypothesisPlus, function);
        (*ps)[2] = HypotheticalQuantity::SP(hypothesisMinus, function);

        return RiskQuantitySP(new RQ(ps));
    }

    int order() const { return 2; }
};

CClassConstSP const SecondDerivative::TYPE = CClass::registerClassLoadMethod(
    "SecondDerivative", typeid(SecondDerivative), load);

typedef SecondDerivative::RQ SecondDerivative_RQ;
CClassConstSP const SecondDerivative_RQ::TYPE = CClass::registerClassLoadMethod(
    "SecondDerivative::RQ", typeid(SecondDerivative_RQ), SecondDerivative_RQ::load);

IScalarDerivativeConstSP IScalarDerivative::second() {
    static IScalarDerivativeConstSP it(new SecondDerivative());
    return it;
}

// 
// ***************
//  underScenario
// ***************
// 

class ConditionedDerivative: public CObject,
                             public virtual IScalarDerivative {

    static void load(CClassSP& clazz) {
        REGISTER(ConditionedDerivative, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScalarDerivative);
        EMPTY_SHELL_METHOD(DefaultConstructor<ConditionedDerivative>::iObject);
        FIELD(derivative, "derivative");
        FIELD(scenario, "scenario");
    }

public:

    static CClassConstSP const TYPE;

    struct RQ: public RiskQuantity {

        static CClassConstSP const TYPE;

        static void load(CClassSP& clazz) {
            REGISTER(RQ, clazz);
            SUPERCLASS(RiskQuantity);
            EMPTY_SHELL_METHOD(DefaultConstructor<RQ>::iObject);
            FIELD(underlying, "underlying");
        }

        RiskQuantityConstSP underlying;

        RQ(RiskQuantityConstSP underlying,
           HypotheticalQuantityArrayConstSP parameters):
            RiskQuantity(TYPE, parameters),
            underlying(underlying)
        {}

        RQ(): RiskQuantity(TYPE) {}

        double value(const DoubleArray& vals, const DoubleArray& dsts) const {
            return underlying->value(vals, dsts);
        }
    };

    IScalarDerivativeConstSP derivative;
    IHypothesisConstSP scenario;

    ConditionedDerivative(IScalarDerivativeConstSP derivative,
                          IHypothesisConstSP scenario):
        CObject(TYPE),
        derivative(derivative),
        scenario(scenario)
    {}

    ConditionedDerivative(): CObject(TYPE) {}

    RiskQuantitySP riskQuantity(IResultsFunctionConstSP function,
                                IRiskAxisConstSP axis,
                                double coeff) const {
        RiskQuantitySP rq = derivative->riskQuantity(function, axis, coeff);
        HypotheticalQuantityArrayConstSP parameters = rq->parameters();

        HypotheticalQuantityArraySP ps(new HypotheticalQuantityArray());

        for (int h = 0; h < parameters->size(); ++h) {
            ps->push_back(HypotheticalQuantity::SP(
                scenario->then((*parameters)[h]->hypothesis()),
                (*parameters)[h]->quantity()));
        }

        return RiskQuantitySP(new RQ(rq, ps));
    }

    int order() const { return derivative->order(); }
};

CClassConstSP const ConditionedDerivative::TYPE = CClass::registerClassLoadMethod(
    "ConditionedDerivative", typeid(ConditionedDerivative), load);

typedef ConditionedDerivative::RQ ConditionedDerivative_RQ;
CClassConstSP const ConditionedDerivative_RQ::TYPE = CClass::registerClassLoadMethod(
    "ConditionedDerivative::RQ", typeid(ConditionedDerivative_RQ), ConditionedDerivative_RQ::load);

IScalarDerivativeConstSP IScalarDerivative::underScenario(
        IHypothesisConstSP scenario) const {
    return IScalarDerivativeConstSP(new ConditionedDerivative(
        IScalarDerivativeConstSP::attachToRef(this), scenario));
}

// 
// ***********************
//  averagedOverScenarios
// ***********************
// 

class ScenarioAveragedDerivative: public CObject,
                                  public virtual IScalarDerivative {

    static void load(CClassSP& clazz) {
        REGISTER(ScenarioAveragedDerivative, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScalarDerivative);
        EMPTY_SHELL_METHOD(DefaultConstructor<ScenarioAveragedDerivative>::iObject);
        FIELD(derivative, "derivative");
        FIELD(scenarioProperty, "scenarioProperty");
        FIELD(scenarioShiftSizes, "scenarioShiftSizes");
    }

public:

    static CClassConstSP const TYPE;

    struct RQ: public RiskQuantity {

        static CClassConstSP const TYPE;

        static void load(CClassSP& clazz) {
            REGISTER(RQ, clazz);
            SUPERCLASS(RiskQuantity);
            EMPTY_SHELL_METHOD(DefaultConstructor<RQ>::iObject);
            FIELD(underlying, "underlying");
            FIELD(numScenarios, "numScenarios");
        }

        RiskQuantityConstSP underlying;
        double numScenarios;

        RQ(RiskQuantityConstSP underlying,
           double numScenarios,
           HypotheticalQuantityArrayConstSP parameters):
            RiskQuantity(TYPE, parameters),
            underlying(underlying),
            numScenarios(numScenarios)
        {}

        RQ(): RiskQuantity(TYPE) {}

        double value(const DoubleArray& vals, const DoubleArray& dsts) const {
            int np = underlying->parameters()->size();

            ASSERT(vals.size() == numScenarios * np &&
                   dsts.size() == numScenarios * np);

            DoubleArray vals_s(np), dsts_s(np);

            double tot = 0;

            for (int s = 0; s < numScenarios; ++s) {
                for (int p = 0; p < np; ++p) {
                    vals_s[p] = vals[p + np*s];
                    dsts_s[p] = dsts[p + np*s];
                }

                tot += underlying->value(vals_s, dsts_s);
            }

            return tot / numScenarios;
        }
    };

    IScalarDerivativeConstSP derivative;
    IScalarRiskPropertyConstSP scenarioProperty;
    CDoubleArraySP scenarioShiftSizes;

    ScenarioAveragedDerivative(IScalarDerivativeConstSP derivative,
                               IScalarRiskPropertyConstSP scenarioProperty,
                               double scenarioShiftSize0,
                               double scenarioShiftSize1 = 0.):
        CObject(TYPE),
        derivative(derivative),
        scenarioProperty(scenarioProperty),
        scenarioShiftSizes(new CDoubleArray(1, scenarioShiftSize0))
    {
        if (scenarioShiftSize1 != 0.) {
            scenarioShiftSizes->push_back(scenarioShiftSize1);
        }
    }

    ScenarioAveragedDerivative(): CObject(TYPE) {}

    RiskQuantitySP riskQuantity(IResultsFunctionConstSP function,
                                IRiskAxisConstSP axis,
                                double coeff) const {
        RiskQuantitySP rq = derivative->riskQuantity(function, axis, coeff);
        HypotheticalQuantityArrayConstSP parameters = rq->parameters();
        IRiskAxisConstSP scenarioAxis =
            scenarioProperty->axisFor(axis->marketDataName());

        HypotheticalQuantityArraySP ps(new HypotheticalQuantityArray());

        for (int s = 0; s < scenarioShiftSizes->size(); ++s) {
            for (int h = 0; h < parameters->size(); ++h) {
                ps->push_back(HypotheticalQuantity::SP(
                    scenarioAxis->hypothesis((*scenarioShiftSizes)[s])->
                        then((*parameters)[h]->hypothesis()),
                    (*parameters)[h]->quantity()));
            }
        }

        return RiskQuantitySP(new RQ(rq, scenarioShiftSizes->size(), ps));
    }

    int order() const { return derivative->order(); }
};

CClassConstSP const ScenarioAveragedDerivative::TYPE = CClass::registerClassLoadMethod(
    "ScenarioAveragedDerivative", typeid(ScenarioAveragedDerivative), load);

typedef ScenarioAveragedDerivative::RQ ScenarioAveragedDerivative_RQ;
CClassConstSP const ScenarioAveragedDerivative_RQ::TYPE = CClass::registerClassLoadMethod(
    "ScenarioAveragedDerivative::RQ", typeid(ScenarioAveragedDerivative_RQ), ScenarioAveragedDerivative_RQ::load);

IScalarDerivativeConstSP IScalarDerivative::averagedOverScenarios(
        IScalarRiskPropertyConstSP scenarioProperty,
        double scenario1Shift,
        double scenario2Shift) const {
    return IScalarDerivativeConstSP(new ScenarioAveragedDerivative(
        IScalarDerivativeConstSP::attachToRef(this), scenarioProperty,
        scenario1Shift, scenario2Shift));
}

IScalarDerivativeConstSP IScalarDerivative::underScenario(
        IScalarRiskPropertyConstSP scenarioProperty,
        double scenarioShift) const {
    return averagedOverScenarios(scenarioProperty, scenarioShift, 0.);
}

CClassConstSP const IScalarDerivative::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IScalarDerivative", typeid(IScalarDerivative), 0);

DEFINE_TEMPLATE_TYPE(IScalarDerivativeArray);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
