/**
 * @file IRiskAxis.cpp
 */

#include "edginc/config.hpp"
#define QLIB_RISKAXIS_CPP
#include "edginc/OutputName.hpp"
#include "edginc/IAbstractRiskProperty.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/RiskAxis.hpp"
#include "edginc/AtomicHypothesis.hpp"
#include "edginc/CompoundHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

IRiskAxis::IRiskAxis() {}

IRiskAxis::~IRiskAxis() {}

CClassConstSP const IRiskAxis::TYPE = CClass::registerInterfaceLoadMethod(
    "IRiskAxis", typeid(IRiskAxis), 0);

DEFINE_TEMPLATE_TYPE(IRiskAxisArray);

// 
// ==========================
//  IRiskAxis::conditioned()
// ==========================
// 

class ConditionedRiskAxis: public CObject,
                           public virtual IRiskAxis,
                           public virtual RiskAxis {

    ConditionedRiskAxis(): CObject(TYPE) {}
    static IObject* emptyShell() {
        return new ConditionedRiskAxis();
    }
    static void load(CClassSP& clazz) {
        REGISTER(ConditionedRiskAxis, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRiskAxis);
        IMPLEMENTS(RiskAxis);
        EMPTY_SHELL_METHOD(emptyShell);
        FIELD(condition, "condition");
        FIELD(axis, "axis");
    }

public:

    static CClassConstSP const TYPE;

private:

    IHypothesisConstSP condition;
    IRiskAxisConstSP axis;

public:

    ConditionedRiskAxis(IHypothesisConstSP condition,
                        IRiskAxisConstSP axis):
        CObject(TYPE),
        condition(condition),
        axis(axis)
    {}

    IAbstractRiskPropertyConstSP abstractProperty() const {
        return IAbstractRiskPropertyConstSP();
#if 0
        IAbstractRiskPropertyConstSP p = axis->abstractProperty();
        return !p ? IAbstractRiskPropertySP() :
                    IAbstractRiskProperty::conditioned(condition, p);
#endif
    }

    OutputNameConstSP marketDataName() const {
        return axis->marketDataName();
    }

    IHypothesisConstSP hypothesis(double coeff) const {
        IHypothesisConstSP hypothesis = axis->hypothesis(coeff);

        AtomicHypothesisArraySP hyps(new AtomicHypothesisArray(
            condition->numAtomics() + hypothesis->numAtomics()));

        for (int h = 0; h < condition->numAtomics(); ++h) {
            (*hyps)[h] = AtomicHypothesisSP::constCast(
                             condition->atomic(h));
        }

        {for (int h = 0; h < hypothesis->numAtomics(); ++h) {
            (*hyps)[h + condition->numAtomics()] =
                AtomicHypothesisSP::constCast(hypothesis->atomic(h));
        }}

        return CompoundHypothesis::SP(hyps);
    }

    RiskAxisConstSP frozen() const {
        return RiskAxisConstSP::attachToRef(this);
    }

    IRiskAxisConstSP thawed() const {
        return IRiskAxisConstSP::attachToRef(this);
    }
};

CClassConstSP const ConditionedRiskAxis::TYPE = CClass::registerClassLoadMethod(
    "ConditionedRiskAxis", typeid(ConditionedRiskAxis), load);

IRiskAxisSP IRiskAxis::conditioned(IHypothesisConstSP condition,
                                   IRiskAxisConstSP axis) {
    return IRiskAxisSP(new ConditionedRiskAxis(condition, axis));
}

// 
// =======================
//  IRiskAxis::compound()
// =======================
// 

class CompoundRiskAxis: public CObject,
                        public virtual IRiskAxis,
                        public virtual RiskAxis {

    CompoundRiskAxis(): CObject(TYPE) {}
    static IObject* emptyShell() {
        return new CompoundRiskAxis();
    }
    static void load(CClassSP& clazz) {
        REGISTER(CompoundRiskAxis, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRiskAxis);
        IMPLEMENTS(RiskAxis);
        EMPTY_SHELL_METHOD(emptyShell);
        FIELD(axes, "axes");
        FIELD(metric, "metric");
        FIELD_MAKE_OPTIONAL(metric);
    }

public:

    static CClassConstSP const TYPE;

private:

    IRiskAxisArrayConstSP axes;
    IHypothesis::IDistanceMetricConstSP metric;

public:

    CompoundRiskAxis(IRiskAxisArrayConstSP axes,
                     IHypothesis::IDistanceMetricConstSP metric):
        CObject(TYPE),
        axes(axes),
        metric(metric)
    {}

    IAbstractRiskPropertyConstSP abstractProperty() const {
        return IAbstractRiskPropertyConstSP();
    }

    OutputNameConstSP marketDataName() const {
        return OutputNameConstSP();
    }

    IHypothesisConstSP hypothesis(double coeff) const {
        IHypothesisArraySP hyps = IHypothesisArray::SP(axes->size());

        for (int a = 0; a < axes->size(); ++a) {
            (*hyps)[a] = IHypothesisSP::constCast((*axes)[a]->hypothesis(coeff));
        }

        // This causes a memory leak with gcc 3.2

        // return CompoundHypothesis::SP(
        //     hyps,
        //     !metric ? IHypothesis::IDistanceMetric::constant(coeff) :
        //               metric);

        IHypothesis::IDistanceMetricConstSP m;
        if (!!metric) m = metric;
        else m = IHypothesis::IDistanceMetric::constant(coeff);

        return CompoundHypothesis::SP(hyps, m);
    }

    RiskAxisConstSP frozen() const {
        return RiskAxisConstSP::attachToRef(this);
    }

    IRiskAxisConstSP thawed() const {
        return IRiskAxisConstSP::attachToRef(this);
    }
};

CClassConstSP const CompoundRiskAxis::TYPE = CClass::registerClassLoadMethod(
    "CompoundRiskAxis", typeid(CompoundRiskAxis), load);

IRiskAxisSP IRiskAxis::compound(IRiskAxisArrayConstSP axes,
                                IHypothesis::IDistanceMetricConstSP metric) {
    return IRiskAxisSP(new CompoundRiskAxis(axes, metric));
}

DRLIB_END_NAMESPACE
