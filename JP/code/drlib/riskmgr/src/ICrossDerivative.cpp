/**
 * @file ICrossDerivative.cpp
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
#include "edginc/ICrossDerivative.hpp"

DRLIB_BEGIN_NAMESPACE

ICrossDerivative::ICrossDerivative() {}
ICrossDerivative::~ICrossDerivative() {}

// 
// *******
//  cross
// *******
// 

class CrossDerivative: public CObject,
                       public virtual ICrossDerivative {

    static void load(CClassSP& clazz) {
        REGISTER(CrossDerivative, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ICrossDerivative);
        EMPTY_SHELL_METHOD(DefaultConstructor<CrossDerivative>::iObject);
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
            if (Maths::isZero(dists[1]) || Maths::isZero(dists[2]) ||
                Maths::isZero(dists[3]) || Maths::isZero(dists[4])) {
                TRACE("Oops, distance is zero");
                throw ModelException("Zero divisor in cross deriv");
            }

            // 00 0- 0+ -0 +0 -- ++

            double g0 = (vals[2] + vals[1] - 2 * vals[0]) / (-dists[1] * dists[2]);
            double g1 = (vals[4] + vals[3] - 2 * vals[0]) / (-dists[3] * dists[4]);

            TRACE_SHOW(g0);
            TRACE_SHOW(g1);

            TRACE_SHOW(
                   (vals[6] + vals[5] - 2 * vals[0] - (g0 * (-dists[1] * dists[2]) +
                                                       g1 * (-dists[3] * dists[4])) /
                    (0.5 * (dists[2] - dists[1]) * (dists[4] - dists[3]))));

            return (vals[6] + vals[5] - 2 * vals[0] - (g0 * (-dists[1] * dists[2]) +
                                                       g1 * (-dists[3] * dists[4]))) /
                    (0.5 * (dists[2] - dists[1]) * (dists[4] - dists[3]));
        }
    };

    static CClassConstSP const TYPE;

    CrossDerivative(): CObject(TYPE) {}

    RiskQuantitySP riskQuantity(IResultsFunctionConstSP derivand,
                                IRiskAxisConstSP axis0,
                                double c0,
                                IRiskAxisConstSP axis1,
                                double c1) const {
        TRACE_METHOD;

        HypotheticalQuantityArraySP ps(new HypotheticalQuantityArray(7));

#       define H(Q) HypotheticalQuantity::SP(Q, derivand)
        (*ps)[0] = H(axis0->hypothesis(  0)->then(axis1->hypothesis(  0)));
        (*ps)[1] = H(axis0->hypothesis(  0)->then(axis1->hypothesis(-c1)));
        (*ps)[2] = H(axis0->hypothesis(  0)->then(axis1->hypothesis(+c1)));
        (*ps)[3] = H(axis0->hypothesis(-c0)->then(axis1->hypothesis(  0)));
        (*ps)[4] = H(axis0->hypothesis(+c0)->then(axis1->hypothesis(  0)));
        (*ps)[5] = H(axis0->hypothesis(-c0)->then(axis1->hypothesis(-c1)));
        (*ps)[6] = H(axis0->hypothesis(+c0)->then(axis1->hypothesis(+c1)));
#       undef H

        return RiskQuantitySP(new RQ(ps));
    }
};

CClassConstSP const CrossDerivative::TYPE = CClass::registerClassLoadMethod(
    "CrossDerivative", typeid(CrossDerivative), load);

typedef CrossDerivative::RQ CrossDerivative_RQ;
CClassConstSP const CrossDerivative_RQ::TYPE = CClass::registerClassLoadMethod(
    "CrossDerivative::RQ", typeid(CrossDerivative_RQ), CrossDerivative_RQ::load);

ICrossDerivativeConstSP ICrossDerivative::cross() {
    static ICrossDerivativeConstSP it(new CrossDerivative());
    return it;
}

CClassConstSP const ICrossDerivative::TYPE =
    CClass::registerInterfaceLoadMethod(
        "ICrossDerivative", typeid(ICrossDerivative), 0);

DEFINE_TEMPLATE_TYPE(ICrossDerivativeArray);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
