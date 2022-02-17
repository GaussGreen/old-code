//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OneToOneMapping.cpp
//
//   Date        : 20 June 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/OneToOneMapping.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE

class IdentityMapping: public OneToOneMapping{
public:
    virtual double operator()(double x) const{    // x in R
        return x;
    }

    virtual double inverse(double y) const{       // y in (lower, +infty)
        return y;
    }

    virtual double derivative(double x) const {    // x in R
        return 1.0;
    }

    IdentityMapping(){}
};

class HalfLineMapping: public OneToOneMapping{
public:
    virtual double operator()(double x) const{    // x in R
        x = (x < max_exp ? x : max_exp);
        return (bound + epsilon * exp(x - 1.0));
    }

    virtual double inverse(double y) const{       // y in (lower, +infty)
        return (1.0 + log(epsilon * (y - bound)));
    }

    virtual double derivative(double x) const {    // x in R
        return (epsilon * exp(x - 1.0));
    }

    HalfLineMapping(bool isUpper, double bound):
    epsilon(isUpper ? +1.0 : -1.0),
    bound(bound){}

private:
    static const double max_exp;
    double epsilon;
    double bound;
};
const double HalfLineMapping::max_exp = DBL_MAX_10_EXP / log(10.0);

class FiniteIntervalMapping: public OneToOneMapping{
public:
    virtual double operator()(double x) const{    // x in R
        return (lower + 0.5 * diff * (tanh(2.0 * x) + 1.0));
    }

    virtual double inverse(double y) const{       // y in (lower, upper)
        return (0.5 * atanh(-1.0 + 2.0 * (y - lower) / diff));
    }

    virtual double derivative(double x) const {    // x in R
        return (diff / Maths::square(cosh(2.0 * x)));
    }

    FiniteIntervalMapping(double lower, double upper):
    lower(lower), upper(upper), diff(upper - lower){}

private:
    double lower;
    double upper;
    double diff;
};

OneToOneMappingConstSP OneToOneMapping::create(const Range& range){
    const string method = "OneToOneMapping::create";
    try{
        Range::checkIsNonEmpty(range);
        Range::checkIsNotSingleton(range);
        Range::checkIsOpen(range);
        const Boundary& lower = range.getLower();
        const Boundary& upper = range.getUpper();
        if (lower.isInfinite() && upper.isInfinite()){    // (-infty, +infty)
            return OneToOneMappingConstSP(new IdentityMapping());
        }
        if (upper.isInfinite()){    // upper half line
            return OneToOneMappingConstSP(new HalfLineMapping(true, lower.getValue()));    // (-infty, upper)
        }
        if (lower.isInfinite()){    // lower half line
            return OneToOneMappingConstSP(new HalfLineMapping(false, upper.getValue()));    // (lower, +infty)
        }
        // (lower, upper)
        return OneToOneMappingConstSP(new FiniteIntervalMapping(lower.getValue(),
                                                                upper.getValue()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

DRLIB_END_NAMESPACE
