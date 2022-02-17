//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LinearInterpolatorMDMD.hpp
//
//   Description : 
//
//   Date        : 08 Oct 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_LINEARINTERPOLATORMD_CPP
#include "edginc/LinearInterpolatorMD.hpp"
#include "edginc/Nrfns.hpp"

DRLIB_BEGIN_NAMESPACE

LinearInterpolatorMD::LinearInterpolatorMD():
CObject(TYPE){}

LinearInterpolatorMD::InterpolantConstSP LinearInterpolatorMD::computeLinearInterp(const DoubleArrayArray& xdata,
                                                                                   const DoubleArrayMD&    fdata) const{
    return LinearInterpolatorMD::InterpolantConstSP(new LinearInterpolatorMD::Interpolant(xdata, fdata));
}

LinearInterpolatorMD::InterpolantConstSP LinearInterpolatorMD::computeLinearInterp(const DoubleArrayArraySP& xdata,
                                                                                   const DoubleArrayMDSP&    fdata) const{
    return LinearInterpolatorMD::InterpolantConstSP(new LinearInterpolatorMD::Interpolant(xdata, fdata));
}

InterpolatorMD::InterpolantConstSP LinearInterpolatorMD::computeInterp(const DoubleArrayArray& xdata,
                                                                       const DoubleArrayMD&    fdata) const{
    return computeLinearInterp(xdata, fdata);
}

InterpolatorMD::InterpolantConstSP LinearInterpolatorMD::computeInterp(const DoubleArrayArraySP& xdata,
                                                                       const DoubleArrayMDSP&    fdata) const{
    return computeLinearInterp(xdata, fdata);
}

class LinearInterpolatorMDHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LinearInterpolatorMD, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(InterpolatorMD);
        EMPTY_SHELL_METHOD(defaultCtor);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new LinearInterpolatorMD();
    }
};

CClassConstSP const LinearInterpolatorMD::TYPE = CClass::registerClassLoadMethod(
    "LinearInterpolatorMD", typeid(LinearInterpolatorMD), LinearInterpolatorMDHelper::load);

LinearInterpolatorMD::Interpolant::Interpolant():
InterpolatorMD::Interpolant(TYPE) {}

LinearInterpolatorMD::Interpolant::Interpolant(const DoubleArrayArray& x, 
                                               const DoubleArrayMD&    y):
InterpolatorMD::Interpolant(TYPE, x, y){}

LinearInterpolatorMD::Interpolant::Interpolant(const DoubleArrayArraySP& x, 
                                               const DoubleArrayMDSP&    y):
InterpolatorMD::Interpolant(TYPE, x, y){}

void LinearInterpolatorMD::Interpolant::lookupValue(const DoubleArray& xx,
                                                    DoubleArrayMD::Index& index) const {
    int i = 0;
    for (; i < dim; ++i){
        unsigned long position;
        // NR searches in an array xx[1,...,n]
        locate(&(*x)[i][0]-1, (*x)[i].size(), xx[i], &position); // position is in [0, n]
        index[i] = position-1;
    }
}

void LinearInterpolatorMD::Interpolant::lookupValueWithGuess(const DoubleArray& xx,
                                                             DoubleArrayMD::Index& index) const {
    int i = 0;
    for (; i < dim; ++i){
        unsigned long position = index[i]+1;
        // NR searches in an array xx[1,...,n]
        hunt(&(*x)[i][0]-1, (*x)[i].size(), xx[i], &position); // position is in [0, n]
        index[i] = position-1;
    }
}

double LinearInterpolatorMD::Interpolant::value(const DoubleArray& xx) const {
    DoubleArrayMD::Index index(y->first());     // no guess => set to start
    lookupValue(xx, index);
    return valueRecurse(xx, index, 0);  // start recursion at 0
}

double LinearInterpolatorMD::Interpolant::valueWithGuess(const DoubleArray& xx,
                                                         DoubleArrayMD::Index& index) const {
    lookupValueWithGuess(xx, index);
    return valueRecurse(xx, index, 0);  // start recursion at 0
}

// Recursive implementation
double LinearInterpolatorMD::Interpolant::valueRecurse(const DoubleArray& thex,
                                                       const DoubleArrayMD::Index& index,
                                                       int i) const {
    // If we're done recursing, simply return y value at given index
    if (i >= dim){
        return (*y)[index];
    }
    // Otherwise, we need to do extrapolation or linear interpolation
    int x_idx = index[i];
    // If index is outside range, extrapolate flat
    if (x_idx < 0){
        DoubleArrayMD::Index index_at_start(index);
        index_at_start[i] = 0;
        return valueRecurse(thex, index_at_start, i+1);
    }
    if (x_idx >= y_size[i]-1){
        DoubleArrayMD::Index index_at_end(index);
        index_at_end[i] = y_size[i]-1;
        return valueRecurse(thex, index_at_end, i+1);
    }
    // If index is inside range, do linear interpolation
    int shifted_x_idx = x_idx+1;
    DoubleArrayMD::Index shifted_index(index);
    shifted_index[i] = shifted_x_idx;
    double xx = (*x)[i][x_idx];
    double shifted_xx = (*x)[i][shifted_x_idx];
    double yy = valueRecurse(thex, index, i+1);
    double shifted_yy = valueRecurse(thex, shifted_index, i+1);
    double slope = (shifted_yy - yy) / (shifted_xx - xx);
    return (yy + slope * (thex[i]-xx));
}

class LinearInterpolatorMD_InterpolantHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LinearInterpolatorMD::Interpolant, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
    }

    static IObject* defaultCtor(){
        return new LinearInterpolatorMD::Interpolant();
    }
};

CClassConstSP const LinearInterpolatorMD::Interpolant::TYPE = CClass::registerClassLoadMethod(
    "LinearInterpolatorMD::Interpolant", typeid(LinearInterpolatorMD::Interpolant), LinearInterpolatorMD_InterpolantHelper::load);

bool LinearInterpolatorMDLoad() {
    return LinearInterpolatorMD::TYPE != NULL;
}

DRLIB_END_NAMESPACE
