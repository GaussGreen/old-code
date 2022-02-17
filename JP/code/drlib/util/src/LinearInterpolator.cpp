//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LinearInterpolator.cpp
//
//   Description : 
//
//   Date        : 06 June 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_LINEARINTERPOLATOR_CPP
#include "edginc/LinearInterpolator.hpp"

DRLIB_BEGIN_NAMESPACE

CDoubleArray computeSlopes(const CDoubleArray& xx, const CDoubleArray& yy) {
    static const string routine("computeSlopes");
    
    try {
        int n = xx.size();
        if(n != yy.size()) {
            throw ModelException("Sizes of x and y arrays are not equal: " +
                                 Format::toString(xx.size()) + ", " + 
                                 Format::toString(yy.size()));

        }
    
        CDoubleArray dydx(n-1);
        int i;
        for(i = 0; i < n-1; i++) {
#if 0
            if(Maths::equals(xx[i+1], xx[i])) {
                throw ModelException("Identical points found in the x-Array " + Format::toString(xx[i]));
            } else {
                dydx[i] = (yy[i+1] - yy[i]) / (xx[i+1] - xx[i]);
            }
#else
                dydx[i] = (yy[i+1] - yy[i]) / (xx[i+1] - xx[i]);
#endif
        }

        return dydx;
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

////////////////////////////////////////////////////////////////////////////////

// EquidistantInterpolantMain
template<> CClassConstSP const EquidistantInterpolantMain<Interpolator::Interpolant>::TYPE = CClass::registerClassLoadMethod(
    "EquidistantInterpolantMain<Interpolator::Interpolant>", 
    typeid(EquidistantInterpolantMain<Interpolator::Interpolant>), 
    EquidistantInterpolantMain<Interpolator::Interpolant>::load);

template<> CClassConstSP const EquidistantInterpolantMain<Interpolator::NullBase>::TYPE = CClass::registerClassLoadMethod(
    "EquidistantInterpolantMain<Interpolator::NullBase>", 
    typeid(EquidistantInterpolantMain<Interpolator::NullBase>), 
    EquidistantInterpolantMain<Interpolator::NullBase>::load);

////////////////////////////////////////////////////////////////////////////////

// LinearInterpolantMain
template<> CClassConstSP const LinearInterpolantMain<Interpolator::Interpolant>::TYPE = CClass::registerClassLoadMethod(
    "LinearInterpolantMain<Interpolator::Interpolant>", 
    typeid(LinearInterpolantMain<Interpolator::Interpolant>), 
    LinearInterpolantMain<Interpolator::Interpolant>::load);

template<> CClassConstSP const LinearInterpolantMain<Interpolator::NullBase>::TYPE = CClass::registerClassLoadMethod(
    "LinearInterpolantMain<Interpolator::NullBase>", 
    typeid(LinearInterpolantMain<Interpolator::NullBase>), 
    LinearInterpolantMain<Interpolator::NullBase>::load);

////////////////////////////////////////////////////////////////////////////////

LinearInterpolator::LinearInterpolator():
CObject(TYPE) {}

Interpolator::InterpolantConstSP LinearInterpolator::computeInterp(
    const CDoubleArray& xdata,
    const CDoubleArray& fdata) const {
    
    static string routine = "LinearInterpolator::computeInterp";

    try {
        return LinearInterpolantConstSP(new LinearInterpolant(xdata, fdata));
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

LinearInterpolantNonVirtualConstSP LinearInterpolator::computeLinearInterpNV(
    const LinearInterpolator& interpolator,
    const CDoubleArray& xdata,
    const CDoubleArray& fdata) {
    
    static string routine = "LinearInterpolator::computeLinearInterpNV";

    try {
        return LinearInterpolantNonVirtualConstSP(new LinearInterpolantNonVirtual(xdata, fdata));    
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

class LinearInterpolatorHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LinearInterpolator, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Interpolator);
        EMPTY_SHELL_METHOD(defaultLinearInterpolator);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultLinearInterpolator(){
        return new LinearInterpolator();
    }
};

CClassConstSP const LinearInterpolator::TYPE = CClass::registerClassLoadMethod(
    "LinearInterpolator", typeid(LinearInterpolator), LinearInterpolatorHelper::load);

bool LinearInterpolatorLoad() {
    return LinearInterpolator::TYPE != NULL;
}

DRLIB_END_NAMESPACE
