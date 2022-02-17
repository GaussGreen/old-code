//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TabulatedFunc.cpp
//
//   Description : 2 arrays of doubles {x} and {f(x)}
//                 Given x, return f(x) with perhaps interpolation
//
//   Date        : May 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_TABULATEDFUNC_CPP
#include "edginc/TabulatedFunc.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Nrfns.hpp"

DRLIB_BEGIN_NAMESPACE

const string TabulatedFunc::INTERP_NONE   = "N";
const string TabulatedFunc::INTERP_LINEAR = "L";
const string TabulatedFunc::INTERP_STAIRS = "S";



TabulatedFunc::TabulatedFunc(
    const DoubleArray&   xArray,
    const DoubleArray&   fArray,
    const string&        interp) : CObject(TYPE),
                                   xArray(xArray),
                                   fArray(fArray),
                                   interp(interp) {
    validatePop2Object();
}

TabulatedFunc::~TabulatedFunc() {}

int TabulatedFunc::length() const {
    return xArray.getLength();
}

/** v. basic interpolator */
double TabulatedFunc::interpolate(double x) const {
    static const string method = "TabulatedFunc::interpolate";
    try {
        if (x < xArray[0] ||
            x > xArray[xArray.size()-1]) {
            throw ModelException(method, "Defined only between " + Format::toString(xArray[0]) +
                                 " and " + Format::toString(xArray[xArray.size()-1]) + " but request for value at " +
                                 Format::toString(x));
        } 
        unsigned long idx;
        // NR. Bisection. 
        locate(&(xArray[0])-1, // NR searches in an array xx[1,...,n]
               xArray.size(), 
               x, 
               &idx);
        idx--; // bring back to [0]-based arrays
        // x lies between xArray[idx] and xArray[idx+1]
        switch (interp[0]) {
        case 'N':
            if (Maths::equals(x, xArray[idx])) {
                return fArray[idx];
            } else if (Maths::equals(x, xArray[idx+1])) {
                return fArray[idx];
            } else {
                throw ModelException(method, "interpNone and ordinate value supplied " + 
                                     Format::toString(x) + " is not found!" );
            }
            break;
        case 'L':
            return (fArray[idx]*(xArray[idx+1]-x) + fArray[idx+1]*(x-xArray[idx]))/(xArray[idx+1]-xArray[idx]);
            break;
        case 'S':
            return fArray[idx];
            break;
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
    return 0; // stop compiler warning
}

void TabulatedFunc::validatePop2Object() {
    static const string method = "TabulatedFunc::validatePop2Object";
    try {

        if ( xArray.size() != fArray.size() ) { 
           throw ModelException(method,
                                "number of ordinates (" + Format::toString(xArray.size()) + 
                                ") must be the same as number of values (" +
                                Format::toString(fArray.size()) + ").");
        }

        for (int i = 1; i < length(); i++) {
            if (xArray[i] <= xArray[i-1]) {
                throw ModelException(method,
                                     "ordinates not in increasing order (" +
                                     Format::toString(xArray[i]) + " <= " +
                                     Format::toString(xArray[i-1]) + ")");
            }
        }

        if (interp != TabulatedFunc::INTERP_NONE &&
            interp != TabulatedFunc::INTERP_LINEAR &&
            interp != TabulatedFunc::INTERP_STAIRS) 
        {
            throw ModelException(method, "invalid interp type " + interp+
                                 ". Must be (N)one, (L)inear or (S)tairs");
        }
    }            
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

TabulatedFunc::TabulatedFunc(): CObject(TYPE) {
    //empty
}

class TabulatedFuncHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(TabulatedFunc, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultTabulatedFunc);
        FIELD(xArray, "ordinates");
        FIELD(fArray, "values");
        FIELD(interp, "interpolation style");
        Addin::registerConstructor("TABULATED_FUNC",
                                   Addin::UTILITIES,
                                   "Creates a TabulatedFunc",
                                   TabulatedFunc::TYPE);
    }

    static IObject* defaultTabulatedFunc(){
        return new TabulatedFunc();
    }
};

CClassConstSP const TabulatedFunc::TYPE = CClass::registerClassLoadMethod(
    "TabulatedFunc", typeid(TabulatedFunc), TabulatedFuncHelper::load);

DRLIB_END_NAMESPACE
