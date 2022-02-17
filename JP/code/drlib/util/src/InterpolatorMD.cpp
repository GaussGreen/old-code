//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InterpolatorMDMD.cpp
//
//   Description : 
//
//   Date        : 08 Oct 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_INTERPOLATORMD_CPP
#include "edginc/InterpolatorMD.hpp"


DRLIB_BEGIN_NAMESPACE

// InterpolatorMD
void InterpolatorMD::load(CClassSP& clazz){
    REGISTER_INTERFACE(InterpolatorMD, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const InterpolatorMD::TYPE = CClass::registerInterfaceLoadMethod(
    "InterpolatorMD", typeid(InterpolatorMD), load);


// InterpolatorMD::Interpolant
InterpolatorMD::Interpolant::Interpolant():
CObject(TYPE),
dim(0){}

InterpolatorMD::Interpolant::Interpolant(const CClassConstSP&    clazz):
CObject(clazz),
dim(0){}

InterpolatorMD::Interpolant::Interpolant(const CClassConstSP&    clazz,
                                         const DoubleArrayArray& xx,
                                         const DoubleArrayMD&    yy):
CObject(clazz),
x(new DoubleArrayArray(xx)),
y(new DoubleArrayMD(yy)),
dim(0){
    validatePop2Object();
}

InterpolatorMD::Interpolant::Interpolant(const CClassConstSP&      clazz,
                                         const DoubleArrayArraySP& xx,
                                         const DoubleArrayMDSP&    yy):
CObject(clazz),
x(xx),
y(yy),
dim(0){
    validatePop2Object();
}

void InterpolatorMD::Interpolant::validatePop2Object(){
    static const string method("InterpolatorMD::Interpolant::validatePop2Object");
    dim = y->dimension();
    y_size = y->size();
    if (dim != x->size()) {
        throw ModelException(method, "Dimensions of x and y arrays are not equal: " +
                                      Format::toString(dim) + ", " + Format::toString(x->size()));
    }

    int i = 0;
    for (; i < dim; ++i){
        if (y_size[i] != (*x)[i].size()) {
            throw ModelException(method, "Sizes of x and y arrays in " + 
                                          Format::toString(i) + 
                                          "-th dimension are not equal: " +
                                          Format::toString((*x)[i].size()) + ", " + Format::toString(y_size[i]));
        }
    }
}

// Interpolant
void InterpolatorMD::Interpolant::load(CClassSP& clazz){
    REGISTER(InterpolatorMD::Interpolant, clazz);
    SUPERCLASS(CObject);
    FIELD(x, "array of X arrays");
    FIELD(y, "Y multidimensional array");
    FIELD(y_size, "");
    FIELD_MAKE_TRANSIENT(y_size);
    FIELD(dim, "");
    FIELD_MAKE_TRANSIENT(dim);
}

const DoubleArrayArray& InterpolatorMD::Interpolant::getXarray() {
    return *x;
}

const DoubleArrayMD& InterpolatorMD::Interpolant::getYarray() {
    return *y;
}

const DoubleArrayMD::Size& InterpolatorMD::Interpolant::size() {
    return y_size;
}

CClassConstSP const InterpolatorMD::Interpolant::TYPE = CClass::registerClassLoadMethod(
    "InterpolatorMD::Interpolant", typeid(InterpolatorMD::Interpolant), load); 

DEFINE_TEMPLATE_TYPE(InterpolatorMD::InterpolantArray);


#if 0
/** ADDIN method for getting an Interpolant */
class GetInterpolantAddin: public CObject{
    static CClassConstSP const TYPE;

    InterpolatorMDSP interpolator;
    DoubleArray   x;
    DoubleArray   y;

    static IObjectSP getInterpolant(GetInterpolantAddin* params){
        static const string routine = "GetInterpolantAddin::getInterpolant";
        try {
            InterpolatorMD::InterpolantConstSP interpolant(params->interpolator->computeInterp(params->x, params->y));
            
            return InterpolatorMD::InterpolantSP::constCast(interpolant);

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    GetInterpolantAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetInterpolantAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetInterpolantAddin);
        FIELD(interpolator, "InterpolatorMD");
        FIELD(x, "x");
        FIELD(y, "y");

        Addin::registerClassObjectMethod("GET_INTERPOLANT",
                                         Addin::RISK,
                                         "Get interpolant",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getInterpolant);

    }

    static IObject* defaultGetInterpolantAddin(){
        return new GetInterpolantAddin();
    }
};

CClassConstSP const GetInterpolantAddin::TYPE = CClass::registerClassLoadMethod(
    "GetInterpolantAddin", typeid(GetInterpolantAddin), load);

/** ADDIN method for interpolating with some interpolant */
class GetInterpolatedValuesAddin: public CObject{
    static CClassConstSP const TYPE;

    InterpolatorMD::InterpolantSP interpolant;
    DoubleArray   x;
    int n;

    static IObjectSP computeValue(GetInterpolatedValuesAddin* params){
        static const string routine = "GetInterpolatedValuesAddin::computeValue";
        try {
            DoubleArray y(params->x.size());
            params->interpolant->value(params->x, params->n, y);
            
            return IObjectSP(y.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    GetInterpolatedValuesAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetInterpolatedValuesAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetInterpolatedValuesAddin);
        FIELD(interpolant, "Interpolant");
        FIELD(x, "x");
        FIELD(n, "n");

        Addin::registerClassObjectMethod("GET_INTERPOLATED_VALUES",
                                         Addin::RISK,
                                         "Get interpolated values",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)computeValue);

    }

    static IObject* defaultGetInterpolatedValuesAddin(){
        return new GetInterpolatedValuesAddin();
    }
};

CClassConstSP const GetInterpolatedValuesAddin::TYPE = CClass::registerClassLoadMethod(
    "GetInterpolatedValuesAddin", typeid(GetInterpolatedValuesAddin), load);
#endif

DRLIB_END_NAMESPACE
