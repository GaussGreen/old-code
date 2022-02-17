//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Spline2D.hpp
//
//   Description : 
//
//   Date        : 24 Oct 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_SPLINE2D_CPP
#include "edginc/Spline2D.hpp"

DRLIB_BEGIN_NAMESPACE

Spline2D::Spline2D():
CObject(TYPE){}

InterpolatorMD::InterpolantConstSP Spline2D::computeInterp(const DoubleArrayArray& xdata,
                                                           const DoubleArrayMD&    fdata) const{
    return Spline2D::InterpolantConstSP(new Spline2D::Interpolant(xdata, fdata));
}

InterpolatorMD::InterpolantConstSP Spline2D::computeInterp(const DoubleArrayArraySP& xdata,
                                                           const DoubleArrayMDSP&    fdata) const{
    return Spline2D::InterpolantConstSP(new Spline2D::Interpolant(xdata, fdata));
}

class Spline2DHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Spline2D, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(InterpolatorMD);
        EMPTY_SHELL_METHOD(defaultCtor);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new Spline2D();
    }
};

CClassConstSP const Spline2D::TYPE = CClass::registerClassLoadMethod(
    "Spline2D", typeid(Spline2D), Spline2DHelper::load);

Spline2D::Interpolant::Interpolant():
InterpolatorMD::Interpolant(TYPE){}

Spline2D::Interpolant::Interpolant(const DoubleArrayArray& x, 
                                   const DoubleArrayMD&    y):
InterpolatorMD::Interpolant(TYPE, x, y){
    validatePop2Object();
}

Spline2D::Interpolant::Interpolant(const DoubleArrayArraySP& x, 
                                   const DoubleArrayMDSP&    y):
InterpolatorMD::Interpolant(TYPE, x, y){
    validatePop2Object();
}

void Spline2D::Interpolant::validatePop2Object(){
    try{
        Interpolant::validatePop2Object();
        if (dim != 2){
            throw ModelException("Spline2D::Interpolant::validatePop2Object",
                                 "dimension should be 2; got"
                                 + Format::toString(dim));
        }
        createSpline();
    }
    catch(exception& e){
        throw ModelException(e, "Spline2D::Interpolant::validatePop2Object");
    }
}

void Spline2D::Interpolant::createSpline(){
    static const string method("Spline2D::Interpolant::createSpline");
    try{
        DoubleArray& xdata = (*x)[0];
        int num_xdata = xdata.size();
        DoubleArray& ydata = (*x)[1];
        int num_ydata = ydata.size();
        // copy y to a contiguous array
        // need to think of speeding that operation up
        vector<double> fdata(num_xdata * num_ydata);
        DoubleArrayMD::Iterator y_it(y->begin());
        vector<double>::iterator fdata_it(fdata.begin());
        for(; fdata_it != fdata.end(); ++y_it, ++fdata_it){
            *fdata_it = *y_it;
        }
        // call imsl
        spline = IMSLSplineSP(new IMSLSpline(imsl_d_spline_2d_interp(num_xdata, &xdata[0],
                                                                     num_ydata, &ydata[0],
                                                                     &fdata[0], 0)));
        IMSLError::throwExceptionIfError();
        if (!spline->sp){
            throw ModelException(method, "imsl_d_spline_2d_interp failed");
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

IObject* Spline2D::Interpolant::clone() const{
    IObjectSP objcopy(CObject::clone());
    Spline2D::Interpolant& thecopy = 
        dynamic_cast<Spline2D::Interpolant&>(*objcopy);
    thecopy.spline = spline;    // not a deep copy
    return &thecopy;
}

double Spline2D::Interpolant::value(const DoubleArray& xx) const {
    try{
        // call imsl
        double rtn = imsl_d_spline_2d_value(xx[0], xx[1], spline->sp, 0);
        IMSLError::throwExceptionIfError();
        return rtn;
    }
    catch(exception& e){
        throw ModelException(e, "Spline2D::Interpolant::value");
    }
}

class Spline2D_InterpolantHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Spline2D::Interpolant, clazz);
        SUPERCLASS(InterpolatorMD::Interpolant);
        EMPTY_SHELL_METHOD(defaultCtor);
    }

    static IObject* defaultCtor(){
        return new Spline2D::Interpolant();
    }
};

CClassConstSP const Spline2D::Interpolant::TYPE = CClass::registerClassLoadMethod(
    "Spline2D::Interpolant", typeid(Spline2D::Interpolant), Spline2D_InterpolantHelper::load);

DRLIB_END_NAMESPACE
