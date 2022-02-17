//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Spline.cpp
//
//   Description : 
//
//   Date        : 05 April 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Spline.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

// IMSL CUBIC SPLINE
class CubicSplineInterpECnd_Ppoly: public InterpolantVirtual {
public:
    static CClassConstSP const TYPE;

    CubicSplineInterpECnd_Ppoly(): 
    InterpolantVirtual(TYPE), ppoly(0), amOwner(false){}

    CubicSplineInterpECnd_Ppoly(Imsl_d_ppoly* ppoly, 
                                const double* xx,
                                const double* yy,
                                int           nbPoints):
    InterpolantVirtual(TYPE), ppoly(ppoly), amOwner(true) {
        n = nbPoints;
        x = CDoubleArray(n);
        y = CDoubleArray(n);
        for (int i = 0; i < n; ++i){
            x[i] = xx[i];
            y[i] = yy[i];
        }
    }

    ~CubicSplineInterpECnd_Ppoly(){
        if (amOwner){
            free(ppoly);
        }
    }

    virtual IObject* clone() const{
        CubicSplineInterpECnd_Ppoly& thisCopy = 
            dynamic_cast<CubicSplineInterpECnd_Ppoly&>(*CObject::clone());
        thisCopy.ppoly = ppoly;
        thisCopy.amOwner = false;
        return &thisCopy;
    }

    virtual double value(double x) const{
        const static string method = "CubicSplineInterpECnd_Ppoly::value";

        double rtn = imsl_d_cub_spline_value(x, ppoly, 0);

        IMSLError::throwExceptionIfError();

        return rtn;
    }

    virtual double value(double x,
                         int    deriv) const{
        const static string method = "CubicSplineInterpECnd_Ppoly::value";

        double rtn = imsl_d_cub_spline_value (x, ppoly,
                                        IMSL_DERIV, deriv,
                                        0);

        IMSLError::throwExceptionIfError();

        return rtn;
    }

    virtual void value(const double* xvec,           
                       int           nvec,
                       int           deriv,
                       double*       valuevec) const{
        const static string method = "CubicSplineInterpECnd_Ppoly::value";

        imsl_d_cub_spline_value (0.0 /* not used */, ppoly,
                                 IMSL_DERIV, deriv,
                                 IMSL_GRID_USER, nvec, xvec, valuevec,
                                 0);

        IMSLError::throwExceptionIfError();
    }

    virtual void value(const CDoubleArray& xvec,
                       int                 deriv,
                       CDoubleArray&       valuevec) const{
        const static string method = "CubicSplineInterpECnd_Ppoly::value";

        try{
            if (xvec.size() != valuevec.size()){
                throw ModelException(method, 
                                     "Size mismatch between xvec ("+ 
                                     Format::toString(xvec.size()) +
                                     ") and valuevec ("+ 
                                     Format::toString(valuevec.size())+ ")");
            }
            value(&xvec[0],
                  xvec.size(),
                  deriv,
                  &valuevec[0]);
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CubicSplineInterpECnd_Ppoly, clazz);
        SUPERCLASS(InterpolantVirtual);
        EMPTY_SHELL_METHOD(defaultCtor);
    }

    static IObject* defaultCtor(){
        return new CubicSplineInterpECnd_Ppoly();
    }

    Imsl_d_ppoly* ppoly; // $unregistered
    bool amOwner; // $unregistered
};

CClassConstSP const CubicSplineInterpECnd_Ppoly::TYPE = CClass::registerClassLoadMethod(
    "CubicSplineInterpECnd_Ppoly", typeid(CubicSplineInterpECnd_Ppoly), CubicSplineInterpECnd_Ppoly::load);

Interpolator::InterpolantConstSP CubicSplineInterpECnd::computeInterp(const CDoubleArray& xdata,
                                                                      const CDoubleArray& fdata) const{
    const static string method = "CubicSplineInterpECnd::computeInterp";

    try{
        if (xdata.size() != fdata.size()){
            throw ModelException(method, 
                                 "Size mismatch between xdata ("+ 
                                 Format::toString(xdata.size()) +
                                 ") and fdata ("+ 
                                 Format::toString(fdata.size())+ ")");
        }

        return computeInterp(&xdata[0],
                             &fdata[0],
                             xdata.size());

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

Interpolator::InterpolantConstSP CubicSplineInterpECnd::computeInterp(const double* inxdata,
                                                                      const double* infdata,
                                                                      int           ndata) const{
    const static string method = "CubicSplineInterpECnd::computeInterp";

    try{
        Imsl_d_ppoly* ppoly = NULL;

        double* xdata = const_cast<double*>(inxdata);
        double* fdata = const_cast<double*>(infdata);

        if (ileft != 0 && iright != 0){
            ppoly = imsl_d_cub_spline_interp_e_cnd (ndata,
                                                    xdata,
                                                    fdata,
                                                    IMSL_LEFT, ileft, left,
                                                    IMSL_RIGHT, iright, right,
                                                    0);
        }
        else if (ileft == 0){
            ppoly = imsl_d_cub_spline_interp_e_cnd (ndata,
                                                    xdata,
                                                    fdata,
                                                    IMSL_RIGHT, iright, right,
                                                    0);
        }
        else if (iright == 0){
            ppoly = imsl_d_cub_spline_interp_e_cnd (ndata,
                                                    xdata,
                                                    fdata,
                                                    IMSL_LEFT, ileft, left,
                                                    0);
        }
        else{
            ppoly = imsl_d_cub_spline_interp_e_cnd (ndata,
                                                    xdata,
                                                    fdata,
                                                    0);
        }

        IMSLError::throwExceptionIfError();

        if (!ppoly){
            throw ModelException(method, "Interpolant cannot be computed");            
        }

        return Interpolator::InterpolantSP(
                    new CubicSplineInterpECnd_Ppoly(ppoly, xdata, fdata, ndata));

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

CubicSplineInterpECnd::CubicSplineInterpECnd():
CObject(TYPE),
ileft(0), left(0.0),
iright(0), right(0.0){}

CubicSplineInterpECnd::CubicSplineInterpECnd(int ileft, double left,
                                             int iright, double right):
CObject(TYPE),
ileft(ileft), left(left),
iright(iright), right(right){
    validatePop2Object();
}

void CubicSplineInterpECnd::validatePop2Object(){
    const static string method = "CubicSplineInterpECnd::validatePop2Object";

    if (!(0 <= ileft && ileft <= 2)){
        throw ModelException(method, 
                             "ileft should be equal to 0 ('not a knot'), 1 ('1st deriv cdn') "
                             "or 2 ('2nd deriv cdn'); got " + Format::toString(ileft) + ".");
    }
    if (!(0 <= iright && iright <= 2)){
        throw ModelException(method, 
                             "iright should be equal to 0 ('not a knot'), 1 ('1st deriv cdn') "
                             "or 2 ('2nd deriv cdn'); got " + Format::toString(iright) + ".");
    }
}

class CubicSplineInterpECndHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CubicSplineInterpECnd, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Interpolator);
        EMPTY_SHELL_METHOD(defaultCubicSplineInterpECnd);
        FIELD(ileft, "ileft");
        FIELD_MAKE_OPTIONAL(ileft);
        FIELD(left, "left");
        FIELD_MAKE_OPTIONAL(left);
        FIELD(iright, "iright");
        FIELD_MAKE_OPTIONAL(iright);
        FIELD(right, "right");
        FIELD_MAKE_OPTIONAL(right);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCubicSplineInterpECnd(){
        return new CubicSplineInterpECnd();
    }
};

CClassConstSP const CubicSplineInterpECnd::TYPE = CClass::registerClassLoadMethod(
    "CubicSplineInterpECnd", typeid(CubicSplineInterpECnd), CubicSplineInterpECndHelper::load);

// IMSL SHAPE PRESERVING CUBIC SPLINE
class CubicShapePresSplineInterpECnd_Ppoly: public InterpolantVirtual {
public:
    static CClassConstSP const TYPE;

    CubicShapePresSplineInterpECnd_Ppoly(): 
    InterpolantVirtual(TYPE), ppoly(0), amOwner(false){}

    CubicShapePresSplineInterpECnd_Ppoly(Imsl_d_ppoly* ppoly,
                                         const double* xx,
                                         const double* yy,
                                         int           nbPoints):
    InterpolantVirtual(TYPE), ppoly(ppoly), amOwner(true) {
        n = nbPoints;
        x = CDoubleArray(n);
        y = CDoubleArray(n);
        for (int i = 0; i < n; ++i){
            x[i] = xx[i];
            y[i] = yy[i];
        }
    }

    ~CubicShapePresSplineInterpECnd_Ppoly(){
        if (amOwner){
            free(ppoly);
        }
    }

    virtual IObject* clone() const{
        CubicShapePresSplineInterpECnd_Ppoly& thisCopy = 
            dynamic_cast<CubicShapePresSplineInterpECnd_Ppoly&>(*CObject::clone());
        thisCopy.ppoly = ppoly;
        thisCopy.amOwner = false;
        return &thisCopy;
    }

    virtual double value(double x) const{
        const static string method =
            "CubicShapePresSplineInterpECnd_Ppoly::value";

        double rtn = imsl_d_cub_spline_value(x, ppoly, 0);

        IMSLError::throwExceptionIfError();

        return rtn;
    }

    virtual double value(double x,
                         int    deriv) const{
        const static string method = 
            "CubicShapePresSplineInterpECnd_Ppoly::value";

        double rtn = imsl_d_cub_spline_value (x, ppoly,
                                        IMSL_DERIV, deriv,
                                        0);

        IMSLError::throwExceptionIfError();

        return rtn;
    }

    virtual void value(const double* xvec,           
                       int           nvec,
                       int           deriv,
                       double*       valuevec) const{
        const static string method =
            "CubicShapePresSplineInterpECnd_Ppoly::value";

        imsl_d_cub_spline_value (0.0 /* not used */, ppoly,
                                 IMSL_DERIV, deriv,
                                 IMSL_GRID_USER, nvec, xvec, valuevec,
                                 0);

        IMSLError::throwExceptionIfError();
    }

    virtual void value(const CDoubleArray& xvec,
                       int                 deriv,
                       CDoubleArray&       valuevec) const{
        const static string method = 
            "CubicShapePresSplineInterpECnd_Ppoly::value";

        try{
            if (xvec.size() != valuevec.size()){
                throw ModelException(method, 
                                     "Size mismatch between xvec ("+ 
                                     Format::toString(xvec.size()) +
                                     ") and valuevec ("+ 
                                     Format::toString(valuevec.size())+ ")");
            }
            value(&xvec[0],
                  xvec.size(),
                  deriv,
                  &valuevec[0]);
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CubicShapePresSplineInterpECnd_Ppoly, clazz);
        SUPERCLASS(InterpolantVirtual);
        EMPTY_SHELL_METHOD(defaultCtor);
    }

    static IObject* defaultCtor(){
        return new CubicShapePresSplineInterpECnd_Ppoly();
    }

    Imsl_d_ppoly* ppoly; // $unregistered
    bool amOwner; // $unregistered
};

CClassConstSP const CubicShapePresSplineInterpECnd_Ppoly::TYPE = CClass::registerClassLoadMethod(
    "CubicShapePresSplineInterpECnd_Ppoly", typeid(CubicShapePresSplineInterpECnd_Ppoly), CubicShapePresSplineInterpECnd_Ppoly::load);

Interpolator::InterpolantConstSP CubicShapePresSplineInterpECnd::computeInterp(const CDoubleArray& xdata,
                                                                               const CDoubleArray& fdata) const{
    const static string method = "CubicShapePresSplineInterpECnd::getPPoly";

    try{
        if (xdata.size() != fdata.size()){
            throw ModelException(method, 
                                 "Size mismatch between xdata ("+ 
                                 Format::toString(xdata.size()) +
                                 ") and fdata ("+ 
                                 Format::toString(fdata.size())+ ")");
        }

        return computeInterp(&xdata[0],
                             &fdata[0],
                             xdata.size());

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

Interpolator::InterpolantConstSP CubicShapePresSplineInterpECnd::computeInterp(const double* inxdata,
                                                                               const double* infdata,
                                                                               int           ndata) const{
    const static string method = "CubicShapePresSplineInterpECnd::computeInterp";
    try{
        Imsl_d_ppoly* ppoly = NULL;

        double* xdata = const_cast<double*>(inxdata);
        double* fdata = const_cast<double*>(infdata);
        if (isConcave){
            ppoly = imsl_d_cub_spline_interp_shape (ndata,
                                                    xdata,
                                                    fdata,
                                                    IMSL_CONCAVE,
                                                    0);
        }
        else{
            ppoly = imsl_d_cub_spline_interp_shape (ndata,
                                                    xdata,
                                                    fdata,
                                                    0);
        }
        IMSLError::throwExceptionIfError();
        if (!ppoly){
            throw ModelException(method, "Interpolant cannot be computed");            
        }
        return Interpolator::InterpolantSP(
                    new CubicShapePresSplineInterpECnd_Ppoly(ppoly, xdata, fdata, ndata));

    } catch (exception& e){
        throw ModelException(e, method);
    }
}

CubicShapePresSplineInterpECnd::CubicShapePresSplineInterpECnd():
CObject(TYPE),
isConcave(false){}

CubicShapePresSplineInterpECnd::CubicShapePresSplineInterpECnd(bool isConcave):
CObject(TYPE),
isConcave(isConcave){}

class CubicShapePresSplineInterpECndHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(CubicShapePresSplineInterpECnd, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Interpolator);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(isConcave, "isConcave")
        clazz->setPublic(); // make visible to EAS/spreadsheet

        Addin::registerConstructor("SHAPE_PRESERVING_SPLINE",
                                   Addin::UTILITIES,
                                   "Creates a handle to a numerical recipe spline",
                                   CubicShapePresSplineInterpECnd::TYPE);
    }

    static IObject* defaultCtor(){
        return new CubicShapePresSplineInterpECnd();
    }
};

CClassConstSP const CubicShapePresSplineInterpECnd::TYPE = CClass::registerClassLoadMethod(
    "CubicShapePresSplineInterpECnd", typeid(CubicShapePresSplineInterpECnd), CubicShapePresSplineInterpECndHelper::load);


// NR CUBIC SPLINE
NRSpline::NRSpline():
CObject(TYPE),
ileft(1), left(0.0),
iright(1), right(0.0){}

NRSpline::NRSpline(int ileft, double left,
                   int iright, double right):
CObject(TYPE),
ileft(ileft), left(left),
iright(iright), right(right){
    validatePop2Object();
}

void NRSpline::validatePop2Object(){
    const static string method = "NRSpline::validatePop2Object";
    if (!(1 <= ileft && ileft <= 2)){
        throw ModelException(method, 
                             "ileft should be equal to 1 ('1st deriv cdn') "
                             "or 2 ('2nd deriv cdn'); got " + Format::toString(ileft) + ".");
    }
    if (ileft == 2 && !Maths::isZero(left)){
        throw ModelException(method, 
                             "left should be equal to 0.0 when ileft == 2; got " + Format::toString(left) + ".");
    }
    ypleft = (ileft == 1 ? left : 1e30);    // see nr's 'spline' for explanation
    if (!(1 <= iright && iright <= 2)){
        throw ModelException(method, 
                             "iright should be equal to 1 ('1st deriv cdn') "
                             "or 2 ('2nd deriv cdn'); got " + Format::toString(iright) + ".");
    }
    if (iright == 2 && !Maths::isZero(right)){
        throw ModelException(method, 
                             "right should be equal to 0.0 when iright == 2; got " + Format::toString(right) + ".");
    }
    ypright = (iright == 1 ? right : 1e30);    // see nr's 'spline' for explanation
}

class NRSplineHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(NRSpline, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Interpolator);
        EMPTY_SHELL_METHOD(defaultNRSpline);
        FIELD(ileft, "ileft");
        FIELD_MAKE_OPTIONAL(ileft);
        FIELD(left, "left");
        FIELD_MAKE_OPTIONAL(left);
        FIELD(iright, "iright");
        FIELD_MAKE_OPTIONAL(iright);
        FIELD(right, "right");
        FIELD_MAKE_OPTIONAL(right);
        // transients
        FIELD(ypleft, "");
        FIELD_MAKE_TRANSIENT(ypleft)
        FIELD(ypright, "");
        FIELD_MAKE_TRANSIENT(ypright)
        clazz->setPublic(); // make visible to EAS/spreadsheet

        Addin::registerConstructor("NRSPLINE",
                                   Addin::UTILITIES,
                                   "Creates a handle to a numerical recipe spline",
                                   NRSpline::TYPE);
    }

    static IObject* defaultNRSpline(){
        return new NRSpline();
    }
};

CClassConstSP const NRSpline::TYPE = CClass::registerClassLoadMethod(
    "NRSpline", typeid(NRSpline), NRSplineHelper::load);

NRSpline::Interpolant::Interpolant(const DoubleArray&  xx,
                                   const DoubleArray&  yy):
InterpolantVirtual(TYPE, xx, yy){
    y2.resize(n);
}

NRSpline::Interpolant::Interpolant(const double* xx,
                                   const double* yy,
                                   int           nn):
InterpolantVirtual(TYPE){
    n = nn;
    x = CDoubleArray(nn);
    y = CDoubleArray(nn);
    for (int i = 0; i < n; ++i){
        x[i] = xx[i];
        y[i] = yy[i];
    }
    y2.resize(n);
}

NRSpline::Interpolant::Interpolant():
InterpolantVirtual(TYPE){}

class NRSpline_InterpolantHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(NRSpline::Interpolant, clazz);
        SUPERCLASS(InterpolantVirtual);
        EMPTY_SHELL_METHOD(defaultCtor);
        // transient
        FIELD(y2, "");
        FIELD_MAKE_TRANSIENT(y2)
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new NRSpline::Interpolant();
    }
};

CClassConstSP const NRSpline::Interpolant::TYPE = CClass::registerClassLoadMethod(
    "NRSpline::Interpolant", typeid(NRSpline::Interpolant), NRSpline_InterpolantHelper::load);

typedef NRSpline::InterpolantArray NRSplineInterpolantArray; // MSVC7 bug
DEFINE_TEMPLATE_TYPE_WITH_NAME("NRSpline::InterpolantArray", NRSplineInterpolantArray);

NRSpline::InterpolantConstSP NRSpline::computeInterp(const NRSpline&     spliner,
                                                     const CDoubleArray& xdata,
                                                     const CDoubleArray& fdata){
    NRSpline::InterpolantSP interpolant(new NRSpline::Interpolant(xdata, fdata));
    int n = xdata.size();
    double* x = const_cast<double*>(&xdata[0]-1);
    double* y = const_cast<double*>(&fdata[0]-1);
    spline(x, y, n, spliner.ypleft, spliner.ypright, &interpolant->y2[0]-1);
    return interpolant;
}

NRSpline::InterpolantConstSP NRSpline::computeInterp(const NRSpline& spliner,
                                                     const double*   xdata,
                                                     const double*   fdata,
                                                     int             ndata){
    NRSpline::InterpolantSP interpolant(new NRSpline::Interpolant(xdata, fdata, ndata));
    int n = ndata;
    double* x = const_cast<double*>(&xdata[0]-1);
    double* y = const_cast<double*>(&fdata[0]-1);
    spline(x, y, n, spliner.ypleft, spliner.ypright, &interpolant->y2[0]-1);
    return interpolant;
}

Interpolator::InterpolantConstSP NRSpline::computeInterp(const CDoubleArray& xdata,
                                                         const CDoubleArray& fdata) const{
    return computeInterp(*this, xdata, fdata);
}

Interpolator::InterpolantConstSP NRSpline::computeInterp(const double* xdata,
                                                         const double* fdata,
                                                         int           ndata) const{
    return computeInterp(*this, xdata, fdata, ndata);
}

/** default implementation - calls nr's 'splint' */
double NRSpline::Interpolant::value(double xx) const{
    double rtn;
    double* xa = const_cast<double*>(&x[0]-1);
    double* ya = const_cast<double*>(&y[0]-1);
    double* y2a = const_cast<double*>(&y2[0]-1);
    splint(xa, ya, y2a, n, xx, &rtn);
    return rtn;
}

#define nrerror(text) throw ModelException("NR error!", text)
#define float double // all NR float will default to double except for ran1 & gasdev(and others rans)
/** Modified version of nr's splint - (i) takes a lookup function
    and an initial guess index; (ii) computes derivatives */
        // modified version of nr's splint
typedef void (*TLookupFunc)(float xx[], unsigned long n, float x, unsigned long *jlo);
void splint(float xa[], float ya[], float y2a[], int n, 
            float x, unsigned long *klo, TLookupFunc lookup, int deriv, float *y){
    int khi;
    float h,b,a;

    lookup(xa, n, x, klo);
    if (*klo==0) *klo=1;
    if ((int)*klo==n) *klo=n-1;
    khi=*klo+1;
    h=xa[khi]-xa[*klo];
    if (h == 0.0) nrerror("Bad xa input to routine splint");
    a=(xa[khi]-x)/h;
    b=(x-xa[*klo])/h;
    if (deriv==0){
        *y=a*ya[*klo]+b*ya[khi]+((a*a*a-a)*y2a[*klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
        return;
    }
    if (deriv==1){
        *y=(ya[khi]-ya[*klo])/h+((3.0*b*b-1.0)*y2a[khi]-(3.0*a*a-1.0)*y2a[*klo])*h/6.0;
        return;
    }
    if (deriv==2){
        *y=b*y2a[khi]+a*y2a[*klo];
        return;
    }
    if (deriv==3){
        *y=(y2a[khi]-y2a[*klo])/h;
        return;
    }
    *y=0.0;
}

double NRSpline::Interpolant::valueWithGuess(double xx,
                                             int&   guess,
                                             int    deriv) const{
    unsigned long jlo = static_cast<unsigned long>(guess+1);
    double rtn;
    double* xa = const_cast<double*>(&x[0]-1);
    double* ya = const_cast<double*>(&y[0]-1);
    double* y2a = const_cast<double*>(&y2[0]-1);
    splint(xa, ya, y2a, n, 
           xx, &jlo, (TLookupFunc)hunt,
           deriv, &rtn);
    guess = jlo - 1;
    return rtn;
}

double NRSpline::Interpolant::value(double xx,
                                    int    deriv) const{
    unsigned long jlo = 0;
    double rtn;
    double* xa = const_cast<double*>(&x[0]-1);
    double* ya = const_cast<double*>(&y[0]-1);
    double* y2a = const_cast<double*>(&y2[0]-1);
    splint(xa, ya, y2a, n, 
           xx, &jlo, (TLookupFunc)locate,
           deriv, &rtn);
    return rtn;
}

void NRSpline::Interpolant::value(const CDoubleArray& xvec,           
                                  int                 deriv,
                                  CDoubleArray&       valuevec) const{
    int guess = -1;
    int i = 0;
    for (; i < xvec.size(); ++i){
        valuevec[i] = valueWithGuess(xvec[i],
                                     guess,
                                     deriv);
    }
}

bool SplineLoad() {
    return CubicShapePresSplineInterpECnd::TYPE != NULL &&
           CubicShapePresSplineInterpECnd_Ppoly::TYPE != NULL &&
           CubicSplineInterpECnd::TYPE != NULL &&
           CubicSplineInterpECnd_Ppoly::TYPE != NULL &&
           NRSpline::TYPE != NULL &&
           NRSpline::Interpolant::TYPE != NULL &&
           NRSplineInterpolantArray::TYPE != NULL;
}

DRLIB_END_NAMESPACE
