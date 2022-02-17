//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Integrator.cpp
//
//   Description : 
//
//   Date        : 23 Nov 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_INTEGRATOR_CPP
#include "edginc/Integrator.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Spline.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

typedef double (IMSLFunction)(double);

/** Class that gives out static function pointers for use with eg IMSL 
    ie use this to turn C++ types with methods into a C stlye function pointer.
    This is not thread-safe but could be made so with some simple locking.
*/
class IMSLFunctionProvider{
    int           index;
public:
    IMSLFunctionProvider(const Function1DDouble* object):
        index(getFirstNullEntry()){
        globalObjects[index] = object; // set global static variable
    }

    ~IMSLFunctionProvider(){
        // clear globalObject
        globalObjects[index] = 0;
    }

    IMSLFunction* function(){
        return globalFuncs[index];
    }

private:
    static int getFirstNullEntry(){
        // search for first null entry
        for (int i = 0; i < numGlobalObjects; i++){
            if (!globalObjects[i]){
                return i;
            }
        }
        throw ModelException(
            "IMSLFunctionProvider (in Integrator.cpp)",
            "Run out of spare static variables");
    }
    
    //// wrapped functions
    static double wrappedMethod1(double x){
        return (*globalObjects[0])(x);
    }
    static double wrappedMethod2(double x){
        return (*globalObjects[1])(x);
    }
    static double wrappedMethod3(double x){
        return (*globalObjects[2])(x);
    }
    static double wrappedMethod4(double x){
        return (*globalObjects[3])(x);
    }

    static IMSLFunction* globalFuncs[];
    static int numGlobalObjects;
    static const Function1DDouble* globalObjects[];
};
//// set up statics
IMSLFunction* IMSLFunctionProvider::globalFuncs[] = 
    {
        wrappedMethod1,
        wrappedMethod2,
        wrappedMethod3,
        wrappedMethod4
    };
int IMSLFunctionProvider::numGlobalObjects =
    sizeof(globalFuncs)/sizeof(IMSLFunction*);
const Function1DDouble* IMSLFunctionProvider::globalObjects[] = {0,0,0,0};

// Integrator1D
void Integrator1D::load(CClassSP &clazz){
    REGISTER_INTERFACE(Integrator1D, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}

CClassConstSP const Integrator1D::TYPE = CClass::registerInterfaceLoadMethod(
    "Integrator1D", typeid(Integrator1D), Integrator1D::load);

Integrator1DSP Integrator1D::createIntegrator(const string& name,
                                              double        absPrecision,
                                              double        relPrecision) {
    static const string method = "createIntegrator";
    try{
        if (CString::equalsIgnoreCase(name, "ClosedRomberg1D")) { // 1. ClosedRomberg1D
            ClosedRomberg1DSP newIntegrator(new ClosedRomberg1D(relPrecision));
            return newIntegrator;
        } else if (CString::equalsIgnoreCase(name, "OpenRomberg1D")) { // 2. OpenRomberg1D
            OpenRomberg1DSP newIntegrator(new OpenRomberg1D(relPrecision));
            return newIntegrator;
        } else if (CString::equalsIgnoreCase(name, "IntFuncInf")) { // 3. IntFuncInf
            IntFuncInfSP newIntegrator(new IntFuncInf(absPrecision, relPrecision));
            return newIntegrator;
        } else if (CString::equalsIgnoreCase(name, "IntFuncGeneral")) { // 4. IntFuncGeneral
            IntFuncGeneralSP newIntegrator(new IntFuncGeneral(absPrecision, relPrecision));
            return newIntegrator;
        } else {
            throw ModelException(method, 
                "Integrator name " + name + " not recognized as integrator.");
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}


///////////////////////////////////////////////////////////////////////


Trapez1DSimple::Trapez1DSimple(int nbSteps): CObject(TYPE), nbSteps(nbSteps) {}


double Trapez1DSimple::integrate(const Function1DDouble& func) const {
    static const string method = "Trapez1DSimple::integrate";
    try{
        const Range& interval = func.getInterval();

        Range::checkIsNonEmpty(interval);

        if (interval.isSingleton()){
            return 0.0;
        }

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        if (!lower.isClosedBracket() || lower.isInfinite()
            || !upper.isClosedBracket() || upper.isInfinite()){
            throw ModelException(method,
                                 "(Semi-)Open and / or infinite intervals are not supported; got " + interval.toString());
        }

        double a = lower.getValue();
        double b = upper.getValue();

        double del = (b - a) / nbSteps;
        
        // Interior Sum
        double intSum = 0.0;
        double x = a;
        for (int j = 0; j < nbSteps - 1; j++) {
            x += del;
            intSum += func(x);
        }
        intSum *= 2.0;
        
        // Endpoints Contribution
        double extSum = func(a) + func(b);
        
        // Integral
        double integral = 0.5 * del * (intSum + extSum);
        return integral;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}


void Trapez1DSimple::load(CClassSP& clazz){
    REGISTER(Trapez1DSimple, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(Integrator1D);
    EMPTY_SHELL_METHOD(defaultTrapez1DSimple);
    FIELD(nbSteps, "Number of steps");
}


IObject* Trapez1DSimple::defaultTrapez1DSimple(){
    return new Trapez1DSimple();
}


Trapez1DSimple::Trapez1DSimple(): CObject(TYPE) {}


Trapez1DSimple::Trapez1DSimple(CClassConstSP clazz):
CObject(clazz) {}


Trapez1DSimple::Trapez1DSimple(CClassConstSP clazz, int nbSteps):
CObject(clazz), nbSteps(nbSteps) {}


CClassConstSP const Trapez1DSimple::TYPE = CClass::registerClassLoadMethod(
    "Trapez1DSimple", typeid(Trapez1DSimple), load);


///////////////////////////////////////////////////////////////////////
                                              
                                              
// Trapez1D
const double Trapez1D::default_accuracy = 1.0e-6;
const int    Trapez1D::default_nbIterMax = 20;

/** Numerical Recipes in C page 137
    qtrap(float(*func)(float), float a, float b)

    Returns the integral of the function func
    from a to b using the trapezoidal rule */
double Trapez1D::integrate(const Function1DDouble& func) const {
    static const string method = "Trapez1D::integrate";
    
    try{
        const Range& interval = func.getInterval();

        Range::checkIsNonEmpty(interval);

        if (interval.isSingleton()){
            return 0.0;
        }

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        if (!lower.isClosedBracket() || lower.isInfinite()
            || !upper.isClosedBracket() || upper.isInfinite()){
            throw ModelException(method,
                                 "(Semi-)Open and / or infinite intervals are not supported; got " + interval.toString());
        }

        double a = lower.getValue();
        double b = upper.getValue();

        double olds, s;
        int j;

        olds = -1.0e30;
        for(j=1;j<=nbIterMax;j++) {
            s=trapzd(&func, a, b, j);
            if(j>5)
                if(Maths::isNegative(fabs(s - olds) - accuracy * fabs(olds))) 
                    return s;
            olds = s;
        }

        // Failure
        throw ModelException(method,
                             "Number of iterations exceeded " + 
                             Format::toString(nbIterMax));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Numerical Recipes in C page 137
    trapzd(float(*func)(float), float a, float b, int n)

    This routine computes the n-th stage of a 
    refinement of an extended trapezoidal rule */
double Trapez1D::trapzd(const Function1DDouble* func, const double a, const double b, const int n) const
{
    static string method = "trapzd"; 

    double x, tnm, sum, del;
    static double s;
    int it, j;

    if (n == 1) {
        return (s = 0.5 * (b-a)* ((*func)(a) + (*func)(b)));
    }
    else {
        for (it=1,j=1;j<n-1;j++) it <<= 1;
        tnm=it;
        del=(b-a)/tnm;
        x=a+0.5*del;
        for (sum=0.0,j=1;j<=it;j++,x+=del) 
            sum += (*func)(x);
        s=0.5*(s+(b-a)*sum/tnm);

        return s;
    }
}

   
Trapez1D::Trapez1D(double accuracy, int nbIterMax):
CObject(TYPE),
accuracy(accuracy),
nbIterMax(nbIterMax){}

Trapez1D::Trapez1D(int notUsed):
CObject(TYPE),
accuracy(default_accuracy),
nbIterMax(default_nbIterMax){}

class Trapez1DHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Trapez1D, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Integrator1D);
        EMPTY_SHELL_METHOD(defaultTrapez1D);
        FIELD(accuracy, "Accuracy");
        FIELD_MAKE_OPTIONAL(accuracy);
        FIELD(nbIterMax, "Maximum Nb of Iterations");
        FIELD_MAKE_OPTIONAL(nbIterMax);
    }

    static IObject* defaultTrapez1D(){
        return new Trapez1D(0);
    }
};

CClassConstSP const Trapez1D::TYPE = CClass::registerClassLoadMethod(
    "Trapez1D", typeid(Trapez1D), Trapez1DHelper::load);

// ClosedRomberg1D
const int ClosedRomberg1D::default_RichardsonPolyOrder = 9;

void ClosedRomberg1D::validatePop2Object() { /* For initializing nbIterMaxP */
    nbIterMaxP = nbIterMax + 1;
}

/** Numerical Recipes in C page 140
    qromb(float(*func)(float), float a, float b)

    Returns the integral of the function func
    from a to b using Romberg integration */
double ClosedRomberg1D::integrate(const Function1DDouble& func) const {
    static const string method = "ClosedRomberg1D::integrate";
    
    try{
        const Range& interval = func.getInterval();

        Range::checkIsNonEmpty(interval);

        if (interval.isSingleton()){
            return 0.0;
        }

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        if (!lower.isClosedBracket() || lower.isInfinite()
            || !upper.isClosedBracket() || upper.isInfinite()){
            throw ModelException(method,
                                 "(Semi-)Open and / or infinite intervals are not supported; got " + interval.toString());
        }

        double a = lower.getValue();
        double b = upper.getValue();

        double ss,dss;
        
        vector<double> s(nbIterMaxP+1);
        vector<double> h(nbIterMaxP+1);

        int j;
        h[1]=1.0;
        for (j=1;j<=nbIterMax;j++) {
            s[j]=trapzd(&func,a,b,j);
            if (j >= RichardsonPolyOrder) {
                polint(&h[j-RichardsonPolyOrder],
                       &s[j-RichardsonPolyOrder],
                       RichardsonPolyOrder,
                       0.0,
                       &ss,
                       &dss);
                if (Maths::isNegative(fabs(dss) - accuracy*fabs(ss)) ||
		    (Maths::isZero(dss) && ((fabs(dss) - accuracy*fabs(ss)) <= 0)))
                {
                    return ss;   
                }
            }
            h[j+1]=0.25*h[j];
        }

        // Failure
        throw ModelException(method,
                            "Number of iterations exceeded " + 
                            Format::toString(nbIterMax));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

ClosedRomberg1D::ClosedRomberg1D(double accuracy, int nbIterMax, int RichardsonPolyOrder):
Trapez1D(accuracy, nbIterMax),
RichardsonPolyOrder(RichardsonPolyOrder),
nbIterMaxP(nbIterMax + 1){ }

ClosedRomberg1D::ClosedRomberg1D(int notUsed):
Trapez1D(0),
RichardsonPolyOrder(default_RichardsonPolyOrder),
nbIterMaxP(default_nbIterMax + 1){}

class ClosedRomberg1DHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ClosedRomberg1D, clazz);
        SUPERCLASS(Trapez1D);
        EMPTY_SHELL_METHOD(defaultClosedRomberg1D);
        FIELD(RichardsonPolyOrder, "Order of Richardson Extrapolation Polynomial");
        FIELD_MAKE_OPTIONAL(RichardsonPolyOrder);
        FIELD(nbIterMaxP, "nbIterMaxP");
        FIELD_MAKE_TRANSIENT(nbIterMaxP);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultClosedRomberg1D(){
        return new ClosedRomberg1D(0);
   }
};

CClassConstSP const ClosedRomberg1D::TYPE = CClass::registerClassLoadMethod(
    "ClosedRomberg1D", typeid(ClosedRomberg1D), ClosedRomberg1DHelper::load);

// OpenRomberg1D
const double OpenRomberg1D::default_accuracy = 1.0e-6;
const int    OpenRomberg1D::default_nbIterMax = 20;
const int    OpenRomberg1D::default_RichardsonPolyOrder = 9;


OpenRomberg1D::MidpntRule::MidpntRule(const Function1DDouble* infunc, 
                                      double                  a, 
                                      double                  b):
func(infunc),
a(a), 
b(b){}

/** Numerical Recipes in C page 142
    midpnt(float(*func)(float), float a, float b, int n)

    Computes the n-th stage refinement of an extended
    midpoint rule. */
double OpenRomberg1D::MidpntRule::operator()(const int n) const
{
    static string method = "OpenRomberg1D::MidpntRule";

    double x,tnm,sum,del,ddel;
    static double s;
    
    int it,j;
    if (n == 1) {
        return (s=(b-a)*(*func)(0.5*(a+b)));
    } 
    else {
        for(it=1,j=1;j<n-1;j++) it *= 3;
        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;
        for (j=1;j<=it;j++) {
            sum += (*func)(x);
            x += ddel;
            sum += (*func)(x);
            x += del;
        }
 
        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}

OpenRomberg1D::ExpRule::ExpRule(const Function1DDouble* infunc,
                                double aa):
func(infunc),
a(0.0), b(exp(-aa)){}

/** Numerical Recipes in C page 146
    midpnt(float(*func)(float), float aa, float bb, int n)

    Computes the n-th stage refinement of an extended
    midpoint rule. */
double OpenRomberg1D::ExpRule::operator()(const int n) const
{
    static string method = "OpenRomberg1D::ExpRule";

    double x,tnm,sum,del,ddel;
    static double s;
    
    int it,j;
    
    if (n == 1) {
        return (s=(b-a)*(*func)(-log(0.5*(a+b)))/(0.5*(a+b)));
    } 
    else {
        for(it=1,j=1;j<n-1;j++) it *= 3;
        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;
        for (j=1;j<=it;j++) {
            sum += (*func)(-log(x))/x;
            x += ddel;
            sum += (*func)(-log(x))/x;
            x += del;
        }
 
        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}


/** Numerical Recipes in C page 143
    qromo(float(*func)(float), float a, float b,
          float(*choose)(float(*)(float), float, float, int))

    Romberg integration on an open interval. Choose is an
    open formula, not evaluating the function at the end
    points a, b */
double OpenRomberg1D::integrate(const Function1DDouble& func) const {
    static const string method = "OpenRomberg1D::integrate";
    
    try{
        const Range& interval = func.getInterval();

        Range::checkIsNonEmpty(interval);

        if (interval.isSingleton()){
            return 0.0;
        }

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        if (lower.isInfinite() && upper.isInfinite()){
            throw ModelException(method,
                                 "Infinite intervals are not supported. Only semi-infinite intervals are supported; got " + interval.toString());
        }

        int j;
        double ss,dss;

        vector<double> h(nbIterMaxP+1);
        vector<double> s(nbIterMaxP+1);
        h[1]=1.0;

        ChooseConstSP choose;

        try{
            if(!lower.isInfinite() && upper.isInfinite()) {
                choose = ChooseConstSP(new ExpRule(&func, lower.getValue()));
            }
            else if(!lower.isInfinite() && !upper.isInfinite()) {
                choose = ChooseConstSP(new MidpntRule(&func, lower.getValue(), upper.getValue()));
            }
            else {
                throw ModelException(method, interval.toString() + "is not a supported interval");
            }
        }
        catch (exception& e){
            throw ModelException(e, method, "Failed to create open rule for function");
        }
        
        for (j=1;j<=nbIterMax;j++) {
            s[j]=(*choose)(j);
            if (j >= RichardsonPolyOrder) {
                polint(&h[j-RichardsonPolyOrder],
                       &s[j-RichardsonPolyOrder],
                       RichardsonPolyOrder,
                       0.0,
                       &ss,
                       &dss);
		 if (Maths::isNegative(fabs(dss) - accuracy*fabs(ss)) ||
		     (Maths::isZero(dss) && ((fabs(dss) - accuracy*fabs(ss)) <= 0))) {
                    return ss;
                }
            }

            h[j+1]=h[j]/9.0; 
        }
              
        throw ModelException(method,
                             "Number of iterations exceeded " + 
                             Format::toString(nbIterMax) + ".");
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void OpenRomberg1D::validatePop2Object() { /* For initializing nbIterMaxP */
    nbIterMaxP = nbIterMax + 1;
}

OpenRomberg1D::OpenRomberg1D(double accuracy, int nbIterMax, int RichardsonPolyOrder):
CObject(TYPE),
accuracy(accuracy),
nbIterMax(nbIterMax),
RichardsonPolyOrder(RichardsonPolyOrder),
nbIterMaxP(nbIterMax + 1){}

OpenRomberg1D::OpenRomberg1D(int notUsed):
CObject(TYPE),
accuracy(default_accuracy),
nbIterMax(default_nbIterMax),
RichardsonPolyOrder(default_RichardsonPolyOrder),
nbIterMaxP(default_nbIterMax + 1){}

class OpenRomberg1DHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(OpenRomberg1D, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Integrator1D);
        EMPTY_SHELL_METHOD(defaultOpenRomberg1D);
        FIELD(accuracy, "Accuracy");
        FIELD_MAKE_OPTIONAL(accuracy);
        FIELD(nbIterMax, "Maximum Nb of Itertions");
        FIELD_MAKE_OPTIONAL(nbIterMax);
        FIELD(RichardsonPolyOrder, "Order of Richardson Extrapolation Polynomial");
        FIELD_MAKE_OPTIONAL(RichardsonPolyOrder);
        FIELD(nbIterMaxP, "nbIterMaxP");
        FIELD_MAKE_TRANSIENT(nbIterMaxP);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultOpenRomberg1D(){
        return new OpenRomberg1D(0);
    }
};

CClassConstSP const OpenRomberg1D::TYPE = CClass::registerClassLoadMethod(
    "OpenRomberg1D", typeid(OpenRomberg1D), OpenRomberg1DHelper::load);

// FFTIntegrator1D::Integral
void FFTIntegrator1D::Integral::load(CClassSP &clazz){
    REGISTER_INTERFACE(FFTIntegrator1D::Integral, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}

CClassConstSP const FFTIntegrator1D::Integral::TYPE = CClass::registerInterfaceLoadMethod(
    "FFTIntegrator1D::Integral", typeid(FFTIntegrator1D::Integral), FFTIntegrator1D::Integral::load);

class FFTIntegrator1D_Integral: public CObject,
                                public FFTIntegrator1D::Integral{
public:
    static CClassConstSP const TYPE;
    friend class FFTIntegrator1D_IntegralHelper;

    virtual double getValue(double weight) const{
        if (!ppoly){
            doTheWork();
        }
        if(weight < lw || weight > uw) {
            throw ModelException(Format::toString("Weight (%f) is outside interpolation range [%f,%f]",
                                                  weight,
                                                  lw,
                                                  uw));
        }
        return ppoly->value(weight);
    }

    virtual void getValue(const CDoubleArray& weights,
                          CDoubleArray&       values) const{
        if (!ppoly){
            doTheWork();
        }
        int iWeight = 0;
        for (; iWeight < weights.size(); ++iWeight){
            if(weights[iWeight] < lw || weights[iWeight] > uw) {
                throw ModelException(Format::toString("Weight (%f) is outside interpolation range [%f,%f]",
                                                      weights[iWeight],
                                                      lw,
                                                      uw));
            }
        }
        ppoly->value(weights, 
                     0,     // not a derivative
                     values);
    }
    
    FFTIntegrator1D_Integral(int nbPoint,
                             double range,
                             const CDoubleArraySP& p):
    CObject(TYPE),
    nbPoint(nbPoint),
    range(range),
    p(p){}
            
private:
    FFTIntegrator1D_Integral(): CObject(TYPE){}

    void doTheWork() const{
        static const string method = "FFTIntegrator1D::doTheWork";
        try{
            /* Call fft routine */
            DoubleArray q(nbPoint);

            imsl_d_fft_real(nbPoint,
                            &(*p)[0],
                            IMSL_BACKWARD,
                            IMSL_RETURN_USER, &q[0],
                            0);
        
            IMSLError::throwExceptionIfError();

            /* Reorder output and create abscissa */
            CDoubleArray qq(nbPoint);
            CDoubleArray ww(nbPoint);
            int j = 0, i = nbPoint / 2;
            int k = (nbPoint + 1) / 2;
//            double dw = range / nbPoint;
            double dw = 2.0 * Maths::PI / range;
            for (; i < nbPoint; ++i, ++j){
                qq[j] = q[i];
                ww[j] = dw * (j - k);
            }
            for (i = 0; i < nbPoint / 2; ++i, ++j){
                qq[j] = q[i];
                ww[j] = dw * (j - k);
            }
            lw = - dw * k;
            uw = dw * (nbPoint - k);

            /* Create spline interpolant.
               NB: Abscissa do not need to be order for imsl */
            CubicSplineInterpECnd spline(0, 0, 0, 0);
            ppoly = Interpolator::InterpolantConstSP(spline.computeInterp(ww, qq));
        }
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    int nbPoint;
    double range;
    CDoubleArraySP p;
    mutable Interpolator::InterpolantConstSP ppoly; // $unregistered
    mutable double lw; // $unregistered
    mutable double uw; // $unregistered
};

class FFTIntegrator1D_IntegralHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FFTIntegrator1D_Integral, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(nbPoint, "");
        FIELD_MAKE_TRANSIENT(nbPoint);
        FIELD(range, "");
        FIELD_MAKE_TRANSIENT(range);
        FIELD(p, "");
        FIELD_MAKE_TRANSIENT(p);
//        FIELD(ppoly, "");
//        FIELD_MAKE_TRANSIENT(ppoly);
    }

    static IObject* defaultCtor(){
        return new FFTIntegrator1D_Integral();
    }
};

CClassConstSP const FFTIntegrator1D_Integral::TYPE = CClass::registerClassLoadMethod(
    "FFTIntegrator1D_Integral", typeid(FFTIntegrator1D_Integral), FFTIntegrator1D_IntegralHelper::load);

// FFTIntegrator1D
const double FFTIntegrator1D::default_range = 20.0;
const int FFTIntegrator1D::default_nbPoint = 1024;

FFTIntegrator1D::IntegralConstSP FFTIntegrator1D::integrate(const Function1DComplex& func) const{
    static const string method = "FFTIntegrator1D::integrate";
    
    try{
        const Range& interval = func.getInterval();

        Range::checkIsNonEmpty(interval);

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        if (!lower.isInfinite() && !upper.isInfinite()){
            throw ModelException(method,
                                 "Need an infinite interval; got " + interval.toString());
        }

//        double du = 2.0 * Maths::PI / range;
        double du = range / nbPoint;
        double mult = du / (2.0 * Maths::PI);
        CDoubleArraySP p(new DoubleArray(nbPoint));
        DoubleArray& pp = *p;   // for ease

        Complex funcval = func(0.0);
        pp[0] = mult * funcval.real();
        double u = du;
        int iPoint = 1;
        for (; iPoint < nbPoint / 2; ++iPoint, u += du){
            funcval = func(u);
            pp[2 * iPoint - 1] = + mult * funcval.real();
            if (iPoint < nbPoint / 2){
                pp[2 * iPoint] = - mult * funcval.imag();
            }
        }
        return FFTIntegrator1D::IntegralSP(new FFTIntegrator1D_Integral(nbPoint, range, p));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void FFTIntegrator1D::validatePop2Object(){
    static const string method = "FFTIntegrator1D::validatePop2Object";
    
    if (nbPoint <= 0){
        throw ModelException(method, "nbPoint should be greater than or equal to 1");
    }
}

FFTIntegrator1D* FFTIntegrator1D::create(double range, int nbPoint){
    return new FFTIntegrator1D(range, nbPoint);
}

FFTIntegrator1D::FFTIntegrator1D(double range, int nbPoint):
CObject(TYPE),
range(range),
nbPoint(nbPoint){}

class FFTIntegrator1DHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FFTIntegrator1D, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(range, "Range of integration");
        FIELD_MAKE_OPTIONAL(range);
        FIELD(nbPoint, "Nb of integration points");
        FIELD_MAKE_OPTIONAL(nbPoint);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new FFTIntegrator1D();
    }
};

CClassConstSP const FFTIntegrator1D::TYPE = CClass::registerClassLoadMethod(
    "FFTIntegrator1D", typeid(FFTIntegrator1D), FFTIntegrator1DHelper::load);

// IntFuncInf
const double IntFuncInf::default_err_abs = sqrt(DBL_EPSILON);
const double IntFuncInf::default_err_rel = sqrt(DBL_EPSILON);
const int    IntFuncInf::default_max_subinter = 500;

//static const Function1DDouble* local_func = 0;
//static double local_fcn(double x) {
//    return (*local_func)(x);
//}

/** imsl's int_func_inf */
double IntFuncInf::integrate(const Function1DDouble& infunc) const {
    static const string method = "IntFuncInf::integrate";
    
    try{
        const Range& interval = infunc.getInterval();

        Range::checkIsNonEmpty(interval);

        if (interval.isSingleton()){
            return 0.0;
        }

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        Imsl_quad imsl_interval;
        double bound = 0.0; // not used if interval == (-infty, + infty)
        if (lower.isInfinite() && upper.isInfinite()) {
            imsl_interval = IMSL_INF_INF;
        }
        else if (lower.isInfinite()) {
            imsl_interval = IMSL_INF_BOUND;
            bound = upper.getValue();
        }
        else if (upper.isInfinite()) {
            imsl_interval = IMSL_BOUND_INF;
            bound = lower.getValue();
        }
        else {
            throw ModelException(method,
                                 "Need an infinite or semi-infinite interval; got " + interval.toString());
        }

        IMSLFunctionProvider funcProvider(&infunc);
        double err_est;
        int n_subinter, n_evals;
        double rtn = imsl_d_int_fcn_inf (funcProvider.function(), 
                                         bound, 
                                         imsl_interval,
                                         IMSL_ERR_ABS, err_abs,
                                         IMSL_ERR_REL, err_rel,
                                         IMSL_ERR_EST, &err_est,
                                         IMSL_MAX_SUBINTER, max_subinter,
                                         IMSL_N_SUBINTER, &n_subinter,
                                         IMSL_N_EVALS, &n_evals,
                                         0);

        IMSLError::throwExceptionIfError();

        return rtn;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}
   
IntFuncInf* IntFuncInf::create(double err_abs,
                               double err_rel,
                               int max_subinter){
    return new IntFuncInf(err_abs, err_rel, max_subinter);
}

IntFuncInf::IntFuncInf(double err_abs,
                       double err_rel,
                       int max_subinter):
CObject(TYPE),
err_abs(err_abs),
err_rel(err_rel),
max_subinter(max_subinter){}

IntFuncInf::IntFuncInf(int notUsed):
CObject(TYPE),
err_abs(default_err_abs),
err_rel(default_err_rel),
max_subinter(default_max_subinter){}

class IntFuncInfHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IntFuncInf, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Integrator1D);
        EMPTY_SHELL_METHOD(defaultIntFuncInf);
        FIELD(err_abs, "Absolute Error");
        FIELD_MAKE_OPTIONAL(err_abs);
        FIELD(err_rel, "Relative Error");
        FIELD_MAKE_OPTIONAL(err_rel);
        FIELD(max_subinter, "Maximum Nb of Sub-intervals");
        FIELD_MAKE_OPTIONAL(max_subinter);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultIntFuncInf(){
        return new IntFuncInf(0);
    }
};

CClassConstSP const IntFuncInf::TYPE = CClass::registerClassLoadMethod(
    "IntFuncInf", typeid(IntFuncInf), IntFuncInfHelper::load);


// IntFuncClosed
const double IntFuncGeneral::default_err_abs = sqrt(DBL_EPSILON);
const double IntFuncGeneral::default_err_rel = sqrt(DBL_EPSILON);
const int    IntFuncGeneral::default_max_subinter = 500;

/** imsl's int_func_inf */
double IntFuncGeneral::integrate(const Function1DDouble& infunc) const {
    static const string method = "IntFuncGeneral::integrate";
    
    try{
        const Range& interval = infunc.getInterval();

        Range::checkIsNonEmpty(interval);

        if (interval.isSingleton()){
            return 0.0;
        }

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        if (!lower.isClosedBracket() || lower.isInfinite()
            || !upper.isClosedBracket() || upper.isInfinite()) {
            throw ModelException(method,
                                 "(Semi-)Open and / or infinite intervals are not supported; got " + interval.toString());
        }

        IMSLFunctionProvider funcProvider(&infunc);
        double err_est;
        int n_subinter, n_evals;
        double rtn = imsl_d_int_fcn(funcProvider.function(), 
                                    lower.getValue(),
                                    upper.getValue(),
                                    IMSL_ERR_ABS, err_abs,
                                    IMSL_ERR_REL, err_rel,
                                    IMSL_ERR_EST, &err_est,
                                    IMSL_MAX_SUBINTER, max_subinter,
                                    IMSL_N_SUBINTER, &n_subinter,
                                    IMSL_N_EVALS, &n_evals,
                                    0);

        IMSLError::throwExceptionIfError();

        return rtn;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

IntFuncGeneral* IntFuncGeneral::create(double err_abs,
                                       double err_rel,
                                       int max_subinter){
    return new IntFuncGeneral(err_abs, err_rel, max_subinter);
}

IntFuncGeneral::IntFuncGeneral(double err_abs,
                               double err_rel,
                               int max_subinter):
CObject(TYPE),
err_abs(err_abs),
err_rel(err_rel),
max_subinter(max_subinter){}

IntFuncGeneral::IntFuncGeneral(int notUsed):
CObject(TYPE),
err_abs(default_err_abs),
err_rel(default_err_rel),
max_subinter(default_max_subinter){}

class IntFuncGeneralHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IntFuncGeneral, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Integrator1D);
        EMPTY_SHELL_METHOD(defaultIntFuncGeneral);
        FIELD(err_abs, "Absolute Error");
        FIELD_MAKE_OPTIONAL(err_abs);
        FIELD(err_rel, "Relative Error");
        FIELD_MAKE_OPTIONAL(err_rel);
        FIELD(max_subinter, "Maximum Nb of Sub-intervals");
        FIELD_MAKE_OPTIONAL(max_subinter);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultIntFuncGeneral(){
        return new IntFuncGeneral(0);
    }
};

CClassConstSP const IntFuncGeneral::TYPE = CClass::registerClassLoadMethod(
    "IntFuncGeneral", typeid(IntFuncGeneral), IntFuncGeneralHelper::load);

// force linker to include this file (avoid having header file) */
bool IntegratorLoad() {
    return (Trapez1D::TYPE != 0
            && ClosedRomberg1D::TYPE != 0
            && OpenRomberg1D::TYPE != 0
            && FFTIntegrator1D::TYPE != 0
            && IntFuncInf::TYPE != 0
            && IntFuncGeneral::TYPE != 0);
}

/*********************************************************************/

/* The ultimate wrapping of Integrator1D, mainly for use in Pyramid 
 */
#define INTEGRATOR1D_TYPE_CLOSEDROMBERG         "ClosedRomberg"
#define INTEGRATOR1D_TYPE_OPENROMBERG           "OpenRomberg"
#define INTEGRATOR1D_TYPE_INTFUNCGENERAL        "IntFuncGeneral"
#define INTEGRATOR1D_TYPE_INTFUNCINF            "IntFuncInf"

class Integrator1DWrapper : public CObject,
                            virtual public ITypeConvert {
public: // how can I have this protected or private?
    string                  integratorType;     // ClosedRomberg, OpenRomberg, IntFuncGeneral, or IntFuncInf
    ClosedRomberg1DSP       closedRomberg;
    OpenRomberg1DSP         openRomberg;
    IntFuncGeneralSP        intFuncGeneral;
    IntFuncInfSP            intFuncInf;

private:
    Integrator1DSP          realIntegrator;

public:
    static CClassConstSP const TYPE;

    // validation
    void validatePop2Object(){
        static const string routine = 
            "Integrator1DWrapper::validatePop2Object";
        try{
            if (integratorType.empty()){
                throw ModelException(routine, 
                                     "Blank integrator specified!");
            }
            if (CString::equalsIgnoreCase(integratorType,INTEGRATOR1D_TYPE_CLOSEDROMBERG)) {
                if (closedRomberg.get()) {
                    realIntegrator = closedRomberg;
                } else {
                    throw ModelException(routine, "Expected ClosedRomberg1D "
                                         "but none supplied!");
                }
            } else if (CString::equalsIgnoreCase(integratorType,INTEGRATOR1D_TYPE_OPENROMBERG)) {
                if (openRomberg.get()) {
                    realIntegrator = openRomberg;
                } else {
                    throw ModelException(routine, "Expected OpenRomberg1D "
                                         "but none supplied!");
                }
            } else if (CString::equalsIgnoreCase(integratorType,INTEGRATOR1D_TYPE_INTFUNCGENERAL)) {
                if (intFuncGeneral.get()) {
                    realIntegrator = intFuncGeneral;
                } else {
                    throw ModelException(routine, "Expected IntFuncGeneral "
                                         "but none supplied!");
                }
            } else if (CString::equalsIgnoreCase(integratorType,INTEGRATOR1D_TYPE_INTFUNCINF)) {
                if (intFuncInf.get()) {
                    realIntegrator = intFuncInf;
                } else {
                    throw ModelException(routine, "Expected IntFuncInf "
                                         "but none supplied!");
                }
            } else {
                throw ModelException(routine, "Unrecognised Intergrator1D "
                                     + integratorType + ". Expected " 
                                     + INTEGRATOR1D_TYPE_CLOSEDROMBERG + ", " 
                                     + INTEGRATOR1D_TYPE_OPENROMBERG + ", " 
                                     + INTEGRATOR1D_TYPE_INTFUNCGENERAL + " or " 
                                     + INTEGRATOR1D_TYPE_INTFUNCINF);
            }
        }  catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** create a proper Integrator1D */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const {
        static const string method = "Integrator1DWrapper::convert";
        try {
            if (requiredType != Integrator1D::TYPE) {
                throw ModelException(method, 
                                     "Cannot convert a Integrator1D into "
                                     "object of type "+requiredType->getName());               
            }
            object = realIntegrator;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Integrator1DWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITypeConvert);
        EMPTY_SHELL_METHOD(defaultIntegrator1DWrapper);
        FIELD(integratorType, "ClosedRomberg, OpenRomberg, IntFuncGeneral or IntFuncInf");
        FIELD(closedRomberg,  "ClosedRomberg1D");
        FIELD_MAKE_OPTIONAL(closedRomberg);
        FIELD(openRomberg,  "OpenRomberg1D");
        FIELD_MAKE_OPTIONAL(openRomberg);
        FIELD(intFuncGeneral,  "IntFuncGeneral");
        FIELD_MAKE_OPTIONAL(intFuncGeneral);
        FIELD(intFuncInf,  "IntFuncInf");
        FIELD_MAKE_OPTIONAL(intFuncInf);
        FIELD(realIntegrator, "real Integrator1D");
        FIELD_MAKE_TRANSIENT(realIntegrator);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    Integrator1DWrapper(): CObject(TYPE){}

    static IObject* defaultIntegrator1DWrapper(){
        return new Integrator1DWrapper();
    }
};

typedef smartPtr<Integrator1DWrapper> Integrator1DWrapperSP;

CClassConstSP const Integrator1DWrapper::TYPE =
CClass::registerClassLoadMethod("Integrator1DWrapper", 
                                typeid(Integrator1DWrapper), load);

DRLIB_END_NAMESPACE
