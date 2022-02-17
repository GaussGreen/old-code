//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LinearInterpolator.hpp
//
//   Description : 
//
//   Date        : 06 June 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Format.hpp"
#include "edginc/smartPtr.hpp"
#include <math.h>
#include "edginc/AtomicArray.hpp"
#include "edginc/Interpolator.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Nrfns.hpp"

#ifndef EDR_LINEARINTERPOLATOR_HPP
#define EDR_LINEARINTERPOLATOR_HPP

DRLIB_BEGIN_NAMESPACE

UTIL_DLL CDoubleArray computeSlopes(const CDoubleArray& xx,
                                    const CDoubleArray& yy);

/** Equidistant interpolant */
template<class _TBase = Interpolator::Interpolant>
class EquidistantInterpolantMain: public InterpolantBase<_TBase> {
public:
    static CClassConstSP const TYPE;
    
    /** Constructor from xx and yy arrays and step size */
    EquidistantInterpolantMain(const CDoubleArray& xx, const CDoubleArray& yy, double step):
    InterpolantBase<_TBase>(TYPE, xx, yy), step(step) {
        static const string routine("EquidistantInterpolantMain::EquidistantInterpolantMain");

        try {
            validatePop2Object();
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Default constructor */
    EquidistantInterpolantMain(): InterpolantBase<_TBase>(TYPE) {}
    
    double value(double xx) const {
        if(xx <= xLow) {
            return this->y[0];
        }

        if(xx >= xHigh) {
            return this->y[this->n-1];
        }

        int position = lookupValue(xx);
        return this->y[position] + dydx[position] * (xx-this->x[position]);
    }

    double value(double xx,
                 int    deriv) const {
        if(deriv == 0) {
            return value(xx);        // y(x)
        } 
    
        if(deriv == 1) {
            // y'(x) = dydx
            if(xx <= xLow || xx >= xHigh) {
                return 0.0;
            }

            int position = lookupValue(xx);
            return dydx[position];
        }
    
        // zero higher derivatives because it's linear
        return 0.0;
    }

    void value(const CDoubleArray& xvec,           
               int                 deriv,
               CDoubleArray&       valuevec) const {
        static const string method("LinearInterpolantMain::value");
        if(xvec.size() != valuevec.size()) {
            throw ModelException(method, "Sizes of x and y arrays are not equal: " +
                                          Format::toString(xvec.size()) + ", " + 
                                          Format::toString(valuevec.size()));

    
        }
        int i;
        for(i = 0; i < xvec.size(); i++) {
            valuevec[i] = value(xvec[i], deriv);
        }
    }

    double invValue(double yy) const {
        static unsigned long position;
        // NR searches in an array xx[1,...,n]
        locate(&this->y[0]-1, this->n, yy, &position); // position is in [0, n]
        
        if(position == 0) {
            return this->x[0];
        } 
    
        if(int(position) == this->n) {
            return this->x[this->n-1];
        }
    
        --position;
        return this->x[position] + (yy-this->y[position]) / dydx[position];
    }

    virtual void validatePop2Object() {
        static const string method = "EquidistantInterpolantMain::validatePop2Object";

        try {
            // Call parent method
            InterpolantBase<_TBase>::validatePop2Object();

            invStep = 1.0 / step;
            xLow = this->x.front();
            xHigh = this->x.back();
            
            // Compute slopes           
            dydx = computeSlopes(this->x, this->y);

        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    /** Optimal search for equidistant interpolant */
    inline int lookupValue(double xx) const {
        return static_cast<int>((xx - xLow) * invStep);
    }
    
    double step;            //!< Step size for equidistant XArray
    double invStep;         //!< Inverse step used in lookUpValue
    double xLow;            //!< Lowest x value. Used for faster access
    double xHigh;           //!< Highest x value. Used for faster access
    bool increasing;        //!< Whether XArray is increasing or not $unregistered
    CDoubleArray dydx;      //!< Slopes

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(EquidistantInterpolantMain<_TBase>, clazz);
        SUPERCLASS(InterpolantBase<_TBase>);
        EMPTY_SHELL_METHOD(defaultEquidistantInterpolantMain);
        FIELD(step, "Step size for XArray");
        FIELD(invStep, "Inverse step size");
        FIELD_MAKE_TRANSIENT(invStep);
        FIELD(xLow, "Lowest x value");
        FIELD_MAKE_TRANSIENT(xLow);
        FIELD(xHigh, "Highest x value");
        FIELD_MAKE_TRANSIENT(xHigh);
        FIELD(dydx, "slopes");
        FIELD_MAKE_TRANSIENT(dydx);
    }

    static IObject* defaultEquidistantInterpolantMain() {
        return new EquidistantInterpolantMain<_TBase>();
    }
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each array template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class _TBase> CClassConstSP const 
EquidistantInterpolantMain<_TBase>::TYPE =
CClass::templateRegisterClass(typeid(EquidistantInterpolantMain<_TBase>));
#endif

typedef EquidistantInterpolantMain<Interpolator::Interpolant> EquidistantLinearInterpolant;
typedef smartPtr<EquidistantLinearInterpolant> EquidistantLinearInterpolantSP;
typedef smartConstPtr<EquidistantLinearInterpolant> EquidistantLinearInterpolantConstSP;
#ifndef QLIB_LINEARINTERPOLATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<EquidistantLinearInterpolant>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<EquidistantLinearInterpolant>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<EquidistantLinearInterpolant>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<EquidistantLinearInterpolant>);
#endif

typedef EquidistantInterpolantMain<Interpolator::NullBase> EquidistantLinearInterpolantNonVirtual;
typedef smartPtr<EquidistantLinearInterpolantNonVirtual> EquidistantLinearInterpolantNonVirtualSP;
typedef smartConstPtr<EquidistantLinearInterpolantNonVirtual> EquidistantLinearInterpolantNonVirtualConstSP;
#ifndef QLIB_LINEARINTERPOLATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<EquidistantLinearInterpolantNonVirtual>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<EquidistantLinearInterpolantNonVirtual>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<EquidistantLinearInterpolantNonVirtual>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<EquidistantLinearInterpolantNonVirtual>);
#endif

////////////////////////////////////////////////////////////////////////////////

/** Linear interpolant */
template<class _TBase = Interpolator::Interpolant>
class LinearInterpolantMain: public InterpolantBase<_TBase> {
public:
    static CClassConstSP const TYPE;
    
    /** Constructor from xx and yy arrays */
    LinearInterpolantMain(const CDoubleArray& xx, const CDoubleArray& yy):
    InterpolantBase<_TBase>(TYPE, xx, yy) {
        static const string routine("EquidistantInterpolantMain::EquidistantInterpolantMain");
        
        try {
            validatePop2Object();
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    LinearInterpolantMain(): InterpolantBase<_TBase>(TYPE) {}
    
    double value(double xx) const {
        int position = lookupValue(xx, this->x);
        if(position == 0) {
            return this->y[0];
        } 
    
        if(position == this->n) {
            return this->y[this->n-1];
        }
    
        --position;
        return this->y[position] + dydx[position] * (xx-this->x[position]);
    }

    double value(double xx,
                 int    deriv) const {
        if(deriv == 0) {
            return value(xx);        // y(x)
        } 
    
        if(deriv == 1) {
            int position = lookupValue(xx, this->x);
            // y'(x) = dydx
            if(position == 0 || position == this->n) {
                return 0.0;
            }

            return dydx[position-1];
        }
    
        // zero higher derivatives because it's linear
        return 0.0;
    }

    void value(const CDoubleArray& xvec,           
               int                 deriv,
               CDoubleArray&       valuevec) const {
        static const string method("LinearInterpolantMain::value");
        if(xvec.size() != valuevec.size()) {
            throw ModelException(method, "Sizes of x and y arrays are not equal: " +
                                          Format::toString(xvec.size()) + ", " + 
                                          Format::toString(valuevec.size()));

    
        }
        int i;
        for(i = 0; i < xvec.size(); i++) {
            valuevec[i] = value(xvec[i], deriv);
        }
    }

    double valueWithGuess(double xx, 
                          int&   position) const {
        
        position = lookupValueWithGuess(xx, this->x, position);

        if(position == 0) {
            return this->y[0];
        } 
    
        if(position == this->n) {
            return this->y[this->n-1];
        } 

        // Do not do --position as position is returned to user as an updated guess
        int pos = position - 1;
        return this->y[pos] + dydx[pos] * (xx-this->x[pos]);
    }

    double valueWithGuess(double xx,
                          int    deriv,
                          int&   position) const {
        if(deriv == 0) {
            return valueWithGuess(xx, position);        // y(x)
        } 
    
        if(deriv == 1) {
            position = lookupValueWithGuess(xx, this->x, position);
            // y'(x) = dydx
            if (position == 0 || position == this->n) {
                return 0.0;
            } 
            return dydx[position-1];
        }

        // zero higher derivatives because it's linear
        return 0.0;
    }

    double invValue(double yy) const {
        int position = lookupValue(yy, this->y);
        if(position == 0) {
            return this->x[0];
        } 
    
        if(position == this->n) {
            return this->x[this->n-1];
        }
    
        --position;
        return this->x[position] + (yy-this->y[position]) / dydx[position];
    }

    double invValue(double yy,
                    int    deriv) const {
        if(deriv == 0) {
            return invValue(yy);        // x(y)
        } 
    
        if(deriv == 1) {
            int position = lookupValue(yy, this->y);
            // flat extrapolation on the inverse as well
            if(position == 0 || position == this->n) {
                return 0.0;
            }
            return 1.0 / dydx[position-1];
        }
    
        return 0.0;
    }

    void invValue(const CDoubleArray& valuevec,           
                  int                 deriv,
                  CDoubleArray&       xvec) const {
        static const string method("LinearInterpolantMain::invValue");
        if(xvec.size() != this->yvec.size()) {
            throw ModelException(method, "Sizes of x and y arrays are not equal: " +
                                          Format::toString(xvec.size()) + ", " + 
                                          Format::toString(this->yvec.size()));

        }
        int i;
        for(i = 0; i < xvec.size(); i++) {
            xvec[i] = invValue(this->yvec[i], deriv);
        }
    }

    virtual void validatePop2Object() {
        static const string method = "LinearInterpolantMain::validatePop2Object";

        try {
            // Compute slopes           
            dydx = computeSlopes(this->x, this->y);

        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    void editY(int idx, double newY) 
	{
	    this->y[idx] = newY;
	    // Recompute slopes           
        dydx = computeSlopes(this->x, this->y);
	}


private:
    CDoubleArray dydx;    //!< Slopes

    /** Searches an ordered table using NR locate */
    inline int lookupValue(double zz, CDoubleArray& z) const {
        static unsigned long position;
        // NR searches in an array xx[1,...,n]
        locate(&z[0]-1, this->n, zz, &position); // position is in [0, n]
        return position;
    }

    /** Searches an ordered table with a guess using NR hunt */
    inline int lookupValueWithGuess(double zz, CDoubleArray& z, int i) const {
        /* static const string method("LinearInterpolant::lookupValueWithGuess");
        What is the speed effect of this?
        if(i<0 || i>n) {
            throw ModelException(method, "Index " + Format::toString(i) + 
                                         " out of range [0," + 
                                         Format::toString(n) + "]");
        }*/

        unsigned long position = i;
        // NR searches in an array xx[1,...,n]
        hunt(&z[0]-1, this->n, zz, &position); // position is in [0, n]
        return position;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LinearInterpolantMain<_TBase>, clazz);
        SUPERCLASS(InterpolantBase<_TBase>);
        EMPTY_SHELL_METHOD(defaultLinearInterpolantMain);
        FIELD(dydx, "slopes");
        FIELD_MAKE_TRANSIENT(dydx);
    }

    static IObject* defaultLinearInterpolantMain() {
        return new LinearInterpolantMain<_TBase>();
    }
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each array template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class _TBase> CClassConstSP const 
LinearInterpolantMain<_TBase>::TYPE =
CClass::templateRegisterClass(typeid(LinearInterpolantMain<_TBase>));
#endif

// Some typedefs
typedef LinearInterpolantMain<Interpolator::Interpolant> LinearInterpolant;
typedef smartPtr<LinearInterpolant> LinearInterpolantSP;
typedef smartConstPtr<LinearInterpolant> LinearInterpolantConstSP;
#ifndef QLIB_LINEARINTERPOLATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<LinearInterpolant>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<LinearInterpolant>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<LinearInterpolant>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<LinearInterpolant>);
#endif

typedef LinearInterpolantMain<Interpolator::NullBase> LinearInterpolantNonVirtual;
typedef smartPtr<LinearInterpolantNonVirtual> LinearInterpolantNonVirtualSP;
typedef smartConstPtr<LinearInterpolantNonVirtual> LinearInterpolantNonVirtualConstSP;
#ifndef QLIB_LINEARINTERPOLATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<LinearInterpolantNonVirtual>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<LinearInterpolantNonVirtual>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<LinearInterpolantNonVirtual>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<LinearInterpolantNonVirtual>);
#endif

////////////////////////////////////////////////////////////////////////////////

/** Produces a linear interpolant. Note due to templated functions within class
    cannot export entire class. */
class LinearInterpolator: public Interpolator, public CObject {
public:
    static UTIL_DLL CClassConstSP const TYPE;
    friend class LinearInterpolatorHelper;
    
    /** Default constructor */
    UTIL_DLL LinearInterpolator();

    /** Produces a linear interpolant */
    virtual UTIL_DLL Interpolator::InterpolantConstSP computeInterp(const CDoubleArray& xdata,
                                                                    const CDoubleArray& fdata) const;

    /** Produces a non virtual casted linear interpolant */
    static UTIL_DLL LinearInterpolantNonVirtualConstSP computeLinearInterpNV(const LinearInterpolator& interpolator,
                                                                             const CDoubleArray& xdata,
                                                                             const CDoubleArray& fdata);
    
    /** Produces an equidistant linear interpolant with virtual base. Func must support () */
    template<class _TFunc>
    Interpolator::InterpolantConstSP computeInterp(double  xlow,
                                                   double  xhigh,
                                                   int     nbSteps,
                                                   _TFunc        func) const {
        static string routine = "LinearInterpolator::computeInterp";

        try {
            DoubleArray x, y;
            double step;
            populateData(xlow, xhigh, nbSteps, func, x, y, step);
            return EquidistantLinearInterpolantConstSP(new EquidistantLinearInterpolant(x, y, step));
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Produces a non virtual equidistant linear interpolant. Func must support () */
    template<class _TFunc> static 
    EquidistantLinearInterpolantNonVirtualConstSP computeEquidInterpNV(const LinearInterpolator& interpolator,
                                                                       double  xlow,
                                                                       double  xhigh,
                                                                       int     nbSteps,
                                                                       _TFunc  func) {
        static string routine = "LinearInterpolator::computeEquidInterpNV";

        try {
            DoubleArray x, y;
            double step;
            populateData(xlow, xhigh, nbSteps, func, x, y, step);
        
            return EquidistantLinearInterpolantNonVirtualConstSP(new EquidistantLinearInterpolantNonVirtual(x, y, step));
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

private:
    /** Helper function that populates x and y arrays and step size for an equidistant interpolant */
    template<class _TFunc> static
    void populateData(double xlow,
                      double xhigh,
                      int    nbSteps,
                      _TFunc func,
                      DoubleArray&  x,
                      DoubleArray&  y,
                      double&       stepSize) {
        static string routine = "LinearInterpolator::populateData";
        
        try {
            if(Maths::isNegative(xhigh - xlow)){
                throw ModelException("Low value " + 
                                     Format::toString(xlow) +
                                     " exceeds high value" +
                                     Format::toString(xhigh));
            }

            if(nbSteps <= 0) {
                throw ModelException("Number of steps must be greater than 0.");
            }

            x.resize(nbSteps);
            y.resize(nbSteps);

            stepSize = nbSteps <= 1? 0.0 : (xhigh - xlow) / (double)(nbSteps - 1);
        
            x[0] = xlow;
            y[0] = func(xlow);
            for(int i = 1; i < nbSteps; i++) {
                x[i] = x[i-1] + stepSize;
                y[i] = func(x[i]);
            }
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }
};

typedef smartPtr<LinearInterpolator> LinearInterpolatorSP;
typedef smartConstPtr<LinearInterpolator> LinearInterpolatorConstSP;    
#ifndef QLIB_LINEARINTERPOLATOR_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<LinearInterpolator>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<LinearInterpolator>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<LinearInterpolator>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<LinearInterpolator>);
#endif

DRLIB_END_NAMESPACE

#endif
