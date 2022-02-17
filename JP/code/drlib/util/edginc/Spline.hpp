//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Spline.hpp
//
//   Description : 
//
//   Date        : 05 April 02
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/smartPtr.hpp"
#include <math.h>
#include "edginc/AtomicArray.hpp"
#include "edginc/Interpolator.hpp"

#ifndef EDR_SPLINE_HPP
#define EDR_SPLINE_HPP

DRLIB_BEGIN_NAMESPACE

/** Wrapper around imsl cubic spline */
class UTIL_DLL CubicSplineInterpECnd: public CObject,
                             public Interpolator {
public:
    static CClassConstSP const TYPE;
    friend class CubicSplineInterpECndHelper;

    CubicSplineInterpECnd();

    CubicSplineInterpECnd(int ileft, double left,
                          int iright, double right);

    virtual void validatePop2Object();

    virtual Interpolator::InterpolantConstSP computeInterp(const CDoubleArray& xdata,
                                                           const CDoubleArray& fdata) const;

    virtual Interpolator::InterpolantConstSP computeInterp(const double* xdata,
                                                           const double* fdata,
                                                           int           ndata) const;

private:
    int ileft;
    double left;
    int iright; 
    double right;
};

/** Wrapper around imsl shape preserving cubic spline */
class UTIL_DLL CubicShapePresSplineInterpECnd: public CObject,
                                      public Interpolator {
public:
    static CClassConstSP const TYPE;
    friend class CubicShapePresSplineInterpECndHelper;

    CubicShapePresSplineInterpECnd(bool isConcave);

    virtual Interpolator::InterpolantConstSP computeInterp(const CDoubleArray& xdata,
                                                           const CDoubleArray& fdata) const;

    virtual Interpolator::InterpolantConstSP computeInterp(const double* xdata,
                                                           const double* fdata,
                                                           int           ndata) const;
private:
    CubicShapePresSplineInterpECnd();

    // registered var
    bool isConcave;
};


/** Wrapper around numerical recipes cubic spline. 
    Modified to allow computation of derivatives
    and to allow fast look up with initial guess */
class UTIL_DLL NRSpline: public CObject,
                public Interpolator {
public:
    static CClassConstSP const TYPE;
    friend class NRSplineHelper;

    class UTIL_DLL Interpolant: public InterpolantVirtual {
    public:
        static CClassConstSP const TYPE;
        friend class NRSpline_InterpolantHelper;
        friend class NRSpline;

        virtual double value(double xx) const;

        virtual double value(double xx,
                             int    deriv) const;

        virtual void value(const CDoubleArray& xvec,           
                           int                 deriv,
                           CDoubleArray&       valuevec) const;

        double valueWithGuess(double xx,
                              int&   guess,
                              int    deriv) const;

    private:
        Interpolant();
        Interpolant(const DoubleArray&  xx,
                    const DoubleArray&  yy);
        Interpolant(const double* xx,
                    const double* yy,
                    int           nn);

        // transient
        DoubleArray y2;
    };
    typedef smartPtr<Interpolant> InterpolantSP;
    typedef smartConstPtr<Interpolant> InterpolantConstSP;
    typedef array<InterpolantSP, Interpolant> InterpolantArray;
    
    NRSpline(int ileft, double left,
             int iright, double right);

    virtual void validatePop2Object();

    virtual Interpolator::InterpolantConstSP computeInterp(const CDoubleArray& xdata,
                                                           const CDoubleArray& fdata) const;

    virtual Interpolator::InterpolantConstSP computeInterp(const double* xdata,
                                                           const double* fdata,
                                                           int           ndata) const;

    static InterpolantConstSP computeInterp(const NRSpline&     spliner,
                                            const CDoubleArray& xdata,
                                            const CDoubleArray& fdata);

    static InterpolantConstSP computeInterp(const NRSpline& spliner,
                                            const double*   xdata,
                                            const double* fdata,
                                            int           ndata);

private:
    NRSpline();

    int ileft;
    double left;
    int iright; 
    double right;

    // transient
    double ypleft;
    double ypright;
};

DRLIB_END_NAMESPACE

#endif




