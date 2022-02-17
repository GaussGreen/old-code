//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Spline2D.hpp
//
//   Description : 
//
//   Date        : 24 Oct 2002
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Format.hpp"
#include <math.h>
#include "edginc/InterpolatorMD.hpp"
#include "edginc/Maths.hpp"
#include "edginc/imslerror.hpp"

#ifndef EDR_SPLINE2D_HPP
#define EDR_SPLINE2D_HPP

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL Spline2D: public InterpolatorMD, public CObject {
public:
    static CClassConstSP const TYPE;
    friend class Spline2DHelper;

    Spline2D();

    virtual InterpolantConstSP computeInterp(const DoubleArrayArray& xdata,
                                             const DoubleArrayMD&    fdata) const;

    virtual InterpolantConstSP computeInterp(const DoubleArrayArraySP& xdata,
                                             const DoubleArrayMDSP&    fdata) const;

    class UTIL_DLL Interpolant: public InterpolatorMD::Interpolant{
    public:
        static CClassConstSP const TYPE;
        friend class Spline2D_InterpolantHelper;
    
        Interpolant();
        Interpolant(const DoubleArrayArray& x,
                    const DoubleArrayMD&    y);
        Interpolant(const DoubleArrayArraySP& x,
                    const DoubleArrayMDSP&    y);

        virtual double value(const DoubleArray& x) const;

        virtual void validatePop2Object();

        virtual IObject* clone() const;

    private:
        void createSpline();

        struct IMSLSpline{
            Imsl_d_spline* sp;

            IMSLSpline(Imsl_d_spline* sp):
            sp(sp){}

            ~IMSLSpline(){
                free(sp);
            }
        };
        typedef refCountPtr<IMSLSpline> IMSLSplineSP;
        IMSLSplineSP spline; // $unregistered
    };
    typedef smartPtr<Interpolant> InterpolantSP;
    typedef smartConstPtr<Interpolant> InterpolantConstSP;
};

typedef smartPtr<Spline2D> Spline2DSP;
typedef smartConstPtr<Spline2D> Spline2DConstSP;
#ifndef QLIB_SPLINE2D_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<Spline2D>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<Spline2D>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<Spline2D>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<Spline2D>);
#endif

DRLIB_END_NAMESPACE
#endif
