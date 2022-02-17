//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LinearInterpolatorMD.hpp
//
//   Description : 
//
//   Date        : 08 Oct 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Format.hpp"
#include "edginc/smartPtr.hpp"
#include <math.h>
#include "edginc/AtomicArray.hpp"
#include "edginc/InterpolatorMD.hpp"
#include "edginc/Maths.hpp"

#ifndef EDR_LINEARINTERPOLATORMD_HPP
#define EDR_LINEARINTERPOLATORMD_HPP

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL LinearInterpolatorMD: public InterpolatorMD, public CObject {
public:
    static CClassConstSP const TYPE;
    friend class LinearInterpolatorMDHelper;

    LinearInterpolatorMD();

    virtual InterpolantConstSP computeInterp(const DoubleArrayArray& xdata,
                                             const DoubleArrayMD&    fdata) const;

    virtual InterpolantConstSP computeInterp(const DoubleArrayArraySP& xdata,
                                             const DoubleArrayMDSP&    fdata) const;

    class UTIL_DLL Interpolant: public InterpolatorMD::Interpolant{
    public:
        static CClassConstSP const TYPE;
        friend class LinearInterpolatorMD_InterpolantHelper;
    
        Interpolant();
        Interpolant(const DoubleArrayArray& x,
                    const DoubleArrayMD&    y);
        Interpolant(const DoubleArrayArraySP& x,
                    const DoubleArrayMDSP&    y);

        virtual double value(const DoubleArray& x) const;

        double valueWithGuess(const DoubleArray&    xx,
                              DoubleArrayMD::Index& index) const;   // index will be modified

    private:
        void lookupValue(const DoubleArray&    xx,
                         DoubleArrayMD::Index& index) const;
        void lookupValueWithGuess(const DoubleArray&    xx,
                                  DoubleArrayMD::Index& index) const;
        double valueRecurse(const DoubleArray&          thex,
                            const DoubleArrayMD::Index& index,
                            int                         i) const;

    };
    typedef smartPtr<Interpolant> InterpolantSP;
    typedef smartConstPtr<Interpolant> InterpolantConstSP;

    InterpolantConstSP computeLinearInterp(const DoubleArrayArray& xdata,
                                           const DoubleArrayMD&    fdata) const;

    InterpolantConstSP computeLinearInterp(const DoubleArrayArraySP& xdata,
                                           const DoubleArrayMDSP&    fdata) const;
};

typedef smartPtr<LinearInterpolatorMD> LinearInterpolatorMDSP;
typedef smartConstPtr<LinearInterpolatorMD> LinearInterpolatorMDConstSP;
#ifndef QLIB_LINEARINTERPOLATORMD_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<LinearInterpolatorMD>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<LinearInterpolatorMD>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<LinearInterpolatorMD>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<LinearInterpolatorMD>);
#endif

DRLIB_END_NAMESPACE
#endif
