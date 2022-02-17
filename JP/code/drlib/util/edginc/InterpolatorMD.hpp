//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InterpolatorMD.hpp
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
#include "edginc/ArrayMD.hpp"

#ifndef EDR_INTERPOLATORMD_HPP
#define EDR_INTERPOLATORMD_HPP

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL InterpolatorMD: virtual public IObject{
public:
    static CClassConstSP const TYPE;

    class UTIL_DLL Interpolant: public CObject{
    public:
        static CClassConstSP const TYPE;
        
        virtual double value(const DoubleArray& x) const = 0;

        // for convenience
        const DoubleArrayArray& getXarray();
        const DoubleArrayMD& getYarray();
        const DoubleArrayMD::Size& size();

        virtual void validatePop2Object();

    protected:
        DoubleArrayArraySP  x;
        DoubleArrayMDSP     y;
        DoubleArrayMD::Size y_size;
        int                 dim;
        
        Interpolant(const CClassConstSP&    clazz);
        Interpolant(const CClassConstSP&    clazz,
                    const DoubleArrayArray& x,
                    const DoubleArrayMD&    y);
        Interpolant(const CClassConstSP&      clazz,
                    const DoubleArrayArraySP& x,
                    const DoubleArrayMDSP&    y);

    private:
        static void load(CClassSP& clazz);

        Interpolant();
    };

    typedef smartPtr<Interpolant> InterpolantSP;
    typedef smartConstPtr<Interpolant> InterpolantConstSP;

    // data will be deep-copied
    virtual InterpolantConstSP computeInterp(const DoubleArrayArray& xdata,
                                             const DoubleArrayMD&    fdata) const = 0;

    // data will not be deep-copied
    virtual InterpolantConstSP computeInterp(const DoubleArrayArraySP& xdata,
                                             const DoubleArrayMDSP&    fdata) const = 0;

    typedef array<InterpolantSP, Interpolant> InterpolantArray;
    typedef smartPtr<InterpolantArray> InterpolantArraySP;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<InterpolatorMD> InterpolatorMDSP;
typedef smartConstPtr<InterpolatorMD> InterpolatorMDConstSP;
#ifndef QLIB_INTERPOLATORMD_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<InterpolatorMD>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<InterpolatorMD>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<InterpolatorMD>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<InterpolatorMD>);
#endif

DRLIB_END_NAMESPACE

#endif
