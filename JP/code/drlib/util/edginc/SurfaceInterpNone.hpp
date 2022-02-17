//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterpNone.hpp
//
//   Description : 'None' interpolator for the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Dec-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SurfaceInterpNone_HPP
#define QLIB_SurfaceInterpNone_HPP
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/SurfaceInterp.hpp"

DRLIB_BEGIN_NAMESPACE

/** None interpolation class.
 */
class UTIL_DLL SurfaceInterpNone : public SurfaceInterp
{
public:
    static CClassConstSP const TYPE;

    /** Default Constructor */
    SurfaceInterpNone();

    /** Constructor */
    SurfaceInterpNone(
        const DoubleArraySP    xArray,
        const DoubleArraySP    yArray,
        const CDoubleMatrixSP  zMatrix,
        ExtrapMethod           extrapMethod = SurfaceInterpNone::EXTRAP_FLAT);

    /** Copy constructor */
    SurfaceInterpNone(const SurfaceInterpNone& rhs);

    /** Destructor */
    virtual ~SurfaceInterpNone();

    /** Get interpolation method */
    virtual const SurfaceInterp::InterpMethod interpMethod() const { return SurfaceInterp::INTERP_NONE; }

    /** Interpolation functions */
    virtual double        operator()(double x, double y) const;
    virtual double        interpolate(double x, double y) const;

    /** Interpolation functions - OVERRIDE THESE IF EFFICIENT VERSIONS REQUIRED*/
    // virtual DoubleArraySP interpolateX(double x) const;
    // virtual DoubleArraySP interpolateY(double y) const;

protected:

private:

    static void load(CClassSP& clazz);
};


// Support for smart pointers, arrays etc
typedef smartConstPtr<SurfaceInterpNone> SurfaceInterpNoneConstSP;
typedef smartPtr<SurfaceInterpNone> SurfaceInterpNoneSP;

DRLIB_END_NAMESPACE

#endif
