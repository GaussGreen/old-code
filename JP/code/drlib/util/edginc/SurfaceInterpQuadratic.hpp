//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterpQuadratic.hpp
//
//   Description : Quadratic interpolator for the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Dec-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SurfaceInterpQuadratic_HPP
#define QLIB_SurfaceInterpQuadratic_HPP
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/SurfaceInterp.hpp"

DRLIB_BEGIN_NAMESPACE

/** Quadratic interpolation class.
 */
class UTIL_DLL SurfaceInterpQuadratic : public SurfaceInterp
{
public:
    static CClassConstSP const TYPE;

    /** Default Constructor */
    SurfaceInterpQuadratic();

    /** Constructor */
    SurfaceInterpQuadratic(
        const DoubleArraySP    xArray,
        const DoubleArraySP    yArray,
        const CDoubleMatrixSP  zMatrix,
        ExtrapMethod           extrapMethod = SurfaceInterpQuadratic::EXTRAP_FLAT);

    /** Copy constructor */
    SurfaceInterpQuadratic(const SurfaceInterpQuadratic& rhs);

    /** Destructor */
    virtual ~SurfaceInterpQuadratic();

    /** Get interpolation method */
    virtual const SurfaceInterp::InterpMethod interpMethod() const { return SurfaceInterp::INTERP_QUADRATIC; }

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
typedef smartConstPtr<SurfaceInterpQuadratic> SurfaceInterpQuadraticConstSP;
typedef smartPtr<SurfaceInterpQuadratic> SurfaceInterpQuadraticSP;

DRLIB_END_NAMESPACE

#endif
