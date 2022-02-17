//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterpBiLinear.hpp
//
//   Description : Linear interpolator for the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 04-Dec-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SurfaceInterpBiLinear_HPP
#define QLIB_SurfaceInterpBiLinear_HPP
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/SurfaceInterp.hpp"

DRLIB_BEGIN_NAMESPACE

/** BiLinear interpolation class.
 */
class UTIL_DLL SurfaceInterpBiLinear : public SurfaceInterp
{
public:
    static CClassConstSP const TYPE;

    /** Default Constructor */
    SurfaceInterpBiLinear();

    /** Constructor */
    SurfaceInterpBiLinear(
        const DoubleArraySP    xArray,
        const DoubleArraySP    yArray,
        const CDoubleMatrixSP  zMatrix,
        ExtrapMethod           extrapMethod = SurfaceInterpBiLinear::EXTRAP_FLAT);

    /** Copy constructor */
    SurfaceInterpBiLinear(const SurfaceInterpBiLinear& rhs);

    /** Destructor */
    virtual ~SurfaceInterpBiLinear();

    /** Get interpolation method */
    virtual const SurfaceInterp::InterpMethod interpMethod() const { return SurfaceInterp::INTERP_BILINEAR; }

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
typedef smartConstPtr<SurfaceInterpBiLinear> SurfaceInterpBiLinearConstSP;
typedef smartPtr<SurfaceInterpBiLinear> SurfaceInterpBiLinearSP;

DRLIB_END_NAMESPACE

#endif
