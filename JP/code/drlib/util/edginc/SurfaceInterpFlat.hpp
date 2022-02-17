//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterpFlat.hpp
//
//   Description : Piecewise Flat interpolator for the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Dec-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SurfaceInterpFlat_HPP
#define QLIB_SurfaceInterpFlat_HPP
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/SurfaceInterp.hpp"

DRLIB_BEGIN_NAMESPACE

/** Flat interpolation class.
 */
class UTIL_DLL SurfaceInterpFlat : public SurfaceInterp
{
public:
    static CClassConstSP const TYPE;

    /** Default Constructor */
    SurfaceInterpFlat();

    /** Constructor */
    SurfaceInterpFlat(
        const DoubleArraySP    xArray,
        const DoubleArraySP    yArray,
        const CDoubleMatrixSP  zMatrix,
        ExtrapMethod           extrapMethod = SurfaceInterpFlat::EXTRAP_FLAT);

    /** Copy constructor */
    SurfaceInterpFlat(const SurfaceInterpFlat& rhs);

    /** Destructor */
    virtual ~SurfaceInterpFlat();

    /** Get interpolation method */
    virtual const SurfaceInterp::InterpMethod interpMethod() const { return SurfaceInterp::INTERP_FLAT; }

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
typedef smartConstPtr<SurfaceInterpFlat> SurfaceInterpFlatConstSP;
typedef smartPtr<SurfaceInterpFlat> SurfaceInterpFlatSP;

DRLIB_END_NAMESPACE

#endif
