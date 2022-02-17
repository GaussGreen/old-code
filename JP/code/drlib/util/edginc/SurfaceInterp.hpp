//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterp.hpp
//
//   Description : Interpolator class associated with the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 04-Dec-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SurfaceInterp_HPP
#define QLIB_SurfaceInterp_HPP
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE


/** Abstract Base class for Surface interpolations.
 */
class UTIL_DLL SurfaceInterp : public CObject
{
public:
    static CClassConstSP const TYPE;

    /** Enum InterpMethod */
    enum InterpMethod
    {
        INTERP_NONE,        // No interpolation - requested coordinates must be node points.
        INTERP_FLAT,        // Piecewise flat.
        INTERP_BILINEAR,    // BiLinear interpolation.
        INTERP_QUADRATIC    // Quadratic interpolation.
    };

    /** Enum ExtrapMethod */
    enum ExtrapMethod
    {
        EXTRAP_ZERO,        // Values 0's outside boundary.
        EXTRAP_FLAT,        // Piecewise flat.
        EXTRAP_BILINEAR     // BiLinear extrapolation (extrapolate using previous two points on boundary).
    };

    /** Default Constructor */
    SurfaceInterp();

    /** Constructor */
    SurfaceInterp(
        const DoubleArraySP    xArray,
        const DoubleArraySP    yArray,
        const CDoubleMatrixSP  zMatrix,
        ExtrapMethod           extrapMethod = SurfaceInterp::EXTRAP_FLAT);

    /** Copy constructors */
    SurfaceInterp(const CClassConstSP& clazz);
    SurfaceInterp(const SurfaceInterp& rhs);
    SurfaceInterp& operator=(const SurfaceInterp& rhs);

    /** Destructor */
    virtual ~SurfaceInterp();

    /** Get surface data */
    DoubleArraySP  Xs()  const { return m_xArray; }
    DoubleArraySP  Ys()  const { return m_yArray; }
    CDoubleMatrixSP Zs() const { return m_zMatrix; }

    /** Set surface data */
    void Xs(const DoubleArraySP  xArray)  { m_xArray = xArray; }
    void Ys(const DoubleArraySP  yArray)  { m_yArray = yArray; }
    void Zs(const CDoubleMatrixSP zMatrix) { m_zMatrix = zMatrix; }
    void setXYZ(const DoubleArraySP xArray, const DoubleArraySP yArray, const CDoubleMatrixSP zMatrix);

    /** Returns the number of x values (abscissas) */
    int numXs() const { return m_xArray.get() ? m_xArray->size() : 0; }

    /** Returns the number of y values (ordinates) */
    int numYs() const { return m_yArray.get() ? m_yArray->size() : 0; }

    /** Set/Get extrapolation method */
    void extrapMethod(const ExtrapMethod extrapMethod) { m_extrapMethod = extrapMethod; }
    ExtrapMethod extrapMethod() const { return m_extrapMethod; }

    /** Get interpolation method */
    virtual const SurfaceInterp::InterpMethod interpMethod() const = 0;

    /** Interpolation functions */
    virtual double        operator()(double x, double y) const = 0;
    virtual double        interpolate(double x, double y) const = 0;
    virtual DoubleArraySP interpolateX(double x) const;
    virtual DoubleArraySP interpolateY(double y) const;

    /** static createSurfaceInterp - returns a new interpolator object */
    static SurfaceInterp* createSurfaceInterp(
        const DoubleArraySP    xArray,
        const DoubleArraySP    yArray,
        const CDoubleMatrixSP  zMatrix,
        InterpMethod           interpMethod = SurfaceInterp::INTERP_BILINEAR,
        ExtrapMethod           extrapMethod = SurfaceInterp::EXTRAP_FLAT);

protected:

    // Protected members
    DoubleArraySP    m_xArray;          // x-axis (abscissas)
    DoubleArraySP    m_yArray;          // y-axis (ordinates)
    CDoubleMatrixSP  m_zMatrix;         // z-axis (2d matrix)
    ExtrapMethod     m_extrapMethod;    // extrapolation method

    /** Return index into array at point to interpolate/extrapolate */
    static int findIndex(const DoubleArray& sortedArray, double dbl, bool& extrap);

private:

    static void load(CClassSP& clazz);

};


// Support for smart pointers, arrays etc.
typedef smartConstPtr<SurfaceInterp> SurfaceInterpConstSP;
typedef smartPtr<SurfaceInterp> SurfaceInterpSP;

DRLIB_END_NAMESPACE

#endif
