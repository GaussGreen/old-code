//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : Surface.hpp
//
//   Description : Class to represent a surface as a matrix (z-axis) with the
//                 additional vector arrays to represent the x and y-axes.
//
//                 NOTE: Contrary to standard graphical conventions, the matrix
//                       rows is referred to as the x-axis here and the columns
//                       as the y-axis.  The origin, zMatrix[0][0] is the top
//                       left corner of the 2D matrix.
//
//   Author      : Anwar E Sidat
//
//   Date        : 03-Dec-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_Surface_HPP
#define QLIB_Surface_HPP
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/SurfaceInterp.hpp"

DRLIB_BEGIN_NAMESPACE

// To avoid redundant file includes.
FORWARD_DECLARE(CDoubleMatrix);

/** Class representing a surface using a 2D matrix and two vectors for the axes.
 */
class UTIL_DLL Surface : public CObject
{
public:
    static CClassConstSP const TYPE;

    /** Constructors */
    Surface();
    Surface(const DoubleArraySP xArray, const DoubleArraySP yArray, const CDoubleMatrixSP zMatrix);

    /** Destructor */
    virtual ~Surface();

    /** Get surface data */
    const DoubleArraySP   Xs() const { return xArray; }
    const DoubleArraySP   Ys() const { return yArray; }
    const CDoubleMatrixSP Zs() const { return zMatrix; }

    /** Set surface data */
    void Xs(const DoubleArraySP   arrayX)  { xArray = arrayX; }
    void Ys(const DoubleArraySP   arrayY)  { yArray = arrayY; }
    void Zs(const CDoubleMatrixSP matrixZ) { zMatrix = matrixZ; }
    void setXYZ(const DoubleArraySP arrayX, const DoubleArraySP arrayY, const CDoubleMatrixSP matrixZ);

    /** Returns the number of x values (abscissas) */
    int numXs() const { return xArray.get() ? xArray->size() : 0; }

    /** Returns the number of y values (ordinates) */
    int numYs() const { return yArray.get() ? yArray->size() : 0; }

    /** Set interpolation method */
    void interpMethod(const SurfaceInterp::InterpMethod interpMethod)
    { m_surfaceInterp.reset(SurfaceInterp::createSurfaceInterp(xArray, yArray, zMatrix, interpMethod)); }

    /** Get interpolation method */
    SurfaceInterp::InterpMethod interpMethod() const
    { return m_surfaceInterp->interpMethod(); }

    /** Set extrapolation method */
    void extrapMethod(const SurfaceInterp::ExtrapMethod extrapMethod)
    { m_surfaceInterp->extrapMethod(extrapMethod); }

    /** Get extrapolation method */
    SurfaceInterp::ExtrapMethod extrapMethod() const
    { return m_surfaceInterp->extrapMethod(); }

    /** Returns the value at the specified point */
    double operator()(double x, double y) const;

    /** Overrides validatePop2Object */
    virtual void validatePop2Object();

    /** Overrides clone */
    IObject* clone() const;

protected:

    /** Copy constructors */
    Surface(const Surface& rhs);
    Surface(const CClassConstSP& clazz);
    const Surface& operator=(const Surface& rhs);

    //Fields
    DoubleArraySP    xArray;    // x-axis (abscissas)
    DoubleArraySP    yArray;    // y-axis (ordinates)
    CDoubleMatrixSP  zMatrix;   // z-axis (2D matrix)

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new Surface(); }

    // Private members
    SurfaceInterpSP  m_surfaceInterp;   // interpolator
};


// Support for smart pointers and arrays
typedef smartConstPtr<Surface> SurfaceConstSP;
typedef smartPtr<Surface> SurfaceSP;
typedef array<SurfaceSP, Surface> SurfaceArray;

DRLIB_END_NAMESPACE

#endif
