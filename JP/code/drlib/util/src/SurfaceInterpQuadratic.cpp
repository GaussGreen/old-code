//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterpQuadratic.cpp
//
//   Description : Quadratic interpolator for the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Dec-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SurfaceInterpQuadratic_CPP
#include "edginc/SurfaceInterpQuadratic.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

/** Default Constructor */
SurfaceInterpQuadratic::SurfaceInterpQuadratic()
{
}


/** Constructor */
SurfaceInterpQuadratic::SurfaceInterpQuadratic(
    const DoubleArraySP   xArray,
    const DoubleArraySP   yArray,
    const CDoubleMatrixSP zMatrix,
    ExtrapMethod          extrapMethod)     // (Optional)
    : 
    SurfaceInterp(xArray, yArray, zMatrix, extrapMethod)
{
}


/** Copy Constructors */
SurfaceInterpQuadratic::SurfaceInterpQuadratic(const SurfaceInterpQuadratic& rhs)
{
    if (this != &rhs)
        *this = rhs;
}


/** Destructor */
SurfaceInterpQuadratic::~SurfaceInterpQuadratic()
{
}


/** SurfaceInterpQuadratic::operator() */
double SurfaceInterpQuadratic::operator()(double x, double y) const
{
    // Get surface data
    const DoubleArray&  arrX = *(m_xArray.get());
    const DoubleArray&  arrY = *(m_yArray.get());
    const DoubleMatrix& matZ = *(m_zMatrix.get());
    const int nX = numXs();
    const int nY = numYs();
    double    z = 0.0;

    // Find nodes
    bool extrapX; // = ((x < arrX[0]) || (x > arrX[nX-1]));
    bool extrapY; // = ((y < arrY[0]) || (y > arrY[nY-1]));
    int i = SurfaceInterp::findIndex(arrX, x, extrapX);
    int j = SurfaceInterp::findIndex(arrY, y, extrapY);

    if ((extrapX || extrapY) && m_extrapMethod == SurfaceInterp::EXTRAP_ZERO)
    {
        // EXTRAP_ZERO    
        z = 0.0;
    }
    else if ((extrapX || extrapY) && m_extrapMethod == SurfaceInterp::EXTRAP_FLAT)
    {
        // EXTRAP_FLAT
        if (x <= arrX[0])    i = 0;
        if (x >= arrX[nX-1]) i = nX - 1;
        if (y <= arrY[0])    j = 0;
        if (y >= arrY[nY-1]) j = nY - 1;
        z = matZ[j][i];
    }
    else
    {
        /**********************************************************************************/
        throw ModelException(__FUNCTION__, "Sorry, Quadratic method not yet implemented!");
        /**********************************************************************************/

        // INTERP_Quadratic and EXTRAP_Quadratic
        // Compute fractions in quadrant (reverts to FLAT if insufficient points for linear)
        double fracX = (i < nX-1) ? (x - arrX[i]) / (arrX[i+1] - arrX[i]) : 0.0;
        double fracY = (j < nY-1) ? (y - arrY[j]) / (arrY[j+1] - arrY[j]) : 0.0;

        z = matZ[j][i] * (1.0 - fracY) * (1.0 - fracX)
          + matZ[j+1][i] * fracY * (1.0 - fracX)
          + matZ[j][i+1] * (1.0 - fracY) * fracX
          + matZ[j+1][i+1] * fracY * fracX;
    }
    
    return z;
}


/** SurfaceInterpQuadratic::interpolate() */
double SurfaceInterpQuadratic::interpolate(double x, double y) const
{
    return (*this)(x,y);  // calls operator()
}


/** Surface::load */
void SurfaceInterpQuadratic::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Quadratic surface interpolation.");
    REGISTER(SurfaceInterpQuadratic, clazz);
    SUPERCLASS(SurfaceInterp);
}


/** Register */
CClassConstSP const SurfaceInterpQuadratic::TYPE = CClass::registerClassLoadMethod(
    "SurfaceInterpQuadratic", typeid(SurfaceInterpQuadratic), SurfaceInterpQuadratic::load);

DRLIB_END_NAMESPACE
