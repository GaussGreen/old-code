//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterpBiLinear.cpp
//
//   Description : Linear interpolator for the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 04-Dec-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SurfaceInterpBiLinear_CPP
#include "edginc/SurfaceInterpBiLinear.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

/** Default Constructor */
SurfaceInterpBiLinear::SurfaceInterpBiLinear()
{
}


/** Constructor */
SurfaceInterpBiLinear::SurfaceInterpBiLinear(
    const DoubleArraySP   xArray,
    const DoubleArraySP   yArray,
    const CDoubleMatrixSP zMatrix,
    ExtrapMethod          extrapMethod)     // (Optional)
    : 
    SurfaceInterp(xArray, yArray, zMatrix, extrapMethod)
{
}


/** Copy Constructors */
SurfaceInterpBiLinear::SurfaceInterpBiLinear(const SurfaceInterpBiLinear& rhs)
{
    if (this != &rhs)
        *this = rhs;
}


/** Destructor */
SurfaceInterpBiLinear::~SurfaceInterpBiLinear()
{
}


/** SurfaceInterpBiLinear::operator() */
double SurfaceInterpBiLinear::operator()(double x, double y) const
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

    // Check for extrapolation and perform interpolation
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
        // INTERP_BILINEAR and EXTRAP_BILINEAR
        // Compute fractions in quadrant
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


/** SurfaceInterpBiLinear::interpolate() */
double SurfaceInterpBiLinear::interpolate(double x, double y) const
{
    return (*this)(x,y);  // calls operator()
}


/** Surface::load */
void SurfaceInterpBiLinear::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Bilinear surface interpolation.");
    REGISTER(SurfaceInterpBiLinear, clazz);
    SUPERCLASS(SurfaceInterp);
}


/** Register */
CClassConstSP const SurfaceInterpBiLinear::TYPE = CClass::registerClassLoadMethod(
    "SurfaceInterpBiLinear", typeid(SurfaceInterpBiLinear), SurfaceInterpBiLinear::load);

DRLIB_END_NAMESPACE
