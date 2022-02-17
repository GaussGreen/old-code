//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterpFlat.cpp
//
//   Description : Piecewise Flat interpolator for the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Dec-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SurfaceInterpFlat_CPP
#include "edginc/SurfaceInterpFlat.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

/** Default Constructor */
SurfaceInterpFlat::SurfaceInterpFlat()
{
}


/** Constructor */
SurfaceInterpFlat::SurfaceInterpFlat(
    const DoubleArraySP   xArray,
    const DoubleArraySP   yArray,
    const CDoubleMatrixSP zMatrix,
    ExtrapMethod          extrapMethod)     // (Optional)
    : 
    SurfaceInterp(xArray, yArray, zMatrix, extrapMethod)
{
}


/** Copy Constructors */
SurfaceInterpFlat::SurfaceInterpFlat(const SurfaceInterpFlat& rhs)
{
    if (this != &rhs)
        *this = rhs;
}


/** Destructor */
SurfaceInterpFlat::~SurfaceInterpFlat()
{
}


/** SurfaceInterpFlat::operator() */
double SurfaceInterpFlat::operator()(double x, double y) const
{
    // Get surface data
    const DoubleArray&  arrX = *(m_xArray.get());
    const DoubleArray&  arrY = *(m_yArray.get());
    const DoubleMatrix& matZ = *(m_zMatrix.get());
    const int nX = numXs();
    const int nY = numYs();
    double z = 0.0;

    // Find nodes
    bool extrapX;
    bool extrapY;
    int i = SurfaceInterp::findIndex(arrX, x, extrapX);
    int j = SurfaceInterp::findIndex(arrY, y, extrapY);

    // Check for extrapolation and perform interpolation
    if ((extrapX || extrapY) && m_extrapMethod == SurfaceInterp::EXTRAP_ZERO)
    {
        // EXTRAP_ZERO    
        z = 0.0;
    }
    else if ((extrapX || extrapY) && m_extrapMethod == SurfaceInterp::EXTRAP_BILINEAR)
    {
        // EXTRAP_BILINEAR
        // Compute fractions in quadrant
        double fracX = (x - arrX[i]) / (arrX[i+1] - arrX[i]);
        double fracY = (y - arrY[j]) / (arrY[j+1] - arrY[j]);
        z = matZ[j][i] * (1.0 - fracY) * (1.0 - fracX)
          + matZ[j+1][i] * fracY * (1.0 - fracX)
          + matZ[j][i+1] * (1.0 - fracY) * fracX
          + matZ[j+1][i+1] * fracY * fracX;
    }
    else
    {
        // INTERP_FLAT and EXTRAP_FLAT
        if (x <= arrX[0])    i = 0;
        if (x >= arrX[nX-1]) i = nX - 1;
        if (y <= arrY[0])    j = 0;
        if (y >= arrY[nY-1]) j = nY - 1;
        z = matZ[j][i];
    }
    
    return z;
}


/** SurfaceInterpFlat::interpolate() */
double SurfaceInterpFlat::interpolate(double x, double y) const
{
    return (*this)(x,y);  // calls operator()
}


/** Surface::load */
void SurfaceInterpFlat::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Flat surface interpolation.");
    REGISTER(SurfaceInterpFlat, clazz);
    SUPERCLASS(SurfaceInterp);
}


/** Register */
CClassConstSP const SurfaceInterpFlat::TYPE = CClass::registerClassLoadMethod(
    "SurfaceInterpFlat", typeid(SurfaceInterpFlat), SurfaceInterpFlat::load);

DRLIB_END_NAMESPACE
