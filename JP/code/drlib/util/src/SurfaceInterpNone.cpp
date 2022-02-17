//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterpNone.cpp
//
//   Description : 'None' interpolator for the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Dec-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SurfaceInterpNone_CPP
#include "edginc/SurfaceInterpNone.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

/** Default Constructor */
SurfaceInterpNone::SurfaceInterpNone()
{
}


/** Constructor */
SurfaceInterpNone::SurfaceInterpNone(
    const DoubleArraySP   xArray,
    const DoubleArraySP   yArray,
    const CDoubleMatrixSP zMatrix,
    ExtrapMethod          extrapMethod)     // (Optional)
    : 
    SurfaceInterp(xArray, yArray, zMatrix, extrapMethod)
{
}


/** Copy Constructors */
SurfaceInterpNone::SurfaceInterpNone(const SurfaceInterpNone& rhs)
{
    if (this != &rhs)
        *this = rhs;
}


/** Destructor */
SurfaceInterpNone::~SurfaceInterpNone()
{
}


/** SurfaceInterpNone::operator() */
double SurfaceInterpNone::operator()(double x, double y) const
{
    static const string method = "SurfaceInterpNone::operator()";

    // Get surface data
    const DoubleArray&  arrX = *(m_xArray.get());
    const DoubleArray&  arrY = *(m_yArray.get());
    const DoubleMatrix& matZ = *(m_zMatrix.get());
    const int nX = numXs();
    const int nY = numYs();
    double    z = 0.0;

    // Find nodes
    bool extrapX;
    bool extrapY;
    int i = SurfaceInterp::findIndex(arrX, x, extrapX);
    int j = SurfaceInterp::findIndex(arrY, y, extrapY);

    // Return point (extrapolation not applicable)
    bool exact = true;
    if (x == arrX[i])
    {
        if (y == arrY[j])
            z = matZ[j][i];
        else if (j+1 < nY && y == arrY[j+1])
            z = matZ[j+1][i];
        else
            exact = false;
    }
    else if (i+1 < nX && x == arrX[i+1])
    {
        if (y == arrY[j])
            z = matZ[j][i+1];
        else if (j+1 < nY && y == arrY[j+1])
            z = matZ[j+1][i+1];
        else
            exact = false;
    }
    
    // Check if exact node point found
    if (!exact)
        throw ModelException(method, "Node point (" + Format::toString(x) + "," + Format::toString(y)
                             + ") node found in zMatrix (interpolation is set to NONE)!");
    return z;
}


/** SurfaceInterpNone::interpolate() */
double SurfaceInterpNone::interpolate(double x, double y) const
{
    return (*this)(x,y);  // calls operator()
}


/** Surface::load */
void SurfaceInterpNone::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("None surface interpolation.");
    REGISTER(SurfaceInterpNone, clazz);
    SUPERCLASS(SurfaceInterp);
}


/** Register */
CClassConstSP const SurfaceInterpNone::TYPE = CClass::registerClassLoadMethod(
    "SurfaceInterpNone", typeid(SurfaceInterpNone), SurfaceInterpNone::load);

DRLIB_END_NAMESPACE
