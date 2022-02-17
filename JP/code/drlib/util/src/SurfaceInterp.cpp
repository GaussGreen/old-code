//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : SurfaceInterp.cpp
//
//   Description : Interpolator class associated with the Surface object.
//
//   Author      : Anwar E Sidat
//
//   Date        : 04-Dec-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_SurfaceInterp_CPP
#include "edginc/SurfaceInterp.hpp"
#include "edginc/BoxedEnum.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/SurfaceInterpNone.hpp"
#include "edginc/SurfaceInterpFlat.hpp"
#include "edginc/SurfaceInterpBiLinear.hpp"
#include "edginc/SurfaceInterpQuadratic.hpp"

DRLIB_BEGIN_NAMESPACE

// Enum InterpMethod
START_PUBLIC_ENUM_DEFINITION(SurfaceInterp::InterpMethod, "Interpolation method applicable to a Surface.");
    ENUM_VALUE_AND_NAME(SurfaceInterp::INTERP_NONE,      "None",      "No interpolation - requested coordinates must be node points.");
    ENUM_VALUE_AND_NAME(SurfaceInterp::INTERP_FLAT,      "Flat",      "Piecewise flat.");
    ENUM_VALUE_AND_NAME(SurfaceInterp::INTERP_BILINEAR,  "BiLinear",  "BiLinear interpolation.");
    ENUM_VALUE_AND_NAME(SurfaceInterp::INTERP_QUADRATIC, "Quadratic", "Quadratic interpolation.");
END_ENUM_DEFINITION(SurfaceInterp::InterpMethod);

// Enum ExtrapMethod
START_PUBLIC_ENUM_DEFINITION(SurfaceInterp::ExtrapMethod, "Extrapolation method applicable to a Surface.");
    ENUM_VALUE_AND_NAME(SurfaceInterp::EXTRAP_ZERO,      "Zero",      "Values 0's outside boundary.");
    ENUM_VALUE_AND_NAME(SurfaceInterp::EXTRAP_FLAT,      "Flat",      "Piecewise flat.");
    ENUM_VALUE_AND_NAME(SurfaceInterp::EXTRAP_BILINEAR,  "BiLinear",  "BiLinear extrapolation (extrapolate using previous two points on boundary).");
END_ENUM_DEFINITION(SurfaceInterp::ExtrapMethod);


/** Default Constructor */
SurfaceInterp::SurfaceInterp()
    :
    CObject(TYPE),
    m_extrapMethod(EXTRAP_FLAT)
{
}


/** Constructor */
SurfaceInterp::SurfaceInterp(
    const DoubleArraySP   xArray,
    const DoubleArraySP   yArray,
    const CDoubleMatrixSP zMatrix,
    ExtrapMethod          extrapMethod)     // (Optional)
    : 
    CObject(TYPE),
    m_xArray(xArray),
    m_yArray(yArray),
    m_zMatrix(zMatrix),
    m_extrapMethod(extrapMethod)
{
}


/** Copy Constructors */
SurfaceInterp::SurfaceInterp(const CClassConstSP& clazz)
    : CObject(clazz)
{
}

SurfaceInterp::SurfaceInterp(const SurfaceInterp& rhs)
    : CObject(rhs.TYPE)
{
    if (this != &rhs)
        *this = rhs;
}


SurfaceInterp& SurfaceInterp::operator=(const SurfaceInterp& rhs)
{
    if (this != &rhs)
    {
        m_xArray = rhs.m_xArray;
        m_yArray = rhs.m_yArray;
        m_zMatrix = rhs.m_zMatrix;
        m_extrapMethod = rhs.m_extrapMethod;
    }
    return *this;
}

/** Destructor */
SurfaceInterp::~SurfaceInterp()
{
}


/** SurfaceInterp::setXYZ */
void SurfaceInterp::setXYZ(
    const DoubleArraySP   xArray,
    const DoubleArraySP   yArray,
    const CDoubleMatrixSP zMatrix)
{
    m_xArray = xArray;
    m_yArray = yArray;
    m_zMatrix = zMatrix;
}


/** SurfaceInterp::interpolateX() */
DoubleArraySP SurfaceInterp::interpolateX(double y) const
{
    // Allocate result
    DoubleArraySP xArray;
    int n = numXs();
    if (n > 0)
        xArray.reset(new DoubleArray(m_xArray->size()));

    // Interpolate column
    for (int i = 0; n; ++i)
        (*xArray)[i] = interpolate((*m_xArray)[i], y);
    
    return xArray;
}


/** SurfaceInterp::interpolateY() */
DoubleArraySP SurfaceInterp::interpolateY(double x) const
{
    // Allocate result
    DoubleArraySP yArray;
    int n = numYs();
    if (n > 0)
        yArray.reset(new DoubleArray(m_yArray->size()));

    // Interpolate row
    for (int j = 0; n; ++j)
        (*yArray)[j] = interpolate(x, (*m_yArray)[j]);
    
    return yArray;
}

/** static Surface::findIndex */
//  Returns value of index for interpolating between index and index+1.
//  Index=0 if extrapolating at start, Index=n-2 if extrapolating at end.
int SurfaceInterp::findIndex(const DoubleArray& sortedArray, double d, bool& extrap)
{
    static const string method = "Surface::findIndex";
    int n = sortedArray.size();
    extrap = true;
    if (d < sortedArray[0])
        return 0;
    else if (d > sortedArray[n-1])
        return (n-2 >= 0 ? n-2 : 0);
    else
    {
        // Binary search
        int i;
        int lo = 0;
        int hi = n - 1;
        while ((hi - lo) > 1)
        {
            i = (lo + hi) >> 1;                 // determine mid point.
            if (d >= sortedArray[i]) 
                lo = i;                         // find in upper half.
            else if (d < sortedArray[i]) 
                hi = i;                         // find in lower half.
        }
        extrap = false;
        return lo;
    }
    
    // Return error
    throw ModelException(method, "Failed to find item in Array - most likely input array was not sorted!");
    return -1.0;
}

/** static createSurfaceInterp */
SurfaceInterp* SurfaceInterp::createSurfaceInterp(
    const DoubleArraySP    xArray,
    const DoubleArraySP    yArray,
    const CDoubleMatrixSP  zMatrix,
    InterpMethod           interpMethod,    // (opt) defaults to INTERP_BILINEAR
    ExtrapMethod           extrapMethod)    // (opt) defaults to EXTRAP_FLAT
{
    static const string method = "Surface::createSurfaceInterp";
    switch (interpMethod)
    {
        case INTERP_NONE:
            return new SurfaceInterpNone(xArray, yArray, zMatrix, extrapMethod);
            break;
        case INTERP_FLAT:
            return new SurfaceInterpFlat(xArray, yArray, zMatrix, extrapMethod);
            break;
        case INTERP_BILINEAR:
            return new SurfaceInterpBiLinear(xArray, yArray, zMatrix, extrapMethod);
            break;
        case INTERP_QUADRATIC:
            return new SurfaceInterpQuadratic(xArray, yArray, zMatrix, extrapMethod);
            break;
        default:
            throw ModelException(method, "Invalid InterpMethod: " + Format::toString(interpMethod));
    }
    return 0;
}


/** Surface::load */
void SurfaceInterp::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Base class for surface interpolation methods.");
    REGISTER(SurfaceInterp, clazz);
    SUPERCLASS(CObject);
}


/** Register */
CClassConstSP const SurfaceInterp::TYPE = CClass::registerClassLoadMethod(
    "SurfaceInterp", typeid(SurfaceInterp), SurfaceInterp::load);

DRLIB_END_NAMESPACE
