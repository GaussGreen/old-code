//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : Surface.cpp
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

#include "edginc/config.hpp"
#define QLIB_Surface_CPP
#include "edginc/Surface.hpp"
#include "edginc/SurfaceInterp.hpp"
#include "edginc/SurfaceInterpBiLinear.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

// Make sure the class links
bool SurfaceLoad() { return (Surface::TYPE != 0); }

/** Constructors */
Surface::Surface()
    :
    CObject(TYPE),
    m_surfaceInterp(new SurfaceInterpBiLinear())
    
{
    m_surfaceInterp->extrapMethod(SurfaceInterp::EXTRAP_FLAT);
}


Surface::Surface(const DoubleArraySP xArray, const DoubleArraySP yArray, const CDoubleMatrixSP zMatrix)
    :
    CObject(TYPE),
    xArray(xArray),
    yArray(yArray),
    zMatrix(zMatrix)
{
    m_surfaceInterp.reset(new SurfaceInterpBiLinear(xArray, yArray, zMatrix, SurfaceInterp::EXTRAP_FLAT));
}


/** Copy Constructors */
Surface::Surface(const Surface& rhs)
    :
    CObject(rhs.TYPE)
{
    if (this != &rhs)
        *this = rhs;
}


Surface::Surface(const CClassConstSP& clazz)
    :
    CObject(clazz),
    m_surfaceInterp(new SurfaceInterpBiLinear())
    
{
    m_surfaceInterp->extrapMethod(SurfaceInterp::EXTRAP_FLAT);
}


const Surface& Surface::operator=(const Surface& rhs)
{
    if (this != &rhs)
    {
        xArray = rhs.xArray;
        yArray = rhs.yArray;
        zMatrix = rhs.zMatrix;
        m_surfaceInterp = rhs.m_surfaceInterp;
    }
    return *this;
}

/** Destructor */
Surface::~Surface()
{
}


/** Surface::validatePop2Object */
void Surface::validatePop2Object()
{
    static const string method = "Surface::validatePop2Object";
    try
    {
        // Check sizes
        int nX = numXs();
        int nY = numYs();
        if (nX <= 0 || nY <= 0)
            throw ModelException(method, "Invalid vector sizes, xArray[" + Format::toString(nX) + "] and yArray[" + Format::toString(nY) + "]");
        if (nX != zMatrix->numRows() || nY != zMatrix->numCols())
            throw ModelException(method, "Invalid matrix size["
                                 + Format::toString(zMatrix->numRows()) + ","
                                 + Format::toString(zMatrix->numCols()) + "]"
                                 + "when vector sizes are xArray[" + Format::toString(nX) + "] and yArray[" + Format::toString(nY) + "]");

        // Check arrays sorted (no duplicates allowed)
        for (int i = 1; i < nX; ++i)
            if ((*xArray)[i-1] >= (*xArray)[i])
                throw ModelException(method, "xArray in Surface must be sorted and duplicates not allowed!");
        for (int j = 1; j < nY; ++j)
            if ((*yArray)[j-1] >= (*yArray)[j])
                throw ModelException(method, "yArray in Surface must be sorted and duplicates not allowed!");
                
        // Set interplator data
        m_surfaceInterp->setXYZ(xArray, yArray, zMatrix);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


/** Surface::clone */
IObject* Surface::clone() const
{
    IObject* pIObject = CObject::clone();
    pIObject->validatePop2Object();
    return pIObject;
}


/** Surface::setXYZ */
void Surface::setXYZ(const DoubleArraySP arrayX, const DoubleArraySP arrayY, const CDoubleMatrixSP matrixZ)
{
    xArray = arrayX;
    yArray = arrayY;
    zMatrix = matrixZ;
}


/** Surface::operator() */
double Surface::operator()(double x, double y) const
{
    return (*m_surfaceInterp.get())(x,y);
}


/** Surface::load */
void Surface::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Surface object.");
    REGISTER(Surface, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(xArray,  "X array (rows)");
    FIELD(yArray,  "Y array (columns)");
    FIELD(zMatrix, "Z matrix (surface data)");

    Addin::registerConstructor("Surface", Addin::RISK, "Creates a surface handle.", Surface::TYPE);
}


/** Register */
CClassConstSP const Surface::TYPE = CClass::registerClassLoadMethod(
    "Surface", typeid(Surface), Surface::load);




/** SurfaceInterpAddin */
class SurfaceInterpAddin: public CObject
{
public:
    static CClassConstSP const TYPE;
    SurfaceSP  surface;
    double     x;
    double     y;
    SurfaceInterp::InterpMethod  interpType;
    SurfaceInterp::ExtrapMethod  extrapType;
    
    SurfaceInterpAddin(): CObject(TYPE) {}
    static IObject* defaultSurfaceInterpAddin()
    {
        return new SurfaceInterpAddin();
    }
    double interpolate()
    {
        surface->interpMethod(interpType);
        surface->extrapMethod(extrapType);
        return (*surface.get())(x,y);
    }
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SurfaceInterpAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSurfaceInterpAddin);
        FIELD(surface, "surface");
        FIELD(x, "row value");
        FIELD(y, "column value");
        FIELD(interpType, "interp");
        FIELD(extrapType, "extrap");

        Addin::registerDoubleMethod("SURFACE_INTERP",
            Addin::UTILITIES,
            "Interpolates on the surface given coordinates (x,y).",
            &SurfaceInterpAddin::interpolate);                                    
    }
};

CClassConstSP const SurfaceInterpAddin::TYPE=CClass::registerClassLoadMethod(
    "SurfaceInterpAddin", typeid(SurfaceInterpAddin), load);

DRLIB_END_NAMESPACE
