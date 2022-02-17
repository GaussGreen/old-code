#ifndef POLY_BASIS_WRAPPER_H
#define POLY_BASIS_WRAPPER_H

#include "edginc/config.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/PolynomialBasis.hpp"
#include "edginc/IFuncBasisWrapper.hpp"

/*This class 'PolyBasisWrapper' acts as a wrapper around 'PolynomialBasis' with all the
reflection information and is the one which is exposed to Excel.  One needs to create, in Excel,
an object of type 'PolyBasisWrapper' and pass it to 'KRadarRepGenerator'.  This class has 
as its variables 'numVar' and 'maxDeg' (same as 'PolynomialBasis').
Author: Kranthi K. Gade/Vladimir A. Grebinskiy
Date: 07/18/2006.
*/

DRLIB_BEGIN_NAMESPACE
class RADAR_DLL PolyBasisWrapper : public CObject, 
                         virtual public IFuncBasisWrapper
{
public:
    static CClassConstSP const TYPE;
    PolyBasisWrapper() : CObject(TYPE){}
    PolyBasisWrapper(CClassConstSP const &type) : CObject(type) {}
    PolyBasisWrapper(int numVars, int maxDeg);

    virtual void validatePop2Object(void);

    virtual IFunctionBasisSP getBasis() const;

private:
    int numVars;
    int maxDeg;

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new PolyBasisWrapper(TYPE); }
};

DECLARE(PolyBasisWrapper);

DRLIB_END_NAMESPACE

#endif
