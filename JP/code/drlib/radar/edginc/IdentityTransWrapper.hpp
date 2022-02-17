#ifndef POLY_BASIS_WRAPPER_H
#define POLY_BASIS_WRAPPER_H

#include "edginc/DECLARE.hpp"
#include "edginc/IdentityTransform.hpp"
#include "edginc/IFittingVarTransWrapper.hpp"

/*This class 'IdentityTransWrapper' acts as a wrapper around 'IdentityTransform' with all the
reflection information and is the one which is exposed to Excel.  One needs to create, in Excel,
an object of type 'IdentityTransWrapper' and pass it to 'KRadarRepGenerator'. 
Author: Kranthi K. Gade/Vladimir A. Grebinskiy
Date: 07/18/2006.
*/

DRLIB_BEGIN_NAMESPACE
class  RADAR_DLL IdentityTransWrapper : public CObject, 
                         virtual public IFittingVarTransWrapper
{
public:
    static CClassConstSP const TYPE;
    IdentityTransWrapper() : CObject(TYPE){}
    IdentityTransWrapper(CClassConstSP const &type) : CObject(type) {}

    virtual void validatePop2Object(void);

    virtual IFittingVarTransformSP getTransform() const;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new IdentityTransWrapper(TYPE); }
};

DECLARE( IdentityTransWrapper);

DRLIB_END_NAMESPACE

#endif
