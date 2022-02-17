#ifndef IFUNCTION_BASIS_WRAPPER_H
#define IFUNCTION_BASIS_WRAPPER_H

/*This class 'IFuncBasisWrapper' is a wrapper around 'IFunctionBasis' with all the
reflection information and is the one which is exposed to Excel.  It has a function called
'getBasis()' which returns an 'IFunctionBasisSP' which in turn is used to create
'RadarRepCreator' and so on.  
Author : Kranthi K. Gade/Vladimir A. Grebinskiy
Date : 07/18/2006
*/

#include "edginc/DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/IFunctionBasis.hpp"

DRLIB_BEGIN_NAMESPACE
class  RADAR_DLL IFuncBasisWrapper : virtual public IObject
{
public:
    static CClassConstSP const TYPE;
    virtual IFunctionBasisSP getBasis() const = 0;
    virtual ~IFuncBasisWrapper() {}

private:
    static void load(CClassSP& clazz);
    //static IObject* defaultConstructor(void);
};

DECLARE(IFuncBasisWrapper);

DRLIB_END_NAMESPACE
#endif
