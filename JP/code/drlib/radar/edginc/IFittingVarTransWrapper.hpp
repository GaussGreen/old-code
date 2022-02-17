#ifndef IFITTING_VAR_TRANS_WRAPPER_H
#define IFITTING_VAR_TRANS_WRAPPER_H

/*This class 'IFittingVarTransWrapper' is a wrapper around 'IFittingVarTransform' with all the
reflection information and is the one which is exposed to Excel.  It has a function called
'getBasis()' which returns an 'IFittingVarTransformSP' which in turn is used to create
'RadarRepCreator' and so on.  
Author : Kranthi K. Gade/Vladimir A. Grebinskiy
Date : 07/18/2006
*/

#include "edginc/config.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/IFittingVarTransform.hpp"

DRLIB_BEGIN_NAMESPACE
class RADAR_DLL  IFittingVarTransWrapper : virtual public IObject
{
public:
    static CClassConstSP const TYPE;
    virtual IFittingVarTransformSP getTransform() const = 0;
    virtual ~IFittingVarTransWrapper() {}

private:
    static void load(CClassSP& clazz);
};

DECLARE(IFittingVarTransWrapper);

DRLIB_END_NAMESPACE
#endif
