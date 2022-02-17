
/** Null exists so we can identify the non-existence of an object when 
converting to/from a DataDictionary
*/
#ifndef EDG_NULL_H
#define EDG_NULL_H

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL CNull: public CObject {
public:
    static CClassConstSP const TYPE;

    //// overridden for performance
    virtual IObject* clone() const;

    /** Create an instance of CNull - forcing clients to go through
        this factory method allows us to reuse the same object for all
        clients. It must not be used before this class has been
        loaded  */
    static IObjectSP create();
private:
    static void load(CClassSP& clazz);
    CNull();
    static IObjectSP commonCNullObject;
};

DRLIB_END_NAMESPACE

#endif
