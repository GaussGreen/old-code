#ifndef EDG_PUBLICOBJECT_H
#define EDG_PUBLICOBJECT_H

DRLIB_BEGIN_NAMESPACE
class IPrivateObject;

/** A public object is a wrapper for a private object (see IPrivateObject).
    It represents the public interface to a private object.

    An object which implements the IPrivateObject will be converted to the
    corresponding IPublicObject before, for example, being written to file or
    being converted to a data dictionary. 

    An object which implements the IPublicObject will be converted to the
    corresponding IPrivateObject after, for example, being read from file or
    being converted from a data dictionary. 
*/
class TOOLKIT_DLL IPublicObject: virtual public IObject{
public:
    /** This method will be used for example when
        serialising in the class contents from XML. The method allows
        the true internal object to be created from the public one */
    virtual IPrivateObject* toPrivateObject() const = 0;
    /** the class representing the IPublicObject interface */
    static CClassConstSP const TYPE; // in Object.cpp
};

DRLIB_END_NAMESPACE
#endif

