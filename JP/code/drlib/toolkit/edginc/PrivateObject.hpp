
#ifndef EDG_PRIVATEOBJECT_H
#define EDG_PRIVATEOBJECT_H

DRLIB_BEGIN_NAMESPACE
class IPublicObject;

/** A private object is one which does not want its fields used as a
    public interface. Such an object needs to implement the toPublicObject
    method below. 

    An object which implements the IPrivateObject will be converted to the
    corresponding IPublicObject before, for example, being written to file or
    being converted to a data dictionary. 

    An object which implements the IPublicObject will be converted to the
    corresponding IPrivateObject after, for example, being read from file or
    being converted from a data dictionary. 
*/
class TOOLKIT_DLL IPrivateObject: virtual public IObject{
public:
    /** This method will be used for example when
        serialising out the class contents via XML */
    virtual IPublicObject* toPublicObject() const = 0;
    /** the class representing the IPrivateObject interface */
    static CClassConstSP const TYPE; // in Object.cpp
};
DRLIB_END_NAMESPACE
#endif

