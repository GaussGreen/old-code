//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : XField.cpp
//
//   Description : Describes fields for classes from other DRI libraries
//
//   Author      : Mark A Robson
//
//   Date        : 27 Feb 2004
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/XField.hpp"
#include "edginc/XClass.hpp"


DRLIB_BEGIN_NAMESPACE

/** creates a field using the data given */
XField* XField::create(XClassSP         declaringClass,
                       const string&    name,
                       XClassSP         type,
                       const string&    description,
                       bool             isOptional){
    return new XField(declaringClass, name, type, description, isOptional);
}
 
/** Returns the name of the field represented by this Field object. */
const string& XField::getName() const throw(){
    return name;
}

/** Returns a Class object that identifies the declared type
        for the field represented by this Field object. */
XClassConstSP XField::getType() const throw(){
    return type;
}

/** Returns the Class object representing the class or
            interface that declares the field represented by this
            Field object. */
XClassConstSP XField::getDeclaringClass() const throw(){
    return declaringClass;
}

/** indicates whether the the type of data held by the field is
        atomic */
bool XField::typeIsPrimitive() const throw(){
    return type->isPrimitive();
}

/** indicates whether the the type of data held by the field is
        an array */
bool XField::typeIsArray() const throw(){
    return type->isArray();
}

/** creates a field using the data given */
XField::XField(XClassSP         declaringClass,
               const string&    name,
               XClassSP         type,
               const string&    description,
               bool             isOptional):
    CObject(TYPE), name(name), declaringClass(declaringClass), 
    type(type), optional(isOptional), description(description){}

/** clones a field */
IObject* XField::clone() const{
    return const_cast<XField*>(this);
}

/* is the field optional */
bool XField::isOptional() const{
    return optional;
}

/* is the field optional */
const string& XField::getDescription() const{
    return description;
}

static void myLoad(CClassSP& clazz){
    REGISTER(XField, clazz);
    SUPERCLASS(CObject);
}

CClassConstSP const XField::TYPE = CClass::registerClassLoadMethod(
    "XField", typeid(XField), myLoad);
    

DRLIB_END_NAMESPACE
