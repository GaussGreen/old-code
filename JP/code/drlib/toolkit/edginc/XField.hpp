//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : XField.hpp
//
//   Description : Describes fields for classes from other DRI libraries
//
//   Author      : Mark A Robson
//
//   Date        : 27 Feb 2004
//
//----------------------------------------------------------------------------

#ifndef EDR_XFIELD_HPP
#define EDR_XFIELD_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

class XField;
typedef const XField* XFieldConstSP;
typedef XField* XFieldSP;
typedef vector<XFieldConstSP> XFieldArray;
class XClass;
typedef const XClass* XClassConstSP;
typedef XClass* XClassSP;


/** A XField provides information about a single field of a class
    belonging to another DRI library */
class TOOLKIT_DLL XField: public CObject{
    friend class ShutTheCompilerUp;

public:
    /** clones a field */
    virtual IObject* clone() const;

    static CClassConstSP const TYPE;

    /** creates a field using the data given */
    static XField* create(XClassSP         declaringClass,
                          const string&    name,
                          XClassSP         type,
                          const string&    description,
                          bool             isOptional);

    /** Returns the name of the field represented by this Field object. */
    const string& getName() const throw();

    /** Returns a Class object that identifies the declared type
        for the field represented by this Field object. */
    XClassConstSP getType() const throw();

    /** Returns the Class object representing the class or
        interface that declares the field represented by this
        Field object. */
    XClassConstSP getDeclaringClass() const throw();

    /** indicates whether the the type of data held by the field is
        atomic */
    bool typeIsPrimitive() const throw();

    /** indicates whether the the type of data held by the field is
        an array */
    bool typeIsArray() const throw();
        
    /** flag the optionality of a field */
    void setOptional(bool isOptional);

    /* is the field optional */
    bool isOptional() const;

    /* is the field optional */
    const string& getDescription() const;

private:
    // Fields must be accessed via pointers
    XField(XClassSP         declaringClass,
           const string&    name,
           XClassSP         type,
           const string&    description,
           bool             isOptional);
    XField();
    XField(const XField &rhs);
    XField& operator=(const XField& rhs);
    
    // **** fields *****
    string        name;     // the name of the field $unregistered
    XClassConstSP declaringClass;     // to which class the field belongs $unregistered
    XClassConstSP type;// the type of the element the field contains $unregistered
    bool          optional; // $unregistered
    string        description; // $unregistered
};

DRLIB_END_NAMESPACE

#endif
