//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : Base class for enums wrapped up as CObjects
//
//   Author      : Mark A Robson
//
//   Date        : 14 Sep 2006
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_ENUM_HPP
#define QLIB_ENUM_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

/** Wrapper for enums ie implements IObject interface. Note that all types
    of enums get wrapped in this class but that we create a different
    CClass for each flavour of enum. */
class TOOLKIT_DLL Enum: public CObject{
public:
    static CClassConstSP const TYPE;
    virtual ~Enum();

    //// this object is immutable - clone() just returns this
    virtual IObject* clone() const;

    /** Returns string corresponding to enum integer */
    string toString() const;

    /** Create a wrapped enum of the specified enumType and with specified
        value specified as a string */
    static IObject* create(CClassConstSP requiredType, 
                           const string& enumValAsString);

#if 0
    /** Create a wrapped enum of the specified enumType and with specified
        enumValue */
    static EnumSP create(CClassConstSP enumType,
                         int           enumValue);

    /** Create a wrapped enum of the specified enumType and with specified
        value specified as a string */
    static EnumSP create(CClassConstSP enumType,
                         const string& enumValAsString);
#endif

    //// returns the wrapped enum's value as an int
    int enumValueAsInt() const;

    //// returns the wrapped enum's value as a string
    const string& enumValueAsString() const;

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;

    /** populate an empty object from description */
    virtual void import(Reader::Node* elem, Reader* reader);

    //// Writes prefix: <type name>::<Enum as string>
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Returns a hash code value for the object */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one */
    virtual bool equalTo(const IObject* obj) const;

protected:
    Enum(CClassConstSP enumType, int enumVal);
    Enum(CClassConstSP enumType);
    /// fields ///
    const int enumVal; // $unregistered

private:
    Enum(const Enum& rhs);
    Enum& operator=(const Enum& rhs);
    
    static void load(CClassSP& classToLoad);
};

typedef smartConstPtr<Enum> EnumConstSP;
typedef smartPtr<Enum> EnumSP;

#ifndef QLIB_ENUM_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<Enum>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<Enum>);
#endif


DRLIB_END_NAMESPACE
#endif

