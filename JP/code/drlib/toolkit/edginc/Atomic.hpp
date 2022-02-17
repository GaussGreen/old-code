
#ifndef EDG_ATOMIC_H
#define EDG_ATOMIC_H
#include "edginc/DRObject.hpp"
#include "edginc/TypeConvert.hpp"
#include "edginc/CombinableResult.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

class CDouble;
typedef smartConstPtr<CDouble> CDoubleConstSP;
typedef smartPtr<CDouble> CDoubleSP;

/** wrapper for double type ie implements IObject interface */
class TOOLKIT_DLL CDouble: public CObject, 
               public virtual CombinableResult,
               public virtual IDRObject{
    friend class DoubleHelper;
public:
    static CClassConstSP const TYPE;
    ~CDouble();

    // inherited 
    virtual IObject* clone() const;
    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    static CDouble* create(double d);
    static CDoubleSP SP(double d);
    double doubleValue() const;
    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object from description */
    virtual void import(Reader::Node* elem, Reader* reader);
    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Returns a hash code value for the object */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one */
    virtual bool equalTo(const IObject* obj) const;

    /** scale by factor x (Implementation of CombinableResult) */
    virtual void scale(double x);
    /** add CDouble object to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** For use in constructing e.g. error messages */
    virtual string toString() const;

    /** Returns a hash code value for the supplied double */
    static int hashCode(double db);

private:
    double db; // $unregistered
    CDouble(const double d);
    CDouble(const CDouble& rhs);
    CDouble& operator=(const CDouble& rhs);
};

#ifndef QLIB_ATOMIC_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<CDouble>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CDouble>);
#endif

class CInt;
typedef smartConstPtr<CInt> CIntConstSP;
typedef smartPtr<CInt> CIntSP;

/** wrapper for int type ie implements IObject interface */
class TOOLKIT_DLL CInt: public CObject,
            public virtual IDRObject{
    friend class IntHelper;
public:
    static CClassConstSP const TYPE;
    ~CInt();

    // inherited 
    virtual IObject* clone() const;
    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    static CInt* create(int i);
    static CIntSP SP(int i);
    int intValue() const;
    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object from description */
    virtual void import(Reader::Node* elem, Reader* reader);
    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Returns a hash code value for the object */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one */
    virtual bool equalTo(const IObject* obj) const;

    /** For use in constructing e.g. error messages */
    virtual string toString() const;

private:
    const int i; // $unregistered
    CInt(const int d);
    CInt(const CInt& rhs);
    CInt& operator=(const CInt& rhs);
};

#ifndef QLIB_ATOMIC_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<CInt>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CInt>);
#endif

/** wrapper for bool type ie implements IObject interface */
class TOOLKIT_DLL CBool: public CObject,
             public virtual IDRObject{
    friend class BoolHelper;
public:
    static CClassConstSP const TYPE;
    ~CBool();

    // inherited 
    virtual IObject* clone() const;
    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    static CBool* create(bool i);
    bool boolValue() const;
    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object from description */
    virtual void import(Reader::Node* elem, Reader* reader);
    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Returns a hash code value for the object */
    virtual int hashCode() const;
    
    /** Returns a hash code value for the supplied boolean */
    static int hashCode(bool b);

    /** Indicates whether some other object is "equal to" this one */
    virtual bool equalTo(const IObject* obj) const;

    /* returns true for strings "Y", "YES", "1", "TRUE" and false for strings
       "N", "NO", "0", "FALSE" else throws an exception. Comparison is
       case independent */
    static bool fromString(const string& value);
    
    // so can convert strings to bools (e.g. StringArray->BoolArray) when 
    // building via DataDictionary
    static IObject* makeFromString(CClassConstSP requiredType, 
                                   const string& data);

    /** For use in constructing e.g. error messages */
    virtual string toString() const;

private:
    const bool b; // $unregistered
    CBool(const bool d);
    CBool(const CBool& rhs);
    CBool& operator=(const CBool& rhs);
};

typedef smartConstPtr<CBool> CBoolConstSP;
typedef smartPtr<CBool> CBoolSP;
#ifndef QLIB_ATOMIC_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<CBool>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CBool>);
#endif

/** wrapper for string type ie implements IObject interface */
class TOOLKIT_DLL CString: public CObject,
               public virtual ITypeConvert,
               public virtual IDRObject{
    friend class StringHelper;
public:
    static CClassConstSP const TYPE;
    ~CString();

    // inherited 
    virtual IObject* clone() const;
    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    static CString* create(const string& s);
    const string& stringValue() const;

    static bool equalsIgnoreCase(const string& s1, const string& s2);
    /** like strncmp - does s2 match first n characters of s1 ? */
    static bool equalsIgnoreCase(const string& s1, const string& s2, int n);

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object from description */
    virtual void import(Reader::Node* elem, Reader* reader);
    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Returns a hash code value for the object */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one */
    virtual bool equalTo(const IObject* obj) const;

    /** remove duplicates and empty strings from an array. Order is otherwise
        unchanged */
    static CStringArraySP trim(const CStringArray& names);

    /** Support for conversion of objects from strings. This defines
        the function that the CString class demands of a type.  */
    typedef IObject* (TObjectFromString)(CClassConstSP requiredType,
                                         const string& data);
    /** Register a conversion method for a given type. Note that the
        method will be invoked for derived classes of the supplied
        type */
    static void registerObjectFromStringMethod(CClassConstSP      targetClass,
                                               TObjectFromString* method);

    /** Converts this object to an instance of the requiredType. Throws an
        exception if a conversion to the required Type is not supported */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const;

    /** For use in constructing e.g. error messages */
    virtual string toString() const;

private:
    const string s; // $unregistered
    CString(const string s);
    CString(const CString& rhs);
    CString& operator=(const CString& rhs);
};

typedef smartConstPtr<CString> CStringConstSP;
typedef smartPtr<CString> CStringSP;
#ifndef QLIB_ATOMIC_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<CString>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CString>);
#endif

DRLIB_END_NAMESPACE
#endif

