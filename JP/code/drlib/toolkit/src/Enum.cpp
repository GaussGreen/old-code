
#include "edginc/config.hpp"
#define QLIB_ENUM_CPP
#include "edginc/Class.hpp"
#include "edginc/Writer.hpp"
#include "edginc/Array.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<Enum>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<Enum>);

Enum::~Enum(){}

/** Returns a hash code value for the object */
int Enum::hashCode() const{
    return enumVal;
}

/** Are these objects equal (ie contain the same data) */
bool Enum::equalTo(const IObject* obj) const{
    if (this == obj){
        return true;
    }
    if (obj->getClass() != getClass()){
        return false;
    }
    const Enum* theEnum = STATIC_CAST(Enum, obj);
    return (theEnum->enumVal == enumVal);
}

//// enums are immutable
IObject* Enum::clone() const{
    return const_cast<Enum*>(this);
}

/** Create a wrapped enum of the specified enumType and with specified
    enumVal */
Enum::Enum(CClassConstSP enumType,
           int           enumVal): CObject(enumType), enumVal(enumVal){}

/** Create a wrapped enum of the specified enumType and with no specified
    enumVal */
Enum::Enum(CClassConstSP enumType): CObject(enumType), enumVal(0){}

//// returns the wrapped enum's value
int Enum::enumValueAsInt() const{
    return enumVal;
}

//// returns the wrapped enum's value as a string
const string& Enum::enumValueAsString() const{
    return getClass()->getEnumValue(enumVal).valueAsString;
}
    
/** write object out to writer */
void Enum::write(const string& tag,  Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            writer->write(enumValueAsString());
        }
        writer->objectEnd(tag, this);
    } catch (exception& e){
        throw ModelException(&e, "Enum::write");
    }
}

/** populate an empty object from reader */
void Enum::import(Reader::Node* elem, Reader* reader){
    // import is privileged
    const_cast<int&>(enumVal) = 
        getClass()->getEnumValue(elem->value()).valueAsInt;
}

//// Writes prefix: <type name>::<Enum as string>
void Enum::outputWrite(const string& linePrefix,
                       const string& prefix, ostream& stream) const{
    const string& s  = enumValueAsString();
    stream << linePrefix << prefix << ": " << getClass()->getName() << 
        "::" << s << endl;
}

//// Conversion method to build appropriate Enum from supplied string
IObject* Enum::create(CClassConstSP requiredType,
                      const string& data){
    // map 'data' to integer enum value
    int enumValue = requiredType->getEnumValue(data).valueAsInt;
    // build empty instance of class
    IObject* enumObject = requiredType->newInstance();
    // set the enum value
    const_cast<int&>(STATIC_CAST(Enum, enumObject)->enumVal) = enumValue;
    return enumObject;
}

/** Returns string corresponding to enum integer */
string Enum::toString() const{
    return enumValueAsString();
}

//// invoked when class is loaded. Records 'reflection' data
void Enum::load(CClassSP& classToLoad){
    classToLoad->setPublic(); // make visible to EAS/spreadsheet
    classToLoad->setNative();
    REGISTER(Enum, classToLoad);
    SUPERCLASS(CObject);
    classToLoad->setIsEnum();
    CString::registerObjectFromStringMethod(TYPE, create);
};

CClassConstSP const Enum::TYPE = CClass::registerClassLoadMethod(
    "Enum", typeid(Enum), load);

//// used by ToolkitLib to force the file to be linked in. 
bool EnumLoad(){
    return true;
}

//// Addin to create an enum of a specific type
class EnumAddin: public CObject{
    string enumType;
    string enumAsString;
public:
    static CClassConstSP const TYPE;

    //// Creates [derived] Enum class with right integer value
    IObjectSP create(){
        CClassConstSP clazz = CClass::forName(enumType);
        return IObjectSP(Enum::create(clazz, enumAsString));
    }

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EnumAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(enumType, "Type of the enum");
        FIELD(enumAsString, "Value for the enum");
        Addin::registerObjectMethod("ENUM",
                                    Addin::UTILITIES,
                                    "Creates an Enum of the specified type and "
                                    "corresponding to the supplied value",
                                    true,
                                    Addin::returnHandle,
                                    &EnumAddin::create);
    }
    EnumAddin(): CObject(TYPE){}
    static IObject* defaultConstructor(){
        return new EnumAddin();
    }
};

CClassConstSP const EnumAddin::TYPE = CClass::registerClassLoadMethod(
    "EnumAddin", typeid(EnumAddin), load);

/** Class illustrating/testing the use of enums with reflection */
class EnumTester: public CObject{
public:
    static CClassConstSP const TYPE;
    //// enums are used in their normal way. Repeated values are allowed.
    //// This definition of the enum should go in the header file if other
    //// files want to use it
    enum MyValues {
        FIRST_VALUE = -2,
        SECOND_VALUE = 5,
        THIRD_VALUE = 3
    };

    //// This line is only required if you want an smartPtr for your enum
    typedef smartPtr<BoxedEnum<MyValues> > MyValuesSP;
    //// This line is only required if you want an array of your enums
    typedef array<MyValues> MyValuesArray;

    //// illustration of being able to switch on the value of an enum, and
    //// the ability to convert the enum to its string value
    void doSomething(){
        // the implementation of this method goes in the source file
        switch (myEnum1){
        case FIRST_VALUE:
            // something
        case SECOND_VALUE:
            // something
        case THIRD_VALUE:
            throw ModelException("EnumTester::doSomething",
                                 BoxedEnum<MyValues>::toString(THIRD_VALUE)+
                                 " is not valid here");
            break;
        default:
            // something else
            break;
        }
    }
    
private:
    //// optional enums should be initialised in the constructor
    EnumTester(): CObject(TYPE), myEnum1(FIRST_VALUE){
        // the implementation of this method goes in the source file
    }

    static IObject* defaultConstructor(){
        // the implementation of this method goes in the source file
        return new EnumTester();
    }

    static void load(CClassSP& clazz){
        // the implementation of this method goes in the source file
        clazz->setPublic();
        REGISTER(EnumTester, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        // enums registered using the FIELD macro ...
        FIELD(myEnum1, "one of the allowed values");
        FIELD(myEnum2, "one of the allowed values or nothing");
        FIELD_MAKE_OPTIONAL(myEnum2);
        FIELD(myEnums, "array of the allowed values");
        //// an addin function that builds this EnumTester object
        Addin::registerConstructor("ENUM_TESTER",
                                   Addin::XL_TESTS,
                                   "Tests QLIB Enum functionality",
                                   TYPE);
    }

    /// fields
    MyValues      myEnum1;
    MyValuesSP    myEnum2; /* you can have an SP if you want but not terribly
                              useful */
    MyValuesArray myEnums;
};

/** This goes in the header file. It is only needed if you want an array of
    your enums. In particular, this follows the usual paradigm of defining
    arrays in that this template specialisation is only need if you don't want
    an array of smart pointers. */
template <> class TOOLKIT_DLL arrayObjectCast<EnumTester::MyValues>{
public:
    /** Casts array element to an IObject */
    static IObjectSP toIObject(EnumTester::MyValues value);

    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    static EnumTester::MyValues fromIObject(IObjectSP value);
};

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<EnumTester::MyValues>::toIObject(
    EnumTester::MyValues value){
    return BoxedEnum<EnumTester::MyValues>::create(value);
}

/** Casts array element from an IObject */
EnumTester::MyValues arrayObjectCast<EnumTester::MyValues>::fromIObject(
    IObjectSP value)
{
    BoxedEnum<EnumTester::MyValues>* ptr =
        DYNAMIC_CAST(BoxedEnum<EnumTester::MyValues>, value.get());
    if (!ptr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an EnumTester::MyValues");
    }
    return ptr->enumValue();
}

//// this goes in the .cpp file
CClassConstSP const EnumTester::TYPE = CClass::registerClassLoadMethod(
    "EnumTester", typeid(EnumTester), load);

//// This goes in the .cpp file too
//// Start by naming the enum type together with a description of what the
//// enum represents
START_PUBLIC_ENUM_DEFINITION(EnumTester::MyValues, "Values for MyValues");
//// then list each enum giving its C++ value, then its string value and
//// finally a description of that value
ENUM_VALUE_AND_NAME(EnumTester::FIRST_VALUE, "FIRST_VALUE", 
                    "Indicates that the 1st value is to be used");
ENUM_VALUE_AND_NAME(EnumTester::SECOND_VALUE, "SECOND_VALUE",
                    "Indicates that the 2nd value is to be used");
ENUM_VALUE_AND_NAME(EnumTester::THIRD_VALUE, "THIRD_VALUE",
                    "Indicates that the 3rd value is to be used");
//// finally complete with definition with this macro
END_ENUM_DEFINITION(EnumTester::MyValues);

//// Only needed if you want an array of your enums. This goes in the .cpp 
//// file too
typedef EnumTester::MyValuesArray EnumTesterMyValuesArray; // MSVC7 bug
DEFINE_TEMPLATE_TYPE_WITH_NAME("EnumTester::MyValuesArray",
                               EnumTesterMyValuesArray);

DRLIB_END_NAMESPACE
