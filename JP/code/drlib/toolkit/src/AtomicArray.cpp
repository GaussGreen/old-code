
#include "edginc/config.hpp"
#define QLIB_ATOMICARRAY_CPP
#include "edginc/Addin.hpp"
#include "edginc/Hashtable.hpp"

DRLIB_BEGIN_NAMESPACE
/** explicity instantiate template if we did an 'extern template' */
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<double _COMMA_ double>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<int _COMMA_ int>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<string _COMMA_ string>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CIntArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CIntArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CStringArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CStringArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<DoubleArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<DoubleArrayArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<StringArraySP _COMMA_ StringArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<IntArraySP _COMMA_ IntArray>);

INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CDoubleArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CDoubleArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CBoolArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CBoolArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<IntArrayArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<IntArrayArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<StringArrayArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<StringArrayArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<DoubleArrayArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<DoubleArrayArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<DoubleArrayArrayArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<DoubleArrayArrayArray>);
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<DoubleArray>(DoubleArray*  t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<IntArray>(IntArray*  t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<StringArray>(StringArray*  t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<BoolArray>(BoolArray*  t));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<DoubleArray>(DoubleArray* t _COMMA_
                                                      IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<IntArray>(IntArray* t
                                                   _COMMA_ IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<StringArray>(StringArray* t  _COMMA_
                                                      IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<BoolArray>(BoolArray* t  _COMMA_ 
                                                    IObjectSP o));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<CDoubleArraySP>(
                         CDoubleArraySP* t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<CIntArraySP>(
                         CIntArraySP* t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<CStringArraySP>(
                         CStringArraySP* t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<CBoolArraySP>(CBoolArraySP* t));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<CDoubleArraySP>(CDoubleArraySP* t,
                                                           IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<CIntArraySP>(CIntArraySP* t,
                                                        IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<CStringArraySP>(CStringArraySP* t,
                                                           IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<CBoolArraySP>(CBoolArraySP* t, 
                                                        IObjectSP o));

//// msvc6 doesn't like it if these are after code which causes the template
//// to be instantiated
//template<> CClassConstSP const CDoubleArray::TYPE = 
//CClass::registerClassLoadMethod("DoubleArray", typeid(CDoubleArray), load);
DEFINE_TEMPLATE_TYPE_WITH_NAME("DoubleArray", CDoubleArray);

//template<> CClassConstSP const DoubleArrayArray::TYPE = 
//CClass::registerClassLoadMethod(
//    "DoubleArrayArray", typeid(DoubleArrayArray), load);
DEFINE_TEMPLATE_TYPE(DoubleArrayArray);

//template<> CClassConstSP const CStringArray::TYPE =
// CClass::registerClassLoadMethod("StringArray", typeid(CStringArray), load);
DEFINE_TEMPLATE_TYPE_WITH_NAME("StringArray", CStringArray);

//template<> CClassConstSP const CIntArray::TYPE = 
//CClass::registerClassLoadMethod(
//    "IntArray", typeid(CIntArray), load);
DEFINE_TEMPLATE_TYPE_WITH_NAME("IntArray", CIntArray);

/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
/** Casts array element to an IObject */
IObjectSP arrayObjectCast<double>::toIObject(double value){
    IObjectSP objValue = IObjectSP(CDouble::create(value));
    return objValue;
}

/** Sets the value of the indexed component of the specified
    array object to the specified new value. */
double arrayObjectCast<double>::fromIObject(IObjectSP value){
    CDouble *dbPtr = DYNAMIC_CAST(CDouble, value.get());
    if (!dbPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a Double");
    }
    return dbPtr->doubleValue();
}

IObject* arrayClone<double>::clone(const CArray* arrayToClone){
    const CDoubleArray& dbArray = 
        static_cast<const CDoubleArray&>(*arrayToClone);
    return new CDoubleArray(dbArray);
}

//// hash code for array of doubles
int arrayCompare<double>::hashCode(const CArray* arrayToHash){
    const CDoubleArray& dbArray =
        static_cast<const CDoubleArray&>(*arrayToHash);
    int hCode = (size_t) CDoubleArray::TYPE;
    for (int i = 0; i < dbArray.size(); i++){
        hCode ^= CDouble::hashCode(dbArray[i]);
    }
    return hCode;
}
//// equals method for array of doubles
bool arrayCompare<double>::equalTo(const CArray*  arrayToCompare, 
                                   const IObject* obj){
    if (arrayToCompare == obj){
        return true;
    }
    if (!obj || obj->getClass() != CDoubleArray::TYPE){
        return false;
    }
    const CDoubleArray& dbArray1 =
        static_cast<const CDoubleArray&>(*arrayToCompare);
    const CDoubleArray& dbArray2 = *STATIC_CAST(CDoubleArray, obj);
    return (dbArray1.equals(dbArray2));
}

/** specialisations of arrayObjectCast for DoubleArrayArray */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<DoubleArray>::toIObject(
    const DoubleArray& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}


/** Casts array element to an IObject */
IObjectSP arrayObjectCast<DoubleArray>::toIObject(DoubleArray& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a DoubleArray */
const DoubleArray& arrayObjectCast<DoubleArray>::fromIObject(IObjectSP& value){
    DoubleArray *dtPtr = DYNAMIC_CAST(DoubleArray, value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a DoubleArray");
    }
    return *dtPtr;
}

// explicit clone for arrays of Double arrays - for performance
IObject* arrayClone<DoubleArray>::clone(const CArray* arrayToClone){
    const DoubleArrayArray& theCluster = 
        static_cast<const DoubleArrayArray&>(*arrayToClone);
    return new DoubleArrayArray(theCluster);
}

/** specialisations of arrayObjectCast for DoubleArrayArrayArray */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<DoubleArrayArray>::toIObject(
    const DoubleArrayArray& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<DoubleArrayArray>::toIObject(DoubleArrayArray& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a DoubleArray */
const DoubleArrayArray& arrayObjectCast<DoubleArrayArray>::fromIObject(IObjectSP& value){
    DoubleArrayArray *dtPtr = DYNAMIC_CAST(DoubleArrayArray, value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a DoubleArrayArray");
    }
    return *dtPtr;
}

// explicit clone for arrays of 3D Double Arrays - for performance
IObject* arrayClone<DoubleArrayArray>::clone(const CArray* arrayToClone){
    const DoubleArrayArrayArray& theCluster = 
        static_cast<const DoubleArrayArrayArray&>(*arrayToClone);
    return new DoubleArrayArrayArray(theCluster);
}

//template<> CClassConstSP const DoubleArrayArrayArray::TYPE = 
//CClass::registerClassLoadMethod(
//    "DoubleArrayArrayArray", typeid(DoubleArrayArrayArray), load);
DEFINE_TEMPLATE_TYPE(DoubleArrayArrayArray);

//template<> CClassConstSP const IntArrayArray::TYPE = 
//CClass::registerClassLoadMethod(
//   "IntArrayArray", typeid(IntArrayArray), load);
DEFINE_TEMPLATE_TYPE(IntArrayArray);

//template<> CClassConstSP const StringArrayArray::TYPE = 
//CClass::registerClassLoadMethod(
//    "StringArrayArray", typeid(StringArrayArray), load);
DEFINE_TEMPLATE_TYPE(StringArrayArray);

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<int>::toIObject(int value){
    IObjectSP objValue = IObjectSP(CInt::create(value));
    return objValue;
}


/** Sets the value of the indexed component of the specified
        array object to the specified new value. */
int arrayObjectCast<int>::fromIObject(IObjectSP value){
    CInt *intPtr = DYNAMIC_CAST(CInt, value.get());
    if (!intPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an Int");
    }
    return intPtr->intValue();
}

IObject* arrayClone<int>::clone(const CArray* arrayToClone){
    const CIntArray& theArray = 
        static_cast<const CIntArray&>(*arrayToClone);
    return new CIntArray(theArray);
}

//// hash code for array of ints
int arrayCompare<int>::hashCode(const CArray* arrayToHash){
    const CIntArray& intArray =
        static_cast<const CIntArray&>(*arrayToHash);
    int hCode = (size_t) CIntArray::TYPE;
    for (int i = 0; i < intArray.size(); i++){
        hCode ^= intArray[i];
    }
    return hCode;
}
//// equals method for array of ints
bool arrayCompare<int>::equalTo(const CArray*  arrayToCompare, 
                                const IObject* obj){
    if (arrayToCompare == obj){
        return true;
    }
    if (!obj || obj->getClass() != CIntArray::TYPE){
        return false;
    }
    const CIntArray& theArray1 =
        static_cast<const CIntArray&>(*arrayToCompare);
    const CIntArray& theArray2 = *STATIC_CAST(CIntArray, obj);
    return (theArray1.equals(theArray2));
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<bool>::toIObject(bool value){
    IObjectSP objValue = IObjectSP(CBool::create(value));
    return objValue;
}

/** Sets the value of the indexed component of the specified
        array object to the specified new value. */
bool arrayObjectCast<bool>::fromIObject(IObjectSP value){
    CBool *boolPtr = DYNAMIC_CAST(CBool, value.get());
    if (!boolPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an Bool");
    }
    return boolPtr->boolValue();
}

//// hash code for array of bools
int arrayCompare<bool>::hashCode(const CArray* arrayToHash){
    const CBoolArray& boolArray =
        static_cast<const CBoolArray&>(*arrayToHash);
    int hCode = (size_t) CBoolArray::TYPE;
    int bitPos = 0; // ranges from 0 to 32 (number of bits in an int)
    for (int i = 0; i < boolArray.size(); i++){
        bool b = boolArray[i];
        if (b){
            hCode ^= 1 << bitPos; // = 2^bitPos
        }
        bitPos++;
        if (bitPos == sizeof(int)*8){ // 8 bits in a byte
            bitPos = 0;
        }
    }
    return hCode;
}

//// equals method for array of bools
bool arrayCompare<bool>::equalTo(const CArray*  arrayToCompare, 
                                 const IObject* obj){
    if (arrayToCompare == obj){
        return true;
    }
    if (!obj || obj->getClass() != CBoolArray::TYPE){
        return false;
    }
    const CBoolArray& theArray1 =
        static_cast<const CBoolArray&>(*arrayToCompare);
    const CBoolArray& theArray2 = *STATIC_CAST(CBoolArray, obj);
    return (theArray1.equals(theArray2));
}

CClassConstSP const CBoolArray::TYPE = 
CClass::registerClassLoadMethod("BoolArray", typeid(CBoolArray), load);

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<string>::toIObject(const string& value){
    IObjectSP objValue = IObjectSP(CString::create(value));
    return objValue;
}

/** Sets the value of the indexed component of the specified
        array object to the specified new value. */
string arrayObjectCast<string>::fromIObject(IObjectSP value){
    CString *sPtr = DYNAMIC_CAST(CString, value.get());
    if (!sPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a string");
    }
    return sPtr->stringValue();
}
IObject* arrayClone<string>::clone(const CArray* arrayToClone){
    const CStringArray& theArray = 
        static_cast<const CStringArray&>(*arrayToClone);
    return new CStringArray(theArray);
}

//// hash code for array of strings
int arrayCompare<string>::hashCode(const CArray* arrayToHash){
    const CStringArray& stringArray =
        static_cast<const CStringArray&>(*arrayToHash);
    int hCode = (size_t) CStringArray::TYPE;
    for (int i = 0; i < stringArray.size(); i++){
        hCode ^= hash_string(stringArray[i]);
    }
    return hCode;
}
//// equals method for array of strings
bool arrayCompare<string>::equalTo(const CArray*  arrayToCompare, 
                                 const IObject* obj){
    if (arrayToCompare == obj){
        return true;
    }
    if (!obj || obj->getClass() != CStringArray::TYPE){
        return false;
    }
    const CStringArray& theArray1 =
        static_cast<const CStringArray&>(*arrayToCompare);
    const CStringArray& theArray2 = *STATIC_CAST(CStringArray, obj);
    return (theArray1.equals(theArray2));
}


/** Addin 'constructor' functions for arrays */

/** Addin for building handles to double arrays */
class DoubleArrayAddin: public CObject{
    static CClassConstSP const TYPE;

    /** the single [optional] parameter that our addin takes */
    CDoubleArraySP  dbArray;

    /** the 'addin function' - returns a DoubleArray (creates empty one if
        no data supplied) */
    static IObjectSP createArray(DoubleArrayAddin* params){
        if (!params->dbArray){
            return IObjectSP(new DoubleArray());
        }
        return params->dbArray;
    }

    /** for reflection */
    DoubleArrayAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DoubleArrayAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDoubleArrayAddin);
        FIELD(dbArray, "Array of doubles");
        FIELD_MAKE_OPTIONAL(dbArray);
        Addin::registerClassObjectMethod("DOUBLE_ARRAY",
                                         Addin::UTILITIES,
                                         "Constructs a handle to an array of "
                                         "doubles",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createArray);
    }

    static IObject* defaultDoubleArrayAddin(){
        return new DoubleArrayAddin();
    }
    
};

CClassConstSP const DoubleArrayAddin::TYPE = CClass::registerClassLoadMethod(
    "DoubleArrayAddin", typeid(DoubleArrayAddin), load);


/** Addin for building handles to int arrays */
class IntArrayAddin: public CObject{
    static CClassConstSP const TYPE;

    /** the single parameter that our addin takes */
    CIntArraySP  intArray;

    /** the 'addin function' - returns an IntArray (creates empty one if
        no data supplied) */
    static IObjectSP createArray(IntArrayAddin* params){
        if (!params->intArray){
            return IObjectSP(new IntArray());
        }
        return params->intArray;
    }

    /** for reflection */
    IntArrayAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IntArrayAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultIntArrayAddin);
        FIELD(intArray, "Array of integers");
        FIELD_MAKE_OPTIONAL(intArray);
        Addin::registerClassObjectMethod("INT_ARRAY",
                                         Addin::UTILITIES,
                                         "Constructs a handle to an array of "
                                         "integers",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createArray);
    }

    static IObject* defaultIntArrayAddin(){
        return new IntArrayAddin();
    }
    
};

CClassConstSP const IntArrayAddin::TYPE = CClass::registerClassLoadMethod(
    "IntArrayAddin", typeid(IntArrayAddin), load);

/** Addin for building handles to bool arrays */
class BoolArrayAddin: public CObject{
    static CClassConstSP const TYPE;

    /** the single parameter that our addin takes */
    CBoolArraySP  boolArray;

    /** the 'addin function' - returns a BoolArray (creates empty one if
        no data supplied) */
    static IObjectSP createArray(BoolArrayAddin* params){
        if (!params->boolArray){
            return IObjectSP(new BoolArray());
        }
        return params->boolArray;
    }

    /** for reflection */
    BoolArrayAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BoolArrayAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBoolArrayAddin);
        FIELD(boolArray, "Array of booleans");
        FIELD_MAKE_OPTIONAL(boolArray);
        Addin::registerClassObjectMethod("BOOLEAN_ARRAY",
                                         Addin::UTILITIES,
                                         "Constructs a handle to an array of "
                                         "booleans",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createArray);
    }

    static IObject* defaultBoolArrayAddin(){
        return new BoolArrayAddin();
    }
    
};

CClassConstSP const BoolArrayAddin::TYPE = CClass::registerClassLoadMethod(
    "BoolArrayAddin", typeid(BoolArrayAddin), load);

/** Addin for building handles to string arrays */
class StringArrayAddin: public CObject{
    static CClassConstSP const TYPE;

    /** the single parameter that our addin takes */
    CStringArraySP  stringArray;

    /** the 'addin function' - returns a StringArray (creates empty one if
        no data supplied) */
    static IObjectSP createArray(StringArrayAddin* params){
        if (!params->stringArray){
            return IObjectSP(new StringArray());
        }
        return params->stringArray;
    }

    /** for reflection */
    StringArrayAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(StringArrayAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultStringArrayAddin);
        FIELD(stringArray, "Array of strings");
        FIELD_MAKE_OPTIONAL(stringArray);
        Addin::registerClassObjectMethod("STRING_ARRAY",
                                         Addin::UTILITIES,
                                         "Constructs a handle to an array of "
                                         "strings",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createArray);
    }

    static IObject* defaultStringArrayAddin(){
        return new StringArrayAddin();
    }
    
};

CClassConstSP const StringArrayAddin::TYPE = CClass::registerClassLoadMethod(
    "StringArrayAddin", typeid(StringArrayAddin), load);


///////////////////////////////////////////////////////////////////////////////
DoubleArraySP operator+(const DoubleArray & g1, const DoubleArray & g2)
{
    static const string routine 
        = "DoubleArraySP operator+(DoubleArray, DoubleArray)";
    try 
    { 
        const int _size = g1.size();
        DoubleArraySP _result(new DoubleArray(_size));
        if (_size != g2.size())
        {
            throw ModelException("Arrays sizes do not match", routine);
        }        
        int s; for(s = 0; s < _size; s++)
        {  
            (*_result)[s] = g1[s] + g2[s];
        }
        return _result;
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
void operator +=(DoubleArray & g1, const DoubleArray & g2)
{
    static const string routine 
        = "void operator+=(DoubleArray, DoubleArray)";
    try 
    { 
        const int _size = g1.size();
        if (_size != g2.size())
        {
            throw ModelException("Arrays sizes do not match", routine);
        }        
        int s; for(s = 0; s < _size; s++)
        {  
            g1[s] += g2[s];
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
DoubleArraySP operator*(double gDouble, const DoubleArray & gArray)
{
    static const string routine="DoubleArraySP operator*(double, DoubleArray)";
    try 
    { 
        const int _size = gArray.size();
        DoubleArraySP _result(new DoubleArray(_size));
        int s; for(s = 0; s < _size; s++)
        {  
            (*_result)[s] = gDouble*gArray[s];
        }
        return _result;
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
DoubleArraySP operator*(const DoubleArray & gArray, double gDouble)
{
    static const string routine="DoubleArraySP operator*(DoubleArray, double)";
    try 
    { 
        const int _size = gArray.size();
        DoubleArraySP _result(new DoubleArray(_size));
        int s; for(s = 0; s < _size; s++)
        {  
            (*_result)[s] = gDouble*gArray[s];
        }
        return _result;
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
void operator *=(DoubleArray & gArray, double gDouble)
{
    static const string routine = "void operator*=(DoubleArray, double)";
    try 
    { 
        const int _size = gArray.size();
        int s; for(s = 0; s < _size; s++)
        {  
            gArray[s] *= gDouble;
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE

