

#include "edginc/config.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Void.hpp"


/** A Field provides information about, and dynamic access to, a
    single field of a class or an interface. The reflected field may
    be a class (static) field or an instance field. */

DRLIB_BEGIN_NAMESPACE
/** creates a field using the data given */
CField* CField::create(CClassSP         declaringClass,
                       const string&    name,
                       const type_info& typeId,  // of field
                       TFieldGetMethod* getMethod,
                       TFieldSetMethod* setMethod,
                       PointerAttribute pointerAttribute,
                       int              offset){ // from start of structure
    CFieldSP field = 0;
    try{
        field = new CField();
        field->declaringClass = declaringClass;
        field->name = name;
        field->type = CClass::forName(typeId);
        field->getMethod = getMethod;
        field->setMethod = setMethod;
        field->typePrimitive = field->type->isPrimitive();
        field->typeArray = field->type->isArray();
        field->offset = offset;
        field->pointerAttribute = pointerAttribute;
        // fields of type 'Void' are artefacts of C++ and templates.
        // We flag them as transient as we don't want them in the public
        // interface
        if (field->type == Void::TYPE){
            field->setTransient(true);
        }
    } catch (exception& e){
        throw ModelException(&e, "CField::create",
                             "Failed to load create field: " +
                             name+" (type: "+typeId.name()+")");
    }
    return field;
}

/** gives a pointer to where the data starts. Routine is used for error
    messages */
void* CField::dataAddress(const IObject* obj,
                          const string&  routine) const{
    if (!obj){
        throw ModelException(routine, "NULL Object. Expected type "+
                             declaringClass->getName());
    }
    if (!declaringClass->isInstance(obj)){
        throw ModelException(routine,
                             "Object ("+obj->getClass()->getName()+
                             ") is not of the right type ("+
                             declaringClass->getName()+
                             ") for field "+getName());
    }
    // static cast as we've just checked the type is ok
    void* structStart = const_cast<void*>(declaringClass->staticCast(obj));
    void *dataStart = (char *)structStart + offset;
    return dataStart;
}

/** Returns the value of the field represented by this Field,
    on the specified object. */
IObjectConstSP CField::constGet(const IObjectConstSP& obj) const{
    return getObj(obj.get());
}

/** Returns the value of the field represented by this Field,
    on the specified object. */
IObjectSP CField::get(const IObjectSP& obj) const{
    return getObj(obj.get());
}

/** Returns the value of the field represented by this Field,
    on the specified object. */
IObjectSP CField::getObj(const IObject* obj) const{
    const static string routine = "CField::getObj";
    void *dataStart = dataAddress(obj, routine);
    return getMethod(dataStart);
}

/** Gets the value of a field as a boolean on the specified object. */
bool CField::getBool(const IObjectConstSP& obj) const{
    const static string routine = "CField::getBool";
    if (type != CBool::TYPE){
        throw ModelException(routine,
                             "Field is a "+type->getName()+", not a bool");
    }
    const void *dataStart = dataAddress(obj.get(), routine);
    if (pointerAttribute == INLINE){
        return *(bool *)dataStart;
    } else if (pointerAttribute == PLAIN_POINTER){
        // field is a pointer to a bool
        bool** b = (bool **)dataStart;
        if (!*b){
            throw ModelException(routine,
                                 "Field contains a bool* which is null");
        }
        return **b;
    } else {
        // field is a smart pointer
        CBoolSP* b = (CBoolSP*)dataStart;
        return (*b)->boolValue();
    }
}

/** Gets the value of a field as a double on the specified object. */
double CField::getDouble(const IObjectConstSP& obj) const{
    const static string routine = "CField::getDouble";
    if (type != CDouble::TYPE){
        throw ModelException(routine,
                             "Field is a "+type->getName()+", not a double");
    }
    const void *dataStart = dataAddress(obj.get(), routine);
    if (pointerAttribute == INLINE){
        return *(double *)dataStart;
    } else if (pointerAttribute == PLAIN_POINTER){
        // field is a pointer to a double
        double** db = (double **)dataStart;
        if (!*db){
            throw ModelException(routine,
                                 "Field contains a double* which is null");
        }
        return **db;
    } else {
        // field is a smart pointer
        CDoubleSP* d = (CDoubleSP*)dataStart;
        return (*d)->doubleValue();
    }
}

/** Gets the value of a field as an int on the specified object. */
int CField::getInt(const IObjectConstSP& obj) const{
    const static string routine = "CField::getInt";
    if (type != CInt::TYPE){
        throw ModelException(routine,
                             "Field is a "+type->getName()+", not an int");
    }
    const void *dataStart = dataAddress(obj.get(), routine);
    if (pointerAttribute == INLINE){
        return *(int *)dataStart;
    } else if (pointerAttribute == PLAIN_POINTER){
        // field is a pointer to an int
        int** i = (int **)dataStart;
        if (!*i){
            throw ModelException(routine,
                                 "Field contains a int* which is null");
        }
        return **i;
    } else {
        // field is a smart pointer
        CIntSP* i = (CIntSP*)dataStart;
        return (*i)->intValue();
    }
}

/** Gets the value of a field as a string on the specified object. */
const string& CField::getString(const IObjectConstSP& obj) const{
    const static string routine = "CField::getString";
    if (type != CString::TYPE){
        throw ModelException(routine,
                             "Field is a "+type->getName()+", not a string");
    }
    const void *dataStart = dataAddress(obj.get(), routine);
    if (pointerAttribute == INLINE){
        return *(string *)dataStart;
    } else if (pointerAttribute == PLAIN_POINTER){
        // field is a pointer to a string
        string** s = (string**)dataStart;
        if (!*s){
            throw ModelException(routine,
                                 "Field contains a string* which is null");
        }
        return **s;
    } else {
        // field is a smart pointer
        CStringSP* s = (CStringSP*)dataStart;
        return (*s)->stringValue();
    }
}

/** Sets the value of the field represented by this Field,
    on the specified object.  */
void CField::set(const IObjectSP& obj, const IObjectSP& valueSP) const{
    const static string routine = "CField::set";
    if ((valueSP.get() && !type->isInstance(valueSP))){
        throw ModelException("CField::set", "Object ("+
                             valueSP->getClass()->getName()+
                             ") not of the right type ("+type->getName()+")");
    }
    void* dataStart = dataAddress(obj.get(), routine);
    setMethod(dataStart, valueSP);
}

/** Sets the value of a field as a boolean on the specified object. */
void CField::setBool(const IObjectSP& obj, bool z) const{
    const static string routine = "CField::setBool";
    const void *dataStart = dataAddress(obj.get(), routine);
    if (type != CBool::TYPE){
        throw ModelException(routine,
                             "Field is a "+type->getName()+", not a bool");
    }
    if (pointerAttribute == INLINE){
        *(bool *)dataStart = z;
    } else if (pointerAttribute == PLAIN_POINTER){
        bool** b = (bool **)dataStart;
        if (!*b){
            *b = new bool;
        }
        **b = z;
    } else {
        // field is a smart pointer
        CBoolSP* b = (CBoolSP*)dataStart;
        b->reset(CBool::create(z));
    }
}

/** Sets the value of a field as a double on the specified object. */
void CField::setDouble(const IObjectSP& obj,  double d) const{
    const static string routine = "CField::setDouble";
    const void *dataStart = dataAddress(obj.get(), routine);
    if (type != CDouble::TYPE){
        throw ModelException(routine,
                             "Field is a "+type->getName()+", not a double");
    }
    if (pointerAttribute == INLINE){
        *(double *)dataStart = d;
    } else if (pointerAttribute == PLAIN_POINTER){
        double** db = (double **)dataStart;
        if (!*db){
            *db = new double;
        }
        **db = d;
    } else {
        // field is a smart pointer
        CDoubleSP* db = (CDoubleSP*)dataStart;
        db->reset(CDouble::create(d));
    }
}

/** Sets the value of a field as an int on the specified object. */
void CField::setInt(const IObjectSP& obj, int i) const{
    const static string routine = "CField::setInt";
    const void *dataStart = dataAddress(obj.get(), routine);
    if (type != CInt::TYPE){
        throw ModelException(routine,
                             "Field is a "+type->getName()+", not an int");
    }
    if (pointerAttribute == INLINE){
        *(int *)dataStart = i;
    } else if (pointerAttribute == PLAIN_POINTER){
        int** intVal = (int **)dataStart;
        if (!*intVal){
            *intVal = new int;
        }
        **intVal = i;
    } else {
        // field is a smart pointer
        CIntSP* intVal = (CIntSP*)dataStart;
        intVal->reset(CInt::create(i));
    }
}

/** Sets the value of a field as a string on the specified object. */
void CField::setString(const IObjectSP& obj, const string& s) const{
    const static string routine = "CField::setString";
    const void *dataStart = dataAddress(obj.get(), routine);
    if (type != CString::TYPE){
        throw ModelException(routine, 
                             "Field is a "+type->getName()+", not a string");
    }
    if (pointerAttribute == INLINE){
        *(string *)dataStart = s;
    } else if (pointerAttribute == PLAIN_POINTER){
        string** strVal = (string **)dataStart;
        if (!*strVal){
            *strVal = new string;
        }
        **strVal = s;
    } else {
        // field is a smart pointer
        CStringSP* strVal = (CStringSP*)dataStart;
        strVal->reset(CString::create(s));
    }
}

/** Returns the name of the field represented by this Field object. */
const string& CField::getName() const throw(){
    return name;
}

/** Returns a Class object that identifies the declared type
        for the field represented by this Field object. */
CClassConstSP CField::getType() const throw(){
    return type;
}

/** Returns the Class object representing the class or
            interface that declares the field represented by this
            Field object. */
CClassConstSP CField::getDeclaringClass() const throw(){
    return declaringClass;
}

/** indicates whether the the type of data held by the field is
        atomic */
bool CField::typeIsPrimitive() const throw(){
    return typePrimitive;
}

/** indicates whether the the type of data held by the field is
        an array */
bool CField::typeIsArray() const throw(){
    return typeArray;
}

/** Returns the pointer attribute for this field ie indicates whether the
    field is a pointer, a smart pointer or 'inline' */
CField::PointerAttribute CField::getPointerAttribute() const{
    return pointerAttribute;
}

CField::CField():CObject(TYPE), modifiers(0), optional(false), 
    transientForIteration(false){}

/** clones a field */
IObject* CField::clone() const{
    return new CField(*this);
}

/** flag the optionality of a field */
void CField::setOptional(bool isOptional){
    this->optional = isOptional;
}

/* is the field optional */
bool CField::isOptional() const{
    return optional;
}

/** flag that the field is transient */
void CField::setTransient(bool isTransient){
    if (isTransient){
        modifiers |= Modifier::TRANSIENT;
    } else {
        modifiers = modifiers & ~Modifier::TRANSIENT;
    }
    transientForIteration = isTransient;
}

/** flag that the field is transient for the purposes of ObjectIteration */
void CField::setTransientForIteration(bool isTransient){
    transientForIteration = isTransient;
}

/** is the field transient for the purposes of ObjectIteration */
bool CField::isTransientForIteration() const{
    return transientForIteration;
}

/** Returns the modifiers for the field represented
    by this Field object, as an integer. (see Modifier class) */
int CField::getModifiers() const{
    return modifiers;
}

/** flag the optionality of a field */
void CField::setDescription(const string& description){
    this->description = description;
}

/* is the field optional */
const string& CField::getDescription() const{
    return description;
}

/** Copies (ie clones) component idenitified by this field from source
    to destination. If component is null no action is
    performed. Method is essntially a get followed by a clone and
    then a set but exists for performance reasons (eg when handling
    native types like doubles) */
void CField::copyComponent(const IObjectConstSP& source, 
                           const IObjectSP&      destination) const{
    bool doneIt;
    if (type->isPrimitive() && pointerAttribute == INLINE /* NB rule if cmpt
                                                             is null */){
        doneIt = true; // most cases
        if (type == CDouble::TYPE){
            setDouble(destination, getDouble(source));
        } else if (type == CInt::TYPE){
            setInt(destination, getInt(source));
        } else if (type == CString::TYPE){
            setString(destination, getString(source));
        } else if (type == CBool::TYPE){
            setBool(destination, getBool(source));
        } else {
            doneIt = false;
        }
    } else {
        doneIt = false;
    }
    
    if (!doneIt){
        // to do at some later stage: see if field is 'in-line' and the
        // type is simple (ie a structure copy = a clone), then can skip clone
        IObjectSP o(getObj(source.get()));
        if (o.get()){
            // can skip clone for inline fields which implicity clone
            bool noClone = pointerAttribute == INLINE &&
                type->doesAssignmentOperatorClone();
            // then clone if required and then set it
            set(destination, noClone? o: IObjectSP(o->clone()));
        }
    }
}

static void myLoad(CClassSP& clazz){
    REGISTER(CField, clazz);
    SUPERCLASS(CObject);
}

CClassConstSP const CField::TYPE = CClass::registerClassLoadMethod(
    "Field", typeid(CField), myLoad);
    
// specialised get/set methods
template <> IObjectSP FieldGetInLine<string>(string* str){
    return IObjectSP(CString::create(*str));
}

template <> void FieldSetInLine<string>(string*   str,
                                        IObjectSP obj){
    if (!obj){
        throw ModelException("FieldSetInLine",
                             "Can't set NULL string");
    }
    const CString* newStr = DYNAMIC_CONST_CAST(CString, obj.get());
    *str = newStr->stringValue();
}

template <> IObjectSP FieldGetInLine<double>(double* db){
    return IObjectSP(CDouble::create(*db));
}

template <> void FieldSetInLine<double>(double* db, IObjectSP  obj){
    if (!obj){
        throw ModelException("FieldSetInLine",
                             "Can't set NULL double");
    }
    const CDouble* newDb = DYNAMIC_CONST_CAST(CDouble, obj.get());
    *db = newDb->doubleValue();
}
template <> IObjectSP FieldGetInLine<int>(int* i){
    return IObjectSP(CInt::create(*i));
}

template <> void FieldSetInLine<int>(int* i, IObjectSP obj){
    if (!obj){
        throw ModelException("FieldSetInLine",
                             "Can't set NULL int");
    }
    const CInt* newI = DYNAMIC_CONST_CAST(CInt, obj.get());
    *i = newI->intValue();
}

template <> IObjectSP FieldGetInLine<bool>(bool* b){
    return IObjectSP(CBool::create(*b));
}

template <> void FieldSetInLine<bool>(bool* b, IObjectSP obj){
    if (!obj){
        throw ModelException("FieldSetInLine",
                             "Can't set NULL bool");
    }
    const CBool* newB = DYNAMIC_CONST_CAST(CBool, obj.get());
    *b = newB->boolValue();
}

template <> IObjectSP FieldGetPlainPtr(string** str){
    if (*str == 0){
        return IObjectSP();
    }
    return IObjectSP(CString::create(**str));
}
    
template <> void FieldSetPlainPtr(string** str, IObjectSP obj){
    if (!obj){
        // set pointer to null - must delete any string there already
        delete *str;
        *str = 0;
    } else {
        // create room for string if not there already
        if (!*str){
            *str = new string;
        }
        const CString* newStr = DYNAMIC_CONST_CAST(CString, obj.get());
        **str = newStr->stringValue();
    }
}

template <> IObjectSP FieldGetPlainPtr(double** db){
    if (*db == 0){
        return IObjectSP();
    }
    return IObjectSP(CDouble::create(**db));
}

template <> void FieldSetPlainPtr(double** db, IObjectSP obj){
    if (!obj){
        // set pointer to null - must delete any double there already
        delete *db;
        *db = 0;
    } else {
        if (!*db){
            *db = new double;
        }
        const CDouble* newDb = DYNAMIC_CONST_CAST(CDouble, obj.get());
        **db = newDb->doubleValue();
    }
}

template <> IObjectSP FieldGetPlainPtr(int** val){
    if (*val == 0){
        return IObjectSP();
    }
    return IObjectSP(CInt::create(**val));
}

template <> void FieldSetPlainPtr(int** val, IObjectSP obj){
    if (!obj){
        // set pointer to null - must delete any double there already
        delete *val;
        *val = 0;
    } else {
        if (!*val){
            *val = new int;
        }
        const CInt* newInt = DYNAMIC_CONST_CAST(CInt, obj.get());
        **val = newInt->intValue();
    }
}

template <> IObjectSP FieldGetPlainPtr(bool** val){
    if (*val == 0){
        return IObjectSP();
    }
    return IObjectSP(CBool::create(**val));
}

template <> void FieldSetPlainPtr(bool** val, IObjectSP obj){
    if (!obj){
        // set pointer to null - must delete any double there already
        delete *val;
        *val = 0;
    } else {
        if (!*val){
            *val = new bool;
        }
        const CBool* newBool = DYNAMIC_CONST_CAST(CBool, obj.get());
        **val = newBool->boolValue();
    }
}


template <> IObjectSP FieldGetInLine(IObjectSP* t){
    return *t;
}

template <> void FieldSetInLine(IObjectSP* t, IObjectSP obj){
    *t = obj;
}    
DRLIB_END_NAMESPACE
