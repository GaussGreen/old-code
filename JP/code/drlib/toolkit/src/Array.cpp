
#include "edginc/config.hpp"
#define QLIB_ARRAY_CPP
#include "edginc/Writer.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Null.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Format.hpp"

//   Provides wrapper around stl vector template to give object like
//   functionality
DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<IArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<IArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<IObjectSP _COMMA_ IObject>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<ObjectArray>);

void IArray::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER_INTERFACE(IArray, clazz);
    EXTENDS(IObject);
}

/** XML tag for recording an array's length */
const string IArray::ARRAY_LENGTH = "length";


CClassConstSP const IArray::TYPE = 
CClass::registerInterfaceLoadMethod("IArray", typeid(IArray), load);


// CArray::CArray(const CClassConstSP& objClass) in header file for performance

/** override default field wise copy */
IObject* CArray::clone() const{
    return copyArray(this);
}

/** generic clone operation (typically used for cloning arrays of pointers) */
IObject* CArray::copyArray(const CArray* arrayToClone){
    try{
        int length = arrayToClone->getLength();
        IArraySP newArray(arrayToClone->getClass()->newArrayInstance(length));
        for (int idx = 0; idx < length; idx++){
            IObjectConstSP  origCmpt(arrayToClone->get(idx));
            // skip over null elements in list
            if (origCmpt.get()){
                IObjectSP  newCmpt(origCmpt.clone());
                newArray->set(idx, newCmpt);
            }
        }
        return newArray.release();
    } catch (exception& e){
        throw ModelException(&e, "CArray::copyArray");
    }
}

/** Returns a hash code value for the object by XORing the hashCode for
    each element in the array */
int CArray::hashCode() const{
    return hashArray(this);
}

/** Indicates whether some other object is "equal to" this one by comparing
    each element in the array */
bool CArray::equalTo(const IObject* obj) const{
    return equalToArray(this, obj);
}

/** default method for computing hashCode for an array */
int CArray::hashArray(const CArray* arrayToHash){
    int hCode = (size_t) arrayToHash->getClass();
    int len = arrayToHash->getLength();
    hCode ^= len;
    for (int i = 0; i < len; i++){
        IObjectConstSP obj(arrayToHash->get(i));
        if (obj.get()){
            hCode ^= obj->hashCode();
        }
    }
    return hCode;
}

/** default method for comparing two arrays */
bool CArray::equalToArray(const CArray* arrayToCompare, const IObject* obj){
    if (arrayToCompare == obj){
        return true; // the obvious comparison
    }
    CClassConstSP  c;
    if (!obj || (c = obj->getClass()) != arrayToCompare->getClass()){
        return false; // if obj null or a different type
    }
    int len1 = arrayToCompare->getLength();
    const CArray* array2 = STATIC_CAST(CArray, obj);
    int len2 = array2->getLength();
    if (len1 != len2){
        return false;
    }
    for (int i = 0; i < len1; i++){
        IObjectConstSP o1(arrayToCompare->get(i));
        IObjectConstSP o2(array2->get(i));
        if ((!o1 && o2.get() != 0) || (!o2 && o1.get())){
            return false;
        }
        // recurse
        if (o1.get() && !o1->equalTo(o2.get())){
            return false;
        }
    }
    return true;
}

/** write object out to Writer */
void CArray::write(const string& tag, Writer* writer) const{
    try{
        int length = getLength();
        char buffer[256];
        sprintf(buffer, "%s='%d'", ARRAY_LENGTH.c_str(), length);
        IObjectConstSP obj(writer->objectStart(tag, buffer, this, true));
        if (obj.get()){
            for (int i = 0; i < length; i++){
                IObjectConstSP elem(get(i));
                sprintf(buffer, "%s%d", "Item", i);
                if (!elem){
                    // need to handle null elements in list
                    writer->writeNull(buffer);
                } else {
                    elem->write(buffer, writer);
                }
            }
        }
        writer->objectEnd(tag, this);
    } 
    catch (exception& e){
        throw ModelException(e, "CArray::write");
    }
}

/** populate an empty object from Reader */
void CArray::import(Reader::Node* elem, Reader* reader){
    static const string method("CArray::xmlImport");
    try {
        Reader::NodeListSP  nl(elem->children());
        int arrayLength = getLength();
        int j = 0;
        for (unsigned int i = 0; i < nl->size(); i++) {
            Reader::Node* child = (*nl)[i].get();
            if (j >= arrayLength){
                throw ModelException(method, "Malformed xml - supplied "
                                     "array has too many elements");
            }
            IObjectSP obj(reader->read(child));
            // skip over null elements
            if (!CNull::TYPE->isInstance(obj)){
                CClassConstSP  desiredType = getClass()->getComponentType();
                if (desiredType != obj->getClass()) {
                    CObject::checkType(obj, desiredType);
                }
                set(j, obj);
            }
            j++;
        }
    }
    catch (exception &e) {
        throw ModelException(&e, method);
    }    
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void CArray::outputWrite(const string& linePrefix,
                         const string& prefix, ostream& stream) const{
    int length = getLength();
    char buffer[16];
    for (int i = 0; i < length; i++){
        IObjectConstSP elem(get(i));
        // skip over null elements
        if (elem.get()){
            sprintf(buffer, "[%d]", i);
            elem->outputWrite(linePrefix, prefix+buffer, stream);
        }
    }
}  


void CArray::load(CClassSP& classToLoad){
    REGISTER(CArray, classToLoad);
    SUPERCLASS(CObject);
    IMPLEMENTS(IArray);
    classToLoad->enableCloneOptimisations(); /* allow particular types of arrays
                                                to switch this on */
    classToLoad->setArrayType(typeid(IObject), 0);
}

CClassConstSP const CArray::TYPE = CClass::registerClassLoadMethod(
    "Array", typeid(CArray), load);

/** specialisation of arrayObjectCast for arrays of IObjectSP */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<IObjectSP>::toIObject(const IObjectSP& value){
    return value;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<IObjectSP>::toIObject(IObjectSP& value){
    return value;
}

/** Sets the value of the indexed component of the specified
    array object to the specified new value. */
IObjectSP arrayObjectCast<IObjectSP>::fromIObject(IObjectSP& value){
    return value;
}

//template<> CClassConstSP const ObjectArray::TYPE =
//CClass::registerClassLoadMethod("IObjectArray", typeid(ObjectArray), load);
DEFINE_TEMPLATE_TYPE_WITH_NAME("IObjectArray", ObjectArray);

/** Addin for building handles to arrays of specific types of objects */
class ObjectArrayAddin: public CObject{
    static CClassConstSP const TYPE;

    /** addin takes two parameters - the name of the type of the components
        and an array of the components themselves */
    string         componentType;
    ObjectArraySP  objectArray; // an array of IObjectSP

    /** the 'addin function' - builds array of correct type */
    static IObjectSP createArray(ObjectArrayAddin* params){
        try{
            CClassConstSP clazz = CClass::forName(params->componentType);
            int modifiers = clazz->getModifiers();
            if (!(Modifier::isPublic(modifiers))){
                throw ModelException("Creation of arrays of non public type "
                                     "("+clazz->getName()+") disallowed");
            }
            int length = (!params->objectArray)? 
                0: params->objectArray->size();
            // need to create array based upon its component type
            IArraySP array(clazz->newArrayInstanceByComponent(length));
            for (int i = 0; i < length; i++){
                IObjectSP cmpt((*params->objectArray)[i]);
                if (cmpt.get()){
                    /* check type (we may get IPublicObjects which need to 
                       be converted before being put in the array) */
                    CObject::checkType(cmpt, clazz);
                    array->set(i, cmpt);
                }
            }
            return array;
        } catch (exception& e){
            throw ModelException(e, "ObjectArrayAddin::createArray", 
                                 "Failed for array of type "+
                                 params->componentType);
        }
    }

    /** for reflection */
    ObjectArrayAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObjectArrayAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObjectArrayAddin);
        FIELD(componentType, "Type of component of array");
        FIELD(objectArray, "Array of objects");
        FIELD_MAKE_OPTIONAL(objectArray);
        Addin::registerClassObjectMethod("ARRAY",
                                         Addin::UTILITIES,
                                         "Constructs a handle to an array of "
                                         "objects",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createArray);
    }

    static IObject* defaultObjectArrayAddin(){
        return new ObjectArrayAddin();
    }
    
};

CClassConstSP const ObjectArrayAddin::TYPE = CClass::registerClassLoadMethod(
    "ObjectArrayAddin", typeid(ObjectArrayAddin), load);

/** Addin for merging n arrays into one big one */
class ArrayMergeAddin: public CObject{
    static CClassConstSP const TYPE;

    /** addin takes one parameter - the array of the arrays to merge.
        Actually components can be objects rather than the array */
    ObjectArraySP  arrayArray;

    /** the 'addin function' - builds array of correct type */
    static IObjectSP createArray(ArrayMergeAddin* params){
        static const string method("ArrayMergeAddin::createArray");
        try{
            ObjectArray&  inputArray = *params->arrayArray;
            // first count total number and check component type
            int total = 0;
            CClassConstSP arrayCmptType = 0;
            for (int i = 0; i < inputArray.size(); i++){
                IObjectSP& elt = inputArray[i];
                if (elt.get()){
                    CClassConstSP eltType = elt->getClass();
                    if (eltType->isArray()){
                        CArray& arrayElt = *DYNAMIC_CAST(CArray,elt.get());
                        total += arrayElt.getLength();
                        eltType=eltType->getComponentType();
                    } else {
                        total++;
                    }
                    if (!arrayCmptType){
                        arrayCmptType = eltType;
                    } else if (!arrayCmptType->isAssignableFrom(eltType)){
                        // require other arrays are of same type
                        throw ModelException(
                            method, "Array element "+Format::toString(i+1)+
                            " contains instances of type "+eltType->getName()+
                            ". Instances of type "+arrayCmptType->getName()+
                            " are required.");
                    }
                }
            }
            if (!arrayCmptType){
                throw ModelException(method, "No arrays supplied");
            }
            // need to create array based upon its component type
            IArraySP myArray(arrayCmptType->
                             newArrayInstanceByComponent(total));
            int index = 0;
            for (int j = 0; j < inputArray.size(); j++){
                IObjectSP elt(inputArray[j]);
                if (elt.get()){
                    CClassConstSP eltType = elt->getClass();
                    if (eltType->isArray()){
                        CArray& arrayElt = *DYNAMIC_CAST(CArray,elt.get());
                        for (int k = 0; k < arrayElt.getLength(); k++){
                            IObjectSP value(arrayElt.get(k));
                            if (value.get()){
                                /* check type (we may get IPublicObjects which
                                   need to be converted before being put in 
                                   the array) */
                                CObject::checkType(value, arrayCmptType);
                                myArray->set(index, value);
                            }
                            index++;
                        }
                    } else {
                        /* check type (we may get IPublicObjects which need
                           to be converted before being put in the array) */
                        CObject::checkType(elt, arrayCmptType);
                        myArray->set(index, elt);
                        index++;
                    }
                }
            }
            return myArray;
        } catch (exception& e){
            throw ModelException(e, "ArrayMergeAddin::createArray");
        }
    }

    /** for reflection */
    ArrayMergeAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ArrayMergeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultArrayMergeAddin);
        FIELD(arrayArray, "Array of arrays");
        Addin::registerClassObjectMethod("ARRAY_MERGE",
                                         Addin::UTILITIES,
                                         "Merges an array of arrays",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createArray);
    }

    static IObject* defaultArrayMergeAddin(){
        return new ArrayMergeAddin();
    }
    
};

CClassConstSP const ArrayMergeAddin::TYPE = CClass::registerClassLoadMethod(
    "ArrayMergeAddin", typeid(ArrayMergeAddin), load);

/** Addin for returning length of an array */
class ArrayLengthAddin: public CObject{
    static CClassConstSP const TYPE;

    IArraySP  theArray; // an array

    /** the 'addin function' - builds array of correct type */
    static int getLength(ArrayLengthAddin* params){
        return params->theArray->getLength();
    }

    /** for reflection */
    ArrayLengthAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ArrayLengthAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultArrayLengthAddin);
        FIELD(theArray, "An Array");
        Addin::registerInstanceIntMethod("ARRAY_LENGTH",
                                         Addin::UTILITIES,
                                         "Returns the length of an array",
                                         TYPE,
                                         (Addin::IntMethod*)getLength);
    }

    static IObject* defaultArrayLengthAddin(){
        return new ArrayLengthAddin();
    }
};

CClassConstSP const ArrayLengthAddin::TYPE = CClass::registerClassLoadMethod(
    "ArrayLengthAddin", typeid(ArrayLengthAddin), load);

DRLIB_END_NAMESPACE
