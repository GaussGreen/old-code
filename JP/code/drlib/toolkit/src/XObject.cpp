//----------------------------------------------------------------------------
//
//   Group       : Global Derivatives Research
//
//   Filename    : XObject.cpp
//
//   Description : Wrapper for external [DR interface] objects
//
//   Author      : Mark A Robson
//
//   Date        : 14 Nov 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_XOBJECT_CPP
#include "edginc/AtomicArray.hpp"
#include "edginc/XMap.hpp"
#include "edginc/GDRDate.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DRLibrary.hpp"
#include "edginc/Handle.hpp"
#include "edginc/Writer.hpp"
#include "edginc/XArray.hpp"
#include "edginc/XMatrix.hpp"
#include "edginc/XClass.hpp"
#include "edginc/DRUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// DRI_CHECK is defined in DRUtil.hpp. This simplifies the invocation of the
// DRI_CHECK macro by fixing some parameters.
#define CHECK(exec, method) \
    DRI_CHECK((exec), (method), 0, driErrorHandler, 0, ;, ;, ;)

// DRI_ERROR_CALLBACK (DRUtil.hpp)
void XObject::driErrorHandler(const char *method,
                              const char *errMsg,
                              void       *cbParam)
{
    throwError(method, errMsg);
}

// xml attribute name for getXClassName()
const string XObject::XTYPE_ATTRIBUTE = "XTYPE";

//// Helper for serialisation
class XObject::XWriter: public Writer {
public:
    virtual ~XWriter(){}
    //// constructor
    XWriter(Writer* writer, const DRLibrarySP& drLibrary):
        writer(writer), drLibrary(drLibrary){}

    /** All redirected to component writer */
    virtual void documentStart(const string& id){
        writer->documentStart(id);
    }
    virtual void documentEnd(const string& id){
        writer->documentEnd(id);
    }
    virtual IObjectConstSP objectStart(const string&  id,
                                       const string&  attributes,
                                       const IObject* object,
                                       bool           appendNewLine){
        return writer->objectStart(id, attributes, object, appendNewLine);
    }
    virtual void objectEnd(const string& id, const IObject* object){
        writer->objectEnd(id, object);
    }
    virtual void commentStart(){
        writer->commentStart();
    }
    virtual void commentEnd(){
        writer->commentEnd();
    }
    virtual void write(const string& data){
        writer->write(data);
    }
    virtual void writeNull(const string& id){
        writer->writeNull(id);
    }
    virtual void outputWrite(const string&  linePrefix,
                             const string&  prefix,
                             const IObject* object) const{
        writer->outputWrite(linePrefix, prefix, object);
    }

    //// returns the wrapped writer
    Writer* getWriter(){
        return writer;
    }
    //// returns the associated DRLibrary object
    DRLibrarySP getDRLibrary(){
        return drLibrary;
    }

private:
    Writer*     writer;
    DRLibrarySP drLibrary;
};

//// Helper for deserialisation
class XObject::XReader: public Reader{
public:
    virtual ~XReader(){}

    //// constructor. Pointer to DRLibrarySP must not go out of scope
    XReader(Reader* reader, const DRLibrarySP& drLibrary):
        reader(reader), drLibrary(drLibrary){}


    /** All Redirected to component reader */
    virtual IObjectSP read(){
        return reader->read();
    }
    virtual Node* root(){
        return reader->root();
    }

    /** Redirected to component reader->read(Node*, this) */
    virtual smartPtr<IObject> read(Node* node, Reader* subNodeReader){
        return reader->read(node, this);
    }

    //// returns the DRLibrary associated with this reader
    DRLibrarySP getDRLibrary(){
        return drLibrary;
    }

    //// returns the wrapped reader
    Reader* getReader(){
        return reader;
    };
private:
    Reader*      reader;
    DRLibrarySP  drLibrary; // points into field in wrapper

};

//// helper class for serialisation/deserialisation. Ensures we record what
//// the corresponding DR libray is
class XObject::Wrapper: public CObject,
               public virtual ITypeConvert,
               public virtual IRegressionTest{
    friend class XObject;
    // fields
    DRLibrarySP      drLibrary;
    XObjectConstSP   xObject;
public:
    ~Wrapper(){}

    //// note: just takes reference
    Wrapper(const DRLibrarySP& drLib, const XObjectConstSP& xObject):
        CObject(WRAPPER_TYPE), drLibrary(drLib), xObject(xObject){}

    /** IRegressionTest interface */
    IObjectSP runTest() const{
        return xObject->execute();
    }

    /** Converts this object to the enclosed XObject */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const{
        if (!requiredType->isInstance(*xObject)){
            throw ModelException("XObject::Wrapper::convert", "Require object "
                                 "of type "+requiredType->getName()+
                                 " but type "+
                                 xObject->getClass()->getName()+" supplied");
        }
        // the problems of const...
        object = XObjectSP::constCast(xObject);
    }

    /** Wraps the supplied writer into a new Writer that allows access
        to DRLibrary */
    virtual void write(const string& tag, Writer* writer) const{
        // create new writer containing DRLibrary
        DRLibrarySP drLib(DRLibrary::getLibrary(xObject->object->svc));
        XWriter xWriter(writer, drLib);
        // call CObject::write with new writer
        CObject::write(tag, &xWriter);
    }

    /** Wraps the supplied reader into a new reader that allows access
        to DRLibrary  */
    virtual void import(Reader::Node* elem, Reader* reader){
        static const string method("XObject::Wrapper::import");
        Reader::NodeListSP nl(elem->children());
        if (nl->empty()){
            throw ModelException(method, "XObject::Wrapper is empty");
        }
        // read drLibrary field in
        Reader::Node* child = nl->front().get();
        IObjectSP obj(reader->read(child));
        if (!DRLibrary::TYPE->isInstance(obj)){
            throw ModelException(method, "First component of XObject::Wrapper "
                                 "must be a DRLibrary");
        }
        DRLibrarySP drLib(DRLibrarySP::dynamicCast(obj));
        // create new XReader
        XReader xReader(reader, drLib);
        // then just call default with new reader
        CObject::import(elem, &xReader);
    }
private:
    Wrapper(): CObject(WRAPPER_TYPE){}
    static IObject* defaultConstructor(){
        return new Wrapper();
    }
    static void load(CClassSP& clazz){
        REGISTER(Wrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITypeConvert);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(drLibrary, "DR Library");
        FIELD(xObject, "External object from DR Library");
    }
};

CClassConstSP const XObject::WRAPPER_TYPE =
CClass::registerClassLoadMethod("XObject::Wrapper", typeid(Wrapper),
                                Wrapper::load);

MyDRObjectInterface::~MyDRObjectInterface(){
    (*svc->fptr->objectFree)(svc, (DRObject)this);
}
void MyDRObjectInterface::operator delete(void *ptr){}

XObject::~XObject(){

}

//// Simple constructor
XObject::XObject(const MyDRObjectInterfaceSP& object):
    CObject(TYPE), object(object.release()) {} // FIXME

/** Simple constructor. Takes ownership of memory */
XObject::XObject(DRObject obj): CObject(TYPE), object((MyDRObject)obj){}

/** Are these objects equal (ie contain the same data) */
bool XObject::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const XObject* theObj = STATIC_CAST(XObject, drObj.get());
    // check they are the same type first
    if (theObj->getXClass() != getXClass()){
        return false;
    }
    // then turn both objects into maps and compare those
    XMapSP xMap1(toMap());
    XMapSP xMap2(theObj->toMap());
    return (xMap1->equals(xMap2));
}

/** If we own the DRObject then we create another XObject pointing to
    the same DRObject. Otherwise copy the object - by turning into map
    and back again */
IObject* XObject::clone() const{
    if (/*object.*/owns()){
        return new XObject(object);
    } else {
        XMapSP xMap(toMap());
        XObjectSP xObj(xMap->toObject());
        return xObj.release();
    }
}

/** Create a clone of this object by recursively cloning every component
    of this object. This is used for regression testing and should not
    be used in general (unnecessarily slow) */
XObjectSP XObject::recursiveClone() const{
    try{
        XMapSP xMap(toMap());
        XObjectSP newMap(xMap->recursiveClone());
        XMap* newxMap = STATIC_CAST(XMap, newMap.get());
        XObjectSP xObj(newxMap->toObject());
        return xObj;
    } catch (exception& e){
        throw ModelException(e, "XObject::recursiveClone", "Failed to "
                             "recursively clone "+getXClassName());
    }
}
/** Override default CObject implementation.  */
void XObject::write(const string& tag, Writer* writer) const{
    // is it our writer?
    XWriter* xWriter = dynamic_cast<XWriter*>(writer);
    if (xWriter){
        // our writer, but is it the same DRLibrary?
        DRLibrarySP drLib(xWriter->getDRLibrary());
        // obvious check to make
        if (object->svc == drLib->getService() ||
            // more subtle check - there may be more than one svc
            // associated with a DRLibrary
            DRLibrary::getLibrary(object->svc)->equals(drLib)){
            // it is the same DRLibrary so just write out contents
            xWrite(tag, writer);
            return;
        }
    }
    // So either top level XObject or we really do have a different DR Library
    // so create Wrapper object and serialise that
    DRLibrarySP drLib(DRLibrary::getLibrary(object->svc));
    Wrapper wrapper(drLib, XObjectConstSP(this));
    wrapper.write(tag, writer);
}

/** Override default CObject implementation. This copes with XObject and
    XMap */
void XObject::xWrite(const string& tag, Writer* writer) const{
    try {
        // write external type as an attribute
        string attribute = XTYPE_ATTRIBUTE+"='" + getXClassName() + "'";
        IObjectConstSP obj(writer->objectStart(tag, attribute, this, true));
        if (obj.get()){ // if not already written out
            XMapConstSP objAsMap(toMap());
            // then write out contents of map
            objAsMap->xWriteElts(writer);
        }
        // all done. Just close object
        writer->objectEnd(tag, this);
    } catch (exception& e){
        throw ModelException(e, "XObject::xWrite");
    }
}

/** Override default CObject implementation */
void XObject::import(Reader::Node* elem, Reader* reader){
    static const string method("XObject::import");
    // make sure it's our Reader
    XReader* xReader = dynamic_cast<XReader*>(reader);
    if (!xReader){
        throw ModelException(method, "Invalid reader. XObject must be wrapped "
                             "in a XObject::Wrapper");
    }
    // retrieve the DRLibrary
    DRLibrarySP drLib(xReader->getDRLibrary());
    xImport(elem, reader, drLib->getService());
}


//// Derived classes should override this. Essentially X version of
//// import()
void XObject::xImport(Reader::Node* elem, Reader* reader, DRService* svc){
    static const string method("XObject::xImport");
    string type;
    try {
        // first get type of the external object
        type = elem->attribute(XTYPE_ATTRIBUTE);
        Reader::NodeListSP nl(elem->children());
        XMapSP xMap(XMap::create(type.c_str(), svc));
        xMap->xImportElts(elem, reader); // import the elements of the map
        // turn map into object
        XObjectSP xObj(xMap->toObject());
        // then copy pointer over
        this->object = xObj->object;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Override default CObject implementation */
void XObject::outputWrite(const string& linePrefix,
                          const string& prefix,
                          ostream&      stream) const{
    XMapConstSP objAsMap(toMap());
    objAsMap->outputWrite(linePrefix, prefix, stream);
}

/** ClientRunnable interface */
IObjectSP XObject::run(){
    return runTest();
}

/** IRegressionTest interface */
IObjectSP XObject::runTest() const{
    return execute();
}

/** CombinableResult interface: scale by factor x */
void XObject::scale(double x){
    // turn into map, scale elements of map, and then rebuild object.
    XMapSP map(toMap());
    map->scale(x);
    XObjectSP xObj(map->toObject());
    // then copy pointer over
    this->object = xObj->object;
}

/** CombinableResult interface:
    add an object (scaled by scaleFactor) to this
    result. Implementations should modify this result. If the x is
    not the same type as this then a [class cast] exception will
    be thrown */
void XObject::add(const CombinableResult& x, double scaleFactor){
    const XObject* objToAdd = DYNAMIC_CAST(XObject, &x);
    // turn into maps, invoke add on map, and then rebuild object.
    XMapSP map(toMap());
    XMapSP mapToAdd(objToAdd->toMap());
    map->add(*mapToAdd, scaleFactor);
    XObjectSP xObj(map->toObject());
    // then copy pointer over
    this->object = xObj->object;
}

/** Creates string from char* which is then DR freed */
string XObject::makeString(const char* str) const{
    DRService* svc = object->svc;
    return makeString(str, svc);
}

//// Creates a DRString from the supplied stl string
DRString XObject::makeDRString(const string& str, DRService* svc){
    const char* s = str.c_str();
    DRString  drString;
    CHECK((*svc->fptr->stringNew)(svc, s, &drString),
          "XObject:::makeDRString");
    return drString;
}

/** Creates string from char* which is then DR freed */
string XObject::makeString(const char* str, DRService* svc){
    string s(str);
    (*svc->fptr->stringFree)(svc, str);
    return s;
}

void XObject::throwError(const char* method, const char *errorMsg) {
    ModelException e(errorMsg);
    e.addMsg(method);
    throw e;
}

/** Returns the type of the external object */
string XObject::getXClassName() const{
    DRString         typeName;
    DRService* svc = object->svc;
    CHECK((*svc->fptr->objectGetType)(svc, (DRObject)object.get(), &typeName),
          "getXClassName");
    return makeString(typeName);
}

/** Returns the type of the external object */
XClassConstSP XObject::getXClass() const{
    DRString         typeName;
    DRService* svc = object->svc;
    CHECK((*svc->fptr->objectGetType)(svc, (DRObject)object.get(), &typeName),
          "getXClass");
    const string& name = makeString(typeName);
    // look up DRLibrary first - this simplifies once we can get service name
    // of the svc (or alternatively guarantee that service ptrs are unique
    // per service
    DRLibrarySP lib(DRLibrary::getLibrary(svc));
    return XClass::forName(lib->getServiceName(), name);
}

/** Turns DR Object into DR Map */
XMapSP XObject::toMap() const{
    DRService* svc = object->svc;
    DRMap      map;
    CHECK((*svc->fptr->objectToMap)(svc, (DRObject)object.get(), &map),
          "XObject::toMap");
    return XMapSP(new XMap(map));
}

/** Turns object into map and retrieves relevant field */
IDRObjectSP XObject::getField(const string& fieldForName) const{
    return getField(fieldForName.c_str());
}

/** Turns object into map and retrieves relevant field */
IDRObjectSP XObject::getField(const char* fieldForName) const{
    XMapSP myMap(toMap());
    return myMap->getItem(fieldForName);
}

/** Turn DRValue into IDRObject. If takeOwnership is true then you
    must not free the DRValue - this is done if needed by this
    function. Alternatively if takeOwnership is false then it is up
    to the caller to free the DRValue (if necessary) */
IDRObjectSP XObject::toIDRObject(const DRValue& value,
                                 bool           takeOwnership,
                                 DRService*     svc){
    static const string method("XObject::toIDRObject");
    IDRObject* obj;
    switch (value.type){
    case DR_UNDEFINED:
        throw ModelException(method, "Undefined type");
        break;
    case DR_BOOL:
        obj = CBool::create(value.value.boolean? true: false);
        break;
    case DR_INT:
        obj = CInt::create(value.value.integer);
        break;
    case DR_DOUBLE:
        obj = CDouble::create(value.value.real);
        break;
    case DR_STRING:
        obj = CString::create(value.value.string);
        if (takeOwnership){
            // now free the string
            (*svc->fptr->stringFree)(svc, value.value.string);
        }
        break;
    case DR_OBJECT:
        return toXObject(value.value.object, takeOwnership);
        break;
    case DR_DATE:
        obj = GDRDate::create(value.value.date);
        break;
    default:
        throw ModelException(method, "Illegal DR type");
    }
    return IDRObjectSP(obj);
}

/** Same as above but for objects */
XObjectSP XObject::toXObject(DRObject       drObj,
                             bool           takeOwnership){
    XObject* obj;
    if (!drObj){
        obj = 0;
    } else {
        DRService* svc = drObj->svc;
        MyDRObject myDrObj = (MyDRObject) drObj;
        MyDRObjectInterfaceSP drObjSP;
        drObjSP = takeOwnership? MyDRObjectInterfaceSP(myDrObj):
                MyDRObjectInterfaceSP(myDrObj /*, NullDeleter()*/); // FIXME: NullDestructor
        if (isArray(drObj)){
            obj = new XArray(drObjSP);
        } else if (isMatrix(drObj, svc)){
            obj = new XMatrix(drObjSP);
        } else {
            obj = new XObject(drObjSP);
        } /*  there doesn't seem to be any provision for passing maps
              through this bit */
        obj->weOwnObject = takeOwnership;

    }
    return XObjectSP(obj);
}


/** Turn IObject into DRValue - only free the DRValue for strings (in which case
    the return value is true (this is a right mess) */
bool XObject::fromIDRObject(
    const IDRObjectSP& theObj, DRService* svc, DRValue& value){
    IDRObject* object = theObj.get();
    if (!object){
        value.type = DR_OBJECT;
        value.value.object = 0;
    } else {
        CClassConstSP clazz = object->getClass();
        if (TYPE->isAssignableFrom(clazz)){
            value.type = DR_OBJECT;
            XObject* xobj = STATIC_CAST(XObject, object);
            value.value.object = (DRObject)xobj->object.get();
        } else if (clazz == CDouble::TYPE){
            value.type = DR_DOUBLE;
            CDouble* obj = STATIC_CAST(CDouble, object);
            value.value.real = obj->doubleValue();
        } else if (clazz == CInt::TYPE){
            value.type = DR_INT;
            CInt* obj = STATIC_CAST(CInt, object);
            value.value.integer = obj->intValue();
        } else if (clazz == CBool::TYPE){
            value.type = DR_BOOL;
            CBool* obj = STATIC_CAST(CBool, object);
            value.value.boolean = obj->boolValue();
        } else if (clazz == CString::TYPE){
            value.type = DR_STRING;
            CString* obj = STATIC_CAST(CString, object);
            // note have to create string through DR Interface - which
            // seems like madness
            value.value.string = makeDRString(obj->stringValue(), svc);
            return true;
        } else if (clazz == GDRDate::TYPE){
            value.type = DR_DATE;
            GDRDate* dt = STATIC_CAST(GDRDate, object);
            value.value.date = dt->dateValue();
        } else {
            throw ModelException("XObject::fromIDRObject",
                                 "Can't turn internal "
                                 "type "+ theObj->getClass()->getName()+" into "
                                 "DRObject");
        }
    }
    return false;
}

/** Try to convert the supplied IObject to an IDRObject. If the
    supplied object is an IDRObject already then this function
    just returns the casted object. If the object is an array of
    IDRObjects (eg array of Strings) then an XArray will be created of
    the type of the array. The service used to create the XArray is
    that in this object (ie only getService() property of this is used) */
IDRObjectSP XObject::toIDRObject(IObjectSP obj) const{
    static const string method("XObject::toIDRObject");
    if (IDRObject::TYPE->isInstance(obj)){
        return IDRObjectSP::dynamicCast(obj);
    }
    // is the object an array
    if (!IArray::TYPE->isInstance(obj)){
        throw ModelException(method, "Cannot turn an object"
                             " of type "+obj->getClass()->getName()+" into "
                             "an IDRObject");
    }
    // get the type of the components
    CClassConstSP cmptType = obj->getClass()->getComponentType();
    if (!IDRObject::TYPE->isAssignableFrom(cmptType)){
        throw ModelException(method, "Cannot convert an array of elements of "
                             "type "+cmptType->getName()+" into an IDRObject");
    }
    const char* xArrayEltType;
    // now switch on the different types possible
    if (cmptType == CString::TYPE){
        xArrayEltType = DRI_TYPE_STRING;
    } else if (cmptType == CDouble::TYPE){
        xArrayEltType = DRI_TYPE_DOUBLE;
    } else if (cmptType == CInt::TYPE){
        xArrayEltType = DRI_TYPE_INT;
    } else if (cmptType == CBool::TYPE){
        xArrayEltType = DRI_TYPE_BOOL;
    } else if (cmptType == GDRDate::TYPE){
        xArrayEltType = DRI_TYPE_DATE;
    } else {
        throw ModelException(method, "Internal error - array of type "+
                             cmptType->getName()+" not catered for");
    }
    IArraySP origArray(IArraySP::dynamicCast(obj));
    int size = origArray->getLength();
    XArraySP theArray(XArray::create(xArrayEltType, size, object->svc));
    for (int i = 0; i < size; i++){
        theArray->set(i, origArray->get(i));
    }
    return theArray;
}

//// const version
IDRObjectConstSP XObject::toIDRObject(IObjectConstSP obj) const{
    return toIDRObject(IObjectSP::constCast(obj));
}

/** Is the DRObject an array */
bool XObject::isArray(DRObject drObj){
    DRService* svc = drObj->svc;
    DRBool  isArray;
    CHECK((*svc->fptr->objectIsArray)(svc, drObj, &isArray),
          "XObject::isArray");
    return isArray? true: false;
}

/** Is the DRObject a matrix */
bool XObject::isMatrix(DRObject drObj, DRService* svc){
    DRBool  isMatrix;
    CHECK((*svc->fptr->objectIsMatrix)(svc, drObj, &isMatrix),
          "XObject::isMatrix");
    return isMatrix? true: false;
}

/** returns the enclosed DRObject which then must be freed. */
DRObject XObject::getDRObject(){
    XMapSP xMap(toMap());
    XObjectSP xObj(xMap->toObject());
    DRObject drObj = (DRObject) xObj->object.release();
    return drObj;
}

//// Returns the DRService from the object
DRService* XObject::getService() const{
    return object->svc;
}

/** Returns the all the types in the specified service */
XArraySP XObject::typeList(DRService* svc) {
    DRArray allTypes;
    CHECK((*svc->fptr->typeList)(svc, &allTypes), "typeList");
    return XArraySP(new XArray(allTypes));
}

/** Invoke 'execute' on this object */
IDRObjectSP XObject::execute() const{
    DRService* svc = object->svc;
    DRValue   item;
    CHECK((*svc->fptr->execute)(svc, (DRObject)object.get(), &item),
          "XObject::execute");
    return toIDRObject(item, true /* we own the object*/, svc);
}

XObject::XObject(const CClassConstSP& clazz): CObject(clazz){}
XObjectSP XObject::create(DRObject obj){
    return XObjectSP(new XObject(obj));
}

XObject::XObject(const CClassConstSP& clazz,
                 DRObject             object):
    CObject(clazz), object((MyDRObject)object){}

XObject::XObject(const CClassConstSP& clazz,
                 const MyDRObjectInterfaceSP& object):
    CObject(clazz), object(object){}

void XObject::load(CClassSP& clazz){
    REGISTER(XObject, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IDRObject);
    IMPLEMENTS(ClientRunnable);
    IMPLEMENTS(CombinableResult);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

XObject::XObject(): CObject(TYPE) {}

IObject* XObject::defaultConstructor(){
    return new XObject();
}

CClassConstSP const XObject::TYPE =
CClass::registerClassLoadMethod("XObject", typeid(XObject), load);

class XObject::Create: public CObject{
    static CClassConstSP const TYPE;

    DRLibrarySP         library;
    Handle::RawStringSP externalHandle;

    static IObjectSP fromHandle(Create* params){
        return params->library->handleToObject(params->
                                               externalHandle->getString());
    }

    Create(): CObject(TYPE){}

    static IObject* defaultConstructor(){
        return new Create();
    }
    static void load(CClassSP& clazz){
        REGISTER(Create, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(library, "The DR Library");
        FIELD(externalHandle, "A handle from the external DR library");
        Addin::registerClassObjectMethod("XOBJECT",
                                         Addin::UTILITIES,
                                         "Reads an object using a handle"
                                         " from an external DR library",
                                         TYPE,
                                         true, // handle name
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)fromHandle);
    }
};

CClassConstSP const XObject::Create::TYPE =
CClass::registerClassLoadMethod("XObject::Create", typeid(Create), load);

/** Addin for building external object from it component parts */
class XObject::New: public CObject{
    static CClassConstSP const TYPE;

    DRLibrarySP              library;
    string                   typeName;
    StringArraySP            fieldNames;
    Handle::RawStringArraySP externalHandles;

    static IObjectSP buildXObject(New* params){
        static const string method("XObject::New::buildXObject");
        // get service ptr
        DRService* svc = params->library->getService();
        try{
            // build up map by first creating empty map
            XMapSP xMap(XMap::create(params->typeName.c_str(), svc));
            if (params->fieldNames.get()){
                // for ease
                CStringArray& fieldNames = *params->fieldNames;
                Handle::RawStringArray&  components = *params->externalHandles;
                // then for each field
                for (int i = 0; i < fieldNames.size(); i++){
                    if (!fieldNames[i].empty()){
                        if (!params->externalHandles || components.size() <= i){
                            throw ModelException(method, "Object for field "+
                                                 fieldNames[i]+" is missing");
                        }
                        if (components[i].get()){
                            const string& handle = components[i]->getString();
                            // skip over empty handles
                            if (!handle.empty()){
                                IDRObjectSP obj;
                                if (Handle::exists(handle)){
                                    // it's one of ours!
                                    IObjectSP theObj(
                                        Handle::fetch(handle, IDRObject::TYPE));
                                    obj = IDRObjectSP::dynamicCast(theObj);
                                } else {
                                    // get object for
                                    // handle (note no need to clone)
                                    obj = params->library->
                                        handleToObject(handle);
                                }
                                xMap->setItem(fieldNames[i], obj);
                            }
                        }
                    }
                }
            }
            return xMap->toObject();
        } catch (exception& e){
            throw ModelException(e, method, "Failed trying to build external"
                                 " object of type "+params->typeName);
        }
    }

    New(): CObject(TYPE){}

    static IObject* defaultConstructor(){
        return new New();
    }
    static void load(CClassSP& clazz){
        REGISTER(New, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(library, "The DR Library");
        FIELD(typeName, "Type of the external object to build");
        FIELD(fieldNames, "Names of the components");
        FIELD_MAKE_OPTIONAL(fieldNames);
        FIELD(externalHandles,"The components themselves");
        FIELD_MAKE_OPTIONAL(externalHandles);
        Addin::registerClassObjectMethod("XOBJECT_CREATE",
                                         Addin::UTILITIES,
                                         "Creates an external DR object using "
                                         "handles from external DR library",
                                         TYPE,
                                         true, // handle name
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)buildXObject);
    }
};

CClassConstSP const XObject::New::TYPE =
CClass::registerClassLoadMethod("XObject::New", typeid(New), load);

DRLIB_END_NAMESPACE
