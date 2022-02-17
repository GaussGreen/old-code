//----------------------------------------------------------------------------
//
//   Group       : Global Derivatives Research
//
//   Filename    : XMap.cpp
//
//   Description : Wrapper for external [DR interface] maps
//
//   Author      : Mark A Robson
//
//   Date        : 14 Nov 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XMap.hpp"
#include "edginc/DRAnalyticsInterface.h"
#include "edginc/Null.hpp"
#include "edginc/Writer.hpp"
#include "edginc/XArray.hpp"
#include "edginc/DRUtil.hpp"
#include "edginc/Hashtable.hpp"
#include ext_hash_set

DRLIB_BEGIN_NAMESPACE

// DRI_CHECK is defined in DRUtil.hpp. This simplifies the invocation of the
// DRI_CHECK macro by fixing some parameters.
#define CHECK(exec, method) \
    DRI_CHECK((exec), (method), 0, driErrorHandler, 0, ;, ;, ;)

XMap::~XMap(){}

XMap::XMap(): XObject(TYPE){}

/** Simple constructor. Takes ownership of memory */
XMap::XMap(DRMap map): XObject(TYPE, map){}

/** Simple constructor. */
XMap::XMap(const MyDRObjectInterfaceSP& obj): XObject(TYPE, obj){}

/** Are these objects equal (ie contain the same data) */
bool XMap::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const XMap* theObj = STATIC_CAST(XMap, drObj.get());
    // check they are the same type first
    if (theObj->getXObjectClass() != getXObjectClass()){
        return false;
    }
    // the tricky bit ...
    /* iterate through one map and get a set of all the fields */
    hash_set<string, Hashtable::StringHash> fields;
    for (Iterator iter(createIterator()); iter.nextElement(); ){
        fields.insert(iter.key());
    }
    // then iterator through other map
    unsigned int count = 0;
    for (Iterator iter2(theObj->createIterator()); iter2.nextElement();
         count++){
        const string& key = iter2.key();
        if (fields.find(iter2.key()) == fields.end()){
            return false;
        }
        IDRObjectSP value(getItem(key));
        IDRObjectSP value2(iter2.value());
        if (!value->equals(value2)){
            return false;
        }
    }
    if (count != fields.size()){
        return false;
    }
    return true;
}
/** If we own the DRObject then nothing is done (NB May need to
    review this since an XMap can be modified). Otherwise copy
    the object - by converting to object and back again (which is
    a bit weak but there is no clone method) */
IObject* XMap::clone() const{
    if (/*object.*/ owns()){
        return const_cast<IObject*>((const IObject*)this);
    }
    XObjectSP xObj(toObject());
    XObjectSP xMap(xObj->toMap());
    return xMap.release();
}

/** Create a clone of this object by recursively cloning every component
    of this object. This is used for regression testing and should not
    be used in general (unnecessarily slow) */
XObjectSP XMap::recursiveClone() const{
    try{
        // create new empty map
        XMapSP newMap(create(getXObjectClass().c_str(), getService()));
        // iterate over elements of map
        for (Iterator iter = createIterator(); iter.nextElement(); ){
            const string& id = iter.key();
            IDRObjectSP obj(iter.value());
            if (XObject::TYPE->isInstance(obj)){
                // for each XObject insert a recursive clone of itself
                XObject* xObj = STATIC_CAST(XObject, obj.get());
                obj = xObj->recursiveClone();
            }
            newMap->setItem(id, obj);
        }
        return newMap;
    } catch (exception& e){
        throw ModelException(e, "XMap::recursiveClone", "Failed to "
                             "recursively clone map of "+getXObjectClass());
    }
}

/** Returns the type of the external object that this map corresponds to */
string XMap::getXObjectClass() const{
    DRString         typeName;
    DRService* svc = object->svc;
    CHECK((*svc->fptr->mapGetTypeName)(svc, (DRObject)object.get(), &typeName),
          "getXObjectClass");
    return makeString(typeName);
}

/** Override default CObject implementation. */
void XMap::xWrite(const string& tag, Writer* writer) const{
    try {
        // write external type as an attribute
        string attribute = XTYPE_ATTRIBUTE+"='" + getXClassName() + "'";
        IObjectConstSP obj(writer->objectStart(tag, attribute, this, true));
        if (obj.get()){ // if not already written out
            xWriteElts(writer);
        }
        // all done. Just close object
        writer->objectEnd(tag, this);
    } catch (exception& e){
        throw ModelException(e, "XMap::write");
    }
}

//// write elements of map to writer
void XMap::xWriteElts(Writer* writer) const{
    for (Iterator iter = createIterator(); iter.nextElement(); ){
        const string& id = iter.key();
        IDRObjectSP obj(iter.value());
        if (!obj){
            writer->writeNull(id);
        } else {
            obj->write(id, writer);
        }
    }
}

/** Override default XObject implementation */
void XMap::xImport(Reader::Node* elem, Reader* reader, DRService* svc){
    static const string method("XMap::xImport");
    string type;
    try {
        // first get type of the external object
        type = elem->attribute(XTYPE_ATTRIBUTE);
        // create empty map
        XMapSP xMap(XMap::create(type.c_str(), svc));
        xMap->xImportElts(elem, reader); // import the elements of the map
        //  copy pointer over
        this->object = xMap->object;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/**  read elements of map from reader */
void XMap::xImportElts(Reader::Node* elem, Reader* reader){
    Reader::NodeListSP nl(elem->children());
    int numElts = nl->size();
    for (int i = 0; i < numElts; i++){
        Reader::Node* child = (*nl)[i].get();
        if (!elem->isNull()) {
            // read in the objects
            IObjectSP obj(reader->read(child));
            if (CNull::TYPE->isInstance(obj)){
                // do nothing
            } else if (!IDRObject::TYPE->isInstance(obj)){
                throw ModelException("XMap::xImportElts",
                                     "Encounterd object of type "+
                                     obj->getClass()->getName()+
                                     "\nNot a DRInterface type!");
            } else {
                IDRObjectSP drObj(IDRObjectSP::dynamicCast(obj));
                setItem(child->name(), drObj);
            }
        }
    }
}

/** Override default CObject implementation */
void XMap::outputWrite(const string& linePrefix,
                       const string& prefix,
                       ostream&      stream) const{
    int totalPrinted = 0;
    for (Iterator iter = createIterator(); iter.nextElement(); ){
        IDRObjectSP elt(iter.value());
        if (elt.get()){
            const string& fieldName = iter.key();
            elt->outputWrite(linePrefix, prefix == ""? fieldName:
                             prefix+"_"+fieldName,
                             stream);
            totalPrinted++;
        }
    }
    if (totalPrinted == 0){
        // print out type of object instead
        stream << linePrefix << prefix << "_XTYPE: " <<
            getXObjectClass() << endl;
    }
}

/** CombinableResult interface: scale by factor x */
void XMap::scale(double x){
    // doubtful whether implementations will support altering a map whilst
    // iterating over it. So create a new empty one
    XMapSP newMap(XMap::create(getXObjectClass().c_str(), object->svc));
    for (Iterator iter(createIterator()); iter.nextElement(); ){
        const string& key = iter.key();
        IDRObjectSP elt(iter.value());
        if (CombinableResult::TYPE->isInstance(elt)){
            CombinableResult& eltCR = dynamic_cast<CombinableResult&>(*elt);
            eltCR.scale(x);
        }
        newMap->setItem(key, elt);
    }
    // then copy pointer over
    this->object = newMap->object;
}

/** CombinableResult interface:
    add an object (scaled by scaleFactor) to this
    result. Implementations should modify this result. If the x is
    not the same type as this then a [class cast] exception will
    be thrown */
void XMap::add(const CombinableResult& x, double scaleFactor){
    const XMap* mapToAdd = DYNAMIC_CAST(XMap, &x);
    // doubtful whether implementations will support altering a map whilst
    // iterating over it. So create a new empty one
    XMapSP newMap(XMap::create(getXObjectClass().c_str(), object->svc));
    for (Iterator iter(createIterator()); iter.nextElement(); ){
        const string& key = iter.key();
        IDRObjectSP elt(iter.value());
        if (CombinableResult::TYPE->isInstance(elt)){
            CombinableResult& eltCR = dynamic_cast<CombinableResult&>(*elt);
            IDRObjectSP eltToAdd(mapToAdd->getItem(key));
            if (!CombinableResult::TYPE->isInstance(eltToAdd)){
                throw ModelException("XMap::add", "Couldn't add field "+
                                     key+" (different types)");
            }
            CombinableResult& eltToAddCR =
                dynamic_cast<CombinableResult&>(*eltToAdd);
            eltCR.add(eltToAddCR, scaleFactor);
        }
        newMap->setItem(key, elt);
    }
    // then copy pointer over
    this->object = newMap->object;
}

/** Creates a map of the specified type using the supplied service */
XMapSP XMap::create(const char* typeName, DRService* svc){
    DRMap   obj;
    CHECK((*svc->fptr->mapNew)(svc, typeName, &obj), "XMap::create");
    return XMapSP(new XMap(obj));
}

/** Turns DR Map into DR Object */
XObjectSP XMap::toObject() const{
    DRService* svc = object->svc;
    DRObject   obj;
    CHECK((*svc->fptr->mapToObject)(svc, (DRObject)object.get(), &obj),
          "XMap::toObject");
    return XObject::create(obj);
}

/** Turns object into map and retrieves relevant field */
IDRObjectSP XMap::getItem(const string& fieldName) const{
    return getItem(fieldName.c_str());
}

/** Turns object into map and retrieves relevant field */
IDRObjectSP XMap::getItem(const char* fieldName) const{
    DRValue value;
    DRService* svc = object->svc;
    CHECK((*svc->fptr->mapGetItem)(svc, (DRObject)object.get(),
                                   fieldName, &value),
          "XMap::getItem");
    return toIDRObject(value, true /* we own the object*/, svc);
}

void XMap::getItem(const char* fieldName, DRValue& drValue, int reqdType) const{
    static const char* drTypes[] = {"Undefined","bool", "int", "double",
                                    "string", "object", "date"};
    DRService* svc = object->svc;
    CHECK((*svc->fptr->mapGetItem)(svc, (DRObject)object.get(),
                                   fieldName, &drValue),
          "XMap::getItem");
    if (drValue.type != reqdType){
        throw ModelException("XMap::getItem", "Field "+string(fieldName)+
                             " is not a "+string(drTypes[reqdType])+
                             " but a "+string(drTypes[drValue.type]));
    }
}

bool XMap::getBool(const string& fieldName) const{
    DRValue drValue;
    getItem(fieldName.c_str(), drValue, DR_BOOL);
    return drValue.value.boolean == 0? false: true;
}

int XMap::getInt(const string& fieldName) const{
    DRValue drValue;
    getItem(fieldName.c_str(), drValue, DR_INT);
    return drValue.value.integer;
}
double XMap::getDouble(const string& fieldName) const{
    DRValue drValue;
    getItem(fieldName.c_str(), drValue, DR_DOUBLE);
    return drValue.value.real;
}
string XMap::getString(const string& fieldName) const{
    DRValue drValue;
    getItem(fieldName.c_str(), drValue, DR_STRING);
    string s(drValue.value.string);
    // now free the string
    (*object->svc->fptr->stringFree)(object->svc, drValue.value.string);
    return s;
}

XObjectSP XMap::getXObject(const string& fieldName) const{
    DRValue drValue;
    getItem(fieldName.c_str(), drValue, DR_OBJECT);
    return toXObject(drValue.value.object, true /* we own the object*/);
}

XArraySP XMap::getXArray(const string& fieldName) const{
    DRValue drValue;
    getItem(fieldName.c_str(), drValue, DR_OBJECT);
    if (!XObject::isArray(drValue.value.object)){
        throw ModelException("XMap::getXArray", "Field "+ fieldName +" is not "
                             "an array");
    }
    return XArraySP(new XArray(drValue.value.object));
}

/** Sets the specified item from the external map */
void XMap::setItem(const string& fieldName, const IDRObjectSP& obj) const{
    setItem(fieldName.c_str(), obj);
}

/** Sets the specified item from the external map */
void XMap::setItem(const char* fieldName, const IDRObjectSP& obj) const{
    DRService* svc = object->svc;
    DRValue value;
    bool freeStr = fromIDRObject(obj, svc, value);

    DRError errMsg = (*svc->fptr->mapAddItem)(svc, (DRObject)object.get(),
                                              fieldName, &value);
    if (freeStr){
        (*svc->fptr->stringFree)(svc, value.value.string);
    }
    CHECK(errMsg, "XMap::setItem");
}

class XMap::Iterator::Imp{
public:
    Imp(const XMap* xMap){
        DRService* svc = xMap->object->svc;
        CHECK((*svc->fptr->mapIteratorGet)(svc, (DRObject)xMap->object.get(),
                                           &drIter),
              "XMap::Iterator::Imp");
    }
    ~Imp(){
        DRService* svc = drIter->svc;
        (*svc->fptr->objectFree)(svc, drIter);
    }
    IDRObjectSP    value;
    string         key;
    DRMapIterator  drIter;
};

//// move to first/next element. Returns false if none available
bool XMap::Iterator::nextElement(){
    DRService* svc = my->drIter->svc;
    DRString  fieldName;
    DRValue   value;
    CHECK((*svc->fptr->mapIteratorNext)(svc, my->drIter, &fieldName, &value),
          "XMap::Iterator::nextElement");
    if (!fieldName){
        return false;
    }
    my->key = makeString(fieldName, svc);
    my->value = toIDRObject(value, true /* we own the object*/, svc);
    return true;
}
XMap::Iterator::~Iterator(){
    delete my;
}

const string& XMap::Iterator::key() const{
    return my->key;
}

IDRObjectSP XMap::Iterator::value() const{
    return my->value;
}

XMap::Iterator::Iterator(const XMap* xMap): my(new Imp(xMap)){}

XMap::Iterator XMap::createIterator() const{
    return Iterator(this);
}

void XMap::load(CClassSP& clazz){
    REGISTER(XMap, clazz);
    SUPERCLASS(XObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

IObject* XMap::defaultConstructor(){
    return new XMap();
}

CClassConstSP const XMap::TYPE =
CClass::registerClassLoadMethod("XMap", typeid(XMap), load);



DRLIB_END_NAMESPACE
