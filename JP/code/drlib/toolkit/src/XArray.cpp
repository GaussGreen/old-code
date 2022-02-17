//----------------------------------------------------------------------------
//
//   Group       : Global Derivatives Research
//
//   Filename    : XArray.cpp
//
//   Description : Wrapper for external [DR interface] objects which are
//                 arrays
//
//   Author      : Mark A Robson
//
//   Date        : 2 Dec 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XArray.hpp"
#include "edginc/DRAnalyticsInterface.h"
#include "edginc/XClass.hpp"
#include "edginc/Format.hpp"
#include "edginc/Writer.hpp"
#include "edginc/Null.hpp"
#include "edginc/DRLibrary.hpp"
#include "edginc/DRUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// DRI_CHECK is defined in DRUtil.hpp. This simplifies the invocation of the
// DRI_CHECK macro by fixing some parameters.
#define CHECK(exec, method) \
    DRI_CHECK((exec), (method), 0, driErrorHandler, 0, ;, ;, ;)

XArray::~XArray(){}

XArray::XArray(): XObject(TYPE){}

/** Simple constructor. */
XArray::XArray(const MyDRObjectInterfaceSP& theArray):
    XObject(TYPE, theArray){}

/** Simple constructor. Takes ownership of memory */
XArray::XArray(DRArray theArray): XObject(TYPE, theArray){}

int XArray::getLength()const{
    return size();
}

/** Returns the value of the indexed component in the
    specified array object. */
IObjectConstSP XArray::get(int index) const{
    return getElt(index);
}

/** Non const version of above */
IObjectSP XArray::get(int index){
    return getElt(index);
}


/** Sets the value of the indexed component of the specified
    array object to the specified new value */
void XArray::set(int index, IObjectSP value){
    IDRObjectSP xValue(IDRObjectSP::dynamicCast(value));
    setElt(index, xValue);
}

/** Appends the supplied object to the end of supplied array. The
    length of the array is increased by 1. (aka push_back).
    Note this implementation is particularly slow due to lack of
    append method in DRInterface */
void XArray::append(IObjectSP value){
    int len = size();
    XArraySP newArray(create(elementType(), len+1, object->svc));
    for (int i = 0; i < len; i++){
        newArray->setElt(len, getElt(i));
    }
    IDRObjectSP xValue(IDRObjectSP::dynamicCast(value));
    newArray->setElt(len, xValue);
    // then copy pointer over
    this->object = newArray->object;
}

/** Are these objects equal (ie contain the same data) */
bool XArray::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const XArray* theObj = STATIC_CAST(XArray, drObj.get());
    // check they are the same type first
    if (theObj->getXClass() != getXClass()){
        return false;
    }
    int mySize = size();
    if (mySize != theObj->size()){
        return false;
    }
    for (int i = 0; i < mySize; i++){
        if (!getElt(i)->equals(theObj->getElt(i))){
            return false;
        }
    }
    return true;
}

/** If we own the DRObject then nothing is done (NB May need to review this
    since an XArray can be modified). Otherwise
    copy the object - by building empty array and getting/setting
    (which is a bit weak but there is no clone method) */
IObject* XArray::clone() const{
    if (owns()){
        return const_cast<IObject*>(/*(const IObject*)*/ static_cast<const IObject*>(this));
    }
    int length = size();
    XArraySP xArray(create(elementType().c_str(), length, object->svc));
    for (int i = 0; i < length; i++){
        xArray->setElt(i, getElt(i));
    }
    return xArray.release();
}

/** Create a clone of this object by recursively cloning every component
    of this object. This is used for regression testing and should not
    be used in general (unnecessarily slow) */
XObjectSP XArray::recursiveClone() const{
    try{
        int length = size();
        XArraySP xArray(create(elementType().c_str(), length, object->svc));
        for (int i = 0; i < length; i++){
            IDRObjectSP elt(getElt(i));
            if (XObject::TYPE->isInstance(elt)){
                // for each XObject insert a recursive clone of itself
                XObject* xObj = STATIC_CAST(XObject, elt.get());
                elt = xObj->recursiveClone();
            }
            xArray->setElt(i, elt);
        }
        return xArray;
    } catch (exception& e){
        throw ModelException(e, "XArray::recursiveClone", "Failed to "
                             "recursively clone "+getXClassName());
    }
}

/** Override default CObject implementation. */
void XArray::xWrite(const string& tag, Writer* writer) const{
    try {
        // write external type + length as an attribute
        int length = size();
        string attribute = XTYPE_ATTRIBUTE+"='"+ elementType() + "' " +
            CArray::ARRAY_LENGTH+"='" + Format::toString(length) + "'";
        // start the object
        IObjectConstSP obj(writer->objectStart(tag, attribute, this, true));
        if (obj.get()){ // if not already written out
            char buffer[20];
            for (int i = 0; i < length; i++){
                IDRObjectSP elem(getElt(i));
                sprintf(buffer, "%s%d", "Item", i);
                if (!elem){
                    // need to handle null elements in list
                    writer->writeNull(buffer);
                } else {
                    elem->write(buffer, writer);
                }
            }
        }
        writer->objectEnd(tag, this); // all done. Just close object
    } catch (exception& e){
        throw ModelException(e, "XArray::write");
    }
}
/** Override default XObject implementation */
void XArray::xImport(Reader::Node* elem, Reader* reader, DRService* svc){
    static const string method("XArray::import");
    string type;
    try {
        // first get type of the external object
        type = elem->attribute(XTYPE_ATTRIBUTE);
        /* a bit dubious as whether the elt is an array or not (as
           far as the reader is concerned) is a property of the
           type */
        int length = elem->arrayLength();
        Reader::NodeListSP nl(elem->children());
        int numElts = nl->size();
        if (numElts > length){
            throw ModelException(method, "Supplied array has too many "
                                 "elements");
        }
        // create empty array
        XArraySP xArray(create(type, length, svc));
        for (int i = 0; i < numElts; i++){
            Reader::Node* child = (*nl)[i].get();
            if (!elem->isNull()) {
                // read in the objects
                IObjectSP obj(reader->read(child));
                if (CNull::TYPE->isInstance(obj)){
                    // do nothing
                } else if (!IDRObject::TYPE->isInstance(obj)){
                    throw ModelException(method, "Encounterd object of type "+
                                         obj->getClass()->getName()+
                                         "\nNot a DRInterface type!");
                } else {
                    IDRObjectSP drObj(IDRObjectSP::dynamicCast(obj));
                    xArray->setElt(i, drObj);
                }
            }
        }
        // then copy pointer over
        this->object = xArray->object;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Override default CObject implementation */
void XArray::outputWrite(const string& linePrefix,
                         const string& prefix,
                         ostream&      stream) const{
    int length = size();
    char buffer[16];
    for (int i = 0; i < length; i++){
        IDRObjectSP elt(getElt(i));
        if (elt.get()){
            sprintf(buffer, "[%d]", i);
            elt->outputWrite(linePrefix, prefix+buffer, stream);
        }
    }
}

/** CombinableResult interface: scale by factor x */
void XArray::scale(double x){
    // just iterate over elements of the array
    int len = size();
    for (int i = 0; i < len; i++){
        IDRObjectSP elt(getElt(i));
        if (CombinableResult::TYPE->isInstance(elt)){
            CombinableResult& eltCR = dynamic_cast<CombinableResult&>(*elt);
            eltCR.scale(x);
            setElt(i, elt);
        }
    }
}

/** CombinableResult interface:
    add an object (scaled by scaleFactor) to this
    result. Implementations should modify this result. If the x is
    not the same type as this then a [class cast] exception will
    be thrown */
void XArray::add(const CombinableResult& x, double scaleFactor){
    const XArray* xToAdd = DYNAMIC_CAST(XArray, &x);
    // just iterate over elements of the array
    int len = size();
    for (int i = 0; i < len; i++){
        IDRObjectSP elt(getElt(i));
        if (CombinableResult::TYPE->isInstance(elt)){
            IDRObjectSP eltToAdd(xToAdd->getElt(i));
            if (!CombinableResult::TYPE->isInstance(eltToAdd)){
                string m("Element "+ Format::toString(i+1)+ " of the two arrays"
                         " cannot be added");
                throw ModelException("XArray::add", m);
            }
            CombinableResult& eltCR = dynamic_cast<CombinableResult&>(*elt);
            CombinableResult& eltToAddCR =
                dynamic_cast<CombinableResult&>(*eltToAdd);
            eltCR.add(eltToAddCR, scaleFactor);
            setElt(i, elt);
        }
    }
}


/** Not sure if this belongs here or in a new class XArray.
    Anyway, provided this object is an array this get's the n'th
    element. */
IDRObjectSP XArray::getElt(int index) const{
    DRService* svc = object->svc;
    DRValue   item;
    CHECK((*svc->fptr->arrayGet)(svc, (DRObject)object.get(), index, &item),
          "XArray::getElt");
    return toIDRObject(item, true /* we own the object*/, svc);
}

/** Not sure if this belongs here or in a new class XArray.
    Anyway, provided this object is an array this set's the n'th
    element. */
void XArray::setElt(int index, IDRObjectSP obj){
    DRService* svc = object->svc;
    DRValue   item;
    bool freeStr = fromIDRObject(obj, svc, item);
    DRError errMsg =
        (*svc->fptr->arraySet)(svc, (DRObject)object.get(), index, &item);
    if (freeStr){
        (*svc->fptr->stringFree)(svc, item.value.string);
    }
    CHECK(errMsg, "XArray::setElt");
}

/** returns its length */
int XArray::size()const{
    DRService* svc = object->svc;
    // get length
    int length;
    CHECK((*svc->fptr->arrayLength)(svc, (DRObject)object.get(), &length),
          "XArray::size");
    return length;
}

/** Returns the type of elements of this array (provided isArray() is
    true). An empty string is returned if the corresponding DR interface
    method returns null */
string XArray::elementType() const{
    DRService* svc = object->svc;
    // get arrayElementType
    DRString type;
    CHECK((*svc->fptr->arrayElementType)(svc, (DRObject)object.get(), &type),
          "XArray::elementType");
    return type? makeString(type): string();
}

/** Returns the XClass representing the elements of this array. Beware
    returns null if the corresponding DRI function returns null (which
    it's not supposed to but Enhanced Magnet still is) */
XClassConstSP XArray::elementClass() const{
    const string& eltType = elementType();
    if (eltType.empty()){
#if 0
        // when enhanced magnet works properly....
        throw ModelException("XArray::elementClass", "Failed to determine "
                             "type of components of array");
#else
        return 0;
#endif
    }
    // then get hold of the service
    DRService* svc = getService();
    // then the corresponding DRLibrary
    DRLibrarySP lib(DRLibrary::getLibrary(svc));
    // then finally the XClass for the component
    return XClass::forName(lib->getServiceName(), eltType);
}

/** Creates an array of specified length with elements of specified type.
    eltTypeName can be empty in which case null is used in the equivalent
    DR interace function */
XArraySP XArray::create(const string& eltTypeName,
                        int           numElements,
                        DRService*    svc){
    const char* arrayEltType = eltTypeName.empty()? 0: eltTypeName.c_str();
    // create Array
    DRArray theArray;
    CHECK((*svc->fptr->arrayNew)(svc, numElements, arrayEltType, &theArray),
          "XArray::createArray");
    return XArraySP(new XArray(MyDRObjectInterfaceSP(
                                   (MyDRObject)theArray)));
}

/** Creates an array from supplied vector of XObjects using specified type.
    The vector must not be empty or contain nulls */
XArraySP XArray::create(const string&            eltTypeName,
                        const vector<XObjectSP>& elts){
    static const string method("XArray::create");
    if (elts.empty() || !elts.front()){
        throw ModelException(method, "Zero/null elements not supported");
    }
    try{
        XArraySP xArr(create(eltTypeName, elts.size(),
                             elts.front()->getService()));
        for (unsigned int i = 0; i < elts.size(); i++){
            XObjectSP elt(elts[i]); // 2 steps: avoid error on gcc 2.95
            xArr->setElt(i, elt);
        }
        return xArr;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

void XArray::load(CClassSP& clazz){
    REGISTER(XArray, clazz);
    SUPERCLASS(XObject);
    IMPLEMENTS(IArray);
    //// we can't instantiate an array unless we have the relevant service ptr
    //// so this is a bit of a hack
    clazz->setArrayType(typeid(XObject), 0);
}

CClassConstSP const XArray::TYPE =
CClass::registerClassLoadMethod("XArray", typeid(XArray), load);



DRLIB_END_NAMESPACE
