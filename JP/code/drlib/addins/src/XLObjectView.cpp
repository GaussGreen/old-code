//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLObjectView.cpp
//
//   Description : Object and data dictionary view addins
//
//   Author      : Stephen Hope
//
//   Date        : 3 Sep 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Null.hpp"
#include "edginc/Handle.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/XArray.hpp"
#include "edginc/XMap.hpp"

DRLIB_BEGIN_NAMESPACE


class ObjectViewer: public CObject{
public:
    static CClassConstSP const TYPE;
    // addin functions just take a single parameter and an optional name
    IObjectSP  object;
    string     fieldName;

    static IObjectSP viewer(ObjectViewer* inputs){
        static const string method = "ObjectViewer::viewer";
        IObjectSP result(0);
        // check if the object implements the IReadableMap interface
        if (IMap::TYPE->isInstance(inputs->object)) {
            IReadableMap& hashLikeObject = 
                dynamic_cast<IReadableMap&>(*inputs->object);
            if (inputs->fieldName.empty()) {
                // return all of the object elements in an array
                ObjectArraySP everything(new ObjectArray(0));

                IMap::IIteratorSP iter(hashLikeObject.createIterator());
                while (iter->hasMoreElements()) {
                    everything->push_back(iter->getElement());
                    iter->increment();
                }

                result = everything;
            } else {
                // return the single element requested
                result = hashLikeObject.get(inputs->fieldName);
                if (!result){
                    result = CNull::create();
                }
            }
        } else if (XObject::TYPE->isInstance(inputs->object)){
            ObjectArraySP elts;
            if (XArray::TYPE->isInstance(inputs->object)){
                XArraySP xArray(XArraySP::dynamicCast(inputs->object));
                elts = ObjectArraySP(new ObjectArray(xArray->getLength()));
                for (int i = 0; i < elts->size(); i++){
                    (*elts)[i] = xArray->getElt(i);
                }
            } else {
                XObjectSP xObj(XObjectSP::dynamicCast(inputs->object));
                XMapSP xMap(xObj->toMap());
                elts = ObjectArraySP(new ObjectArray());
                for (XMap::Iterator iter(xMap->createIterator());
                     iter.nextElement(); /* no op */){
                    elts->push_back(iter.value());
                }
            }
            result =  elts;
        } else {
            if (inputs->fieldName.empty()) // view the whole object
            {
                result = inputs->object;  
            } else {
                // view a particular field within the object
                CFieldConstSP field(0);
                CClassConstSP objClass = inputs->object->getClass();
                // recurse through parent classes if any
                do {
                    field = objClass->hasDeclaredField(inputs->fieldName);
                    if (field) {
                        result = field->get(inputs->object);
                    }
                }while ((objClass = objClass->getSuperClass()) && !field);
                if (!field) {
                    throw ModelException(method,
                                     inputs->fieldName +
                                     " not found in " +
                                     inputs->object->getClass()->getName() +
                                     " or it's parent classes!");
                }
            }
        }
        // ensure we don't return a private representation of an object
        // (In theory the incoming object can't be an internal one either)
        result = CObject::convertToPublicRep(result);
        // be careful we don't return something that might go out of scope
        if (result.get() && result->getRefCount() == 0){
            result = IObjectSP(result->clone());
        }
        return result;
    }

    ObjectViewer():  CObject(TYPE){}
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObjectViewer, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObjectViewer);
        FIELD(object, "Object to view");
        FIELD(fieldName, "Field within Object to view");
        FIELD_MAKE_OPTIONAL(fieldName);
        Addin::registerInstanceObjectMethod("OBJECT_VIEWER",
                                            Addin::XL_TESTS,
                                            "Unpacks object so you can see "
                                            "its component parts",
                                            TYPE,
                                            false,
                                            Addin::expandSimple,
                                            (Addin::ObjMethod*)viewer);

        Addin::registerInstanceObjectMethod("OBJECT_VIEWER_WITH_PREFIX",
                                            Addin::XL_TESTS,
                                            "Unpacks object so you can see "
                                            "its component parts",
                                            TYPE,
                                            true,
                                            Addin::expandSimple,
                                            (Addin::ObjMethod*)viewer);

        Addin::registerInstanceObjectMethod("OBJ_VIEW",
                                            Addin::XL_TESTS,
                                            "Unpacks object so you can see "
                                            "its component parts",
                                            TYPE,
                                            false,
                                            Addin::expandMulti,
                                            (Addin::ObjMethod*)viewer);

        Addin::registerInstanceObjectMethod("OBJ_VIEW_WITH_PREFIX",
                                            Addin::XL_TESTS,
                                            "Unpacks object so you can see "
                                            "its component parts",
                                            TYPE,
                                            true,
                                            Addin::expandMulti,
                                            (Addin::ObjMethod*)viewer);
}

    static IObject* defaultObjectViewer(){
        return new ObjectViewer();
    }
};

CClassConstSP const ObjectViewer::TYPE = CClass::registerClassLoadMethod(
    "ObjectViewer", typeid(ObjectViewer), load);

class FieldGet: public CObject,
                public virtual Handle::IName{
public:
    static CClassConstSP const TYPE;

    /** Return a default handle name when this class is used for an
        addin function */
    virtual string defaultHandleName() const{
        // return className.fieldName
        string name(object->getClass()->getName()+"."+fieldName);
        return name;
    }
    
    // take a handle and a component name
    IObjectSP  object;
    string     fieldName;

    static IObjectSP fieldViewer(FieldGet* inputs){
        static const string method = "FieldGet::fieldViewer";
        IObjectSP subObject;
        // check if the object implements the IReadableMap interface
        if (IReadableMap::TYPE->isInstance(inputs->object)) {
            IReadableMap& hashLikeObject
                = dynamic_cast<IReadableMap&>(*inputs->object);
            subObject = hashLikeObject.get(inputs->fieldName);
            if (!subObject){
                subObject = CNull::create();
            }
        } else {
            CFieldConstSP field(0);
            CClassConstSP objClass = inputs->object->getClass();
            // recurse through parent classes if any
            do {
                field = objClass->hasDeclaredField(inputs->fieldName);
                if (field) {
                    subObject = field->get(inputs->object);
                }
            }while ((objClass = objClass->getSuperClass()) && !field);
            
            if (!field) {
                throw ModelException(method,
                                     inputs->fieldName +
                                     " not found in " +
                                     inputs->object->getClass()->getName() +
                                     " or it's parent classes!");
            }
        }
        // ensure we don't return a private representation of an object
        // (In theory the incoming object can't be an internal one either)
        subObject = CObject::convertToPublicRep(subObject);
        // be careful we don't return something that might go out of scope
        if (subObject.get() && subObject->getRefCount() == 0){
            subObject = IObjectSP(subObject->clone());
        }
        return subObject;
    }
    
    FieldGet():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FieldGet, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFieldGet);
        FIELD(object, "Object containing component");
        FIELD(fieldName, "Component to view");
        Addin::registerInstanceObjectMethod("FIELD_GET",
                                            Addin::XL_TESTS,
                                            "returns a handle to a "
                                            "component part of the object",
                                            TYPE,
                                            false,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)fieldViewer);

        Addin::registerInstanceObjectMethod("FIELD_GET_WITH_PREFIX",
                                            Addin::XL_TESTS,
                                            "returns a handle to a "
                                            "component part of the object",
                                            TYPE,
                                            true ,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)fieldViewer);

        
    }

    static IObject* defaultFieldGet(){
        return new FieldGet();
    }
};

CClassConstSP const FieldGet::TYPE = CClass::registerClassLoadMethod(
    "FieldGet", typeid(FieldGet), load);


class XLFile_Load: public CObject{
public:
    static CClassConstSP const TYPE;
    string filename;
    
    static IObjectSP xmlFileLoad(XLFile_Load *params) {
        static const string method = "xmlFileLoad";        
        
        XMLReader reader(params->filename, true);
        // get the root of the document
        Reader::NodeSP root(reader.root());
        string rootname = root->name();
        if (rootname == "MULTIPERM") {
            // if multiple tests, just output first one.
            Reader::NodeListSP nl(root->children());
            for (unsigned int i = 0; i < nl->size(); i++) {
                Reader::Node* child = (*nl)[i].get();
                try {
                    return IObjectSP(reader.read(child));
                }
                catch (exception& e) {
                    throw ModelException(e, method);
                }
            }
            // didn't find ELEMENT_NODE
            throw ModelException(method,
                                 "couldn't read in " + params->filename);
        }
        else {
            // just one test
            try {
                return reader.read(root.get());
            }
            catch (exception& e) {
                throw ModelException(e, method, 
                                     "couldn't read in " + params->filename);
            }
        }
    }
    
    XLFile_Load(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLFile_Load, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLFile_Load);
        FIELD(filename, "filename");
        Addin::registerInstanceObjectMethod("XML_FILE_LOAD",
                                            Addin::XL_TESTS,
                                            "Returns handle to object "
                                            "contained in xml file",
                                            TYPE,
                                            true,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)xmlFileLoad);
    }
    
    static IObject* defaultXLFile_Load(){
        return new XLFile_Load();
    }
};

CClassConstSP const XLFile_Load::TYPE = CClass::registerClassLoadMethod(
    "XLFile_Load", typeid(XLFile_Load), load);


/* non static variable to force linkage of this file. Outside of class to
   avoid necessity of header file */
bool XLObjectViewRegistered = true;


DRLIB_END_NAMESPACE
