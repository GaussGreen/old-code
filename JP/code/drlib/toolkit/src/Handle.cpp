//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Handle.hpp
//
//   Description : Class for holding handles - which are string references
//                 to objects
//
//   Author      : Mark A Robson
//
//   Date        : 21 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_HANDLE_CPP
#include "edginc/Addin.hpp"
#include "edginc/Handle.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/PrivateObject.hpp"


DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<Handle>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<Handle>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<Handle::RawString>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<Handle::RawStringSP _COMMA_ 
                     Handle::RawString>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<Handle::RawStringArray>);

/** Holds all handles in use. Handles are essentially string
    references to objects. They are primarily used so that pointers
    can be 'returned' to the spreadsheet - and provide for type safety */

HashtableSP Handle::handles(new Hashtable());

Handle::SpreadSheetCellNameMethod* Handle::sheetNameMethod = 0;
const int Handle::MAX_ID = 10000; // max handle
const int Handle::MAX_ID_LENGTH = 4; /* number of character needed
                                        for MAX_ID */
const int Handle::MIN_HANDLE_LENGTH = 6; /* 4 for handle, 1 for user name
                                            1 for '.' */
const int Handle::MAX_HANDLE_LENGTH = 256; // fairly arbitrary

/** Sets the method used for finding the current cell in the spreadsheet.
    Method can be null (default), in which case "NULL" is used */
void Handle::setSpreadSheetCellNameMethod(
    SpreadSheetCellNameMethod* sSheetMethod){
    sheetNameMethod = sSheetMethod;
}


Handle::Handle(const IObjectSP& object): 
    CObject(TYPE), id(1), object(object){
}

/** create and store a handle to an object. Returns complete handle name */
string Handle::create(const string&     userName,
                      const IObjectSP&  object){
    static const string routine = "Handle::create";
    if (userName.empty()){
        throw ModelException(routine, "Handle name is empty");
    }
    if (!isalpha(userName[0])) {
        throw ModelException(routine,"handle "+userName+
                             " must start with a character [a-zA-Z]");
    }
    if (!object){
        throw ModelException(routine,"Object for handle "+userName+" is null");
    }

    string hashName = userName + "." +
        (sheetNameMethod? sheetNameMethod(): "NULL");
    if ((int)hashName.size() + MIN_HANDLE_LENGTH > MAX_HANDLE_LENGTH){
        throw ModelException("Handle name is too long: "+hashName);
    }
    IObjectSP objectToStore;
    // convert to public rep if needed (eg when creating a handle to a 
    // component in an object or an array)
    if (IPrivateObject::TYPE->isInstance(object)){
        objectToStore = convertToPublicRep(object);
    } else {
        objectToStore = object;
    }
    HandleSP handle(HandleSP::dynamicCast(handles->get(hashName)));
    if (!handle){
        handle = HandleSP(new Handle(objectToStore));
        handles->put(hashName, handle);
    } else {
        // increment handle count and update object
        handle->id = (handle->id + 1) % MAX_ID;
        handle->object = objectToStore;
    }
    sprintf(handle->idAsString, "%4.4d", handle->id); 
    return (string(handle->idAsString)+hashName);
}

/** Returns true if given string has the right format to be a handle. */
bool Handle::valid(const string& name){
    if (name.empty() || (int)name.size() < MIN_HANDLE_LENGTH ||
        !isalpha(name[MAX_ID_LENGTH]) ||
        (int)name.size() > MAX_HANDLE_LENGTH){
        return false;
    }
    for (int idx = 0; idx < MAX_ID_LENGTH; idx++){
        if (!isdigit(name[idx])){
            return false;
        }
    }
    return true;
}

/** Returns the user part of a full handle name eg for 0001Divs.[Sheet1!A1]
    returns Divs */
string Handle::getUserHandleName(const string& fullHandle){
    int length = fullHandle.size();
    if (length < MAX_ID_LENGTH){
        throw ModelException("Handle::getUserHandleName", "Invalid handle");
    }
    int i;
    for (i = MAX_ID_LENGTH; i < length && fullHandle[i] != '.'; i++){
        // empty loop
    }
    return fullHandle.substr(MAX_ID_LENGTH, i-MAX_ID_LENGTH);
}

/** Returns true if given string has the right format to be a handle. Note
    here countedStringL is a 'counted string' - it is not null terminated 
    but rather the first byte gives the number of characters */
bool Handle::valid(const char* countedStringL){
    if (!countedStringL || *countedStringL < MIN_HANDLE_LENGTH  || 
        !isalpha(countedStringL[MAX_ID_LENGTH+1])){
        return false;
    }
    for (int idx = 0; idx < MAX_ID_LENGTH; idx++){
        if (!isdigit(countedStringL[idx+1])){ /* characters offset by 1 
                                               as countedString */
            return false;
        }
    }
    for (int i = MAX_ID_LENGTH; i < *countedStringL; i++){
        if (countedStringL[i+1] == '.'){
            // now check remaining string is something that Excel/VB gives us
            // Just check length for now...
            return (i < (*countedStringL) - 1);
        }
    }
    return false;
}

/** Returns true if the handle with the given name exists */
bool Handle::exists(const string& name){
    if (valid(name)){
        // need to extract bit after id
        string hashName = name.substr(MAX_ID_LENGTH);
        return ((!handles->get(hashName))? false: true);
    }
    return false;
}
        
IObjectSP Handle::fetch(const string& name,
                        CClassConstSP objectClass){ // can be null
    static const string routine = "Handle::fetch";
    IObjectSP           obj; // return value
    try{
        // check that we've got a handle
        if (!valid(name)){
            throw ModelException(routine, "Invalid handle name: "+name+
                                 "\nHandle names have the format: "
                                 "0001<user handle name>.<cell Reference>");
        }
        // need to extract bit after id
        string hashName = name.substr(MAX_ID_LENGTH);
        HandleSP handle(HandleSP::dynamicCast(handles->get(hashName)));
        
        if (!handle){
            throw ModelException(routine, "Unknown handle: "+name);
        }

        // check to see that handle is not stale
        if (strncmp(handle->idAsString, name.c_str(), MAX_ID_LENGTH)){
            throw ModelException(routine, "Stale handle ("+name+") used. "
                                 "The handle has "
                                 "been updated elsewhere on the spreadsheet");
        }
        
        obj = handle->object;
        if (objectClass){
            CClassConstSP clazz = obj->getClass(); // save before any conversion
            // be flexible for data dictionaries etc
            try{
                CObject::checkType(obj, objectClass);
            } catch (exception& e){
                throw ModelException(e,
                                     routine, "Type mismatch. Expected "
                                     "handle of type "+objectClass->getName()+
                                     " but got handle of type "+
                                     clazz->getName());
            }
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
    return obj;
}

/** Removes the specified handle from storage */
void Handle::destroy(const string&  name){
    // check that we've got a handle
    if (valid(name)){
        // need to extract bit after id
        string hashName = name.substr(MAX_ID_LENGTH);
        handles->remove(hashName);
    }        
}

/** Removes all handles from storage */
void Handle::destroyAll(){
    handles->clear();
}

Handle::IName::~IName(){}

static void myLoad(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(Handle, clazz);
    SUPERCLASS(CObject);
}


CClassConstSP const Handle::TYPE = CClass::registerClassLoadMethod(
    "Handle", typeid(Handle), myLoad);

/** returns the string - which could be a handle name */
const string& Handle::RawString::getString() const{
    return rawString;
}

void Handle::RawString::load(CClassSP& clazz){
    REGISTER(Handle::RawString, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultCreate);
    FIELD(rawString, "unprocessed string");
}
IObject* Handle::RawString::defaultCreate(){
    return new RawString();
}

Handle::RawString::RawString(): CObject(TYPE){}

/** constructor */
Handle::RawString::RawString(const string& rawString): 
    CObject(TYPE), rawString(rawString){}

CClassConstSP const Handle::RawString::TYPE = CClass::registerClassLoadMethod(
    "Handle::RawString", typeid(Handle::RawString), load);

typedef Handle::RawStringArray HandleRawStringArray;
//template<> CClassConstSP const HandleRawStringArray::TYPE =
//CClass::registerClassLoadMethod("Handle::RawStringArray",
//                                typeid(Handle::RawStringArray), load);
DEFINE_TEMPLATE_TYPE_WITH_NAME("Handle::RawStringArray", HandleRawStringArray);

//// 'actions' for manipulating handles through DR Interface
class Handle::FromObject: public CObject,
              public virtual ClientRunnable{
    static CClassConstSP const TYPE;

    //// parameters
    string    userName; // 'prefix' for handle
    IObjectSP object;   // the object to create the handle to

    IObjectSP run(){
        return IObjectSP(CString::create(Handle::create(userName, object)));
    }
    /** for reflection */
    FromObject():  CObject(TYPE){}
    static IObject* defaultConstructor(){
        return new FromObject();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(FromObject, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(userName, "'prefix' for handle");
        FIELD(object, "the object to create the handle to");
    }
};
CClassConstSP const Handle::FromObject::TYPE = 
CClass::registerClassLoadMethod(
    "Handle::FromObject", typeid(FromObject), load);

//// 'actions' for manipulating handles through DR Interface
class Handle::ToObject: public CObject,
              public virtual ClientRunnable{
    static CClassConstSP const TYPE;

    //// parameters
    string    name; // handle name as returned by FromObject

    IObjectSP run(){
        return Handle::fetch(name, 0);
    }
    /** for reflection */
    ToObject():  CObject(TYPE){}
    static IObject* defaultConstructor(){
        return new ToObject();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(ToObject, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(name, "handle name as returned by FromObject");
    }
};
CClassConstSP const Handle::ToObject::TYPE = 
CClass::registerClassLoadMethod(
    "Handle::ToObject", typeid(ToObject), load);

/** class for taking an array of 'raw' strings */
class DeleteAllHandles: public CObject{
public:
    static CClassConstSP const TYPE;

    DeleteAllHandles(): CObject(TYPE){}

    static bool run(DeleteAllHandles* params){
        Handle::destroyAll();
        return true;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DeleteAllHandles, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDeleteAllHandles);
        Addin::registerClassBoolMethod("DELETE_ALL_HANDLES",
                                       Addin::UTILITIES,
                                       "Deletes all handles",
                                       TYPE,
                                       (Addin::BoolMethod*)run);

    }

    static IObject* defaultDeleteAllHandles(){
        return new DeleteAllHandles();
    }
};

CClassConstSP const DeleteAllHandles::TYPE = CClass::registerClassLoadMethod(
    "DeleteAllHandles", typeid(DeleteAllHandles), load);

DRLIB_END_NAMESPACE
