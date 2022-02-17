//////////////////////////////////////////////////////////////
// Equity Derivatives Research: Models Library Interface
// Originally from EDRInterface.cpp
//
//
//////////////////////////////////////////////////////////////

#include "edginc/config.hpp"
#include "edginc/DataDictionary.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Library.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Null.hpp"
#include "edginc/Version.hpp"
#include "edginc/GDRDate.hpp"
#include "edginc/XObject.hpp"
#include "edginc/PrivateObject.hpp"
#include "edginc/Map.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/EDGServices.h"
#include "edginc/DRUtil.hpp"
#include "edginc/Lock.hpp"
#include <algorithm>
#include <sstream>

DRLIB_BEGIN_NAMESPACE
#define EDG_ASSERT(x, errmsg)                                     \
    {                                                             \
        if (!(x)) {                                               \
            if (errmsg) {                                         \
                throwError(method, errmsg, service);              \
            }                                                     \
            else {                                                \
                throwError(method, #x" not specified", service);  \
            }                                                     \
        }                                                         \
    }

// DRI_CHECK is defined in DRUtil.hpp. This simplifies the invocation of the
// DRI_CHECK macro by fixing some parameters.
#define CHECK(svc, exec, userMsg)                                       \
    DRI_CHECK(exec, method.c_str(), string(userMsg).c_str(), handleError, \
              svc, ;, ;, ;)

// DRI_ERROR_CALLBACK (DRUtil.hpp)
static void handleError(const char *method,
                        const char *errMsg,
                        void       *cbParam);

static void throwError(const string &method,
                       const string &errMsg,
                       void         *cbParam);

// Return the stack trace of a ModelException with out-of-memory provisions. 
static DRError EDGSExceptionStackTrace(const ModelException& exception){
    char *err = exception.stackTrace();
    return err ? err : driOutOfMemoryError();
}

// this is the internal representation of DRObject
// has a service pointer so we can identify EDG flavour DRObjects
class ObjectHandle {
  private:
    static size_t       idCounter;
    static const string thisClassName;
    
  public:
    // constructor from smart pointer
    ObjectHandle(DRService *svc, const IObjectSP& object): 
        svc(svc), object(object){
        if (!svc) {
            throw ModelException(thisClassName, "NULL service supplied.");
        }
        if (!object){
            throw ModelException(thisClassName, "NULL pointer supplied.");
        }
    }

    // constructor from pointer
    ObjectHandle(DRService *svc, IObject *object)
    : svc(svc)
    , object(object)
    , objectId(idCounter++)
    {
        if (!svc) {
            throw ModelException(thisClassName, "NULL service supplied.");
        } 
        if (!object){
            throw ModelException(thisClassName, "NULL pointer supplied.");
        }
    }

    // fields
    DRService*    svc; 
    IObjectSP     object;
    size_t        objectId;
};

size_t ObjectHandle::idCounter = 0;
const string ObjectHandle::thisClassName = "ObjectHandle";

// this is the internal representation of DRService
class EDGService {
  public:
    // the function pointers that make up the interface
    DRServiceInterface_  *fptr;

    // service is a "singleton" - i.e. "there can be only one", whether it be
    // the thread-safe version or not.
    static EDGService    *svc;

  private:
    static DRServiceInterface_ REGULAR_DRFUNCS;

  protected:
    // name of EDGService
    string                      svcName;

    static int                  svcCount;    // ref count the singleton

  public:
    EDGService(DRServiceInterface_ *funcTable = 0, 
               const char          *serviceName = 0)
    : fptr(funcTable ? funcTable : &REGULAR_DRFUNCS)
    , svcName(serviceName ? serviceName : Library::SERVICE_NAME.c_str())
    {
        // Do not try to update svcCount here. EDGService must be
        // constructable and deletable as a baseclass for ThreadSafeEDGService.
    }

    // This is a hack.  To support typecasting of both EDGService and
    // ThreadSafeEDGService to DRService, fptr needs to be the first data
    // member in objects of each class.  While it is carefully maintained that
    // the derived class ThreadSafeEDGService does not have any data members, a
    // vtable pointer will still offset fptr.  So, we made the destructor of
    // EDGService non-virtual, and use our own lookup to determine if the
    // destructor of the derived object should be called.  All destruction of
    // EDGService must therefore go through serviceFree to avoid memory leaks.
    // In any case, EDGService is usually killed at the end of the program.
    // (See serviceFree.)
    ~EDGService() {}

    // Methods added to support ThreadSafeEDGService
    const string& serviceName() const { return svcName; }

    // Singleton 
    static EDGService* getService() {
        if (svcCount == 0) {
            svc = new EDGService();
        }
        else if (isThreadSafe(svc)) {
            // Only allow either thread-safe or non-thread-safe service.
            return 0;
        }
        ++svcCount;
        return svc;
    }

    // Decrement reference count and return true if the service gets deleted.
    // Note that deleteService is not virtual because we are using function
    // hiding.
    static bool deleteService() {
        if (--svcCount == 0) {
            delete svc;
            return true;
        }
        return false;
    }

    static bool isThreadSafe(EDGService *svc) {
        return svc->serviceName() == Library::THREADSAFE_SERVICE_NAME;
    }

    // Is this our service?
    static DRError checkEdgService(DRService *service)
    {
        // Only compare base class because (for now), the objects are agnostic
        // to the type of EDGService they are passed into. Also, as a bonus of
        // not making ~EDGService virtual, typeid will always return the
        // base-class typeinfo.
        if (!service || service != (DRService*) svc) {
            return "service is not EDG";
        }
        
        return 0;
    }

    // was the DRObject created by the EDG service
    static DRError checkEdgObject(DRObject object)
    {
        if (!object) {
            //return "object not specified";
            return 0;   // from the old implementation
        }
        return checkEdgService(object->svc);
    }

    // is this DRObject a DRMapIterator?
    static bool isDRMapIterator(DRObject object)
    {
        IObject* iterObj = EDGService::drObject2Object(object);
        return IMap::IIterator::TYPE->isInstance(iterObj);
    }

    // convert a DRValue to an IObject
    static IObjectSP drValue2Object(const DRValue *drValue) 
    {
        static const string method = "drValue2Object";
        try {
            if (!drValue) {
                throw ModelException(method, "DRValue is null.");
            }
                
            switch (drValue->type)
            {
              case DR_BOOL:
                return IObjectSP(CBool::create(
                                     drValue->value.boolean == DR_TRUE));
              case DR_INT:
                return IObjectSP(CInt::create(drValue->value.integer));
              case DR_DOUBLE:
                return IObjectSP(CDouble::create(drValue->value.real));
              case DR_STRING:
                return IObjectSP(CString::create(drValue->value.string));
              case DR_OBJECT:
                if (!drValue->value.object){
                    return CNull::create();
                }
                return IObjectSP(drObject2Object(drValue->value.object));
              case DR_DATE:
                //throw ModelException(
                //method, "DR_DATE type is obsoleted in DRI 1.2 for QLIB. "
                //"Please use the DateTime object instead.");

                // For backward compatibility with Kapital and other apps,
                // we still need to support GDRDate.
                return IObjectSP(GDRDate::create(drValue->value.date));
              case DR_UNDEFINED:
                // DRI 1.2
                return IObjectSP(CNull::create());
                //throw ModelException(method, 
                //"DR_VALUE type is DR_UNDEFINED");
              default:
              {
                  throw ModelException(method,
                                       "Unrecognised DR_VALUE type (" + 
                                       Format::toString(drValue->type) + ").");
             
              }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // convert an IObject to a DRValue
    void object2DRValue(IObject *object, DRValue *drValue)
    {
        static const string method = "object2DRValue";
        try {
            if (!object){
                drValue->type = DR_UNDEFINED;  // used to be DR_OBJECT;
                drValue->value.object = 0;
            }
            else {                
                CClassConstSP clazz = object->getClass();
                if (clazz->isPrimitive()){
                    if (clazz == CDouble::TYPE){
                        drValue->type = DR_DOUBLE;
                        CDouble* obj = STATIC_CAST(CDouble, object);
                        drValue->value.real = obj->doubleValue();
                    } 
                    else if (clazz == CInt::TYPE){
                        drValue->type = DR_INT;
                        CInt* obj = STATIC_CAST(CInt, object);
                        drValue->value.integer = obj->intValue();
                    } 
                    else if (clazz == CBool::TYPE){
                        drValue->type = DR_BOOL;
                        CBool* obj = STATIC_CAST(CBool, object);
                        drValue->value.boolean = obj->boolValue();
                    } 
                    else if (clazz == CString::TYPE){
                        drValue->type = DR_STRING;
                        CString* sp = STATIC_CAST(CString, object);
                        drValue->value.string = 
                            driStringCreate(sp->stringValue().c_str());
                    } 
                    else if (clazz == GDRDate::TYPE){
                        drValue->type = DR_DATE;
                        GDRDate* dt = STATIC_CAST(GDRDate, object);
                        drValue->value.date = dt->dateValue();
                    } 
                    else if (clazz->isEnum()){
                        const Enum* theEnum = STATIC_CAST(Enum, object);
                        const string& enumAsString=theEnum->enumValueAsString();
                        drValue->type = DR_STRING;
                        drValue->value.string =
                            driStringCreate(enumAsString.c_str());
                    } else {
                        throw ModelException(method, "Internal error.");
                    }
                } 
                else {
                    drValue->type = DR_OBJECT;
                    if (XObject::TYPE->isAssignableFrom(clazz)){
                        // return the external object that's being wrapped
                        XObject* xObj = STATIC_CAST(XObject, object);
                        drValue->value.object = xObj->getDRObject();
                    } 
                    else if (clazz == CNull::TYPE){
                        drValue->type = DR_UNDEFINED; 
                        drValue->value.object = 0;
                    }
                    else {
                        drValue->value.object = (DRObject)new 
                            ObjectHandle((DRService*) EDGService::svc, object);
                    }
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // convert an IObject to a DRObject
    DRObject object2DRObject(IObject *object)
    {
        return (DRObject)new ObjectHandle((DRService*)svc, object);
    }      
    
    // convert a DRObject to an IObject
    static IObject* drObject2Object(DRObject drObj)
    {
         if (!drObj){
            return 0;
        }
        bool isEDR = !checkEdgObject(drObj);  // no errmsg - must be EDR
        if (isEDR){
            ObjectHandle* handle = (ObjectHandle*)drObj;
            return handle->object.get();
        }
        else { 
            // An object from another library(!). Wrap up in an XObject
            // we get a reference to objects - therefore we do not own
            // the memory - so pass false to toIDRObject()
            XObjectSP theObject(XObject::toXObject(drObj, false));
            // must clone because we have a reference
            return theObject->clone();
        }
    }
};

// Do not use inheritance because the vtable pointer will offset fptr which
// needs to come first to enable the typecasting of EDGService to DRService,
// and vice versa.  Also, do not add any data members to this derived class.
class ThreadSafeEDGService : public EDGService {
  private:
    static DRServiceInterface_ THREADSAFE_DRFUNCS;

  public:
    static Lock        lock;        // one global lock for every function call

    ThreadSafeEDGService(DRServiceInterface_ *funcTable = 0,
                         const char          *serviceName = 0)
    : EDGService(funcTable ? funcTable : &THREADSAFE_DRFUNCS,
                 serviceName ? 
                 serviceName : Library::THREADSAFE_SERVICE_NAME.c_str())
    {}

    // Singleton 
    static EDGService* getService() {
        if (svcCount == 0) {
            svc = new ThreadSafeEDGService();
        }
        else if (!isThreadSafe(svc)) {
            // Only allow either thread-safe or non-thread-safe service.
            return 0;
        }
        ++svcCount;
        return svc;
    }
};

EDGService* EDGService::svc = 0;
int EDGService::svcCount = 0;
Lock ThreadSafeEDGService::lock = Lock();

//----------------------------------------------------------------------------
// These functions actually implement the DR Interface.
//
// As a C interface, no exceptions are thrown in any function. 
//----------------------------------------------------------------------------

// Extract the state of the specified DRService by means of having this method
// invoke the specified createServiceArgs callback, given the pointer to this
// DRService instance.  This method is useful for storing and varying the
// parameters used to create a service (without ownership issues because of the
// callback), cloning a service with these (varied) parameters, or recreating a
// service when the state of the given service is an unknown (e.g., in the
// situation where the service is to be recreated remotely).
//
// Example createServiceArgs callbacks that client could put in his/her code in
// preparation for calling serviceCreateArgs. 
//
// DRServiceInitArgs initArgs;
// static void createServiceArgs(void              *userData,
//                               DRServiceInitArgs *args)
// {
//     *initArgs = *args;
// }
//
static DRError serviceCreateArgs(DRService           *service,
                                 DRCreateServiceArgs  createServiceArgs,
                                 void                *userData)
{
    static const string method = "serviceCreateArgs";

    try {
        CHECK(service, EDGService::checkEdgService(service), "");
        EDG_ASSERT(createServiceArgs, 0);

        // Initialize parameters of DRServiceInitArgs pertaining to this
        // service.
        DRServiceInitArgs args;
        EDGSInitServiceCreateArgs(&args);

        // Pass this DRServiceInitArgs to the createServiceArgs callback, which
        // will copy this structure to the output container and add
        // initialization of parameters that only the client can specify, e.g.,
        // errorHandler and its first parameter.
        (*createServiceArgs)(userData, &args);

        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// The EDR service description
static DRError serviceDescription(DRService *service, 
                                  DRString  *description) 
{
    static const string method = "serviceDescription";

    try {
        CHECK(service, EDGService::checkEdgService(service), "");
        EDG_ASSERT(description, 0);
        *description = 0;

        string qlibVer = CVersion::DRLibVersion();
        return service->fptr->stringNew(service,
                                        qlibVer.c_str(),
                                        description);
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Kill service 
static DRError serviceFree(DRService *service)
{
    // This is a hack.  To support typecasting of both EDGService and
    // ThreadSafeEDGService to DRService, fptr needs to be the first data
    // member in objects of each class.  While it is carefully maintained that
    // the derived class ThreadSafeEDGService does not have any data members, a
    // vtable pointer will still offset fptr.  So, we made the destructor of
    // EDGService non-virtual, and use our own lookup to determine if the
    // destructor of the derived object should be called.  All destruction of
    // EDGService must therefore go through serviceFree to avoid memory leaks.
    // In any case, EDGService is usually killed at the end of the program.
    //
    // (The lookup isn't so bad anyway as we have done the same name check
    // when DRCreateService decides which service to create.)

    static const string method = "serviceFree";

    try {
        CHECK(service, EDGService::checkEdgService(service), "");

        EDGService *svc = (EDGService*) service;
        svc->deleteService();
        Library::shutdown();
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Run an executable DRObject i.e 'Client Runnable' 
static DRError execute(DRService *service,
                       DRObject   input,
                       DRValue   *output) 
{
    static const string method = "execute";

    try {
        EDG_ASSERT(input, 0);
        EDG_ASSERT(output, 0);
        DRI_VALUE_CLEAR(*output);

        IObject* object = EDGService::drObject2Object(input);
        if (!object) {
            throwError(method, "cannot convert input into an object.", 
                       service);
        }
            
        ClientRunnable* runner = dynamic_cast<ClientRunnable*>(object);

        if (!runner){
            throwError(method, 
                       "Input is of incorrect type (" +
                       object->getClass()->getName() + ")" +
                       " Object must be derived from ClientRunnable.",
                       service);
        }

        IObjectSP result(runner->run());

        ((EDGService*) service)->object2DRValue(result.get(), output);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Get the data type of an object 
static DRError objectGetType(DRService *service,
                             DRObject   object,
                             DRString  *typeName)
{
    static const string method = "objectGetType";

    try {
        EDG_ASSERT(typeName, 0);
        *typeName = 0;
        EDG_ASSERT(object, 0);

        // all QLib objects have a type
        IObject* obj = EDGService::drObject2Object(object);
        *typeName = driStringCreate(obj->getClass()->getName().c_str());
        EDG_ASSERT(*typeName, "string copy failed.");

        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Get the object ID
static DRError objectGetId(DRService  *service,
                           DRObject    object,
                           void      **id)
{
    static const string method = "objectGetId";

    try  {
        EDG_ASSERT(id, 0);
        *id = 0;
        CHECK(service, EDGService::checkEdgObject(object), "");
        EDG_ASSERT(!EDGService::isDRMapIterator(object), 
                   "DRMapIterator does not have an ID.");

        ObjectHandle* handle = (ObjectHandle*)object;
        *id = &handle->objectId;
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}    

// Free a DRObject 
static DRError objectFree(DRService *service, DRObject object) 
{
    static const string method = "objectFree";

    try{
        EDG_ASSERT(object, 0);
        CHECK(service, EDGService::checkEdgObject(object), "");
        ObjectHandle* handle = (ObjectHandle*)object;
        delete handle;
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Does an object represent an array ? 
static DRError objectIsArray(DRService *service,
                             DRObject   object,
                             DRBool    *isArray) 
{
    static const string method = "objectIsArray";
    
    try {
        EDG_ASSERT(object, 0);
        EDG_ASSERT(isArray, 0);

        IObject* obj = EDGService::drObject2Object(object);
        *isArray = obj->getClass()->isArray();
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Does an object represent a matrix of doubles ? 
static DRError objectIsMatrix(DRService *service,
                              DRObject   object,
                              DRBool    *isMatrix) 
{
    static const string method = "objectIsMatrix";
    
    try {
        EDG_ASSERT(object, 0);
        EDG_ASSERT(isMatrix, 0);

        IObject* obj = EDGService::drObject2Object(object);
        *isMatrix = CDoubleMatrix::TYPE->isInstance(obj);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Can execute() be called on the object ? Is it ClientRunnable ? 
static DRError objectIsExecutable(DRService *service,
                                  DRObject   object,
                                  DRBool    *isExecutable) 
{
    static const string method = "objectIsExecutable";
    
    try {
        EDG_ASSERT(object, 0);
        EDG_ASSERT(isExecutable, 0);

        IObject* obj = EDGService::drObject2Object(object);
        *isExecutable = ClientRunnable::TYPE->isInstance(obj);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

/* Convert an object into a map where its component parts are
   accessible as key/value pairs */
static DRError objectToMap(DRService *service,
                           DRObject   object,
                           DRMap     *drmap)
{
    static const string method = "objectToMap";

    try {
        EDG_ASSERT(object, 0);
        EDG_ASSERT(drmap, 0);
        *drmap = 0;

        IObjectSP obj(EDGService::drObject2Object(object));
        obj = CObject::toProxy(obj);
        // see if it's already a map in DR Interface terms
        IMap* edrmap = dynamic_cast<IMap*>(obj.get());
        if (edrmap && edrmap->isTrueMap()) {
            // just return what came in - except we clone it in case anyone
            // modifies it. Ideally should just do a clone of the map without
            // cloning the contents
            IObjectSP clone(obj->clone());
            *drmap = ((EDGService*) service)->object2DRObject(clone.get());
        }
        else {
            CDataDictionarySP dd(CDataDictionary::pop2DataDict(obj));
            *drmap = ((EDGService*) service)->object2DRObject(dd.release());
        }
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Convert a map of key/value pairs into an object 
static DRError mapToObject(DRService *service,
                           DRMap      drmap,
                           DRObject  *object) 
{
    static const string method = "mapToObject";

    try {
        CHECK(service, EDGService::checkEdgObject(drmap), "");
        EDG_ASSERT(object, 0);
        *object = 0;

        IObjectSP mapobj(EDGService::drObject2Object(drmap));
        if (CDataDictionary::TYPE->isInstance(mapobj)) {
            CDataDictionary* dict = STATIC_CAST(CDataDictionary, mapobj.get());
            IObjectSP obj(dict->pop2Object());
            obj = CObject::fromProxy(obj);
            *object = ((EDGService*) service)->object2DRObject(obj.release());
        }
        else if (!IMap::TYPE->isInstance(mapobj)) {
            throwError(method, "input map is not an EDG map.", service);
        }
        else {
            // must be a true map (hopefully by definition)
            // Just return what came in - except we clone it in case anyone
            // modifies it. Ideally should just do a clone of the map without
            // cloning the contents
            IObjectSP clone(mapobj->clone());
            *object = ((EDGService*) service)->object2DRObject(clone.get());
        }
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

/* Returns some version of an object if that version is supported
   this function is called by DR implementation 
   it should not be used directly */
static DRError objectVersionGet(DRService *service,
                                DRObject   object,
                                int       *devKitId,
                                int       *devKitVersion) 
{
    static const string method = "objectVersionGet";

    try {
        EDG_ASSERT(object, 0);
        EDG_ASSERT(devKitId, 0);
        EDG_ASSERT(devKitVersion, 0);

        *devKitId = DRI_LIBRARY_ID_EDG;
        *devKitVersion = 1;
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Create a DRString - i.e. a string copy
static DRError stringNew(DRService  *service,
                         const char *inString,
                         DRString   *outString) 
{
    static const string method = "stringNew";

    try {
        CHECK(service, EDGService::checkEdgService(service), "");
        EDG_ASSERT(inString, 0);
        EDG_ASSERT(outString, 0);

        *outString = driStringCreate(inString);
        EDG_ASSERT(*outString, "string copy failed.");

        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Free a DRString 
static DRError stringFree(DRService *service, DRString str) 
{
    static const string method = "stringFree";

    try {
        EDG_ASSERT(str, 0);

        // driStringFree contains out-of-memory provisions that the system
        // free does not have.
        driStringFree((char*) str);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Free any DRString or DRObject in a DRValue 
static DRError valueClear(DRService *service, DRValue *value)
{
    static const string method = "valueClear";

    try {
        EDG_ASSERT(value, 0);

        // No need to free the atomic types
        if (value->type == DR_OBJECT && value->value.object != 0) {
            CHECK(service, 
                  service->fptr->objectFree(service, value->value.object), 
                  "");
        }
        else if (value->type == DR_STRING && value->value.string != 0) {
            CHECK(service, 
                  service->fptr->stringFree(service, value->value.string), 
                  "");
        }
    
        DRI_VALUE_CLEAR(*value);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }    
}

/* Create an array of size numElements */
static DRError arrayNew(DRService  *service,
                        int         numElements,
                        const char *typeName,
                        DRArray    *drarray)
{
    static const string method = "arrayNew";

    try {
        EDG_ASSERT(typeName, 0);
        EDG_ASSERT(drarray, 0);
        *drarray = 0;

        // look up class of componentType
        CClassConstSP componentClass= CClass::forName(string(typeName));
        if (!Modifier::isPublic(componentClass->getModifiers())){
            throwError(method, string(typeName) +
                       " is not a public class - instantiation disallowed.",
                       service);
        }
        IArraySP edrarray(componentClass->
                          newArrayInstanceByComponent(numElements));

        *drarray = ((EDGService*) service)->object2DRObject(edrarray.release());
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// what's the type of an element of an array?
static DRError arrayElementType(DRService *service,
                                DRArray    array,
                                DRString  *typeName) 
{
    static const string method = "arrayElementType";

    try {
        EDG_ASSERT(array, 0);
        EDG_ASSERT(typeName, 0);
        *typeName = 0;

        IObject* object = EDGService::drObject2Object(array);
        CClassConstSP objClazz   = object->getClass();
        CClassConstSP arrayClazz = objClazz->getComponentType();

        string className = arrayClazz->getName();
        *typeName = driStringCreate(className.c_str());
        if (!*typeName) {
            throwError(method, "string copy failed.", service);
        }        
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// set the i'th item from an array (starts at 0) 
static DRError arraySet(DRService     *service,
                        DRArray        drarray,
                        int            index,
                        const DRValue *valueToSet)
{
    static const string method = "arraySet";

    try {
        EDG_ASSERT(drarray, 0);
        EDG_ASSERT(valueToSet, 0);

        IObjectSP object(EDGService::drObject2Object(drarray));
        IArraySP array(IArraySP::dynamicCast(object));

        IObjectSP elem(EDGService::drValue2Object(valueToSet));
        //// convert any public objects to private objects etc if needed
        CObject::checkType(elem, array->getClass()->getComponentType());

        if (!CNull::TYPE->isInstance(elem)){
            CClassConstSP desiredType = array->getClass()->getComponentType();
            if (desiredType != elem->getClass() ) {
                CObject::checkType(elem, desiredType);
            }
        }

        array->set(index, elem);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// get the i'th item from an array (starts at 0) 
static DRError arrayGet(DRService *service,
                        DRArray    drarray,
                        int        index,
                        DRValue   *itemHolder) 
{
    static const string method = "arrayGet";

    try {
        EDG_ASSERT(drarray, 0);
        EDG_ASSERT(itemHolder, 0);
        DRI_VALUE_CLEAR(*itemHolder);

        IObjectSP object(EDGService::drObject2Object(drarray));
        IArraySP array(IArraySP::dynamicCast(object));
        IObjectSP theElem(array->get(index));
        // make sure we don't give out private objects
        if (IPrivateObject::TYPE->isInstance(theElem)){
            theElem = CObject::convertToPublicRep(theElem);
        }
        IObjectSP elem(!theElem? 0: theElem->clone());
        ((EDGService*) service)->object2DRValue(elem.get(), itemHolder);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// get the number of elements in the array  
static DRError arrayLength(DRService *service, DRArray drarray, int *length) 
{
    static const string method = "arrayLength";

    try {
        EDG_ASSERT(drarray, 0);
        EDG_ASSERT(length, 0);
        *length = -1;

        IObjectSP object(EDGService::drObject2Object(drarray));
        IArraySP array(IArraySP::dynamicCast(object));
        *length = array->getLength();
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// create a new, empty matrix of doubles 
static DRError matrixNew(DRService *service,
                         int        numCols,
                         int        numRows,
                         DRMatrix  *matrix) 
{
    static const string method = "matrixNew";

    try {
        EDG_ASSERT(matrix, 0);
        *matrix = 0;

        // DRI 1.2
        if (numCols == 0 && numRows != 0) {
            throwError(method, "Construction of an nx0 matrix is denied.",
                       service);
        }

        if (numCols != 0 && numRows == 0) {
            throwError(method, "Construction of a 0xn matrix is denied.",
                       service);
        }

        *matrix = ((EDGService*) service)->object2DRObject(
            new CDoubleMatrix(numCols,numRows));
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// set the [i][j]th element of a matrix 
static DRError matrixSet(DRService *service,
                         DRMatrix   matrix,
                         int        col,
                         int        row,
                         double     value) 
{
    static const string method = "matrixSet";

    try {
        EDG_ASSERT(matrix, 0);

        IObjectSP object(EDGService::drObject2Object(matrix));
        CDoubleMatrixSP matrixSP(
            CDoubleMatrixSP::dynamicCast(object));
        (*matrixSP)[col][row] = value;
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// get the [i][j]th element of a matrix 
static DRError matrixGet(DRService *service,
                         DRMatrix   matrix,
                         int        col,
                         int        row,
                         double    *valueHolder) 
{
    static const string method = "matrixGet";

    try {
        EDG_ASSERT(matrix, 0);
        EDG_ASSERT(valueHolder, 0);
        *valueHolder = 0;

        IObjectSP object(EDGService::drObject2Object(matrix));
        CDoubleMatrixSP matrixSP(
            CDoubleMatrixSP::dynamicCast(object));
        *valueHolder = (*matrixSP)[col][row];
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// get the dimensions of a matrix  
static DRError matrixSize(DRService *service,
                          DRMatrix   matrix,
                          int       *numCol,
                          int       *numRow) 
{
    static const string method = "matrixSize";

    try {
        EDG_ASSERT(matrix, 0);
        EDG_ASSERT(numCol, 0);
        *numCol = -1;
        EDG_ASSERT(numRow, 0);
        *numRow = -1;

        IObjectSP object(EDGService::drObject2Object(matrix));
        CDoubleMatrixSP matrixSP(
            CDoubleMatrixSP::dynamicCast(object));
        *numCol = matrixSP->numCols();
        *numRow = matrixSP->numRows();
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// Create a new (empty) map to represent the "typename" class 
static DRError mapNewFromId(DRService *service,
                            void      *typeId, // this is a CClass*
                            DRMap     *drmap) 
{
    static const string method = "mapNewFromId";

    try {
        CHECK(service, EDGService::checkEdgService(service), "");
        EDG_ASSERT(typeId, 0);
        EDG_ASSERT(drmap, 0);
        *drmap = 0;

        const CClass* clazz = (const CClass*) typeId;
        const CClass* proxy = CObject::proxyClass(clazz);
        if (proxy) {
            clazz = proxy;
        }
        if (IMap::TYPE->isAssignableFrom(clazz)){
            // have a map
            IObjectSP map(clazz->newInstance());
            /* is it a true map? This is weak. Looks like the trueMap 
               attribute should be a property of the type. Perhaps some sort
               of marker interface */
            IMap* map_ptr = dynamic_cast<IMap*>(map.get());
            if (map_ptr && map_ptr->isTrueMap()){
                *drmap = 
                    ((EDGService*) service)->object2DRObject(map.release());
                return 0;
            } 
            else {
                // ... do nothing fall into code below (map is deleted)
            }
        }
        *drmap = ((EDGService*) service)->object2DRObject(
            CDataDictionary::create(clazz->getName().c_str()));      
        return 0;
    } 
    catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}

// Create a new (empty) map to represent the "typename" class 
static DRError mapNew(DRService  *service,
                      const char *typeName,
                      DRMap      *drmap)
{
    static const string method = "mapNew";
    try {
        EDG_ASSERT(typeName, 0);
        EDG_ASSERT(drmap, 0);
        *drmap = 0;

        const CClass* clazz = CClass::forName(typeName);
        return mapNewFromId(service, (void*) clazz, drmap);
    } 
    catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}

/* Same as mapNewFromValues but type provided by void* */
static DRError mapNewFromIdAndValues(DRService *service,
                                     void      *typeId,
                                     DRValue   *items,
                                     int        numItems,
                                     DRMap     *drmap)
{
    static const string method = "mapNewFromIdAndValues";
    
    try {
        CHECK(service, EDGService::checkEdgService(service), "");
        EDG_ASSERT(typeId, 0);
        EDG_ASSERT(drmap, 0);
        *drmap = 0;

        // conciously don't worry about proxies here - if something has
        // an exported interface then you are saying you are happy with
        // your interface for all time
        CClassConstSP clazz = (CClassConstSP) typeId;
        const string& typeName = clazz->getName();
        if (Modifier::isExport(clazz->getModifiers())) {
            // Do a bit of validation first
            CFieldArray fields(Addin::getDataClassFields(clazz));
            int numFields = fields.size();

            if (numFields != numItems) {
                throwError(method,
                           "number of fields for type " +
                           string(typeName) + " (" + 
                           Format::toString(numFields) + ")" + 
                           " does not equal " +
                           "numItems in parameter list!",
                           service);
            }

            CDataDictionarySP dict(CDataDictionary::create(typeName));

            for (int i = 0; i < numItems; i++) {
                const string fieldName = fields[i]->getName();
                IObjectSP object(EDGService::drValue2Object(&items[i]));
                dict->put(fieldName, object);
            }
   
            *drmap = ((EDGService*) service)->object2DRObject(dict.release()); 
            
        } else {
            throwError(method, "type " + typeName + 
                       " is NOT exported.", service);
        }
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

/* Calling mapNewFromValues should have the same effect as
   calling mapNew and then calling mapAddItem, where 
   "fieldName" argument  takes first "numItems" field names
   of the type that is being created.
   This function is supported only by types marked as exported 
*/
static DRError mapNewFromValues(DRService  *service,
                                const char *typeName,
                                DRValue    *items,
                                int         numItems,
                                DRMap      *drmap) 
{
    static const string method = "mapNewFromValues";
    try {
        EDG_ASSERT(typeName, 0);
        EDG_ASSERT(drmap, 0);
        *drmap = 0;

        const CClass* clazz = CClass::forName(typeName);
        return mapNewFromIdAndValues(service, (void*)clazz, items,
                                     numItems, drmap);
    } 
    catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}
 
//// this function lets user save time on type name lookups 
static DRError idFromName(DRService   *service,
                          const char  *typeName,
                          void       **id) 
{
    static const string method = "idFromName";
    
    try {
        EDG_ASSERT(typeName, 0);
        EDG_ASSERT(id, 0);
        *id = 0;
        *id = (void*) CClass::forName(typeName);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// add a value to a map 
static DRError mapAddItem(DRService     *service,
                          DRMap          drmap,
                          const char    *fieldName,
                          const DRValue *valueToAdd) 
{
    static const string method = "mapAddItem";

    try {
        CHECK(service, EDGService::checkEdgObject(drmap), "");
        EDG_ASSERT(fieldName, 0);
        EDG_ASSERT(valueToAdd, 0);

        IObjectSP dict(EDGService::drObject2Object(drmap));
        IObjectSP object(EDGService::drValue2Object(valueToAdd));
        IWriteableMap* dd = dynamic_cast<IWriteableMap*>(dict.get());
        if (!dd) {
            throwError(method, "Input map (" + 
                       dict->getClass()->getName() + 
                       ") is not a writeable map.", 
                       service);
        }

        dd->put(fieldName, object.get() ? object: IObjectSP(   ));
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

// get a value out of a map 
static DRError mapGetItem(DRService  *service,
                          DRMap       drmap,
                          const char *fieldName,
                          DRValue    *itemHolder)
{
    static const string method = "mapGetItem";

    try {
        CHECK(service, EDGService::checkEdgObject(drmap), "");
        EDG_ASSERT(itemHolder, 0);
        DRI_VALUE_CLEAR(*itemHolder);

        EDG_ASSERT(fieldName, 0);

        IObjectSP dict(EDGService::drObject2Object(drmap));

        IReadableMap* dd = dynamic_cast<IReadableMap*>(dict.get());
        if (!dd) {
            throwError(method, "Input map (" + 
                       dict->getClass()->getName() + 
                       ") is not a readable map.", 
                       service);
        }

        // dd->get will throw an exception only if fieldName is invalid, and
        // output a CNull object when the value for fieldName is missing so
        // that object2DRValue can convert it into a DR_UNDEFINED DRValue.
        IObjectSP object(dd->get(fieldName));

        if (object.get() && object->getRefCount() != 1) {
            // held elsewhere - must clone
            object = IObjectSP(object->clone());
        }

        ((EDGService*) service)->object2DRValue(object.get(), itemHolder);
        return 0;
    }
    catch (exception& e) {
        itemHolder->type         = DR_UNDEFINED;
        itemHolder->value.object = 0;
        return EDGSExceptionStackTrace(e);
    }
}

// Returns the type of the map
static DRError mapGetTypeName(DRService *service,
                              DRMap      drmap,
                              DRString  *typeNameHolder)
{
    static const string method = "mapGetTypeName";

    try {
        CHECK(service, EDGService::checkEdgObject(drmap), "");
        EDG_ASSERT(typeNameHolder, 0);
        *typeNameHolder = 0;

        IObjectSP mapobj(EDGService::drObject2Object(drmap));
        CClassConstSP ddtype;
        if (CDataDictionary::TYPE->isInstance(mapobj)) {
            CDataDictionary* dict = STATIC_CAST(CDataDictionary, mapobj.get());
            // it's possible that the map represents a proxy for what the
            // client thinks they're building, but we want to hide this from
            // the client
            ddtype = dict->getType();
        } 
        else if (!IMap::TYPE->isInstance(mapobj)) {
            throwError(method, "input map is not an EDG map.", service);
        }
        else {
            // must be a true map (hopefully by definition)
            ddtype = mapobj->getClass();
        }
        CClassConstSP actual = CObject::actualClass(ddtype);
        string typeName = actual ? actual->getName() : ddtype->getName();
        if (!(*typeNameHolder = driStringCreate(typeName.c_str()))) {
            throwError(method, "string copy failed.", service);
        }
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

/** get array of values from the map.
    Values are returned in the same order as corresponding 
    data members in the type represnted by the map.
    Omitted values have type DR_UNDEFINED.
    This function is supported only by types marked as exported.
    Calling this function for other types will set *array to 0.  
*/
// guess what - CMLib !
static DRError mapGetValues(DRService *service,
                            DRMap      drmap,
                            DRArray   *drarray) 
{
    static const string method = "mapGetValues";

    try {
        CHECK(service, EDGService::checkEdgObject(drmap), "");
        EDG_ASSERT(drarray, 0);
        *drarray = 0;

        // handy if you really call functions 
        // rather than build up classes via reflection
        IObjectSP mapobj(EDGService::drObject2Object(drmap));
        if (CDataDictionary::TYPE->isInstance(mapobj)){
            // find out what this map is meant to be
            // again, no need to worry about proxies as this is a guaranteed
            // interface that you're happy with for ever
            CDataDictionary* dict = STATIC_CAST(CDataDictionary, mapobj.get());
            CClassConstSP clazz = dict->getType();
       
            if (Modifier::isExport(clazz->getModifiers())) {
                CFieldArray fields(Addin::getDataClassFields(clazz));
                int numFields = fields.size();
                
                // make 'variant' array - it's going to contain arbitrary data
                CHECK(service, service->fptr->arrayNew(service, numFields,
                                              DRI_TYPE_VARIANT, drarray),
                      "array build failed.");
                
                // can't use an iterator - the order 
                // we get the fields is important
                for (int i = 0; i < numFields; i++) {
                    DRValue obj;
                    const string fieldName = fields[i]->getName();
                    
                    CHECK(service, service->fptr->mapGetItem(service, 
                                                    drmap, 
                                                    fieldName.c_str(),
                                                    &obj), "");
                    
                    CHECK(service, 
                          service->fptr->arraySet(service, *drarray, i, &obj), 
                          "");
                    service->fptr->valueClear(service, &obj);
                }
            }
            else {
                *drarray = 0; // apparently we don't fail here
            }
        }
        else if (!IMap::TYPE->isInstance(mapobj)) {
            throwError(method, "input map is not an EDG map.", service);
        } 
        else {
            *drarray = 0; // apparently we don't fail here
        }
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

/* Creates a map iterator (which must eventually be freed)
   This object can then be used to explore the contents of a map */
static DRError mapIteratorGet(DRService     *service,
                              DRMap          drmap,
                              DRMapIterator *mapIterator)
{
    static const string method = "mapIteratorGet";

    try {
        CHECK(service, EDGService::checkEdgObject(drmap), "");
        EDG_ASSERT(mapIterator, 0);
        *mapIterator = 0;

        IObjectSP dict(EDGService::drObject2Object(drmap));
        IMap* dd = dynamic_cast<IMap*>(dict.get());
        if (!dd) {
            throwError(method, 
                        "Input (" + dict->getClass()->getName() + 
                        ") is not a map.",
                        service);
        }
       
        *mapIterator = ((EDGService*) service)->object2DRObject(
            dd->createIterator().get());
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

/* Returns next <field, object> pair in iteration 
   field and object will be null when there are no more elements 
   in the iteration. */
static DRError mapIteratorNext(DRService     *service,
                               DRMapIterator  mapIterator,
                               DRString      *fieldName,
                               DRValue       *objectHolder)
{
    static const string method = "mapIteratorNext";

    try
    {
        CHECK(service, EDGService::checkEdgObject(mapIterator), "");
        EDG_ASSERT(objectHolder, 0);
        DRI_VALUE_CLEAR(*objectHolder);

        EDG_ASSERT(fieldName, 0);

        IObject* iterobj = EDGService::drObject2Object(mapIterator);
        
        IMap::IIterator* dditer = dynamic_cast<IMap::IIterator*>(iterobj);
        if (!dditer) {
            throwError(method, "invalid iterator.", service);
        }

        if (!dditer->hasMoreElements()){
            *fieldName = 0;
            objectHolder->type = DR_UNDEFINED;
        } 
        else {
            // don't like this being dynamic
            *fieldName = driStringCreate(dditer->getKey().c_str());
            IObject* theObj = dditer->getElement().get();
            IObjectSP item(theObj ? theObj->clone(): 0); // to do: review
            ((EDGService*) service)->object2DRValue(item.get(), objectHolder);
            dditer->increment();
        }
        return 0;
    } 
    catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}

// This returns an object representing the type description 
// of the "typename" class.
static DRError typeDescription(DRService  *service,
                               const char *typeName,
                               DRObject   *typeObjectHolder) 
{
    static const string method = "typeDescription";

    try {
        EDG_ASSERT(typeObjectHolder, 0);
        *typeObjectHolder = 0;

        EDG_ASSERT(typeName, 0);

        // look up class
        CClassConstSP clazz = CClass::forName(string(typeName));
        // see if it's a proxy
        if (CObject::actualClass(clazz) != 0){
            // proxies do not exist as far as the client is concerned
            throwError(method, 
                       "Class not found: " + string(typeName), 
                       service);
        }
        IObjectSP driType(clazz->getDRIType()); // get hold of DRIType
        *typeObjectHolder = 
            (DRObject)new ObjectHandle(service, driType);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}


/** Return a list of all data types that can be created via the service
    i.e. classes built via maps */
static DRError typeList(DRService *service, DRArray *typeArray) 
{
    static const string method = "typeList";

    try {
        EDG_ASSERT(typeArray, 0);
        *typeArray = 0;

        const CClassVec& allClasses = CClass::allClasses();
        unsigned int numTypes = allClasses.size();
        unsigned int i;
        // return a list of strings that are the class names
        StringArraySP array(new StringArray(0));
        
        for (i = 0; i < numTypes; i++) {
            // hide private types and proxies from the Client
            if (!Modifier::isPrivate(allClasses[i]->getModifiers()) &&
                CObject::actualClass(allClasses[i]) == 0) {
                array->push_back(allClasses[i]->getName());
            }
        }
        sort(array->begin(), array->end());

        // apparently we don't return type names, but an array of
        // all the reflection information!
        CHECK(service, service->fptr->arrayNew(service, array->size(),
                                               DRI_TYPE_VARIANT, typeArray), 
              "couldn't build output array");
        
        DRValue vinfo;
        for (i = 0; i < (unsigned int)array->size(); i++) {
            DRObject info = 0;
            CHECK(service, 
                  service->fptr->typeDescription(service, (*array)[i].c_str(),
                                                 &info), 
                  "couldn't get info for " + (*array)[i]);

            vinfo.type         = DR_OBJECT;
            vinfo.value.object = info;
            CHECK(service, 
                  service->fptr->arraySet(service, *typeArray, i, &vinfo), 
                  "couldn't add info for " + (*array)[i]);
            valueClear(service, &vinfo);
        }

        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}



// The structure of service function pointers
DRServiceInterface_ EDGService::REGULAR_DRFUNCS = 
{
    DRLIB_VERSION(1,2,0,0),
    DRI_LIBRARY_ID_QLIB,
    &serviceCreateArgs,
    &serviceDescription,
    &serviceFree,
    &execute,
    &objectGetType,
    &objectGetId,
    &objectFree,
    &objectIsArray,
    &objectIsMatrix,
    &objectIsExecutable,
    &objectToMap,
    &mapToObject,
    &objectVersionGet,
    &stringNew,
    &stringFree,
    &valueClear,
    &arrayNew,
    &arrayElementType,
    &arraySet,
    &arrayGet,
    &arrayLength,
    &matrixNew,
    &matrixSet,
    &matrixGet,
    &matrixSize,
    &mapNew,
    &mapNewFromValues,
    &idFromName,
    &mapNewFromId,
    &mapNewFromIdAndValues,
    &mapAddItem,
    &mapGetItem,
    &mapGetTypeName,
    &mapGetValues,
    &mapIteratorGet,
    &mapIteratorNext,
    &typeList,
    &typeDescription 
};

//
// Thread-safe version
//
#define MAKE_THREAD_SAFE(f, arg)                          \
    {                                                     \
        Guard lockGuard(&ThreadSafeEDGService::lock);     \
        return ((f)arg);                                  \
    }

static DRError ts_serviceCreateArgs(DRService           *service,
                                    DRCreateServiceArgs  createServiceArgs,
                                    void                *userData)
{
    MAKE_THREAD_SAFE(serviceCreateArgs, 
                     (service, createServiceArgs, userData));
}

// The EDR service description
static DRError ts_serviceDescription(DRService *service, 
                                     DRString  *description) 
{
    MAKE_THREAD_SAFE(serviceDescription, (service, description));
}

// kill service 
static DRError ts_serviceFree(DRService *service) 
{
    MAKE_THREAD_SAFE(serviceFree, (service));
}

// Run an executable DRObject i.e 'Client Runnable' 
static DRError ts_execute(DRService *service,
                          DRObject   input,
                          DRValue   *output) 
{
    MAKE_THREAD_SAFE(execute, (service, input, output));
}

// Get the data type of an object 
static DRError ts_objectGetType(DRService *service,
                                DRObject   object,
                                DRString  *typeName) 
{
    MAKE_THREAD_SAFE(objectGetType, (service, object, typeName));
}

// Get the object ID
static DRError ts_objectGetId(DRService  *service,
                              DRObject    object,
                              void      **id)
{
    MAKE_THREAD_SAFE(objectGetId, (service, object, id));
}    

// Free a DRObject 
static DRError ts_objectFree(DRService *service, DRObject object) 
{
    MAKE_THREAD_SAFE(objectFree, (service, object));
}

// Does an object represent an array ? 
static DRError ts_objectIsArray(DRService *service,
                                DRObject   object,
                                DRBool    *isArray) 
{
    MAKE_THREAD_SAFE(objectIsArray, (service, object, isArray));
}

// Does an object represent a matrix of doubles ? 
static DRError ts_objectIsMatrix(DRService *service,
                              DRObject   object,
                              DRBool    *isMatrix) 
{
    MAKE_THREAD_SAFE(objectIsMatrix, (service, object, isMatrix));
}

// Can execute() be called on the object ? Is it ClientRunnable ? 
static DRError ts_objectIsExecutable(DRService *service,
                                     DRObject   object,
                                     DRBool    *isExecutable) 
{
    MAKE_THREAD_SAFE(objectIsExecutable, (service, object, isExecutable));
}

/* Convert an object into a map where its component parts are
   accessible as key/value pairs */
static DRError ts_objectToMap(DRService *service,
                              DRObject   object,
                              DRMap     *drmap)
{
    MAKE_THREAD_SAFE(objectToMap, (service, object, drmap));
}

// Convert a map of key/value pairs into an object 
static DRError ts_mapToObject(DRService *service,
                              DRMap      drmap,
                              DRObject  *object) 
{
    MAKE_THREAD_SAFE(mapToObject, (service, drmap, object));
}

/* Returns some version of an object if that version is supported
   this function is called by DR implementation 
   it should not be used directly */
static DRError ts_objectVersionGet(DRService *service,
                                   DRObject   object,
                                   int       *devKitId,
                                   int       *devKitVersion) 
{
    MAKE_THREAD_SAFE(objectVersionGet, 
                     (service, object, devKitId, devKitVersion));
}

// Create a DRString - i.e. a string copy
static DRError ts_stringNew(DRService  *service,
                            const char *inString,
                            DRString   *outString) 
{
    MAKE_THREAD_SAFE(stringNew, (service, inString, outString));
}

// Free a DRString 
static DRError ts_stringFree(DRService *service, DRString str) 
{
    MAKE_THREAD_SAFE(stringFree, (service, str));
}

// Free any DRString or DRObject in a DRValue 
static DRError ts_valueClear(DRService *service, DRValue *value)
{
    MAKE_THREAD_SAFE(valueClear, (service, value));
}

/* Create an array of size numElements */
static DRError ts_arrayNew(DRService  *service,
                           int         numElements,
                           const char *typeName,
                           DRArray    *drarray)
{
    MAKE_THREAD_SAFE(arrayNew, (service, numElements, typeName, drarray));
}

// what's the type of an element of an array?
static DRError ts_arrayElementType(DRService *service,
                                   DRArray    array,
                                   DRString  *typeName) 
{
    MAKE_THREAD_SAFE(arrayElementType, (service, array, typeName));
}


// set the i'th item from an array (starts at 0) 
static DRError ts_arraySet(DRService     *service,
                           DRArray        drarray,
                           int            index,
                           const DRValue *valueToSet)
{
    MAKE_THREAD_SAFE(arraySet, (service, drarray, index, valueToSet));
}


// get the i'th item from an array (starts at 0) 
static DRError ts_arrayGet(DRService *service,
                           DRArray    drarray,
                           int        index,
                           DRValue   *itemHolder) 
{
    MAKE_THREAD_SAFE(arrayGet, (service, drarray, index, itemHolder));
}

// get the number of elements in the array  
static DRError ts_arrayLength(DRService *service, DRArray drarray, int *length)
{
    MAKE_THREAD_SAFE(arrayLength, (service, drarray, length));
}


// create a new, empty matrix of doubles 
static DRError ts_matrixNew(DRService *service,
                            int        numCols,
                            int        numRows,
                            DRMatrix  *matrix) 
{
    MAKE_THREAD_SAFE(matrixNew, (service, numCols, numRows, matrix));
}

// set the [i][j]th element of a matrix 
static DRError ts_matrixSet(DRService *service,
                            DRMatrix   matrix,
                            int        col,
                            int        row,
                            double     value) 
{
    MAKE_THREAD_SAFE(matrixSet, (service, matrix, col, row, value));
}

// get the [i][j]th element of a matrix 
static DRError ts_matrixGet(DRService *service,
                            DRMatrix   matrix,
                            int        col,
                            int        row,
                            double    *valueHolder) 
{
    MAKE_THREAD_SAFE(matrixGet, (service, matrix, col, row, valueHolder));
}


// get the dimensions of a matrix  
static DRError ts_matrixSize(DRService *service,
                             DRMatrix   matrix,
                             int       *numCol,
                             int       *numRow) 
{
    MAKE_THREAD_SAFE(matrixSize, (service, matrix, numCol, numRow));
}

// Create a new (empty) map to represent the "typename" class 
static DRError ts_mapNewFromId(DRService *service,
                               void      *typeId, // this is a CClass*
                               DRMap     *drmap) 
{
    MAKE_THREAD_SAFE(mapNewFromId, (service, typeId, drmap));
}

// Create a new (empty) map to represent the "typename" class 
static DRError ts_mapNew(DRService  *service,
                         const char *typeName,
                         DRMap      *drmap)
{
    MAKE_THREAD_SAFE(mapNew, (service, typeName, drmap));
}

/* Same as mapNewFromValues but type provided by void* */
static DRError ts_mapNewFromIdAndValues(DRService *service,
                                        void      *typeId,
                                        DRValue   *items,
                                        int        numItems,
                                        DRMap     *drmap)
{
    MAKE_THREAD_SAFE(mapNewFromIdAndValues, 
                     (service, typeId, items, numItems, drmap));
}

/* Calling mapNewFromValues should have the same effect as
   calling mapNew and then calling mapAddItem, where 
   "fieldName" argument  takes first "numItems" field names
   of the type that is being created.
   This function is supported only by types marked as exported 
*/
static DRError ts_mapNewFromValues(DRService  *service,
                                   const char *typeName,
                                   DRValue    *items,
                                   int         numItems,
                                   DRMap      *drmap) 
{
    MAKE_THREAD_SAFE(mapNewFromValues,
                     (service,
                      typeName,
                      items,
                      numItems,
                      drmap));
}
 
//// this function lets user save time on type name lookups 
static DRError ts_idFromName(DRService   *service,
                             const char  *typeName,
                             void       **id) 
{
    MAKE_THREAD_SAFE(idFromName, (service, typeName, id));
}

// add a value to a map 
static DRError ts_mapAddItem(DRService     *service,
                             DRMap          drmap,
                             const char    *fieldName,
                             const DRValue *valueToAdd) 
{
    MAKE_THREAD_SAFE(mapAddItem, (service, drmap, fieldName, valueToAdd));
}

// get a value out of a map 
static DRError ts_mapGetItem(DRService  *service,
                             DRMap       drmap,
                             const char *fieldName,
                             DRValue    *itemHolder)
{
    MAKE_THREAD_SAFE(mapGetItem, (service, drmap, fieldName, itemHolder));
}

// Returns the type of the map
static DRError ts_mapGetTypeName(DRService *service,
                                 DRMap      drmap,
                                 DRString  *typeNameHolder)
{
    MAKE_THREAD_SAFE(mapGetTypeName, (service, drmap, typeNameHolder));
}

/** get array of values from the map.
    Values are returned in the same order as corresponding 
    data members in the type represnted by the map.
    Omitted values have type DR_UNDEFINED.
    This function is supported only by types marked as exported.
    Calling this function for other types will set *array to 0.  
*/
// guess what - CMLib !
static DRError ts_mapGetValues(DRService *service,
                               DRMap      drmap,
                               DRArray   *drarray) 
{
    MAKE_THREAD_SAFE(mapGetValues, (service, drmap, drarray));
}

/* Creates a map iterator (which must eventually be freed)
   This object can then be used to explore the contents of a map */
static DRError ts_mapIteratorGet(DRService     *service,
                                 DRMap          drmap,
                                 DRMapIterator *mapIterator)
{
    MAKE_THREAD_SAFE(mapIteratorGet, (service, drmap, mapIterator));
}

/* Returns next <field, object> pair in iteration 
   field and object will be null when there are no more elements 
   in the iteration. */
static DRError ts_mapIteratorNext(DRService     *service,
                                  DRMapIterator  mapIterator,
                                  DRString      *fieldName,
                                  DRValue       *objectHolder)
{
    MAKE_THREAD_SAFE(mapIteratorNext, 
                     (service, mapIterator, fieldName, objectHolder));
}

// This returns an object representing the type description 
// of the "typename" class.
static DRError ts_typeDescription(DRService  *service,
                                  const char *typeName,
                                  DRObject   *typeObjectHolder) 
{
    MAKE_THREAD_SAFE(typeDescription, (service, typeName, typeObjectHolder));
}

/** Return a list of all data types that can be created via the service
    i.e. classes built via maps */
static DRError ts_typeList(DRService *service, DRArray *typeArray) 
{
    MAKE_THREAD_SAFE(typeList, (service, typeArray));
}

// The structure of service function pointers
DRServiceInterface_ ThreadSafeEDGService::THREADSAFE_DRFUNCS = 
{
    DRLIB_VERSION(1,2,0,0),
    DRI_LIBRARY_ID_TS_QLIB,
    &ts_serviceCreateArgs,
    &ts_serviceDescription,
    &ts_serviceFree,
    &ts_execute,
    &ts_objectGetType,
    &ts_objectGetId,
    &ts_objectFree,
    &ts_objectIsArray,
    &ts_objectIsMatrix,
    &ts_objectIsExecutable,
    &ts_objectToMap,
    &ts_mapToObject,
    &ts_objectVersionGet,
    &ts_stringNew,
    &ts_stringFree,
    &ts_valueClear,
    &ts_arrayNew,
    &ts_arrayElementType,
    &ts_arraySet,
    &ts_arrayGet,
    &ts_arrayLength,
    &ts_matrixNew,
    &ts_matrixSet,
    &ts_matrixGet,
    &ts_matrixSize,
    &ts_mapNew,
    &ts_mapNewFromValues,
    &ts_idFromName,
    &ts_mapNewFromId,
    &ts_mapNewFromIdAndValues,
    &ts_mapAddItem,
    &ts_mapGetItem,
    &ts_mapGetTypeName,
    &ts_mapGetValues,
    &ts_mapIteratorGet,
    &ts_mapIteratorNext,
    &ts_typeList,
    &ts_typeDescription 
};

// DRI_ERROR_CALLBACK (DRUtil.hpp)
static void handleError(const char *method,
                        const char *errMsg,
                        void       *cbParam)
{
    EDGService *svc = reinterpret_cast<EDGService*>(cbParam);
    throw ModelException(svc->serviceName() + " Service:" + method, errMsg);
}

static void throwError(const string &method,
                       const string &errMsg,
                       void         *cbParam)
{
    EDGService *svc = reinterpret_cast<EDGService*>(cbParam);
    throw ModelException(svc->serviceName() + " Service:" + method, errMsg);
}

DRLIB_END_NAMESPACE

USING_DRLIB_NAMESPACE

// Initialize an error string for each of the service invocation error.
static bool initErrorStrings(vector<string>&  errStrs)
{
    static const char method[] = "DRCreateService";

    // Use a generic name "QLib" - let the caller determine which specific
    // QLib library should be used for output.
    const char *libName = "QLib"; 
    
    char buf[256];

    sprintf(buf, "%s cannot create %s service - library startup failure "
            "(error code = %d)", method, libName, DR_FAILED);
    errStrs[DR_FAILED] = buf;

    sprintf(buf, "%s cannot create %s service - invalid parameter "
            "(error code = %d)", method, libName, DR_INVALID_PARAMETER);
    errStrs[DR_INVALID_PARAMETER] = buf;
    
    sprintf(buf, "%s cannot create %s service - unknown service name "
            "(error code = %d)", method, libName, DR_UNKNOWN_SERVICE_NAME);
    errStrs[DR_UNKNOWN_SERVICE_NAME] = buf;

    sprintf(buf, "%s cannot create %s service - invalid input service "
            "version or internal service version cannot be parsed "
            "(error code = %d)", method, libName, DR_INVALID_SERVICE_VERSION);
    errStrs[DR_INVALID_SERVICE_VERSION] = buf;

    sprintf(buf, "%s cannot create %s service - mismatched service "
            "version: %u.%u expected (error code = %d)",
            method, libName, 
            DRI_VERSION >> 24, (DRI_VERSION & 0x00FF0000) >> 16,
            //(DRI_VERSION & 0x0000FF00) >> 8, DRI_VERSION & 0x000000FF, 
            DR_INVALID_INTERFACE_VERSION);
    errStrs[DR_INVALID_INTERFACE_VERSION] = buf;

    sprintf(buf, "%s cannot create %s service - missing startup arguments "
            "(error code = %d)", method, libName, DR_MISSING_OPTION);
    errStrs[DR_MISSING_OPTION] = buf;

    return true;
}

// Define and return static error strings (so that the errorHandler if invoked
// in DRCreateService won't have to delete the error strings).
static const string& errorString(DR_ERROR errorCode)
{
    // Use static strings for return because errorHandler is not permitted
    // to free error strings.
    static vector<string> errStrs(DR_MAX_ERROR + 1);
    static const bool initStrs = initErrorStrings(errStrs);
    return errStrs[errorCode];
}

// Exported function that actually starts up the EDG DRService
DRI_DLL_EXPORT DRBool DRCreateService(DRServiceInitArgs  *startupArgs, 
                                      DRService         **service)
{
    DR_ERROR    errorCode = DR_SUCCEEDED;
    string      libName = Library::SERVICE_NAME;
    EDGService *svc = 0;
    
    *service = 0;

    if (startupArgs->interfaceVersion != DRI_VERSION) {
        errorCode = DR_INVALID_INTERFACE_VERSION;
        goto invocationError;
    }
    
    if (!startupArgs) {
        errorCode = DR_MISSING_OPTION;
        goto invocationError;
    }

    // turn library on
    try {
        Library::startup();
    }
    catch (...) {
        errorCode = DR_FAILED;
        goto invocationError;
    }

    if (startupArgs->serviceVersion != Library::SERVICE_VERSION) {
        errorCode = DR_INVALID_SERVICE_VERSION;
        goto invocationError;
    }

    if (CString::equalsIgnoreCase(startupArgs->serviceName, 
                                  Library::SERVICE_NAME)) {
        libName = Library::SERVICE_NAME;
        if (!(svc = EDGService::getService())) {
            errorCode = DR_INVALID_PARAMETER;
            goto invocationError;
        }
    }
    else if (CString::equalsIgnoreCase(
                 startupArgs->serviceName, 
                 Library::THREADSAFE_SERVICE_NAME)) {
        libName = Library::THREADSAFE_SERVICE_NAME;
        if (!(svc = ThreadSafeEDGService::getService())) {
            errorCode = DR_INVALID_PARAMETER;
            goto invocationError;
        }
    }
    else {
        errorCode = DR_UNKNOWN_SERVICE_NAME;
        goto invocationError;
    }
    
    *service = (DRService*) svc;
    return true;

  invocationError:
    if (startupArgs->errorHandler) {
        (*startupArgs->errorHandler)(
            startupArgs->userData,
            (libName + ": " + errorString(errorCode)).c_str(),
            errorCode);
    }

    return false;
}

//// Allow EDRInterface access to these functions. A bit hacky but saves
//// a lot of work and (in theory....) the need for these will go soon
// convert an IObject to a DRValue
void object2DRValue(IObject *object, DRValue *drValue, DRService *service)
{
    ((EDGService*) service)->object2DRValue(object, drValue);
}

// convert a DRValue to an IObject
IObjectSP drValue2Object(const DRValue *drValue) 
{
    return EDGService::drValue2Object(drValue);
}

// convert a DRObject to an IObject
QLIB_DLL_EXPORT IObject* drObject2Object(DRObject drObj)
{
    return EDGService::drObject2Object(drObj);
}

// convert an IObject to a DRObject
QLIB_DLL_EXPORT DRObject object2DRObject(IObject *object, DRService *service)
{
    return ((EDGService*) service)->object2DRObject(object);
}

void isEdrObject(DRObject object)
{
    static const string method = "isEdrObject";
    CHECK(object->svc, EDGService::checkEdgObject(object), "");
}

DRLIB_BEGIN_NAMESPACE
// symbol (referenced by AddinLib.cpp) to ensure file gets linked in
bool DRInterfaceLink = true;
DRLIB_END_NAMESPACE
