//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : DRLibrary.cpp
//
//   Description : Class representing a reference to an external DR Library
//
//   Author      : Mark A Robson
//
//   Date        : 27 Nov 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_DRLIBRARY_CPP
#include "edginc/DRLibrary.hpp"
#include "edginc/DRAnalyticsInterface.h"
#include "edginc/Format.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Handle.hpp"
#include "edginc/XMap.hpp"
#include "edginc/CommandLineParams.hpp"
#include ext_hash_map
#if defined(_MSC_VER)
#include "windows.h"
#else
#include <dlfcn.h>
#if defined(sun)
#include <link.h>
#endif
#endif

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<DRLibrary>);

struct DRLibrary::LibHandle{
    // simple constructor
    LibHandle(): svc(0), osHandle(0){}
    // unload library when this handle goes out of scope
    ~LibHandle(){
        try{
            unloadLibrary();
        } catch (exception&){} // ignore
    }
    /** Force library to be unloaded */
    bool unloadLibrary(){
        if (osHandle){
            if (svc){
                /* free all handles - should stop crashes (at least on the
                   spreadsheet) */
                Handle::destroyAll();
                svc->fptr->serviceFree(svc);
                svc = 0;
            }
            if (
#if defined(_MSC_VER)
                FreeLibrary((HINSTANCE)osHandle)
#else
                dlclose(osHandle)
#endif
                == 0){
                throw ModelException("DRLibrary::unloadLibrary",
                                     "Unable to unload library (error code "+
#if defined(_MSC_VER)
                                     Format::toString((int)(GetLastError()))
#else
                                     string(dlerror())
#endif
                    );
            }
            osHandle = 0;
        }
        return true;
    }

    // fields
    DRService* svc;
    void*      osHandle;
};

/** Utility class for use in hashmap templates */
struct DRLibrary::Hash{
    //// equals method
    bool operator()(const DRLibrarySP& lib1, 
                    const DRLibrarySP& lib2) const{
        return (lib1->equals(lib2));
    }
    //// hash method
    size_t operator()(const DRLibrarySP& lib) const{
        return (hash_string(lib->libraryName) ^ hash_string(lib->serviceName) ^
                hash_string(lib->serviceVersion));
    }
    typedef hash_map<DRLibrarySP, refCountPtr<LibHandle>, Hash,Hash> LibraryMap;
};
/** Utility class for use in hashmap templates */
struct DRLibrary::ServiceHash{
    size_t operator()(DRService* svc) const {return (size_t)svc;}
    typedef hash_map<DRService*, DRLibrarySP, ServiceHash> Table;
    typedef hash_map<string, DRLibrarySP, Hashtable::StringHash> StringTable;
};


struct DRLibrary::Cache{
    static DRLibrary::Hash::LibraryMap          loadedLibraries;
    static DRLibrary::ServiceHash::Table        librariesByService;
    static DRLibrary::ServiceHash::StringTable  librariesByServiceName;
};

DRLibrary::Hash::LibraryMap DRLibrary::Cache::loadedLibraries;
DRLibrary::ServiceHash::Table DRLibrary::Cache::librariesByService;
DRLibrary::ServiceHash::StringTable DRLibrary::Cache::librariesByServiceName;

DRLibrary::~DRLibrary(){}

/** Loads the specified library */
DRLibrary::DRLibrary(const string& libraryName, 
                     const string& serviceName,
                     const string& serviceVersion):
    CObject(TYPE), libraryName(libraryName), serviceName(serviceName),
    serviceVersion(serviceVersion) {
    validatePop2Object();
}
/** Loads the specified DRI library */
DRLibrarySP DRLibrary::create(const string& libraryName, 
                              const string& serviceName,
                              const string& serviceVersion){
    return DRLibrarySP(new DRLibrary(libraryName, serviceName, serviceVersion));
}

DRLibrary::DRLibrary(DRService*    svc, 
                     const string& libraryName, 
                     const string& serviceName, 
                     const string& serviceVersion):
    CObject(TYPE), libraryName(libraryName), serviceName(serviceName),
    serviceVersion(serviceVersion) {
    handle = refCountPtr<LibHandle>(new LibHandle());
    handle->osHandle = 0; // means we didn't load it
    handle->svc = svc;
    // now be brave and load the types
    // XClass::loadAllClasses(handle->svc); - not just yet!
    // record in our various caches
    Cache::librariesByService[handle->svc] = DRLibrarySP(this);
    Cache::librariesByServiceName[serviceName] = DRLibrarySP(this);
    Cache::loadedLibraries[DRLibrarySP(this)] = handle;
}
    
/** Override default CObject implementation */
IObject* DRLibrary::clone() const{
    // think we want single reference to library so that we can unload it
    return const_cast<IObject*>((const IObject*)this);
}

/** Returns true if this DRLibrary specifies the same external DR
    Library as the supplied drLib */
bool DRLibrary::equals(const DRLibrarySP& drLib){
    if (this == drLib.get()){
        return true;
    }
    return (libraryName == drLib->libraryName &&
            serviceName == drLib->serviceName &&
            serviceVersion == drLib->serviceVersion);
}

/** Returns the name of the service (eg EDG). This should be unique
    per library and OS independent */
const string& DRLibrary::getServiceName() const{
    return serviceName;
}

typedef DRBool (DRCreateServiceFunc)(DRServiceInitArgs* args, 
                                     DRService**        service);

typedef DRBool (DRHandleToDRValue)(const char* handle,
                                   DRValue*    value);
#define HANDLE_TO_DRVALUE "HandleToDRValue"

//// returns true if the file exists
static bool fileExists(const string& file){
    FILE* filePtr = fopen(file.c_str(), "r");
    if (filePtr){
        fclose(filePtr);
        return true;
    }
    return false;
}

//// uses the CommandLineParams::DLLPath to try and find library
static string findDLL(const string& lib){
    string libName = lib;
#ifdef _MSC_VER
    const char delim = ';';
    // turn any .so into .dll
    if (libName.substr(lib.size()-3) == ".so"){
        libName.replace(lib.size()-3, 3, ".dll");
    }
#else
    const char delim = ':';
    // turn any .dll or .xll into .so
    if (libName.substr(lib.size()-4) == ".dll"){
        libName.replace(lib.size()-4, 4, ".so");
    } else if (libName.substr(lib.size()-4) == ".xll"){
        libName.replace(lib.size()-4, 4, ".so");
    }
#endif
    if (fileExists(lib)){
        return lib;
    }
    // strip path off
    size_t lastSlash = lib.find_last_of('\\');
    if (lastSlash == string::npos){
        // not found
        lastSlash = lib.find_last_of('/');
    }
    if (lastSlash == string::npos){
        libName = lib;
    } else {
        libName = string(lib, lastSlash+1);
    }
    IObjectSP searchPath(
        CommandLineParams::getParameterArgs(CommandLineParams::DLLPath));
    bool found = false;
    if (searchPath.get()){
        CStringSP str(CStringSP::dynamicCast(searchPath));
        const string& path = str->stringValue();
        string thePath;
        size_t currPos = 0;
        do {
            // pull of first part
            size_t delimPos = path.find_first_of(delim, currPos);
            if (delimPos == string::npos){
                currPos = path.size();
            }
            thePath = string(path, 0, delimPos-1);
            thePath += '/';
            thePath += libName;
            found = fileExists(thePath);
        } while (!found && currPos < path.size());
        if (found){
            return thePath;
        }
    }
    throw ModelException("DRLibrary::findDLL", "Unable to locate "+lib);
}

/** Called after an instance of a class (which implements this
    interface) is constructed via 'pop2Object' ie reflection is
    used to fill the internal fields of the object directly */
void DRLibrary::validatePop2Object(){
    static const string method("DRLibrary::validatePop2Object");
    static const char* CREATE_SERVICE = "DRCreateService";
    // change forward slashes into backslashes.
    // This is crucial to avoid loading multiple instances of the same xll
    string libName = libraryName;
    for (size_t i = 0; i < libName.size(); i++){
        if (libName[i] == '/'){
            libName[i] = '\\';
        }
    }
    if (!handle || !handle->svc){
        // is this in the cache
        Hash::LibraryMap::iterator iter = 
            Cache::loadedLibraries.find(DRLibrarySP(this));
        if (iter != Cache::loadedLibraries.end()){
            // copy over handle SP
            handle = iter->second;
        } else {
            handle = refCountPtr<LibHandle>(new LibHandle());
            libName = findDLL(libName); // search presupplied path
            // try and load library
#if defined(_MSC_VER)
            HINSTANCE dllHandle = LoadLibrary(libName.c_str());
#else
            void*     dllHandle = dlopen(libName.c_str(), RTLD_LAZY);
#endif
            if (!dllHandle) {
                throw ModelException(method, "Failed to load " + libName + "("+
#if defined(_MSC_VER)
                                     Format::toString((int) GetLastError())
#else
                                     string(dlerror())
#endif
                                     +")"
                    );
            }
            handle->osHandle = (void*) dllHandle; // save
            try{
                DRCreateServiceFunc* funcPtr = (DRCreateServiceFunc*)
#if defined(_MSC_VER)
                    GetProcAddress(dllHandle, CREATE_SERVICE);
#else
                dlsym(dllHandle, CREATE_SERVICE);
#endif
                if (!funcPtr) {
                    throw ModelException(
                        method, "Failed to lookup " +
                        string(CREATE_SERVICE) +
#if defined(_MSC_VER)
                        Format::toString((int)(GetLastError()))
#else
                        string(dlerror())
#endif
                        );
                }
                DRServiceInitArgs initArgs;
                initArgs.serviceName = serviceName.c_str();
                if (!serviceVersion.empty()){
                    throw ModelException(method, "Service Version argument not"
                                         " currently supported");
                }
                initArgs.serviceVersion = 0;
                initArgs.interfaceVersion = DRI_VERSION;
                initArgs.nOptions = 0;
                initArgs.options = 0;
                initArgs.ignoreUnrecognized = DR_TRUE;
                // start up service
                if (!funcPtr(&initArgs, &handle->svc) || !handle->svc){
                    throw ModelException(method, 
                                         "Call to CreateService failed for "+
                                         libraryName);
                }
                // now be brave and load the types
                // XClass::loadAllClasses(handle->svc); - not just yet!
            } catch (exception&){
                try{
                    handle->unloadLibrary();
                } catch (exception&){} //ignore
                throw;
            }
            Cache::librariesByService[handle->svc] = DRLibrarySP(this);
            Cache::librariesByServiceName[serviceName] = DRLibrarySP(this);
            Cache::loadedLibraries[DRLibrarySP(this)] = handle;
        }
    }
}

/** Returns the DRService pointer for this external DR Library */
DRService* DRLibrary::getService() {
    if (!handle || !handle->svc){
        throw ModelException("DRLibrary::getService",libraryName+"["+
                             serviceName+"] has been unloaded");
    }
    return handle->svc;
}

/** A record is kept of all libraries that have been loaded. This function
    returns the associated DRLibrary for the supplied service */
DRLibrarySP DRLibrary::getLibrary(DRService* svc){
    ServiceHash::Table::iterator iter =
        Cache::librariesByService.find(svc);
    if (iter == Cache::librariesByService.end()){
        throw ModelException("DRLibrary::getLibrary", 
                             svc? "No matching library":
                             "Service pointer is null");
    }
    return iter->second;
}

/** A record is kept of all libraries that have been loaded. This function
    returns the associated DRLibrary for the supplied service name */
DRLibrarySP DRLibrary::getLibrary(const string& svcName){
    ServiceHash::StringTable::iterator iter =
        Cache::librariesByServiceName.find(svcName);
    if (iter == Cache::librariesByServiceName.end()){
        throw ModelException("DRLibrary::getLibrary", 
                             "No library loaded for service with name "+
                             svcName);
    }
    return iter->second;
}
/** For within a spreadsheet environment, turns supplied handle into
    corresponding object */
IDRObjectSP DRLibrary::handleToObject(const string& handle) {
    static const string method("DRLibrary::handleToObject");
#if defined(_MSC_VER)
    // hack it for now. Long term have an object to do this. Object created
    // by factory using a user supplied string
    try{
        if (serviceName == "EDG"){
            // get the service
            DRService* svc = getService();
            // create the right action map
            XMapSP xMap(XMap::create("Handle::ToObject", svc));
            // pass in our handle name
            xMap->setItem("name",  CStringSP(CString::create(handle)));
            // turn into object
            XObjectSP xObj(xMap->toObject());
            // execute
            IDRObjectSP theObj(xObj->execute());
            return theObj;
        } else if (true /*serviceName == "IRLib"*/){
            HMODULE dllHandle = GetModuleHandle("DRLoaderAddin.xll");
            if (!dllHandle){
                throw ModelException(method,
                                     "Couldn't attach to DRLoaderAddin.xll");
            }
            DRHandleToDRValue* funcPtr = 
                (DRHandleToDRValue*)GetProcAddress(dllHandle,
                                                   HANDLE_TO_DRVALUE);
            if (!funcPtr){
                throw ModelException(method,
                                     "Couldn't find function "
                                     HANDLE_TO_DRVALUE "in DRLoaderAddin.xll");
            }
            DRValue drValue;
            if (!funcPtr(handle.c_str(), &drValue)){
                throw ModelException(method, "Couldn't retrieve handle "+
                                     handle);
            }
            /* It is possible that the service pointer in this DRLibrary is
               different to the svc pointer in objects returned
               from this function (depending on how the createService works).
               So if different we just add to our hash of librariesByService */
            DRService* svc = getService();
            bool isObj = drValue.type == DR_OBJECT;
            if (isObj && svc != drValue.value.object->svc){
                Cache::librariesByService[drValue.value.object->svc] = 
                    DRLibrarySP(this);
            }
            // we get a reference to objects - therefore we do not own
            // the memory - so pass false to toIDRObject()
            IDRObjectSP theObject(XObject::toIDRObject(drValue, false, svc));
            if (isObj){
                // again, must clone because we have a reference
                theObject = IDRObjectSP(theObject.clone());
            }
            return theObject;
        } else {
            throw ModelException(method, "Don't know how to retrieve handles"
                                 " for "+serviceName);
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
#else
    throw ModelException(method, "Only supported on Windows");
#endif
}

/** Registers the previously loaded DRI library */
DRLibrarySP DRLibrary::registerLib(DRService*    svc,
                                   const string& libraryName, 
                                   const string& serviceName,
                                   const string& serviceVersion){
    return DRLibrarySP(new DRLibrary(svc, libraryName,
                                     serviceName, serviceVersion));
}

/** Unloads the DRI library associated with this DRLibrary object. 
    Operations on this object which require the library to be loaded will
    subsequently fail  */
void DRLibrary::unload(){
    DRLibrarySP library(this);
    DRService* svc = handle->svc;
    handle->unloadLibrary();
    Cache::loadedLibraries.erase(library);
    Cache::librariesByService.erase(svc);
    // then delete any other matching entries in librariesByService
    Hash hashFnc;
    for (ServiceHash::Table::iterator iter = 
             Cache::librariesByService.begin(); 
         iter != Cache::librariesByService.end(); /* ++ in loop */){
        if (hashFnc(iter->second, library)){ // if equals
            Cache::librariesByService.erase(iter++);
        } else {
            ++iter;
        }
    }
    // similarly for librariesByServiceName
    { // work around MSVC
        for (ServiceHash::StringTable::iterator iter = 
                 Cache::librariesByServiceName.begin(); 
             iter != Cache::librariesByServiceName.end(); /* ++ in loop */){
            if (hashFnc(iter->second, library)){ // if equals
                Cache::librariesByServiceName.erase(iter++);
            } else {
                ++iter;
            }
        }
    }
}

DRLibrary::DRLibrary(): CObject(TYPE){}

IObject* DRLibrary::defaultConstructor(){
    return new DRLibrary();
}

void DRLibrary::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(DRLibrary, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(libraryName, "Name, including path, of DR Library");
    FIELD(serviceName, "Name of service to invoke");
    FIELD(serviceVersion, "Version to use");
    FIELD_MAKE_OPTIONAL(serviceVersion);
    Addin::registerConstructor("LOAD_LIBRARY",
                               Addin::UTILITIES,
                               "Loads a DR library",
                               TYPE);
}

CClassConstSP const DRLibrary::TYPE = 
CClass::registerClassLoadMethod("DRLibrary", typeid(DRLibrary), load);

class DRLibrary::Unload: public CObject{
    static CClassConstSP const TYPE;

    DRLibrarySP library;

    static bool unload(Unload* params){
        DRLibrarySP library(params->library);
        library->unload();
        return true;
    }

    Unload(): CObject(TYPE){}

    static IObject* defaultConstructor(){
        return new Unload();
    }
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(Unload, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(library, "The DR Library");
        Addin::registerInstanceBoolMethod("UNLOAD_LIBRARY", 
                                          Addin::UTILITIES,
                                          "Unloads a DR library",
                                          TYPE,
                                          (Addin::BoolMethod*)unload);
    }
};

CClassConstSP const DRLibrary::Unload::TYPE = 
CClass::registerClassLoadMethod("DRLibrary::Unload", typeid(Unload), load);
    

DRLIB_END_NAMESPACE

