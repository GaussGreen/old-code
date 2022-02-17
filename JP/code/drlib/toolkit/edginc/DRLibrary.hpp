//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : DRLibrary.hpp
//
//   Description : Class representing a reference to an external DR Library
//
//   Author      : Mark A Robson
//
//   Date        : 27 Nov 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_DRLIBRARY_HPP
#define EDR_DRLIBRARY_HPP
#include "edginc/Object.hpp"
#include "edginc/DRObject.hpp"
typedef struct DRService_ DRService;

DRLIB_BEGIN_NAMESPACE
class DRLibrary;
typedef smartPtr<DRLibrary> DRLibrarySP;
#ifndef QLIB_DRLIBRARY_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DRLibrary>);
#endif

/** Class representing a reference to an external DR Library */
class TOOLKIT_DLL DRLibrary: public CObject{
public:
    static CClassConstSP const TYPE; // defined in Perturbation.cpp

    virtual ~DRLibrary();

    /** Loads the specified DRI library */
    static DRLibrarySP create(const string& libraryName, 
                              const string& serviceName,
                              const string& serviceVersion);

    /** registers the already loaded DRI library */
    static DRLibrarySP registerLib(DRService*    svc,
                                   const string& libraryName, 
                                   const string& serviceName,
                                   const string& serviceVersion);

    /** Override default CObject implementation */
    virtual IObject* clone() const;

    /** Called after an instance of a class (which implements this
        interface) is constructed via 'pop2Object' ie reflection is
        used to fill the internal fields of the object directly */
    virtual void validatePop2Object();

    /** Returns the DRService pointer for this external DR Library */
    DRService* getService();

    /** Returns the name of the service (eg EDG). This should be unique
        per library and OS independent */
    const string& getServiceName() const;

    /** For within a spreadsheet environment, turns supplied handle into
        corresponding object */
    IDRObjectSP handleToObject(const string& handle);

    /** Returns true if this DRLibrary specifies the same external DR
        Library as the supplied drLib */
    bool equals(const DRLibrarySP& drLib);

    /** A record is kept of all libraries that have been loaded. This function
        returns the associated DRLibrary for the supplied service */
    static DRLibrarySP getLibrary(DRService* svc);

    /** A record is kept of all libraries that have been loaded. This function
        returns the associated DRLibrary for the supplied service name */
    static DRLibrarySP getLibrary(const string& svcName);

    
    /** Unloads the DRI library associated with this DRLibrary object. 
        Operations on this object which require the library to be loaded will
        subsequently fail  */
    void unload();

    struct Hash;
    struct ServiceHash;
    struct LibHandle;
private:
    class Unload;
    struct Cache;
    friend class Unload;
    friend struct Hash;
    DRLibrary(const DRLibrary& rhs); // don't use
    DRLibrary& operator=(const DRLibrary& rhs); // don't use
    DRLibrary();
    DRLibrary(const string& libraryName, const string& serviceName,
              const string& serviceVersion);
    DRLibrary(DRService*    svc, const string& libraryName, 
              const string& serviceName, const string& serviceVersion);

    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);

    //// fields
    const string     libraryName; // includes path
    const string     serviceName; // eg fxCom
    const string     serviceVersion; // optional
    refCountPtr<LibHandle> handle; // not registered $unregistered
};
DRLIB_END_NAMESPACE
#endif
