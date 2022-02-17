/*  ---------------------------------------------------------------------------
 *  Equity Derivatives Research: EDGServices.cpp
 *
 *  Functions to supplement the EDG implementation of the DR Interface
 *
 *
*/

#include "edginc/config.hpp"
#include "edginc/EDGServices.h"
#include "edginc/DRLibrary.hpp"
#include "edginc/DRAnalyticsInterface.h"
#include "edginc/Format.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/Library.hpp"
#include "edginc/Version.hpp"
#include "edginc/DRUtil.hpp"

#ifdef _MSC_VER
/* only for NT of course */
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

USING_DRLIB_NAMESPACE

//----------------------------------------------------------------------------
// These functions implement the additional EDG services
//----------------------------------------------------------------------------

//// Allow EDRInterface access to these functions. A bit hacky but saves
//// a lot of work and (in theory....) the need for these will go soon
//// convert an IObject to a DRValue
//// They are defined in DRInterface.cpp
extern void object2DRValue(IObject* object, DRValue* drValue, DRService *svc);
extern IObjectSP drValue2Object(const DRValue* drValue);
extern IObject* drObject2Object(DRObject drObj);
extern DRObject object2DRObject(IObject* object, DRService *svc);
extern void isEdrObject(DRObject object);

static char* stringCopy(const char* in) {
    char *out = (char*)malloc(sizeof(char) * strlen(in) + 1);
    strcpy(out, in);
    return out;
}

// Return the stack trace of a ModelException with out-of-memory provisions. 
static DRError EDGSExceptionStackTrace(const ModelException& exception){
    char *err = exception.stackTrace();
    return err ? err : driOutOfMemoryError();
}

/** Populate a DRValue e.g. EDGSDRValue(&value, DR_INT, 7) */
DLL_EXPORT DRError EDGSDRValue(DRValue* drValue, int type, ...) {
    static const char method[] = "EDGSDrValue";
    try {
        va_list args;
        va_start(args, type);
        if (!drValue) {
            throw ModelException(method, "DRValue is null");
        }
        switch (type)
        {
        case DR_BOOL:
            drValue->value.boolean = va_arg(args, DRBool);
            break;
        case DR_INT:
            drValue->value.integer = va_arg(args, int);
            break;
        case DR_DOUBLE:
            drValue->value.real = va_arg(args, double);
            break;
        case DR_STRING:
            if (!(drValue->value.string = stringCopy(va_arg(args, char*)))) {
                throw ModelException(method, "string copy failed");
            }
            break;
        case DR_OBJECT:
        {
            drValue->value.object = va_arg(args, DRObject);
            isEdrObject(drValue->value.object);
            break;
        }
        case DR_DATE:
            throw ModelException(method,
                                 "EDGService does not support DR_DATE type");
        case DR_UNDEFINED:
            break;  // undefined is OK - maybe it's an initialisation
        default:
        {
            throw ModelException(method,
                                 "Unrecognised DR_VALUE type (" +
                                 Format::toString(type) + ")");

        }
        }
        // set that type
        drValue->type = type;

        va_end(args);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

/** Create a copy of the supplied object - plugs gap prior to DRI v2 */
DLL_EXPORT DRObject EDGSObjectClone(DRService* service,
                                    DRObject   object) {
    try {
        IObjectSP obj(drObject2Object(object));
        IObjectSP copy(obj->clone());
        return object2DRObject(copy.get(), service);
    }
    catch (exception&){
        return 0;
    }
}

/** Appends the supplied object to the end of supplied array. The
    length of the array is increased by 1. (aka push_back).
    The supplied object is not deep copied when it is appended.
    Returns DR_FALSE on failure*/
DLL_EXPORT DRError EDGSArrayAppend(DRService* service,
                                   DRArray    drarray,
                                   DRValue*   drValue){
    try {
        IObjectSP object(drObject2Object(drarray));
        IArraySP array(IArraySP::dynamicCast(object));

        IObjectSP elem(drValue2Object(drValue));
        array->append(elem);
        return 0;
    }
    catch (exception& e) {
        return EDGSExceptionStackTrace(e);
    }
}

/** what is the type of element in an array ?
    do NOT free the output of this function */
DLL_EXPORT const char* EDGSArrayElemType(DRService* service,
                                         const char* arrayType) {
    const static char method[] = "EDGSArrayElemType";
    try {
        if (!arrayType){
            throw ModelException(method, "NULL type supplied");
        }
        CClassConstSP clazz = CClass::forName(string(arrayType));
        CClassConstSP arrayClazz = clazz->getComponentType();

        return (arrayClazz->getName().c_str());
    }
    catch (exception&) {
        return 0;
    }
}

/** return wrapper errors from file as a dynamic C style string -
    use stringFree() to cleanup */
DLL_EXPORT DRString EDGSDRWrapperError(DRService* service, const char* filename) {
    try {
        string err(OutputFile::xmlReadErrors(filename));
        return (stringCopy(err.c_str()));
    }
    catch (exception& /*e*/){
        return 0;
    }
}

// read wrapper results back from file
DLL_EXPORT DRObject EDGSDRWrapperRead(DRService* service, const char* filename) {
    try {
        IObjectSP o(OutputFile::xmlReadResults(filename));
        return object2DRObject(o.release(), service);
    }
    catch (exception& /*e*/){
        return 0;
    }
}

// read a DRValue from a file
DLL_EXPORT DRError EDGSValueRead(DRService*  service,
                                 const char* fileName,
                                 DRValue*    object) {
    try {
        XMLReader reader(fileName, true);
        IObjectSP o(reader.read());
        object2DRValue(o.get(), object, service);
        return 0;
    }
    catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}

// write a DRValue to a file
DLL_EXPORT DRError EDGSValueWrite(DRService*  service,
                                  const char* fileName,
                                  DRValue*    object) {
    try {
        IObjectSP o(drValue2Object(object));
        XMLWriter xml(fileName);
        o->write("OBJECT", &xml);
        return 0;
    }
    catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}

// read a DRValue from an XML buffer
DLL_EXPORT DRError EDGSValueXMLRead(DRService*  service,
                                    const char* buffer,
                                    DRValue*    object) {
    static const char method[] = "EDGSValueXMLRead";
    try {
        if (!buffer){
            throw ModelException(method, "buffer is null");
        }
        // convert char* into string
        string xmlin(buffer);
        XMLReader reader(xmlin, false);
        IObjectSP o(reader.read());
        object2DRValue(o.get(), object, service);
        return 0;
    }
    catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}

// write a DRValue to an XML buffer
DLL_EXPORT DRString EDGSValueXMLWrite(DRService* service, DRValue* object) {
    try {
        IObjectSP o(drValue2Object(object));
        // create xml output stream
        string buffer;
        XMLWriter xml(buffer, false);
        o->write("OBJECT", &xml);
        return stringCopy(buffer.c_str());
    }
    catch (exception& /*e*/){
        return 0;
    }
}

/** Writes the supplied results object to file (includes a 'summary' of the
    results at the top of the file) */
DLL_EXPORT DRError EDGSWriteResultsToFile(DRService*  service,
                                          const char* filename,
                                          DRObject    results) {
    static const char method[] = "EDGSWriteResultsToFile";
    try {
        if (!filename){
            throw ModelException(method, "filename is null");
        }
        IObject* o = drObject2Object(results);
        OutputFile  output(filename);
        output.write(o);
        return 0;
    }
    catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}
/** Dynamically loads the specified DRI library and returns a service
    pointer */
DLL_EXPORT DRService* EDGSLoadLibrary(const char* libName,
                                      const char* serviceName,
                                      const char* serviceVersion){
    static const char method[] = "EDGSLoadLibrary";
    try{
        if (!libName){
            throw ModelException(method, "libName is null");
        }
        DRLibrarySP drLib(
            DRLibrary::create(libName,
                              !serviceName? "": serviceName,
                              !serviceVersion? "": serviceVersion));
        return drLib->getService();
    } catch (exception& /*e*/){
        return 0;
    }
}

/** Unloads the specified DRI library */
DLL_EXPORT DRError EDGSUnloadLibrary(DRService* drService){
    try{
        DRLibrarySP drLib(DRLibrary::getLibrary(drService));
        drLib->unload();
        drLib.reset();
        return 0;
    } catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}

/** Registers the previously loaded DRI library with this library */
DLL_EXPORT DRError EDGSRegisterLibrary(DRService* drService,
                                       const char* libName,
                                       const char* serviceName,
                                       const char* serviceVersion){
    static const char method[] = "EDGSRegisterLibrary";
    try{
        if (!libName){
            throw ModelException(method, "libName is null");
        }
        if (!drService){
            throw ModelException(method, "drService is null");
        }
        DRLibrary::registerLib(drService, libName,
                               !serviceName? "": serviceName,
                               !serviceVersion? "": serviceVersion);
        return 0;
    } catch (exception& e){
        return EDGSExceptionStackTrace(e);
    }
}

// Initialize DRCreateInitArgs for invoking EDG service.
DLL_EXPORT void EDGSInitServiceCreateArgs(DRServiceInitArgs *args)
{
    if (!args) {
        return;
    }

    args->serviceName         =  Library::SERVICE_NAME.c_str();
    args->serviceVersion      =  DRLIB_VERSION( 1, 0, 0, 0 );
    args->interfaceVersion    =  DRI_VERSION;
    args->nOptions            =  0;
    args->options             =  0;
    args->ignoreUnrecognized  =  DR_TRUE;
    args->userData            =  0;

    driServiceInitSetDefaultErrorHandler(args);
}

// Return default error structure written in support of EDG service.
// Do not free errorString.
DLL_EXPORT DR_ERROR EDGSGetServiceInvocationError(const char **errorString)
{
    return driGetServiceInvocationError(errorString);
}

DRLIB_BEGIN_NAMESPACE
// symbol (referenced by AddinLib.cpp) to ensure file gets linked in
bool EDGServicesLink = true;
DRLIB_END_NAMESPACE
