/*  ---------------------------------------------------------------------------
 *  Equity Derivatives Research: EDGServices.h
 *
 *  Functions to supplement the EDG implementation of the DR Interface
 *  Focus around market data cache/DR Wrapper/results handling
 *
*/


#ifndef EDG_SERVICES_H
#define EDG_SERVICES_H

#include "edginc/DRAnalyticsInterface.h"

#ifdef __cplusplus
extern "C" {
#endif

    /************* DRValue Support *********************/
    /** Populate a DRValue e.g. EDGSDRValue(&value, DR_INT, 7) */
    DRI_DLL_EXPORT DRError EDGSDRValue(DRValue* drValue, int type, ...);

    /************* Supplementary object support  *********************/
    /** Create a copy of the supplied object - plugs gap prior to DRI v2 */
    DRI_DLL_EXPORT DRObject EDGSObjectClone(DRService* service,
                                            DRObject   object);

    /************* Supplementary array support  *********************/

    /** Appends the supplied object to the end of supplied array. The
        length of the array is increased by 1. (aka push_back).
        The supplied object is not deep copied when it is appended.
        Returns DR_FALSE on failure*/
    DRI_DLL_EXPORT DRError EDGSArrayAppend(DRService *service,
                                           DRArray    array,
                                           DRValue   *drValue);

    /** what is the type of element in an array ?
        do NOT free the output of this function */
    DRI_DLL_EXPORT const char* EDGSArrayElemType(DRService*  service,
                                                 const char* arrayType);


    /************* DR Wrapper/Regression Test Creation support **************/

    /** return wrapper errors from file as a dynamic C style string -
        use stringFree() to cleanup */
    DRI_DLL_EXPORT DRString EDGSDRWrapperError(DRService*  service,
                                           const char* filename);

    /** read wrapper results back from file. */
    DRI_DLL_EXPORT DRObject EDGSDRWrapperRead(DRService*  service,
                                          const char* filename);

    /** read a DRValue from a file */
    DRI_DLL_EXPORT DRError EDGSValueRead(DRService  *service,
                                         const char *fileName,
                                         DRValue    *object);

    /** write a DRValue to a file */
    DRI_DLL_EXPORT DRError EDGSValueWrite(DRService  *service,
                                          const char *fileName,
                                          DRValue    *object);

    /** read a DRValue from a buffer (XML string) */
    DRI_DLL_EXPORT DRError EDGSValueXMLRead(DRService  *service,
                                            const char *buffer,
                                            DRValue    *object);

    /** write a DRValue to an XML buffer */
    DRI_DLL_EXPORT DRString EDGSValueXMLWrite(DRService* service, DRValue* object);

    /** Writes the supplied results object to file (includes a 'summary' of the
        results at the top of the file) */
    DRI_DLL_EXPORT DRError EDGSWriteResultsToFile(DRService  *service,
                                                  const char *filename,
                                                  DRObject    results);

    /** Dynamically loads the specified DRI library and returns a
     * service pointer */
    DRI_DLL_EXPORT DRService* EDGSLoadLibrary(const char* libName,
                                              const char* serviceName,
                                              const char* serviceVersion);

    /** Unloads the specified DRI library */
    DRI_DLL_EXPORT DRError EDGSUnloadLibrary(DRService* drService);

    /** Registers the previously loaded DRI library with this library */
    DRI_DLL_EXPORT DRError EDGSRegisterLibrary(DRService  *drService,
                                               const char *libName,
                                               const char *serviceName,
                                               const char *serviceVersion);

    // Initialize DRCreateInitArgs for invoking EDG service.
    DRI_DLL_EXPORT void EDGSInitServiceCreateArgs(DRServiceInitArgs *args);

    // Return default error structure written in support of EDG service.
    DRI_DLL_EXPORT DR_ERROR EDGSGetServiceInvocationError(
        const char **errorString);

}

#endif
