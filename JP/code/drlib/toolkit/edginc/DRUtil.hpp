/*----------------------------------------------------------------------------

   Group       : Equity Derivatives Research

   Filename    : DRUtil.hpp

   Description : A collection of utilities to facilitate the use of the
                 C interface of DRI 1.2 (in applications written in C/C++)
                 The following are some cumbersome tasks made easier by this
                 utility:
                   - Error handling during DRService invocation 
                   - DRString (de-)allocation
                   - DRError handling after invocating a method in DRService

   Author      : Alex W. Fung

   Date        : Mar 24, 2006

----------------------------------------------------------------------------*/

#ifndef EDR_DRUTIL_HPP
#define EDR_DRUTIL_HPP
#include "edginc/DRAnalyticsInterface.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Error handling during service invocation
 */

/* Initialize default error handler. Not thread-safe. */
TOOLKIT_DLL void driServiceInitSetDefaultErrorHandler(DRServiceInitArgs *args);

/* Return default error structure written in support of EDG service.
   Do not free errorString. Not thread-safe. */
TOOLKIT_DLL DR_ERROR driGetServiceInvocationError(const char **errorString);

/*
 * DRValue utilities
 */
#define DRI_VALUE_CLEAR(v) \
    (v).type = DR_UNDEFINED; \
    (v).value.integer = 0; 

/*
 * DRString utilities
 */
TOOLKIT_DLL DRString driStringCreate(const char* in);

TOOLKIT_DLL void driStringFree(DRString str);

/*
 * Memory utilities
 */
#define DRI_FREE(svc, o)                                        \
    {                                                           \
        DRError err = 0;                                        \
        if (o) {                                                \
            if ((err = (svc)->fptr->objectFree((svc), (o)))) {  \
                driStringFree(err);                             \
            }                                                   \
        }                                                       \
    }

#define DRI_STRING_FREE(svc, s) { if (s) (svc)->fptr->stringFree((svc), (s)); }

#define DRI_ERROR_FREE(svc, e)  DRI_STRING_FREE((svc), (e))

/*
 * DRError utililities
 */ 

/* Client can always define another macro to fix usermsg and callback. */
#define DRI_CHECK(driFunc,                                              \
                  method,                                               \
                  usermsg,                                              \
                  callback,                                             \
                  cbParam,                                              \
                  onError,                                              \
                  postError,                                            \
                  onNoError)                                            \
    {                                                                   \
        DRError errmsg = (driFunc);                                     \
        if (errmsg) {                                                   \
            onError;                                                    \
            driHandleError((method), errmsg, (usermsg),                 \
                           callback, cbParam);                          \
            postError;                                                  \
        }                                                               \
        else {                                                          \
            onNoError;                                                  \
        }                                                               \
    }

/* Users can throw exceptions or log errors. */
typedef void (*DRI_ERROR_CALLBACK)(const char *method, 
                                   const char *errMsg,
                                   void       *callbackParam);

/* Compose the full error message and invoke callback. */
TOOLKIT_DLL void driHandleError(const char         *method,
                                DRError             errmsg,
                                const char         *userMsg,
                                DRI_ERROR_CALLBACK  callback,
                                void               *callbackParam);

TOOLKIT_DLL DRError driErrorCreate(const char* in);

TOOLKIT_DLL DRError driOutOfMemoryError();

#ifdef __cplusplus
}
#endif

#endif
