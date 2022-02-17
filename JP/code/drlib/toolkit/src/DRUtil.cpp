/*----------------------------------------------------------------------------

   Group       : Equity Derivatives Research

   Filename    : DRUtil.c

   Description : A collection of utilities to facilitate the use of the
                 C interface of DRI 1.2 (in applications written in C/C++)

   Author      : Alex W. Fung

   Date        : Mar 24, 2006

----------------------------------------------------------------------------*/

#include "edginc/config.hpp"
#include "edginc/DRUtil.hpp"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * Error handling during service invocation
 */

/* Example errorHandler and createServiceArgs callbacks. */
typedef struct ServiceInvocationError_
{
    const char *errorString;
    DR_ERROR    errorCode;
} ServiceInvocationError;

/* Default error holder (to be used by localErrorHandler) */
static ServiceInvocationError s_errorHandler;

/* Default error handler written in support of EDG service */
static void localErrorHandler(void       *userData,
                              const char *errorString,
                              DR_ERROR    errorCode)
{
    if (userData) {
        ServiceInvocationError *err = (ServiceInvocationError*) userData;
        err->errorString = errorString;
        err->errorCode   = errorCode;
    }
    else {
        s_errorHandler.errorString = errorString;
        s_errorHandler.errorCode   = errorCode;
    }
}

/* Initialize error handler. */
void driServiceInitSetDefaultErrorHandler(DRServiceInitArgs *args)
{
    if (!args) {
        return;
    }

    args->errorHandler = localErrorHandler;
}

/* Return default error structure written in support of EDG service.
   Do not free errorString. Not thread-safe. */
DR_ERROR driGetServiceInvocationError(const char **errorString)
{
    *errorString = s_errorHandler.errorString;
    return s_errorHandler.errorCode;
}

/*
 * DRString/DRError handling with provisions for out-of-memory errors
 */
static const char OUT_OF_MEMORY[] = "Out of memory when creating DRError";

DRString driStringCreate(const char* in)
{
    char *out;
    if ((out = (char*) malloc(sizeof(char) * strlen(in) + 1))) {
        strcpy(out, in);
    }
    return out;
}

DRError driErrorCreate(const char* in)
{
    char *out;
    if ((out = (char*) malloc(sizeof(char) * strlen(in) + 1))) {
        strcpy(out, in);
    }
    return OUT_OF_MEMORY;
}

DRError driOutOfMemoryError()
{
    return OUT_OF_MEMORY;
}

void driStringFree(DRString str)
{
    /* Assume that the caller has provisions to handle out-of-memory errors
       (e.g., by freeing unused memory) and want to proceed and free the
       DRError collected. */
    if (str && str != OUT_OF_MEMORY) {
        free((void*) str);
    }
}

/*
 * DRError handling (without printing line number for now)
 */ 
void driHandleError(const char         *method,
                    DRError             errMsg,
                    const char         *userMsg,
                    DRI_ERROR_CALLBACK  callback,
                    void               *callbackParam)
{
    enum { 
        METHOD_LEN    = 32, 
        USR_MSG_LEN   = 128,
        ERR_MSG_LEN   = 128,
        MSG_LEN       = USR_MSG_LEN + ERR_MSG_LEN + 3
    };

    char wholeMsg[MSG_LEN + METHOD_LEN + 1];

    if (userMsg && *userMsg) {
        /* snprintf is not enabled. */
        strncpy(wholeMsg, userMsg, USR_MSG_LEN);
        strcat(wholeMsg, " - ");
        strncat(wholeMsg, errMsg, ERR_MSG_LEN);
    }
    else {
        strncpy(wholeMsg, errMsg, ERR_MSG_LEN);
    }

    driStringFree((DRString) errMsg);

    if (callback) {
        callback(method, wholeMsg, callbackParam);
    }
    else {
        fprintf(stderr, "%s: %s\n", method, wholeMsg);
    }
}

