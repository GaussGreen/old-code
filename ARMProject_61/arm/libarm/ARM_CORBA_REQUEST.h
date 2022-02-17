/*
 * $Log: ARM_CORBA_REQUEST.h,v $
 * Revision 1.1  2000/04/04 12:07:51  mab
 * Initial revision
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Set of functions headers setting ARM_CORBA_REQUEST parameters              */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#ifndef _ARM_CORBA_REQUEST_H
#define _ARM_CORBA_REQUEST_H


#include "ARM_CORBA_PARAM.h" 


extern ARM_CorbaRequestModule::ARM_CORBA_REQUEST* REQ_NewInitializedRequest(
                                                     long reqId, long nbPar);

/* OLD CODE
extern void REQ_Set(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                    long reqId,
                    long nbPar, ARM_CorbaRequestModule::ARM_CORBA_PARAM* data);
*/


extern void REQ_Print(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req);

extern void REQ_SetLong(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                        long val, long rank);

extern long REQ_GetLong(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                        long rank);

extern void REQ_SetDouble(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                          double val, long rank);

extern double REQ_GetDouble(
           const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, long rank);

extern void REQ_SetString(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                          const char* val, long rank);

extern char* REQ_GetString(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                           long rank);

extern void REQ_SetLongArray(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                             long* arVals, long nbVal, long rank);

extern void REQ_GetLongArray(
                        const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                        long** arVals, long* nbVal, long rank);

extern void REQ_SetDoubleArray(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                               double* arVals, long nbVal, long rank);

extern void REQ_GetDoubleArray(
       const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, double** arVals,
       long* nbVal, long rank);

extern void REQ_SetStringArray(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                        char** arVals, long nbVal, long rank);

extern void REQ_GetStringArray(
                           ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                               char**& arVals, long* nbVal, long rank);

extern ARM_CorbaRequestModule::ARM_CORBA_PARAM* REQ_GetNthParam(
               ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, long rank); 

extern void REQ_SetNthParam(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                 ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                 long rank);



#endif 	// _ARM_CORBA_REQUEST_H
/*----------------------------------------------------------------------------*/
/*---- End of File ----*/


// EOF %M%
