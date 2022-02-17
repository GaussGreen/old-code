/*
 * $Log: ARM_CORBA_REQUEST.cpp,v $
 * Revision 1.3  2000/04/05 14:36:39  mab
 * Version initiale
 *
 * Revision 1.2  2000/04/04 12:17:36  mab
 * version initiale
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Set of functions headers setting ARM_CORBA_REQUEST parameters              */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#include <stdio.h>


#include "req.h"
#include "ARM_CORBA_REQUEST.h" 






ARM_CorbaRequestModule::ARM_CORBA_REQUEST* REQ_NewInitializedRequest(long reqId,
                                                                     long nbPar)
{
    ARM_CorbaRequestModule::ARM_CORBA_REQUEST* newReq 
                  = new ARM_CorbaRequestModule::ARM_CORBA_REQUEST();

/* OLD CODE
    ARM_CorbaRequestModule::ARM_CORBA_PARAM* data 
 = ARM_CorbaRequestModule::ARM_CORBA_REQUEST::paramsList::allocbuf(nbPar); 

    REQ_Set(newReq, reqId, nbPar, data);
*/

    IT_Cxx_USeq< ARM_CorbaRequestModule::ARM_CORBA_PARAM > tmpParamList;

    ARM_CorbaRequestModule::ARM_CORBA_PARAM* data;
    data = tmpParamList.allocbuf(nbPar);

    IT_Cxx_USeq< ARM_CorbaRequestModule::ARM_CORBA_PARAM > inParamList(nbPar,
                                                           nbPar, data, 1);

    newReq->paramsList = inParamList;

    newReq->IsRequest  = ARM_CorbaRequestModule::OBJ_REQUEST;
    newReq->reqId      = reqId;                             
    newReq->nbParams   = nbPar;                             

    return(newReq);
}


/* OLD CODE
void REQ_Set(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, long reqId,
             long nbPar,
             ARM_CorbaRequestModule::ARM_CORBA_PARAM* data)
{
    req->IsRequest  = ARM_CorbaRequestModule::OBJ_REQUEST;
    req->reqId      = reqId;                                    
    req->nbParams   = nbPar;                                    

    req->paramsList.assignbuf(nbPar, nbPar, data, 1);
}
*/


void REQ_Print(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req)
{
    long i;

    fprintf(stderr, "\n===============> REQ <===============\n");

    fprintf(stderr, "\n ReqID = %ld\n", req->reqId);
    fprintf(stderr, "\n NbPar = %ld\n", req->nbParams);

    for (i = 0; i < req->nbParams; i++)
    {
        PARAM_Print(&req->paramsList[i], i);
    }
 
    fprintf(stderr, "\n<=============== REQ ===============>\n");
}


void REQ_SetLong(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                 long val, long rank)
{
    PARAM_SetLong(&req->paramsList[rank], val);
}


long REQ_GetLong(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                 long rank)
{
    return(PARAM_GetLong(&req->paramsList[rank]));
}


void REQ_SetDouble(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                   double val, long rank)
{
    PARAM_SetDouble(&req->paramsList[rank], val);
}


double REQ_GetDouble(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                     long rank)
{
    return(PARAM_GetDouble(&req->paramsList[rank]));
}


void REQ_SetString(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                   const char* val, long rank)
{
    PARAM_SetString(&req->paramsList[rank], val);
}


char* REQ_GetString(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                    long rank)
{
    return(PARAM_GetString(&req->paramsList[rank]));
} 


void REQ_SetLongArray(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                      long* arVals, long nbVal, long rank)
{
    PARAM_SetLongArray(&req->paramsList[rank], arVals, nbVal); 
}


void REQ_GetLongArray(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                      long** arVals, long* nbVal, long rank)
{
    PARAM_GetLongArray(&req->paramsList[rank], arVals, nbVal); 
}


void REQ_SetDoubleArray(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                        double* arVals, long nbVal, long rank)
{
    PARAM_SetDoubleArray(&req->paramsList[rank], arVals, nbVal); 
}


void REQ_GetDoubleArray(const ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                        double** arVals, long* nbVal, long rank)
{
    PARAM_GetDoubleArray(&req->paramsList[rank], arVals, nbVal); 
}


void REQ_SetStringArray(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                        char** arVals, long nbVal, long rank)
{
    PARAM_SetStringArray(&req->paramsList[rank], arVals, nbVal); 
}


void REQ_GetStringArray(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                        char**& arVals, long* nbVal, long rank)
{
    PARAM_GetStringArray(&req->paramsList[rank], arVals, nbVal);
}


ARM_CorbaRequestModule::ARM_CORBA_PARAM* 
REQ_GetNthParam(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req, long rank)
{
    return(&req->paramsList[rank]);
}


void REQ_SetNthParam(ARM_CorbaRequestModule::ARM_CORBA_REQUEST* req,
                   ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long rank)
{
    req->paramsList[rank] = *param;
}


/*----------------------------------------------------------------------------*/
/*---- End of File ----*/
