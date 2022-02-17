/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Set of functions for ARM_CorbaRequestModule::ARM_CORBA_PARAM               */
/*                                                                            */
/*----------------------------------------------------------------------------*/
 

#ifndef ARM_CORBA_PARAM_H 
#define ARM_CORBA_PARAM_H


#include "req.h"


extern void PARAM_Print(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                        long idx = 0); 

 
extern void PARAM_SetLong(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                          long lVal);
extern CORBA::Long PARAM_GetLong(
                 const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param); 

 
extern void PARAM_SetDouble(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                            double dVal);
extern CORBA::Double PARAM_GetDouble(
                  const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param); 

 
extern void PARAM_SetString(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, 
                            const char* cVal);
extern char* PARAM_GetString(
                      const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param); 

 
extern void PARAM_SetLongArray(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                               long* arVals, long nbVal);
extern void PARAM_GetLongArray(
                 const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                               long** arVals, long* nbVal); 

 
extern void PARAM_GetDoubleArray(
                       const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                                 double** arVals, long* nbVal);
extern void PARAM_SetDoubleArray(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                                 double* arVals, long nbVal);

 
extern void PARAM_GetStringArray(
                      ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                                 char**& arVals, long* nbVal);
extern void PARAM_SetStringArray(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                                 char** arVals, long nbVal); 



extern CORBA::Long PARAM_GetLongFromArray(
           const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long index);

extern CORBA::Double PARAM_GetDoubleFromArray(
           const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long index);

extern void PARAM_GetStringFromArray(
   ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long index, char* value);
 
extern CORBA::Long PARAM_GetNbValues(
           const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param);

extern void PARAM_SetNbValues(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                              long nbVals);


extern ARM_CorbaRequestModule::ARM_CORBA_VALUE_TYPE PARAM_GetValueType(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param);

extern void PARAM_SetValueType(
            ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, 
            ARM_CorbaRequestModule::ARM_CORBA_VALUE_TYPE vType); 


extern long PARAM_GetNthLongValue(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long rank);

extern double PARAM_GetNthDoubleValue(
   const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long rank);

extern char* PARAM_GetNthStringValue(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                                     long rank);



#endif 	// ARM_CORBA_PARAM_H
/*----------------------------------------------------------------------------*/
/*---- End of File ----*/

// EOF %M%
