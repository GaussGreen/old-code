/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Set of functions setting ARM_CORBA_PARAM parameters                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <string.h>



#include "ARM_CORBA_PARAM.h"



void PARAM_Print(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long idx)
{
    fprintf(stderr, "\n====> ARM_CORBA_PARAM : %ld \n", idx);

    fprintf(stderr, "\n nbValues = %d", PARAM_GetNbValues(param));

    CORBA::ULong i = 0;

    switch(PARAM_GetValueType(param))
    {
        case ARM_CorbaRequestModule::ARM_CORBA_DOUBLE:
        {
            fprintf(stderr, "\n val = %lf\n", param->val.doubleVal());
        };
        break;

        case ARM_CorbaRequestModule::ARM_CORBA_LONG:
        {
            fprintf(stderr, "\n val = %ld\n", param->val.longVal());
        };
        break;

        case ARM_CorbaRequestModule::ARM_CORBA_STRING:
        {
            fprintf(stderr, "\n val = %s\n", param->val.stringVal ());
        };
        break;

        case ARM_CorbaRequestModule::ARM_CORBA_LONG_ARRAY:
        {
            for (i = 0; i < param->val.longArray().length(); i++)
            {
                fprintf(stderr, "\n val[%d] = %ld\n", i, 
                                param->val.longArray()[i]);
            }
        };
        break;

        case ARM_CorbaRequestModule::ARM_CORBA_DOUBLE_ARRAY:
        {
            for (i = 0; i < param->val.doubleArray().length (); i++)
            {
                fprintf(stderr, "\n val[%d] = %lf\n", i, 
                                param->val.doubleArray()[i]);
            }
        };
        break;

        case ARM_CorbaRequestModule::ARM_CORBA_STRING_ARRAY:
        {
            for (i = 0; i < param->val.stringArray().length(); i++)
            {
                fprintf(stderr, " val[%d] = %s\n", i, 
                        (const char*)(param->val.stringArray()[i]));
            }
        };
        break;

        default:
            break;
    }

    fprintf(stderr, "\n<==== ARM_CORBA_PARAM\n");
}


void PARAM_SetLong(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long lVal)
{
    param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

    param->nbValues = 1;

    param->val.longVal((CORBA::Long) lVal);
}


CORBA::Long PARAM_GetLong(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param)
{
    return(param->val.longVal());
}


void PARAM_SetDouble(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                     double dVal)
{ 
    param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

    param->nbValues = 1;

    param->val.doubleVal((CORBA::Double) dVal);
} 


CORBA::Double PARAM_GetDouble(
                     const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param)
{
    return(param->val.doubleVal());
}


void PARAM_SetString(
          ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, const char* cVal)  
{  
    param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

    param->nbValues = 1;

    param->val.stringVal(cVal);
}


char* PARAM_GetString(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param)
{ 
    return((char *) param->val.stringVal());
}


long PARAM_GetNthLongValue(
            const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, 
            long rank)
{
    return((param->val.longArray())[rank]);
}


double PARAM_GetNthDoubleValue(
            const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long rank)
{
    return((param->val.doubleArray())[rank]);
}


char* PARAM_GetNthStringValue(
        ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long rank)
{
    return((char *) (param->val.stringArray())[rank]);
}


void PARAM_GetLongArray(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                        long** arVals, long* nbVal)
{
    long* values;
    long  nb;
    int   i; 
    
    nb = param->val.longArray().length();

    values = new long [nb];

    for (i = 0; i < nb; i++)
    {
        values[i] = (param->val.longArray())[i];
    }

    *nbVal = nb;

    *arVals = values;
}


void PARAM_SetLongArray(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                        long* arVals, long nbVal)
{  
    param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;;

    param->nbValues = nbVal;

    ARM_CorbaRequestModule::ARM_CORBA_LONG_ARRAY_TYPE arrayLong(nbVal, 
                                                         nbVal, arVals, 0);

    param->val.longArray(arrayLong);
}


void PARAM_GetDoubleArray(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, 
                          double** arVals, long* nbVal)
{
    double* values;
    long    nb;
    long    i; 
 
    nb = param->val.doubleArray().length();
 
    values = new double [nb];
 
    for (i = 0; i < nb; i++)
    {
        values[i] = (param->val.doubleArray())[i];
    }
 
    *nbVal = nb;
 
    *arVals = values;
}


void PARAM_SetDoubleArray(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                          double* arVals, long nbVal)  
{
    param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

    param->nbValues = nbVal;

    ARM_CorbaRequestModule::ARM_CORBA_DOUBLE_ARRAY_TYPE arrayDouble(nbVal, 
                                                 nbVal, arVals, 0);
 
    param->val.doubleArray(arrayDouble);
}


void PARAM_GetStringArray(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, 
                          char**& arVals, long* nbVal)
{
    char** values;
    long   nb;
    int    i;
 
    nb = param->val.stringArray().length();
 
    values = new char*[nb];
 
    for (i = 0; i < nb; i++)
    {
        char* curVal;
        long  strLen;

        curVal = (char *) (param->val.stringArray())[i];
        strLen = strlen(curVal);

        values[i] = new char[strLen+1];
        strcpy(values[i], curVal);
    }
 
    *nbVal = nb;
 
    arVals = values;
}


void PARAM_SetStringArray(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                          char** arVals, long nbVal)  
{
    param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

    param->nbValues = nbVal;
 
    ARM_CorbaRequestModule::ARM_CORBA_STRING_ARRAY_TYPE  as(nbVal, nbVal, arVals, 0);

    param->val.stringArray(as);
}


CORBA::Long PARAM_GetLongFromArray(
         const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long index)
{
    long* arVals = NULL;
    long  nbVal  = 0;
    long  value;


    PARAM_GetLongArray(param, &arVals, &nbVal);

    value = arVals[index];

    delete [] arVals;

    return(value);   
}

 
CORBA::Double PARAM_GetDoubleFromArray(
        const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long index)
{
    double* arVals = NULL;
    long  nbVal  = 0;
    double  value;
 
 
    PARAM_GetDoubleArray(param, &arVals, &nbVal);
 
    value = arVals[index];   
 
    delete [] arVals;
 
    return(value);
}

 
void PARAM_GetStringFromArray(
             ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
             long index, char* value )
{
    char** arVals = NULL;
    long  nbVal  = 0;

 
 
    PARAM_GetStringArray(param, arVals, &nbVal);
 
    strcpy(value, arVals[index]);   
 
    for (long i = 0; i < nbVal; i++)
    {
        delete arVals[i];
    }

    delete [] arVals;
}


CORBA::Long PARAM_GetNbValues(
                   const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param)
{
    return(param->nbValues);
}


void PARAM_SetNbValues(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                       long nbVals)
{
    param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

    param->nbValues = nbVals;
}


ARM_CorbaRequestModule::ARM_CORBA_VALUE_TYPE PARAM_GetValueType(
                const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param)
{
    return(param->val._d());
}


void PARAM_SetValueType(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param,
                        ARM_CorbaRequestModule::ARM_CORBA_VALUE_TYPE vType)
{
    param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

    param->val._d(vType);
}




/*----------------------------------------------------------------------------*/
/*---- End of File ----*/
// EOF %M%
