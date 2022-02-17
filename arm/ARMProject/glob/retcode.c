/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : retcode.c                                                    */
/*                                                                            */
/* DESCRIPTION : API Return codes Management                                  */
/*                                                                            */
/* DATE        : Tue Oct 29 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
 
 



#include <string.h>
#include <stdio.h>



#include "retcode.h"





/*----------------------------------------------------------------------------*/


void RC_Print(ARM_RET_CODE rc)
{
}



void SetDoubleValue(ARM_RET_CODE* retCode, double val)
{
    memset(retCode, 0, sizeof(ARM_RET_CODE));

    retCode->valueType = ARM_DOUBLE;

    retCode->value.doubleVal = val;

    retCode->code = 0;
}



ARM_RET_CODE RC_GetDoubleValue(double val)
{
    ARM_RET_CODE rc;
 
 
 
    SetDoubleValue(&rc, val);
 
    return(rc);
}



void SetIntValue(ARM_RET_CODE* retCode, int val)
{
    memset(retCode, 0, sizeof(ARM_RET_CODE));

    retCode->valueType = ARM_INT;
 
    retCode->value.intVal = val;
 
    retCode->code = 0;
}



ARM_RET_CODE RC_GetIntValue(int val)
{
    ARM_RET_CODE rc;

 

    SetIntValue(&rc, val);

    return(rc);
}



void SetStringValue(ARM_RET_CODE* retCode, char* val)
{
    memset(retCode, 0, sizeof(ARM_RET_CODE));

    retCode->valueType = ARM_STRING;
 
    strcpy(retCode->value.stringVal, val);
 
    retCode->code = 0;
}



ARM_RET_CODE RC_GetStringValue(char* val)
{
    ARM_RET_CODE rc;
 


    SetStringValue(&rc, val);

    return(rc);
}



ARM_RET_CODE RC_GetErr(long code, char* msg,
                       char* file, char* func, int line)
{
    ARM_RET_CODE rc;
 
    

    SetErr(&rc, code, msg, file, func, line);

    fprintf(stderr, "\n");
    fprintf(stderr, rc.msg);
    fprintf(stderr, "\n");

    return(rc);
}



void SetErr(ARM_RET_CODE* retCode, long code, char* msg, 
            char* file, char* func, int line)
{
    char buf[500];


    memset(retCode, 0, sizeof(ARM_RET_CODE));

    sprintf(buf, "?-----> ERROR : %ld : %s, FILE : %s, FUNC : %s, LINE : %d",
                   code, msg, file, func, line);

    retCode->valueType = ARM_ERR;

    strcpy(retCode->msg, buf);

    retCode->code = code;
}
              


ARM_RET_CODE RC_GetOK(void)
{
    ARM_RET_CODE rc;
 
   
    memset(&rc, 0, sizeof(ARM_RET_CODE));
 
    rc.valueType = ARM_INT;
    rc.value.intVal = 0;

    return(rc);
}








/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
