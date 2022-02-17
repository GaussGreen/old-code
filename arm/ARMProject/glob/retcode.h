/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: retcode.h,v $
 * Revision 1.5  2003/10/20 15:23:44  mab
 * Comments in C-Style
 *
 * Revision 1.4  2003/10/09 18:50:25  ebenhamou
 * remove ,
 *
 * Revision 1.3  2003/10/07 15:16:54  ebenhamou
 * added some more type
 *
 * Revision 1.2  2003/10/02 17:51:19  ebenhamou
 * added new retcode
 *
 *
 */


/*! \file retcode.h
 *
 *  \brief defined some general type 
 *
 *	\version revisited for the generic pricer
 *	\date October 2003
 */

 
#ifndef _RETCODE_H
#define _RETCODE_H



/* this is an enumeration of possible type
 * names are explicit!
 * string are string or char*
 */

typedef enum
{
    ARM_ERR = -1,
	ARM_UNKNOWN = 0,
    ARM_DOUBLE =1,		/* defined with value to make the two enums equal */
	ARM_DOUBLE_TYPE = 1,/* defined with value to make the two enums equal */
    ARM_INT = 2,		/* defined with value to make the two enums equal */
	ARM_INT_TYPE = 2,	/* defined with value to make the two enums equal */
    ARM_STRING = 3,		/* defined with value to make the two enums equal */
	ARM_STRING_TYPE = 3,/* for type compliance                            */
	ARM_DATE_TYPE,
	ARM_TOBECOMPUTED,
	ARM_BOOL_TYPE,
	ARM_MISSING_TYPE
		
} ARM_VALUE_TYPE;



typedef union _ARM_VALUE
{
    double doubleVal;
    int    intVal;
    char   stringVal[254];
}
ARM_VALUE;


#define	ARM_EXCEPTION_MSG_MAX_SIZE	1024 


typedef struct _ARM_RET_CODE
{
    long           code;
    char           msg[ARM_EXCEPTION_MSG_MAX_SIZE];
    ARM_VALUE_TYPE valueType;
    ARM_VALUE      value;
}
ARM_RET_CODE;



#define GET_ERR(code, msg) \
    RC_GetErr(code, msg, __FILE__, __ARM_FUNC__, __LINE__)

#define SET_ERR(rc, code, msg) \
    SetErr(rc, code, msg, __FILE__, __ARM_FUNC__, __LINE__)


extern void SetDoubleValue(ARM_RET_CODE* retCode, double val);
extern ARM_RET_CODE RC_GetDoubleValue(double val);

extern void SetIntValue(ARM_RET_CODE* retCode, int val);
extern ARM_RET_CODE RC_GetIntValue(int val);

extern void SetStringValue(ARM_RET_CODE* retCode, char* val);
extern ARM_RET_CODE RC_GetStringValue(char* val);

extern void SetErr(ARM_RET_CODE* retCode, long code, char* msg, 
            char* file, char* func, int line);

extern ARM_RET_CODE RC_GetErr(long code, char* msg,
                       char* file, char* func, int line);
              

extern ARM_RET_CODE RC_GetOK(void);

extern void RC_Print(ARM_RET_CODE rc);




#endif 
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
