// ARMModule.h: Definition of the ARMModule class
//
//////////////////////////////////////////////////////////////////////

#ifndef _ACTIVEX_MODULE_H_ 
#define _ACTIVEX_MODULE_H_ 

#include <comdef.h>

//ERROR TYPE
#define ARM_ERROR_IRCURVE 1
#define ARM_ERROR_DFCURVE 2
#define ARM_ERROR_DIFSZCURVE 3
#define ARM_ERROR_SZMATRIX 4
#define ARM_ERROR_PRICER 5
#define ARM_ERROR_FREQ 6
#define ARM_ERROR_DAYCOUNT 7
#define ARM_ERROR_ACCONDEF 8
#define ARM_ERROR_AMONDEF 9
#define ARM_ERROR_INTONDEF 10
#define ARM_ERROR_CORRMATRIX 11
#define ARM_ERROR_SECURITY 12
#define ARM_ERROR_MODEL 13


// this function to be moved. 
class CCString; 
#include <string>
void ERROR_MSG(CCString Err_mess,VARIANT* pRet,long noerr = 0);

class ActiveXModule 
{
public:
	#include "ActiveXLib.h"
public:
	void init(); 
	void release(); 
	STDMETHOD(createErrorInfo)(const std::string& src,const std::string& desc) ;
};

#endif // _ACTIVEX_MODULE_H_ 


