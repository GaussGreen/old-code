#ifndef _ARM_LOCAL_INTERFACE_H
#define _ARM_LOCAL_INTERFACE_H


#include <libCCxll\CCxll.h>

//#include "ARM_init.h"

							

#define ARM_ERR()				XL_result.xltype = xltypeStr;\
								XL_result.val.str = "\007ARM_ERR";\
								SetCurCellErrValue (C_result.getMsg ());


#define PXL_ARM_ERR()			XL_result.xltype = xltypeStr;\
								XL_result.val.str = "\007ARM_ERR";

//MSG_printf_message (MSG_TRACE, "===>IN: %s, %d", file, line);


#define ARM_NOCALCIFWIZ()		if(XL_IsCalledByFuncWiz ())\
								{\
									XL_result.xltype = xltypeStr;\
									XL_result.val.str = "";\
									return (LPXLOPER)&XL_result;\
								}


#endif

/*---- End Of File ----*/