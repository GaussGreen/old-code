// kexception.cpp: interface for the kexception class.
//
//////////////////////////////////////////////////////////////////////

#include "kexception.h"


void HYErrLog(const char* msg)
{
	std::ofstream file;
	if(ErrLog_OnOff == HY_ERR_ON)
	{
		file.open("hy_err.log",std::ios_base::app);
		file <<msg;
		file.close();
	}
}


