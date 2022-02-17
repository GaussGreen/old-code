#include "firsttobeincluded.h"
#include <ARM\libarm_frometk\VariantTools.h>
#include "ActiveXModule.h"


#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\ARM_local_license.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARMKernel\glob\dates.h>
#include <ARMKernel\glob\calend.h>


//	--------------------------------------------------------------
void 
ActiveXModule::init()
{
	ARM_Date TMP_DATE ;
	int k = isAuthorized();
	CREATE_GLOBAL_OBJECT();
	LOCALARM_IniFileRead();
	ARM_Calendar* cal;
	int rc;
	char calenDarFileName[200];
	rc = ARM_GetCalendarFile(calenDarFileName);
	if ( rc == RET_OK )
	{
	   cal = new ARM_Calendar(calenDarFileName);
	   TMP_DATE.SetCalendar(cal);
	}
}
//	--------------------------------------------------------------
void 
ActiveXModule::release()
{
	if (LOCAL_PERSISTENT_OBJECTS)
		LOCAL_PERSISTENT_OBJECTS->FreeAllObjects();
	deconnection_etoolkit();
}
STDMETHODIMP ActiveXModule::createErrorInfo(const std::string& src_,const std::string& desc_)
{
	ICreateErrorInfo *pcerrinfo;
	IErrorInfo *perrinfo;
	HRESULT hr;

	hr = CreateErrorInfo(&pcerrinfo);
	// set fields here by calling ICreateErrorInfo methods on pcerrinfo
	_bstr_t desc; VariantTools::convert(desc_,desc); 
	_bstr_t src; VariantTools::convert(src_,src); 
	pcerrinfo->SetDescription(desc) ; 
	pcerrinfo->SetSource(src); 
	// pcerrinfo->SetHelpContext(dwhelpcontext);  // and so forth

	hr = pcerrinfo->
		 QueryInterface(IID_IErrorInfo, (LPVOID FAR*) &perrinfo);
	if (SUCCEEDED(hr))
	{
		SetErrorInfo(0, perrinfo);
		perrinfo->Release();
	}
	pcerrinfo->Release();
	// then, eventually...
	return E_FAIL;  // E_FAIL or other appropriate failure code
}
