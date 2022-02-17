// local_DLLARM.cpp : Implementation of DLL Exports.


// Note: Proxy/Stub Information
//      To build a separate proxy/stub DLL, 
//      run nmake -f local_DLLARMps.mk in the project directory.

#include "stdafx.h"
#include "resource.h"
#include <initguid.h>
#include "local_DLLARM.h"

#include "local_DLLARM_i.c"
#include "ARMModule.h"

#include <comdef.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\ARM_local_license.h>

#include <ARM\libarm_frometk\ARM_local_etoolkit.h>

#include <glob\dates.h>
#include <glob\calend.h>

CComModule _Module;

ARM_Date TMP_DATE;


BEGIN_OBJECT_MAP(ObjectMap)
OBJECT_ENTRY(CLSID_ARMModule, ARMCommonModule)
END_OBJECT_MAP()

/////////////////////////////////////////////////////////////////////////////
// DLL Entry Point

extern "C"
BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID /*lpReserved*/)
{
    if (dwReason == DLL_PROCESS_ATTACH)
    {
        _Module.Init(ObjectMap, hInstance, &LIBID_LOCAL_DLLARMLib);
        DisableThreadLibraryCalls(hInstance);

		setlocale(LC_ALL,"English");

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
    else if (dwReason == DLL_PROCESS_DETACH)
    {
		_Module.Term();

		setlocale(LC_ALL,"C");

		if (LOCAL_PERSISTENT_OBJECTS)
			LOCAL_PERSISTENT_OBJECTS->FreeAllObjects();

		deconnection_etoolkit();
	}

    return TRUE;    // ok
}

/////////////////////////////////////////////////////////////////////////////
// Used to determine whether the DLL can be unloaded by OLE

STDAPI DllCanUnloadNow(void)
{
    return (_Module.GetLockCount()==0) ? S_OK : S_FALSE;
}

/////////////////////////////////////////////////////////////////////////////
// Returns a class factory to create an object of the requested type

STDAPI DllGetClassObject(REFCLSID rclsid, REFIID riid, LPVOID* ppv)
{
    return _Module.GetClassObject(rclsid, riid, ppv);
}

/////////////////////////////////////////////////////////////////////////////
// DllRegisterServer - Adds entries to the system registry

STDAPI DllRegisterServer(void)
{
    // registers object, typelib and all interfaces in typelib
    return _Module.RegisterServer(TRUE);
}

/////////////////////////////////////////////////////////////////////////////
// DllUnregisterServer - Removes entries from the system registry

STDAPI DllUnregisterServer(void)
{
    return _Module.UnregisterServer(TRUE);
}


