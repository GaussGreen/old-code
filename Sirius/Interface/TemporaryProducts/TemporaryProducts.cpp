// TemporaryProducts.cpp : Implementation of DLL Exports.
//
//
// Note: Proxy/Stub Information
//      To build a separate proxy/stub DLL, 
//      run nmake -f QuantLibps.mk in the project directory.
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "resource.h"
#include <initguid.h>
#include "TemporaryProducts.h"
#include "TemporaryProducts_i.c"
#include "sirius_i.c"
#include "product_put.h"
#include "product_mattoption.h"
#include "product_iguana.h"
#include "montecarlo.h"
#include "mleqparameterlist.h"
#include "product_funding.h"
#include "sirius.h"

CTemporaryProductsComModule					_Module;

BEGIN_OBJECT_MAP(ObjectMap)	
	OBJECT_ENTRY(CLSID_ProductPut, CProductPut)
	OBJECT_ENTRY(CLSID_MattOption, CMattOption)
	OBJECT_ENTRY(CLSID_Iguana, CIguana)
	OBJECT_ENTRY(CLSID_Funding, CFunding)
END_OBJECT_MAP() 

/////////////////////////////////////////////////////////////////////////////
// DLL Entry Point

extern "C"
BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID /*lpReserved*/)
{
    if (dwReason == DLL_PROCESS_ATTACH){
		_Module.Init(ObjectMap, hInstance, &LIBID_TemporaryProductsLib);
        DisableThreadLibraryCalls(hInstance);
    } else if (dwReason == DLL_PROCESS_DETACH){
		_Module.Term();
	}
    return TRUE;
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


