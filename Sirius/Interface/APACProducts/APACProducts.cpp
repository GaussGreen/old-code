// APACProducts.cpp : Implementation of DLL Exports.
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
#include "APACProducts.h"
#include "APACProducts_i.c"
#include "sirius_i.c"
#include "sirius.h"
#include "Puff.h"

CAPACProductsComModule					_Module;

BEGIN_OBJECT_MAP(ObjectMap)	
OBJECT_ENTRY(CLSID_Puff, CPuff)
END_OBJECT_MAP() 







/////////////////////////////////////////////////////////////////////////////
// Global Functions
//



HRESULT __stdcall APAddTwoNumbers(VARIANT Number1, VARIANT Number2, VARIANT* pResult)
{
	begin_function
	
	map_parameter(Number1, double, f1);
	map_parameter(Number2, double, f2);
	return CComVariant(f1 + f2).Detach(pResult);

	end_function
}


HRESULT __stdcall APAddTwoSpotsFromTwoAssets(VARIANT Asset1, VARIANT Asset2, VARIANT* pResult)
{
	begin_function

	map_object_parameter(Asset1, Asset, h1);
	map_object_parameter(Asset2, Asset, h2);
	double f1 = h1->GetNaturalSpot(h1->GetDateHandle()->GetDate());
	double f2 = h2->GetNaturalSpot(h2->GetDateHandle()->GetDate());
	
	// Return in an array handle, f1, f2 and f1 + f2

	MlEqArrayHandle h = new MlEqArray;

	h->resize(3);
	(*h)[0] = f1;
	(*h)[1] = f2;
	(*h)[2] = f1 + f2;

	map_analytic_to_com(h, Array, sp);

	CComBSTR sHandle;
	_Module.GetSiriusApplication()->InsertObject(sp, VARIANT_TRUE, &sHandle);

	
	return CComVariant(sHandle).Detach(pResult);
	end_function
}







/////////////////////////////////////////////////////////////////////////////
// DLL Entry Point

extern "C"
BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID /*lpReserved*/)
{
    if (dwReason == DLL_PROCESS_ATTACH){
		_Module.Init(ObjectMap, hInstance, &LIBID_APACProducts);
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









