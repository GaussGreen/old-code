//	SiriusRisk.cpp : Implementation of DLL Exports.
//
//	Note: Proxy/Stub Information
//        To build a separate proxy/stub DLL, 
//        run nmake -f QuantLibps.mk in the project directory.
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "resource.h"
#include <initguid.h>
#include "SiriusRisk.h"
#include "sirius_i.c"
#include "SiriusRisk_i.c"
#include "Engine.h"
#include "JetAggregator.h"
#include "Aggregators.h"
#include "risk_delta.h"
#include "risk_price.h"
#include "risk_gamma.h"
#include "risk_vega.h"
#include "risk_Theta.h"
#include "risk_FixedDelta.h"
#include "risk_SwimDelta.h"
#include "risk_gamma5.h"
#include "risk_VegaSkew.h"
#include "risk_DivDelta.h"
#include "risk_rho.h"

CSiriusRiskComModule					_Module;

BEGIN_OBJECT_MAP(ObjectMap)
	OBJECT_ENTRY(CLSID_Engine, CEngine)
	OBJECT_ENTRY(CLSID_JetAggregator, CJetAggregator)
	OBJECT_ENTRY(CLSID_Aggregators, CAggregators)
	OBJECT_ENTRY(CLSID_Delta, CDelta)
	OBJECT_ENTRY(CLSID_Price, CPrice)
	OBJECT_ENTRY(CLSID_Theta, CTheta)
	OBJECT_ENTRY(CLSID_Gamma, CGamma)
	OBJECT_ENTRY(CLSID_Vega, CVega)
	OBJECT_ENTRY(CLSID_FixedDelta, CFixedDelta)
	OBJECT_ENTRY(CLSID_SwimDelta, CSwimDelta)
	OBJECT_ENTRY(CLSID_Gamma5, CGamma5)
	OBJECT_ENTRY(CLSID_VegaSkew, CVegaSkew)
	OBJECT_ENTRY(CLSID_DivDelta, CDivDelta)
	OBJECT_ENTRY(CLSID_Rho, CRho)
END_OBJECT_MAP()

class CSiriusRiskApp : public CWinApp
{
public:
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSiriusRiskApp)
	public:
    virtual BOOL InitInstance();
    virtual int ExitInstance();
	//}}AFX_VIRTUAL

	//{{AFX_MSG(CSiriusRiskApp)		
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

BEGIN_MESSAGE_MAP(CSiriusRiskApp, CWinApp)
	//{{AFX_MSG_MAP(CSiriusRiskApp)		
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

CSiriusRiskApp theApp;

BOOL CSiriusRiskApp::InitInstance()
{
    _Module.Init(ObjectMap, m_hInstance, &LIBID_SiriusRisk);
    return CWinApp::InitInstance();
}

int CSiriusRiskApp::ExitInstance()
{
    _Module.Term();
    return CWinApp::ExitInstance();
}

/////////////////////////////////////////////////////////////////////////////
// Used to determine whether the DLL can be unloaded by OLE

STDAPI DllCanUnloadNow(void)
{
    AFX_MANAGE_STATE(AfxGetStaticModuleState());
    return (AfxDllCanUnloadNow()==S_OK && _Module.GetLockCount()==0) ? S_OK : S_FALSE;
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


