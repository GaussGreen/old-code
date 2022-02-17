// ARMCCSModule.h : Declaration of the ARMCCSModule

#ifndef __ARMCCSMODULE_H_
#define __ARMCCSMODULE_H_

#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// ARMCCSModule
class ATL_NO_VTABLE ARMCommonModule : 
	public CComObjectRootEx<CComMultiThreadModel>,
	public CComCoClass<ARMCommonModule, &CLSID_ARMCCSModule>,
	public ISupportErrorInfo,
	public IDispatchImpl<IARMCCSModule, &IID_IARMCCSModule, &LIBID_ARMCCSLIBLib>
{
public:
	ARMCommonModule()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_ARMCCSMODULE)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(ARMCommonModule)
	COM_INTERFACE_ENTRY(IARMCCSModule)
	COM_INTERFACE_ENTRY(IDispatch)
	COM_INTERFACE_ENTRY(ISupportErrorInfo)
END_COM_MAP()

// ISupportsErrorInfo
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

// IARMCCSModule
public:
	#include "..\ActiveX_ARM\ActiveXLib.h" 
};

#endif //__ARMCCSMODULE_H_
