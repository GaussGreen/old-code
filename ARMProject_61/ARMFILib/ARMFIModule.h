// ARMFIModule.h : Declaration of the ARMFIModule

#ifndef __ARMFIMODULE_H_
#define __ARMFIMODULE_H_

#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// ARMFIModule
class ATL_NO_VTABLE ARMCommonModule : 
	public CComObjectRootEx<CComMultiThreadModel>,
	public CComCoClass<ARMCommonModule , &CLSID_ARMFIModule>,
	public ISupportErrorInfo,
	public IDispatchImpl<IARMFIModule, &IID_IARMFIModule, &LIBID_ARMFILib>
{
public:
	ARMCommonModule ()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_ARMFIMODULE)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(ARMCommonModule)
	COM_INTERFACE_ENTRY(IARMFIModule)
	COM_INTERFACE_ENTRY(IDispatch)
	COM_INTERFACE_ENTRY(ISupportErrorInfo)
END_COM_MAP()

// ISupportsErrorInfo
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

// IARMFIModule
public:
	#include "..\ActiveX_ARM\ActiveXLib.h" 

};

#endif //__ARMFIMODULE_H_
