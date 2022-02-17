// product_cppi.h : Declaration of the CCppi

#ifndef __CPPI_H_
#define __CPPI_H_

#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// CCppi
class ATL_NO_VTABLE CCppi : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CCppi, &CLSID_Cppi>,
	public IDispatchImpl<ICppi, &IID_ICppi, &LIBID_Products>
{
public:
	CCppi()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_CPPI)
DECLARE_NOT_AGGREGATABLE(CCppi)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(CCppi)
	COM_INTERFACE_ENTRY(ICppi)
	COM_INTERFACE_ENTRY(IDispatch)
END_COM_MAP()

// ICppi
public:
};

#endif //__CPPI_H_
