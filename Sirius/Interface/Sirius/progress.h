//	progress.h : Declaration of the CProgress
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PROGRESS_H_
#define __PROGRESS_H_

#include "resource.h"

class ATL_NO_VTABLE CProgress : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CProgress, &CLSID_Progress>,
	public ISupportErrorInfo,
	public IDispatchImpl<IProgress, &IID_IProgress, &LIBID_Sirius>
{
public:
	CProgress(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_PROGRESS)
	DECLARE_NOT_AGGREGATABLE(CProgress)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CProgress)
		COM_INTERFACE_ENTRY(IProgress)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	HRESULT FinalConstruct();
	HRESULT	FinalRelease();
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_PercentCompleted)(/*[out, retval]*/ short *pVal);
	STDMETHOD(put_PercentCompleted)(/*[in]*/ short newVal);

protected:			
	HRESULT								GetWindowText(std::string* psz) const;	
	estring								m_szInitial;
	HWND								m_hWnd;
	short int							m_nPercent;
};

#endif
