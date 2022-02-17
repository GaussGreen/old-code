//	handle_matrix.h : Declaration of the CInterfaceMatrix
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MATRIX_H_
#define __MATRIX_H_

#include "resource.h"
#include "mleqmatrix.h"					// for m_h

class ATL_NO_VTABLE CInterfaceMatrix : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CInterfaceMatrix, &CLSID_Matrix>,
	public ISupportErrorInfo,
	public IDispatchImpl<IMatrix, &IID_IMatrix, &LIBID_Sirius>
{
public:
	CInterfaceMatrix() : m_h(new MlEqMatrix){}
	DECLARE_REGISTRY_RESOURCEID(IDR_MATRIX)
	DECLARE_NOT_AGGREGATABLE(CInterfaceMatrix)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CInterfaceMatrix)
		COM_INTERFACE_ENTRY(IMatrix)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CInterfaceMatrix, MlEqMatrix);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
};

#endif
