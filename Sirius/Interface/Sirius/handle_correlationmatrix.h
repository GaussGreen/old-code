//	handle_correlationmatrix.h : Declaration of the CCorrelationMatrix
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __CORRELATIONMATRIX_H_
#define __CORRELATIONMATRIX_H_

#include "resource.h"
#include "mleqcorrelationmatrix.h"		// for m_h

#include "MlEqHandles.h"

class ATL_NO_VTABLE CCorrelationMatrix : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CCorrelationMatrix, &CLSID_CorrelationMatrix>,
	public ISupportErrorInfo,
	public IDispatchImpl<ICorrelationMatrix, &IID_ICorrelationMatrix, &LIBID_Sirius>
{
public:
	static estring						s_szCorrelationMatrixName;		// the interface sets all correlation matrix names to this value
	
	CCorrelationMatrix() : m_h(new MlEqCorrelationMatrix) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_CORRELATIONMATRIX)
	DECLARE_NOT_AGGREGATABLE(CCorrelationMatrix)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CCorrelationMatrix)
		COM_INTERFACE_ENTRY(ICorrelationMatrix)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CCorrelationMatrix, MlEqCorrelationMatrix);
	declare_serialisable;
	
	HRESULT								FinalConstruct(void);
	double								GetCorrelation(const std::string& sz1, const std::string& sz2);
	static HRESULT						Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<ICorrelationMatrix>& spCorrelationMatrix);

	STDMETHOD(Add)(/*[in]*/ ICorrelationMatrix* CorrelationMatrix, VARIANT_BOOL Replace);
	STDMETHOD(Clear)();
	STDMETHOD(CopyToClipboard)();
	STDMETHOD(GetCorrelation)(/*[in]*/ BSTR Name1, /*[in]*/ BSTR Name2, /*[out, retval]*/ double* pVal);
	STDMETHOD(get_Matrix)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(get_Name)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(put_Matrix)(/*[in]*/ VARIANT newVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);		
	STDMETHOD(Save)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(SetCorrelation)(/*[in]*/ BSTR Name1, /*[in]*/ BSTR Name2, /*[in]*/ double Correlation);
	
protected:
	declare_member_variable_rekey(DATE, Date);
	declare_member_variable_rekey(DataSourceEnum, DataSource);
	void								GetParameterMap(CParameterMap* ppmValue);
	MlEqCorrelationMatrixHandle			PutMap(CParameterMap* ppm, const std::string& szName, DataSourceEnum ds, long nDate) const;
};
 
#endif
