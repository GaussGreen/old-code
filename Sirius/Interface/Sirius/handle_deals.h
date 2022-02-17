//	handle_deals.h : Declaration of the CDeals
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DEALS_H_
#define __DEALS_H_

#include "resource.h"
#include "comobjectcollectionserialisable.h"		// serialisable template collection class

class ATL_NO_VTABLE CDeals : 
	public CComObjectCollectionSerialisable<IDeal, CDeals, IDeals, &CLSID_Deals, &IID_IDeals, &LIBID_Sirius, &CLSID_Deal, IDR_DEALS>,
	public IDispatchImpl<IEvaluatable, &IID_IEvaluatable, &LIBID_Sirius>
{
public:
	static HRESULT						Load(const std::string& szBook, DataSourceEnum ds, long nDate, CComPtr<IDeals>& spDeals);

	BEGIN_COM_MAP(CDeals)
		COM_INTERFACE_ENTRY(IDeals)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IDeals)
		COM_INTERFACE_ENTRY(IEvaluatable)
	END_COM_MAP()

protected:
	STDMETHOD(Evaluate)(/*[in, defaultvalue("Price")]*/ BSTR Calculate, /*[out, retval]*/ IResult** pVal);
	STDMETHOD(GetFxUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(GetUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(get_Date)(/*[out, retval]*/ DATE* pVal);
	STDMETHOD(GetVolatilityStructures)(/*[out, retval]*/ IVolatilityStructures** pVal);
	STDMETHOD(GetZeroCurves)(/*[out, retval]*/ IZeroCurves** pVal);
	STDMETHOD(PutDataSource)(/*[in]*/ DataSourceEnum newVal);
	STDMETHOD(PutDate)(/*[in]*/ DATE newVal);
};

#endif
