//	handle_positions.h : Declaration of the CPositions
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __POSITIONS_H_
#define __POSITIONS_H_

#include "resource.h"
#include "comobjectcollectionserialisable.h"		// serialisable template collection class

class ATL_NO_VTABLE CPositions : 
	public CComObjectCollectionSerialisable<IPosition, CPositions, IPositions, &CLSID_Positions, &IID_IPositions, &LIBID_Sirius, &CLSID_Position, IDR_POSITIONS>,
	public IDispatchImpl<IEvaluatable, &IID_IEvaluatable, &LIBID_Sirius>
{	
public:
	static HRESULT						Load(const std::string& szBook, DataSourceEnum ds, long nDate, CComPtr<IPositions>& spPositions);

	BEGIN_COM_MAP(CPositions)
		COM_INTERFACE_ENTRY(IPositions)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IPositions)
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
