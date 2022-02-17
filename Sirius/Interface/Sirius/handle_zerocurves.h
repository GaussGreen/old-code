//	handle_zerocurves.h : Declaration of the CZeroCurves
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __ZEROCURVES_H_
#define __ZEROCURVES_H_

#include "resource.h"
#include "comobjectcollectionserialisable.h"

class ATL_NO_VTABLE CZeroCurves : public CComObjectCollectionSerialisable<IZeroCurve, CZeroCurves, IZeroCurves, &CLSID_ZeroCurves, &IID_IZeroCurves, &LIBID_Sirius, &CLSID_ZeroCurve, IDR_ZEROCURVES>
{
public:
	static HRESULT						Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<IZeroCurves>& spZeroCurves);
	STDMETHOD(Reset)(void);
	STDMETHOD(Shift)(/*[in]*/ double Amount);
	STDMETHOD(Stick)(void);

protected:
	HRESULT								LoadGDAObject(const std::string& szName, long nDate, CComPtr<IZeroCurve>& spObject);
};

#endif
