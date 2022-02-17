//	handle_volatilitystructures.h : Declaration of the CVolatilityStructures
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __VOLATILITYSTRUCTURES_H_
#define __VOLATILITYSTRUCTURES_H_

#include "resource.h"
#include "comobjectcollectionserialisable.h"

class ATL_NO_VTABLE CVolatilityStructures : public CComObjectCollectionSerialisable<IVolatilityStructure, CVolatilityStructures, IVolatilityStructures, &CLSID_VolatilityStructures, &IID_IVolatilityStructures, &LIBID_Sirius, &CLSID_VolatilityStructure, IDR_VOLATILITYSTRUCTURES>
{
public:
	static HRESULT						Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<IVolatilityStructures>& spVolatilityStructures);
	STDMETHOD(Reset)(void);
	STDMETHOD(Shift)(/*[in]*/ double Amount);
	STDMETHOD(Stick)(void);
};

#endif
