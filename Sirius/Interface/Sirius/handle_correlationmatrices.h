//	handle_correlationmatrices.h : Declaration of the CCorrelationMatrices
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __CORRELATIONMATRICES_H_
#define __CORRELATIONMATRICES_H_

#include "resource.h"
#include "comobjectcollectionserialisable.h"

class ATL_NO_VTABLE CCorrelationMatrices : public CComObjectCollectionSerialisable<ICorrelationMatrix, CCorrelationMatrices, ICorrelationMatrices, &CLSID_CorrelationMatrices, &IID_ICorrelationMatrices, &LIBID_Sirius, &CLSID_CorrelationMatrix, IDR_CORRELATIONMATRICES> 
{
public:
	static HRESULT								Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<ICorrelationMatrices>& spCorrelationMatrices);

protected:
	void										CheckAddHandle(std::string* pszHandle);
	HRESULT										SpecialAdd(CComPtr<ICorrelationMatrix> spSingular, const CComBSTR& sKey);
};

#endif
