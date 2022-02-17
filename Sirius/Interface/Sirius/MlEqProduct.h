//	MlEqProduct.h :			 Product handle base class.
//							 This is an interface class.
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQPRODUCT_H
#define _MLEQPRODUCT_H

#pragma once

class CParameters;
class MlEqProduct;
typedef RCPtr<MlEqProduct>						MlEqProductHandle;
#include "smart.h"

class MlEqProduct : public RCObject
{
public:	
	MlEqProduct(void);	
	void										CopyToClipboard(void);	
	void										Evaluate(const BSTR& Calculate, CComPtr<IResult>& spResult) const;
	DataSourceEnum								GetDataSource(void) const;
	long										GetDate(void) const;	
	void										GetObject(IDispatch** pVal) const;
	CComPtr<IParameters>						GetParameters(void) const;	
	estring										GetProductType(void) const;				
	HRESULT										GetUnderlyings(CComPtr<IDispatch> spObject, CComPtr<IAssets>& spAssets) const;
	CParameters*								Parameters(void) const;
	void										PutDataSource(DataSourceEnum ds);
	void										PutDate(long nDate);
	void										PutParameters(CComPtr<IParameters> spParameters);
	void										PutProductType(const std::string& szProductType);

protected:	
	void										GetBaseUnderlyings(CComPtr<IDispatch> spObject, CComPtr<IAssets>& spAssets) const;
	std::string									m_szProductType;
	CComPtr<IParameters>						m_spParameters;
	long										m_nDate;
	DataSourceEnum								m_ds;	
};

#endif
