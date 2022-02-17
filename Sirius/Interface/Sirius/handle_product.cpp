//	handle_product.cpp :  Implementation of CProduct
//
//	author:				  David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_product.h"
#include "handle_assets.h"
#include "handle_asset.h"
#include "MlEqObjects.h"
#include "MlEqDictionary.h"
#include "parameters.h"
#include "siriusapplication.h"
#include "comobjectcollectionserialisabledefaulter.h"
#include "comobjectcollectionfunctions.h"
#include "handle_zerocurves.h"


STDMETHODIMP CProduct::CopyToClipboard()
{
	begin_function	
	CParameterMap pmValue;
	GetParameterMap(&pmValue);
	return pmValue.CopyToClipboard();	
	end_function
}

STDMETHODIMP CProduct::Evaluate(BSTR Calculate, IResult** pVal)
{
	CComPtr<IResult>					spResult;	

	begin_function
	m_h->Evaluate(Calculate, spResult);
	return spResult.CopyTo(pVal);
	end_function
}

//	returns a collection of foreign exchange assets on which the product depends
STDMETHODIMP CProduct::GetFxUnderlyings(IAssets** pVal)
{
	HRESULT								hr;	
	CComPtr<IDispatch>					spObject;
	CComPtr<IAssets>					spAssets;
	CComPtr<IAssets>					spFxAssets;
	CAssets*							pAssets;
	long								nCount = 0L;
			
	begin_function
	m_h->GetObject(&spObject);
	m_h->GetUnderlyings(spObject, spAssets);
	if ((hr = spFxAssets.CoCreateInstance(CLSID_Assets)) || (hr = spAssets->get_Count(&nCount))) return hr;
	if (!(pAssets = dynamic_cast<CAssets*>(spFxAssets.p))) return E_FAIL;
	
	for (long n = 1; n <= nCount; n++){
		CComPtr<IAsset>					spAsset, spCurrencyAsset, spCompositeCurrencyAsset, spPayCurrencyAsset;				
		if (spAssets->get_Item(CComVariant(n), &spAsset)) continue;
		spAsset->get_CurrencyAsset(&spCurrencyAsset);
		spAsset->get_PayCurrencyAsset(&spPayCurrencyAsset);
		spAsset->get_CompositeCurrencyAsset(&spCompositeCurrencyAsset);
		pAssets->MergeSingular(spCurrencyAsset);
		pAssets->MergeSingular(spPayCurrencyAsset);
		pAssets->MergeSingular(spCompositeCurrencyAsset);
	}
	return spFxAssets.CopyTo(pVal);
	end_function
}

STDMETHODIMP CProduct::GetObject(IDispatch** pVal)
{
	CComPtr<IDispatch>					spObject;
	
	begin_function	
	m_h->GetObject(&spObject);
	return spObject.CopyTo(pVal);
	end_function
}

void CProduct::GetParameterMap(CParameterMap* ppm) const
{
	MlEqDictionary<CComVariant>			dict;
	
	m_h->Parameters()->GetDictionary(&dict);
	dict.GetValue(ppm);
	
	if (m_h->GetDataSource() != NoDataSource){
		ppm->InsertRow(0);
		ppm->SetValue(0, 0, CStringResource(IDS_HEADER_DATA_SOURCE));
		ppm->SetValue(0, 1, CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()));
	}	
	if (m_h->GetDate()){
		ppm->InsertRow(0);
		ppm->SetValue(0, 0, CStringResource(IDS_HEADER_DATE));
		ppm->SetValue(0, 1, MlEqDate(m_h->GetDate()).GetString());
	}
	ppm->InsertRow(0);
	ppm->SetValue(0, 0, CStringResource(IDS_HEADER_PRODUCT_TYPE));
	ppm->SetValue(0, 1, m_h->GetProductType());	
}

STDMETHODIMP CProduct::get_Parameters(IParameters** pVal)
{
	begin_function
	return m_h->GetParameters().CopyTo(pVal);
	end_function
}

STDMETHODIMP CProduct::GetPayCurrency(CurrencyEnum* pVal)
{
	begin_function	
	CComPtr<IAssets>					spUnderlyings;
	CurrencyEnum						ceRet = NoCurrency;
	long								nCount;
	
	if (GetUnderlyings(&spUnderlyings)) propagate_error;
	if (spUnderlyings->get_Count(&nCount)) propagate_error;
	if (!nCount) throw "The product does not have any underlyings";
	for (CComVariant vItem = 1L; vItem.lVal <= nCount; vItem.lVal++){
		CurrencyEnum		ce;
		CComPtr<IAsset>		spAsset;
		if (spUnderlyings->get_Item(vItem, &spAsset) || spAsset->get_PayCurrency(&ce)) propagate_error;
		if (ceRet == NoCurrency){
			ceRet = ce;
		} else if (ceRet != ce){
			throw "The product has more than one pay currency";
		}
	}
	*pVal = ceRet;
	end_function
}

//	returns a list of parameter names that the product type implementation object expects
STDMETHODIMP CProduct::GetRequiredParameters(VARIANT *pVal)
{	
	HRESULT								hr;
	CParameterMap						pmRet, pm;

	pmRet.SetSize(3, 1);
	pmRet.SetValue(0, 0, CStringResource(IDS_HEADER_PRODUCT_TYPE));
	pmRet.SetValue(1, 0, CStringResource(IDS_HEADER_DATE));
	pmRet.SetValue(2, 0, CStringResource(IDS_HEADER_DATA_SOURCE));
	if (hr = g_pApplication->GetObjectManager().GetObjectProperties(INVOKE_PROPERTYPUT, m_h->GetProductType(), NULL, &pm)) return hr;
	return pmRet.AddToEnd(pm).GetValue(pVal);
}

STDMETHODIMP CProduct::GetUnderlyings(IAssets** pVal)
{				
 	CComPtr<IDispatch>					spObject;
	CComPtr<IAssets>					spUnderlyings;
		
	begin_function
	m_h->GetObject(&spObject);
	m_h->GetUnderlyings(spObject, spUnderlyings);
	return spUnderlyings.CopyTo(pVal);
	end_function
}

STDMETHODIMP CProduct::GetVolatilityStructures(IVolatilityStructures** pVal)
{
	HRESULT								hr;
	CComPtr<IAssets>					spUnderlyings;
	CComPtr<IVolatilityStructures>		spVolatilityStructures;

	if (hr = GetUnderlyings(&spUnderlyings)) return hr;
	if (hr = spUnderlyings->GetVolatilityStructures(&spVolatilityStructures)) return hr;
	return spVolatilityStructures.CopyTo(pVal);
}

//	returns the value stored in a given property name
HRESULT CProduct::GetValue(const std::string& szName, VARIANT *pVal)
{
	begin_function
	HRESULT								hr;	
	CComPtr<IParameter>					spParameter;
	
	if (!(hr = m_h->GetParameters()->get_Item(estring::GetValue(szName), &spParameter))){
		return spParameter->get_Value(pVal);
	}
	// We need special cases for ProductType, Date and DataSource
	if (!CStringResource(IDS_HEADER_PRODUCT_TYPE).CompareNoCase(szName)){
		return m_h->GetProductType().GetValue(pVal);
	} else if (!CStringResource(IDS_HEADER_DATE).CompareNoCase(szName)){
		return CComVariant(m_h->GetDate()).Detach(pVal);
	} else if (!CStringResource(IDS_HEADER_DATA_SOURCE).CompareNoCase(szName)){		
		return CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()).GetValue(pVal);
	}
	return hr;		
	end_function
}
 
STDMETHODIMP CProduct::get_Value(VARIANT* pVal)
{
	// This implementation returns a flat matrix for ease of presentation.
	begin_function
	CParameterMap						pmValue;
	GetParameterMap(&pmValue);
	return pmValue.GetValue(pVal);	
	end_function
}

STDMETHODIMP CProduct::GetZeroCurves(IZeroCurves** pVal)
{	
	return CComObjectCollectionFunctions<IProduct>(CComQIPtr<IEvaluatable>(this)).GetZeroCurves(pVal);	
}

STDMETHODIMP CProduct::HasParameter(BSTR Name, VARIANT_BOOL* pVal)
{
	bool b = m_h->Parameters()->IsInCollection(Name);
	*pVal = b ? VARIANT_TRUE : VARIANT_FALSE;
	return S_OK;
}

STDMETHODIMP CProduct::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IProduct};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CProduct::PutDataSource(/*[in]*/ DataSourceEnum newVal)
{
	return put_DataSource(newVal);
}

STDMETHODIMP CProduct::PutDate(DATE newVal)
{
	return put_Date(newVal);
}

//	Set the parameters and product type of this object from an input Sirius product.
//  That it, this function performs the opposite of GetObject.
STDMETHODIMP CProduct::PutObject(IDispatch* Product)
{
	begin_function
	CComPtr<IDispatch>					spProduct(Product);
	std::string							szProductType;		
	CParameterMap						pmProductData;
	CComPtr<IParameters>				spParameters;
	
	if (CParameterMap::GetObjectName(spProduct, &szProductType)) propagate_error;

	// Can we actually create that product?
	if (!g_pApplication->GetObjectManager().IsSiriusProduct(szProductType)){
		throw "You cannot attach an object with name '" + szProductType + "' to the product";
	} else {	
		CComVariant vProductData;
		g_pApplication->GetObjectManager().ImplementGetValue(spProduct, &vProductData);		
		pmProductData.SetValue(vProductData);
		if (pmProductData.GetCols() != 2) throw "The input object's state is invalid for this operation";
	}		
	
	// Create the parameters collection
	spParameters.CoCreateInstance(CLSID_Parameters);	
	for (long nRow = 0; nRow < pmProductData.GetRows(); nRow++){
		CComPtr<IParameter> spParameter;				
		CComBSTR			sName;
		CComVariant			vValue;
		
		spParameter.CoCreateInstance(CLSID_Parameter);
		pmProductData.GetValue(nRow, 0, &sName);
		pmProductData.GetValue(nRow, 1, &vValue);								
		spParameter->put_Name(sName);
		spParameter->put_Value(vValue);
		if (spParameters->Add(CComVariant(sName), spParameter)) propagate_error;
	}

	m_h->PutParameters(spParameters);
	m_h->PutProductType(szProductType);	
	end_function
}

STDMETHODIMP CProduct::put_Parameters(IParameters* newVal)
{
	begin_function
	m_h->PutParameters(newVal);
	end_function
}

HRESULT CProduct::PutValue(const std::string& szName, const CComVariant& Value)
{	
	begin_function		
	
	// We need special cases for ProductType, Date and DataSource since
	// we don't want to add them to the dictionary.
	bool bProductType = CStringResource(IDS_HEADER_PRODUCT_TYPE).CompareNoCaseAndSpace(szName) ? false : true;
	bool bDate = CStringResource(IDS_HEADER_DATE).CompareNoCaseAndSpace(szName) ? false : true;
	bool bDataSource = CStringResource(IDS_HEADER_DATA_SOURCE).CompareNoCaseAndSpace(szName) ? false : true;
	if (bProductType || bDate || bDataSource){			
		if (bProductType){
			map_parameter(Value, estring, szProductType);
			m_h->PutProductType(szProductType);			
			return S_OK;
		} else if (bDate){
			map_parameter(Value, long, nDate);
			m_h->PutDate(nDate);
			return S_OK;
		} else if (bDataSource){
			map_enum_parameter(Value, DataSourceEnum, ds);
			m_h->PutDataSource(ds);
			return S_OK;
		} else {
			ATLASSERT(false);
		}
	}

	// General case from here		
	CComPtr<IParameter> spParameter;	
	spParameter.CoCreateInstance(CLSID_Parameter);
	spParameter->put_Name(estring::GetBSTR(szName));	
	spParameter->put_Value(Value);
	m_h->GetParameters()->Add(estring::GetValue(szName), spParameter);	
	end_function
}

STDMETHODIMP CProduct::put_Value(VARIANT newVal)
{	
	/*parameter list is 0 - Product Type
						1 - DataSourceOpt
						2 - DateOpt
						3 - Vector of parameter names
						4 - Vector of paraméter values*/
		
	begin_function	
	HRESULT								hr;
	std::vector<CParameterMap>			vpm;	
	MlEqProductHandle					h = new MlEqProduct;
			
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	if (vpm.size() == 1){
		// We need to support this case since get_Value returns a flat case.			
		MlEqDictionary<CComVariant>		dict;		
		dict.PutValue(vpm[0]);
				
		// Abstract the product type, data source and date from the list
		CComVariant ProductType, DataSource, Date;
		dict.GetValue(IDS_HEADER_PRODUCT_TYPE, &ProductType, true);
		map_parameter(ProductType, std::string, szProductType);
		dict.Delete(IDS_HEADER_PRODUCT_TYPE, true);
		dict.GetValue(IDS_HEADER_DATA_SOURCE, &DataSource, false);
		map_optional_enum_parameter(DataSource, DataSourceEnum, ds, Last);
		dict.Delete(IDS_HEADER_DATA_SOURCE, false);
		dict.GetValue(IDS_HEADER_DATE, &Date, false);
		map_optional_parameter(Date, long, nDate, 0);
		dict.Delete(IDS_HEADER_DATE, false);
		
		// The remaining elements in the dictionary form the parameters of the product
		CComPtr<IParameters>		spParameters;
		spParameters.CoCreateInstance(CLSID_Parameters);
		for (std::map<std::string, CComVariant>::const_iterator it = dict.GetMap().begin(); it != dict.GetMap().end(); it++){
			CComPtr<IParameter> spParameter;			
			spParameter.CoCreateInstance(CLSID_Parameter);
			spParameter->put_Name(estring::GetBSTR(dict.GetLCMap().find(it->first)->second));
			spParameter->put_Value(it->second);
			spParameters->Add(estring(dict.GetLCMap().find(it->first)->second), spParameter);
		}
		
		h->PutProductType(szProductType);
		h->PutDataSource(ds);
		h->PutDate(nDate);
		h->PutParameters(spParameters);
	} else if (vpm.size() == 5){
		map_parameter(vpm[0], estring, ProductType);
		map_optional_enum_parameter(vpm[1], DataSourceEnum, DataSourceOpt, Last);
		map_optional_parameter(vpm[2], long, Date, 0);

		if (!vpm[3].IsVector()) return CParameterMap::ReturnErrorR(IDS_VECTOR_NAMES_INVALID);
		if (!vpm[4].IsVector()) return CParameterMap::ReturnErrorR(IDS_VECTOR_VALUES_INVALID);
		if (vpm[3].GetCols() != 1 && vpm[3].Transpose()) ATLASSERT(false);
		if (vpm[4].GetCols() != 1 && vpm[4].Transpose()) ATLASSERT(false);
		if (vpm[3].GetRows() < vpm[4].GetRows()) return CParameterMap::ReturnErrorR(IDS_VECTOR_VALUES_INVALID);
		vpm[4].SetRows(vpm[3].GetRows());	// this allows for the case where the values parameter map has blanks at the end
						
		CComPtr<IParameters>		spParameters;
		spParameters.CoCreateInstance(CLSID_Parameters);
		for (long nRow = 0; nRow < vpm[3].GetRows(); nRow++){
			estring					szName;
			CComVariant				v;
			CComPtr<IParameter>		spParameter;
			if (vpm[3].GetValue(nRow, 0, &szName)) return CParameterMap::ReturnErrorR(IDS_VECTOR_NAMES_INVALID);
			if (vpm[4].GetValue(nRow, 0, &v)) return CParameterMap::ReturnErrorR(IDS_VECTOR_VALUES_INVALID);
			spParameter.CoCreateInstance(CLSID_Parameter);
			spParameter->put_Name(szName.GetBSTR());
			spParameter->put_Value(v);
			spParameters->Add(szName, spParameter);
		}
		
		h->PutProductType(ProductType);
		h->PutDataSource(DataSourceOpt);
		h->PutDate(Date);
		h->PutParameters(spParameters);
	} else {
		return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IProduct);
	}

	m_h = h;
	end_function
}