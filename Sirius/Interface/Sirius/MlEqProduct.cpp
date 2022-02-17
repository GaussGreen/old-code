//	MlEqProduct.cpp :		 Implementation of MlEqProduct
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqProduct.h"
#include "MlEqDictionary.h"
#include "parameters.h"
#include "siriusapplication.h"
#include "handle_assets.h"
#include "handle_result.h"
#include "handle_correlationmatrix.h"
#include "MlEqObjects.h"
#include "comobjectcollectionserialisabledefaulter.h"

// Note. If you write a new function here, you need to carefully consider whether you need to
//       use the macro "pin_date". If you have any doubt, then use it; overuse has no effect. 
//       Consider yourself warned!

MlEqProduct::MlEqProduct(void) : m_nDate(0), m_ds(NoDataSource)
{
}

void MlEqProduct::CopyToClipboard(void)
{
	MlEqDictionary<CComVariant>			dict;
	Parameters()->GetDictionary(&dict);
	dict.CopyToClipboard();
}

void MlEqProduct::Evaluate(const BSTR& Calculate, CComPtr<IResult>& spResult) const
{			
	pin_date(m_ds, m_nDate)
	
	// Perform the evaluation based on the data in this object
	// We need to set up the result object and call the Evaluate method on the product type interface.
	CComPtr<IDispatch>				spObject;		
			
	// I have made a rather ostensible decision here concerning object creation. I am aware
	// that there is a performance hit in creating a new product type object every time we
	// need one. I think that expecting the quant to program in such a way that allows you
	// to simply change various object properties - thereby allowing us to keep an object
	// in memory - and then evaluate it correctly is asking too much. The quant, may, for
	// example place code in FinalConstruct() which will obviously only be called once if
	// the object is cached.	
	
	GetObject(&spObject);
	
	CComDispatchDriverEx			ddObject(spObject);
	
	CResult::CreateInstance(&spResult);
	if (spResult->SetCalculate(Calculate)) propagate_error;
	
	// ensure that all correlations between underlying assets are defined
	CComPtr<IAssets>				spUnderlyings;
	
	GetBaseUnderlyings(spObject, spUnderlyings);		
	CAssets::InsertCorrelations(spUnderlyings );
	
	// now find the evaluate method on the product and call it
	DISPID				dispid;		
	if (ddObject.GetIDOfName(L"Evaluate", &dispid)) throw "The product '" + m_szProductType + "' does not have an evaluate function";
	if (ddObject.Invoke1(dispid, &CComVariant(spResult))){
		propagate_error;
	}		
}

void MlEqProduct::GetBaseUnderlyings(CComPtr<IDispatch> spObject, CComPtr<IAssets>& spAssets) const
{	
	pin_date(m_ds, m_nDate);
	CComPtr<IAssets>					spUnderlyings;
	
	GetUnderlyings(spObject, spUnderlyings);
	spUnderlyings->GetBaseUnderlyings(&spAssets);
}

DataSourceEnum MlEqProduct::GetDataSource(void) const
{
	return m_ds;
}

long MlEqProduct::GetDate(void) const
{
	return m_nDate;
}

//	Creates and sets up the parameters of the underlying object assoicated with the product.
void MlEqProduct::GetObject(IDispatch** pVal) const
{		
	pin_date(m_ds, m_nDate);
	HRESULT								hr;
	MlEqDictionary<CComVariant>			dict;
	std::vector<CParameterMap>			vpm;
	CComPtr<IDispatch>					spObject;						// this is the created product	
	
	Parameters()->GetDictionary(&dict);
	if (dict.GetValue(&vpm)) propagate_error;
	if (hr = g_pApplication->GetObjectManager().CreateObject(m_szProductType, spObject)) propagate_error_ex(hr);
	if (g_pApplication->GetObjectManager().ImplementPutValue(vpm, spObject)) propagate_error;
	spObject.CopyTo(pVal);
}

CComPtr<IParameters> MlEqProduct::GetParameters(void) const
{
	return m_spParameters;
}

//	returns a collection of assets on which the product depends
HRESULT MlEqProduct::GetUnderlyings(CComPtr<IDispatch> spObject, CComPtr<IAssets>& spAssets) const
{
	pin_date(m_ds, m_nDate);
	
	// The strategy here is to create an assets collection and loop through all read properties
	// associated with the product's object. If the property type is an Asset then we add this
	// to the collection.If it's a collection of assets then we add each to the collection. 
	// Else we ignore that property.

	HRESULT								hr;	
	
	if (hr = spAssets.CoCreateInstance(CLSID_Assets)) return hr;	

	CAssets*							pAssets = dynamic_cast<CAssets*>(spAssets.p);	
	CComPtr<ITypeInfo>					pti;
	TYPEATTR*							pTypeAttr;
	FUNCDESC*							pFuncDesc;
	
	if (hr = spObject->GetTypeInfo(0, LOCALE_USER_DEFAULT, &pti)) return hr;
	if (pti->GetTypeAttr(&pTypeAttr)) return hr;	
	for (unsigned int nIndex = 0; nIndex < pTypeAttr->cFuncs && !hr; nIndex++){
		if (!(pti->GetFuncDesc(nIndex, &pFuncDesc))){
			if (pFuncDesc->invkind & DISPATCH_PROPERTYGET){
				// we expect one level of indirection for a property get on either asset or assets
				VARTYPE vt = pFuncDesc->elemdescFunc.tdesc.vt;
				if (vt == VT_PTR){
					// one indirection level
					ELEMDESC* ped = (ELEMDESC*)pFuncDesc->elemdescFunc.tdesc.lptdesc;
					if (ped->tdesc.vt == VT_USERDEFINED){
						CComPtr<ITypeInfo>	pti_ref;
						if (!(hr = pti->GetRefTypeInfo(ped->tdesc.hreftype, &pti_ref))){
							TYPEATTR*	pTypeAttrRef;
							if (!(hr = pti_ref->GetTypeAttr(&pTypeAttrRef))){
								GUID guid = pTypeAttrRef->guid;
								// consider the asset and assets case
								if (guid == IID_IAsset){										
									// get the property
									CComVariant vObject;
									if (!(hr = CComDispatchDriverEx(spObject).GetProperty(pFuncDesc->memid, &vObject))){
										// add the object to the collection																				
										if (vObject.pdispVal) hr = spAssets->Add(CComVariant(), dynamic_cast<IAsset*>(vObject.pdispVal));
									}										
								} else if (guid == IID_IAssets){
									// we need to get the property - which is a collection - and add this collection to spAssets
									CComVariant vCollection;
									if (!(hr = CComDispatchDriverEx(spObject).GetProperty(pFuncDesc->memid, &vCollection))){
										CComPtr<IAssets> spAssetsCollection = dynamic_cast<IAssets*>(vCollection.pdispVal);
										try {
											pAssets->MergeCollection(spAssetsCollection);
										} catch (...){
											hr = E_FAIL;
										}
									}
								}								
								pti_ref->ReleaseTypeAttr(pTypeAttrRef);
							}									
						}
					}					
				}
			}
			pti->ReleaseFuncDesc(pFuncDesc);
		}	
	}	
	pti->ReleaseTypeAttr(pTypeAttr);
	return hr;	
}

CParameters* MlEqProduct::Parameters(void) const
{
	return dynamic_cast<CParameters*>(m_spParameters.p);
}

estring MlEqProduct::GetProductType(void) const
{
	return m_szProductType;
}

void MlEqProduct::PutDate(long nDate)
{	
	m_nDate = nDate;
}

void MlEqProduct::PutDataSource(DataSourceEnum ds)
{	
	m_ds = ds;	
}

void MlEqProduct::PutParameters(CComPtr<IParameters> spParameters)
{
	m_spParameters = spParameters;
}

void MlEqProduct::PutProductType(const std::string& szProductType)
{
	m_szProductType.assign(szProductType);
}
