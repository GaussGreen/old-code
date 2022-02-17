//	handle_array.cpp : Implementation of CArray
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_array.h"


STDMETHODIMP CArray::get_Element(long Index, double *pVal)
{
	begin_function
	if (Index < 1 || Index > m_h->getsize()) throw CStringResource(IDS_ARRAY_INDEX);
	*pVal = (*m_h)[Index - 1];
	end_function
}

STDMETHODIMP CArray::get_Size(long *pVal)
{
	begin_function
	*pVal = m_h->getsize();
	end_function
}

STDMETHODIMP CArray::get_Value(VARIANT *pVal)
{	
	begin_function
	unmap_parameter(*m_h, pm)
	return pm.GetValue(pVal);
	end_function
}

STDMETHODIMP CArray::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IArray };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CArray::put_Element(long Index, double newVal)
{
	begin_function
	if (Index < 1 || Index > m_h->getsize()) throw CStringResource(IDS_ARRAY_INDEX);
	(*m_h)[Index - 1] = newVal;
	end_function
}

STDMETHODIMP CArray::put_Size(long newVal)
{
	begin_function
	m_h->PutSize(newVal);
	end_function
}

STDMETHODIMP CArray::put_Value(VARIANT newVal)
{
	begin_function
	map_optional_parameter(newVal, CVector, vector, CVector())		// a vector with no elements in it is OK
	*m_h = vector;	
	end_function
}