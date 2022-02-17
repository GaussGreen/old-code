//	handle_matrix.cpp : Implementation of CInterfaceMatrix
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Sirius.h"
#include "handle_matrix.h"


STDMETHODIMP CInterfaceMatrix::get_Value(VARIANT *pVal)
{	
	begin_function
	unmap_parameter(*m_h, pm)
	return pm.GetValue(pVal);
	end_function
}

STDMETHODIMP CInterfaceMatrix::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IMatrix };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CInterfaceMatrix::put_Value(VARIANT newVal)
{
	begin_function
	map_optional_parameter(newVal, CMatrix, matrix, CMatrix())		// a matrix with no elements in it is OK
	*m_h = matrix;
	end_function
}
