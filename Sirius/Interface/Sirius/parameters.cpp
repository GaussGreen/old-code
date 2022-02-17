//	parameters.cpp : Implementation of CParameters
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "parameters.h"
#include "MlEqDictionary.h"


void CParameters::GetDictionary(MlEqDictionary<CComVariant>* pdict)
{		
	long								nCount = 0L;
		
	get_Count(&nCount);
	pdict->Clear();			
	for (long n = 1; n <= nCount; n++){
		CComPtr<IParameter>				spParameter;
		CComBSTR						sName;
		CComVariant						vValue;

		if (get_Item(CComVariant(n), &spParameter)) continue;
		if (spParameter->get_Name(&sName)) continue;
		if (spParameter->get_Value(&vValue)) continue;
		pdict->AddValue(estring(sName), vValue);
	}
}





bool CParameters::IsInCollection(const CComBSTR& s)
{	
	return m_coll.find(s) != m_coll.end() ? true : false;
}