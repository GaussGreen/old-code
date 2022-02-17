//	handle_results.cpp : Implementation of CResults
//
//  Author :			 David Cuin
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_results.h"
#include "handle_result.h"

// Returns a 3 column list consisting of Field Name, Union Data Type, Field Size.
// I anticipate this function being used when building some kind of serialised
// representation of a results set.
STDMETHODIMP CResults::GetFields(VARIANT* pVal)
{
	begin_function

	std::map<std::string, VARTYPE>		mapNameToVarType;
	std::map<std::string, std::string>	mapLowerToProperCase;			// Lower case to proper case map
	std::map<std::string, long>			mapNameToSize;
	
	for (ContainerType::const_iterator it = m_coll.begin(); it != m_coll.end(); it++){
		CResult* pResult = dynamic_cast<CResult*>(it->second.pdispVal);		
		ATLASSERT(pResult);
		for (std::map<std::string, double>::const_iterator it = pResult->GetStringToDoubleMap().begin(); it != pResult->GetStringToDoubleMap().end(); it++){
			std::string szProperCase = pResult->GetLowerCaseToProperCaseMap().find(it->first)->second;
			GetFieldsAddRecord(it->first, szProperCase, &mapNameToVarType, &mapLowerToProperCase, &mapNameToSize, VT_R8, sizeof(double));
		}
		for (std::map<std::string, CComVariant>::const_iterator it = pResult->GetStringToVariantMap().begin(); it != pResult->GetStringToVariantMap().end(); it++){
			std::string szProperCase = pResult->GetLowerCaseToProperCaseMap().find(it->first)->second;									
			GetFieldsAddRecord(it->first, szProperCase, &mapNameToVarType, &mapLowerToProperCase, &mapNameToSize, it->second.vt, CParameterMap::GetVariantSize(it->second));
		}
	}
	
	CParameterMap						pmRet;
	long								nRow(0L);
	pmRet.SetSize(mapNameToVarType.size(), 3);	
	for (std::map<std::string, VARTYPE>::const_iterator itNameToVarType = mapNameToVarType.begin(); itNameToVarType != mapNameToVarType.end(); itNameToVarType++){
		pmRet.SetValue(nRow, 0, itNameToVarType->first);
		pmRet.SetValue(nRow, 1, itNameToVarType->second);		
		pmRet.SetValue(nRow, 2, mapNameToSize[itNameToVarType->first]);
		++nRow;
	}
	return pmRet.GetValue(pVal);
	end_function
}


void CResults::GetFieldsAddRecord(const std::string& szLowerCase, const std::string& szProperCase, std::map<std::string, VARTYPE>* pmapNameToVarType, std::map<std::string, std::string>* pmapLowerToProperCase, std::map<std::string, long>* pmapNameToSize, VARTYPE vt, long nSize)
//	vt - type of field under consideration.
//	nSize - size of the current record (bytes).
{
	std::map<std::string, std::string>::iterator itLowerToProperCase = pmapLowerToProperCase->find(szLowerCase);
	if (itLowerToProperCase == pmapLowerToProperCase->end()){
		// Insert new				
		(*pmapNameToVarType)[szProperCase] = vt;
		(*pmapLowerToProperCase)[szLowerCase] = szProperCase;
		(*pmapNameToSize)[szProperCase] = nSize;
	} else {
		// Update existing, taking care of a field that has changed from a number to a string				
		VARTYPE vtCurrent = (*pmapNameToVarType)[szProperCase];
		if (vt == VT_BSTR && CParameterMap::IsVariantTypeNumeric(vtCurrent)){
			nSize = 32;	// reasonable length for a string representation of a number
		}
		nSize = std::max(nSize, (*pmapNameToSize)[szProperCase]);
		(*pmapNameToSize)[szProperCase] = nSize;
		(*pmapNameToVarType)[szProperCase] = CParameterMap::GetUnionVariantType(vtCurrent, vt);
	}
}