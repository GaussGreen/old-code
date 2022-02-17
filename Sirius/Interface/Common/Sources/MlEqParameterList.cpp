//	MlEqParameterList.cpp :	 Implementation of MlEqParameterList
//							 THIS SHOULD NOT BE IN THE ANALYTICS LIBRARY
//                           SINCE IT CONTAINS VARIANTS ETC.
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqParameterList.h"

//  Attempts to return a long value associated with the parameter name
void MlEqParameterList::GetValue(const std::string& szName, long* pn) const
{
	CComVariant							v;

	MlEqDictionary<CComVariant>::GetValue(szName, &v, true);
	if (v.ChangeType(VT_I4)){
		throw "Invalid parameter for parameter '" + szName + "'";
	}
	*pn = v.lVal;
}

//	Attempts to return a double value associated with the parameter name.
void MlEqParameterList::GetValue(const std::string& szName, double* pf) const
{
	CComVariant							v;

	MlEqDictionary<CComVariant>::GetValue(szName, &v, true);
	if (v.ChangeType(VT_R8)){
		throw "Invalid value for parameter '" + szName + "'";
	}
	*pf = v.dblVal;
}
//	Attempts to return a double value associated with the parameter name.
//	We support a default value.
void MlEqParameterList::GetValue(const std::string& szName, double* pf, double fDefault) const
{
	CComVariant							v;

	try {
		MlEqDictionary<CComVariant>::GetValue(szName, &v, true);
	} catch (...){
		*pf = fDefault;
	}
	if (v.ChangeType(VT_R8)){
		throw "Invalid value for parameter '" + szName + "'";
	}
	*pf = v.dblVal;
}

//	Attempts to set an (already allocated) vector of doubles in the order
//	specified by an input delimited string.
void MlEqParameterList::GetValue(const std::string& szNames, const std::string& szDelimit, double* af) const
{		
	std::vector<std::string>			vector;						// vector of names
	
	estring::Split(szNames, szDelimit, &vector);
	for (long nElement = 0; nElement < vector.size(); nElement++){
		GetValue(vector[nElement], af + nElement);
	}
}		
//	Attempts to set an (already allocated) vector of doubles in the order
//	specified by an input delimited string. We support default values
//  by use of a second delimiter which separates the variable name from
//  the stringified version of the default value. I know this not the most
//  computationally efficient way of doing things but bear in mind that this
//  simplifies the work of the quant and the performance hit is neglibable
//  when you consider the time taken to set up and run the quant evaluation.
void MlEqParameterList::GetValue(const std::string& szNames, const std::string& szDelimit, const std::string& szDefaultDelimit, double* af) const
{
	std::vector<std::string>			vector;						// vector of names (may also contain default values)

	estring::Split(szNames, szDelimit, &vector);
	for (long nElement = 0; nElement < vector.size(); nElement++){
		if (vector[nElement].find(szDefaultDelimit) == std::string::npos){
			// no default value
			GetValue(vector[nElement], af + nElement);
		} else {
			// default value
			std::string szName = estring::left(&vector[nElement], "=");
			estring::trim(&szName);
			estring szDefaultValue = estring::mid(&vector[nElement], "=");
			szDefaultValue.trim();
			double fDefault;
			if (!szDefaultValue.isdouble(&fDefault)) throw "'" + szDefaultValue + "' is an invalid default value for parameter list element '" + szName + "'";
			GetValue(szName, af + nElement, fDefault);
		}
	}
}