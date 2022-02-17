//	MlEqParameterList.h :	 Parameter list handle class
//							 THIS SHOULD NOT BE IN THE ANALYTICS LIBRARY
//                           SINCE IT CONTAINS VARIANTS ETC.
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQPARAMETERLIST_H_
#define _MLEQPARAMETERLIST_H_

#pragma once
#include "mleqdictionary.h"

/* Unused Code -> // ´returns an object, if any, associated with a dictionary elemenet
template<class ISingular> HRESULT GetValue(const std::string& szName, CComPtr<ISingular>& sp)
{
	HRESULT							hr;
	CComVariant						v;

	if (hr = MlEqDictionary<CComVariant, &IID_IParameterList>::GetValue(szName, &v)) return hr;
	return g_pApplication->GetObjectManager().GetObject(szName, sp);
}*/

class MlEqParameterList : public MlEqDictionary<CComVariant>
{
public:
	void								GetValue(const std::string& szName, long* pn) const;
	void								GetValue(const std::string& szName, double* pf) const;
	void								GetValue(const std::string& szName, double* pf, double fDefault) const;
	void								GetValue(const std::string& szNames, const std::string& szDelimit, double* af) const;
	void								GetValue(const std::string& szNames, const std::string& szDelimit, const std::string& szDefaultDelimit, double* af) const;
};

#endif