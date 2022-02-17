//	mleqdictionary.h:		Interface for the MlEqDictionary class.
//							This is a case-insensitive dictionary.
//
//							I also provide a mechanism (via the
//							BeginExhaust and CheckUnexhausted functions)
//							to test whether all the elements in a
//							dictionary have been used after a series
//							of GetValue calls.
//
//							This is not to be regarded as an analytics-level
//							class since it makes reference to BSTRs and
//							CParameterMap classes
//							and BSTR.
//
//	author:					David Cuin
//
//////////////////////////////////////////////////////////////////////

#ifndef _DICTIONARY_H
#define _DICTIONARY_H

#pragma once
#include "smart.h"
#include "parametermap.h"

template<class T>
class MlEqDictionary : public RCObject
{
private:		
	// we never want inherited classes to access variables
	// or the exhaust feature may unwittingly broken.
	std::map<std::string, T>					m_map;							// parameter list map
	std::map<std::string, std::string>			m_mapLowerCaseToProperCase;		// transforms the internal lower case representation of the input values to the originally specified string values
	mutable std::map<std::string, std::string>	m_mapExhaust;					// used to check parameter exhaustion (map is as for m_mapLowerCaseToProperCase)

public:
	void AddValue(const BSTR& sName, const T& Value)
	{
		AddValue(estring(sName), Value);
	}
	
	void AddValue(UINT nResourceID, const T& Value)
	{
		AddValue(CStringResource(nResourceID), Value);
	}
	
	void AddValue(const std::string& szName, const T& Value)
	{
		estring							szName_lc;						// lower case form of szName
		szName_lc.assign(szName);
		szName_lc.StripWhiteSpace();
		szName_lc.lc();		
		if (!szName.size()) throw "Invalid name";
		m_mapLowerCaseToProperCase[szName_lc] = szName;
		m_map[szName_lc] = Value;
	}
	
	void BeginExhaust(void) const
	{
		m_mapExhaust = m_mapLowerCaseToProperCase;
	}

	// Returns the number of unused elements. If a string pointer is passed
	// then it is blanked if there are no unused elements or set to the 
	// first unused non-blank element. If a vector of string pointer is passed then
	// it is also blanked if there are no unused elements or it is set to
	// the proper case representation of the unused elements.
	long CheckUnexhausted(std::string* psz, std::vector<std::string>* pasz) const
	{
		if (psz) psz->clear();
		if (pasz) pasz->clear();
		for (std::map<std::string, std::string>::const_iterator it = m_mapExhaust.begin(); it != m_mapExhaust.end(); it++){
			if (psz && !psz->size()) psz->assign(it->second);
			if (pasz) pasz->push_back(it->second);
		}
		return m_mapExhaust.size();
	}
	// As above but we throw a string exception if there are unused elements
	void CheckUnexhausted(void) const
	{
		std::string						sz;
		if (CheckUnexhausted(&sz, NULL)) throw "Unused parameter '" + sz + "' in the parameter list";
	}

	void Clear()
	{
		m_map.clear();
		m_mapLowerCaseToProperCase.clear();
	}
	
	HRESULT CopyToClipboard(void) const
	{
		CParameterMap						pm;
		HRESULT								hr;

		if (hr = GetValue(&pm)) return hr;
		return pm.CopyToClipboard();
	}

	void Delete(UINT nResourceID, bool bThrow)
	{
		Delete(CStringResource(nResourceID), bThrow);
	}

	void Delete(const std::string& szName, bool bThrow)
	{
		estring											szName_lc;						// lower case form of szName
		std::map<std::string, T>::iterator				it; 
		std::map<std::string, std::string>::iterator	itCase;
		szName_lc.assign(szName);
		szName_lc.StripWhiteSpace();
		szName_lc.lc();
		if (!szName.size()){
			if (bThrow) throw "Invalid name";
			return;
		}
		it = m_map.find(szName_lc);
		itCase = m_mapLowerCaseToProperCase.find(szName_lc);
		if (it == m_map.end() || itCase == m_mapLowerCaseToProperCase.end()){
			if (bThrow) throw "'" + szName + "' does not appear in the list";
			return;
		}
		m_map.erase(it);
		m_mapLowerCaseToProperCase.erase(itCase);
	}

	long GetCount(void) const
	{
		return m_map.size();
	}	

	const std::map<std::string, T>&				GetMap(void) const {return m_map;}
	const std::map<std::string, std::string>&	GetLCMap(void) const {return m_mapLowerCaseToProperCase;}
	
	// Exhaustive GetValue functions. These return the entire dictionary
	// in various forms. Notice that these functions return HRESULT.
	HRESULT GetValue(CParameterMap* ppm) const
	{
		ppm->Clear();
		if (m_map.size()){
			ppm->SetSize(m_map.size(), 2);
			long							nRow(0);	
			HRESULT							hr;
			for (std::map<std::string, T>::const_iterator itCount = m_map.begin(); itCount != m_map.end(); ++itCount){
				std::map<std::string, std::string>::const_iterator itString = m_mapLowerCaseToProperCase.find(itCount->first);
				ATLASSERT(itString != m_mapLowerCaseToProperCase.end());								
				if (hr = ppm->SetValue(nRow, 0, itString->second)) return hr;
				if (hr = ppm->SetValue(nRow, 1, itCount->second)) return hr;
				nRow++;		
			}				
		}
		return S_OK;
	}	
	// sets two parameter maps in a vector - the first one holds the dictionary headers, the second holds the dictionary values
	HRESULT GetValue(std::vector<CParameterMap>* pvpm) const
	{		
		pvpm->clear();
		pvpm->resize(2);
		if (m_map.size()){
			(*pvpm)[0].SetSize(m_map.size(), 1);	// header
			(*pvpm)[1].SetSize(m_map.size(), 1);	// values
			long							nRow(0);
			HRESULT							hr;
			for (std::map<std::string, T>::const_iterator itCount = m_map.begin(); itCount != m_map.end(); ++itCount){
				std::map<std::string, std::string>::const_iterator itString = m_mapLowerCaseToProperCase.find(itCount->first);
				ATLASSERT(itString != m_mapLowerCaseToProperCase.end());								
				if (hr = (*pvpm)[0].SetValue(nRow, 0, itString->second)) return hr;
				if (hr = (*pvpm)[1].SetValue(nRow, 0, itCount->second)) return hr;
				nRow++;		
			}				
		}
		return S_OK;
	}
	// sets two parameter maps - the first one holds the dictionary headers, the second holds the dictionary values
	HRESULT GetValue(CParameterMap* ppmHeaders, CParameterMap* ppmValues)
	{
		HRESULT							hr;
		std::vector<CParameterMap>		vpm;

		if (hr = GetValue(&vpm)) return hr;
		if (hr = ppmHeaders->Attach(&vpm[0])) return hr;
		if (hr = ppmValues->Attach(&vpm[1])) return hr;
		return S_OK;
	}
	// returns a flat form of the dictionary (therefore compatible with MLGetValue calls into spreadsheets)
	HRESULT	GetValue(CComVariant* pValue) const
	{
		CParameterMap						pm;
		HRESULT								hr;

		if (hr = GetValue(&pm)) return hr;
		return pm.GetValue(pValue);
	}
	void GetValue(VARIANT* pValue) const
	{
		GetValue((CComVariant*)pValue);
	}
	
	// Non-exhaustive GetValue functions. These return a subset of the dictionary
	// (probably only one element thereof). Notice that these functions return void.	
	void GetValue(UINT nResourceID, T* Value, bool bThrow) const
	{
		GetValue(CStringResource(nResourceID), Value, bThrow);
	}

	void GetValue(const BSTR& sName, T* pValue, bool bThrow) const
	{
		GetValue(estring(sName), pValue, bThrow);
	}
	
	// This is the lowest-level GetValue function. All non-exhastive GetValue
	// functions must call into this one.
	void GetValue(const std::string& szName, T* pValue, bool bThrow) const
	{
		estring								szName_lc;					// lower case form of szName
		szName_lc.assign(szName);
		szName_lc.StripWhiteSpace();
		szName_lc.lc();		
		std::map<std::string, T>::const_iterator it = m_map.find(szName_lc);
		if (it == m_map.end()){
			if (bThrow) throw "Cannot find element '" + szName + "'";
			return;
		}
		
		// implement exhaustion
		std::map<std::string, std::string>::iterator itExhaust = m_mapExhaust.find(szName_lc);
		if (itExhaust != m_mapExhaust.end()) m_mapExhaust.erase(itExhaust);
	
		// retun the value
		*pValue = it->second;
	}
	
	//	Returns true if the dictionary contains the supplied key
	bool IsDefined(const std::string& szName) const
	{
		estring								szName_lc;					// lower case form of szName
		szName_lc.assign(szName);
		szName_lc.StripWhiteSpace();
		szName_lc.lc();
		std::map<std::string, T>::const_iterator it = m_map.find(szName_lc);
		return (it != m_map.end());
	}

	void PutValue(const CParameterMap& pm)
	{
		estring								szName;
		estring								szName_lc;
		T									Element;
		std::map<std::string, T>			map;							// temporary for m_map
		std::map<std::string, std::string>	mapLowerCaseToProperCase;		// temporary for m_mapLowerCaseToProperCase
		
		if (pm.IsBlank()){
			m_map.clear();
			m_mapLowerCaseToProperCase.clear();
			return;
		}

		if (pm.GetCols() != 2) throw CStringResource(IDS_COLUMNS_INVALID);
		if (!pm.GetRows()) throw CStringResource(IDS_ROWS_INVALID);
		
		for (long nRow = 0; nRow < pm.GetRows(); nRow++){
			if (pm.GetValue(nRow, 0, &szName)) throw CStringResource(IDS_INVALID_NAME);
			szName_lc.assign(szName);
			szName_lc.StripWhiteSpace();
			szName_lc.lc();
			if (!szName_lc.size()){				
				if (pm.IsBlank(nRow)) continue;		// don't return an error if the entire row is blank
				throw CStringResource(IDS_INVALID_NAME);
			}
			if (mapLowerCaseToProperCase.find(szName_lc) != mapLowerCaseToProperCase.end()){				
				throw "The name '" + szName_lc + "' is already in use";
			}
			if (pm.GetValue(nRow, 1, &Element)) throw CStringResource(IDS_INVALID_PARAMETER_VALUE);
			mapLowerCaseToProperCase[szName_lc] = szName;
			map[szName_lc] = Element;
		}

		m_map = map;
		m_mapLowerCaseToProperCase = mapLowerCaseToProperCase;
	}
	
	void PutValue(const CComVariant& Value)
	//	Value - Data to pass into the dictionary. Blank variants are supported.
	{
		CParameterMap						pm;
		if (pm.SetValue(Value)) throw "Unhandled exception in MlEqDictionary::PutValue";
		PutValue(pm);
	}
};

#endif
