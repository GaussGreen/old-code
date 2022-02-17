//	xmlstreamer.h: interface for the CXmlStreamer class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _XMLSTREAMER_H 
#define _XMLSTREAMER_H

#pragma once

#include "FastStream.h"
class CSiriusApplication;
class CParameterMap;


//////////////////////////////////////////////////////////////////////
//	CXmlStreamer
//
//	Implements VARIANT <-> Xml
//
class CXmlStreamer
{
protected:
	class xml_attribute
	{
	public:
		const char*							sz_name_start;
		const char*							sz_name_end;
		const char*							sz_value_start;
		const char*							sz_value_end;

		xml_attribute()
		{
			clear();
		}
		void clear(void)
		{
			sz_name_start = sz_name_end = sz_value_start = sz_value_end = NULL;
		}
	};

	class xml_attributes : public std::vector<xml_attribute>		
	{
	public:
		HRESULT								GetValue(const char* szName, std::string* pszValue);
		HRESULT								GetValue(const char* szName, long* pnValue);
		HRESULT								GetValue(const char* szName, unsigned short* pnValue);
	};
	
public:
	CXmlStreamer();
	virtual									~CXmlStreamer();
	
	template<class ISerialisable> static HRESULT GetObject(char* szIn, DataSourceEnum ds, CComPtr<ISerialisable>& spObject)
	{
		HRESULT								hr;
		CComVariant							v;
				
		if (hr = GetVariant(szIn, &ds, v)) return hr;				
		ATLASSERT(v.vt == VT_DISPATCH);
		spObject = dynamic_cast<ISerialisable*>(v.pdispVal);
		ATLASSERT(spObject);
		// The XML data is invariant w.r.t. the data source so we have to manually set the data source.
		spObject->put_DataSource(ds);
		return S_OK;
	}
	static HRESULT							GetVariant(const char* szIn, CComVariant& vOut);
	
	template<class ISerialisable> static HRESULT GetXML(CComPtr<ISerialisable> spObject, xmlstreamer& ssOut)
	{
		// We have to make the XML invariant w.r.t. the data source or the PL process breaks.
		DataSourceEnum						ds;
		HRESULT								hr;
		if (spObject->get_DataSource(&ds)) propagate_error;
		if (spObject->put_DataSource(NoDataSource)) propagate_error;
		try {			
			hr = GetXML(CComVariant(spObject.p), ssOut, true, 0, NULL);
		} catch (...){
			hr = E_FAIL;
		}
		if (spObject->put_DataSource(ds)) propagate_error;
		return hr;
	}
	static HRESULT							GetXML(const CComVariant& vIn, xmlstreamer& ssOut);

protected:	
	static HRESULT							GetVariant(const char* szIn, DataSourceEnum* pds, CComVariant& vOut);	// Don't ever make this public.
	static HRESULT							GetVariantProcessFreeText(const char* szTextStart, const char* szTextEnd, std::vector<std::string>& vector_tag_stack, std::vector<void*>& vector_object_stack, DataSourceEnum* pds);
	static HRESULT							GetVariantProcessTag(const char* szTagStart, const char* szTagEnd, xml_attributes& attributes, std::vector<std::string>& vector_tag_stack, std::vector<void*>& vector_object_stack);
	static HRESULT							GetXML(const CComVariant& vIn, xmlstreamer& ssOut, bool bNullifyDataSource, long nRecursionCount, CComPtr<IDispatch> spParent);	// Don't ever make this public.	
};

#endif