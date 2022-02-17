//	ComObjectCollectionSerialisablekey.h : This is used to key serialisable objects by Name, Date and DataSource
//
//  Author :							   David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _COMOBJECTCOLLECTIONSERIALISABLEKEY_H
#define _COMOBJECTCOLLECTIONSERIALISABLEKEY_H

#include "comobjectcollectionserialisabledefaulter.h"

class CComObjectCollectionSerialisableKey
{
public:	
	template<class I> explicit CComObjectCollectionSerialisableKey(CComPtr<I> spObject) : m_bObjectNameKnown(false)
	{		
		CComBSTR							sName;		
		DataSourceEnum						ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
		DATE								date = CComObjectCollectionSerialisableDefaulter::GetDate();
		
		spObject->get_Name(&sName);
		spObject->get_DataSource(&ds);
		spObject->get_Date(&date);						
		Set(estring(sName), date, ds);
		CParameterMap::GetObjectName(spObject, &m_szObjectName);
		m_bObjectNameKnown = m_szObjectName.size() ? true : false;
	}	
	explicit CComObjectCollectionSerialisableKey(IDispatch* pObject);	// Using the templatised constructor is preferable.
	explicit CComObjectCollectionSerialisableKey(const CComBSTR& sName, long nDate, DataSourceEnum ds);
	explicit CComObjectCollectionSerialisableKey(const std::string& szName, long nDate, DataSourceEnum ds);
	explicit CComObjectCollectionSerialisableKey(const CComBSTR& sAssetName, CurrencyEnum ceComposite, CurrencyEnum cePay, long nDate, DataSourceEnum ds);
	operator CComBSTR(void) const;	
	operator CComVariant(void) const;
	static bool Respell(const std::string& szObjectName, const std::string& szIndex, std::string* pszKeyOut, std::string* pszShortKeyOut, std::string* pszNameOut, long* pnDateOut, DataSourceEnum* pdsOut, bool* pbDateDefaulted, bool* pbDataSourceDefaulted);
	operator std::string(void) const;
	estring				GetKeyAndObjectName(bool bNullifyDataSource, bool bExplicitDate) const;

protected:
	std::string			Get(bool bNullifyDataSource, bool bExplicitDate) const;
	void				Set(const std::string& szName, long nDate, DataSourceEnum ds);
	
	std::string			m_szName;
	long				m_nDate;
	DataSourceEnum		m_ds;
	std::string			m_szDate;			// this is a function of m_nDate but is included here for speed
	std::string			m_szDataSource;		// this is a function of m_ds included to prevent many string copy operations and calls to CEnumMap::GetString
	std::string			m_szObjectName;		// if the spObject constructor is called then we can set this
	bool				m_bObjectNameKnown;	// true if m_szObjectName is valid
};

#endif
