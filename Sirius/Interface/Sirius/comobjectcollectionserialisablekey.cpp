//	comobjectcollectionserialisablekey.cpp: implementation of CComObjectCollectionDefaulterKey
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionserialisablekey.h"
#include "MlEqDate.h"


//////////////////////////////////////////////////////////////////////////////
//	CComObjectCollectionSerialisableKey
//
/*explicit*/ CComObjectCollectionSerialisableKey::CComObjectCollectionSerialisableKey(IDispatch* pObject) : m_bObjectNameKnown(false)
{		
	CComPtr<IDispatch>					spObject(pObject);
	CComVariant							vName, vDataSource, vDate;
	DataSourceEnum						ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	long								nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
	
	CComDispatchDriverEx(spObject).GetPropertyByName(L"Name", &vName);
	CComDispatchDriverEx(spObject).GetPropertyByName(L"DataSource", &vDataSource);
	CComDispatchDriverEx(spObject).GetPropertyByName(L"Date", &vDate);
	CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, vDataSource, &ds);
	if (!vDate.ChangeType(VT_I4)){
		nDate = vDate.lVal;
	}
	Set(estring(vName), nDate, ds);
	CParameterMap::GetObjectName(spObject, &m_szObjectName);
	m_bObjectNameKnown = m_szObjectName.size() ? true : false;
}
/*explicit*/ CComObjectCollectionSerialisableKey::CComObjectCollectionSerialisableKey(const CComBSTR& sName, long nDate, DataSourceEnum ds) : m_bObjectNameKnown(false)
{
	Set(estring(sName), nDate, ds);
}
/*explicit*/ CComObjectCollectionSerialisableKey::CComObjectCollectionSerialisableKey(const std::string& szName, long nDate, DataSourceEnum ds) : m_bObjectNameKnown(false)
{
	Set(szName, nDate, ds);
}
/*explicit*/ CComObjectCollectionSerialisableKey::CComObjectCollectionSerialisableKey(const CComBSTR& sAssetName, CurrencyEnum ceComposite, CurrencyEnum cePay, long nDate, DataSourceEnum ds)
{
	estring szName(estring(sAssetName) + " (" + CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, ceComposite) + ") > " + CEnumMap::GetString("CurrencyEnum", LIBID_Sirius, cePay));
	Set(szName, nDate, ds);
}
CComObjectCollectionSerialisableKey::operator CComBSTR(void) const
{	
	return estring::GetBSTR(Get(false, false));
}
CComObjectCollectionSerialisableKey::operator CComVariant(void) const
{
	return estring::GetValue(Get(false, false));
}
CComObjectCollectionSerialisableKey::operator std::string(void) const
{
	return Get(false, false);
}

std::string CComObjectCollectionSerialisableKey::Get(bool bNullifyDataSource, bool bExplicitDate) const
//  bNullifyDataSource - true if we set all object keys to have Name.Date.NoDataSource.
//	bExplicitDate - true if we return the current date if m_nDate is zero.
{
	if (bNullifyDataSource){
		if (m_nDate){
			return m_szName + "@" + m_szDate + ".NoDataSource";
		} else {
			if (bExplicitDate){
				return m_szName + "@" + MlEqDate(MlEqDate::GetCurrentDate()).GetString() + ".NoDataSource";
			} else {
				return m_szName + "@NoDataSource";
			}
		}
	} else {
		if (m_nDate){
			return m_szName + "@" + m_szDate + "." + m_szDataSource;
		} else {
			if (bExplicitDate){
				return m_szName + "@" + MlEqDate(MlEqDate::GetCurrentDate()).GetString() + "." + m_szDataSource;
			} else {
				return m_szName + "@" + m_szDataSource;
			}
		}
	}
}

estring CComObjectCollectionSerialisableKey::GetKeyAndObjectName(bool bNullifyDataSource, bool bExplicitDate) const
//  bNullifyDataSource - true if we set all object keys to have Name.Date.NoDataSource
//	bExplicitDate - true if we return the current date if m_nDate is zero.
{
	if (!m_bObjectNameKnown) throw "Error in CComObjectCollectionSerialisableKey::GetKeyAndObjectName(): Object name not known!";
	return m_szObjectName + "::" + Get(bNullifyDataSource, bExplicitDate);
}

/*static*/ bool CComObjectCollectionSerialisableKey::Respell(const std::string& szObjectName, const std::string& szIndex, std::string* pszKeyOut, std::string* pszShortKeyOut, std::string* pszNameOut, long* pnDateOut, DataSourceEnum* pdsOut, bool* pbDateDefaulted, bool* pbDataSourceDefaulted)
//  szObjectName - e.g. "Asset", "ZeroCurve"
//	szIndex - object key to respell
//  pszKeyOut, pszNameOut, pnDateOut, pdsOut - (returned, nullable)
//  pszShortKeyOut - this is the key returned when the output date and data source are not defaulted (nullable).
//	pbDateDefaulted - true if the output date has been defaulted (nullable)
//  pbDataSourceDefaulted - true if the output data source has been defaulted (nullable)
{
	estring								szName;
	long								nDate = 0L;
	DataSourceEnum						ds = NoDataSource;
	std::string::size_type				stName, stDate;	
						
	if ((stName = szIndex.find('@')) == szIndex.npos){
		szName.assign(szIndex);
	} else {
		// contins date and / or data source information
		szName.assign(szIndex.substr(0, stName));
		if ((stDate = szIndex.find('.', stName + 1)) == szIndex.npos){			
			estring szDateOrDataSource = szIndex.substr(stName + 1);
			if (CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szDateOrDataSource, &ds) && CParameterMap::StringToDate(szDateOrDataSource, &nDate)){
				// failed
				return false;
			}
		} else {
			// date and data source given
			std::string szL = szIndex.substr(stName + 1, stDate - stName - 1);
			std::string szR = szIndex.substr(stDate + 1);
			if (!CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szL, &ds) && !CParameterMap::StringToDate(szR, &nDate)){
				// sz is of the form Name@DataSource.Date
			} else if (!CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szR, &ds) && !CParameterMap::StringToDate(szL, &nDate)){
				// sz is of the form Name@Date.DataSource
			} else {
				return false;
			}
		}
	}

	// The correctly spelt name should never contain the name of the object.
	if (!szName.left("::").CompareNoCase(szObjectName)){
		szName.assign(szName.right("::"));
	}

	// Success if this point is reached.
	// Set the short key before defaulting the date and data source!
	if (pszShortKeyOut) pszShortKeyOut->assign(CComObjectCollectionSerialisableKey(szName, nDate, ds));	
	if (nDate){
		if (pbDateDefaulted) *pbDateDefaulted = false;
	} else {
		nDate = CComObjectCollectionSerialisableDefaulter::GetDate();
		if (pbDateDefaulted) *pbDateDefaulted = true;
	}
	if (ds != NoDataSource){
		if (pbDataSourceDefaulted) *pbDataSourceDefaulted = false;
	} else {
		ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
		if (pbDataSourceDefaulted) *pbDataSourceDefaulted = true;
	}

	if (pszKeyOut) pszKeyOut->assign(CComObjectCollectionSerialisableKey(szName, nDate, ds));
	if (pszNameOut) pszNameOut->assign(szName);
	if (pnDateOut) *pnDateOut = nDate;
	if (pdsOut) *pdsOut = ds;
	return true;
}

void CComObjectCollectionSerialisableKey::Set(const std::string& szName, long nDate, DataSourceEnum ds)
{
	m_szName.assign(szName);
	m_nDate = nDate;
	m_ds = ds;	
	m_szDate.assign(MlEqDate(nDate).GetString());			
	if (CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_ds, &m_szDataSource)) ATLASSERT(false);
}

