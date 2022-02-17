//	comobjectcollectionserialisabledefaulter.cpp
//
//	Author : David Cuin
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionserialisabledefaulter.h"
#include "MlEqDate.h"

/*static*/ std::stack<CComObjectCollectionSerialisableDefaulter::key>	CComObjectCollectionSerialisableDefaulter::s_keys;

/*static*/ const std::string& CComObjectCollectionSerialisableDefaulter::GetDataSourceStr(void)
{
	if (!s_keys.size()) return _Module.GetDefaultDataSourceStr();
	return s_keys.top().m_szDataSource;
}
/*static*/ DataSourceEnum CComObjectCollectionSerialisableDefaulter::GetDataSource(void)
{
	if (!s_keys.size() || !s_keys.top().m_bDataSourceValid) return _Module.GetDefaultDataSource();
	return s_keys.top().m_ds != NoDataSource ? s_keys.top().m_ds : _Module.GetDefaultDataSource();
}
/*static*/ long CComObjectCollectionSerialisableDefaulter::GetDate(void)
{
	if (!s_keys.size()) return MlEqDate::GetCurrentDate();	
	return s_keys.top().m_nDate ? s_keys.top().m_nDate : MlEqDate::GetCurrentDate();
}
/*explicit*/ CComObjectCollectionSerialisableDefaulter::CComObjectCollectionSerialisableDefaulter(const std::string& szDataSource, long nDate)
{
	key									k;
	DataSourceEnum						ds;
	
	if (!nDate && s_keys.size()){
		nDate = s_keys.top().m_nDate;	// Don't use GetDate function since we don't want any defaulting here.
	}
	if (!szDataSource.size() && s_keys.size()){
		k.m_szDataSource = s_keys.top().m_szDataSource;
	} else {
		k.m_szDataSource = szDataSource;
	}
			
	if (CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, k.m_szDataSource, &ds)){
		// failure
		k.m_ds = NoDataSource;
		k.m_bDataSourceValid = false;
	} else {
		k.m_ds = ds;
		k.m_bDataSourceValid = true;		
	}	
	k.m_nDate = nDate;
	s_keys.push(k);
}
/*explicit*/ CComObjectCollectionSerialisableDefaulter::CComObjectCollectionSerialisableDefaulter(DataSourceEnum ds, long nDate)
{
	if (!nDate && s_keys.size()) nDate = s_keys.top().m_nDate;			// Don't use GetDate function since we don't want any defaulting here.
	if (ds == NoDataSource && s_keys.size()) ds = s_keys.top().m_ds;
			
	key k;
	k.m_ds = ds;
	k.m_szDataSource = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds);
	k.m_nDate = nDate;	
	k.m_bDataSourceValid = true;
	s_keys.push(k);
}
/*virtual*/ CComObjectCollectionSerialisableDefaulter::~CComObjectCollectionSerialisableDefaulter(void)
{
	ATLASSERT(s_keys.size());
	s_keys.pop();
}