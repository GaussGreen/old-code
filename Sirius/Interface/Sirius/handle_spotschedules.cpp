//	handle_spotschedules.cpp : Implementation of CSpotSchedules
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_spotschedules.h"
#include "handle_spotschedule.h"
#include "comobjectcollectionserialisabledefaulter.h"

long CSpotSchedules::GetDate(void) const
{
	if (_Module.GetFileSystem() == fsNTFS){
		return 0;	// spot schedules are not available historically
	} else {
		return CComObjectCollectionSerialisableDefaulter::GetDate();
	}
}

/*static*/ HRESULT CSpotSchedules::Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<ISpotSchedules>& spSpotSchedules)
{
	HRESULT								hr;	
	FileSystemEnum						fs;
	
	if (szDummy.size()) throw "The identifier value '" + szDummy + "' is invalid for the load spot schedules request";
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (hr = spSpotSchedules.CoCreateInstance(CLSID_SpotSchedules)) return hr;
	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		throw "Loading spot schedules from a SQL Server database is not supported";
	} else if (fs == fsNTFS){										
		return CComObjectCollectionFunctions<ISpotSchedule>().ImplementLoadCollection_NTFS(spSpotSchedules, 0, "", _Module.GetLocation(), "", ds, nDate, CSpotSchedule::Load);
	}
	return E_FAIL;
}
