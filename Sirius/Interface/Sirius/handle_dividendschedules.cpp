//	handle_dividendschedules.cpp : Implementation of CDividendSchedules
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "siriusapplication.h"
#include "handle_dividendschedules.h"
#include "handle_dividendschedule.h"
#include "MlEqZeroCurve.h"
#include "MlEqDate.h"

/*static*/ HRESULT CDividendSchedules::Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<IDividendSchedules>& spDividendSchedules)
{
	HRESULT								hr;	
	FileSystemEnum						fs;
	
	if (szDummy.size()) throw "The identifier value '" + szDummy + "' is invalid for the load dividend schedules request";
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();	
	if (hr = spDividendSchedules.CoCreateInstance(CLSID_DividendSchedules)) return hr;
	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		std::stringstream				ssQuery;
		std::string						szError = "No dividend schedules found on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		ssQuery << "sp_user_get_dividendschedules " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
		return CComObjectCollectionFunctions<IDividendSchedule>(&szError).ImplementLoadCollection_SQL(ssQuery, ds, spDividendSchedules);
	} else if (fs == fsNTFS){										
		return CComObjectCollectionFunctions<IDividendSchedule>().ImplementLoadCollection_NTFS(spDividendSchedules, 0, "", _Module.GetLocation(), "", ds, nDate, CDividendSchedule::Load);
	}
	return E_FAIL;
}

STDMETHODIMP CDividendSchedules::Reset()
{
	begin_function
	CComObjectCollectionFunctions<IDividendSchedule>(&m_coll).Reset();
	end_function
}

STDMETHODIMP CDividendSchedules::Shift(double Yield)
{
	begin_function
	CComObjectCollectionFunctions<IDividendSchedule>(&m_coll).Shift(Yield);
	end_function
}

STDMETHODIMP CDividendSchedules::Stick(void)
{
	begin_function
	CComObjectCollectionFunctions<IDividendSchedule>(&m_coll).Stick();
	end_function
}