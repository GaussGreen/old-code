//	handle_volatilitystructures.cpp : Implementation of CVolatilityStructures
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "siriusapplication.h"
#include "handle_volatilitystructures.h"
#include "handle_volatilitystructure.h"
#include "MlEqDate.h"

/*static*/ HRESULT CVolatilityStructures::Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<IVolatilityStructures>& spVolatilityStructures)
{
	HRESULT								hr;	
	FileSystemEnum						fs;
	
	if (szDummy.size()) throw "The identifier value '" + szDummy + "' is invalid for the load volatility structures request";
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();	
	if (hr = spVolatilityStructures.CoCreateInstance(CLSID_VolatilityStructures)) return hr;
	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		std::stringstream				ssQuery;
		std::string						szError;
		
		if (!_Module.GetLocationFirstTryEnabled() || !_Module.GetLocationFirstTry().size() || !estring::CompareNoCase(_Module.GetLocationFirstTry(), _Module.GetLocation())){
			// not using trial location
			ssQuery << "sp_user_get_volatilitystructures " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
			szError = "No volatility structures found on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		} else {
			// are using trial location
			ssQuery << "sp_user_get_volatilitystructures " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocationFirstTry() << "', '" << _Module.GetLocation() << "'";
			szError = "No volatility structures found on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocationFirstTry() + "' or '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		}
		return CComObjectCollectionFunctions<IVolatilityStructure>(&szError).ImplementLoadCollection_SQL(ssQuery, ds, spVolatilityStructures);
	} else if (fs == fsNTFS){
		if (!_Module.GetLocationFirstTryEnabled() || !_Module.GetLocationFirstTry().size() || !estring::CompareNoCase(_Module.GetLocationFirstTry(), _Module.GetLocation())){
			// not using trial location
			return CComObjectCollectionFunctions<IVolatilityStructure>().ImplementLoadCollection_NTFS(spVolatilityStructures, 0, "", _Module.GetLocation(), "", ds, nDate, CVolatilityStructure::Load);
		} else {
			// are using trial location
			return CComObjectCollectionFunctions<IVolatilityStructure>().ImplementLoadCollection_NTFS(spVolatilityStructures, 0, "", _Module.GetLocationFirstTry(), _Module.GetLocation(), ds, nDate, CVolatilityStructure::Load);
		}
	}
	return E_FAIL;
}

STDMETHODIMP CVolatilityStructures::Reset()
{
	begin_function
	CComObjectCollectionFunctions<IVolatilityStructure>(&m_coll).Reset();
	end_function
}

STDMETHODIMP CVolatilityStructures::Shift(double Amount)
{
	begin_function
	CComObjectCollectionFunctions<IVolatilityStructure>(&m_coll).Shift(Amount);
	end_function
}

STDMETHODIMP CVolatilityStructures::Stick(void)
{
	begin_function
	CComObjectCollectionFunctions<IVolatilityStructure>(&m_coll).Stick();
	end_function
}