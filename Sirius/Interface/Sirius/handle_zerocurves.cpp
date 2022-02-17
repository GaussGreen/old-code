//	handle_zerocurves.cpp : Implementation of CZeroCurves
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "siriusapplication.h"
#include "handle_zerocurves.h"
#include "handle_zerocurve.h"
#include "MlEqDate.h"


/*static*/ HRESULT CZeroCurves::Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<IZeroCurves>& spZeroCurves)
{
	HRESULT								hr;	
	FileSystemEnum						fs;
	
	if (szDummy.size()) throw "The identifier value '" + szDummy + "' is invalid for the load zero curves request";
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();	
	if (hr = spZeroCurves.CoCreateInstance(CLSID_ZeroCurves)) return hr;
	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		std::stringstream				ssQuery;
		std::string						szError = "No zero curves found on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		ssQuery << "sp_user_get_zerocurves " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
		hr = CComObjectCollectionFunctions<IZeroCurve>(&szError).ImplementLoadCollection_SQL(ssQuery, ds, spZeroCurves);
		return hr;
	} else if (fs == fsNTFS){										
		return CComObjectCollectionFunctions<IZeroCurve>().ImplementLoadCollection_NTFS(spZeroCurves, 0, "", _Module.GetLocation(), "", ds, nDate, CZeroCurve::Load);
	}
	return E_FAIL;
}

HRESULT CZeroCurves::LoadGDAObject(const std::string& szName, long nDate, CComPtr<IZeroCurve>& spObject)
{
	// We don't use the check_load_enabled enabled here since any GDA objects are effectively loaded.
	try {
		estring szIndex;
		szIndex = ("YC." + szName);					
		GDA::Context_ref					ctx = GDA::GetLibraryInstance()->getDefaultNamespace()->getContext();
		bool								bHasCurve = ctx->contains(szIndex.c_str());		
		GDA::HDElement						hdeYieldCurve = ctx->resolve(szIndex.c_str());

		if (hdeYieldCurve.type() == GDA::HDElement::ConstObject){		
			GDA::Functor_const_ref			hYieldCurve(hdeYieldCurve.asConstObject());
			MlEqZeroCurveHandle				hZeroCurve = new MlEqZeroCurve;
			hZeroCurve->PutYieldCurve(hYieldCurve);

			// Only add this curve if the date matches.
			if (hZeroCurve->GetReferenceDate() == nDate){				
				map_analytic_to_com(hZeroCurve, ZeroCurve, spZeroCurve);
				spObject = spZeroCurve;
				spZeroCurve->put_DataSource(Last);			
				return S_OK;
			}
		}
		return E_FAIL;
	} catch (...){
		// All errors originate from GDA. We don't care about them as it corresponds to an object
		// not being loaded.
		return E_FAIL;
	}		
}

STDMETHODIMP CZeroCurves::Reset()
{
	begin_function
	CComObjectCollectionFunctions<IZeroCurve>(&m_coll).Reset();
	end_function
}

STDMETHODIMP CZeroCurves::Shift(double Amount)
{
	begin_function
	CComObjectCollectionFunctions<IZeroCurve>(&m_coll).Shift(Amount);
	end_function
}

STDMETHODIMP CZeroCurves::Stick(void)
{
	begin_function
	CComObjectCollectionFunctions<IZeroCurve>(&m_coll).Stick();
	end_function
}