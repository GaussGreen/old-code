//	handle_svlvoldata.cpp : Implementation of CSVLVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_svlvoldata.h"
#include "handle_strikes.h"
#include "siriusapplication.h"

STDMETHODIMP CSVLVolData::get_BasicSvlParameters(VARIANT* pVal)
{
	CParameterMap						pm;
	
	begin_function
	GetFromOld(m_h->getSvlParameters(), &pm);
	return pm.GetValue(pVal);
	end_function
}

STDMETHODIMP CSVLVolData::get_BasicValue(VARIANT *pVal)
{
	begin_function
	return GetValue(pVal, &GetFromOld);
	end_function
}

/*static*/ bool CSVLVolData::GetFromNew(const std::vector<SvlCoeffs>& svl, CParameterMap* ppm)
{
	ppm->SetSize(svl.size(), 6);
	for (long n = 0 ; n < svl.size(); n++){
		// Get the 'new' members of the SvlCoeffs structure		
		if (!svl[n].bNewValid) throw "The svi parameterisation values are not valid";
		ppm->SetValue(n, 0, svl[n].ATMVolF);	 
		ppm->SetValue(n, 1, svl[n].Skew);			
		ppm->SetValue(n, 2, svl[n].CallWing);
		ppm->SetValue(n, 3, svl[n].PutWing);
		ppm->SetValue(n, 4, svl[n].VolMin);
		ppm->SetValue(n, 5, svl[n].t);
	}
	return true;
}

/*static*/ bool CSVLVolData::GetFromOld(const std::vector<SvlCoeffs>& svl, CParameterMap* ppm)
{
	ppm->SetSize(svl.size(), 6);
	for (int n = 0; n < svl.size(); n++){
		// Get the 'old' members of the SvlCoeffs structure
		if (!svl[n].bOldValid) throw "The svl parameterisation values are not valid";
		ppm->SetValue(n, 0, svl[n].a);	 
		ppm->SetValue(n, 1, svl[n].b);			
		ppm->SetValue(n, 2, svl[n].rho);
		ppm->SetValue(n, 3, svl[n].m);
		ppm->SetValue(n, 4, svl[n].sig);
		ppm->SetValue(n, 5, svl[n].t);
	}
	return false;
}

STDMETHODIMP CSVLVolData::get_SvlParameters(VARIANT* pVal)
{
	CParameterMap						pm;
	
	begin_function
	GetFromNew(m_h->getSvlParameters(), &pm);
	return pm.GetValue(pVal);
	end_function
}

HRESULT CSVLVolData::GetValue(VARIANT* pVal, bool (*Get)(const std::vector<SvlCoeffs>& , CParameterMap*))
{
	HRESULT								hr;	
	const vector<SvlCoeffs>				svl = m_h->getSvlParameters();	
	CParameterMap						pmStrikes, pmCoeffs, pmSviParametersGivenOpt;
	
	if (hr = pmStrikes.SetValue(m_aspStrike)) return hr;
	bool b = Get(svl, &pmCoeffs);	// Will return true if this function sets from Svl and false if set from Svi.
	pmSviParametersGivenOpt.SetValue(!b);
	return CParameterMap::VariableArgumentListToArray(pVal, 3, pmStrikes, pmCoeffs, pmSviParametersGivenOpt);
}

STDMETHODIMP CSVLVolData::get_Value(VARIANT *pVal)
{
	begin_function	
	// Attempt to use the svl parameters. However, if these are invalid then use svi.
	bool bNewValid = true;
	for (long n = 0; n < m_h->getSvlParameters().size(); ++n){
		if (!m_h->getSvlParameters()[n].bNewValid){
			bNewValid = false;
			break;
		}
	}
	return GetValue(pVal, bNewValid ? &GetFromNew : &GetFromOld);
	end_function
}

STDMETHODIMP CSVLVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_ISVLVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CSVLVolData::put_BasicSvlParameters(VARIANT newVal)
{
	std::vector<SvlCoeffs>				svl;
	
	begin_function
	map_parameter(newVal, CParameterMap, pm);
	if (pm.GetCols() != 6) throw CStringResource(IDS_COLUMNS_INVALID);
	PutFromOld(pm, &svl);
	m_h->initialize(m_h->getStrikes(), svl);
	return S_OK;
	end_function
}

STDMETHODIMP CSVLVolData::put_BasicValue(VARIANT newVal)
{
	begin_function
	HRESULT											hr;
	std::vector<CComVariant>						vv;	

	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;		
	PutValue(vv, &PutFromOld);
	end_function
}

/*static*/ bool CSVLVolData::PutFromNew(const CParameterMap& pm, std::vector<SvlCoeffs>* psvl)
{
	psvl->resize(pm.GetRows());
	for (long n = 0; n < pm.GetRows(); n++){
		pm.GetValue(n, 0, &(*psvl)[n].ATMVolF);
		pm.GetValue(n, 1, &(*psvl)[n].Skew);			
		pm.GetValue(n, 2, &(*psvl)[n].CallWing);
		pm.GetValue(n, 3, &(*psvl)[n].PutWing);
		pm.GetValue(n, 4, &(*psvl)[n].VolMin);
		pm.GetValue(n, 5, &(*psvl)[n].t);
		(*psvl)[n].bNewValid = true;
		MlSvlVolData::initOldSvl((*psvl)[n]);
	}
	return true;
}

/*static*/ bool CSVLVolData::PutFromOld(const CParameterMap& pm, std::vector<SvlCoeffs>* psvl)
{
	psvl->resize(pm.GetRows());
	for (long n = 0; n < pm.GetRows(); n++){
		pm.GetValue(n, 0, &(*psvl)[n].a);
		pm.GetValue(n, 1, &(*psvl)[n].b);					
		pm.GetValue(n, 2, &(*psvl)[n].rho);
		pm.GetValue(n, 3, &(*psvl)[n].m);
		pm.GetValue(n, 4, &(*psvl)[n].sig);
		pm.GetValue(n, 5, &(*psvl)[n].t);
		(*psvl)[n].bOldValid = true;
		MlSvlVolData::initNewSvl((*psvl)[n], false);		
	}
	return false;
}

STDMETHODIMP CSVLVolData::put_SvlParameters(VARIANT newVal)
{
	std::vector<SvlCoeffs>				svl;
	
	begin_function
	map_parameter(newVal, CParameterMap, pm);
	if (pm.GetCols() != 6) throw CStringResource(IDS_COLUMNS_INVALID);
	PutFromNew(pm, &svl);
	m_h->initialize(m_h->getStrikes(), svl);
	return S_OK;
	end_function
}

void CSVLVolData::PutValue(const std::vector<CComVariant>& vv, bool (*Put)(const CParameterMap&, std::vector<SvlCoeffs>*))
{		
	std::vector<std::vector<MlEqStrikeHandle> >				pStrikes;			//[date][strike]
	vector<SvlCoeffs>										svl;
	
	if (vv.size() != 2) throw CStringResource(IDS_NUMBER_PARAMETERS);
	map_strikes_vector_parameter(vv[0], hStrikes);
	map_parameter(vv[1], CParameterMap, Parameters);
	if (Parameters.GetCols() != 6) throw CStringResource(IDS_COLUMNS_INVALID);
	Put(Parameters, &svl);
	m_h->initialize(hStrikes, svl);
	m_aspStrike = _shStrikes;	
}

STDMETHODIMP CSVLVolData::put_Value(VARIANT newVal)
{
	/*parameter list is 0 - Strikes
						1 - Parameters
						2 - SviParametersGivenOpt*/		

	begin_function
	HRESULT											hr;
	std::vector<CComVariant>						vv;
	bool											bSviParametersGivenOpt = false;

	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() == 3){
		map_optional_parameter(vv[2], bool, SviParametersGivenOpt, false);
		bSviParametersGivenOpt = SviParametersGivenOpt;
		vv.resize(2);
	}

	PutValue(vv, bSviParametersGivenOpt ? &PutFromOld : &PutFromNew);
	end_function
}