//	handle_strikes.cpp : Implementation of CStrikes
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_strikes.h"

STDMETHODIMP CStrikes::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IStrikes };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CStrikes::get_Value(VARIANT *pVal)
{
	CParameterMap						pmStartDate;
	
	begin_function
	if (!m_h->m_hStrikes.size()){
		// no strikes defined
		return ::VariantClear(pVal);
	}
	// Strike Type (we take the strike type from the first element in m_h->m_hStrikes)
	unmap_parameter(m_h->m_hStrikes[0]->m_strikeType, pmStrikeType);
	
	// StartDate (not always used)
	pmStartDate.SetToObjectPointer(m_spStartDate);

	// Strikes is always the zeroth parameter
	std::vector<double> vStrikes(m_h->m_hStrikes.size());
	for (long n = 0; n < m_h->m_hStrikes.size(); n++){
		if (m_h->m_hStrikes[n]->m_strikeType != m_h->m_hStrikes[0]->m_strikeType) throw "Inconsistent strike types encountered";
		vStrikes[n] = m_h->m_hStrikes[n]->m_strike;
	}
	unmap_parameter(vStrikes, pmStrikes);	

	// The remaining parameters depend on the strikes type
	MlEqStrike* pStrikeRef = &*(m_h->m_hStrikes[0]);		
	switch (m_h->m_hStrikes[0]->m_strikeType){
	case Fixed:
		return CParameterMap::VariableArgumentListToArray(pVal, 2, pmStrikeType, pmStrikes);
	case Normalised:
		{
			MlEqNormalizedStrike* p = dynamic_cast<MlEqNormalizedStrike*>(pStrikeRef);
			unmap_parameter(p->m_fwdparameter, pmForward);
			unmap_parameter(p->m_vol, pmVol);
			unmap_parameter(p->GetEndDate(), pmEndDate);
			unmap_parameter(p->GetIsFloating(), pmIsFloating);
			return CParameterMap::VariableArgumentListToArray(pVal, 7, pmStrikeType, pmStrikes, pmForward, pmVol, pmEndDate, pmStartDate, pmIsFloating);
		}
	case SpotBased:
		{
			MlEqSpotBasedStrike* p = dynamic_cast<MlEqSpotBasedStrike*>(pStrikeRef);
			unmap_parameter(p->m_spotparameter, pmSpot);
			unmap_parameter(p->GetIsFloating(), pmIsFloating);
			return CParameterMap::VariableArgumentListToArray(pVal, 5, pmStrikeType, pmStrikes, pmSpot, pmStartDate, pmIsFloating);
		}		
	case ForwardBased:
		{
			MlEqForwardBasedStrike* p = dynamic_cast<MlEqForwardBasedStrike*>(pStrikeRef);
			unmap_parameter(p->m_fwdparameter, pmForward);
			unmap_parameter(p->GetEndDate(), pmMaturity);
			unmap_parameter(p->GetIsFloating(), pmIsFloating);
			return CParameterMap::VariableArgumentListToArray(pVal, 6, pmStrikeType, pmStrikes, pmForward, pmMaturity, pmStartDate, pmIsFloating);
		}
	case LogBased:
		{
			MlEqLogStrike* p = dynamic_cast<MlEqLogStrike*>(pStrikeRef);
			unmap_parameter(p->m_fwdparameter, pmLog);
			unmap_parameter(p->GetEndDate(), pmMaturity);
			unmap_parameter(p->GetIsFloating(), pmIsFloating);
			return CParameterMap::VariableArgumentListToArray(pVal, 6, pmStrikeType, pmStrikes, pmLog, pmMaturity, pmStartDate, pmIsFloating);
		}
	}
	throw "Unhandled exception in CStrikes::get_Value";
	end_function
}

STDMETHODIMP CStrikes::put_Value(VARIANT newVal)
{
	HRESULT								hr;	
	std::vector<CParameterMap>			vpm;
	CComPtr<IDate>						spStartDate;
		
	begin_function		
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	if (vpm.size() == 1){
		// in this case we create a fixed strike handle
		CParameterMap pm;
		if (pm.SetValue(Fixed)) ATLASSERT(false);
		vpm.insert(vpm.begin(), pm);
	}
		
	map_parameter(vpm[0], long, StrikeType);		
	switch (StrikeType){
	case Fixed:
		if (vpm.size() != 2){
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			map_parameter(vpm[1], CVector, Strikes);
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqStrike(Strikes[n]);
			}
		}
		break;
	case Normalised:
		if (vpm.size() < 4 || vpm.size() > 7){
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			vpm.resize(7);
			map_parameter(vpm[1], CVector, Strikes);
			map_parameter(vpm[2], double, ReferenceForward);
			map_parameter(vpm[3], double, ReferenceVolatility);
			map_parameter(vpm[4], long, nMaturity);
			map_optional_object_parameter(vpm[5], Date, StartDateHandleOpt);
			map_optional_parameter(vpm[6], bool, IsFloatingOpt, true);
			if (!StartDateHandleOpt) StartDateHandleOpt = new MlEqDate(MlEqDate::GetCurrentDate());
			spStartDate = _sStartDateHandleOpt;
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqNormalizedStrike(Strikes[n], nMaturity, ReferenceVolatility, ReferenceForward, IsFloatingOpt, StartDateHandleOpt);
			}
		}
		break;
	case StrikesTypeFromAsset::NormalisedFromAsset:
		// create normalised strike from asset
		if (vpm.size() < 4 || vpm.size() > 6){
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			vpm.resize(6);
			map_parameter(vpm[1], CVector, Strikes);
			map_object_parameter(vpm[2], Asset, hAsset);
			map_parameter(vpm[3], long, nMaturity);
			map_optional_object_parameter(vpm[4], Date, StartDateHandleOpt);
			map_optional_parameter(vpm[5], bool, IsFloatingOpt, true);
			if (!StartDateHandleOpt) StartDateHandleOpt = new MlEqDate(*hAsset->GetDateHandle());
			spStartDate = _sStartDateHandleOpt;
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqNormalizedStrike(*hAsset, Strikes[n], nMaturity, IsFloatingOpt, StartDateHandleOpt);
			}
		}
		break;
	case SpotBased:
		if (vpm.size() < 3 || vpm.size() > 5){
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			vpm.resize(5);
			map_parameter(vpm[1], CVector, Strikes);
			map_parameter(vpm[2], double, ReferenceSpot);
			map_optional_object_parameter(vpm[3], Date, StartDateHandleOpt);
			map_optional_parameter(vpm[4], bool, IsFloatingOpt, true);
			if (!StartDateHandleOpt) StartDateHandleOpt = new MlEqDate(MlEqDate::GetCurrentDate());
			spStartDate = _sStartDateHandleOpt;
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqSpotBasedStrike(Strikes[n], ReferenceSpot, StartDateHandleOpt, IsFloatingOpt);
			}
		}
		break;
	case StrikesTypeFromAsset::SpotBasedFromAsset:
		if (vpm.size() < 3 || vpm.size() > 5){
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			vpm.resize(5);
			map_parameter(vpm[1], CVector, Strikes);
			map_object_parameter(vpm[2], Asset, hAsset);
			map_optional_object_parameter(vpm[3],Date, StartDateHandleOpt);			
			map_optional_parameter(vpm[4], bool, IsFloatingOpt, true);
			if (!StartDateHandleOpt) StartDateHandleOpt = new MlEqDate(*hAsset->GetDateHandle());
			spStartDate = _sStartDateHandleOpt;
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);			
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqSpotBasedStrike(*hAsset, Strikes[n], StartDateHandleOpt, IsFloatingOpt);
			}
		}
		break;
	case ForwardBased:		
		if (vpm.size() < 4 || vpm.size() > 6){
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			vpm.resize(6);
			map_parameter(vpm[1], CVector, Strikes);
			map_parameter(vpm[2], double, ReferenceForward);	
			map_parameter(vpm[3], long, nMaturity);
			map_optional_object_parameter(vpm[4], Date, StartDateHandleOpt);
			map_optional_parameter(vpm[5], bool, IsFloatingOpt, true);
			if (!StartDateHandleOpt) StartDateHandleOpt = new MlEqDate(MlEqDate::GetCurrentDate());
			spStartDate = _sStartDateHandleOpt;
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);			
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqForwardBasedStrike(Strikes[n], ReferenceForward, nMaturity, StartDateHandleOpt, IsFloatingOpt);
			}
		}
		break;
	case StrikesTypeFromAsset::ForwardBasedFromAsset:		
		if (vpm.size() < 4 || vpm.size() > 6){
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			vpm.resize(6);
			map_parameter(vpm[1], CVector, Strikes);
			map_object_parameter(vpm[2], Asset, hAsset);
			map_parameter(vpm[3], long,nMaturity);
			map_object_parameter(vpm[4],Date, StartDateHandleOpt);			
			map_parameter(vpm[5], bool, IsFloatingOpt);			
			if (!StartDateHandleOpt) StartDateHandleOpt = new MlEqDate(*hAsset->GetDateHandle());
			spStartDate = _sStartDateHandleOpt;
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);			
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqForwardBasedStrike(*hAsset, Strikes[n], nMaturity, StartDateHandleOpt, IsFloatingOpt);
			}
		}
		break;
	case LogBased:
		if (vpm.size() < 4 || vpm.size() > 6){
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			vpm.resize(6);
			map_parameter(vpm[1], CVector, Strikes);
			map_parameter(vpm[2], double, ForwardValue);			
			map_parameter(vpm[3], long, ReferenceMaturityDate);			
			map_optional_object_parameter(vpm[4], Date, StartDateHandleOpt);			
			map_optional_parameter(vpm[5], bool, IsFloatingOpt, true);
			if (!StartDateHandleOpt) StartDateHandleOpt = new MlEqDate(MlEqDate::GetCurrentDate());
			spStartDate = _sStartDateHandleOpt;
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);			
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqLogStrike(Strikes[n], ForwardValue, ReferenceMaturityDate, StartDateHandleOpt, IsFloatingOpt);
			}
		}	
		break;
	case StrikesTypeFromAsset::LogBasedFromAsset:		
		if (vpm.size() < 4 || vpm.size() > 6) {
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IStrikes);
		} else {
			map_parameter(vpm[1], CVector, Strikes);
			map_object_parameter(vpm[2], Asset, hAsset);
			map_parameter(vpm[3], long,nMaturity);
			map_object_parameter(vpm[4],Date, StartDateHandleOpt);
			map_parameter(vpm[5], bool, IsFloatingOpt);
			if (!StartDateHandleOpt) StartDateHandleOpt = new MlEqDate(*hAsset->GetDateHandle());
			spStartDate = _sStartDateHandleOpt;
			int nStrikes = Strikes.getsize();
			m_h->m_hStrikes.resize(nStrikes);			
			for (int n = 0; n < nStrikes; n++){
				m_h->m_hStrikes[n] = new MlEqLogStrike(*hAsset, Strikes[n], nMaturity, StartDateHandleOpt, IsFloatingOpt);
			}
		}
		break;
	default:
		return CParameterMap::ReturnErrorR(IDS_STRIKE_TYPE, IID_IStrikes);
	}
	
	// members for get_Value - ToDo - remove
	m_spStartDate = spStartDate;
	end_function
}

STDMETHODIMP CStrikes::get_StrikesType(StrikesTypeEnum* pVal)
{
	*pVal = (StrikesTypeEnum)m_h->m_hStrikes[0]->m_strikeType;
	return S_OK;
}