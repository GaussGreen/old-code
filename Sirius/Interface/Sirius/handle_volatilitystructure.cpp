//	handle_volatilitystructure.cpp : Implementation of CVolatilityStructure
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "comobjectcollectionfunctions.h"
#include "handle_volatilitystructure.h"
#include "handle_asset.h"
#include "handle_interpolator.h"
#include "handle_strikes.h"
#include "handle_date.h"
#include "handle_voldata.h"
#include "handle_jumpwingvoldata.h"
#include "handle_svlvoldata.h"
#include "handle_svlfitvoldata.h"
#include "handle_hullvoldata.h"
#include "siriusapplication.h"
#include "handle_dateschedule.h"
#include "xmlstreamer.h"
#include "ComObjectCollectionSerialisableKey.h"

implement_member_variable_rekey(VolatilityStructure, VolatilityStructures, DATE, Date, long, DateNumber);
implement_member_variable_rekey(VolatilityStructure, VolatilityStructures, DataSourceEnum, DataSource, DataSourceEnum, DataSource);
implement_serialisable(VolatilityStructure);

STDMETHODIMP CVolatilityStructure::get_Forwards(IArray** pVal)
{		
	begin_function
	MlEqArrayHandle h = new MlEqArray(m_h->getForwards());
	map_analytic_to_com(h, Array, spArray);
	return spArray.CopyTo(pVal);
	end_function
}

STDMETHODIMP CVolatilityStructure::get_DiscountFactors(IArray** pVal)
{
	begin_function
	MlEqArrayHandle h = new MlEqArray(m_h->getDiscounts());
	map_analytic_to_com(h, Array, spArray);
	return spArray.CopyTo(pVal);
	end_function
}

//	Returns the index of m_aspSchedule which corresponds to an input interval.
//  (We use the first two dates in the schedule to deciper this - OK due
//  to input validation.)
long CVolatilityStructure::FindSchedule(const estring& szName) const
{
	std::string							szInterval;
	if (!szName.CompareNoCaseAndSpace("Annual")){
		szInterval = "1y";
	} else if (!szName.CompareNoCaseAndSpace("SemiAnnual")){
		szInterval = "6m";
	} else if (!szName.CompareNoCaseAndSpace("Quarterly")){
		szInterval = "3m";
	} else if (!szName.CompareNoCaseAndSpace("Monthly")){
		szInterval = "1m";
	}
	if (szInterval.size()){
		for (long n = 0; n < m_aspSchedules.size(); n++){
			map_com_to_analytic(m_aspSchedules[n].m_T, DateSchedule, hDateSchedule);
			if (hDateSchedule->size() < 2) continue;
			std::map<long, std::vector<double> >::const_iterator it = hDateSchedule->begin();
			long nFrom = it++->first;
			long nTo = it->first;
			if (MlEqDate(nFrom).GetTenor(nTo, NoBusinessDayConvention, "", false) == szInterval){
				return n;
			}
		}
	}
	throw "Schedule '" + szName + "' not found.";
	return -1;
}

STDMETHODIMP CVolatilityStructure::get_AsymptoticVolatilityStructure(IVolatilityStructure** pVal)
{
	if (m_spAsymptoticVolatilityStructure){		
		return m_spAsymptoticVolatilityStructure.CopyTo(pVal);
	} else {
		return E_POINTER;
	}	
}

STDMETHODIMP CVolatilityStructure::get_DateHandle(IDate** pVal)
{
	begin_function
	return m_spDateHandle.CopyTo(pVal);
	end_function	
}

STDMETHODIMP CVolatilityStructure::get_Dates(VARIANT* pVal)
{
	begin_function
	unmap_parameter(m_h->getDates(), pm);
	pm.GetColumn(0, VT_DATE, pVal);
	end_function
}

STDMETHODIMP CVolatilityStructure::get_InterpolateType(InterpolateTypeEnum* pVal)
{
	begin_function
	*pVal = (InterpolateTypeEnum)m_h->getInterpAcrossTimeFlag();
	end_function
}

STDMETHODIMP CVolatilityStructure::get_Name(BSTR* pVal)
{	
	begin_function	
	return estring(m_szName).GetBSTR(pVal);	
	end_function
}

estring CVolatilityStructure::GetLocation(void) const
{
	if (m_szLocation.size()){
		return m_szLocation;
	} else {
		return _Module.GetLocation();
	}
}
STDMETHODIMP CVolatilityStructure::get_Location(BSTR* pVal)
{
	return estring::GetBSTR(GetLocation(), pVal);
}

STDMETHODIMP CVolatilityStructure::get_Parameter(BSTR Name, VARIANT Index1Opt, VARIANT Index2Opt, VARIANT Index3Opt, VARIANT *pVal)
{	
	estring								szName(Name);
	MlEqStochBetaVolatilityStructure*	pvs;
	double								f;
	
	begin_function		
	if (!szName.CompareNoCaseAndSpace("Level")){
		map_parameter(Index1Opt, long, n1);
		map_parameter(Index2Opt, long, n2);
		if (!(pvs = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_h))) return CParameterMap::ReturnErrorRS(IDS_CANNOT_FIND_PARAMETER, Name, IID_IVolatilityStructure);
		f = pvs->getLevel(n1, n2);
		return CComVariant(f).Detach(pVal);	
	} else if (!szName.CompareNoCaseAndSpace("Beta")){
		map_parameter(Index1Opt, long, n1);
		map_parameter(Index2Opt, long, n2);
		if (!(pvs = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_h))) return CParameterMap::ReturnErrorRS(IDS_CANNOT_FIND_PARAMETER, Name, IID_IVolatilityStructure);
		f = pvs->getBeta(n1, n2);
		return CComVariant(f).Detach(pVal);
	} else if (!szName.CompareNoCaseAndSpace("VolVol")){
		map_parameter(Index1Opt, long, n1);
		map_parameter(Index2Opt, long, n2);
		if (!(pvs = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_h))) return CParameterMap::ReturnErrorRS(IDS_CANNOT_FIND_PARAMETER, Name, IID_IVolatilityStructure);
		f = pvs->getVolVol(n1, n2);
		return CComVariant(f).Detach(pVal);
	} else if (!szName.CompareNoCaseAndSpace("MeanReversion")){
		map_parameter(Index1Opt, long, n1);
		map_parameter(Index2Opt, long, n2);
		if (!(pvs = dynamic_cast<MlEqStochBetaVolatilityStructure*>(&*m_h))) return CParameterMap::ReturnErrorRS(IDS_CANNOT_FIND_PARAMETER, Name, IID_IVolatilityStructure);
		f = pvs->getMeanReversion(n1, n2);
		return CComVariant(f).Detach(pVal);
	}

	// Try extracting a beta schedule with the name Name last.
	long nSchedule;
	try {
		nSchedule = FindSchedule(szName);
	} catch (const std::string&){
		throw "'" + szName + "' is not a parameter associated with the volatility structure '" + (std::string)CComObjectCollectionSerialisableKey(CComPtr<IVolatilityStructure>(this)) + "'";
	}
	return CComVariant(m_aspSchedules[nSchedule].m_T).Detach(pVal);
	end_function
}

STDMETHODIMP CVolatilityStructure::get_Spot(double* pVal)
{
	begin_function
	*pVal = m_h->getSpot();
	end_function
}

STDMETHODIMP CVolatilityStructure::get_Value(VARIANT *pVal)
{			
	CParameterMap		pmDate;							
	
	pmDate.SetToObjectPointer(m_spDateHandle);
	if (m_vpm.size() > 6) m_vpm[6] = pmDate;
	
	if (m_vpm.size() > 2){
		std::string szDataSource = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource());
		m_vpm[1].SetValue(szDataSource);
	}

	// The asymptotic volatility structure needs to be fully identifiable
	if (m_vpm.size() > 13){		
		m_vpm[13].SetValue(CComObjectCollectionSerialisableKey(m_spAsymptoticVolatilityStructure).GetKeyAndObjectName(false, false));
	}
	
	
	return CParameterMap::VectorToArray(m_vpm, pVal);

	// Note that we don't return any shift data. This is by design.	
	/*CParameterMap		pmName;							pmName.SetValue(m_szName);
	std::string							szDataSource = CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource());
	unmap_parameter(szDataSource, pmDataSource);
	 CParameterMap		pmSpot;							pmSpot.SetValue(m_h->getSpot());
	CParameterMap		pmForwards;						pmForwards.SetValue(m_h->getForwards());
	CParameterMap		pmDiscountFactors;				pmDiscountFactors.SetValue(m_h->getDiscounts());	
	CParameterMap		pmTimes;						pmTimes.SetValue(m_h->getDates());
	CParameterMap		pmDate;							pmDate.SetToObjectPointer(m_spDate);
	CParameterMap		pmInterpolator;					pmInterpolator.SetToObjectPointer(m_spInterpolator);	
	CParameterMap		pmBid;							pmBid.SetValue(m_aspBid);
	CParameterMap		pmMid;							pmMid.SetValue(m_aspMid);
	CParameterMap		pmAsk;							pmAsk.SetValue(m_aspAsk);
	CParameterMap		pmInterpolateAcrossTimeFlag;	pmInterpolateAcrossTimeFlag.SetValue(m_h->getInterpAcrossTimeFlag());		
	unmap_parameter(m_h->m_volType, pmVolType);

						
	if (m_aspSchedules.size()){
		// Stochastic Beta volatility structure
		ATLASSERT(false);
		CParameterMap	pmStochasticBeta;				pmStochasticBeta.SetValue(m_aspSchedules);		
		CParameterMap		pmAsymptoticVolatilityStructure;		pmAsymptoticVolatilityStructure.SetToObjectPointer(m_spAsymptoticVolatilityStructure);
		return CParameterMap::VariableArgumentListToArray(pVal, 15, pmName, pmDataSource, pmSpot, pmForwards, pmDiscountFactors, pmTimes,	pmDate,	pmInterpolator, pmBid, pmMid, pmAsk, pmInterpolateAcrossTimeFlag, pmVolType, pmAsymptoticVolatilityStructure, pmStochasticBeta);
	} else if (m_spAsymptoticVolatilityStructure){
		ATLASSERT(false);
		CParameterMap		pmAsymptoticVolatilityStructure;		pmAsymptoticVolatilityStructure.SetToObjectPointer(m_spAsymptoticVolatilityStructure);
		return CParameterMap::VariableArgumentListToArray(pVal, 14, pmName, pmDataSource, pmSpot, pmForwards, pmDiscountFactors, pmTimes,	pmDate,	pmInterpolator, pmBid, pmMid, pmAsk, pmInterpolateAcrossTimeFlag, pmVolType, pmAsymptoticVolatilityStructure);
	} else {
		// normal volatility structure
		return CParameterMap::VariableArgumentListToArray(pVal, 13, pmName, pmDataSource, pmSpot, pmForwards, pmDiscountFactors, pmTimes,	pmDate,	pmInterpolator, pmBid, pmMid, pmAsk, pmInterpolateAcrossTimeFlag, pmVolType);
	}*/
}

STDMETHODIMP CVolatilityStructure::GetVolatility(VARIANT Strikes, DATE Maturity, BidAskMidEnum BidAskMid, VARIANT *pResult)
{
	begin_function
	long								nStrikes;		
		
	map_object_parameter(Strikes, Strikes, hStrike);
	map_parameter(Maturity, DATE, date);
				
	// get the volatility for each strike	
	nStrikes = hStrike->m_hStrikes.size();
	CVector vector(nStrikes);
	for (int i = 0; i < nStrikes; i++ ){		
		vector[i] = m_h->getVol(*hStrike->m_hStrikes[i], date, BidAskMid);
	}
	if (nStrikes == 1){
		return CComVariant(vector[0]).Detach(pResult);
	} else {
		CParameterMap pm;
		pm.SetValue(vector);
		return pm.GetValue(pResult);
	}
	end_function
}

STDMETHODIMP CVolatilityStructure::get_VolatilityDataType(VolatilityDataTypeEnum* pVal)
{	
	begin_function
	*pVal = m_h->m_volType;
	end_function
}

STDMETHODIMP CVolatilityStructure::get_VolatilityData(BidAskMidEnum BidAskMid, VARIANT* pVal)
{
	begin_function		
	CParameterMap pm(GetVolatilityData(BidAskMid));
	if (pm.IsBlank()){
		throw CStringResource(IDS_NO_DATA_FOUND);
	} else if (pm.IsScalar()){
		return pm.GetValue(0, 0, pVal);
	} else if (pm.GetCols() == 1){
		// This is the normal case
		pm.GetColumn(0, VT_DISPATCH, pVal);		
	} else {
		// General case (should never be used since m_aspBid etc. is always a vector)
		ATLASSERT(false);	// Why did we end up here?
		return pm.GetValue(pVal);
	}
	end_function
}

CParameterMap CVolatilityStructure::GetVolatilityData(BidAskMidEnum BidAskMid) const
{
	CParameterMap pm;
	switch (BidAskMid){
	case Bid:		
		pm.SetValue(m_aspBid);
		break;
	case Middle:
		pm.SetValue(m_aspMid);		
		break;
	case Ask:
		pm.SetValue(m_aspAsk);
		break;
	default:
		throw "Unknown or handled value of BidAskMid";
	}
	return pm;
}

STDMETHODIMP CVolatilityStructure::get_VolatilityDataName(BidAskMidEnum BidAskMid, VARIANT* pVal)
{
	begin_function
	CParameterMap pm(GetVolatilityData(BidAskMid)), pmOut;
	pmOut.SetSize(pm.GetRows(), pm.GetCols());
	
	for (long nRow = 0; nRow < pm.GetRows(); nRow++){
		for (long nCol = 0; nCol < pm.GetCols(); nCol++){
			CComPtr<IDispatch>	sp;
			std::string			sz;			
			pm.GetValue(nRow, nCol, sp);
			CParameterMap::GetObjectName(sp, &sz);
			pmOut.SetValue(nRow, nCol, sz);
		}
	}
	
	if (pmOut.IsBlank()){
		throw CStringResource(IDS_NO_DATA_FOUND);
	} else if (pmOut.IsScalar()){
		return pmOut.GetValue(0, 0, pVal);
	} else if (pmOut.GetCols() == 1){
		// This is the normal case		
		pmOut.GetColumn(0, VT_BSTR, pVal);
	} else {
		// General case (should never be used since m_aspBid etc. is always a vector)
		ATLASSERT(false);	// Why did we end up here?
		return pmOut.GetValue(pVal);
	}
	end_function
}

STDMETHODIMP CVolatilityStructure::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IVolatilityStructure };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

/*static*/ HRESULT CVolatilityStructure::Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<IVolatilityStructure>& spVolatilityStructure)
{
	FileSystemEnum						fs;	
			
	if (ds == NoDataSource) ds = CComObjectCollectionSerialisableDefaulter::GetDataSource();
	if (!nDate) nDate = CComObjectCollectionSerialisableDefaulter::GetDate();

	// load the XML for the volatility structure
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute
		_RecordsetPtr					prs;							// ADO recordset
		CComVariant						vField;							// volatility structure field value
		estring							szError;
		
		if (!_Module.GetLocationFirstTryEnabled() || !_Module.GetLocationFirstTry().size() || !estring::CompareNoCase(_Module.GetLocationFirstTry(), _Module.GetLocation())){
			// not using trial location			
			ssQuery << "sp_user_get_volatilitystructure '" << szIdentifier << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocation() << "'";
			szError = "Cannot find volatility structure '" + szIdentifier + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		} else {
			// are using trial location
			ssQuery << "sp_user_get_volatilitystructure '" << szIdentifier << "', " << nDate << ", '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) << "', '" << _Module.GetLocationFirstTry() << "', '" << _Module.GetLocation() << "'";
			szError = "Cannot find volatility structure '" + szIdentifier + "' on date '" + MlEqDate(nDate).GetString() + "' in location '" + _Module.GetLocationFirstTry() + "' or '" + _Module.GetLocation() + "' with data source '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, ds) + "'";
		}
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
			if (prs->adoEOF) throw szError;
			vField = prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue();
		} catch (_com_error& e){
			throw estring(e);
		}
		if (vField.vt != VT_BSTR) return E_FAIL;
		
		// parse the XML into a variant and set this volatility structure object to that variant
		return CXmlStreamer::GetObject((char*)_bstr_t(vField.bstrVal), ds, spVolatilityStructure);
	} else if (fs == fsNTFS){
		if (!_Module.GetLocationFirstTryEnabled() || !_Module.GetLocationFirstTry().size() || !estring::CompareNoCase(_Module.GetLocationFirstTry(), _Module.GetLocation())){
			// not using trial location			
			CComObjectCollectionFunctions<IVolatilityStructure>::ImplementLoad_NTFS(IID_IVolatilityStructure, _Module.GetLocation(), "", szIdentifier, ds, nDate, spVolatilityStructure);
		} else {
			// are using trial location
			CComObjectCollectionFunctions<IVolatilityStructure>::ImplementLoad_NTFS(IID_IVolatilityStructure, _Module.GetLocationFirstTry(), _Module.GetLocation(), szIdentifier, ds, nDate, spVolatilityStructure);
		}
		return S_OK;		
	} else {
		return E_FAIL;
	}		
}

void CVolatilityStructure::Make2DInterpolator(const std::map<long, std::map<long, double> > map, MlEq2DInterpolatorHandle& hOut) const
{
	MlEqInterpolatorHandle											hiy;
	std::vector<MlEqInterpolatorHandle>								ahix;
	std::map<long, double>::const_iterator							itx;
	std::map<long, std::map<long, double> >::const_iterator			ity;
	CVector															ax, ay;
	int																n;
	
	// get the 'y' interpolator
	ay.resize(map.size());
	n = 0;	
	for (ity = map.begin(); ity != map.end(); ity++){
		ay[n] = (double)ity->first;
		n++;
	}	
	hiy = new MlEqInterpolator(ay, ay);

	// get the 'x' interpolators
	for (ity = map.begin(); ity != map.end(); ity++){				
		ax.resize(ity->second.size());
		ay.resize(ity->second.size());
		n = 0;						
		for (itx = ity->second.begin(); itx != ity->second.end(); itx++){
			ax[n] = (double)itx->first;
			ay[n] = itx->second;
			n++;
		}
		MlEqInterpolatorHandle h = new MlEqInterpolator(ax, ay);
		ahix.push_back(h);
	}

	// create the 2D interpolator
	hOut = new MlEq2DInterpolator(hiy, ahix);
}

STDMETHODIMP CVolatilityStructure::put_Location(BSTR newVal)
{
	m_szLocation.Set(newVal);
	m_szLocation.trim();
	return S_OK;
}

STDMETHODIMP CVolatilityStructure::put_Value(VARIANT newVal)
{	
	/*parameter list is 0 - Name
						1 - DataSourceAndOrLocation	(can be of the form DataSource (Location), Location (DataSource), DataSource, Location)
						2 - Spot
						3 - Forwards
						4 - Discounts
						5 - Times
						6 - DateHandle
						7 - Interpolator
						8 - BidVolData
						9 - MidVolData
						10 - AskVolData
						11 - InterpolateType
						12 - VolatilityDataTypeEnum
						13 - AsymptoticVolatilityStructure (optional only if the next two parameters are blank)
						14 - DecayFactor (optional only if the previous parameter is blank)
	  for a stochastic beta volatility structure, the list continues:
						15 - StochasticBetaSchedules
						16 - ExcessVolVol2DInterpolatorOpt 	
						17 - HermiteElasticityOpt 
	*/			
	HRESULT											hr;
	std::vector<CParameterMap>						vpm;	
	
	// temporary assignments	
	std::vector<CAdapt<CComQIPtr<IDispatch> > >		aspBid;
	std::vector<CAdapt<CComQIPtr<IDispatch> > >		aspAsk;
	std::vector<CAdapt<CComQIPtr<IDispatch> > >		aspMid;
	std::vector<MlEqVolDataHandle>					volBidHandles;	
	std::vector<MlEqVolDataHandle>					volMidHandles;	
	std::vector<MlEqVolDataHandle>					volAskHandles;
	CComPtr<IVolatilityStructure>					spAsymptotic;
	CComPtr<IInterpolator>							spExcessVolVolInterpolator;
	CComPtr<IInterpolator>							spElasticityInterpolator;
	std::vector<CAdapt<CComQIPtr<IDateSchedule> > >	aspSchedules;			
	MlEqVolatilityStructureHandle					hAsymptoticVolatilityStructure = NULL;
	double											fDecayFactor = 0.0;	
	DataSourceEnum									ds = _Module.GetDefaultDataSource();	// set to default initially
	estring											szLocation = _Module.GetLocation();		// set to defaults
	
	// map the inputs			
	begin_function	
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	if (vpm.size() != 13 && vpm.size() != 15 && vpm.size() != 16 && vpm.size() != 17 && vpm.size() != 18) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IVolatilityStructure);
	map_parameter(vpm[0], estring, Name);
	
	// map vpm[1] to a data source DataSource and an estring called Location
	{
		estring	szA, szB;
		map_optional_parameter(vpm[1], estring, DataSource, "");
		estring::SplitBracket(DataSource, &szA, &szB);
		if (!szA.size()) szA.assign(szB);		
		if (szA.size()){
			// overrides have been specified					
			if (!CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szA, &ds)){
				if (szB.size()) szLocation.assign(szB);
			} else if (szB.size()){
				if (!CEnumMap::GetEnum("DataSourceEnum", LIBID_Sirius, szB, &ds)){
					szLocation.assign(szA);
				} else {
					// invalid if neither non-blank string is a data source
					throw "Invalid data source '" + DataSource + "'";
				}
			}
		}
	}
	
	map_parameter(vpm[2], double, Spot);
	map_parameter(vpm[3], CVector, Forwards);
	map_parameter(vpm[4], CVector, DiscountFactors);
	map_parameter(vpm[5], vector<DATE>, Times);
	map_object_parameter(vpm[6], Date, DateHandle);
	map_object_parameter(vpm[7], Interpolator, Interpolator);
	CObjectManager::Cast(vpm[8], aspBid, volBidHandles);
	CObjectManager::Cast(vpm[9], aspMid, volMidHandles);
	CObjectManager::Cast(vpm[10], aspAsk, volAskHandles);			
	map_enum_parameter(vpm[11], InterpolateTypeEnum, InterpAcrossTimeFlag);
	map_enum_parameter(vpm[12], VolatilityDataTypeEnum, voltype);
	if (vpm.size() > 13){
		// An asymptotic volatility structure has been passed in
		map_object_parameter(vpm[13], VolatilityStructure, AsymptoticVol);
		if (spAsymptotic = _sAsymptoticVol) hAsymptoticVolatilityStructure = AsymptoticVol;
		map_parameter(vpm[14], double, DecayFactor);
		fDecayFactor = DecayFactor;
	}

	// set the MlEqVolatilityStructure members		
	if (vpm.size() == 13){
		// Plain volatility structure case. No asymptotic volatility structure has been given.
		ATLASSERT(!spAsymptotic);
		MlEqVolatilityStructureHandle				h(new MlEqVolatilityStructure);
		h->initialize(Spot, DiscountFactors, Forwards, Times, DateHandle, volBidHandles, volMidHandles, volAskHandles, Interpolator, InterpAcrossTimeFlag, voltype);
		m_h = h;
	} else if (vpm.size() == 15){
		// The user is setting a non-stochastic beta volatiltiy structure with the final parameter being a handle
		// to an asymptotic volatility structure.
		ATLASSERT(spAsymptotic);
		MlEqVolatilityStructureHandle				h(new MlEqVolatilityStructure);
		h->initialize(Spot, DiscountFactors, Forwards, Times, DateHandle, volBidHandles, volMidHandles, volAskHandles, Interpolator, InterpAcrossTimeFlag, voltype, hAsymptoticVolatilityStructure, fDecayFactor);
		m_h = h;
	} else if (vpm.size() == 16 || vpm.size() == 17 || vpm.size() == 18){
		// The user is creating a Stochastic Beta volatility structure. The antepenultimate parameter is a handle
		// to an asymptotic volatility structure. The penultimate final parameter is a set of				
		// stochastic beta handles.
		ATLASSERT(spAsymptotic);		
		MlEqStochBetaVolatilityStructureHandle		h(new MlEqStochBetaVolatilityStructure);
		
		h->initialize(Spot, DiscountFactors, Forwards, Times, DateHandle, volBidHandles, volMidHandles, volAskHandles, Interpolator, InterpAcrossTimeFlag, voltype, hAsymptoticVolatilityStructure, fDecayFactor);
		m_h = &*h;
		
		// Set the stochastic beta schedule parts
		if (vpm[15].GetRows() == 1) vpm[15].Transpose();
		vpm[15].RemoveBlankRows();
		
		
		map_com_object_vector_parameter(vpm[15].GetValue(), DateSchedule, aspDateSchedules);

		aspSchedules = aspDateSchedules;
		
		std::map<long, std::map<long, double> >		mapLevel;						// outer map is by time interval (days); inner map is by date
		std::map<long, std::map<long, double> >		mapBeta;
		std::map<long, std::map<long, double> >		mapVolVol;
		std::map<long, std::map<long, double> >		mapMeanReversion;					
		
			// convert the schedule data into the above maps
		for (std::vector<CAdapt<CComQIPtr<IDateSchedule> > >::const_iterator it = aspSchedules.begin(); it != aspSchedules.end(); it++){
			CComVariant				vSchedule;
			CParameterMap			pmSchedule;
			long					nInterval;
			std::string				szInterval;
			
			it->m_T->get_Value(&vSchedule);
			pmSchedule.SetValue(vSchedule);								
			
			// get the average interval (in whole days)
			if (pmSchedule.GetRows() < 2){
				throw CStringResource(IDS_ROWS_INVALID);
			} else {								
		/*		long nFirst, nLast;
				if (pmSchedule.GetValue(0, 0, &nFirst)) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
				if (pmSchedule.GetValue(pmSchedule.GetRows() - 1, 0, &nLast)) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
				nInterval = (nLast - nFirst) / (pmSchedule.GetRows() - 1);
		*/
				// this allows only monthly, quarterly, semi or annually frequencies
				long nFirst, nSecond;
				if (pmSchedule.GetValue(0, 0, &nFirst)) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
				if (pmSchedule.GetValue(1, 0, &nSecond)) throw CStringResource(IDS_UNHANDLED_EXCEPTION);
				std::string szLocalInterval = MlEqDate(nFirst).GetTenor(nSecond, NoBusinessDayConvention, "", false);

				if	   ( szLocalInterval == "1m" ) nInterval = 1;
				else if( szLocalInterval == "3m" ) nInterval = 3;
				else if( szLocalInterval == "6m" ) nInterval = 6;
				else if( szLocalInterval == "1y" ) nInterval = 12;
				else 
					throw "This frequancy is not supported in the release";

			}
			
			// ToDo - check all the elements in the interval are constant w.r.t. MLEDateInv output
			// It doesn't matter if we have already processed a schedule with this interval - any common values will be overwritten.
							
			for (long nRow = 0; nRow < pmSchedule.GetRows() - 1; nRow++){
				long				nFrom, nTo;
				double				fValue;
				std::string			szLocalInterval;
				
				pmSchedule.GetValue(nRow, 0, &nFrom);
				pmSchedule.GetValue(nRow + 1, 0, &nTo);
				if (!nRow){
					szInterval = MlEqDate(nFrom).GetTenor(nTo, NoBusinessDayConvention, "", false);
				} else if (szInterval != (szLocalInterval = MlEqDate(nFrom).GetTenor(nTo, NoBusinessDayConvention, "", false))){
					throw "Non-constant interval(s) were encountered in the stochastic beta schedules";
				}
				
				if (!pmSchedule.GetValue(nRow + 1, 1, &fValue)){
					mapLevel[nInterval][nFrom] = fValue;
				}
				if (!pmSchedule.GetValue(nRow + 1, 2, &fValue)){
					mapBeta[nInterval][nFrom] = fValue;
				}
				if (!pmSchedule.GetValue(nRow + 1, 3, &fValue)){
					mapVolVol[nInterval][nFrom] = fValue;
				}
				if (!pmSchedule.GetValue(nRow + 1, 4, &fValue)){
					mapMeanReversion[nInterval][nFrom] = fValue;
				}
			}
		}
		
		Make2DInterpolator(mapLevel, h->m_hbiLevel);
		Make2DInterpolator(mapBeta, h->m_hbiBeta);
		Make2DInterpolator(mapVolVol, h->m_hbiVolVol);
		Make2DInterpolator(mapMeanReversion, h->m_hbiMeanReversion);
	
		// set the ExcessVolVolInterpolator (if supplied)
		if (vpm.size() >= 17){
			map_object_parameter(vpm[16], Interpolator, ExcessVolVol2DInterpolatorOpt);
			spExcessVolVolInterpolator = _sExcessVolVol2DInterpolatorOpt;
			MlEq2DInterpolatorHandle h2D;			
			if (!(h2D = dynamic_cast<MlEq2DInterpolator*>(&*ExcessVolVol2DInterpolatorOpt))){
				throw "You must pass in a 2 dimensional interpolator for parameter 'ExcessVolVol2DInterpolatorOpt'";
			}
			h->m_excessVolOfVol = h2D;
		}

		if (vpm.size() == 18){
			map_object_parameter(vpm[17], Interpolator, HermiteElasticityOpt);
			spElasticityInterpolator = _sHermiteElasticityOpt;
			h->m_hermiteElasticity = HermiteElasticityOpt;
		}
	} else {
		ATLASSERT(false);
	}
		
	// if all is well assign the member variables
	m_h->PutDataSource(ds);
	m_szName = Name;
	m_szLocation = szLocation;
	m_spDateHandle = _sDateHandle;
	m_spInterpolator = _sInterpolator;
	m_aspBid = aspBid;
	m_aspMid = aspMid;
	m_aspAsk = aspAsk;		
	m_aspSchedules = aspSchedules;
	m_spAsymptoticVolatilityStructure = spAsymptotic;
	m_spExcessVolVolInterpolator = spExcessVolVolInterpolator;
	m_spElasticityInterpolator = spElasticityInterpolator;
	m_vpm = vpm;
	end_function
}

STDMETHODIMP CVolatilityStructure::ReinitialiseStrike(IAsset* Asset)
{
	begin_function
	map_bare_com_to_analytic(Asset, Asset, hAsset);		
	m_h->reinitializeStrikes(*hAsset);
	end_function
}

STDMETHODIMP CVolatilityStructure::Save(BSTR *pVal)
//	pVal - returned, nullable
{
	HRESULT								hr;	
	xmlstreamer							ssXML;							// XML representation of vData			
	FileSystemEnum						fs;
	
	begin_function
	check_publishing_enabled
	if (!m_h->GetDateNumber()) throw "You cannot save the volatility structure '" + m_szName + "' since it does not have a date.";

	hr = CXmlStreamer::GetXML(CComPtr<IVolatilityStructure>(this), ssXML);	
	if (hr) return hr;

	// save it using the appropriate stored procedure	
	if ((fs = _Module.GetFileSystem()) == fsSQLServer){
		_ConnectionPtr					spConnection;
		std::stringstream				ssQuery;						// SQL query to execute	
		_RecordsetPtr					prs;							// ADO recordset
		
		if (m_h->GetDataSource() != Last) throw "Volatility structures with data sources other than 'Last' cannot be saved. The volatiltiy structure '" + m_szName + "' has the source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "'";
		ssQuery << "sp_user_add_volatilitystructure '" << m_szName << "', '" << m_h->getDateToDouble()->GetDate() << "', '" << CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) << "', '" << _Module.GetLocation() << "', '" << (char*)ssXML << "'";
		try {
			_Module.GetConnection(spConnection);
			prs = spConnection->Execute(_bstr_t(ssQuery.str().c_str()), NULL, -1);
		} catch (_com_error& e){
			throw estring(e);
		}
		// set the return type the returned value of the stored procedure
		CComVariant vValue(prs->GetFields()->GetItem(CComVariant(0, VT_I4))->GetValue());
		if (vValue.vt != VT_BSTR){
			ATLASSERT(false);
			return CParameterMap::ReturnErrorR(IDS_UNHANDLED_EXCEPTION);
		}		
		if (pVal) return CComBSTR(vValue.bstrVal).CopyTo(pVal);
		return S_OK;
	} else if (fs == fsNTFS){
		// we relax the condition that only 'Last' data sources can be saved
		if (m_h->GetDataSource() == NoDataSource) throw "The volatility structure '" + m_szName + "' has the data source set to '" + CEnumMap::GetString("DataSourceEnum", LIBID_Sirius, m_h->GetDataSource()) + "' and so cannot be saved";
		return CComObjectCollectionFunctions<IVolatilityStructure>::ImplementSave_NTFS(ssXML, IID_IVolatilityStructure, _Module.GetLocation(), m_szName, m_h->GetDataSource(), m_h->getDateToDouble()->GetDate(), pVal);
	}	
	end_function
}


/////////////////////////////////////////////////////////////////////////////
//	Scenario Functions
//

STDMETHODIMP CVolatilityStructure::FsaLiquidityShift(double AverageVolume, double CoreVolatility, double NetGamma, double PositionDirection, VARIANT_BOOL IsStock)
{
	begin_function
	m_h->FsaLiquidityShift(AverageVolume, CoreVolatility, NetGamma, PositionDirection, IsStock ? true : false);
	end_function
}

STDMETHODIMP CVolatilityStructure::Reset()
{
	begin_function
	m_h->Reset();
	end_function
}

STDMETHODIMP CVolatilityStructure::SetVol(double vol)
{	
	begin_function
	m_h->SetVol(vol);
	end_function
}

STDMETHODIMP CVolatilityStructure::Shift(double Amount, VARIANT_BOOL IsRelative)
{	
	begin_function
	m_h->parallelShiftVol(Amount, IsRelative ? true : false);
	end_function
}

STDMETHODIMP CVolatilityStructure::ShiftBeta(double Amount, BSTR Tenor)
{
	begin_function		
	m_h->setUpBetaShift(Amount, estring(Tenor));
	end_function
}

STDMETHODIMP CVolatilityStructure::ShiftForwardSkew(double shiftAmount, DATE WindowStartDate, DATE WindowEndDate, long PeriodStartDays, long PeriodEndDays)
{
	begin_function
	m_h->setUpFwdSkewShift(shiftAmount, WindowStartDate, WindowEndDate, PeriodStartDays, PeriodEndDays);
	end_function
}

STDMETHODIMP CVolatilityStructure::ShiftSkew(IStrikes* Strikes, IArray* Maturities, IMatrix* ShiftAcrossMaturity, IAsset* Asset)
{
	begin_function
	map_bare_com_to_analytic(Strikes, Strikes, hStrikes);		
	map_bare_com_to_analytic(Maturities, Array, hMaturities);
	map_bare_com_to_analytic(ShiftAcrossMaturity, Matrix, hShiftAcrossMaturity);
	map_bare_com_to_analytic(Asset, Asset, hAsset);
	m_h->shiftSkew(hStrikes->m_hStrikes, *hMaturities, *hShiftAcrossMaturity, hAsset);
	end_function
}

STDMETHODIMP CVolatilityStructure::ShiftStrikeSkew(IAsset* Asset, IInterpolator* VolShiftInterpolator, IInterpolator* VolSkewShiftInterpolator)
{
	begin_function
	map_bare_com_to_analytic(Asset, Asset, hAsset);
	map_bare_com_to_analytic(VolShiftInterpolator, Interpolator, hvolshift);
	map_bare_com_to_analytic(VolSkewShiftInterpolator, Interpolator, hskewshift);
	m_h->shiftStrikeSkew(hvolshift,hskewshift,hAsset);
	end_function
}

STDMETHODIMP CVolatilityStructure::ShiftStrikeInvariantSkew(IAsset* Asset, IInterpolator* VolShiftInterpolator, IInterpolator* VolSkewShiftInterpolator, VARIANT_BOOL IsRelative)
{
	begin_function
	map_bare_com_to_analytic(Asset, Asset, hAsset);
	map_bare_com_to_analytic(VolShiftInterpolator, Interpolator, hvolshift);
	map_bare_com_to_analytic(VolSkewShiftInterpolator, Interpolator, hskewshift);
	m_h->shiftStrikeInvariantSkew(hvolshift, hskewshift, hAsset, IsRelative ? true : false);
	end_function
}

STDMETHODIMP CVolatilityStructure::Stick(void)
{
	begin_function
	m_h->Stick();
	end_function
}