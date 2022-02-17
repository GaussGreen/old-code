//	handle_montecarlo.cpp : Implementation of CMonteCarlo
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "siriusapplication.h"
#include "handle_montecarlo.h"
#include "handle_volatilitystructure.h"
#include "handle_strikes.h"
#include "handle_asset.h"
#include "handle_correlationmatrix.h"
#include "handle_dividendschedule.h"
#include "handle_interpolator.h"
#include "handle_parameterlist.h"
#include "handle_spotschedule.h"
#include "handle_zerocurve.h"
#include "handle_assets.h"
#include "MlEqArray.h"
#include "MlEqPde.h"

/*static*/ const std::string					CMonteCarlo::s_szTimes = "Times";
/*static*/ const std::string					CMonteCarlo::s_szDiscountFactors = "Discount Factors";
/*static*/ const std::string					CMonteCarlo::s_szForwards = "Forwards";
/*static*/ const std::string					CMonteCarlo::s_szVolatilities = "Volatilities";
/*static*/ const std::string					CMonteCarlo::s_szBetas = "Betas";
/*static*/ const std::string					CMonteCarlo::s_szVolOfVols = "Volatility of Volatilities";
/*static*/ const std::string					CMonteCarlo::s_szMeanReversionRates = "Mean Reversion Rates";
/*static*/ const std::string					CMonteCarlo::s_szMeanReversionLevels = "Mean Reversion Levels";
/*static*/ const std::string					CMonteCarlo::s_szModelInfo = "Model Info";
/*static*/ const std::string					CMonteCarlo::s_sz_cL_Header = "cL";
/*static*/ const std::string					CMonteCarlo::s_sz_cR_Header = "cR";
/*static*/ const std::string					CMonteCarlo::s_sz_addTanhWings_Header = "addTanhWing";
/*static*/ const std::string					CMonteCarlo::s_sz_yPower_Header = "yPower";
/*static*/ const std::string					CMonteCarlo::s_sz_seed_Header = "seed";
/*static*/ const std::string					CMonteCarlo::s_sz_npaths_Header = "npaths";
/*static*/ const std::string					CMonteCarlo::s_sz_calibflag_Header = "calibflag";
/*static*/ const std::string					CMonteCarlo::s_sz_numberVolStates_Header = "numberVolStates";
/*static*/ const std::string					CMonteCarlo::s_sz_localVolFlag_Header = "localVolFlag";
/*static*/ const std::string					CMonteCarlo::s_sz_saveVolVolInfo_Header = "savevolvolInfo";
/*static*/ const std::string					CMonteCarlo::s_sz_numberGridPoints_Header = "numberGridPoints";
/*static*/ const std::string					CMonteCarlo::s_sz_saveVolGrid_Header = "saveVolGridFlag";
/*static*/ const std::string					CMonteCarlo::s_sz_randomNumberFlag_Header = "randomNumberFlag";
/*static*/ const std::string					CMonteCarlo::s_sz_contolVariateFlag_Header = "controlVariate";
/*static*/ const std::string					CMonteCarlo::s_sz_globalCalibFlag_Header = "globalCalib";
/*static*/ const std::string					CMonteCarlo::s_sz_UseHermite_Header = "useHermite";

STDMETHODIMP CMonteCarlo::get_Dates(long *pVal)
{
	*pVal = m_h->GetNumberOfFutureDates() ;
	return S_OK;
}

STDMETHODIMP CMonteCarlo::GetTime(long Date, double* pVal)
{
	if ( Date < 0 || Date > m_h->GetNumberOfFutureDates() ) throw CStringResource(IDS_ARRAY_INDEX);	

	*pVal = m_h->GetTimeAtPayoffDate(Date);
	if ( Date ){
		*pVal -= m_h->GetTimeAtPayoffDate(Date-1);
	}


//	*pVal = m_h->Getdt(Date);	
	return S_OK;
}


HRESULT CMonteCarlo::GetModelInfoParameter(CParameterMap& pm, const std::string& szHeading, long nRow, double fValue) const
{
	pm.SetValue(nRow, 0, szHeading);
	pm.SetValue(nRow, 1, fValue);
	return S_OK;
}

STDMETHODIMP CMonteCarlo::get_Paths(long *pVal)
{
	begin_function	
	*pVal = m_h->m_nPaths ;
	end_function
}

STDMETHODIMP CMonteCarlo::get_Assets(long *pVal)
{
	begin_function
	*pVal = m_h->m_nAssets;
	end_function
}
	
STDMETHODIMP CMonteCarlo::get_Value(VARIANT *pVal)
{	
	HRESULT								hr;
	MonteCarloTypeEnum					mct;

	// ToDo - finish
	if (hr = get_MonteCarloType(&mct)) return hr;
	if (mct == ForwardSkew){
		CForwardSkewMC* p = dynamic_cast<CForwardSkewMC*>(&*m_h);
						
		// ToDo - the first parameter should be the Monte Carlo type
		unmap_parameter(std::string("ToDo - Forward Volatility"), pmForwardVolatility);						// ToDo
		unmap_parameter(p->GetBetas(), pmBetas);
		unmap_parameter(p->GetVolVols(), pmVolVols);
		unmap_parameter(p->GetMeanReversionRates(), pmMeanReversionRates);
		unmap_parameter(p->GetMeanReversionLevels(), pmMeanReversionLevels);

		CParameterMap	pmModelParameters(14, 2);
		GetModelInfoParameter(pmModelParameters, s_sz_cL_Header, 0, -123456);								// ToDo
		GetModelInfoParameter(pmModelParameters, s_sz_cR_Header, 1, -123456);								// ToDo
		GetModelInfoParameter(pmModelParameters, s_sz_addTanhWings_Header, 2, -123456);						// ToDo
		GetModelInfoParameter(pmModelParameters, s_sz_yPower_Header, 3, -123456);							// ToDo
		GetModelInfoParameter(pmModelParameters, s_sz_seed_Header, 4, p->GetSeed());
		GetModelInfoParameter(pmModelParameters, s_sz_npaths_Header, 5, p->m_nPaths);
		GetModelInfoParameter(pmModelParameters, s_sz_calibflag_Header, 6, -123456);						// ToDo
		GetModelInfoParameter(pmModelParameters, s_sz_numberVolStates_Header, 7, p->GetNumberVolStates());	
		GetModelInfoParameter(pmModelParameters, s_sz_localVolFlag_Header, 8, p->GetLocalVolFlag());
		GetModelInfoParameter(pmModelParameters, s_sz_saveVolVolInfo_Header, 9, p->GetSaveVolofVolInfo());
		GetModelInfoParameter(pmModelParameters, s_sz_numberGridPoints_Header, 10, p->GetNumberGridPoints());
		GetModelInfoParameter(pmModelParameters, s_sz_saveVolGrid_Header, 11, p->GetSaveVolGrid());
		GetModelInfoParameter(pmModelParameters, s_sz_randomNumberFlag_Header, 12, p->GetRandomNumberFlag());
		GetModelInfoParameter(pmModelParameters, s_sz_contolVariateFlag_Header, 13, p->GetControlVariate());
		GetModelInfoParameter(pmModelParameters, s_sz_globalCalibFlag_Header, 14, p->GetGlobalCalib());
		GetModelInfoParameter(pmModelParameters, s_sz_UseHermite_Header, 15, -123456);						// ToDo

		unmap_parameter(std::string("ToDo - Strikes"), pmStrikes);											// ToDo	
		unmap_parameter(p->GetCalibTimeIndex(), pmCalibTimeIndex);
		unmap_parameter(p->GetCalibVols(), pmCalibVols);
		unmap_parameter(std::string("ToDo - CalibStrikes"), pmCalibStrikes);								// ToDo
		unmap_parameter(p->GetControlVariatePrices(), pmControlVariatePrices);		
		unmap_parameter(std::string("ToDo - ExcessVolVolInterpolatorOpt"), pmExcessVolVolInterpolator);		// ToDo
									
		if (dynamic_cast<CStochIRForwardSkewMC*>(p)){
			// i.e. StochasticIRMonteCarlo		
			unmap_parameter(std::string("ToDo - ZeroCurveOpt (only include this parameter if we have a stochastic IR Monte Carlo)"), pmZeroCurve);
			return CParameterMap::VariableArgumentListToArray(pVal, 13, pmForwardVolatility, pmBetas, pmVolVols, pmMeanReversionRates, pmMeanReversionLevels, pmModelParameters, pmStrikes, pmCalibTimeIndex, pmCalibVols, pmCalibStrikes, pmControlVariatePrices, pmExcessVolVolInterpolator, pmZeroCurve);
		} else {
			return CParameterMap::VariableArgumentListToArray(pVal, 12, pmForwardVolatility, pmBetas, pmVolVols, pmMeanReversionRates, pmMeanReversionLevels, pmModelParameters, pmStrikes, pmCalibTimeIndex, pmCalibVols, pmCalibStrikes, pmControlVariatePrices, pmExcessVolVolInterpolator);
		}
	} else {
		return CParameterMap::VectorToArray(m_vpm, pVal);
	}
}

STDMETHODIMP CMonteCarlo::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IMonteCarlo};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}
	
STDMETHODIMP CMonteCarlo::put_Value(VARIANT newVal)
{
	HRESULT								hr;
	std::vector<CParameterMap>			vpm;
	std::vector<CComVariant>			vv;
	
	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, &vpm)) return hr;
	map_parameter(vpm[0], long, MonteCarloType);
	
	switch (MonteCarloType){
	case ForwardSkew:
		if (hr = put_ValueForwardSkew(vpm, vv)) return hr;
		break;
	case General:
		if (hr = put_ValueGeneral(vpm, vv)) return hr;
		break;
	case Quasi:
		if (hr = put_ValueQuasi(vpm, vv)) return hr;
		break;
	case HybridForwardSkew:
		if (hr = put_ValueHybridForwardSkew(vpm, vv)) return hr;
		break;
	case CalibrateForwardSkew:
		if (hr = put_CalibrateMC(vpm, vv)) return hr;
		break;
	case CalibrateHybridModel:
		if (hr = put_CalibrateHybridMC(vpm,vv)) return hr;
		break;
	case LocalVolMonteCarlo:
		if (hr = put_ValueLocalVolMC(vpm,vv)) return hr;
		break;		
	default:
		return CParameterMap::ReturnErrorR(IDS_MONTE_CARLO_TYPE, IID_IMonteCarlo);
	}
	m_vpm = vpm;
	end_function
}

//	this helps set up the model information parameter used by CForwardSkewMC::Initialize
HRESULT CMonteCarlo::SetModelInfoParameter(const CParameterMap& pmIn, CVector* pv) const
{
	HRESULT								hr;
	CParameterMap						pmOut;

	if (hr = pmOut.SetSize(16, 1)) return hr;	
	if (hr = SetModelInfoParameter(pmIn, s_sz_cL_Header, 0, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_cR_Header, 1, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_addTanhWings_Header , 2, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_yPower_Header , 3, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_seed_Header, 4, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_npaths_Header, 5, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_calibflag_Header, 6, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_numberVolStates_Header, 7, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_localVolFlag_Header, 8, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_saveVolVolInfo_Header, 9, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_numberGridPoints_Header , 10, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_saveVolGrid_Header , 11, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_randomNumberFlag_Header , 12, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_contolVariateFlag_Header , 13, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_globalCalibFlag_Header, 14, &pmOut)) return hr;
	if (hr = SetModelInfoParameter(pmIn, s_sz_UseHermite_Header, 15, &pmOut)) return hr;
	return pmOut.GetValue(pv);	
}
HRESULT CMonteCarlo::SetModelInfoParameter(const CParameterMap& pm, const std::string& szHeading, long nRow, CParameterMap* ppm) const
{
	double								fVal;

	if (pm.GetValueFromRowHeader(szHeading, &fVal)) return CParameterMap::ReturnErrorRS(IDS_DATA_ELEMENT_NOT_FOUND, szHeading, IID_IMonteCarlo);
	return ppm->SetValue(nRow, 0, fVal);
}


STDMETHODIMP CMonteCarlo::Average(long Path, long DateFrom,long DateUpTo,long Asset, double* pVal)
{
	begin_function
	if (DateFrom > DateUpTo || Path <= 0 || Path > m_h->m_nPaths || DateUpTo < 0 || DateUpTo > m_h->GetNumberOfFutureDates() ||DateFrom < 0 || DateFrom > m_h->GetNumberOfFutureDates() || Asset <= 0 || Asset > m_h->m_nAssets ) throw CStringResource(IDS_ARRAY_INDEX);	
	double res(0);
	
	for(long i=DateFrom;i<=DateUpTo;i++)
	{
		res+=m_h->GetPathValue(Path - 1, i, Asset - 1);	
	}
	
	res	/= (DateUpTo-DateFrom+1);

	*pVal = res;

	end_function

}


STDMETHODIMP CMonteCarlo::Max(long Path, long DateFrom,long DateUpTo,long Asset, double* pVal)
{
	begin_function
	if (DateFrom > DateUpTo || Path <= 0 || Path > m_h->m_nPaths || DateUpTo < 0 || DateUpTo > m_h->GetNumberOfFutureDates() ||DateFrom < 0 || DateFrom > m_h->GetNumberOfFutureDates() || Asset <= 0 || Asset > m_h->m_nAssets ) throw CStringResource(IDS_ARRAY_INDEX);	
	double res(0);
	
	for(long i=DateFrom;i<=DateUpTo;i++)
		if(res<m_h->GetPathValue(Path - 1, i, Asset - 1))
			res=m_h->GetPathValue(Path - 1, i, Asset - 1);	
	
	
	*pVal = res;

	end_function

}

STDMETHODIMP CMonteCarlo::BasketMax(long Path, long DateFrom,long DateUpTo,IArray* Weights, double* pVal)
{
	begin_function


	// map Weights to a cvector

	
	map_bare_com_to_analytic(Weights, Array, hWeights)


	if (!hWeights) throw "You have not passed in an array of weights";



	if (DateFrom > DateUpTo || Path <= 0 || Path > m_h->m_nPaths || DateUpTo < 0 || DateUpTo > m_h->GetNumberOfFutureDates() ||DateFrom < 0 || DateFrom > m_h->GetNumberOfFutureDates()) throw CStringResource(IDS_ARRAY_INDEX);	
	double res(0);
	
	for(long i=DateFrom;i<=DateUpTo;i++)
	{
		
		double basketval(0);

		for(long iasset=0;iasset< m_h->m_nAssets;iasset++)
			basketval+=m_h->GetPathValue(Path - 1, i, iasset)* (*hWeights)[(int)iasset];


		if(res<basketval)
			res=basketval;
	
	}
	
	*pVal = res;

	end_function

}


STDMETHODIMP CMonteCarlo::Basket(long Path, long Date,IArray* Weights, double* pVal)
{
	begin_function


	// map Weights to a cvector
	
	map_bare_com_to_analytic(Weights, Array, hWeights)


	if (!hWeights) throw "You have not passed in an array of weights";

	if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 || Date > m_h->GetNumberOfFutureDates()) throw CStringResource(IDS_ARRAY_INDEX);	
	
	
	double basketval(0);

		for(long iasset=0;iasset< m_h->m_nAssets;iasset++)
			basketval+=m_h->GetPathValue(Path - 1, Date, iasset)* (*hWeights)[(int)iasset];
	
	
	*pVal = basketval;

	end_function

}

STDMETHODIMP CMonteCarlo::GetSpotValue(long Path, long Date, long Asset, double* pVal)
{
	begin_function		
	if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 || Date > m_h->GetNumberOfFutureDates() || Asset <= 0 || Asset > m_h->m_nAssets ) throw CStringResource(IDS_ARRAY_INDEX);	
	*pVal = m_h->GetPathValue(Path - 1, Date, Asset - 1);	
	end_function
}

STDMETHODIMP CMonteCarlo::PutSpotValue(long Path, long Date, long Asset, double newVal)
{
	begin_function		
	if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 || Date > m_h->GetNumberOfFutureDates() || Asset <= 0 || Asset > m_h->m_nAssets ) throw CStringResource(IDS_ARRAY_INDEX);
	m_h->PutPathValue(Path - 1, Date, Asset - 1, newVal);
	end_function
}

STDMETHODIMP CMonteCarlo::GetForwardValue(long Date, long Asset, double* pVal)
{
	begin_function
	if (Date < 0 || Date > m_h->GetNumberOfFutureDates() || Asset <= 0 || Asset > m_h->m_nAssets ) throw CStringResource(IDS_ARRAY_INDEX);	
	*pVal = m_h->GetForwardValue(Date, Asset - 1);	
	end_function
}


STDMETHODIMP CMonteCarlo::GetVolatility(long Path, long Date, double strike,long Asset, double* pVal)
{

	if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 || Date > m_h->GetNumberOfFutureDates() || Asset > m_h->m_nAssets || Asset <= 0 ) return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);


	try {
		*pVal = m_h->GetBridgeVolatility(Path -1, Date,strike,Asset-1 );		// this transforms 1-based to 0-based
	} catch (...){
		return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);
	}
	return S_OK;
}

STDMETHODIMP CMonteCarlo::GetImpliedLogContractVol( long Path, long Date,long Asset,long ngauss,double lowerStdev,double upperStdev, double* pVal)
{

	if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 || Date > m_h->GetNumberOfFutureDates() || Asset > m_h->m_nAssets || Asset <= 0 ) return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);


	try {
		*pVal = m_h->GetImpliedLogContractVol(Path -1, Date,Asset-1,ngauss,lowerStdev,upperStdev );		// this transforms 1-based to 0-based

	} catch (...){
		return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);
	}
	return S_OK;
}


STDMETHODIMP CMonteCarlo::GetDiscountFactor(long Path, long Date, double* pVal)
{
	begin_function			

	if ( Date <= 500 )
	{
//		Date is interpreted as a Monte Carlo time index

		if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 || Date > m_h->GetNumberOfFutureDates() ){
			throw CStringResource(IDS_ARRAY_INDEX);
		}
	}
	else
	{
//		Date is interpreted as a excel date

		if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 ){
			throw CStringResource(IDS_ARRAY_INDEX);
		}

	}


	*pVal = m_h->GetStochasticDiscount(Path - 1, Date);		// this transforms 1-based to 0-based

	end_function
}



STDMETHODIMP CMonteCarlo::GetZeroBond(long Path, long Date, long AsOfDate, double* pVal)
{
	begin_function			

	if ( Date <= 500 )
	{
//		Date is interpreted as a Monte Carlo time index

		if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 || Date > m_h->GetNumberOfFutureDates() ){
			throw CStringResource(IDS_ARRAY_INDEX);
		}
	}
	else
	{
//		Date is interpreted as a excel date

		if (Path <= 0 || Path > m_h->m_nPaths || Date < 0 ){
			throw CStringResource(IDS_ARRAY_INDEX);
		}

	}


	if (AsOfDate < 0 || AsOfDate >= m_h->GetNumberOfFutureDates()){
		return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);
	}

//		function is currently disabled
//		return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);

	if ( Date > 500 )
	{
		// we assume now that Date is an excel date
		long excelDate = Date;
		*pVal = m_h->GetDiscountToDate(Path-1,Date,AsOfDate);
	}
	else
	{
		// we assume now that Date is a Monte Carlo index that refers to a simulation date
		*pVal = m_h->GetDiscount(Path-1,AsOfDate,Date,AsOfDate);
	}


	end_function
}


STDMETHODIMP CMonteCarlo::GetDate(long Date, long* pVal)
{
	if ( Date < 0 || Date > m_h->GetNumberOfFutureDates() ) return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);
	try {
//		*pVal = CForwardSkewMC::GetPathValue(Path - 1, Date - 1);		// this transforms 1-based to 0-based
		*pVal = m_h->GetMCDate(Date);		// this transforms 1-based to 0-based

	} catch (...){
		return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);
	}
	return S_OK;
}


STDMETHODIMP CMonteCarlo::GetIndexFromDate(long Date, long* pVal)
{
	//if ( Date < 0 || Date > m_h->m_nDates ) return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);
	try {

		for(long i=0;i<=m_h->GetNumberOfFutureDates();i++)
		{
			if(m_h->GetMCDate(i)==Date)
			{
				*pVal=i; break;
			}
		
		}

//		*pVal = m_h->GetMCDate(Date);		// this transforms 1-based to 0-based

	} catch (...){
		return CParameterMap::ReturnErrorR(IDS_ARRAY_INDEX, IID_IMonteCarlo);
	}
	return S_OK;
}


STDMETHODIMP CMonteCarlo::get_MonteCarloType(MonteCarloTypeEnum *pVal)
{
	begin_function
	*pVal = NoMonteCarloType;

	std::string szName = m_h->GetName();

	if (szName == "CHermiteMC"){
		*pVal = Hermite;
	} else if (szName == "CQuasiMC"){
		*pVal = Quasi;
	} else if (szName == "CLocalVolMC" || szName == "CMultiLocalVolMC"){
		*pVal = LocalVolMonteCarlo;
	} else if (szName == "CStochIRForwardSkewMC"){
		*pVal = HybridForwardSkew;
	} else if (szName == "CForwardSkewMC"){
		*pVal = ForwardSkew;
	} else if (szName == "CMultiAssetForwardSkewMC"){
		*pVal = General;
	} else {
		throw "No Monte Carlo type associated with class '" + szName + "'";
	}	
	end_function
}

HRESULT CMonteCarlo::put_ValueForwardSkew(std::vector<CParameterMap>& vpm, std::vector<CComVariant>&	vv)
{
	/*parameter list is	vpm[0] - MonteCarloType (ForwardSkew in this case)
						vpm[1] - ForwardVolatility
						vpm[2] - Betas
						vpm[3] - VolVols
						vpm[4] - MeanReversionRates
						vpm[5] - MeanReversionLevels
						vpm[6] - ModelInfo
						vpm[7] - Strikes
						vpm[8] - CalibTimeIndex
						vpm[9] - CalibVolatilities
						vpm[10] - CalibStrikes
						vpm[11] - Spot
						vpm[12] - ControlVariatePricesOpt
						vpm[13] - ExcessVolVolInterpolatorOpt
						vpm[14] - ZeroCurveOpt*/
	
	begin_function
	
	HRESULT								hr;
	CVector								vModelInfo;	

	if (vpm.size() < 12 || vpm.size() > 15) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IMonteCarlo);
	
	map_object_parameter(vpm[1], VolatilityStructure, hVol);
	if (hVol->m_volType != ForwardVolatilities){
		throw "ForwardVolatility must be of type Forwardvol";
	}
	map_parameter(vpm[2], CVector, Betas);
	map_parameter(vpm[3], CVector, VolVols);
	map_parameter(vpm[4], CVector, MeanReversionRates);
	map_parameter(vpm[5], CVector, MeanReversionLevels);

	// model info
	if (hr = SetModelInfoParameter(vpm[6], &vModelInfo)) return hr;
	
	// strikes	
	map_object_vector_parameter(vv[7], Strikes, apStrike)

	std::vector<std::vector<MlEqStrikeHandle> >		pStrikes;			// [idate][istrike]
	pStrikes.resize(apStrike.size());
	for (int i = 0; i < apStrike.size(); i++){
		pStrikes[i] = apStrike[i]->m_hStrikes;
	}	

	const vector<DATE>& simDates = hVol->getDates();
	if (!simDates.size()) throw "No simulation dates";

	CMatrix vols;
	int which;
	if ( pStrikes[0][0]->m_strikeType != Normalised) return E_FAIL;
	
	int nstrikes = apStrike[0]->m_hStrikes.size(); 
	if ( nstrikes == 0 ) return E_FAIL;

	vols.resize(simDates.size()+1,nstrikes);
	for (int j = 0 ; j < nstrikes; j++ ){
		vols[0][j] = pStrikes[0][j]->m_strike;
	}

	for ( int i = 0 ; i < simDates.size(); i++ )
	{
		if (  apStrike.size() > 1 )
		 {
			int istrike = apStrike[i]->m_hStrikes.size();
			if ( istrike != nstrikes ){
				return E_FAIL;
			}
		 }

		 if (  apStrike.size() > 1 ){
			  which = i;}else{
			  which = 0;
		 }

		 for ( int j = 0 ; j < nstrikes; j++ )
		 {
			 if ( pStrikes[which][0]->m_strikeType != Normalised) return E_FAIL;
			 if ( ( vols[i+1][j] = hVol->getVol(*pStrikes[which][j],i, Middle) )   < 0 ){
					return E_FAIL;
			 }
		  }
	  }
													
//	retrieve timeindex for calibration	
	map_parameter(vpm[8], CMatrix, calibTimeIndexRaw);
	int ndim=0;
	for ( int i = 0 ; i < calibTimeIndexRaw.rows(); i++ )
	{
		if ( calibTimeIndexRaw[i][0] ){
			ndim++;
		}
	}

	CMatrix calibTimeIndex(ndim,calibTimeIndexRaw.cols());
	int k=0;
	for (int i = 0 ; i < calibTimeIndexRaw.rows(); i++ )
	{
		if ( calibTimeIndexRaw[i][0] )
		{
			calibTimeIndex[k][0] = i+1;
			for ( int j = 1 ; j < calibTimeIndexRaw.cols(); j++ ){
				calibTimeIndex[k][j] = calibTimeIndexRaw[i][j];
			}
			k++;
		}
	}


		//		retrieve marketvols
	map_object_parameter(vpm[9], VolatilityStructure, hCalibVol);

//	retrieve calibStrikes and vols
	map_object_vector_parameter(vv[10], Strikes, apCalibStrike)

	std::vector<std::vector<MlEqStrikeHandle> >		pCalibStrikes;			//[idate][istrike]
	pCalibStrikes.resize(apCalibStrike.size());
	for (int i = 0; i < apCalibStrike.size(); i++){
		pCalibStrikes[i] = apCalibStrike[i]->m_hStrikes;
	}	

	CMatrix calibVols,calibStrikes;	
	if (pCalibStrikes[0][0]->m_strikeType != Normalised) return E_FAIL;
		
	nstrikes = apCalibStrike[0]->m_hStrikes.size(); 
	if ( nstrikes == 0 ) return E_FAIL;

	calibStrikes.resize(calibTimeIndex.rows(),nstrikes);
	calibVols.resize(calibTimeIndex.rows(),nstrikes);

	for (int i = 0 ; i < calibTimeIndex.rows(); i++ )
	 {
		if (  apCalibStrike.size() > 1 )
		{
			int istrike = apCalibStrike[i]->m_hStrikes.size();
			if ( istrike != nstrikes ){
				return E_FAIL;
			}
		}

		if (apCalibStrike.size() > 1 ){
			which = calibTimeIndex[i][0]-1;
		}else{
			which = 0;
		}

		for (int j = 0 ; j < nstrikes; j++ )
		{
			if ( 	calibVols.rows() > i ){
				calibStrikes[i][j] = pCalibStrikes[which][j]->m_strike;
			}

			if ( pCalibStrikes[which][0]->m_strikeType != Normalised) return E_FAIL;
			if ( ( calibVols[i][j] = hCalibVol->getVol(*pCalibStrikes[which][j],(int)(which+1e-11), Middle) )   < 0 ){
				return E_FAIL;
			}
		}
	}

	map_optional_parameter(vpm.size() > 12 ? vv[12] : CComVariant(), CMatrix, controlVarPrices, CMatrix());
	
	
	map_optional_object_parameter(vpm.size() > 13 ? vv[13] : CComVariant(), Interpolator, hInterpolator);	
	MlEq2DInterpolatorHandle hJim;
	
	
	if (!!hInterpolator){
		if (!(hJim = dynamic_cast<MlEq2DInterpolator*>(&*hInterpolator))){
			throw "You must pass in a 2 dimensional interpolator for the excess vol of vol";
		}
	}


	map_parameter(vpm[11], double, spot);
//	this is really good



	vector<long> simdates(simDates.size());
	for (int i = 0 ; i < simdates.size(); i++ ){
		simdates[i] = simDates[i];
	}
		
	map_optional_object_parameter(vv.size() > 13 ? vv[14] : CComVariant(), ZeroCurve, hZeroCurve);

	if (!hZeroCurve){				
		CForwardSkewMC* pMC;		
		m_h = pMC = new CForwardSkewMC(hVol->getDateToDouble());				
		pMC->Initialize(spot, simdates, 
				hVol->getForwards(), 
				hVol->getDiscounts(),
				vols, 
				Betas, 
				VolVols, 
				MeanReversionRates, 
				MeanReversionLevels, 
				vModelInfo, 
				calibTimeIndex,calibVols,pCalibStrikes,controlVarPrices,
				hJim);
	} else {
		m_h = new CStochIRForwardSkewMC(hVol->getDateToDouble(),spot,
				simdates, 
				hVol->getForwards(), 
				hVol->getDiscounts(),
				vols, 
				Betas, 
				VolVols, 
				MeanReversionRates, 
				MeanReversionLevels, 
				vModelInfo, 
				calibTimeIndex,
				calibVols,
				pCalibStrikes,
				controlVarPrices,				
				hZeroCurve,
				hJim);
	}

	m_h->generatePaths();	
	end_function
}







HRESULT CMonteCarlo::put_ValueGeneral(std::vector<CParameterMap>& vpm, std::vector<CComVariant>&	vv)
{
	/*parameter list is	vpm[0] - MonteCarloType (General in this case)
						vpm[1] - Assets
						vpm[2] - Dates
						vpm[3] - ModelParameters
						vpm[4] - HullAndWhiteModels*/

	HRESULT								hr;	
	CVector								vModelInfo;
	
	if (vpm.size() != 5 && vpm.size() != 4) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IMonteCarlo);
	vpm.resize(5);
	map_object_vector_parameter(vv[1], Asset, Assets);
	map_parameter(vpm[2], std::vector<long>, Dates);
	if( Dates.size() == 0 ){
		throw "You must pass dates in";
	}			
	if (hr = SetModelInfoParameter(vpm[3], &vModelInfo)) return hr;		
	map_optional_object_parameter(vpm[4], ParameterList, hParameterList);

	GVector<int> rateIsStochastic(Assets.size());
	for (int n = 0; n < Assets.size(); n++){
		rateIsStochastic[n] = 0;
	}

	// In setting up the Monte Carlo, the analytics will require mutual asset correlations.
	// We need to pass these to the analytics. Take the date and data source from the first
	// asset in the input array.
	CAssets::InsertCorrelations(Assets);

	if ( Assets.size() > 1){
		m_h = new CMultiAssetForwardSkewMC(Assets[0]->GetDateHandle(), Assets, Dates, rateIsStochastic, vModelInfo, CMatrix(), CMatrix(), std::vector<std::vector<MlEqStrikeHandle> >());
	} else {
		m_h = new CForwardSkewMC(Assets[0]->GetDateHandle(), *Assets[0], Dates, vModelInfo, CMatrix(), CMatrix(), std::vector<std::vector<MlEqStrikeHandle> >());		
	}
	m_h->generatePaths();	
	return S_OK;
}


HRESULT CMonteCarlo::put_ValueLocalVolMC(std::vector<CParameterMap>& vpm, std::vector<CComVariant>&	vv)
{
	/*parameter list is	vpm[0] - MonteCarloType (LocalVol in this case)
						vpm[1] - Asset
						vpm[2] - Dates
						vpm[3] - ModelParameters, that is number of path, steps and rng flag 
	*/

	CVector			vModelInfo(4);

	if (vpm.size() != 4) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IMonteCarlo);
	vpm.resize(4);
	map_object_vector_parameter(vv[1], Asset, Assets);
	map_parameter(vpm[2], std::vector<long>, Dates);

	if( Dates.size() == 0 ){
		throw "You must pass dates in";
	}

	map_object_parameter(vpm[3], ParameterList, ModelParameters);
	
	ModelParameters->BeginExhaust();
	ModelParameters->GetValue("npaths = 50000; numberStepsAYear = 50; randomNumberFlag = 6", ";", "=", &vModelInfo[0]);
	ModelParameters->CheckUnexhausted();

	CAssets::InsertCorrelations(Assets);
	
	int npaths	= vModelInfo[0];
	int nSteps	= vModelInfo[1];
	int rngFlag = vModelInfo[2];
		
	if ( Assets.size() == 1 && !Assets[0]->IsBasket() )
	{
		CLocalVolMC* pMC;
		m_h = pMC = new CLocalVolMC(Assets[0]->GetDateHandle(), Assets[0]);
		pMC->initialize(Dates, npaths, rngFlag, nSteps);		
	}
	else
	{
		CLocalVolMC* pMC;
		if (Assets.size() > 1){
			m_h = pMC = new CMultiLocalVolMC(Assets[0]->GetDateHandle(), Assets);
		} else if (Assets[0]->IsBasket()){
			m_h = pMC = new CMultiLocalVolMC(Assets[0]->GetDateHandle(), Assets[0]);
		}

		pMC->CLocalVolMC::initialize(Dates, npaths, rngFlag, nSteps);
	} 
	
	m_h->generatePaths();	

	return S_OK;
}

HRESULT CMonteCarlo::put_ValueQuasi(std::vector<CParameterMap>& vpm, std::vector<CComVariant>& vv)
{
	/*parameter list is	0 - MonteCarloType (General in this case)
						1 - SpotSchedulesArr
						2 - Dates*/

	if (vpm.size() != 3) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IMonteCarlo);
	map_object_vector_parameter(vv[1], SpotSchedule, SpotSchedulesArr);
	map_parameter(vpm[2], std::vector<long>, Dates);	
	if( Dates.size() == 0 ){
		throw "You must pass dates in";
	}
	
	CQuasiMC* pMC;	
	m_h = pMC = new CQuasiMC(new MlEqDate(MlEqDate::GetCurrentDate()));
	pMC->Initialize(SpotSchedulesArr, Dates);
	return S_OK;
}

HRESULT CMonteCarlo::put_CalibrateMC(std::vector<CParameterMap>& vpm, std::vector<CComVariant>&	vv)
{

	/*parameter list is	vpm[0] - MonteCarloType (General in this case)
						vpm[1] - Assets
						vpm[2] - Dates
						vpm[3] - ModelParameters
						vpm[4] - calibTimeIndexRaw
						vpm[5] - CalibStrike
	*/


	HRESULT								hr;
	CVector								vModelInfo;
	
	if (vpm.size() != 6 ) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IMonteCarlo);

//	vpm.resize(6);

	map_object_parameter(vv[1], Asset, Assets);
	map_parameter(vpm[2], GVector<long>, Dates);
	if (hr = SetModelInfoParameter(vpm[3], &vModelInfo)) return hr;
	map_parameter(vpm[4], CMatrix, calibTimeIndexRaw);
	
	std::vector<long> mcDates(Dates.getsize());
	for ( int i= 0 ; i < Dates.getsize(); i++ ){
		mcDates[i] = Dates[i];
	}

	int ndim=0;
	for ( int i = 0 ; i < calibTimeIndexRaw.rows(); i++ )
	{
		if ( calibTimeIndexRaw[i][0] ){
			ndim++;
		}
	}

	CMatrix calibTimeIndex(ndim,calibTimeIndexRaw.cols());
	int k=0;
	for (int i = 0 ; i < calibTimeIndexRaw.rows(); i++ )
	{
		if ( calibTimeIndexRaw[i][0] )
		{
			calibTimeIndex[k][0] = i+1;
			for ( int j = 1 ; j < calibTimeIndexRaw.cols(); j++ ){
				calibTimeIndex[k][j] = calibTimeIndexRaw[i][j];
			}
			k++;
		}
	}



//	retrieve calibStrikes and vols
	map_object_vector_parameter(vv[5], Strikes, apCalibStrike)

	std::vector<std::vector<MlEqStrikeHandle> >		pCalibStrikes;			//[idate][istrike]
	pCalibStrikes.resize(apCalibStrike.size());
	for (int i = 0; i < apCalibStrike.size(); i++){
		pCalibStrikes[i] = apCalibStrike[i]->m_hStrikes;
	}	

	CMatrix calibVols,calibStrikes;	
	if ( pCalibStrikes[0][0]->m_strikeType != Normalised) return E_FAIL;
		
	int nstrikes = apCalibStrike[0]->m_hStrikes.size(); 
	if ( nstrikes == 0 ) return E_FAIL;

	calibStrikes.resize(calibTimeIndex.rows(),nstrikes);
	calibVols.resize(calibTimeIndex.rows(),nstrikes);

	int which;
	for (int i = 0 ; i < calibTimeIndex.rows(); i++ )
	 {
		if (  apCalibStrike.size() > 1 )
		{
			int istrike = apCalibStrike[i]->m_hStrikes.size();
			if ( istrike != nstrikes ){
				return E_FAIL;
			}
		}

		if (apCalibStrike.size() > 1 ){
			which = calibTimeIndex[i][0]-1;
		}else{
			which = 0;
		}

		for (int j = 0 ; j < nstrikes; j++ )
		{
			if ( 	calibVols.rows() > i ){
				calibStrikes[i][j] = pCalibStrikes[which][j]->m_strike;
			}

			if ( pCalibStrikes[which][0]->m_strikeType != Normalised) return E_FAIL;
			if ( ( calibVols[i][j] = Assets->GetCompositeVolatility(*pCalibStrikes[which][j],mcDates[0],(double)mcDates[which+1], Middle) ) < 0 ){
				return E_FAIL;
			}
		}
	}

	m_h = new CForwardSkewMC(Assets->GetDateHandle(), *Assets, mcDates, vModelInfo, calibTimeIndex, calibVols, pCalibStrikes);	
	m_h->generatePaths();
	return S_OK;
}


HRESULT CMonteCarlo::put_CalibrateHybridMC(std::vector<CParameterMap>& vpm, std::vector<CComVariant>& vv)
{

	/*parameter list is	vpm[0] - MonteCarloType (HybridForwardSkew in this case)
						vpm[1] - Asset
						vpm[2] - Dates
						vpm[3] - ModelParameters
						vpm[4] - HullAndWhiteModels
						vpm[5] - Correlation (scalar)
						vpm[6] - isMinimalHybrid (scalar)
						vpm[7] - calibTimeIndexRaw
						vpm[8] - CalibStrike
						*/

	HRESULT								hr;
	CVector								vModelInfo;

	if (vpm.size() != 9) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IMonteCarlo);

	vpm.resize(9);
	map_object_parameter(vv[1], Asset, Asset);
	map_parameter(vpm[2], GVector<long>, Dates);
	if (hr = SetModelInfoParameter(vpm[3], &vModelInfo)) return hr;
	map_object_vector_parameter(vv[4], HullAndWhite, hHullAndWhites);
	map_parameter(vpm[5], double, fCorrelation);
	map_parameter(vpm[6], bool, isMinimalHybrid);

	if ( hHullAndWhites.size() == 0 ){
		throw("no HullWhite models have been entered");
	}

	// hHullAndWhites is std::vectorMlEqHullAndWhiteHandle>

	std::vector<long> mcDates(Dates.getsize());
	for ( int i= 0 ; i < Dates.getsize(); i++ ){
		mcDates[i] = Dates[i];
	}


	std::string	szPayCcy = Asset->GetPayZeroCurve(true)->GetName();
	std::string	szNaturalCcy = Asset->GetNaturalZeroCurve(true)->GetName();

	MlEqHullAndWhiteHandle spHWNatural, spHWPay;

	if ( hHullAndWhites.size() >= 1 )
	{
		MlEqZeroCurveHandle pCurve = hHullAndWhites[0]->getYieldCurve();
		std::string curr = 	pCurve->GetName();

		if ( ( curr != szPayCcy ) && ( curr != szNaturalCcy ) ){
			throw("first HullWhite model entered is inconsistent with pay and natural currency of the asset");
		}

		if ( curr == szPayCcy ){
			spHWPay = hHullAndWhites[0];
		}

		if ( curr == szNaturalCcy ){
			spHWNatural = hHullAndWhites[0];
		}


		if ( hHullAndWhites.size() == 2 )
		{
			pCurve = hHullAndWhites[1]->getYieldCurve();
			curr = 	pCurve->GetName();

			if ( ( curr != szPayCcy ) && ( curr != szNaturalCcy ) ){
				throw("second HullWhite model entered is inconsistent with pay and natural currency of the asset");
			}

			if ( curr == szPayCcy ){
				spHWPay = hHullAndWhites[1];
			}
			if ( curr == szNaturalCcy ){
				spHWNatural = hHullAndWhites[1];
			}

		}
		else if ( hHullAndWhites.size() > 2 ) {
			throw("too many HullWhite models have been entered");
		}

	}


	map_parameter(vpm[7], CMatrix, calibTimeIndexRaw);
	

	int ndim=0;
	for ( int i = 0 ; i < calibTimeIndexRaw.rows(); i++ )
	{
		if ( calibTimeIndexRaw[i][0] ){
			ndim++;
		}
	}

	CMatrix calibTimeIndex(ndim,calibTimeIndexRaw.cols());
	int k=0;
	for (int i = 0 ; i < calibTimeIndexRaw.rows(); i++ )
	{
		if ( calibTimeIndexRaw[i][0] )
		{
			calibTimeIndex[k][0] = i+1;
			for ( int j = 1 ; j < calibTimeIndexRaw.cols(); j++ ){
				calibTimeIndex[k][j] = calibTimeIndexRaw[i][j];
			}
			k++;
		}
	}


//	retrieve calibStrikes and vols
	map_object_vector_parameter(vv[8], Strikes, apCalibStrike)

	std::vector<std::vector<MlEqStrikeHandle> >		pCalibStrikes;			//[idate][istrike]
	pCalibStrikes.resize(apCalibStrike.size());
	for (int i = 0; i < apCalibStrike.size(); i++){
		pCalibStrikes[i] = apCalibStrike[i]->m_hStrikes;
	}	

	CMatrix calibVols,calibStrikes;	
	if ( pCalibStrikes[0][0]->m_strikeType != Normalised) return E_FAIL;
		
	int nstrikes = apCalibStrike[0]->m_hStrikes.size(); 
	if ( nstrikes == 0 ) return E_FAIL;

	calibStrikes.resize(calibTimeIndex.rows(),nstrikes);
	calibVols.resize(calibTimeIndex.rows(),nstrikes);

	int which;
	for (int i = 0 ; i < calibTimeIndex.rows(); i++ )
	 {
		if (  apCalibStrike.size() > 1 )
		{
			int istrike = apCalibStrike[i]->m_hStrikes.size();
			if ( istrike != nstrikes ){
				return E_FAIL;
			}
		}

		if (apCalibStrike.size() > 1 ){
			which = calibTimeIndex[i][0]-1;
		}else{
			which = 0;
		}

		for (int j = 0 ; j < nstrikes; j++ )
		{
			if ( 	calibVols.rows() > i ){
				calibStrikes[i][j] = pCalibStrikes[which][j]->m_strike;
			}

			if ( pCalibStrikes[which][0]->m_strikeType != Normalised) return E_FAIL;
			if ( ( calibVols[i][j] = Asset->GetCompositeVolatility(*pCalibStrikes[which][j],mcDates[0],(double)mcDates[which+1], Middle) )   < 0 ){
				return E_FAIL;
			}
		}
	}


	CStochIRForwardSkewMC* pMC;
	m_h = pMC = new CStochIRForwardSkewMC(Asset->GetDateHandle());
	pMC->Initialize(*Asset,mcDates,vModelInfo,calibTimeIndex,calibVols,pCalibStrikes, spHWNatural, spHWPay, fCorrelation,isMinimalHybrid);	
	m_h->generatePaths();		
	return S_OK;
}


HRESULT CMonteCarlo::put_ValueHybridForwardSkew(std::vector<CParameterMap>& vpm, std::vector<CComVariant>& vv)
{
	/*parameter list is	vpm[0] - MonteCarloType (HybridForwardSkew in this case)
						vpm[1] - Asset
						vpm[2] - Dates
						vpm[3] - ModelParameters
						vpm[4] - HullAndWhiteModels
						vpm[5] - Correlation (scalar)
						vpm[5] - isMinimalHybrid (scalar)*/

	HRESULT								hr;
	CVector								vModelInfo;

	if (vpm.size() != 7) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IMonteCarlo);

	vpm.resize(7);
	map_object_parameter(vv[1], Asset, Asset);
	map_parameter(vpm[2], GVector<long>, Dates);
	if (hr = SetModelInfoParameter(vpm[3], &vModelInfo)) return hr;
	map_object_vector_parameter(vv[4], HullAndWhite, hHullAndWhites);
	map_parameter(vpm[5], double, fCorrelation);
	map_parameter(vpm[6], bool, isMinimalHybrid);

	if ( hHullAndWhites.size() == 0 ){
		throw("no HullWhite models have been entered");
	}

	// hHullAndWhites is std::vector<MlEqHullAndWhiteHandle>

	std::vector<long> mcDates(Dates.getsize());
	for ( int i= 0 ; i < Dates.getsize(); i++ ){
		mcDates[i] = Dates[i];
	}

	CMatrix calibTimeIndex,calibVols;
	std::vector<std::vector<MlEqStrikeHandle> > pCalibStrikes;

	std::string	szPayCcy = Asset->GetPayZeroCurve(true)->GetName();
	std::string	szNaturalCcy = Asset->GetNaturalZeroCurve(true)->GetName();

	MlEqHullAndWhiteHandle spHWNatural, spHWPay;

	if ( hHullAndWhites.size() >= 1 )
	{
		MlEqZeroCurveHandle pCurve = hHullAndWhites[0]->getYieldCurve();
		std::string curr = 	pCurve->GetName();

		if ( ( curr != szPayCcy ) && ( curr != szNaturalCcy ) ){
			throw("first HullWhite model entered is inconsistent with pay and natural currency of the asset");
		}

		if ( curr == szPayCcy ){
			spHWPay = hHullAndWhites[0];
		}

		if ( curr == szNaturalCcy ){
			spHWNatural = hHullAndWhites[0];
		}


		if ( hHullAndWhites.size() == 2 )
		{
			pCurve = hHullAndWhites[1]->getYieldCurve();
			curr = 	pCurve->GetName();

			if ( ( curr != szPayCcy ) && ( curr != szNaturalCcy ) ){
				throw("second HullWhite model entered is inconsistent with pay and natural currency of the asset");
			}

			if ( curr == szPayCcy ){
				spHWPay = hHullAndWhites[1];
			}
			if ( curr == szNaturalCcy ){
				spHWNatural = hHullAndWhites[1];
			}

		}
		else if ( hHullAndWhites.size() > 2 ) {
			throw("too many HullWhite models have been entered");
		}

	}
	
	CStochIRForwardSkewMC* pMC;
	m_h = pMC = new CStochIRForwardSkewMC(Asset->GetDateHandle());
	pMC->Initialize(*Asset,mcDates, vModelInfo,calibTimeIndex,calibVols,pCalibStrikes, spHWNatural, spHWPay, fCorrelation,isMinimalHybrid);	
	m_h->generatePaths();	

/*
//  start tests:

	int ndates = m_h->m_spIRmc->m_nDates;

	CVector discounts(ndates);

	for ( int ipath = 0 ; ipath < m_h->m_spIRmc->m_nPaths; ipath++ )
	{
		double dsk;

		for ( int idate = 0 ; idate < ndates; idate++ ){

			dsk = m_h->m_spIRmc->GetStochasticDiscount(ipath,idate);
			discounts[idate] += dsk;
		}
	}


	for ( int idate = 0 ; idate < ndates; idate++ ){
		discounts[idate] /= m_h->m_spIRmc->m_nPaths;
	}
*/
	
	return S_OK;
}


