//	CRisk.cpp: implementation of the CRisk class.
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "siriusrisk.h"
#include "risk.h"


CComPtr<IResult> CRisk::GetResultObject(IResult*** pppResult) const
{
	CComPtr<IResult>					spResult;
	
	if (!**pppResult){
		// Create a result object
		spResult.CoCreateInstance(CLSID_Result);
		spResult.CopyTo(*pppResult);
	} else {
		spResult = **pppResult;		
	}
	return spResult;
}

bool CRisk::IsCached(const estring& szTag, CComPtr<IResult> spResult, double* pfOut) const
{				
	CComVariant							v;
		
	if (!spResult->GetValue(szTag.GetBSTR(), &v) && !v.ChangeType(VT_R8)){
		*pfOut = v.dblVal;
		return true;
	}
	return false;
}


double CRisk::Delta(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult, bool bReinitialiseAll, std::string szUnderlyingFilter  ) const
{
	double fUnderPrice =  UnderlyingPrice(Portfolio, spResult);
	double f = Price(Portfolio, spResult);
	double fHigh = PriceSpot(fShift, Portfolio, spResult,bReinitialiseAll, szUnderlyingFilter);	
	double fDelta =  (fHigh - f) / ( fShift *  fUnderPrice );

	estring szDeltaKey("Delta");

	if (szUnderlyingFilter.size()) {
		szDeltaKey +="=" +szUnderlyingFilter; 
	}

	spResult->AddValue(szDeltaKey.GetBSTR(), CComVariant(fDelta));
	return fDelta;

}

double CRisk::PriceDiv(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult) const
{

	HRESULT								hr = S_OK;
	estring								szKey;
	double								fPriceDiv= 0.0;
	szKey = estring("SpotDiv:") + estring(fShift) ;

	if (!IsCached(szKey , spResult, &fPriceDiv )){
		CComPtr<IAssets>				spAssets;
		CComPtr<IAssets>				spBaseAssets;
		DATE							date;
		std::vector<double>				afCurrentSpots;
		long							nAssets;
		
		// Perturb the market data
		hr = Portfolio->GetUnderlyings(&spAssets);
		hr = spAssets->GetBaseUnderlyings(&spBaseAssets);

		std::map<estring,CComVariant> mapDividend;

		if (Portfolio->get_Date(&date)) propagate_error;	
		spBaseAssets->get_Count(&nAssets);
		afCurrentSpots.resize(nAssets);
		for (long nAsset = 1; nAsset <= nAssets; nAsset++){
			CComPtr<IAsset>						spAsset;
			CComBSTR							szDividendName;
			CComVariant							vDividendData;
			CComPtr<IDividendSchedule>			spDividendSchedule;
		
			hr = spBaseAssets->get_Item(CComVariant(nAsset), &spAsset);
			hr = spAsset->get_DividendSchedule(&spDividendSchedule);

			if (spDividendSchedule ) {
				// deep copy 
				hr = spDividendSchedule->get_Name(&szDividendName);
				hr = spDividendSchedule->get_Value(&vDividendData);
				mapDividend.insert(std::pair<estring,CComVariant>(estring(szDividendName),vDividendData));
				hr = spDividendSchedule->DiscreteShift(fShift ,L"1Y");		
			}
			
		}

		// Perform the evaluation
		CComPtr<IResult>					sp;
		sp.CoCreateInstance(CLSID_Result);
		if (Portfolio->Evaluate(L"Price", &sp) || sp->GetPrice(&fPriceDiv)) propagate_error;
		spResult->AddValue(szKey.GetBSTR(), CComVariant(fPriceDiv));
				
		// Reset the market data.
		for (long nAsset = 1; nAsset <= nAssets; nAsset++){
			CComPtr<IAsset>			spAsset;					
			CComVariant							vDividendData;
			CComPtr<IDividendSchedule>			spDividendSchedule;
			CComBSTR							szDividendName;
			estring								szName;

			hr = spBaseAssets->get_Item(CComVariant(nAsset), &spAsset);
			hr = spAsset->get_DividendSchedule(&spDividendSchedule);

			if (spDividendSchedule  ==  NULL )  continue; 
			hr = spDividendSchedule->get_Name(&szDividendName);
			szName = estring(szDividendName);
			std::map<estring,CComVariant>::iterator it =  mapDividend.find(szName)  ;
			if ( it != mapDividend.end() ) {
				hr = spDividendSchedule->put_Value(it->second);

				hr = spAsset->Refresh(VARIANT_FALSE);


			} else {
				throw "unable to reset " + szName +" Dividend Schedule";
			}						
		}
				
	}
	
	return fPriceDiv;


}



double CRisk::PriceIR(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult) const
{

	HRESULT								hr = S_OK;
	estring								szKey;
	double								fPriceIR = 0.0;
	szKey = estring("SpotIR:") + estring(fShift) ;

	if (!IsCached(szKey , spResult, &fPriceIR )){
		// Perform the evaluation
		CComPtr<IResult>					sp;
		long								nZerocurves = 0L;
		CComPtr<IZeroCurves>					spZeroCurves;
	
		sp.CoCreateInstance(CLSID_Result);
		hr  = Portfolio->GetZeroCurves(&spZeroCurves);

		//ZeroCurves shift 
		hr = spZeroCurves->get_Count(&nZerocurves);
		hr = spZeroCurves->Shift(fShift);
		

		if (Portfolio->Evaluate(L"Price", &sp) || sp->GetPrice(&fPriceIR)) propagate_error;

		spResult->AddValue(szKey.GetBSTR(), CComVariant(fPriceIR));
		//result 

		//Reset you need to shift  back down because you might be in a sceario where
		//a shift has already been applied
		hr = spZeroCurves->Shift(-fShift);;
		
	}
	
	return fPriceIR ;



}



double CRisk::DivDelta(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult ) const
{
	double f = Price(Portfolio, spResult);
	double fDivShift = PriceDiv(fShift , Portfolio, spResult); ;
	double fDivDelta =fDivShift - f ;

	spResult->AddValue(L"DivDelta", CComVariant(fDivDelta));

	return fDivDelta ;
}

double CRisk::Rho(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult) const
{
	double f = Price(Portfolio, spResult);
	double fIRShift = PriceIR(fShift , Portfolio, spResult); ;
	double fRho =fIRShift - f ;

	spResult->AddValue(L"Rho", CComVariant(fRho));

	return fRho ;
}



double CRisk::Price(IEvaluatable* Portfolio, CComPtr<IResult> spResult) const
{
	double								fPrice;		

	if (!IsCached("Price", spResult, &fPrice)){
		// Perform the evaluation
		CComPtr<IResult>					sp;
		sp.CoCreateInstance(CLSID_Result);
		if (Portfolio->Evaluate(L"Price", &sp) || sp->GetPrice(&fPrice)) propagate_error;
		spResult->AddValue(L"Price", CComVariant(fPrice));
	}
	
	return fPrice;
}

double CRisk::PriceSpot(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult, bool bReinitialiseAll, std::string szUnderlyingFilter) const
//	fShift - additive shift such that the spot (S) is transformed to S * (1.0 + fShift)
{
	double								fPrice;		
	HRESULT								hr = S_OK;
	estring								szKey;

	szKey = estring("Spot:") + estring(fShift) + estring(":") + estring(bReinitialiseAll) + estring(":");
	if(szUnderlyingFilter.size()) {
		szKey +=":"+ szUnderlyingFilter;

	}

		
	 
	if (!IsCached( szKey , spResult, &fPrice)){
		CComPtr<IAssets>				spAssets;
		CComPtr<IAssets>				spBaseAssets;
		DATE							date;
		std::vector<double>				afCurrentSpots;
		long							nAssets;
		
		// Perturb the market data
		hr = Portfolio->GetUnderlyings(&spAssets);
		hr = spAssets->GetBaseUnderlyings(&spBaseAssets);

		if (Portfolio->get_Date(&date)) propagate_error;	
		spBaseAssets->get_Count(&nAssets);
		afCurrentSpots.resize(nAssets);
		for (long nAsset = 1; nAsset <= nAssets; nAsset++){
			CComPtr<IAsset>			spAsset;
			double					fSpot;			
			CComBSTR				sName;
			spBaseAssets->get_Item(CComVariant(nAsset), &spAsset);
			spAsset->get_Name(&sName);
			spAsset->get_NaturalSpot(date, &fSpot);
			afCurrentSpots[nAsset - 1] = fSpot;		
			if (szUnderlyingFilter.size() == 0  || estring::CompareNoCase(szUnderlyingFilter, sName) == 0 )  { 
				spAsset->put_NaturalSpot(date, fSpot * (1.0 + fShift));
				ReinitialiseStrike(spAsset,bReinitialiseAll);
			}
		}

		// Perform the evaluation
		CComPtr<IResult>					sp;
		sp.CoCreateInstance(CLSID_Result);
		if (Portfolio->Evaluate(L"Price", &sp) || sp->GetPrice(&fPrice)) propagate_error;
		spResult->AddValue(szKey.GetBSTR(), CComVariant(fPrice));
				
		// Reset the market data.
		for (long nAsset = 1; nAsset <= nAssets; nAsset++){
			CComPtr<IAsset>			spAsset;						
			spBaseAssets->get_Item(CComVariant(nAsset), &spAsset);
			spAsset->put_NaturalSpot(date, afCurrentSpots[nAsset - 1]);
			ReinitialiseStrike(spAsset,bReinitialiseAll);
		}
	}
	
	return fPrice;
}

double CRisk::PriceVolatility(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult, std::string szUnderlyingFilter) const
//	fShift - additive shift such that the volatility (V) is transformed to V * (1.0 + fShift)
{
	double								fPrice;
	HRESULT								hr =S_OK;
	estring								sKey;


	sKey = "Volatility" + estring(fShift);
	if(szUnderlyingFilter.size()) {
		sKey +=":"+ szUnderlyingFilter;

	}

	if (!IsCached(sKey, spResult, &fPrice)){						
		CComPtr<IVolatilityStructures>	spVolatilityStructures;
		
		// Perturb the volatility structures
		Portfolio->GetVolatilityStructures(&spVolatilityStructures);

		if (szUnderlyingFilter.size()) {
			CComPtr<IVolatilityStructure> spVol; 
			long nVols = 0 ; 
			hr = spVolatilityStructures->get_Count(&nVols);

			for (long vIndx = 1  ; vIndx  <= nVols ; vIndx++ ) {
				CComBSTR sVol;
				if ( hr = spVolatilityStructures->get_Item(_variant_t(vIndx), &spVol)) propagate_error_ex(hr);
				spVol->get_Name(&sVol);
				if ( estring::CompareNoCase(szUnderlyingFilter, sVol) == 0 ) {
					spVol->Shift(fShift);
					break;
				}
			}
		}
		else {
			spVolatilityStructures->Shift(fShift);
		}
		
		// Perform the evaluation
		CComPtr<IResult>					sp;
		sp.CoCreateInstance(CLSID_Result);
		if (Portfolio->Evaluate(L"Price", &sp) || sp->GetPrice(&fPrice)) propagate_error;

		
		spResult->AddValue(sKey.GetBSTR(), CComVariant(fPrice));
		
		// Reset the volatility structures.
		spVolatilityStructures->Reset();		
	}

	return fPrice;
}


double CRisk::VegaSkew(IEvaluatable* Portfolio, CComPtr<IResult> spResult, DATE dtExpiry) const
{
	double								fVegaSkew;
	HRESULT								hr = S_OK; 

	if (!IsCached("VegaSkew" , spResult, &fVegaSkew)){		
		
		CComPtr<IAssets>				spAssets;
		CComPtr<IAsset>					spAsset;
		CComPtr<IStrikes>				spStrikes;
		CComPtr<IAssets>				spBaseAssets;
		DATE							date;
		long							dateOneYear;
		long							nAssets = 0 ;

		CParameterMap					pmStrikes;
		CParameterMap					pmDates;
		CParameterMap					pmShifts;

		double fPrice = Price(Portfolio, spResult);



		pmStrikes.SetSize(1,3);
		pmStrikes.SetValue(0,0, 0.9);
		pmStrikes.SetValue(0,1, 1.0);
		pmStrikes.SetValue(0,2, 1.1);


		pmShifts.SetSize(2,3);

		Portfolio->get_Date(&date);

		// Perturb the market data
		hr = Portfolio->GetUnderlyings(&spAssets);
		hr = spAssets->GetBaseUnderlyings(&spBaseAssets);

		CComPtr<IVolatilityStructures>	spVolatilityStructures;
		CComPtr<IVolatilityStructure> spVolatilityStructure;

		hr = spBaseAssets->get_Count(&nAssets);

		CComPtr<IDate> spValDate;
		CComPtr<IDate> spDateOneYear;
		CComPtr<IDate> spFarAwayDate;

		double fTimeToExpFactor ;
		long nFarDate;
		hr = spValDate.CoCreateInstance(L"Sirius.Date");
		if (hr = spValDate->put_SerialNumber (date) ) propagate_error_ex(hr);
		hr = spValDate->GetYearFraction( dtExpiry,  &fTimeToExpFactor );
		if (hr = spValDate->AddTenor(L"50Y" ,NoChange, vtMissing,VARIANT_TRUE,&spFarAwayDate)) propagate_error_ex(hr);
		if (hr = spValDate->AddTenor(L"1Y" ,NoChange, vtMissing,VARIANT_TRUE,&spDateOneYear)) propagate_error_ex(hr);

		hr =  spFarAwayDate->get_SerialNumber(&nFarDate);
		hr =  spDateOneYear->get_SerialNumber(&dateOneYear);

		pmDates.SetSize(2,1);
		pmDates.SetValue(0,0, date);
		pmDates.SetValue(1,0, nFarDate);

		fTimeToExpFactor = 1/pow(fTimeToExpFactor,0.5);


		CParameterMap pmStrikeData; 

		pmStrikeData.SetSize(1,3);
		pmStrikeData.SetValue(0,0,SpotBased);
		pmStrikeData.SetValue(0,1,pmStrikes.GetValue());


		for (long idx = 1 ; idx <= nAssets; ++idx ) {

			hr = spBaseAssets->get_Item (_variant_t(idx), &spAsset);			
			
			
			double fSpot = 0.0; 

			hr = spAsset->get_NaturalSpot(date , &fSpot );

			CComVariant vATMVol; 
			hr = MLGetNaturalATMVolatilityFromAsset(CComVariant(CComQIPtr<IDispatch>(spAsset)), CComVariant(date), CComVariant(dateOneYear),CComVariant(Spot), &vATMVol);
			hr = vATMVol.ChangeType(VT_R8);
			hr = spAsset->get_VolatilityStructure(&spVolatilityStructure);

			fTimeToExpFactor  *= vATMVol.dblVal  ; 
			fTimeToExpFactor  *= 0.01  ; 

			pmStrikeData.SetValue(0,2, fSpot);

			pmShifts.SetValue(-1,0,  fTimeToExpFactor );
			pmShifts.SetValue(-1,1,  0);
			pmShifts.SetValue(-1,2,  -fTimeToExpFactor );
			hr = spStrikes.CoCreateInstance(L"Sirius.Strikes");
			if (hr = spStrikes->put_Value(pmStrikeData.GetValue()))   propagate_error_ex(hr);


			CComVariant vDates  = pmDates.GetValue();
			CComVariant vShifts = pmShifts.GetValue();
			CComPtr<IArray> spDates;
			CComPtr<IMatrix> spShifts;

			//map_com_object_parameter(vDates, Array, spDates);
			hr = spDates.CoCreateInstance(L"Sirius.Array");
			hr = spDates->put_Value(vDates);

			//map_com_object_parameter(vShifts, Matrix, spShifts);
			hr = spShifts.CoCreateInstance(L"Sirius.Matrix");
			hr = spShifts->put_Value(vShifts);



			hr = spVolatilityStructure->ShiftSkew( spStrikes ,spDates, spShifts , spAsset);

		}


		// Perform the evaluation
		CComPtr<IResult>					sp;
		sp.CoCreateInstance(CLSID_Result);
		if (Portfolio->Evaluate(L"Price", &sp) || sp->GetPrice(&fVegaSkew)) propagate_error;


		fVegaSkew -=  fPrice; 

		spResult->AddValue(estring("VegaSkew").GetBSTR(), CComVariant(fVegaSkew));


		


		// Perturb the volatility structures
		hr = Portfolio->GetVolatilityStructures(&spVolatilityStructures);
		// Reset the volatility structures.
		spVolatilityStructures->Reset();		
	}

	return fVegaSkew ;
}

double CRisk::UnderlyingPrice(IEvaluatable* Portfolio, CComPtr<IResult> spResult) const
{
	double								fUnderlyingPrice = 0.0;
	HRESULT								hr = S_OK;
	CComPtr<IAssets>					spAssets;
	
	
	

	
	if (!IsCached("UnderlyingPrice", spResult, &fUnderlyingPrice)){	
		
		DATE		date = 0.0;
		hr = Portfolio->GetUnderlyings(&spAssets);
		hr = Portfolio->get_Date(&date);

		fUnderlyingPrice = GetUnderlyingPrice(spAssets, date);
	}

	spResult->AddValue(L"UnderlyingPrice", CComVariant(fUnderlyingPrice));

	return fUnderlyingPrice;

}

void ReinitialiseStrike( IAsset *Asset , bool bReinitialiseAll ) 
{

	HRESULT					hr = S_OK;
	CComPtr<IVolatilityStructure>	spVolatilityStructure;
	CComPtr<IVolatilityStructure>	spAsymVolatilityStructure;

	hr  = Asset->get_VolatilityStructure(&spVolatilityStructure);

	if (spVolatilityStructure) {
		if (bReinitialiseAll) {
			hr = spVolatilityStructure->ReinitialiseStrike(Asset);
		}

		hr  = spVolatilityStructure->get_AsymptoticVolatilityStructure(&spAsymVolatilityStructure);
	}


	if (SUCCEEDED(hr) && spAsymVolatilityStructure) {
		hr = spAsymVolatilityStructure->ReinitialiseStrike(Asset);
	}
	
}

bool CRisk::IsBasket(IAsset* pAsset) const 
{
	HRESULT					hr = S_OK;
	CComPtr<IAsset>			spAsset(pAsset);
	CComPtr<IAssets>		spAssets;
	long					nAssets;

	if ( hr = spAsset->get_Assets(&spAssets)) propagate_error ;

	if ( spAssets  == NULL ) return false ;

	if ( hr = spAssets->get_Count( &nAssets) ) propagate_error ;

	return  ( nAssets > 1  ) ;
}

double CRisk::GetUnderlyingPrice(IAssets *pAsset, DATE date ) const
{

	double								fUnderlyingPrice = 0.0;
	long								nAssets = 0L;
	HRESULT								hr = S_OK;
	CComPtr<IAssets>					spAssets(pAsset);
	double								fWeightTotal = 0.0;	


	hr = spAssets->get_Count(&nAssets);
	
	for (long idx = 1 ; idx  <= nAssets ; idx++) {
		CComPtr<IAsset>			spAsset;
		double					fSpot = 0.0;
		double					fWeight = 0.0;
		
		hr = spAssets->get_Item(CComVariant(idx), &spAsset);
		if (IsBasket(spAsset)) {

			CComPtr<IAssets>					spBskComp;
			hr = spAsset->get_Assets(&spBskComp);
			fSpot = GetUnderlyingPrice(spBskComp, date);
		}
		else 
		{
			hr = spAsset->get_NaturalSpot(date, &fSpot);
		}	

		hr = spAssets->get_Notional( _variant_t(CComQIPtr<IDispatch>(spAsset)), &fWeight);
		fUnderlyingPrice += fSpot * fWeight;
		fWeightTotal += fWeight;
	}

	return fUnderlyingPrice / fWeightTotal;

}




double CRisk::PriceTheta(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult) const
{

	double								fTheta =  0.0;		
	double								fPrice =0.0;		
	HRESULT								hr = S_OK;
	std::string							szError;

	if (!IsCached("Theta" + estring(fShift), spResult, &fPrice)) {


		CComPtr<IResult>					sp;
		sp.CoCreateInstance(CLSID_Result);

	
		

		//for now assume that the IEvaluatable* Portfolio is a product 
		CComQIPtr<IProduct> spProduct (Portfolio) ;
		if ( spProduct == NULL ) throw "currently can only calculate theta for a product";

		DATE	date;
		spProduct->get_Date(&date);

		spProduct->put_Date(date + fShift);

		//need to advance the date back
		AdvanceMarketDataDate(date , date +fShift );
		try
		{
			if (spProduct->Evaluate(L"Price", &sp) || sp->GetPrice(&fTheta)) propagate_error;
		} catch( LPCTSTR szErr ) {
				szError = szErr;
		} catch( std::string szErr ) {
				szError = szErr;
		} catch(...) {
			szError = "unknow error in calculating theta"; 
		}


		AdvanceMarketDataDate(date +fShift , date );
		spProduct->put_Date(date );


		if(szError.size()) throw szError;

		fTheta -= fPrice ; 

		spResult->AddValue(estring("Theta"+ estring(fShift)).GetBSTR(), CComVariant(fTheta));

		}

	return  fTheta;
}

HRESULT CRisk::AdvanceMarketDataDate(DATE dateFrom, DATE dateTo) const 
{
	HRESULT		hr =  S_OK ;

	AdvanceMarketDataDateProcessor processor (dateFrom,dateTo ) ;
	if (hr = ProcessMarketData("ZeroCurves", &processor ) ) return hr ;
	if (hr = ProcessMarketData("CorrelationMatrices", &processor ) ) return hr ;
	if (hr = ProcessMarketData("DividendSchedules", &processor ) ) return hr ;

	VolAdvanceMarketDataDateProcessor volProcessor(dateFrom,dateTo );
	if (hr = ProcessMarketData("VolatilityStructures", &volProcessor) ) return hr ;
	
	SpotAdvanceMarketDataDateProcessor spotProcessor(dateFrom,dateTo ); 
	if (hr = ProcessMarketData("SpotSchedules", &spotProcessor ) ) return hr ;

	AssetAdvanceMarketDataDateProcessor assetProcessor(dateFrom,dateTo ); 
	if (hr = ProcessMarketData("Assets", &assetProcessor ) ) return hr ;

	return hr;
}



HRESULT ProcessMarketData(std::string szCollectionName, MarketDataProcessor* pMarketData)
{

	HRESULT							hr =  S_OK ;
	CComQIPtr<IMarketData>			spMarketData;
	CComVariant						vCollection;
	CComVariant						vObject; 
	CComVariant						vDate=  0.0;
	CComVariant						vCount; 


	if (hr = _Module.GetSiriusApplication()->get_MarketData(&spMarketData) ) return hr;

	if (hr = CComDispatchDriverEx(spMarketData).GetPropertyByName(estring(szCollectionName).GetBSTR(), &vCollection)) return hr ;
	
	if (hr = CComDispatchDriverEx(vCollection.pdispVal).GetPropertyByName(L"Count", &vCount)) return hr;
	
	for (CComVariant vIdx = 1L ; vIdx.lVal <= vCount.lVal ; vIdx.lVal++ ) {

		if (hr = CComDispatchDriverEx(vCollection.pdispVal).GetPropertyByName(L"Item", &vIdx, 1,   &vObject)) return hr; 

		pMarketData->Apply( vObject.pdispVal);
	}

	return hr;
}





HRESULT VolAdvanceMarketDataDateProcessor::Apply( IDispatch* pObject)   
{
	HRESULT hr = S_OK;
	if (hr  =  AdvanceMarketDataDateProcessor::Apply(pObject) ) return hr ;

	CComQIPtr<IVolatilityStructure> spVol(pObject) ; 
	
	if (spVol == NULL) return E_POINTER;


	VolatilityDataTypeEnum vtEnum;
	if (spVol->get_VolatilityDataType(&vtEnum) ) return hr; 

	if ( vtEnum == TermVolatilities )  {

		CComPtr<IVolatilityStructure>		spAsymVol; 

		if( hr = spVol->get_AsymptoticVolatilityStructure(&spAsymVol) ) return hr ;
		if (spAsymVol == NULL ) return E_POINTER; 
		if( hr = spVol->put_Date(m_dateTo) ) return hr ;

	}

	return hr ; 
}


HRESULT SpotAdvanceMarketDataDateProcessor::Apply( IDispatch* pObject)   
{
	HRESULT hr = S_OK;
	if (hr  =  AdvanceMarketDataDateProcessor::Apply(pObject) ) return hr ;

	double fSpot;

	CComQIPtr<ISpotSchedule> spSpot(pObject) ; 
	
	if (spSpot == NULL) return E_POINTER;

	if (hr = spSpot->get_ValueAt(m_dateFrom , &fSpot) ) return hr ; 

	//ATLASSERT(m_dateFrom <  m_dateTo);

	for (long dt = m_dateFrom ; dt <= m_dateTo ; dt++   ) 
	{
		if (hr = spSpot->put_ValueAt(dt, fSpot) ) return hr ; 
	}

	return hr ; 
}


HRESULT AdvanceMarketDataDateProcessor::Apply( IDispatch* pObject)    
{
	HRESULT hr = S_OK;
	CComVariant				vDate; 
	CComPtr<IDispatch> spDispatch(pObject);

	if (hr = CComDispatchDriverEx(spDispatch).GetPropertyByName(L"Date",  &vDate)) return false; 

	vDate.ChangeType(VT_R8);
	if ( (		(long) (m_dateFrom)	    ==  (long)(vDate.dblVal) )
		  ||	((long)(m_dateTo)		==  (long)(vDate.dblVal) )
	   ) {

		CComVariant vToDate(m_dateTo);

		if (hr = CComDispatchDriverEx(spDispatch).PutPropertyByName(L"Date", &vToDate)) return hr; 

	}
	return hr ;
}

HRESULT AssetAdvanceMarketDataDateProcessor::Apply( IDispatch* pObject)    
{
	HRESULT							hr = S_OK;
	CComQIPtr<IAsset>				spAsset(pObject) ; 
	CComPtr<IAssets>				spAssets ; 
	
	CComQIPtr<IMarketData>			spMarketData;
	CComVariant						vCollection;

	if (hr  =  AdvanceMarketDataDateProcessor::Apply(pObject) ) return hr ;


	CComPtr<IZeroCurve> spCurrencyCurve ; 
	
	if (spAsset == NULL) return E_POINTER;

	if (hr = spAsset->get_CurrencyCurve( &spCurrencyCurve ))  return hr ;

	if (spCurrencyCurve ) {

		if (hr = spCurrencyCurve->put_Date(m_dateTo) ) return hr ;
	}


	ReinitialiseStrike(spAsset, true);

	//for baskets 
	hr = spAsset->get_RealAssets(&spAssets);

	long nComponents; 
	if (spAssets) {
		hr = spAssets->get_Count(&nComponents);
		for ( long idx = 1 ; idx <= nComponents; idx ++ ) {
			CComPtr<IAsset>			spAsset;

			hr = spAssets->get_Item(CComVariant(idx), &spAsset);
			if (spAsset) 
			{
				hr = Apply(spAsset);

			}

		}
	}

	if (hr = spAsset->Refresh( VARIANT_FALSE) ) return hr ;


	// if (hr = _Module.GetSiriusApplication()->get_MarketData(&spMarketData) ) return hr;
	// if (hr = CComDispatchDriverEx(spMarketData).GetPropertyByName(L"Assets", &vCollection)) return hr ;
	// CComQIPtr<IAssets>	spAssets(vCollection.pdispVal) ; 
	// if (hr =  spAssets->Refresh(VARIANT_FALSE) ) return hr ;

	return hr ;
}



