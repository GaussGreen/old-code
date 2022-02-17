//	risk.h: interface for the CRisk class.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _RISK_H
#define _RISK_H


#pragma once

class CRisk
{
public:
	
	HRESULT								AdvanceMarketDataDate( DATE nFromDate , DATE  nToDate ) const ;
	double								GetUnderlyingPrice (IAssets* pAsset, DATE date) const;
	bool								IsBasket( IAsset* pAsset) const ;
	
	CComPtr<IResult>					GetResultObject(IResult*** pppResult) const;
	double								Price(IEvaluatable* Portfolio, CComPtr<IResult> spResult) const;
	double								PriceSpot(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult, bool bReinitialiseAll, std::string szUnderlyingFilter) const;
	double								PriceVolatility(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult,std::string szUnderlyingFilter) const;
	double								PriceIR(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult) const;
	double								PriceDiv(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult) const;


	double								UnderlyingPrice(IEvaluatable* Portfolio, CComPtr<IResult> spResult) const; 
	double								PriceTheta(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult) const;

	double								Delta(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult, bool bReinitialiseAll, std::string szUnderlingFilter  = "") const;
	double								Rho(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult ) const;
	double								DivDelta(double fShift, IEvaluatable* Portfolio, CComPtr<IResult> spResult ) const;
	double								VegaSkew(IEvaluatable* Portfolio, CComPtr<IResult> spResult, DATE dtExpiry) const;
	

	
protected:
	bool								IsCached(const estring& szTag, CComPtr<IResult> spResult, double* pfOut) const;
};

void									ReinitialiseStrike( IAsset* Asset ,  bool bReinitialiseAll ) ;

/*
template <class ICollection , class  ISingular > 
HRESULT AdvanceMarketDataDate(std::string szCollectionName, DATE  nFromDate, DATE nToDate)
{

	HRESULT						hr =  S_OK ;
	CComPtr<IMarketData>		spMarketData;
	CComVariant					vCollection;
	CComPtr<ICollection>		spCollection;
	long						nCount;
	CComPtr<ISingular>			spSingular; 
	DATE						dtDate=  0.0;


	if (hr = _Module.GetSiriusApplication()->get_MarketData(&spMarketData) ) return hr;
	
	if (hr = CComDispatchDriverEx(spMarketData).GetPropertyByName(estring(szCollectionName).GetBSTR(), &vCollection)) return hr ;
	
	spCollection = vCollection.pdispVal;

	if (hr = spMarketData->get_ZeroCurves(&spZeroCurves) ) return hr;

	if (hr = spCollection->get_Count(& nCount ) ) return hr; 

	for (long idx = 0 ; idx < nCount ; idx++ ) {
		if ( hr = spCollection->get_Item(_variant_t(idx) , &spSingular ) ) return hr ; 
		if ( hr = spSingular->get_Date(&dtDate ) ) return hr; 

		if ( (		(long)(nFromDate)	==  (long)(dtDate) )
			  ||	(long)(nToDate)		==  (long)(dtDate) )
		   ) {

			if (hr = spSingular->put_Date(nFromDate) ) return hr;
		}
	}

	return hr;
}

*/

/*
HRESULT CRisk::AdvanceZeroCurveDate(DATE  nFromDate, DATE nToDate)
{

	HRESULT						hr =  S_OK ;
	CComPtr<IMarketData>		spMarketData;
	CComPtr<IZeroCurves>		spZeroCurves;
	long						nCount;
	CComPtr<IZeroCurve>			spZeroCurve; 
	DATE						dtDate=  0.0;


	if ( hr = _Module.GetSiriusApplication()->get_MarketData(&spMarketData) ) return hr;
	if ( hr = spMarketData->get_ZeroCurves(&spZeroCurves) ) return hr;
	if ( hr = spZeroCurves->get_Count(& nCount ) ) return hr; 

	for (long idx = 0 ; idx < nCount ; idx++ ) {
		if ( hr = spZeroCurves->get_Item(_variant_t(idx) , &spZeroCurve ) ) return hr ; 
		if ( hr = spZeroCurve->get_Date(&dtDate ) ) return hr; 

		if ( (		static_cast<long>(nFromDate)	==  static_cast<long>(dtDate) )
			  ||	(static_cast<long>(nToDate)		==  static_cast<long>(dtDate) )
		   ) {

			if (hr = spZeroCurve->put_Date(nFromDate) ) return hr;
		}
	}

	return hr;
}

*/

class MarketDataProcessor 
{
public:
	virtual HRESULT Apply( IDispatch* pObject) = 0 ;
};


class AdvanceMarketDataDateProcessor : public MarketDataProcessor 
{
protected:
	DATE  m_dateFrom;
	DATE  m_dateTo;
public:

	AdvanceMarketDataDateProcessor(DATE  dateFrom, DATE dateTo)
		:m_dateFrom(dateFrom)
		,m_dateTo(dateTo)
	{
	}
	
	virtual HRESULT Apply( IDispatch* pObject)  ;
};


class VolAdvanceMarketDataDateProcessor : public AdvanceMarketDataDateProcessor  
{
public:

	VolAdvanceMarketDataDateProcessor  (DATE  dateFrom, DATE dateTo)
		: AdvanceMarketDataDateProcessor(dateFrom, dateTo)
	{

	}

	virtual HRESULT Apply( IDispatch* pObject)   ;

} ;


class SpotAdvanceMarketDataDateProcessor : public AdvanceMarketDataDateProcessor  
{
public:

	SpotAdvanceMarketDataDateProcessor  (DATE  dateFrom, DATE dateTo)
		: AdvanceMarketDataDateProcessor(dateFrom, dateTo)
	{

	}

	virtual HRESULT Apply( IDispatch* pObject)   ;
} ;


class AssetAdvanceMarketDataDateProcessor : public AdvanceMarketDataDateProcessor  
{
public:

	AssetAdvanceMarketDataDateProcessor  (DATE  dateFrom, DATE dateTo)
		: AdvanceMarketDataDateProcessor(dateFrom, dateTo)
	{

	}

	virtual HRESULT Apply( IDispatch* pObject)   ;
} ;



HRESULT ProcessMarketData(std::string szCollectionName, MarketDataProcessor* pMarketData);




#endif
