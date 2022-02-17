// product_funding.cpp : Implementation of CFunding
#include "stdafx.h"
#include "TemporaryProducts.h"
#include "product_funding.h"
#include "MlEqZeroCurve.h"

/////////////////////////////////////////////////////////////////////////////
// CFunding



STDMETHODIMP CFunding::Evaluate(IResult *pVal)
{
	begin_function

	HRESULT					hr;	
								
	if (!m_hFixingSchedule)	throw "No fixing dates defined" ;
	if (!m_hCoupons)		throw "No coupons defined";	
	if (!m_hZeroCurve)		throw "No zero curve defined";
					
	begin_calculate(pVal, L"Price")

	double discountedCouponSum = 0.0;

	long nToday = m_hZeroCurve->GetReferenceDate();
	GVector<long>	afDates;
	m_hFixingSchedule->GetDates(afDates);

	int sz = afDates.getsize();
	if(m_hCoupons->getsize() != sz)	throw "Date and coupon schedule sizes mismatch";

	for(int i=0; i<sz;i++)
	{
		long date = afDates[i];
		if(date >= nToday) // >=
		{
			double df = m_hZeroCurve->GetDiscountFactor(nToday, date);
			discountedCouponSum += 	df * (*m_hCoupons)[i];
		}
	}
	
	
	if (hr = pVal->AddValue(L"Price", CComVariant(discountedCouponSum))) return hr;
	end_calculate


	end_function
}

HRESULT CFunding::FinalConstruct()
{
	return S_OK;
}



STDMETHODIMP CFunding::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFunding
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}
