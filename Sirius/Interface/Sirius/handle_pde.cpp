//	handle_pde.cpp : Implementation of CPde
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_pde.h"
#include "MlEqDate.h"
#include "MlEqShortModels.h"
//#include "MlEqPde.h"


STDMETHODIMP CPde::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IPde };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CPde::get_Value(VARIANT *pVal)
{		
	begin_function
	throw "CPde::get_Value is not yet implemented.";
	end_function
}

STDMETHODIMP CPde::put_Value(VARIANT newVal)
{
	HRESULT								hr;	
	std::vector<CParameterMap>			vpm;
				
	begin_function		
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	
	map_parameter(vpm[0], long, PdeType);
	switch (PdeType){
	case BlackScholes:
		if (vpm.size() !=5){
			throw "Invalid number of parameters";
		} else {					
			RCPtr<pdeBlackScholes>			h;
			
			map_parameter(vpm[1], double, Vol);
			map_parameter(vpm[2], double, LoanSpread);
			map_parameter(vpm[3], double, DiscountRate);
			map_parameter(vpm[4], double, Spot);
		
			h = new pdeBlackScholes;
			h->m_vol = Vol;
			h->m_loanspread = LoanSpread;
			h->m_discount_rate = DiscountRate;
			h->m_spot = Spot;
		
			// ToDo - if all is well then assign to m_h
			m_h = &*h;
		}
		break;
	case LocalVol:
		if (vpm.size() != 6){
			throw "Invalid number of parameters";
		} else {

			RCPtr<pdeLocalVol> h = new pdeLocalVol;

			DupireLocalVolHandle lv;

			map_object_parameter(vpm[1], Asset, hAsset);
			map_parameter(vpm[1], double, lowSpot);
			map_parameter(vpm[2], double, highSpot);
			map_parameter(vpm[3], long, maturityDate);
			map_parameter(vpm[4], long, nt);
			map_parameter(vpm[5], long, nx);

			lv->initialize(hAsset,lowSpot,highSpot,maturityDate,nt,nx);
			//h->init(lv);

			m_h = &*h;

			throw "I'm trying to create a Horsey type pde here";
		}
		break;
	default:
		throw "Invalid or unsupported pde type";
	}
		
	end_function
}
