//	Aggregators.cpp : Implementation of CAggregators
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SiriusRisk.h"
#include "Aggregators.h"


STDMETHODIMP CAggregators::Aggregate(/*[in]*/ BSTR Name, /*[in]*/ IResults* pResults)
{
	begin_function
	HRESULT hr;	
	for (ContainerType::const_iterator it = m_coll.begin(); it != m_coll.end(); it++){
		CComQIPtr<IAggregatable> spAggregator = it->second.pdispVal;		
		if (spAggregator){
			if (hr = spAggregator->Aggregate(Name, pResults)) return hr;
		}
	}
	end_function
}