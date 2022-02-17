//	handle_results.h : Declaration of the CResults
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __RESULTS_H_
#define __RESULTS_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CResults : public CSiriusComObjectCollection<IResult, CResults, IResults, &CLSID_Results, &IID_IResults, &LIBID_Sirius, &CLSID_Result, IDR_RESULTS>
{
protected:
	STDMETHOD(GetFields)(/*[out, retval]*/VARIANT* pVal);
	void								GetFieldsAddRecord(const std::string& szLowerCase, const std::string& szProperCase, std::map<std::string, VARTYPE>* pmapNameToVarType, std::map<std::string, std::string>* pmapLowerToProperCase, std::map<std::string, long>* pmapNameToSize, VARTYPE vt, long nSize);
};

#endif
