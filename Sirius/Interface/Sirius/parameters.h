//	parameters.h : Declaration of the CParameters
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PARAMETERS_H_
#define __PARAMETERS_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"
template<class T> class MlEqDictionary;

class ATL_NO_VTABLE CParameters : public CSiriusComObjectCollection<IParameter, CParameters, IParameters, &CLSID_Parameters, &IID_IParameters, &LIBID_Sirius, &CLSID_Parameter, IDR_PARAMETERS>
{
public:
	bool IsInCollection(const CComBSTR& s);
	void GetDictionary(MlEqDictionary<CComVariant>* pdict);
};

#endif
