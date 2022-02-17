//	handle_parameterlists.h : Declaration of the CParameterLists
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PARAMETERLISTS_H_
#define __PARAMETERLISTS_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CParameterLists : public CSiriusComObjectCollection<IParameterList, CParameterLists, IParameterLists, &CLSID_ParameterLists, &IID_IParameterLists, &LIBID_Sirius, &CLSID_ParameterList, IDR_PARAMETERLISTS>
{
};

#endif
