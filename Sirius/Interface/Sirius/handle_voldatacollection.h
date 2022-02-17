//	handle_voldatacollection.h : Declaration of the CVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __VOLDATACOLLECTION_H_
#define __VOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CVolDataCollection : public CSiriusComObjectCollection<IVolData, CVolDataCollection, IVolDataCollection, &CLSID_VolDataCollection, &IID_IVolDataCollection, &LIBID_Sirius, &CLSID_VolData, IDR_VOLDATACOLLECTION>
{
};

#endif
