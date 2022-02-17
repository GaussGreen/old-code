//	handle_hullvoldatas.h : Declaration of the CHullVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __HULLVOLDATACOLLECTION_H_
#define __HULLVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CHullVolDataCollection : public CSiriusComObjectCollection<IHullVolData, CHullVolDataCollection, IHullVolDataCollection, &CLSID_HullVolDataCollection, &IID_IHullVolDataCollection, &LIBID_Sirius, &CLSID_HullVolData, IDR_HULLVOLDATACOLLECTION>
{
};

#endif
