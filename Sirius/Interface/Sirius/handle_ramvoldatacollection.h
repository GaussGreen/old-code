//	handle_ramvoldatacollection.h : Declaration of the CRamVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __RAMVOLDATACOLLECTION_H_
#define __RAMVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CRamVolDataCollection : public CSiriusComObjectCollection<IRamVolData, CRamVolDataCollection, IRamVolDataCollection, &CLSID_RamVolDataCollection, &IID_IRamVolDataCollection, &LIBID_Sirius, &CLSID_RamVolData, IDR_RAMVOLDATACOLLECTION>
{
};

#endif
