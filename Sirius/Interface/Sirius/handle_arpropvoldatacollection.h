//	handle_arpropvoldatacollection.h : Declaration of the CARPropVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __ARPROPVOLDATACOLLECTION_H_
#define __ARPROPVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CARPropVolDataCollection : public CSiriusComObjectCollection<IARPropVolData, CARPropVolDataCollection, IARPropVolDataCollection, &CLSID_ARPropVolDataCollection, &IID_IARPropVolDataCollection, &LIBID_Sirius, &CLSID_ARPropVolData, IDR_ARPROPVOLDATACOLLECTION>
{
};

#endif
