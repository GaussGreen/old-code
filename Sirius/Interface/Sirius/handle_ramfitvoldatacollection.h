//	handle_ramfitvoldatacollection.h : Declaration of the CRamFitVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __RAMFITVOLDATACOLLECTION_H_
#define __RAMFITVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CRamFitVolDataCollection : public CSiriusComObjectCollection<IRamFitVolData, CRamFitVolDataCollection, IRamFitVolDataCollection, &CLSID_RamFitVolDataCollection, &IID_IRamFitVolDataCollection, &LIBID_Sirius, &CLSID_RamFitVolData, IDR_RAMFITVOLDATACOLLECTION>
{
};

#endif
