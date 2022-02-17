//	handle_svlvoldatacollection.h : Declaration of the CSVLVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SVLVOLDATACOLLECTION_H_
#define __SVLVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CSVLVolDataCollection : public CSiriusComObjectCollection<ISVLVolData, CSVLVolDataCollection, ISVLVolDataCollection, &CLSID_SVLVolDataCollection, &IID_ISVLVolDataCollection, &LIBID_Sirius, &CLSID_SVLVolData, IDR_SVLVOLDATACOLLECTION>
{
};

#endif
