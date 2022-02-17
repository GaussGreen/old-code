//	handle_arpropfitvoldatacollection.h : Declaration of the CARPropFitVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __ARPROPFITVOLDATACOLLECTION_H_
#define __ARPROPFITVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CARPropFitVolDataCollection : public CSiriusComObjectCollection<IARPropFitVolData, CARPropFitVolDataCollection, IARPropFitVolDataCollection, &CLSID_ARPropFitVolDataCollection, &IID_IARPropFitVolDataCollection, &LIBID_Sirius, &CLSID_ARPropFitVolData, IDR_ARPROPFITVOLDATACOLLECTION>
{
};


#endif
