//	handle_hermitevoldatacollection.h : Declaration of the CHermiteVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __HERMITEVOLDATACOLLECTION_H_
#define __HERMITEVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CHermiteVolDataCollection : public CSiriusComObjectCollection<IHermiteVolData, CHermiteVolDataCollection, IHermiteVolDataCollection, &CLSID_HermiteVolDataCollection, &IID_IHermiteVolDataCollection, &LIBID_Sirius, &CLSID_HermiteVolData, IDR_HERMITEVOLDATACOLLECTION>
{
};

#endif
