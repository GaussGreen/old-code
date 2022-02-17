//	handle_pdecollection.h : Declaration of the CPdeCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PDECOLLECTION_H_
#define __PDECOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CPdeCollection : public CSiriusComObjectCollection<IPde, CPdeCollection, IPdeCollection, &CLSID_PdeCollection, &IID_IPdeCollection, &LIBID_Sirius, &CLSID_Pde, IDR_PDECOLLECTION>
{
};

#endif
