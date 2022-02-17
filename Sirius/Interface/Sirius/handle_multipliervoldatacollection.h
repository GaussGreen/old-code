//	handle_multipliervoldatacollection.h : Declaration of the CMultiplierVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MULTIPLIERVOLDATACOLLECTION_H_
#define __MULTIPLIERVOLDATACOLLECTION_H_

#include "resource.h"       
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CMultiplierVolDataCollection : public CSiriusComObjectCollection<IMultiplierVolData, CMultiplierVolDataCollection, IMultiplierVolDataCollection, &CLSID_MultiplierVolDataCollection, &IID_IMultiplierVolDataCollection, &LIBID_Sirius, &CLSID_MultiplierVolDataCollection,IDR_MULTIPLIERVOLDATACOLLECTION>
{
};

#endif 
