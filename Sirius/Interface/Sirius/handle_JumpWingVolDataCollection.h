//	handle_jumpwingvoldatacollection.h : Declaration of the CJumpWingVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __JUMPWINGVOLDATACOLLECTION_H_
#define __JUMPWINGVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CJumpWingVolDataCollection : public CSiriusComObjectCollection<IJumpWingVolData, CJumpWingVolDataCollection, IJumpWingVolDataCollection, &CLSID_JumpWingVolDataCollection, &IID_IJumpWingVolDataCollection, &LIBID_Sirius, &CLSID_JumpWingVolData, IDR_JUMPWINGVOLDATACOLLECTION>
{
};

#endif
