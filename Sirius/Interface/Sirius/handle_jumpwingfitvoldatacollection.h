//	handle_jumpwingfitvoldatacollection.h : Declaration of the CJumpWingFitVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __JUMPWINGFITVOLDATACOLLECTION_H_
#define __JUMPWINGFITVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CJumpWingFitVolDataCollection : public CSiriusComObjectCollection<IJumpWingFitVolData, CJumpWingFitVolDataCollection, IJumpWingFitVolDataCollection, &CLSID_JumpWingFitVolDataCollection, &IID_IJumpWingFitVolDataCollection, &LIBID_Sirius, &CLSID_JumpWingFitVolData, IDR_JUMPWINGFITVOLDATACOLLECTION>
{
};

#endif
