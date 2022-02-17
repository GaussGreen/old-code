//	handle_strikescollection.h : Declaration of the CStrikesCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __STRIKESCOLLECTION_H_
#define __STRIKESCOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CStrikesCollection : public CSiriusComObjectCollection<IStrikes, CStrikesCollection, IStrikesCollection, &CLSID_StrikesCollection, &IID_IStrikesCollection, &LIBID_Sirius, &CLSID_Strikes, IDR_STRIKESCOLLECTION>
{
};

#endif
