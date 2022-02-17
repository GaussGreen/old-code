//  handle_hullandwhitecollection.h : Declaration of the CHullAndWhiteCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __HULLANDWHITECOLLECTION_H_
#define __HULLANDWHITECOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CHullAndWhiteCollection : public CSiriusComObjectCollection<IHullAndWhite, CHullAndWhiteCollection, IHullAndWhiteCollection, &CLSID_HullAndWhiteCollection, &IID_IHullAndWhiteCollection, &LIBID_Sirius, &CLSID_HullAndWhite, IDR_HULLANDWHITECOLLECTION>
{
};

#endif
