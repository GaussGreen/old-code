//	handle_dateschedules.h : Declaration of the CDateSchedules
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DATESCHEDULES_H_
#define __DATESCHEDULES_H_

#include "resource.h"					// main symbols
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CDateSchedules : public CSiriusComObjectCollection<IDateSchedule, CDateSchedules, IDateSchedules, &CLSID_DateSchedules, &IID_IDateSchedules, &LIBID_Sirius, &CLSID_DateSchedule, IDR_DATESCHEDULES>
{
};

#endif
