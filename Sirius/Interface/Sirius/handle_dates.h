//	handle_dates.h : Declaration of the CDates
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DATES_H_
#define __DATES_H_

#include "resource.h"					// main symbols
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CDates : public CSiriusComObjectCollection<IDate, CDates, IDates, &CLSID_Dates, &IID_IDates, &LIBID_Sirius, &CLSID_Date, IDR_DATES>
{
};

#endif
