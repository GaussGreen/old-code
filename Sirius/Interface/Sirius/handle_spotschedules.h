//	handle_spotschedules.h : Declaration of the CSpotSchedules
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SPOTSCHEDULES_H_
#define __SPOTSCHEDULES_H_

#include "resource.h"
#include "comobjectcollectionserialisable.h"		// template collection class

class ATL_NO_VTABLE CSpotSchedules : public CComObjectCollectionSerialisable<ISpotSchedule, CSpotSchedules, ISpotSchedules, &CLSID_SpotSchedules, &IID_ISpotSchedules, &LIBID_Sirius, &CLSID_SpotSchedule, IDR_SPOTSCHEDULES>
{
public:
	static HRESULT						Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<ISpotSchedules>& spSpotSchedules);
protected:
	long								GetDate(void) const;	// override
};

#endif
