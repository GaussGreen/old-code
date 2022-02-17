//	handle_montecarlocollection.h : Declaration of the CMonteCarloCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MONTECARLOCOLLECTION_H_
#define __MONTECARLOCOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CMonteCarloCollection : public CSiriusComObjectCollection<IMonteCarlo, CMonteCarloCollection, IMonteCarloCollection, &CLSID_MonteCarloCollection, &IID_IMonteCarloCollection, &LIBID_Sirius, &CLSID_MonteCarlo, IDR_MONTECARLOCOLLECTION>
{
};

#endif
