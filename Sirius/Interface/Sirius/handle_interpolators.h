//	handle_interpolators.h : Declaration of the CInterpolators
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __INTERPOLATORS_H_
#define __INTERPOLATORS_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CInterpolators : public CSiriusComObjectCollection<IInterpolator, CInterpolators, IInterpolators, &CLSID_Interpolators, &IID_IInterpolators, &LIBID_Sirius, &CLSID_Interpolator, IDR_INTERPOLATORS>
{
};

#endif
