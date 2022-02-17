//	handle_arrays.h : Declaration of the CArrays
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __ARRAYS_H_
#define __ARRAYS_H_

#include "resource.h"					// main symbols
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CArrays : public CSiriusComObjectCollection<IArray, CArrays, IArrays, &CLSID_Arrays, &IID_IArrays, &LIBID_Sirius, &CLSID_Array, IDR_ARRAYS>
{	
};

#endif
