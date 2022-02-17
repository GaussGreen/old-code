//	handle_matrices.h : Declaration of the CInterfaceMatrices
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MATRICES_H_
#define __MATRICES_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CInterfaceMatrices: public CSiriusComObjectCollection<IMatrix, CInterfaceMatrices, IMatrices, &CLSID_Matrices, &IID_IMatrices, &LIBID_Sirius, &CLSID_Matrix, IDR_MATRICES>
{
};

#endif
