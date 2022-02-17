//	handle_svlfitvoldatacollection.h : Declaration of the CSVLFitVolDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SVLFITVOLDATACOLLECTION_H_
#define __SVLFITVOLDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CSVLFitVolDataCollection : public CSiriusComObjectCollection<ISVLFitVolData, CSVLFitVolDataCollection, ISVLFitVolDataCollection, &CLSID_SVLFitVolDataCollection, &IID_ISVLFitVolDataCollection, &LIBID_Sirius, &CLSID_SVLFitVolData, IDR_SVLFITVOLDATACOLLECTION>
{
};

#endif
