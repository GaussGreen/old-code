//	ResultsCollection.h : Declaration of the CResultsCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __RESULTSCOLLECTION_H_
#define __RESULTSCOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CResultsCollection : public CSiriusComObjectCollection<IResults, CResultsCollection, IResultsCollection, &CLSID_ResultsCollection, &IID_IResultsCollection, &LIBID_Sirius, &CLSID_Results, IDR_RESULTSCOLLECTION>
{
};

#endif
