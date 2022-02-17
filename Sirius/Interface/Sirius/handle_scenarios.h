//	handle_scenarios.h : Declaration of the CScenarios
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SCENARIOS_H_
#define __SCENARIOS_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CScenarios : public CSiriusComObjectCollection<IScenario, CScenarios, IScenarios, &CLSID_Scenarios, &IID_IScenarios, &LIBID_Sirius, &CLSID_Scenario, IDR_SCENARIOS>
{
};

#endif
