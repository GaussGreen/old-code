

//
// The construction of the cb contracts is
// done here to make the rest of the code
// clearer
//
#ifndef _UTILS_
#define _UTILS_

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/convertiblebond.h"

ito33::shared_ptr<ito33::finance::ConvertibleBond> InitCB(
  const ito33::shared_ptr<ito33::finance::SessionData>& pSessionData,
  const ito33::Date issueDate,const ito33::Date maturityDate,
  const ito33::shared_ptr<ito33::finance::ConversionSchedule>& pconversionSchedule);

ito33::shared_ptr<ito33::finance::ConvertibleBond> InitCB(
  const ito33::shared_ptr<ito33::finance::SessionData>& pSessionData,
  const ito33::Date issueDate,const ito33::Date maturityDate);

ito33::shared_ptr<ito33::finance::SessionData> 
InitSessionData(const ito33::Date issueDate);

#endif
