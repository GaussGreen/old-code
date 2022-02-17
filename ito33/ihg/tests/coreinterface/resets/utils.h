

//
// The construction of the cb contracts is
// done here to make the rest of the code
// clearer
//
#ifndef _UTILS_
#define _UTILS_

#include "ito33/date.h"
#include "ito33/sharedptr.h"

namespace ito33
{
  namespace finance
  {
    class SessionData;
    class ConvertibleBond;
    class Reset;
    class ResetConversionSchedule;
    class ConversionSchedule;
  }

shared_ptr<finance::ConvertibleBond> InitCB(
  const shared_ptr<finance::SessionData>& pSessionData,
  Date issueDate, Date maturityDate,
  const shared_ptr<finance::ConversionSchedule>& 
  pconversionSchedule, double dParValue);

shared_ptr<finance::ConvertibleBond> InitCB(
  const shared_ptr<finance::SessionData>& pSessionData,
  Date issueDate, Date maturityDate, double dParValue);

shared_ptr<finance::Reset> InitReset(
  const shared_ptr<finance::SessionData>& pSessionData,
  Date issueDate,
  Date maturityDate,
  const shared_ptr<finance::ResetConversionSchedule>& 
        pResetConversionSchedule, double dParValue);

shared_ptr<finance::SessionData> 
   InitSessionData(Date issueDate, bool bFakeDividend, Date valuationDate = Date());
}

#endif
