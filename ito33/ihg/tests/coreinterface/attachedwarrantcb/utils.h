// Some helper functions to construct contracts
#ifndef _UTILS_
#define _UTILS_

#include "ito33/sharedptr.h"

namespace ito33
{
  class Date;

  namespace finance
  {
    class SessionData;
    class ConvertibleBond;
    class BondTerms;
    class AttachedWarrantConvertibleBond;
    class CallSchedule;
  }

shared_ptr<finance::AttachedWarrantConvertibleBond> 
InitAttachedWarrantCB(const shared_ptr<finance::SessionData>& pSessionData,
                      double dShareFactor = 1.0,
                      double dStrike = -1.0,
                      double dCap = 2.0);

shared_ptr<finance::AttachedWarrantConvertibleBond> 
InitAttachedWarrantCBNoReset(double dShareFactor, double dCap);

shared_ptr<finance::AttachedWarrantConvertibleBond> 
  InitCarnivalTest(bool bShiftCallDate = false);

shared_ptr<finance::AttachedWarrantConvertibleBond> 
   InitAmericanExpressTest(bool bShiftCallDate = false);

shared_ptr<finance::AttachedWarrantConvertibleBond> 
  InitAutobacsTest();

shared_ptr<finance::ConvertibleBond> 
InitCB(const shared_ptr<finance::SessionData>& pSessionData);

shared_ptr<finance::BondTerms> InitBondTerms();

shared_ptr<finance::SessionData> InitSessionData();

shared_ptr<finance::CallSchedule> GetCallSchedule(Date callStart,Date callEnd);

shared_ptr<finance::AttachedWarrantConvertibleBond> InitGettyTest();
}

#endif
