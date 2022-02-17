/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/callschedule.h
// Purpose:     call schedule (common for CB) 
// Author:      Wang
// Created:     2004/08/17 
// RCS-ID:      $Id: callschedule.h,v 1.33 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/finance/bondlike/callschedule.h
   @brief declaration of the CallSchedule class for CB-like    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_CALLSCHEDULE_H_
#define _ITO33_FINANCE_BONDLIKE_CALLSCHEDULE_H_

#include "ito33/list.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/call.h"
#include "ito33/finance/bondlike/callperiod.h"
#include "ito33/finance/bondlike/triggeraspercentageof.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace XML
{
  class Tag;
}

namespace finance
{

/**
   CallSchedule class, it is composed of several periods of call.

   @iterator GetAll
 */
class ITO33_DLLDECL CallSchedule : public Call
{
public:

  /**
     Ctor calls the ctor of the base class. 
   */
  CallSchedule() : Call() {}

  /// type of the call data structure
  typedef std::list< shared_ptr<CallPeriod> > Elements;
   
  ///@name methods for initializing CallSchedule
  //@{

  /**
     Adds a call period to the call schedule.

     @param pCallPeriod A shared pointer to a call period. 
   */
  void AddCallPeriod(const shared_ptr<CallPeriod>& pCallPeriod);

  //@}  // name methods for initializing CallSchedule

  ///@name methods for accessing  CallSchedule
  //@{
  
  /**
     @brief Gets the list of the call periods.

     @return The list of the call periods.

     @noexport COM
   */
  const Elements& GetAll() const { return m_calls; }

  //@}  // name methods for accessing CallSchedule 

  /**
     @internal
     @brief Indicates if the callschedule contains a peiod with guaranteed yield.

     @return true if the callschedule contains a peiod with guaranteed yield, 
             false otherwise

     @noexport
   */
  bool HasYield() const;

  /**
     @internal
     @brief Indicates if the callschedule contains a soft call peiod.

     @return true if the callschedule contains a soft call peiod, 
             false otherwise

     @noexport
   */
  bool HasSoftCall() const;

  // implement virtual function  
  virtual bool HasMakeWhole() const;

  /// Validates the call schedule and throws exception if validation fails
  void Validate() const;
     
  /**
     @internal
     
     @noexport
   */
  XML::Tag Dump(XML::Tag& tagParent) const;


private:

  /// call data structure
  Elements m_calls;

}; // class CallSchedule


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_CALLSCHEDULE_H_
