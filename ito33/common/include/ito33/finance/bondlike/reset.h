/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/reset.h
// Purpose:     resettable convertible bond class
// Author:      ITO 33
// Created:     2004/10/05
// RCS-ID:      $Id: reset.h,v 1.11 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/reset.h
    @brief declaration of the financial resettable convertible bond class
    
 */

#ifndef _ITO33_FINANCE_BONDLIKE_RESET_H_
#define _ITO33_FINANCE_BONDLIKE_RESET_H_

#include "ito33/finance/bondlike/cb_base.h"

#include "ito33/dlldecl.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL ResetConversionSchedule;

/**
   A resettable convertible bond.
 */
class ITO33_DLLDECL Reset : public CBBase
{
public:
 
  /**
     Creates a reset bond using BondTerms and ResetConversionSchedule.

     @param pBondTerms the principle characteristics of the bond
     @param pResetConversionSchedule reset parameters
   */
  Reset(const shared_ptr<BondTerms>& pBondTerms,
        const shared_ptr<ResetConversionSchedule>& pResetConversionSchedule);

  /// virtual dtor for base class
  virtual ~Reset() { }

  /**
      @name Methods for accessing convertible bond specific date.
   */
  //@{

  /**
     Gets the reset schedule.

     @return reset schedule
   */
  const shared_ptr<ResetConversionSchedule>& GetResetConversionSchedule() const
  {
    return m_pResetConversionSchedule;
  }

  //@}

  void SetCallSchedule(const shared_ptr<CallSchedule>& pCallSchedule);

  virtual void Validate() const;
  virtual void ValidateWith(const SessionData& sessionData) const;

  virtual void Visit(DerivativeVisitor& visitor) const;
  virtual void Visit(DerivativeModifyingVisitor& visitor);
  virtual XML::Tag Dump(XML::Tag& tagParent) const;

protected:

  /**
      Checks that the call has no trigger period since we don't support it.

      Throws exception if the call has trigger period.
   */
  void CheckTriggerPeriod() const;

  /// reset schedule (dates, caps, floors, etc)
  shared_ptr<ResetConversionSchedule> m_pResetConversionSchedule;

}; // class Reset


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_RESET_H_
