/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/onetouch.h
// Purpose:     Contracts class for OneTouch
// Created:     2005/07/04
// RCS-ID:      $Id: onetouch.h,v 1.1 2005/07/05 10:03:23 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/pricing/onetouch.h
   @brief The declaration of the OneTouch contract class.
 */

#ifndef _ITO33_PRICING_ONETOUCH_H_
#define _ITO33_PRICING_ONETOUCH_H_

#include "ito33/finance/exoticoption/onetouch.h"

#include "ito33/pricing/contract.h"

namespace ito33
{

namespace pricing
{


/// The declaration of the OneTouch contract class.
class OneTouch : public Contract 
{
public:
  
  OneTouch(const finance::OneTouch& oneTouch);

  // Default dtor is ok 
  double GetBarrier() const { return m_dBarrier; }

  finance::BarrierType GetBarrierType() const { return m_barrierType; }

  bool IsImmediateRebate() const { return m_bImmmediate; }


private:

  double m_dBarrier;

  finance::BarrierType m_barrierType;

  bool m_bImmmediate;

  NO_COPY_CLASS(OneTouch);

}; // class OneTouch;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_ONETOUCH_H_
