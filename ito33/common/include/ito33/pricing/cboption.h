/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/cboption.h
// Purpose:     contract class for CBOption (backward)
// Author:      Nabil
// Created:     2005/10/13
// RCS-ID:      $Id: cboption.h,v 1.4 2006/04/04 16:29:46 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/cboption.h
    @brief The declaration of the CBOption contract class.

    The base class for CBOption contracts.  
 */

#ifndef _ITO33_PRICING_CBOPTION_H_
#define _ITO33_PRICING_CBOPTION_H_

#include "ito33/autoptr.h"
#include "ito33/vector.h"
#include "ito33/sharedptr.h"

#include "ito33/pricing/cb.h"
#include "ito33/pricing/cboptiondata.h"
#include "ito33/pricing/contract.h"

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL CBOption;
  class ITO33_DLLDECL YieldCurve;
}

namespace pricing
{

/// The declaration of the (backward) CBOption contract class.
class CBOption : public Contract 
{

public:

  /**
     Creates a CB by financial CBOption object.

     @param cboption financial CBOption
   */
  CBOption(const finance::CBOption& cboption);

  /// virtual dtor for base class
  virtual ~CBOption() {}

  CB& GetCB() { return m_cb; }
  
  /**
     Gets the cb option data.

     @return a pointer to the cb option data.
   */
  CBOptionData* GetCBOptionData() const { return m_pCBOptionData.get(); }

protected:

  /// cb of the cb option
  CB m_cb;

  /// cb option data (for the cb option contract)
  AutoPtr<CBOptionData> m_pCBOptionData;

}; // class CBOption;


} // namespace pricing

} // namespace ito33

#endif // #ifndef _ITO33_PRICING_CBOPTION_H_

