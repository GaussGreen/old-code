/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/basket_goodtype.h
// Purpose:     Calibration basket containing the good type of the instruments
// Created:     2006/06/26
// RCS-ID:      $Id: basket_goodtype.h,v 1.2 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/basket_goodtype.h
    @brief Calibration basket containing the good type of the instruments.

    This is the combination of derivative visitor and termstructure
    enumator.
 */

#ifndef _ITO33_FINANCE_BASKET_GOODTYPE_H_
#define _ITO33_FINANCE_BASKET_GOODTYPE_H_

#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivatives.h"

namespace ito33
{

namespace finance
{

/**
    Basket containing the right type of the object.

    This is a temporary solution since DerivativeVisitorGoodType copies the
    object instead of taking the pointer which should/can probably be avoided
    (maybe just keep a weak pointer instead of a copy).

    Currently, it relies on the validation previously done elsewhere(for 
    example, by parametrization info).

    Please assert before using the returned pointer.
 */
class BasketGoodType : public DerivativeVisitorGoodType, 
                       public TermStructureEnumerator
{
public:
  
  /**
      Gets the derivatives (if any) that can be used for minimization.

      @return The derivatives that can be used for minimization.
   */
  const shared_ptr<finance::Derivatives>& GetDerivatives() const
  {
    return m_pDerivatives;
  }

  /**
      Sets the derivatives (if any) that can be used for minimization.

      @param pDerivatives The derivatives that can be used for minimization.
   */
  void SetDerivatives(const shared_ptr<finance::Derivatives>& pDerivatives)
  {
    m_pDerivatives = pDerivatives;
  }

private:

  shared_ptr<finance::Derivatives> m_pDerivatives;

}; // class BasketGoodType


} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BASKET_GOODTYPE_H_
