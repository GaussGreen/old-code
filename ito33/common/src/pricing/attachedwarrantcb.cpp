/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/attachedwarrantcb.cpp
// Author:      Ito33
// Created:     2005/01/13
// RCS-ID:      $Id: attachedwarrantcb.cpp,v 1.9 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"

#include "ito33/finance/bondlike/attachedwarrantconvertiblebond.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"

#include "ito33/finance/sessiondata.h"

#include "ito33/pricing/attachedwarrantcb.h"

namespace ito33
{

  ITO33_IMPLEMENT_AUTOPTR(pricing::AttachedWarrantConvertibleBond);

namespace pricing
{

AttachedWarrantConvertibleBond::AttachedWarrantConvertibleBond(
     const finance::AttachedWarrantConvertibleBond & warrant) : CB(),
     m_dCurrentConversionRatio(0.),
     m_conversionsShareDependent( warrant.GetShareDependentConversion(),
                                  warrant.GetSessionData()->GetValuationDate() )
{
  const shared_ptr<finance::ShareDependentConversion>&
    pShareDepConv( warrant.GetShareDependentConversion() );

  // note: must construct m_conversions before calling
  GetCBBaseData(warrant);

  if ( pShareDepConv->HasResetDate() )
  {
    // Check if the current conversion ratio has been set  
    m_dCurrentConversionRatio = pShareDepConv->GetCurrentRatio();

    // Note that if the current conversion ratio has been set but the
    // valuation date is before the begining of the reset period
    // the current conversion ratio does not affect the results. The current
    // conversion ratio set by the user if any should simply be ignored.
    if (    warrant.GetSessionData()->GetValuationDate()
         <= pShareDepConv->GetResetDate() )     
      m_dCurrentConversionRatio = 0.;
    else 
    {
      //update the trigger to still use the original base ratio
      m_conversionsShareDependent.ChangeTriggerRates(m_dCurrentConversionRatio);

      //current conversion is set
      m_conversionsShareDependent.SetRatios(m_dCurrentConversionRatio);
    }
  }
}


} //namespace pricing

} //namespace ito33
