/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/spotfxrates.cpp
// Purpose:     implementation of SpotFXRates
// Author:      Wang
// Created:     2004/09/02
// RCS-ID:      $Id: spotfxrates.cpp,v 1.12 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/numeraire.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/spotfxrates.h"
#include "ito33/xml/finance/common.h"

extern const ito33::finance::Error ITO33_SPOTFXRATES_NEGATIVERATE,
                          ITO33_SPOTFXRATES_MISSING_RATE,
                          ITO33_NUMERAIRE_INVALID;

namespace ito33
{
 
namespace finance
{

void SpotFXRates::SetFXRate(const shared_ptr<Numeraire>& pNumeraire1, 
                            const shared_ptr<Numeraire>& pNumeraire2, 
                            double dRate)
{

  // standard checks
  CHECK_COND( dRate > 0.0, ITO33_SPOTFXRATES_NEGATIVERATE );
  CHECK_COND( pNumeraire1, ITO33_NUMERAIRE_INVALID );
  CHECK_COND( pNumeraire2, ITO33_NUMERAIRE_INVALID );

  // Compare the objects
  if ( *pNumeraire1 == *pNumeraire2)
    return;

  Elements::iterator pRate;
  
  // update if the rate (numeraire2, numeraire1) is already set
  pRate = m_FXRates.find
          ( NumerairePair(pNumeraire2->GetCode(), pNumeraire1->GetCode()) );

  if ( pRate != m_FXRates.end() )
  {
    m_FXRates[pRate->first] = 1. / dRate;

    return;
  }   
 
  m_FXRates[NumerairePair(pNumeraire1->GetCode(), pNumeraire2->GetCode())] 
    = dRate;
}


double SpotFXRates::GetFXRate(const shared_ptr<Numeraire>& pNumeraire1, 
                              const shared_ptr<Numeraire>& pNumeraire2)
{
  // standard checks
  CHECK_COND( pNumeraire1, ITO33_NUMERAIRE_INVALID );
  CHECK_COND( pNumeraire2, ITO33_NUMERAIRE_INVALID );

  // Compare objects
  if ( *pNumeraire1 == *pNumeraire2 )
    return 1.;

  // need for check for (num1, num2) and (num2, num1)
  Elements::iterator pRate;

  pRate = m_FXRates.find
          ( NumerairePair(pNumeraire1->GetCode(), pNumeraire2->GetCode()) );

  if ( pRate != m_FXRates.end() )
    return pRate->second;

  pRate = m_FXRates.find
          ( NumerairePair(pNumeraire2->GetCode(), pNumeraire1->GetCode()) );

  if ( pRate != m_FXRates.end() )
    return 1. / pRate->second;

  // rate not found.
  throw EXCEPTION( ITO33_SPOTFXRATES_MISSING_RATE );

  // Just to fix compiler warning. Should never get here.
  return 0.;
}


XML::Tag SpotFXRates::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagFX(XML_TAG_SPOT_FX_RATES_ROOT, tagParent);

  Elements::const_iterator iter;

  for(iter = m_FXRates.begin(); iter != m_FXRates.end(); iter++)
  {
    XML::Tag tagElement(XML_TAG_SPOT_FX_ROOT, tagFX);

    tagElement.Element
                (XML_TAG_SPOT_FX_FOREIGN_CURRENCY)
                (iter->first.first);
    tagElement.Element
                (XML_TAG_SPOT_FX_BASE_CURRENCY)
                (iter->first.second);
 
    tagElement.Element(XML_TAG_FINANCE_RATE)(iter->second);
  }

  return tagFX;
}


} // namespace finance

} // namespace ito33
