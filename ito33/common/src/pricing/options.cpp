/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/options.cpp
// Created:     2004/02/10
// RCS-ID:      $Id: options.cpp,v 1.16 2006/08/09 13:45:19 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/list.h"
#include "ito33/dateutils.h"

#include "ito33/finance/forwardoption.h"

#include "ito33/pricing/options.h"

using namespace ito33::finance;
using namespace ito33::pricing;
using namespace ito33::numeric::mesh;

Options::Options(const ForwardOption& forwardOption)
{
  ASSERT( !forwardOption.IsEmpty() );

  const ForwardOption::Elements& options( forwardOption.GetAll() );

  ForwardOption::Elements::const_iterator iter;

  for (iter = options.begin(); iter != options.end(); ++iter)
  {
    const finance::Option& option = *iter->m_pOption;

    m_pMaturityDates.push_back( option.GetMaturityDate() );
    m_pdStrikes.push_back( option.GetStrike() );
    m_pOptionTypes.push_back( option.GetOptionType() );
    m_pdWeights.push_back( iter->m_dWeight );

    if ( iter->m_pOption->HasMarketPrice() )
      m_pdMarketPrices.push_back( iter->m_pOption->GetMarketPrice() );
    else
       m_pdMarketPrices.push_back(0);
  }
}

double Options::GetMaturityTime() const 
{ 
  return GetDoubleFrom( m_pMaturityDates.back() ); 
}

void Options::GetSpecialTimes(SpecialTimes& specialTimes) const
{
  specialTimes.clear();

  std::vector<Date>::const_iterator pMaturity;

  for (pMaturity = m_pMaturityDates.begin(); 
       pMaturity != m_pMaturityDates.end();
       ++pMaturity)
  {
    specialTimes.push_back( SpecialTime(GetDoubleFrom(*pMaturity)) );
  }
}
