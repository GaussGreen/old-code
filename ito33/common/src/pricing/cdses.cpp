/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cdses.cpp
// Author:      David
// Created:     2004/03/31
// RCS-ID:      $Id: cdses.cpp,v 1.14 2006/08/21 16:10:00 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/list.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/constants.h"
#include "ito33/useexception.h"

#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/cdslike.h"

#include "ito33/pricing/cdses.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::pricing;

using namespace ito33::numeric::mesh;

CDSes::CDSes(const std::list< shared_ptr<ito33::finance::CDSLike> >& pCDSes)
: m_pSpreadStream( pCDSes.back()->GetSpreadStream() )
{
  m_pMaturityDates.clear();

  std::list< shared_ptr<ito33::finance::CDSLike> >::const_iterator pCDS;

  pCDS = pCDSes.begin();

  for (; pCDS != pCDSes.end(); ++pCDS)
  { 
    CashFlowStream::const_iterator 
      pPayment = (**pCDS).GetSpreadStream()->begin();

    for ( ; pPayment != (**pCDS).GetSpreadStream()->end(); ++pPayment)
    {
      m_pMaturityDates.push_back( pPayment->first );
    }

    //m_pMaturityDates.push_back( (*pCDS)->GetMaturityDate() );
  }

  m_dRecoveryRate = pCDSes.front()->GetRecoveryRate(); 
}


CDSes::CDSes(const std::list< ito33::finance::CDSLike* >& pCDSes)
: m_pSpreadStream( pCDSes.back()->GetSpreadStream() )
{
  m_pMaturityDates.clear();

  std::list< ito33::finance::CDSLike* >::const_iterator pCDS;

  pCDS = pCDSes.begin();

  for (; pCDS != pCDSes.end(); ++pCDS)
  { 
    CashFlowStream::const_iterator 
      pPayment = (**pCDS).GetSpreadStream()->begin();

    for ( ; pPayment != (**pCDS).GetSpreadStream()->end(); ++pPayment)
    {
      m_pMaturityDates.push_back( pPayment->first );
    }

    //m_pMaturityDates.push_back( (*pCDS)->GetMaturityDate() );
  }

  m_dRecoveryRate = pCDSes.front()->GetRecoveryRate();
}


void CDSes::GetSpecialTimes(SpecialTimes& specialTimes) const
{
  std::list<Date>::const_iterator pMaturity;

  for (pMaturity = m_pMaturityDates.begin(); 
       pMaturity != m_pMaturityDates.end();
       ++pMaturity)
  {
    specialTimes.push_back( SpecialTime(GetDoubleFrom(*pMaturity)) );
  }
}
