/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/forwardoption.cpp
// Purpose:     Implementation of ForwardOption class
// Author:      ITO33
// Created:     2005/05/05
// RCS-ID:      $Id: forwardoption.cpp,v 1.10 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/debug.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/option.h"
#include "ito33/finance/forwardoption.h"

namespace ito33
{

namespace finance
{

void ForwardOption::Add(const shared_ptr<Option>& pOption, double dWeight)
{
  ASSERT( pOption && pOption->GetExerciseType() == ExerciseType_European );

  m_elements.push_back( Element(pOption, dWeight) );
}

Date ForwardOption::GetMaturityDate() const
{
  ASSERT( !m_elements.empty() );

  Date maturityDate = m_elements.front().m_pOption->GetMaturityDate();
  Elements::const_iterator iter;
  for (iter = m_elements.begin(); iter != m_elements.end(); ++iter)
  {
    if ( iter->m_pOption->GetMaturityDate() > maturityDate)
      maturityDate = iter->m_pOption->GetMaturityDate();
  }

  return maturityDate;
}

const shared_ptr<SessionData>& ForwardOption::GetSessionData() const
{
  ASSERT( !m_elements.empty() );

  return m_elements.front().m_pOption->GetSessionData();
}


} // namespace finance

} // namespace ito33
