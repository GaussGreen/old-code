/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/option/testsuite/asianinterface.cpp
// Purpose:     Base class for testing IHG projects
// Author:      Ito33 Canada
// Created:     2005/06/13
// RCS-ID:      $Id: asianoptioninterface.cpp,v 1.9 2006/08/22 10:10:43 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"
#include "ito33/useexception.h"

#include "asianoptioninterface.h"
#include "testparam.h"
#include "utiltest.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33 
{

namespace ihg
{

namespace test
{

void AsianOptionInterface::SetAverageStartDate(Date averageStartDate)
{
  if ( m_pAsianOption->GetFixedStrike() > 0 )
  {
    m_pAsianOption = make_ptr( new finance::AsianOption
        ( m_pAsianOption->GetFixedStrike(),
          m_pAsianOption->GetMaturityDate(),
          m_pAsianOption->GetOptionType(), m_pAsianOption->GetExerciseType(),
          averageStartDate, m_pAsianOption->GetNbSamplingAverages() ) );
  }
  else
  {     
    m_pAsianOption = make_ptr( new finance::AsianOption
        ( m_pAsianOption->GetMaturityDate(),
          m_pAsianOption->GetOptionType(), m_pAsianOption->GetExerciseType(),
          averageStartDate, m_pAsianOption->GetNbSamplingAverages() ) );
  }
}

void AsianOptionInterface::Solve()
{
  m_pAsianOption->SetSessionData(m_pSessionData);

  shared_ptr<finance::ModelOutput> 
    pOutput =  m_pModel->Compute( *m_pAsianOption );

  m_dPrice = pOutput->GetPrice();
  m_dGamma = pOutput->GetGamma();
  m_dDelta = pOutput->GetDelta();
  m_dTheta = pOutput->GetTheta();
}

} //namespace test

} // namespace ihg

} // namespace ito33

