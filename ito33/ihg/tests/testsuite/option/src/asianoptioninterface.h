/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/testsuite/option/src/asianoptioninterface.h
// Purpose:     Base class for testing asian option
// Author:      Ito33Canada 
// Created:     2005/06/09
// RCS-ID:      $Id: asianoptioninterface.h,v 1.6 2006/08/22 10:10:43 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////


#ifndef _IHG_TESTS_TESTSUITE_OPTION_SRC_ASIANOPTIONINTERFACE_H_
#define _IHG_TESTS_TESTSUITE_OPTION_SRC_ASIANOPTIONINTERFACE_H_


#include "ito33/sharedptr.h"

#include "ito33/finance/exoticoption/asianoption.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "testparam.h"
#include "optioninterface.h"

namespace ito33 
{

namespace ihg 
{

namespace test
{


class AsianOptionInterface: public OptionInterface
{

protected:
  shared_ptr<finance::AsianOption>  m_pAsianOption;


public:
  AsianOptionInterface(shared_ptr<finance::SessionData> pSessionData,
    shared_ptr<finance::AsianOption> pAsianOption,
    shared_ptr<TheoreticalModel> pModel)
    : OptionInterface(pSessionData, pModel),
      m_pAsianOption(pAsianOption)
  {
    double dStrike = 1.0;
    if ( m_pAsianOption->GetFixedStrike() > 0 )
      dStrike = m_pAsianOption->GetFixedStrike();
     
    m_pOption = make_ptr( new finance::Option
        (dStrike, m_pAsianOption->GetMaturityDate(),
         m_pAsianOption->GetOptionType(), m_pAsianOption->GetExerciseType() ) );
  }

  virtual ~AsianOptionInterface() {}

  
  double GetStrike() const
  {
    return m_pAsianOption->GetFixedStrike();
  }

  void SetAverageStartDate(Date averageStartDate);

  shared_ptr<finance::Derivative> GetDerivative() const
  {
    return m_pAsianOption;
  }

  void Solve();

};



} //namespace test

} // namespace ihg

} // namespace ito33



#endif
