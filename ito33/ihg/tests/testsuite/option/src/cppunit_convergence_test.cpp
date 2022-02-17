
/////////////////////////////////////////////////////////////////////////////
// Name:       ihg/tests/testsuite/option/src/cppunit_option_test.cpp
// Purpose:     Base class for testing IHG projects
// Author:      Ito33Canada
// Created:     2005/06/13
// RCS-ID:      $Id: cppunit_convergence_test.cpp,v 1.4 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"

#include "ihg/xml/pricingreader.h"

#include "cppunit_convergence_test.h"
#include "convergencetests.h"
#include "testparam.h"

extern ito33::XML::RootTag convergenceRootTag;
extern ito33::ihg::test::TestParam testParam;

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33 
{
  

namespace ihg
{

namespace test
{

void CppUnitConvergenceTest::setup()
{
  XML::PricingReader reader( testParam.GetFileName() );
  shared_ptr<finance::SessionData> pSessionData(reader.ReadSessionData());

  finance::DerivativeVisitorGoodType visitor;
  reader.ReadDerivatives(visitor);

  if ( visitor.GetOption() )
  {
    m_nNbTests = 5;
    m_pDerivative = visitor.GetOption();
  }

  if ( visitor.GetAsianOption() )
  {
    m_nNbTests = 3;
    m_pDerivative = visitor.GetAsianOption();
  }
  
  if( !m_pDerivative )
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "Derivative type no recognized.");

  m_pModel = make_ptr( new ito33::ihg::TheoreticalModel );
  reader.ReadTheoreticalModel(m_pModel);

  m_pDerivative->SetSessionData( pSessionData );
}
  
void CppUnitConvergenceTest::CheckConvergence()
{
  setup();

  bool bResult = DerivativeCheckConvergence(m_pModel, m_pDerivative,
    testParam, convergenceRootTag, m_nNbTests);

  CPPUNIT_ASSERT( bResult );
}


} //end namespace test
}//end namespace ihg
}//end namespace ito33

