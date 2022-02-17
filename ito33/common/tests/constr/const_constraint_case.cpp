#include "ito33/pricing/constconstraint.h"

#include "ito33/tests/const_constraint_case.h"
#include "ito33/cppunit.h"

using namespace ito33;

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( ConstConstraintCase, "ConstConstraintCase" );

void ConstConstraintCase::Get()
{
  size_t nNb = 10;

  double pdP[50];
  int piFlag[50];
  double dConstraint = 100;
  size_t n;

  for(n = 0; n < nNb; n++)
    piFlag[n] = 1;

  // test Get() function when Constraint is On
  m_cc.Update(dConstraint);

  piFlag[0] = 0;

  m_cc.Get(pdP, piFlag, nNb);

  for(n = 1; n < nNb; n++)
  {
    CPPUNIT_ASSERT(piFlag[n] != 0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(pdP[n], dConstraint, 0 );
  }

  // test Get() function when Constraint is Off
  m_cc.TurnOff();
  m_cc.Get(pdP, piFlag, nNb);
  for(n = 1; n < nNb; n++)
  {
    CPPUNIT_ASSERT(piFlag[n] == 0);
  }
}

void ConstConstraintCase::Apply()
{
  size_t nNb = 50;

  double pdP[50];
  int piFlag[50];
  double dConstraint = 100;
  size_t n;

  for(n = 0; n < nNb; n++)
    piFlag[n] = 1;

  // test Apply() function when Constraint is On
  m_cc.Update(dConstraint);

  piFlag[0] = 0;

  nNb = 2;
  pdP[0] = 101;
  pdP[1] = 50;

  m_cc.Apply(pdP, piFlag, nNb);

  CPPUNIT_ASSERT(piFlag[0] == 0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pdP[0], 101, 0 );
  CPPUNIT_ASSERT(piFlag[1] != 0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pdP[1], dConstraint, 0 );

  // test Apply() function when Constraint is Off
  m_cc.TurnOff();
  nNb = 2;
  pdP[0] = 101;
  pdP[1] = 50;

  m_cc.Apply(pdP, piFlag, nNb);
  for(n = 0; n < nNb; n++)
  {
    CPPUNIT_ASSERT(piFlag[n] == 0);
  }
}


void ConstConstraintCase::ApplyGet()
{
  size_t nNb = 50;

  double pdP[50], pdNew[50];
  int piFlag[50];
  double dConstraint = 75.;
  size_t n;

  for(n = 0; n < nNb; n++)
    pdP[n] = 50. + n;

  m_cc.Update(dConstraint);
  m_cc.Apply(pdP, piFlag, nNb);
  m_cc.Get(pdNew, piFlag, nNb);

  for(n = 0; n < nNb; n++)
    if(piFlag[n])
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pdP[n], pdNew[n], 0);

  //m_cc.Update(dConstraint - 10);
  //m_cc.Get(pdNew, piFlag, nNb);

  //for(n = 0; n < nNb; n++)
  //  CPPUNIT_ASSERT(pdNew[n] <= pdP[n]);
}
