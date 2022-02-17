#ifndef _CONST_CONSTRAINT_CASE_H_
#define _CONST_CONSTRAINT_CASE_H_

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"

#include "ito33/pricing/minconstconstraint.h"
#include "const_constraint_case.h"

#include "ito33/cppunit.h"

// normally there should be no "using" in header but this is just a test...
using ito33::pricing::MinConstConstraint;

class ConstConstraintCase : public CppUnit::TestCase
{
public:
  ConstConstraintCase() { }

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( ConstConstraintCase );
    CPPUNIT_TEST( Get );
    CPPUNIT_TEST( Apply );
    CPPUNIT_TEST( ApplyGet );
  CPPUNIT_TEST_SUITE_END();

  void Get();
  void Apply();
  void ApplyGet();

  MinConstConstraint m_cc;

  NO_COPY_CLASS(ConstConstraintCase);
};


#endif
