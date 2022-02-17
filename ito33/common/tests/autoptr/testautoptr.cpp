/////////////////////////////////////////////////////////////////////////////
// Name:        tests/autoptr/testautoptr.cpp
// Purpose:     implementation file of AutoPtr test program
// Author:      Vadim Zeitlin
// Created:     10.09.03
// RCS-ID:      $Id: testautoptr.cpp,v 1.2 2004/10/05 09:13:48 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"
#include "ito33/autoptr.h"

#include "ito33/list.h"
#include "ito33/vector.h"

#include "ito33/cppunit.h"

#include "ito33/tests/testautoptr.h"

using namespace ito33;

size_t Pointee::ms_nObjects = 0;


// check that we may declare an auto pointer for a forward declared class
class AutoPtrForwardTest;

AutoPtr<AutoPtrForwardTest> g_ptr;

class AutoPtrForwardTest
{
public:
  int n;
};

namespace ito33
{    
    ITO33_IMPLEMENT_AUTOPTR(Pointee);
    ITO33_IMPLEMENT_AUTOPTR(AutoPtrForwardTest);
}

// ----------------------------------------------------------------------------
// test class
// ----------------------------------------------------------------------------

void AutoPtrTestCase::Ctor()
{
  AutoPtr<Pointee> ptr(new Pointee);
}

void AutoPtrTestCase::CopyCtor()
{
  AutoPtr<Pointee> ptr1(new Pointee(17));
  AutoPtr<Pointee> ptr2(ptr1);

  CPPUNIT_ASSERT( ptr2 && !ptr1 );
  CPPUNIT_ASSERT( ptr2->GetValue() == 17 );
}

void AutoPtrTestCase::AssignmentOperator()
{
  AutoPtr<Pointee> ptr1(new Pointee(17));
  AutoPtr<Pointee> ptr2;
  ptr2 = ptr1;

  CPPUNIT_ASSERT( ptr2 && !ptr1 );
  CPPUNIT_ASSERT( ptr2->GetValue() == 17 );
}

void AutoPtrTestCase::Get()
{
  AutoPtr<Pointee> ptr(new Pointee(17));

  CPPUNIT_ASSERT( ptr->GetValue() == ptr.get()->GetValue() );

  ptr->SetValue(7);

  CPPUNIT_ASSERT( ptr->GetValue() == ptr.get()->GetValue() );

  ptr.get()->SetValue(9);

  CPPUNIT_ASSERT( ptr->GetValue() == ptr.get()->GetValue() );
}

void AutoPtrTestCase::OperatorBool()
{
  AutoPtr<Pointee> ptrInvalid;

  CPPUNIT_ASSERT( !ptrInvalid );

  AutoPtr<Pointee> ptrValid(new Pointee);

  CPPUNIT_ASSERT( ptrValid );
}

void AutoPtrTestCase::OperatorsStarAndArrow()
{
  AutoPtr<Pointee> ptr(new Pointee(17));

  CPPUNIT_ASSERT( (*ptr).GetValue() == ptr->GetValue() );

  (*ptr).SetValue(7);

  CPPUNIT_ASSERT( (*ptr).GetValue() == ptr->GetValue() );
}
