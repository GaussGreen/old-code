/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testfloatingrates.h
// Purpose:     header file for floating rates test
// Author:      Nabil
// Created:     2005/09/08
// RCS-ID:      $Id: testfloatingrates.h,v 1.3 2006/08/19 22:28:08 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/common.h"
#include "ito33/cppunit.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"
#include "ito33/exception.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/floatingrates.h"

namespace ito33
{

class FloatingRatesTest : public CppUnit::TestFixture 
{

public:

  FloatingRatesTest() 
  { 
    InitMainData(); 
    m_pFloatingRatesInit = InitFloatingRates(); 
  }

  void tearDown() {}

private:

  CPPUNIT_TEST_SUITE( FloatingRatesTest );

    CPPUNIT_TEST_EXCEPTION( TooLowCap, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( TooGreatCap, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( NegativeFloor, ito33::Exception );

    CPPUNIT_TEST_EXCEPTION( FirstUnknownPaymentDateBeforeStartOfAccruedOne, 
      ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( FirstUnknownPaymentDateAfterLastOne, 
      ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( UndefinedFrequency, ito33::Exception );

    CPPUNIT_TEST_EXCEPTION( EmptyArrayOfKnownPaymentRates, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( DecreasingDatesArrayOfKnownPaymentRates, 
      ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( NegativeKnownPaymentRates, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( DifferentSizeForTheDataArrays, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( KnownPaymentAfterFirstUnknownOne, 
      ito33::Exception );

    CPPUNIT_TEST_EXCEPTION( EmptyFloatingRates_GetFirstUnknown, 
      ito33::Exception );
    CPPUNIT_TEST_EXCEPTION( EmptyFloatingRates_GetLastUnknown, 
      ito33::Exception );

    CPPUNIT_TEST( FloatingRatesWithJustOneUnknownPayment );

  CPPUNIT_TEST_SUITE_END();

  void TooLowCap();
  void TooGreatCap();
  void NegativeFloor();

  // Test of the function Validate()

  void FirstUnknownPaymentDateBeforeStartOfAccruedOne();
  void FirstUnknownPaymentDateAfterLastOne();
  void UndefinedFrequency();

  // Test of the function SetKnownPaymentRates(...)

  void EmptyArrayOfKnownPaymentRates();
  void DecreasingDatesArrayOfKnownPaymentRates();
  void NegativeKnownPaymentRates();
  void DifferentSizeForTheDataArrays();
  void KnownPaymentAfterFirstUnknownOne();

  // Tests for an empty floating rates
  void EmptyFloatingRates_GetFirstUnknown();
  void EmptyFloatingRates_GetLastUnknown();

  // Create a floating rates with only one unkown payment must be possible
  void FloatingRatesWithJustOneUnknownPayment();

  // useful functions
  void InitMainData()
  {
    m_dMargin = 0.;

    m_startOfAccruedDate = Date(2002, Date::Jan, 14);    
    m_lastUnknownPaymentDate = Date(2004, Date::Jan, 14);
    m_firstUnknownPaymentDate = m_lastUnknownPaymentDate;
    
    m_paymentFrequency = finance::Frequency_Quarterly;    
    
    const int iInterval = 12 / m_paymentFrequency; 

    m_firstUnknownPaymentDate = 
      AddMonthsAdjustedForEndOfMonth(m_firstUnknownPaymentDate, 
        - 3 * iInterval);
  }

  shared_ptr<finance::FloatingRates> InitFloatingRates() const
  {      
    Date::DayCountConvention dcc = Date::DayCountConvention_30360;
    
    shared_ptr<finance::FloatingRates> pFloatingRates
           ( new finance::FloatingRates(m_dMargin, m_startOfAccruedDate, 
                   m_firstUnknownPaymentDate, m_lastUnknownPaymentDate, 
                   m_paymentFrequency)
           );

    pFloatingRates->SetDayCountConvention(dcc);

    return pFloatingRates;

  }

  shared_ptr<finance::FloatingRates> InitEmptyFloatingRates() const
  {          
    shared_ptr<finance::FloatingRates> 
      pFloatingRates
           ( new finance::FloatingRates(m_startOfAccruedDate,
                                        m_paymentFrequency)
           );

    return pFloatingRates;
  }

  double m_dMargin;

  Date m_startOfAccruedDate;    
  Date m_lastUnknownPaymentDate;
  Date m_firstUnknownPaymentDate;
  
  finance::Frequency m_paymentFrequency;
  
  shared_ptr<finance::FloatingRates> m_pFloatingRatesInit;

  NO_COPY_CLASS( FloatingRatesTest );

}; // Class FloatingRatesTest

} // namespace ito33

