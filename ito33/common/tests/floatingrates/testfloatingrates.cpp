/////////////////////////////////////////////////////////////////////////////
// Name:        tests/floatingrates/testfloatingrates.cpp
// Purpose:     test file for floating rates
// Author:      Nabil
// Created:     2005/09/08
// RCS-ID:      $Id: testfloatingrates.cpp,v 1.3 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"

#include "ito33/finance/floatingrates.h"

#include "ito33/array.h"
#include "ito33/list.h"
#include "ito33/vector.h"

#include "ito33/arraycheckers.h"

#include "ito33/cppunit.h"

#include "ito33/tests/testfloatingrates.h"

using namespace ito33;

using namespace ito33::finance;

void FloatingRatesTest::TooLowCap()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (m_pFloatingRatesInit);

  pFloatingRates->SetCap(0.);
}

void FloatingRatesTest::TooGreatCap()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (m_pFloatingRatesInit);

  pFloatingRates->SetCap(1.1);

}

void FloatingRatesTest::NegativeFloor()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (m_pFloatingRatesInit);

  pFloatingRates->SetFloor(-0.1);
}

// Test of the function validate()

void FloatingRatesTest::FirstUnknownPaymentDateAfterLastOne()
{
  InitMainData();

  m_firstUnknownPaymentDate = m_lastUnknownPaymentDate;
  m_firstUnknownPaymentDate.AddDays(1);

  shared_ptr<finance::FloatingRates> 
    pFloatingRates = InitFloatingRates();
}

void FloatingRatesTest::FirstUnknownPaymentDateBeforeStartOfAccruedOne()
{
  InitMainData();

  m_firstUnknownPaymentDate = m_startOfAccruedDate;
  m_firstUnknownPaymentDate.AddDays(-1);

  shared_ptr<finance::FloatingRates> 
    pFloatingRates = InitFloatingRates();
}

void FloatingRatesTest::UndefinedFrequency()
{
  InitMainData();

  m_paymentFrequency = finance::Frequency_Undefined;

  shared_ptr<finance::FloatingRates> 
    pFloatingRates = InitFloatingRates();
}

// Test of the function SetKnownPaymentRates(..)

void FloatingRatesTest::EmptyArrayOfKnownPaymentRates()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (m_pFloatingRatesInit);

  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  pFloatingRates->SetKnownPaymentStream(paymentDates, paymentRates);
}

void FloatingRatesTest::DecreasingDatesArrayOfKnownPaymentRates()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (m_pFloatingRatesInit);

  // Known payments added
  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  Date KnownPaymentDate;  
  double dAmount;
   
  KnownPaymentDate = 
    AddMonthsAdjustedForEndOfMonth(m_firstUnknownPaymentDate, - 1);
  dAmount = 0.02;
  paymentDates.push_back(KnownPaymentDate);
  paymentRates.push_back(dAmount);
   
  KnownPaymentDate = 
    AddMonthsAdjustedForEndOfMonth(m_firstUnknownPaymentDate, - 2);
  dAmount = 0.02;
  paymentDates.push_back(KnownPaymentDate);
  paymentRates.push_back(dAmount);

  pFloatingRates->SetKnownPaymentStream(paymentDates, paymentRates);
}

void FloatingRatesTest::NegativeKnownPaymentRates()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (m_pFloatingRatesInit);

  // Known payments added
  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  Date KnownPaymentDate;  
  double dAmount;
   
  KnownPaymentDate = 
    AddMonthsAdjustedForEndOfMonth(m_firstUnknownPaymentDate, - 1);
  dAmount = -0.001;
  paymentDates.push_back(KnownPaymentDate);
  paymentRates.push_back(dAmount);

  pFloatingRates->SetKnownPaymentStream(paymentDates, paymentRates);
}

void FloatingRatesTest::DifferentSizeForTheDataArrays()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (m_pFloatingRatesInit);

  // Known payments added
  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  Date KnownPaymentDate;  
  double dAmount;
   
  KnownPaymentDate = 
    AddMonthsAdjustedForEndOfMonth(m_firstUnknownPaymentDate, - 2);
  dAmount = 0.02;
  paymentDates.push_back(KnownPaymentDate);
  paymentRates.push_back(dAmount);
   
  KnownPaymentDate = 
    AddMonthsAdjustedForEndOfMonth(m_firstUnknownPaymentDate, - 1);
  dAmount = 0.02;
  paymentDates.push_back(KnownPaymentDate);
  paymentRates.push_back(dAmount);
  
  // without date !!!!!!
  paymentRates.push_back(dAmount);

  pFloatingRates->SetKnownPaymentStream(paymentDates, paymentRates);
}

void FloatingRatesTest::KnownPaymentAfterFirstUnknownOne()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (m_pFloatingRatesInit);

  // Known payments added
  std::vector<Date> paymentDates;
  std::vector<double> paymentRates;

  Date KnownPaymentDate;  
  double dAmount;
   
  KnownPaymentDate = m_firstUnknownPaymentDate;
  dAmount = 0.02;
  paymentDates.push_back(KnownPaymentDate);
  paymentRates.push_back(dAmount);

  pFloatingRates->SetKnownPaymentStream(paymentDates, paymentRates);
}

void FloatingRatesTest::EmptyFloatingRates_GetFirstUnknown()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (InitEmptyFloatingRates());

  pFloatingRates->GetFirstUnknownPaymentDate();
}

void FloatingRatesTest::EmptyFloatingRates_GetLastUnknown()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates (InitEmptyFloatingRates());

  pFloatingRates->GetLastUnknownPaymentDate();
}

void FloatingRatesTest::FloatingRatesWithJustOneUnknownPayment()
{
  shared_ptr<finance::FloatingRates> 
    pFloatingRates
    ( new finance::FloatingRates
          (m_dMargin, m_startOfAccruedDate, 
           m_firstUnknownPaymentDate, 
           m_firstUnknownPaymentDate,// the same as first one 
           m_paymentFrequency)
    );

  CPPUNIT_ASSERT_EQUAL (pFloatingRates->GetFirstUnknownPaymentDate(), 
                        pFloatingRates->GetLastUnknownPaymentDate());
}
