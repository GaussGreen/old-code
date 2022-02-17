/////////////////////////////////////////////////////////////////////////////
// Name:        tests/cashflowstream/testcashflowstream.cpp
// Purpose:     test file for cashflowstream
// Author:      Zhang (converted to cppunit by David)
// Modified:    Pedro (test for first and last coupons, minor changes)
// Created:     24.06.04
// RCS-ID:      $Id: testcashflowstream.cpp,v 1.22 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include <iostream>
#include <cmath>

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"

#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/cashflowstream_general.h"
#include "ito33/finance/cashflowstream_uniform.h"
#include "ito33/finance/frequency.h"

#include "ito33/list.h"
#include "ito33/vector.h"
#include "ito33/array.h"

#include "ito33/cppunit.h"

#include "ito33/tests/testcashflowstream.h"
#include "ito33/tests/utilexml.h"
#include "ito33/arraycheckers.h"

using namespace ito33;

using namespace ito33::finance;

void CashFlowStreamTest::BadPaymentDates()
{
  std::vector<Date> dates;
  std::vector<double> rates;

  dates.push_back(Date(200));
  dates.push_back(Date(200));
  
  rates.push_back(5);
  rates.push_back(5);

  CashFlowStreamGeneral c(Date(100), dates, rates, 
     Date::DayCountConvention_Act365,Frequency_BiMonthly);
}

void CashFlowStreamTest::NegativeAmount()
{
  double dAnnualAmt = -.25;

  CashFlowStreamUniform 
    csu(
         Date(2003, Date::Feb, 1),
         Date(2003, Date::May, 1),
         Date(2011, Date::May, 1),
         dAnnualAmt,
         Date::DayCountConvention_Act365,
         finance::Frequency_BiMonthly
       );
}

void CashFlowStreamTest::UndefinedFrequency()
{
  CashFlowStreamUniform 
    csu(
         Date(2003, Date::Feb, 1),
         Date(2003, Date::May, 1),
         Date(2011, Date::May, 1),
         0.5,
         Date::DayCountConvention_Act365,
         finance::Frequency_Undefined
       );
}

void CashFlowStreamTest::LastDateBeforeContractingDate()
{
  CashFlowStreamUniform 
    csu( 
         Date(2003, Date::Feb, 1),
         Date(2003, Date::May, 1),
         Date(2000, Date::May, 1),
         0.5,
         Date::DayCountConvention_Act365,
				 finance::Frequency_BiMonthly
       );
}
 
void CashFlowStreamTest::FirstDateBeforeContractingDate()
{                       
  CashFlowStreamUniform 
    csu(
         Date(2003, Date::Feb, 1),
         Date(2000, Date::May, 1),
         Date(2011, Date::May, 1),
         0.5,
         Date::DayCountConvention_Act365,
				 finance::Frequency_BiMonthly
       );
}

void CashFlowStreamTest::Uniform_Accrued()
{
  CashFlowStreamUniform c(Date(2000, Date::Feb, 5),
                          Date(2000, Date::May, 10),
                          Date(2010, Date::May, 10),
                          0.05,
                          Date::DayCountConvention_Act365,
                          Frequency_SemiAnnual
                          );

  double dAcc;
  dAcc = c.GetAccrued(Date(2001, Date::Nov, 10));
  CPPUNIT_ASSERT_EQUAL(dAcc, 0.);

  dAcc = c.GetAccrued(Date(1999, Date::Nov, 10));
  CPPUNIT_ASSERT_EQUAL(dAcc, 0.);

  dAcc = c.GetAccrued(Date(2010, Date::May, 10));
  CPPUNIT_ASSERT_EQUAL(dAcc, 0.);
}

void CashFlowStreamTest::Uniform_Accrued_From_SWX()
{
  // This function contains the examples of the SWX doc.

  Date accDate;

  double
    dExpectedAcc,
    dAcc;

  // Example of regular interest period with semi-annual coupon payments.
  CashFlowStreamUniform c(Date(1996, Date::Apr, 22),
                          Date(1996, Date::Oct, 22),
                          Date(1997, Date::Oct, 22),
                          2.75,
                          Date::DayCountConvention_ActAct_NO_EOM,
                          Frequency_SemiAnnual
                          );
  
  accDate = Date(1997, Date::Jan, 15);
  
  dExpectedAcc = (2.75 / 2.) * (85. / 182.);

  dAcc = c.GetAccrued(accDate);
  CPPUNIT_ASSERT_EQUAL(dExpectedAcc, dAcc);

  // Example of irregular (long) initial interest period with semi-annual 
  // coupon payments.
  c = CashFlowStreamUniform(Date(1995, Date::Aug, 15),
                            Date(1996, Date::Jul, 15),
                            Date(1997, Date::Jan, 15),
                            4.,
                            Date::DayCountConvention_ActAct_NO_EOM,
                            Frequency_SemiAnnual
                            );
  
  accDate = Date(1996, Date::Mar, 15);
  
  dExpectedAcc = (4. / 2) * (153. / 184. + 60. / 182.);

  dAcc = c.GetAccrued(accDate);
  CPPUNIT_ASSERT_EQUAL(dExpectedAcc, dAcc);

  // Example of irregular (short) initial interest period with annual 
  // coupon payments.
  c = CashFlowStreamUniform(Date(1999, Date::Feb, 1),
                            Date(1999, Date::Jul, 1),
                            Date(2000, Date::Jul, 1),
                            8.,
                            Date::DayCountConvention_ActAct_NO_EOM,
                            Frequency_Annual
                            );
  
  accDate = Date(1999, Date::May, 5);
  
  dExpectedAcc = (8. / 1.) * (93. / 365.);

  dAcc = c.GetAccrued(accDate);
  CPPUNIT_ASSERT_EQUAL(dExpectedAcc, dAcc);

  // Example of irregular (long) final interest period with quarterly 
  // coupon payments.
  c = CashFlowStreamUniform(Date(1998, Date::Aug, 31),
                            Date(1998, Date::Nov, 30),
                            Date(2000, Date::Apr, 30),
                            5.,
                            Date::DayCountConvention_ActAct,
                            Frequency_Quarterly,
                            LastPaymentType_Long
                            );
  
  accDate = Date(2000, Date::Apr, 3);
  
  dExpectedAcc = (5. / 4) * (91. / 91. + 34. / 92.);

  dAcc = c.GetAccrued(accDate);
  CPPUNIT_ASSERT_EQUAL(dExpectedAcc, dAcc);

  // Example of irregular (short) final interest period with semi-annual 
  // coupon payments.
  c = CashFlowStreamUniform(Date(1997, Date::Jun, 30),
                            Date(1998, Date::Jan, 30),
                            Date(2000, Date::Jun, 30),
                            6.,
                            Date::DayCountConvention_ActAct_NO_EOM,
                            Frequency_SemiAnnual,
                            LastPaymentType_Short
                            );
  
  accDate = Date(2000, Date::May, 15);
  
  dExpectedAcc = (6. / 2) * (106. / 182.);

  dAcc = c.GetAccrued(accDate);
  CPPUNIT_ASSERT_EQUAL(dExpectedAcc, dAcc);
}

void CashFlowStreamTest::Uniform_Irregular_Extreme_Period_From_ISMA()
{
  // This function contains the examples of irregular extreme periods 
  // from the ISMA doc.

  Date accDate;

  double
    dIrregCoupon,
    dRegCoupon;

  // Example of short first calculation period
  CashFlowStreamUniform c(Date(1999, Date::Feb, 1),
                          Date(1999, Date::Jul, 1),
                          Date(2001, Date::Jul, 1),
                          10,
                          Date::DayCountConvention_ActAct,
                          Frequency_Annual
                          );
  
  dIrregCoupon = 10. * (150. / (365. * 1.));

  dRegCoupon = 10.;

  CashFlowStream::const_iterator 
    iter,
    before_end_iter;

  for (iter = c.begin(); iter != c.end(); ++iter)
  { 
    if ( iter == c.begin() )
      CPPUNIT_ASSERT_DOUBLES_EQUAL (dIrregCoupon, iter->second, 1e-4);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL (dRegCoupon, iter->second, 1e-4);
  }

  // Example of long first calculation period
  c = CashFlowStreamUniform(Date(2002, Date::Aug, 15),
                            Date(2003, Date::Jul, 15),
                            Date(2005, Date::Jan, 15),
                            10,
                            Date::DayCountConvention_ActAct,
                            Frequency_SemiAnnual
                            );
  
  dIrregCoupon = 10. * ( 181. / (181. * 2.) + 153. / (184. * 2.) );

  dRegCoupon = 5;

  for (iter = c.begin(); iter != c.end(); ++iter)
  { 
    if ( iter == c.begin() )
      CPPUNIT_ASSERT_DOUBLES_EQUAL (dIrregCoupon, iter->second, 1e-4);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL (dRegCoupon, iter->second, 1e-4);
  }

  // Example of short final calculation period
  c = CashFlowStreamUniform(Date(1999, Date::Jul, 30),
                            Date(2000, Date::Jan, 30),
                            Date(2004, Date::Jun, 30),
                            10,
                            Date::DayCountConvention_ActAct,
                            Frequency_SemiAnnual,
                            LastPaymentType_Short
                            );
  
  dIrregCoupon = 10. * ( 152. / (182. * 2.)  );

  dRegCoupon = 5;

  before_end_iter = (--c.end());
  for (iter = c.begin(); iter != c.end(); ++iter)
  { 
    if ( iter == before_end_iter )
      CPPUNIT_ASSERT_DOUBLES_EQUAL (dIrregCoupon, iter->second, 1e-4);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL (dRegCoupon, iter->second, 1e-4);
  }

  // Example of long final calculation period
  c = CashFlowStreamUniform(Date(1999, Date::Aug, 31),
                            Date(1999, Date::Nov, 30),
                            Date(2004, Date::Apr, 30),
                            10,
                            Date::DayCountConvention_ActAct,
                            Frequency_Quarterly,
                            LastPaymentType_Long
                            );
  
  dIrregCoupon = 10. * ( 91 / (91. * 4.) + 61. / (92. * 4.) );

  dRegCoupon = 2.5;

  before_end_iter = (--c.end());
  for (iter = c.begin(); iter != c.end(); ++iter)
  { 
    if ( iter == before_end_iter )
      CPPUNIT_ASSERT_DOUBLES_EQUAL (dIrregCoupon, iter->second, 1e-4);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL (dRegCoupon, iter->second, 1e-4);
  }    
}

void CashFlowStreamTest::Uniform_LastRegular()
{
  CashFlowStreamUniform uniformcfs
    (
     Date(2003,Date::Dec,15),
     Date(2004,Date::Jan,1),
     Date(2010,Date::Jan,1),
     10.0,
     Date::DayCountConvention_ActAct,
     Frequency_SemiAnnual
    );
  
  int i=0;
  int year=2004;
  double amt=0.0;
  for (CashFlowStream::const_iterator pt=uniformcfs.begin(); pt != uniformcfs.end(); ++pt)
  {
    Date::Month month[2]={Date::Jan,Date::Jul};
    
    Date dt(year,month[i%2],1);
    amt=pt->second;
    Date dtcds=pt->first;
    if ( i > 0 )
      CPPUNIT_ASSERT_DOUBLES_EQUAL (amt,5.0,1e-4);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL (amt,0.461957,1e-4);

    CPPUNIT_ASSERT( dt == dtcds );

    year += i%2;
    i++;
  } 

  CPPUNIT_ASSERT( i == 13 );
}

void CashFlowStreamTest::Uniform_LastIrregular()
{
  CashFlowStreamUniform uniformcfs
    (
     Date(2003,Date::Dec,15),
     Date(2004,Date::Jan,1),
     Date(2010,Date::Jan,31),
     10.0,
     Date::DayCountConvention_Act365,
     Frequency_SemiAnnual,
     LastPaymentType_Short
    );

  int i=0;
  int year=2004;
  double amt=0.0;
  for (CashFlowStream::const_iterator pt=uniformcfs.begin(); pt != uniformcfs.end(); ++pt)
  {
    Date::Month month[2]={Date::Jan,Date::Jul};
    
    Date dt(year,month[i%2],1);
    Date dtcds=pt->first;
    amt=pt->second;
    if ( i > 0  && i != 13 )
      CPPUNIT_ASSERT_DOUBLES_EQUAL (amt,5.0,1e-4);
    else if ( i == 0)
      CPPUNIT_ASSERT_DOUBLES_EQUAL (amt,0.465753,1e-4);
    else
      CPPUNIT_ASSERT_DOUBLES_EQUAL (amt,0.821918,1e-4);

    if ( i != 13) 
      CPPUNIT_ASSERT( dt == dtcds );
    else
      CPPUNIT_ASSERT( Date(2010,Date::Jan,31) == pt->first );
    
    year += i%2;
    i++;
  } 

  CPPUNIT_ASSERT( i == 14 );
}

void CashFlowStreamTest::Uniform_EndOfMonth()
{
  CashFlowStreamUniform c(Date(2004, Date::Aug, 31),
                          Date(2004, Date::Nov, 30),
                          Date(2005, Date::Nov, 30),
                          0.04,
                          Date::DayCountConvention_Act365,
                          Frequency_Quarterly      
                          );

  Date dates[] = { Date(2004, Date::Nov, 30),
                   Date(2005, Date::Feb, 28),
                   Date(2005, Date::May, 31),
                   Date(2005, Date::Aug, 31),
                   Date(2005, Date::Nov, 30) };

  CashFlowStream::const_iterator pt; 
  size_t nIdx;
  for (pt = c.begin(), nIdx = 0; pt != c.end(); ++pt, ++nIdx)
  {
    CPPUNIT_ASSERT_EQUAL(dates[nIdx], pt->first);
    CPPUNIT_ASSERT_EQUAL(0.01, pt->second);
  }
}

void CashFlowStreamTest::Uniform_NoEndOfMonth()
{
  CashFlowStreamUniform c(Date(2004, Date::Aug, 30),
                          Date(2004, Date::Nov, 30),
                          Date(2005, Date::Nov, 30),
                          0.04,
                          Date::DayCountConvention_Act365_NO_EOM,
                          Frequency_Quarterly      
                          );

  Date dates[] = { Date(2004, Date::Nov, 30),
                   Date(2005, Date::Feb, 28),
                   Date(2005, Date::May, 30),
                   Date(2005, Date::Aug, 30),
                   Date(2005, Date::Nov, 30) };

  CashFlowStream::const_iterator pt; 
  size_t nIdx;
  for (pt = c.begin(), nIdx = 0; pt != c.end(); ++pt, ++nIdx)
  {
    CPPUNIT_ASSERT_EQUAL(dates[nIdx], pt->first);
    CPPUNIT_ASSERT_EQUAL(0.01, pt->second);
  }
}

void CashFlowStreamTest::Uniform_InMonths()
{
  CashFlowStreamUniform c(Date(2004, Date::Aug, 31),
                          15,
                          Frequency_Quarterly,
                          0.04,
                          Date::DayCountConvention_Act365);

  Date dates[] = { Date(2004, Date::Nov, 30),
                   Date(2005, Date::Feb, 28),
                   Date(2005, Date::May, 31),
                   Date(2005, Date::Aug, 31),
                   Date(2005, Date::Nov, 30) };

  CashFlowStream::const_iterator pt; 
  size_t nIdx;
  for (pt = c.begin(), nIdx = 0; pt != c.end(); ++pt, ++nIdx)
  {
    CPPUNIT_ASSERT_EQUAL(dates[nIdx], pt->first);
    CPPUNIT_ASSERT_EQUAL(0.01, pt->second);
  }
}


void CashFlowStreamTest::Uniform_AlwaysRegularButNotEndOfMonth()
{
  // Should not throw exception
  CashFlowStreamUniform c(Date(2004, Date::Dec, 30),
                          Date(2005, Date::Jun, 30),
                          Date(2010, Date::Dec, 30),
                          0.04,
                          Date::DayCountConvention_Act365_NO_EOM,
                          Frequency_SemiAnnual);
}

void CashFlowStreamTest::General_Accrued()
{
  // FIXME: I comment out this code (PF)
  //        as the ASSERTS are missing
  /*
  std::vector<Date> dates;
  std::vector<double> rates;

  dates.push_back(Date(2001, Date::May, 5));
  dates.push_back(Date(2001, Date::Nov, 5));
  dates.push_back(Date(2002, Date::May, 7));

  rates.push_back(.01);
  rates.push_back(0.02);
  rates.push_back(0.02);

  CashFlowStreamGeneral c(Date(2000, Date::Feb, 5),
                          dates, rates, Date::DayCountConvention_Act365);

  double t =  c.GetAccrued(Date(2001, Date::Jun, 1)); 
  */
}

void CashFlowStreamTest::FirstAndLastCoupon()
{
  double dAnnualAmt = .12;

  std::vector<finance::Frequency> freq;
  freq.push_back(Frequency_Annual);
  freq.push_back(Frequency_SemiAnnual);
  freq.push_back(Frequency_Quarterly);
  freq.push_back(Frequency_BiMonthly);
  freq.push_back(Frequency_Monthly);

  std::vector<Date::DayCountConvention> dcc;
  for (int i=0; i<Date::DayCountConvention_Max; i++)
    dcc.push_back((Date::DayCountConvention) i); 

  Date contracting=Date(2000,Date::Jan,1);
  Date lastDate=Date(2011,Date::Jan,1);

  for (std::vector<finance::Frequency>::const_iterator pf=freq.begin();
       pf != freq.end();
       ++pf)
  {
    for (std::vector<Date::DayCountConvention>::const_iterator pd=dcc.begin();
         pd != dcc.end();
         ++pd)
    {
      Date first=contracting;
      first.AddMonths(12 / ((int)(*pf)));
      Date last=lastDate;
      last.AddMonths(12 / ((int)(*pf)));
      
      CashFlowStreamUniform 
        csu(
            contracting,
            first,
            last,
            dAnnualAmt,
            *pd,
            *pf
           );
      CashFlow cf=csu.GetAll().front();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(cf.GetAmount(),dAnnualAmt/((int) *pf),1e-5);
      
      cf=csu.GetAll().back();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(cf.GetAmount(),dAnnualAmt/((int) *pf),1e-5);

      // TODO: Add a test for last coupon shorter and longer than the period
    }
  }
}

// This test is to check the leap year adjustment when the cash flow starts 
// at a 29. Data provided by W O'Connor.
// PF
// TODO: Add other frequencies
void CashFlowStreamTest::LeapYearAdjust()
{
  // This should work
  CashFlowStreamUniform cf(Date(2005, Date::Aug, 29),
                           Date(2006, Date::Feb, 28),
                           Date(2008, Date::Feb, 29),
                           0.04,
                           Date::DayCountConvention_Act360,
                           Frequency_SemiAnnual      
                          );

  CashFlow c=cf.GetAll().back();
  CPPUNIT_ASSERT_EQUAL(c.GetDate(),Date(2008, Date::Feb, 29));

  // Tests for all frequencies
  std::vector<finance::Frequency> freq;
  freq.push_back(Frequency_Annual);
  freq.push_back(Frequency_SemiAnnual);
  freq.push_back(Frequency_Quarterly);
  freq.push_back(Frequency_BiMonthly);
  freq.push_back(Frequency_Monthly);

  for (std::vector<finance::Frequency>::const_iterator pf=freq.begin();
       pf != freq.end();
       ++pf)
  {
    CashFlowStreamUniform cf(Date(2005, Date::Feb, 28),
                             Date(2006, Date::Feb, 28),
                             Date(2008, Date::Feb, 29),
                             0.04,
                             Date::DayCountConvention_Act360,
                             *pf      
                            );
    CashFlow c=cf.GetAll().back();
    CPPUNIT_ASSERT_EQUAL(c.GetDate(),Date(2008, Date::Feb, 29));
  }
}

void CashFlowStreamTest::LeapYearAdjustException()
{
  // This should throw
  CashFlowStreamUniform cf(Date(2005, Date::Aug, 29),
                           Date(2006, Date::Feb, 28),
                           Date(2008, Date::Feb, 28),
                           0.04,
                           Date::DayCountConvention_Act360,
                           Frequency_SemiAnnual      
                          );
}

void CashFlowStreamTest::Uniform_ShortLastPayment()
{
  Date 
    contractingDate(2005, Date::Jan, 01),
    firstPaymentDate(2006, Date::Jan, 01),
    regularLastPaymentDate(2016, Date::Jan, 01),
    irregularLastPaymentDate(2016, Date::May, 15);
  
  double dAnnualAmt = 0.04;

  std::vector<finance::Frequency> pFreq;
  pFreq.push_back(Frequency_Annual);
  pFreq.push_back(Frequency_SemiAnnual);
  pFreq.push_back(Frequency_Quarterly);
  pFreq.push_back(Frequency_BiMonthly);
  pFreq.push_back(Frequency_Monthly);

  std::vector<Date::DayCountConvention> pDCC;
  for (int i = 0; i<Date::DayCountConvention_Max; ++i)
    pDCC.push_back((Date::DayCountConvention) i); 

  std::vector<finance::Frequency>::const_iterator pf;
  
  std::vector<Date::DayCountConvention>::const_iterator pdcc;

  std::vector<CashFlow> 
    pCFS_short,
    pCFS_long,
    pCFS_uniform;

  size_t
    nIdx,
    nIdxFreq,
    nNbCFS_short,
    nNbCFS_long,
    nNbCFS_uniform;

  // Regular last payment date test
  for (pf = pFreq.begin(); pf != pFreq.end(); ++pf)
  {
    for (pdcc = pDCC.begin(); pdcc != pDCC.end(); ++pdcc)
    {      
      CashFlowStreamUniform 
        csu_short(
                   contractingDate,
                   firstPaymentDate,
                   regularLastPaymentDate,
                   dAnnualAmt,
                   *pdcc,
                   *pf,
                   LastPaymentType_Short
                 );

      CashFlowStreamUniform 
        csu_long(
                  contractingDate,
                  firstPaymentDate,
                  regularLastPaymentDate,
                  dAnnualAmt,
                  *pdcc,
                  *pf,
                  LastPaymentType_Long
                );

      CashFlowStreamUniform 
        csu(
             contractingDate,
             firstPaymentDate,
             regularLastPaymentDate,
             dAnnualAmt,
             *pdcc,
             *pf
            );
      
      // REMARK: The three cfs must be the same.

      pCFS_short = csu_short.GetAll();
      pCFS_long = csu_long.GetAll();
      pCFS_uniform = csu.GetAll();

      nNbCFS_short = pCFS_short.size();
      nNbCFS_long = pCFS_long.size();
      nNbCFS_uniform = pCFS_uniform.size();

      // Same size
      CPPUNIT_ASSERT_EQUAL(nNbCFS_short, nNbCFS_long);
      CPPUNIT_ASSERT_EQUAL(nNbCFS_long, nNbCFS_uniform);
      
      // Same cash flows
      for (nIdx = 0; nIdx < nNbCFS_uniform; ++nIdx)
      {
        // Same date
        CPPUNIT_ASSERT_EQUAL(pCFS_short[nIdx].GetDate(), 
                             pCFS_long[nIdx].GetDate());
        CPPUNIT_ASSERT_EQUAL(pCFS_long[nIdx].GetDate(), 
                             pCFS_uniform[nIdx].GetDate());
        // Same amount
        CPPUNIT_ASSERT_DOUBLES_EQUAL(pCFS_short[nIdx].GetAmount(), 
                                     pCFS_long[nIdx].GetAmount(),
                                     1.e-5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(pCFS_long[nIdx].GetAmount(), 
                                     pCFS_uniform[nIdx].GetAmount(),
                                     1.e-5);
      }
    }
  }

  // Irregular last payment date test (short and long)
  
  Date lastButOne_short[5] = {"2016/01/01", "2016/01/01", "2016/04/01", 
                              "2016/05/01", "2016/05/01"};

  Date lastButOne_long[5] = {"2015/01/01", "2015/07/01", "2016/01/01", 
                             "2016/03/01", "2016/04/01"};

  for (pf = pFreq.begin(), nIdxFreq = 0; pf != pFreq.end(); ++pf, ++nIdxFreq)
  {
    for (pdcc = pDCC.begin(); pdcc != pDCC.end(); ++pdcc)
    { 
      // Short case

      CashFlowStreamUniform 
        csu_short(
                   contractingDate,
                   firstPaymentDate,
                   irregularLastPaymentDate,
                   dAnnualAmt,
                   *pdcc,
                   *pf,
                   LastPaymentType_Short 
                 );

      // Long case

      CashFlowStreamUniform 
        csu_long(
                  contractingDate,
                  firstPaymentDate,
                  irregularLastPaymentDate,
                  dAnnualAmt,
                  *pdcc,
                  *pf,
                  LastPaymentType_Long
                );

      pCFS_short = csu_short.GetAll();
      pCFS_long = csu_long.GetAll();

      nNbCFS_short = pCFS_short.size();
      nNbCFS_long = pCFS_long.size();

      // Size test
      CPPUNIT_ASSERT_EQUAL(nNbCFS_short, nNbCFS_long + 1);

      // Same cash flows for the "long" case and the "short" case until 
      // the prelast one of the "long" case
      for (nIdx = 0; nIdx < nNbCFS_long - 1; ++nIdx)
      {
        // Same date
        CPPUNIT_ASSERT_EQUAL(pCFS_long[nIdx].GetDate(), 
                             pCFS_short[nIdx].GetDate());
        // Same amount
        CPPUNIT_ASSERT_DOUBLES_EQUAL(pCFS_long[nIdx].GetAmount(), 
                                     pCFS_short[nIdx].GetAmount(),
                                     1.e-5);
      }

      // Test for the last but one payment dates
      
      // short case
      CPPUNIT_ASSERT_EQUAL(pCFS_short[nNbCFS_short - 2].GetDate(), 
                           lastButOne_short[nIdxFreq]);
      
      // long case
      CPPUNIT_ASSERT_EQUAL(pCFS_long[nNbCFS_long - 2].GetDate(), 
                           lastButOne_long[nIdxFreq]);
    }
  }
}

void CashFlowStreamTest::Uniform_OnePayment()
{
  // This has only one coupon
  CashFlowStreamUniform cf(Date(2005, Date::Aug, 29),
                           Date(2006, Date::Jan, 28),
                           Date(2006, Date::Jan, 28),
                           0.04,
                           Date::DayCountConvention_Act360,
                           Frequency_Quarterly     
                          );

  CPPUNIT_ASSERT(cf.GetAll().size() == 1);
}

void CashFlowStreamTest::Uniform_TwoPayments()
{
  // This has two coupons
  CashFlowStreamUniform cf(Date(2005, Date::Aug, 29),
                           Date(2006, Date::Jan, 28),
                           Date(2006, Date::Feb, 28),
                           0.04,
                           Date::DayCountConvention_Act360,
                           Frequency_Quarterly,
                           LastPaymentType_Short
                          );

  CPPUNIT_ASSERT(cf.GetAll().size() == 2);

  // This has two coupons
  CashFlowStreamUniform cf2(Date(2005, Date::Aug, 29),
                            Date(2006, Date::Jan, 28),
                            Date(2006, Date::May, 28),
                            0.04,
                            Date::DayCountConvention_Act360,
                            Frequency_Quarterly,
                            LastPaymentType_Long
                           );

  CPPUNIT_ASSERT(cf2.GetAll().size() == 2);
}

void CashFlowStreamTest::Uniform_Dump()
{
  double 
    dAnnualAmt = 0.04;

  Date 
    contractingDate(2005, Date::Jan, 01),
    firstPaymentDate(2006, Date::Jan, 01),
    regularLastPaymentDate(2016, Date::Jan, 01),
    irregularLastPaymentDate(2016, Date::May, 15);

  Date::DayCountConvention dcc = Date::DayCountConvention_ActAct;

  Frequency freq = Frequency_Annual;

  std::ostringstream 
    oss_short,
    oss_long,
    oss_regular;

  // Short case
  ExpectedXML expected_short(oss_short,
    "<?xml version=\"1.0\"?>\n"
    "<root>\n"
    "<cash_distribution>\n"
    "<cash_flow_stream_uniform>\n"
    "<contracting_date>2005-01-01</contracting_date>\n"
    "<day_count_convention>actact</day_count_convention>\n"
    "<payment_frequency>annual</payment_frequency>\n"
    "<first_date>2006-01-01</first_date>\n"
    "<last_payment_type>short</last_payment_type>\n"
    "<last_date>2016-05-15</last_date>\n"
    "<annual_amount>0.04</annual_amount>\n"      
    "</cash_flow_stream_uniform>\n" 
    "</cash_distribution>\n"
    "</root>\n");
  
  shared_ptr<CashFlowStream> pcsu_short( new CashFlowStreamUniform(
                                              contractingDate,
                                              firstPaymentDate,
                                              irregularLastPaymentDate,
                                              dAnnualAmt,
                                              dcc,
                                              freq,
                                              LastPaymentType_Short 
                                            )
                                       );


  ito33::XML::RootTag root_short("root",oss_short);

  pcsu_short->Dump("cash_distribution", root_short);

  // Long case
  ExpectedXML expected_long(oss_long,
    "<?xml version=\"1.0\"?>\n"
    "<root>\n"
    "<cash_distribution>\n"
    "<cash_flow_stream_uniform>\n"
    "<contracting_date>2005-01-01</contracting_date>\n"
    "<day_count_convention>actact</day_count_convention>\n"
    "<payment_frequency>annual</payment_frequency>\n"
    "<first_date>2006-01-01</first_date>\n"
    "<last_payment_type>long</last_payment_type>\n"
    "<last_date>2016-05-15</last_date>\n"
    "<annual_amount>0.04</annual_amount>\n"      
    "</cash_flow_stream_uniform>\n" 
    "</cash_distribution>\n"
    "</root>\n");
  
  shared_ptr<CashFlowStream> pcsu_long( new CashFlowStreamUniform(
                                              contractingDate,
                                              firstPaymentDate,
                                              irregularLastPaymentDate,
                                              dAnnualAmt,
                                              dcc,
                                              freq,
                                              LastPaymentType_Long 
                                            )
                                       );


  ito33::XML::RootTag root_long("root",oss_long);

  pcsu_long->Dump("cash_distribution", root_long);

  // Regular case
  ExpectedXML expected_regular(oss_regular,
    "<?xml version=\"1.0\"?>\n"
    "<root>\n"
    "<cash_distribution>\n"
    "<cash_flow_stream_uniform>\n"
    "<contracting_date>2005-01-01</contracting_date>\n"
    "<day_count_convention>actact</day_count_convention>\n"
    "<payment_frequency>annual</payment_frequency>\n"
    "<first_date>2006-01-01</first_date>\n"
    "<last_date>2016-01-01</last_date>\n"
    "<annual_amount>0.04</annual_amount>\n"      
    "</cash_flow_stream_uniform>\n" 
    "</cash_distribution>\n"
    "</root>\n");
  
  shared_ptr<CashFlowStream> pcsu_regular( new CashFlowStreamUniform(
                                              contractingDate,
                                              firstPaymentDate,
                                              regularLastPaymentDate,
                                              dAnnualAmt,
                                              dcc,
                                              freq
                                            )
                                       );


  ito33::XML::RootTag root_regular("root",oss_regular);

  pcsu_regular->Dump("cash_distribution", root_regular);
}
