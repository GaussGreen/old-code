/////////////////////////////////////////////////////////////////////////////
// Purpose:     main file of domain/surface test programs
// Author:      David Pooley
// Created:     17.05.04
// RCS-ID:      $Id: testsurfaces.cpp,v 1.9 2006/08/19 23:22:41 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"

#include "ito33/dateutils.h"

#include "ito33/numeric/surfacegeneral.h"
#include "ito33/numeric/domain_general.h"
#include "ito33/numeric/domain_fixedspacemesh.h"

#include "ito33/tests/testsurfaces.h"

using namespace ito33;
using namespace ito33::numeric;

const bool g_bEndOfGrid = true;

// ----------------------------------------------------------------------------
// domain_fixedspacemesh test class
// ----------------------------------------------------------------------------

void DomainFixedSpaceMeshTest::SinglePoint()
{
  double dS = 100.0;
  DomainFixedSpaceMesh domain(dS);

  finance::Domain::Spots spots;

  // basic test
  spots = domain.GetSpaceMesh();
  CPPUNIT_ASSERT(spots.size() == 1);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS);

  // make sure all get mesh functions work
  size_t nIdx = 0;

  spots = domain.GetFirstSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == 1);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS);

  spots = domain.GetLastSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == 1);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS);

  spots = domain.GetOutputSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == 1);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS);

  // make sure the index does not matter for fixed grids
  nIdx = 1000;
  spots = domain.GetFirstSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == 1);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS);

  spots = domain.GetLastSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == 1);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS);

  spots = domain.GetOutputSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == 1);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS);

}

void DomainFixedSpaceMeshTest::MultiplePoints()
{
  size_t dSize = 5;
  double dS[] = {0.0, 50.0, 100.0, 200.0, 500.0};
  DomainFixedSpaceMesh domain(dS, dSize);

  finance::Domain::Spots spots;

  // basic test
  spots = domain.GetSpaceMesh();
  CPPUNIT_ASSERT(spots.size() == dSize);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[2], dS[2]);

  // make sure all get mesh functions work
  size_t nIdx = 0;

  spots = domain.GetFirstSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == dSize);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS[0]);

  spots = domain.GetLastSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == dSize);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[1], dS[1]);

  spots = domain.GetOutputSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == dSize);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[3], dS[3]);

  // make sure the index does not matter for fixed grids
  nIdx = 1000;
  spots = domain.GetFirstSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == dSize);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[4], dS[4]);

  spots = domain.GetLastSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == dSize);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[2], dS[2]);

  spots = domain.GetOutputSpaceMeshAt(nIdx);
  CPPUNIT_ASSERT(spots.size() == dSize);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[2], dS[2]);
  
}

void DomainFixedSpaceMeshTest::GetSetSpots()
{
  double dS = 100.0;
  DomainFixedSpaceMesh domain(dS);

  finance::Domain::Spots spotsOut;

  spotsOut = domain.GetSpaceMesh();
  CPPUNIT_ASSERT(spotsOut.size() == 1);

  double spotsIn[5];
  spotsIn[0] = 0.0;
  spotsIn[1] = 50.0;
  spotsIn[2] = 100.0;
  spotsIn[3] = 200.0;
  spotsIn[4] = 500.0;

  DomainFixedSpaceMesh domain2(spotsIn,5);

  spotsOut = domain2.GetSpaceMesh();
  CPPUNIT_ASSERT(spotsOut.size() == 5);
  ITO33_ASSERT_DOUBLES_EQUAL(spotsOut[2], spotsIn[2]);

}

void DomainFixedSpaceMeshTest::AddTimes()
{
  double dS = 100.0;
  DomainFixedSpaceMesh domain(dS);

  double dTime;

  domain.AddTime(100.0);
  dTime = domain.GetTimeAt(0);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 100.0);

  domain.AddTime(200.0);
  dTime = domain.GetTimeAt(0);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 100.0);
  dTime = domain.GetTimeAt(1);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 200.0);

  domain.AddTime(300.0);
  dTime = domain.GetTimeAt(0);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 100.0);
  dTime = domain.GetTimeAt(1);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 200.0);
  dTime = domain.GetTimeAt(2);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 300.0);

  // Check that no dates have been generated
  finance::Domain::Dates dates;
  dates = domain.GetDates();

  CPPUNIT_ASSERT(dates.size() == 0);

  // Generate dates, and check size
  domain.GenerateOutputDates();
  dates = domain.GetDates();
  CPPUNIT_ASSERT(dates.size() == 3);

  // Check that the dates are correct
  ITO33_ASSERT_DOUBLES_EQUAL(GetDoubleFrom(dates[0]), 100.0);
  ITO33_ASSERT_DOUBLES_EQUAL(GetDoubleFrom(dates[1]), 200.0);
  ITO33_ASSERT_DOUBLES_EQUAL(GetDoubleFrom(dates[2]), 300.0);

}


void DomainFixedSpaceMeshTest::DuplicatedDates()
{
  double dS = 100.0;
  DomainFixedSpaceMesh domain(dS);

  Date date1(2004, Date::Jan, 1);
  Date date2(2004, Date::Jan, 2);
  Date date3(2004, Date::Jan, 10);

  // Define 3 dates, but add 4 dates
  double dTime1 = GetDoubleFrom(date1);
  double dTime2 = GetDoubleFrom(date2);
  double dTime3 = GetDoubleFrom(date3);

  domain.AddTime(dTime1);
  domain.AddTime( 0.5 * (dTime1 + dTime2) );
  domain.AddTime(dTime2);
  domain.AddTime(dTime3);

  // make sure 4 times were added
  double dTime;
  dTime = domain.GetTimeAt(0);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, dTime1);
  dTime = domain.GetTimeAt(1);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 0.5*(dTime1+dTime2));
  dTime = domain.GetTimeAt(2);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, dTime2);
  dTime = domain.GetTimeAt(3);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, dTime3);

  // Generate dates, and check size
  domain.GenerateOutputDates();  
  finance::Domain::Dates dates = domain.GetDates();
  CPPUNIT_ASSERT(dates.size() == 3);

  // Check that the dates are correct
  CPPUNIT_ASSERT(dates[0] == date1);
  CPPUNIT_ASSERT(dates[1] == date2);
  CPPUNIT_ASSERT(dates[2] == date3);

}

void DomainFixedSpaceMeshTest::DateToTimeMapping()
{
  double dS = 100.0;
  DomainFixedSpaceMesh domain(dS);

  Date date1(2004, Date::Jan, 1);
  Date date2(2004, Date::Jan, 2);
  Date date3(2004, Date::Jan, 10);

  // Define 3 dates, but add 4 dates
  double dTime1 = GetDoubleFrom(date1);
  double dTime2 = GetDoubleFrom(date2);
  double dTime3 = GetDoubleFrom(date3);

  double dIncrement = dTime2 - dTime1;
  domain.AddTime(dTime1);
  domain.AddTime( dTime1 + 0.2 * dIncrement );
  domain.AddTime( dTime1 + 0.4 * dIncrement );
  domain.AddTime( dTime1 + 0.8 * dIncrement );
  domain.AddTime(dTime2);
  domain.AddTime( dTime2 + 0.5 * dIncrement );
  domain.AddTime(dTime3);

  // Verify output size
  domain.GenerateOutputDates();  
  finance::Domain::Dates dates = domain.GetDates();
  CPPUNIT_ASSERT(dates.size() == 3);

  // check mapping. Gives the last time at that date, not the first
  domain.GenerateOutputDates();
  size_t nIdx = domain.GetTimeIndexFromDateIndex(0);
  CPPUNIT_ASSERT(nIdx == 3);

  nIdx = domain.GetTimeIndexFromDateIndex(1);
  CPPUNIT_ASSERT(nIdx == 5);

  nIdx = domain.GetTimeIndexFromDateIndex(2);
  CPPUNIT_ASSERT(nIdx == 6);

}



// ----------------------------------------------------------------------------
// SurfaceGeneral test class
// ----------------------------------------------------------------------------

void SurfaceGeneralTest::Empty()
{
  double dS = 100.0;
  shared_ptr<Domain> domain(new DomainFixedSpaceMesh(dS));
  
  SurfaceGeneral surface(domain);  
}


void SurfaceGeneralTest::BasicAdd()
{

  Domain::Spots spots(3);
  spots[0] = 0.0;
  spots[1] = 100.0;
  spots[2] = 500.0;
  
  shared_ptr<DomainFixedSpaceMesh> 
    domain(new DomainFixedSpaceMesh(&spots[0], 3));  
  
  SurfaceGeneral surface(domain);  

  domain->AddTime(100.0);
  SurfaceDouble::Doubles values(3);
  values[0] = 1.0;
  values[1] = 2.0;
  values[2] = 6.0;

  surface.Add(values);

  // Make sure if gives back the exact values at the spots
  SurfaceDouble::Doubles valuesOut(3);
  surface.GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(values[0], valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(values[1], valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(values[2], valuesOut[2]);

  // Check interpolation
  spots[0] = 50.0;
  spots[1] = 200.0;
  spots[2] = 600.0;
  surface.GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.5, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(3.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(7.0, valuesOut[2]);

  // Make sure size of input does not matter
  // test smaller size than the grid
  spots.resize(2);
  spots[0] = 90.0;
  spots[1] = 150.0;
  valuesOut.resize(2);
  surface.GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.9, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(2.5, valuesOut[1]);

  // test larger size than the grid
  spots.resize(4);
  spots[0] = 80.0;
  spots[1] = 120.0;
  spots[2] = 330.0;
  spots[3] = 720.0;
  valuesOut.resize(4);
  surface.GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.8, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(2.2, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(4.3, valuesOut[2]);
  ITO33_ASSERT_DOUBLES_EQUAL(8.2, valuesOut[3]);

}


void SurfaceGeneralTest::MultipleAdds()
{
  Domain::Spots spots(3);
  SurfaceDouble::Doubles values1(3);
  SurfaceDouble::Doubles values2(3);
  SurfaceDouble::Doubles values3(3);

  // Setup domain
  spots[0] = 100.0; spots[1] = 200.0; spots[2] = 300.0;
  
  shared_ptr<DomainFixedSpaceMesh> 
    domain(new DomainFixedSpaceMesh(&spots[0], 3));  
  
  SurfaceGeneral surface(domain);  

  // Add one data member at first time
  domain->AddTime(100.0);
  values1[0] = 1.0; values1[1] = 2.0; values1[2] = 3.0;
  surface.Add(values1);

  // Add two data members at 2nd time
  domain->AddTime(200.0);
  values1[0] = 10.0; values1[1] = 20.0; values1[2] = 30.0;
  values2[0] = 11.0; values2[1] = 21.0; values2[2] = 31.0;
  surface.Add(values1);
  surface.Add(values2);

  // Add three data members at time 3
  domain->AddTime(300.0);
  values1[0] = 100.0; values1[1] = 200.0; values1[2] = 300.0;
  values2[0] = 110.0; values2[1] = 210.0; values2[2] = 310.0;
  values3[0] = 120.0; values3[1] = 220.0; values3[2] = 320.0;
  surface.Add(values1);
  surface.Add(values2);
  surface.Add(values3, g_bEndOfGrid);
  
  // Add two data members at 4th time
  domain->AddTime(400.0);
  values1[0] = 10.0; values1[1] = 20.0; values1[2] = 30.0;
  values2[0] = 11.0; values2[1] = 21.0; values2[2] = 31.0;
  surface.Add(values1);
  surface.Add(values2, g_bEndOfGrid);

  // Make sure everything can be extracted

  // Test all ways of getting the first set of data
  SurfaceDouble::Doubles valuesOut(3);

  surface.GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(3.0, valuesOut[2]);

  surface.GetFirstValuesAt( (size_t) 0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(3.0, valuesOut[2]);

  surface.GetLastValuesAt( (size_t) 0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(3.0, valuesOut[2]);

  // Repeat for the 2nd set of data
  surface.GetValuesAt(1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(11.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(21.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(31.0, valuesOut[2]);

  surface.GetFirstValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(10.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(20.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(30.0, valuesOut[2]);

  surface.GetLastValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(11.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(21.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(31.0, valuesOut[2]);

  // the third set of data, with three data values
  surface.GetValuesAt(2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(110.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(210.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(310.0, valuesOut[2]);

  surface.GetFirstValuesAt( (size_t) 2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(100.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(200.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(300.0, valuesOut[2]);

  surface.GetLastValuesAt( (size_t) 2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(120.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(220.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(320.0, valuesOut[2]);

  // Finally the 4th set of data
  surface.GetValuesAt(3, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(10.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(20.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(30.0, valuesOut[2]);

  surface.GetFirstValuesAt( (size_t) 3, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(10.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(20.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(30.0, valuesOut[2]);

  surface.GetLastValuesAt( (size_t) 3, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(11.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(21.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(31.0, valuesOut[2]);

}


void SurfaceGeneralTest::Append()
{
  Domain::Spots spots(3);
  SurfaceDouble::Doubles values1(3);
  SurfaceDouble::Doubles values2(3);

  // Setup domain
  spots[0] = 100.0; spots[1] = 200.0; spots[2] = 300.0;
  
  shared_ptr<DomainFixedSpaceMesh> 
    domain(new DomainFixedSpaceMesh(&spots[0], 3));  
  
  SurfaceGeneral surface(domain);  

  // Advance time in the domain, but add twice to the surface.  The second
  // should be an append
  domain->AddTime(100.0);
  values1[0] = 1.0; values1[1] = 2.0; values1[2] = 3.0;
  surface.Add(values1);
  values1[0] = 1.1; values1[1] = 2.1; values1[2] = 3.1;
  surface.Add(values1);

  // Advance time again, and use different add functions
  domain->AddTime(200.0);
  values1[0] = 10.0; values1[1] = 20.0; values1[2] = 30.0;
  values2[0] = 11.0; values2[1] = 21.0; values2[2] = 31.0;
  surface.Add(values1);
  surface.Add(values2);
  values1[0] = 10.1; values1[1] = 20.1; values1[2] = 30.1;
  surface.Add(values1, g_bEndOfGrid);

  // Repeat, using different add functions yet again
  domain->AddTime(300.0);
  values1[0] = 100.0; values1[1] = 200.0; values1[2] = 300.0;
  surface.Add(values1);
  values1[0] = 100.1; values1[1] = 200.1; values1[2] = 300.1;
  values2[0] = 110.1; values2[1] = 210.1; values2[2] = 310.1;
  surface.Add(values1);
  surface.Add(values2, g_bEndOfGrid);

  // Make sure everything can be extracted

  // Test all ways of getting the first set of data
  SurfaceDouble::Doubles valuesOut(3);

  surface.GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.1, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(2.1, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(3.1, valuesOut[2]);

  surface.GetFirstValuesAt( (size_t) 0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(3.0, valuesOut[2]);

  surface.GetLastValuesAt( (size_t) 0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.1, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(2.1, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(3.1, valuesOut[2]);

  // Repeat for the 2nd set of data
  surface.GetValuesAt(1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(11.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(21.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(31.0, valuesOut[2]);

  surface.GetFirstValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(10.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(20.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(30.0, valuesOut[2]);

  surface.GetLastValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(10.1, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(20.1, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(30.1, valuesOut[2]);

  // Finally, the third set of data, with three data values
  surface.GetValuesAt(2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(100.1, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(200.1, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(300.1, valuesOut[2]);

  surface.GetFirstValuesAt( (size_t) 2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(100.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(200.0, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(300.0, valuesOut[2]);

  surface.GetLastValuesAt( (size_t) 2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(110.1, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(210.1, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(310.1, valuesOut[2]);

}


void SurfaceGeneralTest::DeltaGamma()
{
  Domain::Spots spots(3);
  SurfaceDouble::Doubles values1(3);
  SurfaceDouble::Doubles values2(3);
  SurfaceDouble::Doubles values3(3);

  // Setup domain
  spots[0] = 0.0; spots[1] = 1.0; spots[2] = 5.0;
 
  shared_ptr<DomainFixedSpaceMesh> 
    domain(new DomainFixedSpaceMesh(&spots[0], 3));  

  SurfaceGeneral surface(domain);  

  // Setup test surface. Have one data set at time 1, 3 at time 2, 
  // and 2 at time 3
  domain->AddTime(100.0);
  values1[0] = 0.0; values1[1] = 1.0; values1[2] = 25.0; // x^2
  surface.Add(values1);

  domain->AddTime(101.0);
  values1[0] = 0.0; values1[1] = -1.0; values1[2] = -25.0; // -x^2
  values2[0] = 0.0; values2[1] = 1.0; values2[2] = 5.0;  // linear
  values3[0] = 2.0; values3[1] = 2.0; values3[2] = 2.0;  // const
  surface.Add(values1);
  surface.Add(values2);
  surface.Add(values3, g_bEndOfGrid);

  domain->AddTime(102.0);
  values1[0] = 1.0; values1[1] = 8.0; values1[2] = 214.0;  // x^3
  values2[0] = 0.0; values2[1] = 3.0; values2[2] = 35.0;   // x^2+2x
  surface.Add(values1);
  surface.Add(values2, g_bEndOfGrid);

  // compute delta and gamma
  shared_ptr<SurfaceDouble> surfaceDelta(new SurfaceGeneral(domain) );  
  shared_ptr<SurfaceDouble> surfaceGamma(new SurfaceGeneral(domain) );

  surface.GetDeltaAndGamma(surfaceDelta, surfaceGamma);

  // make sure it was calculated correctly (making some assumptions
  // about how delta and gamma are calculated...if those calculations
  // change, this function may have to be updated)
  SurfaceDouble::Doubles valuesOut(3);

  surfaceDelta->GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(6.0, valuesOut[2]); // backward

  surfaceDelta->GetFirstValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-1.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(-6.0, valuesOut[2]); // backward

  surfaceDelta->GetValuesAt(1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[2]); // backward

  surfaceDelta->GetLastValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[2]); // backward


  surfaceDelta->GetFirstValuesAt( (size_t) 2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(7.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL( (51.5 + 7.0*4.0)/5.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(206./4., valuesOut[2]); // backward

  surfaceDelta->GetLastValuesAt( (size_t) 2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(3.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(4.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(8.0, valuesOut[2]); // backward


  // now check gamma
  // end points are set to adjacent values, in this case value at 1
  surfaceGamma->GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[2]); 

  surfaceGamma->GetFirstValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[2]); 

  surfaceGamma->GetValuesAt(1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[2]); 

  surfaceGamma->GetLastValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[2]); 


  surfaceGamma->GetFirstValuesAt( (size_t) 2, spots, valuesOut);
  double dTmp = ((214. - 8.)/4.0 - (8. - 1.)/1.0)/2.5;
  ITO33_ASSERT_DOUBLES_EQUAL(dTmp, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(dTmp, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(dTmp, valuesOut[2]); 

  surfaceGamma->GetLastValuesAt( (size_t) 2, spots, valuesOut);  
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[2]); 

}


void SurfaceGeneralTest::Theta()
{
  Domain::Spots spots(3);
  SurfaceDouble::Doubles values1(3);
  SurfaceDouble::Doubles values2(3);
  SurfaceDouble::Doubles values3(3);

  // Setup domain
  spots[0] = 0.0; spots[1] = 1.0; spots[2] = 5.0;
  
  shared_ptr<DomainFixedSpaceMesh> 
    domain(new DomainFixedSpaceMesh(&spots[0], 3));  
  
  SurfaceGeneral surface(domain);  

  // Setup test surface. Have one data set at time 1, 3 at time 2, 
  // and 2 at time 3
  domain->AddTime(-100.0);
  values1[0] = 0.0; values1[1] = 1.0; values1[2] = 2.0;
  surface.Add(values1);

  domain->AddTime(-101.0);
  values1[0] = 1.0; values1[1] = 3.0; values1[2] = 5.0; 
  values2[0] = 3.0; values2[1] = 3.0; values2[2] = 3.0; 
  values3[0] = 2.0; values3[1] = 2.0; values3[2] = 2.0; 
  surface.Add(values1);
  surface.Add(values2);
  surface.Add(values3, g_bEndOfGrid);

  domain->AddTime(-104.0);
  values1[0] = 4.0; values1[1] = 6.0; values1[2] = 8.0; 
  values2[0] = 4.0; values2[1] = 4.0; values2[2] = 4.0;  
  surface.Add(values1);
  surface.Add(values2);

  // compute theta
  shared_ptr<SurfaceDouble> surfaceTheta(new SurfaceGeneral(domain) );  

  surface.GetThetaBackwardOnly(surfaceTheta);

  // make sure theta was calculated correctly (assuming finite 
  // difference calculations)

  // Only GetValuesAt() gives correct theta.
  // we don't need to check "first" and "last" values.
  SurfaceDouble::Doubles valuesOut(3);

  surfaceTheta->GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-1.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-3.0, valuesOut[2]); 

  surfaceTheta->GetValuesAt(1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-1.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-3.0, valuesOut[2]); 

  // last time, so backward differencing
  surfaceTheta->GetValuesAt(2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-2./3., valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-4./3., valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-6./3., valuesOut[2]);
}


void SurfaceGeneralTest::FiniteDifference()
{
  Domain::Spots spots(3);
  SurfaceDouble::Doubles values1(3);
  SurfaceDouble::Doubles values2(3);
  SurfaceDouble::Doubles values3(3);

  // Setup domain
  spots[0] = 0.0; spots[1] = 1.0; spots[2] = 5.0;
  
  shared_ptr<DomainFixedSpaceMesh>
    domain(new DomainFixedSpaceMesh(&spots[0], 3));  
  
  SurfaceGeneral surface(domain);  

  // Setup test surface. Have one data set at time 1, 3 at time 2, 
  // and 2 at time 3
  domain->AddTime(100.0);
  values1[0] = 0.0; values1[1] = 1.0; values1[2] = 2.0;
  surface.Add(values1);

  domain->AddTime(101.0);
  values1[0] = 1.0; values1[1] = 3.0; values1[2] = 5.0; 
  values2[0] = 3.0; values2[1] = 3.0; values2[2] = 3.0; 
  values3[0] = 2.0; values3[1] = 2.0; values3[2] = 2.0; 
  surface.Add(values1);
  surface.Add(values2);
  surface.Add(values3, g_bEndOfGrid);

  domain->AddTime(104.0);
  values1[0] = 4.0; values1[1] = 6.0; values1[2] = 8.0; 
  values2[0] = 4.0; values2[1] = 4.0; values2[2] = 4.0;  
  surface.Add(values1);
  surface.Add(values2);

  // Make a shifted surface
  double dShift = 0.111;
  
  shared_ptr<DomainFixedSpaceMesh> 
    shiftedDomain(new DomainFixedSpaceMesh(&spots[0], 3));  
  
  SurfaceGeneral shiftedSurface(shiftedDomain);  

  shiftedDomain->AddTime(100.0);
  values1[0] = 0.0 + dShift; values1[1] = 1.0 + dShift; values1[2] = 2.0 + dShift;
  shiftedSurface.Add(values1);

  shiftedDomain->AddTime(101.0);
  values1[0] = 1.0 + 2.0*dShift; values1[1] = 3.0 + dShift*dShift; values1[2] = 5.0 + 0.0; 
  values2[0] = 3.0 + 3.0*dShift; values2[1] = 3.0 + dShift/2.0;    values2[2] = 3.0 + 1.0; 
  values3[0] = 2.0 + 4.0*dShift; values3[1] = 2.0 + dShift*dShift; values3[2] = 2.0 + 2.0; 
  shiftedSurface.Add(values1);
  shiftedSurface.Add(values2);
  shiftedSurface.Add(values3, g_bEndOfGrid);

  shiftedDomain->AddTime(104.0);
  values1[0] = 4.0 + dShift*dShift; values1[1] = 6.0 + dShift/3.0; values1[2] = 8.0 + dShift; 
  values2[0] = 4.0 + dShift*dShift; values2[1] = 4.0 + dShift/2.0; values2[2] = 4.0 + dShift;  
  shiftedSurface.Add(values1);
  shiftedSurface.Add(values2);

  // Do the test
  shared_ptr<SurfaceDouble> 
    surfaceDerivs = surface.ComputeFiniteDifference(shiftedSurface, 1.0/dShift);  

  SurfaceDouble::Doubles valuesOut(3);

  surfaceDerivs->GetValuesAt(0, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[2]); 

  // 1st array at second time
  surfaceDerivs->GetFirstValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[0]);
  ITO33_ASSERT_DOUBLES_EQUAL(dShift, valuesOut[1]);
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[2]);

  // 2nd array at second time
  surfaceDerivs->GetValuesAt(1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(3.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(0.5, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL((4.0-3.0)/dShift, valuesOut[2]); 

  // 3rd array at second time
  surfaceDerivs->GetLastValuesAt( (size_t) 1, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(4.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(dShift, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL((4.0-2.0)/dShift, valuesOut[2]); 

  // 1st array at 3rd (last) time
  surfaceDerivs->GetFirstValuesAt( (size_t) 2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(dShift, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(1./3., valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[2]);

  // 2nd array at 3rd (last) time
  surfaceDerivs->GetLastValuesAt( (size_t) 2, spots, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(dShift, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(0.5, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[2]); 

  // Check the domain
  shared_ptr<Domain> derivDomain = surfaceDerivs->GetDomain();
  double dTime = domain->GetTimeAt(0);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 100.0);
  dTime = domain->GetTimeAt(1);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 101.0);
  dTime = domain->GetTimeAt(2);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 104.0);
}

// ----------------------------------------------------------------------------
// domain_fixedspacemesh test class
// ----------------------------------------------------------------------------
void DomainGeneralWithSurfaceTest::DomainGetSetSpots()
{
  DomainGeneral domain;
  size_t nIdx;
  finance::Domain::Spots spotsOut;

  finance::Domain::Spots spots;

  // time 100, add 1 spots
  const size_t n1 = 5;
  spots.resize(n1);
  for(nIdx = 0; nIdx < n1; nIdx++)
    spots[nIdx] = nIdx;
  domain.AddSpotsAtTime(spots, 100);

  // time 101, add 3 spots
  const size_t n20 = 5;
  spots.resize(n20);
  for(nIdx = 0; nIdx < n20; nIdx++)
    spots[nIdx] = nIdx + 1;
  domain.AddSpotsAtTime(spots, 101);
  const size_t n21 = 6;
  spots.resize(n21);
  for(nIdx = 0; nIdx < n21; nIdx++)
    spots[nIdx] = nIdx + 1.1;
  domain.AddSpotsAtTime(spots, 101);
  const size_t n22 = 7;
  spots.resize(n22);
  for(nIdx = 0; nIdx < n22; nIdx++)
    spots[nIdx] = nIdx + 1.2;
  domain.AddSpotsAtTime(spots, 101, g_bEndOfGrid);

  // time 102, add 2 spots
  const size_t n30 = 5;
  spots.resize(n30);
  for(nIdx = 0; nIdx < n30; nIdx++)
    spots[nIdx] = nIdx + 2;
  domain.AddSpotsAtTime(spots, 102);
  const size_t n31 = 6;
  spots.resize(n31);
  for(nIdx = 0; nIdx < n31; nIdx++)
    spots[nIdx] = nIdx + 2.1;
  domain.AddSpotsAtTime(spots, 102);

  // time 100
  CPPUNIT_ASSERT(domain.GetTimeAt(0) == 100);

  CPPUNIT_ASSERT(domain.GetFirstSpaceMeshAt(0).size() == n1);
  CPPUNIT_ASSERT(domain.GetOutputSpaceMeshAt(0).size() == n1);
  CPPUNIT_ASSERT(domain.GetLastSpaceMeshAt(0).size() == n1);

  CPPUNIT_ASSERT(domain.GetFirstSpaceMeshAt(0)[2] == 2);
  CPPUNIT_ASSERT(domain.GetOutputSpaceMeshAt(0)[3] == 3);
  CPPUNIT_ASSERT(domain.GetLastSpaceMeshAt(0)[4] == 4);
  
  // time 101
  CPPUNIT_ASSERT(domain.GetTimeAt(1) == 101);

  CPPUNIT_ASSERT(domain.GetFirstSpaceMeshAt(1).size() == n20);
  CPPUNIT_ASSERT(domain.GetOutputSpaceMeshAt(1).size() == n21);
  CPPUNIT_ASSERT(domain.GetLastSpaceMeshAt(1).size() == n22);
  
  CPPUNIT_ASSERT(domain.GetFirstSpaceMeshAt(1)[2] == 3);
  CPPUNIT_ASSERT(domain.GetOutputSpaceMeshAt(1)[3] == 4.1);
  CPPUNIT_ASSERT(domain.GetLastSpaceMeshAt(1)[4] == 5.2);

  // time 102
  CPPUNIT_ASSERT(domain.GetTimeAt(2) == 102);
  
  CPPUNIT_ASSERT(domain.GetFirstSpaceMeshAt(2).size() == n30);
  CPPUNIT_ASSERT(domain.GetOutputSpaceMeshAt(2).size() == n31);
  CPPUNIT_ASSERT(domain.GetLastSpaceMeshAt(2).size() == n31);
  
  CPPUNIT_ASSERT(domain.GetFirstSpaceMeshAt(2)[2] == 4);
  CPPUNIT_ASSERT(domain.GetOutputSpaceMeshAt(2)[3] == 5.1);
  CPPUNIT_ASSERT(domain.GetLastSpaceMeshAt(2)[4] == 6.1);
}

void DomainGeneralWithSurfaceTest::DomainAddTimes()
{
  finance::Domain::Spots spotsIn(5, 1);
  DomainGeneral domain;

  double dTime;

  // time 100
  domain.AddSpotsAtTime(spotsIn, 100.0);

  dTime = domain.GetTimeAt(0);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 100.0);

  // time 200
  domain.AddSpotsAtTime(spotsIn, 200.0);

  dTime = domain.GetTimeAt(0);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 100.0);
  dTime = domain.GetTimeAt(1);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 200.0);

  // continue
  domain.AddSpotsAtTime(spotsIn, 200.0);
  domain.AddSpotsAtTime(spotsIn, 200.0, g_bEndOfGrid);
  domain.AddSpotsAtTime(spotsIn, 300.0);

  dTime = domain.GetTimeAt(0);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 100.0);
  dTime = domain.GetTimeAt(1);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 200.0);
  dTime = domain.GetTimeAt(2);
  ITO33_ASSERT_DOUBLES_EQUAL(dTime, 300.0);

  // Check that no dates have been generated
  finance::Domain::Dates dates;
  dates = domain.GetDates();

  CPPUNIT_ASSERT(dates.size() == 0);

  // Generate dates, and check size
  domain.GenerateOutputDates();
  dates = domain.GetDates();
  CPPUNIT_ASSERT(dates.size() == 3);

  // Check that the dates are correct
  ITO33_ASSERT_DOUBLES_EQUAL(GetDoubleFrom(dates[0]), 100.0);
  ITO33_ASSERT_DOUBLES_EQUAL(GetDoubleFrom(dates[1]), 200.0);
  ITO33_ASSERT_DOUBLES_EQUAL(GetDoubleFrom(dates[2]), 300.0);
}


void DomainGeneralWithSurfaceTest::Empty()
{
  shared_ptr<Domain> domain(new DomainGeneral());
  
  SurfaceGeneral surface(domain);  
}


void DomainGeneralWithSurfaceTest::DeltaGamma()
{
  Domain::Spots
    spots1(3),
    spots2(3),
    spots3(4);
  SurfaceDouble::Doubles
    values1(3),
    values2(3),
    values3(4);

  // Setup spot mesh
  spots1[0] = 0.0; spots1[1] = 1.0; spots1[2] = 5.0;
  spots2[0] = 1.0; spots2[1] = 2.0; spots2[2] = 6.0;
  spots3[0] = 0.0; spots3[1] = 1.0; spots3[2] = 2.0; spots3[3] = 3.0;
 
  shared_ptr<DomainGeneral> domain(new DomainGeneral);  

  SurfaceGeneral surface(domain);  

  // Setup test surface. Have one data set at time 1, 3 at time 2, 
  // and 2 at time 3

  // time 100
  values1[0] = 0.0; values1[1] = 1.0; values1[2] = 25.0; // x^2
  domain->AddSpotsAtTime(spots1, 100.0);
  surface.Add(values1);

  // time 101
  values1[0] = 0.0; values1[1] = -1.0; values1[2] = -25.0; // -x^2
  values2[0] = 0.0; values2[1] = 1.0; values2[2] = 5.0;  // linear
  values3[0] = 0.0; values3[1] = 1.0; values3[2] = 4;  values3[3] = 9; // x^2
  domain->AddSpotsAtTime(spots1, 101.0);
  surface.Add(values1);
  domain->AddSpotsAtTime(spots2, 101.0);
  surface.Add(values2);
  domain->AddSpotsAtTime(spots3, 101.0, g_bEndOfGrid);
  surface.Add(values3, g_bEndOfGrid);

  values1[0] = 1.0; values1[1] = 8.0; values1[2] = 214.0;  // x^3
  values2[0] = 0.0; values2[1] = 3.0; values2[2] = 35.0;   // x^2+2x
  domain->AddSpotsAtTime(spots1, 102.0);
  surface.Add(values1);
  domain->AddSpotsAtTime(spots2, 102.0);
  surface.Add(values2);

  // compute delta and gamma
  shared_ptr<SurfaceDouble> surfaceDelta(new SurfaceGeneral(domain) );  
  shared_ptr<SurfaceDouble> surfaceGamma(new SurfaceGeneral(domain) );

  surface.GetDeltaAndGamma(surfaceDelta, surfaceGamma);

  // make sure it was calculated correctly (making some assumptions
  // about how delta and gamma are calculated...if those calculations
  // change, this function may have to be updated)
  SurfaceDouble::Doubles valuesOut(4);

  surfaceDelta->GetValuesAt(0, spots1, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(6.0, valuesOut[2]); // backward

  surfaceDelta->GetFirstValuesAt( (size_t) 1, spots1, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-1.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(-6.0, valuesOut[2]); // backward

  surfaceDelta->GetValuesAt(1, spots2, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[2]); // backward

  surfaceDelta->GetLastValuesAt( (size_t) 1, spots3, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(4.0, valuesOut[2]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(5.0, valuesOut[3]); // backward


  surfaceDelta->GetFirstValuesAt( (size_t) 2, spots1, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(7.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL( (51.5 + 7.0*4.0)/5.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(206./4., valuesOut[2]); // backward

  surfaceDelta->GetLastValuesAt( (size_t) 2, spots2, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(3.0, valuesOut[0]); // forward difference
  ITO33_ASSERT_DOUBLES_EQUAL(4.0, valuesOut[1]); // central
  ITO33_ASSERT_DOUBLES_EQUAL(8.0, valuesOut[2]); // backward


  // now check gamma
  // end points are set to adjacent values, in this case value at 1
  surfaceGamma->GetValuesAt(0, spots1, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[2]); 

  surfaceGamma->GetFirstValuesAt( (size_t) 1, spots1, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[2]); 

  surfaceGamma->GetValuesAt(1, spots2, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(0.0, valuesOut[2]); 

  surfaceGamma->GetLastValuesAt( (size_t) 1, spots3, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(2., valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2., valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2., valuesOut[2]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2., valuesOut[3]); 


  surfaceGamma->GetFirstValuesAt( (size_t) 2, spots1, valuesOut);
  double dTmp = ((214. - 8.)/4.0 - (8. - 1.)/1.0)/2.5;
  ITO33_ASSERT_DOUBLES_EQUAL(dTmp, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(dTmp, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(dTmp, valuesOut[2]); 

  surfaceGamma->GetLastValuesAt( (size_t) 2, spots2, valuesOut);  
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(2.0, valuesOut[2]); 
}


void DomainGeneralWithSurfaceTest::Theta()
{
  Domain::Spots spotsInit(3);
  SurfaceDouble::Doubles values1(3);
  SurfaceDouble::Doubles values2(3);
  SurfaceDouble::Doubles values3(3);

  // Setup domain
  spotsInit[0] = 0.0; spotsInit[1] = 1.0; spotsInit[2] = 5.0;
  
  shared_ptr<DomainGeneral> 
    domain(new DomainGeneral);  
  
  SurfaceGeneral surface(domain);  

  // Setup test surface. Have one data set at time 1, 3 at time 2, 
  // and 2 at time 3
  domain->AddSpotsAtTime(spotsInit, -100.0);
  values1[0] = 0.0; values1[1] = 1.0; values1[2] = 2.0;
  surface.Add(values1);

  values1[0] = 1.0; values1[1] = 3.0; values1[2] = 5.0; 
  domain->AddSpotsAtTime(spotsInit, -101.0);
  surface.Add(values1);

  // add a middle value
  // in fact, they have no effect to the result
  values2[0] = 3.0; values2[1] = 3.0; values2[2] = 3.0; 
  domain->AddSpotsAtTime(spotsInit, -101.0);
  surface.Add(values2);

  Domain::Spots spots2(3);
  spots2[0] = 1; spots2[1] = 2; spots2[2] = 3;
 
  values3[0] = 4.0; values3[1] = 5.0; values3[2] = 6.0;  // linear x + 3
  domain->AddSpotsAtTime(spots2, -101.0, g_bEndOfGrid);
  surface.Add(values3, g_bEndOfGrid);

  spots2.resize(2);
  spots2[0] = 1; spots2[1] = 2;
  values1[0] = 1; values1[1] = 2;  // linear  f(x) = x
  domain->AddSpotsAtTime(spots2, -104.0);
  surface.Add(values1);

  // compute theta
  shared_ptr<SurfaceDouble> surfaceTheta(new SurfaceGeneral(domain) );  

  surface.GetThetaBackwardOnly(surfaceTheta);

  // make sure theta was calculated correctly (assuming finite 
  // difference calculations)

  // Only GetValuesAt() gives correct theta.
  // we don't need to check "first" and "last" values.
  SurfaceDouble::Doubles valuesOut(3);

  surfaceTheta->GetValuesAt(0, spotsInit, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-1.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-3.0, valuesOut[2]); 

  surfaceTheta->GetValuesAt(1, spotsInit, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(-1.0, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-2.0, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(-3.0, valuesOut[2]); 

  surfaceTheta->GetValuesAt( 2, spotsInit, valuesOut);
  ITO33_ASSERT_DOUBLES_EQUAL(1, valuesOut[0]); 
  ITO33_ASSERT_DOUBLES_EQUAL(1, valuesOut[1]); 
  ITO33_ASSERT_DOUBLES_EQUAL(1, valuesOut[2]); 
}

