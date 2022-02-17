/////////////////////////////////////////////////////////////////////////////
// Purpose:     testing the time mesh generation
// Author:      David
// Created:     25.05.04
// RCS-ID:      $Id: testtimemesh.cpp,v 1.5 2006/01/05 14:19:58 yann Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>
#include <vector>
#include <math.h>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"
#include "ito33/dateutils.h"

#include "ito33/pricing/optionmeshmanager.h"

#include "ito33/numeric/mesh/specialtimes.h"

#include "ito33/tests/testtimemesh.h"

using namespace ito33;
using namespace ito33::pricing;
using namespace ito33::numeric::mesh;

extern void GeneralTimeMesh(SpecialTimes &pdSpecialTimes,
                     int iDirection,
                     size_t nMinNbTimeSteps,
                     std::vector<double> &pdTimes);


void TimeMeshTest::NoEvents()
{

  SpecialTimes specialTimes;
  specialTimes.clear();
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Jan, 1)) );
  specialTimes.push_back( GetDoubleFrom(Date(2005, Date::Jan, 1)) );

  std::vector<double> pdTimes;
  GeneralTimeMesh(specialTimes, 1, 50, pdTimes);

  CheckMesh(specialTimes, 1, pdTimes);

  /*
  CPPUNIT_ASSERT(spots.size() == 1);
  ITO33_ASSERT_DOUBLES_EQUAL(spots[0], dS);
  */
}

void TimeMeshTest::SingleEvent()
{

  // Test forward direction
  SpecialTimes specialTimes;
  specialTimes.clear();
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Jan, 1)) );
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Mar, 1)) );
  specialTimes.push_back( GetDoubleFrom(Date(2005, Date::Jan, 1)) );

  std::vector<double> pdTimes;
  GeneralTimeMesh(specialTimes, 1, 50, pdTimes);

  CheckMesh(specialTimes, 1, pdTimes);

  // Test backward direction
  specialTimes.reverse();
  GeneralTimeMesh(specialTimes, 0, 50, pdTimes);
  CheckMesh(specialTimes, 0, pdTimes);

}

void TimeMeshTest::TwoEvents()
{

  SpecialTimes specialTimes;
  specialTimes.clear();
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Jan, 1)) );
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Mar, 1)) );
  specialTimes.push_back( GetDoubleFrom(Date(2005, Date::Feb, 1)) );
  specialTimes.push_back( GetDoubleFrom(Date(2005, Date::Jul, 1)) );

  std::vector<double> pdTimes;
  GeneralTimeMesh(specialTimes, 1, 50, pdTimes);

  CheckMesh(specialTimes, 1, pdTimes);

  // Test backward direction
  specialTimes.reverse();
  GeneralTimeMesh(specialTimes, 0, 50, pdTimes);
  CheckMesh(specialTimes, 0, pdTimes);

}

void TimeMeshTest::FiveEvents()
{

  SpecialTimes specialTimes;
  specialTimes.clear();
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Jan, 1)) );
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Mar, 1)) );
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Jun, 2)) );
  specialTimes.push_back( GetDoubleFrom(Date(2004, Date::Nov, 3)) );
  specialTimes.push_back( GetDoubleFrom(Date(2005, Date::Jan, 4)) );
  specialTimes.push_back( GetDoubleFrom(Date(2005, Date::Feb, 5)) );
  specialTimes.push_back( GetDoubleFrom(Date(2005, Date::Jul, 1)) );

  std::vector<double> pdTimes;
  GeneralTimeMesh(specialTimes, 1, 50, pdTimes);

  CheckMesh(specialTimes, 1, pdTimes);

  // Test backward direction
  specialTimes.reverse();
  GeneralTimeMesh(specialTimes, 0, 50, pdTimes);
  CheckMesh(specialTimes, 0, pdTimes);
}


void TimeMeshTest::MoreEventsThanTimes()
{
  SpecialTimes specialTimes;
  specialTimes.clear();

  // Create 100 special points
  Date date(2000, Date::Jan, 1);
  for (size_t nIdx = 0; nIdx < 100; nIdx++)
    specialTimes.push_back( GetDoubleFrom(date.AddDays(2)) );

  // Request 50 points
  std::vector<double> pdTimes;
  GeneralTimeMesh(specialTimes, 1, 50, pdTimes);
  CheckMesh(specialTimes, 1, pdTimes);

}


void TimeMeshTest::CloseEvents()
{
  SpecialTimes specialTimes;
  specialTimes.clear();

  // Use a start time typical for the code
  double dStart = GetDoubleFrom( Date(2000, Date::Jan, 1) );
  double dOneDay = GetDoubleFrom( Date(2000, Date::Jan, 2) )
                 - dStart;
  specialTimes.push_back(dStart);
  specialTimes.push_back(dStart + 50.0*dOneDay);
  specialTimes.push_back(dStart + 50.0*dOneDay + 1.e-5);
  specialTimes.push_back(dStart + 50.0*dOneDay + 2.e-5);
  specialTimes.push_back(dStart + 120.0*dOneDay);
  specialTimes.push_back(dStart + 120.0*dOneDay + 1.e-5);
  specialTimes.push_back(dStart + 365.0*dOneDay);

  // Request 50 points
  std::vector<double> pdTimes;
  GeneralTimeMesh(specialTimes, 1, 50, pdTimes);
  CheckMesh(specialTimes, 1, pdTimes);

  // Test backward direction
  specialTimes.reverse();
  GeneralTimeMesh(specialTimes, 0, 50, pdTimes);
  CheckMesh(specialTimes, 0, pdTimes);
}


void TimeMeshTest::CheckMesh(SpecialTimes& pdSpecialTimes,
                             int iDirection,
                             std::vector<double>& pdTimes)
{
  // Must have at least pricing date and maturity date
  size_t nNbTimes = pdTimes.size();
  CPPUNIT_ASSERT(nNbTimes >= 2);
  CPPUNIT_ASSERT(pdSpecialTimes.size() >= 2);
  
  // Make sure the times are in order, and all special times are
  // in the mesh
  size_t nIdx;
  SpecialTimes::iterator iterSpecial = pdSpecialTimes.begin();
  for (nIdx = 0; nIdx < nNbTimes; nIdx++)
  {
    if ( fabs( (*iterSpecial).GetTime() - pdTimes[nIdx] ) < TIMETOLERANCE )
      ++iterSpecial;

    if (nIdx > 0)
      if (iDirection)
        CPPUNIT_ASSERT( pdTimes[nIdx] > pdTimes[nIdx-1] );
      else
        CPPUNIT_ASSERT( pdTimes[nIdx] < pdTimes[nIdx-1] );
  }

  CPPUNIT_ASSERT( iterSpecial == pdSpecialTimes.end() );

}


void TimeMeshTest::PrintMesh(SpecialTimes &pdSpecialTimes,
                             std::vector<double> &pdTimes)
{
  std::cout.precision(10);
  size_t nIdx;
  SpecialTimes::iterator iterSpecial = pdSpecialTimes.begin();
  for (nIdx = 0; nIdx < pdTimes.size(); nIdx++)
  {
    if ( fabs( (*iterSpecial).GetTime() - pdTimes[nIdx] ) < TIMETOLERANCE )
    {
      std::cout << "* " << pdTimes[nIdx];
      ++iterSpecial;
    }
    else
    {
      std::cout << "  " << pdTimes[nIdx];
    }
    if (nIdx > 0)
      std::cout << "  " << pdTimes[nIdx] - pdTimes[nIdx-1] << std::endl;
    else
      std::cout << std::endl;
  }
  std::cout << std::endl << "num points = " << pdTimes.size() << std::endl;

}
