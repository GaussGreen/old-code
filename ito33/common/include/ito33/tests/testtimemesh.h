/////////////////////////////////////////////////////////////////////////////
// Purpose:     header file for testing the time mesh generation
// Author:      David
// Created:     25.05.04
// RCS-ID:      $Id: testtimemesh.h,v 1.5 2004/11/23 14:22:16 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_TIMEMESH_H_
#define _ITO33_TEST_TIMEMESH_H_

#include "ito33/common.h"
#include "ito33/cppunit.h"
#include "ito33/numeric/mesh/specialtimes.h"

class TimeMeshTest : public CppUnit::TestCase
{
public:
  TimeMeshTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( TimeMeshTest );
    CPPUNIT_TEST( NoEvents );
    CPPUNIT_TEST( SingleEvent );
    CPPUNIT_TEST( TwoEvents );
    CPPUNIT_TEST( FiveEvents );
    CPPUNIT_TEST( MoreEventsThanTimes );
    CPPUNIT_TEST( CloseEvents );
  CPPUNIT_TEST_SUITE_END();

  // Test different numbers of events
  void NoEvents();
  void SingleEvent();
  void TwoEvents();
  void FiveEvents();

  // Check what happens when there are more events than times
  void MoreEventsThanTimes();

  // have events very close to each other
  void CloseEvents();

  // Helper function that checks the generated mesh
  void CheckMesh(ito33::numeric::mesh::SpecialTimes &pdSpecialTimes,
                 int iDirection,
                 std::vector<double> &pdTimes);

  // For debugging
  void PrintMesh(ito33::numeric::mesh::SpecialTimes &pdSpecialTimes,
                 std::vector<double> &pdTimes);

  NO_COPY_CLASS(TimeMeshTest);
};

#endif
