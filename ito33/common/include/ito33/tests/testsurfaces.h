/////////////////////////////////////////////////////////////////////////////
// Purpose:     main file of domain/surface test programs
// Author:      David Pooley
// Created:     17.05.04
// RCS-ID:      $Id: testsurfaces.h,v 1.4 2004/10/14 16:03:56 zhang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_TEST_SURFACE_H_
#define _ITO33_TEST_SURFACE_H_

#include "ito33/common.h"
#include "ito33/cppunit.h"

// ----------------------------------------------------------------------------
// domain_fixedspacemesh test class
// ----------------------------------------------------------------------------

class DomainFixedSpaceMeshTest : public CppUnit::TestCase
{
public:
  DomainFixedSpaceMeshTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( DomainFixedSpaceMeshTest );
    CPPUNIT_TEST( SinglePoint );
    CPPUNIT_TEST( MultiplePoints );
    CPPUNIT_TEST( GetSetSpots );
    CPPUNIT_TEST( AddTimes );
    CPPUNIT_TEST( DuplicatedDates );
    CPPUNIT_TEST( DateToTimeMapping );
  CPPUNIT_TEST_SUITE_END();

  // Single grid point
  void SinglePoint();

  // Full grid
  void MultiplePoints();

  // Setting/Getting spots at finance level
  void GetSetSpots();

  // Adding and getting times
  void AddTimes();

  // add 2 times that map to same date
  void DuplicatedDates();

  // check mapping from dates to times
  void DateToTimeMapping();

  NO_COPY_CLASS(DomainFixedSpaceMeshTest);
};

// ----------------------------------------------------------------------------
// SurfaceGeneral test class
// ----------------------------------------------------------------------------

class SurfaceGeneralTest : public CppUnit::TestCase
{
public:
  SurfaceGeneralTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( SurfaceGeneralTest );
 
    CPPUNIT_TEST( Empty );
    CPPUNIT_TEST( BasicAdd );
    CPPUNIT_TEST( MultipleAdds );
    CPPUNIT_TEST( Append );
    CPPUNIT_TEST( DeltaGamma );
   
    CPPUNIT_TEST( Theta ); 
    CPPUNIT_TEST( FiniteDifference );
  CPPUNIT_TEST_SUITE_END();

  // test empty surface
  void Empty();

  // Add one set of data
  void BasicAdd();

  // Add multiple data at one time
  void MultipleAdds();

  // Use adds to append data to the surface
  void Append();

  // Check delta and gamma calculations
  void DeltaGamma();

  // Check theta calculation
  void Theta();

  // Check finite difference calculation
  void FiniteDifference();


  NO_COPY_CLASS(SurfaceGeneralTest);
};


// ----------------------------------------------------------------------------
// DomainGeneral and SurfaceGeneral with DomainGeneral test class
// ----------------------------------------------------------------------------

class DomainGeneralWithSurfaceTest : public CppUnit::TestCase
{
public:

  DomainGeneralWithSurfaceTest() { }

  virtual void tearDown() { }

private:
  CPPUNIT_TEST_SUITE( DomainGeneralWithSurfaceTest );
 
    CPPUNIT_TEST( DomainGetSetSpots );
    CPPUNIT_TEST( DomainAddTimes );

    CPPUNIT_TEST( Empty );
    CPPUNIT_TEST( DeltaGamma );
   
    CPPUNIT_TEST( Theta ); 

    // we don't check FiniteDifference as it is independant of domain
    // CPPUNIT_TEST( FiniteDifference );

  CPPUNIT_TEST_SUITE_END();

  // Setting/Getting spots at finance level
  void DomainGetSetSpots();

  // Adding and getting times
  void DomainAddTimes();

  // test empty surface
  void Empty();

  // Check delta and gamma calculations
  void DeltaGamma();

  // Check theta calculation
  void Theta();



  NO_COPY_CLASS(DomainGeneralWithSurfaceTest);
};

#endif
