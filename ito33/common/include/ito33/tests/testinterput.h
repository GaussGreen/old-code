#ifndef _ITO33_TEST_INTERPUT_H_
#define _ITO33_TEST_INTERPUT_H_

#include "ito33/common.h"
#include "ito33/cppunit.h"


// Test class
class InterpUtTest : public CppUnit::TestCase { 
public: 
  InterpUtTest( ) {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( InterpUtTest );
    CPPUNIT_TEST ( LinearSearch );
    CPPUNIT_TEST ( LinearInterp );
    CPPUNIT_TEST ( QuadraInt );
    CPPUNIT_TEST ( Spline );
    CPPUNIT_TEST ( Splint );

    CPPUNIT_TEST( LinearInterpolateUsingMatrix );
    CPPUNIT_TEST( QuadraticInterpolateUsingMatrix );

  CPPUNIT_TEST_SUITE_END();
  
  void LinearSearch();
  void LinearInterp();
  void QuadraInt();
  void Spline();
  void Splint();

  void LinearInterpolateUsingMatrix();
  void QuadraticInterpolateUsingMatrix();

  NO_COPY_CLASS( InterpUtTest );
};

#endif
