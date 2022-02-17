
#ifndef _ITO33_TEST_BISECNEWT_H_
#define _ITO33_TEST_BISECNEWT_H_

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"

// Test class
class BisecNewtTest : public CppUnit::TestCase { 
public: 
  BisecNewtTest( ) {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( BisecNewtTest );
    CPPUNIT_TEST ( Sin );
    CPPUNIT_TEST ( Polynom );
    CPPUNIT_TEST_EXCEPTION ( NegativeTolerance, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION ( InitEqual, ito33::Exception );
    CPPUNIT_TEST_EXCEPTION ( MaxIterations , ito33::Exception );
  CPPUNIT_TEST_SUITE_END();
  
  void Sin();
  void Polynom();
  void NegativeTolerance();
  void InitEqual();
  void MaxIterations();

  NO_COPY_CLASS( BisecNewtTest );
};

#endif
