/////////////////////////////////////////////////////////////////////////////
// Name:        tests/testxmlreader_bondlike.h
// Purpose:     header file for bondlike test
// Author:      Zhang (converted to cppunit by David)
// Created:     09.03.2004
// RCS-ID:      $Id: testxmlreader_bondlike.h,v 1.7 2006/05/24 12:35:37 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/common.h"
#include "ito33/cppunit.h"
#include "ito33/xml/read.h"


class XMLReaderBondLikeTest : public CppUnit::TestCase
{ 

public: 
  XMLReaderBondLikeTest() {}

  void tearDown() {}

private:
  CPPUNIT_TEST_SUITE( XMLReaderBondLikeTest );

    //-----------------//
    // bond like terms //
    //-----------------//
    CPPUNIT_TEST           ( CheckBondLikeTerms );
    CPPUNIT_TEST           ( CheckBondTerms );
    CPPUNIT_TEST_EXCEPTION ( CheckBondTermsWrongRedemption,           
                                      ito33::XML::Exception );

    //---------------//
    // Call Schedule //
    //---------------//
    CPPUNIT_TEST           ( CheckCallSchedule );

    //---------------------//
    // Conversion Schedule //
    //---------------------//
    CPPUNIT_TEST           ( CheckConversionSchedule );

  CPPUNIT_TEST_SUITE_END();

//_____________________________________________________________________________
    void CheckBondLikeTerms();
    void CheckBondTerms();
    void CheckBondTermsWrongRedemption();

    void CheckCallSchedule();

    void CheckConversionSchedule();

  NO_COPY_CLASS( XMLReaderBondLikeTest );
}; // Class XMLReaderBondLikeTest

