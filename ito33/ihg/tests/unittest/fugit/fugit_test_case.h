#include "ito33/common.h"
#include "ito33/cppunit.h"

class FugitTestCase : public CppUnit::TestCase
{
public:
  FugitTestCase() {}

private:
    CPPUNIT_TEST_SUITE( FugitTestCase );
        CPPUNIT_TEST( ZeroBeingExercisedCase1 );
    CPPUNIT_TEST_SUITE_END();

    void ZeroBeingExercisedCase1();

    NO_COPY_CLASS(FugitTestCase);
};

