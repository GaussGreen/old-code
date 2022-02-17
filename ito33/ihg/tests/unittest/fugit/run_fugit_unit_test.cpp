#include "fugit_test_case.h"
#include "ihg/tests/run_fugit_unit_test.h"

int RunFugitUnitTest()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(FugitTestCase::suite());

    return runner.run("") ? 0 : 1;
}
