#include "edginc/config.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
using boost::unit_test::test_suite;
USING_DRLIB_NAMESPACE

#include "tests/SRMTests.hpp"

test_suite*  
init_unit_test_suite( int argc, char* argv[] )
{
    test_suite* test = BOOST_TEST_SUITE( "Master test suite" );
    test->add(new HyperTrigTestSuite);
    return test;
}                      