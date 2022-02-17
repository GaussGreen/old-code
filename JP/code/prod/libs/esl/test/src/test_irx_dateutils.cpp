#if 0
#include "boost/test/included/unit_test_framework.hpp"
#endif

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

using boost::unit_test::test_suite;

#include "irx/dateutils.h"

// Tolerance is specified in percent
#define TOL_PCT 1.0E-13

static void test_irxDayCountFractionActAct() {
    IrxTDate    date1;
    IrxTDate    date2;
    int rc;
    double dcf;

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2001, 1, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 366.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2000, 1, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 0.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2000, 1, 2);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2000, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 28);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/366.0, TOL_PCT);

    date1 = irxDate(2001, 2, 28);
    date2 = irxDate(2001, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/365.0, TOL_PCT);

    date1 = irxDate(2003, 11, 1);
    date2 = irxDate(2004, 5, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 61.0/365.0 + 121.0/366.0, TOL_PCT);

    date1 = irxDate(1999, 2, 1);
    date2 = irxDate(1999, 7, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 150.0/365.0, TOL_PCT);

    date1 = irxDate(1999, 7, 1);
    date2 = irxDate(2000, 7, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 184.0/365.0 + 182.0/366.0, TOL_PCT);

    date1 = irxDate(2002, 8, 15);
    date2 = irxDate(2003, 7, 15);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 139.0/365.0 + 195.0/365.0, TOL_PCT);

    date1 = irxDate(2003, 7, 15);
    date2 = irxDate(2004, 1, 15);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 170.0/365.0 + 14.0/366.0, TOL_PCT);

    date1 = irxDate(1999, 7, 30);
    date2 = irxDate(2000, 1, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 155.0/365.0 + 29.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 1, 30);
    date2 = irxDate(2000, 6, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 152.0/366.0, TOL_PCT);

    date1 = irxDate(1999, 11, 30);
    date2 = irxDate(2000, 4, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 32.0/365.0 + 120.0/366.0, TOL_PCT);

    date1 = irxDate(2001, 3, 1);
    date2 = irxDate(2001, 2, 28);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -1.0/365.0, TOL_PCT);

    date1 = irxDate(2000, 3, 1);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -1.0/366.0, TOL_PCT);

    date1 = irxDate(2001, 3, 1);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -307.0/366.0 + -59.0/365.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2001, 2, 28);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 307.0/366.0 + 58.0/365.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2004, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 4.0, TOL_PCT);

    date1 = irxDate(1999, 2, 28);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 307.0/365.0 + 59.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2004, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 3.0 + 307.0/366.0 + 60.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2003, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 307.0/366.0 + 789.0/365.0, TOL_PCT);

    date1 = irxDate(2000, 2, 28);
    date2 = irxDate(2004, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 3.0 + 308.0/366.0 + 60.0/366.0, TOL_PCT);
}

static void test_irxDayCountFractionActActAFB() {
    IrxTDate    date1;
    IrxTDate    date2;
    int rc;
    double dcf;

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2001, 1, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0, TOL_PCT);

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2000, 1, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 0.0, TOL_PCT);

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2000, 1, 2);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/365.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2000, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 28);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/365.0, TOL_PCT);

    date1 = irxDate(2001, 2, 28);
    date2 = irxDate(2001, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/365.0, TOL_PCT);

    date1 = irxDate(2003, 11, 1);
    date2 = irxDate(2004, 5, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 182.0/366.0, TOL_PCT);

    date1 = irxDate(1999, 2, 1);
    date2 = irxDate(1999, 7, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 150.0/365.0, TOL_PCT);

    date1 = irxDate(1999, 7, 1);
    date2 = irxDate(2000, 7, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 366.0/366.0, TOL_PCT);

    date1 = irxDate(2002, 8, 15);
    date2 = irxDate(2003, 7, 15);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 334.0/365.0, TOL_PCT);

    date1 = irxDate(2003, 7, 15);
    date2 = irxDate(2004, 1, 15);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 184.0/365.0, TOL_PCT);

    date1 = irxDate(1999, 7, 30);
    date2 = irxDate(2000, 1, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 184.0/365.0, TOL_PCT);

    date1 = irxDate(2000, 1, 30);
    date2 = irxDate(2000, 6, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 152.0/366.0, TOL_PCT);

    date1 = irxDate(1999, 11, 30);
    date2 = irxDate(2000, 4, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 152.0/366.0, TOL_PCT);

    date1 = irxDate(2001, 3, 1);
    date2 = irxDate(2001, 2, 28);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -1.0/365.0, TOL_PCT);

    date1 = irxDate(2000, 3, 1);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -1.0/366.0, TOL_PCT);

    date1 = irxDate(2001, 3, 1);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -1.0 - 1.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2001, 2, 28);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 366.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2004, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 4.0, TOL_PCT);

    date1 = irxDate(1999, 2, 28);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 365.0/365.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2004, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 4.0 + 1.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2003, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 3.0 + 1.0/366.0, TOL_PCT);

    date1 = irxDate(2000, 2, 28);
    date2 = irxDate(2004, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_ACT_AFB, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 4.0 + 2.0/366.0, TOL_PCT);
}

static void test_irxDayCountFractionAct360() {
    IrxTDate    date1;
    IrxTDate    date2;
    int rc;
    double dcf;

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2001, 1, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 366.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2000, 1, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 0.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 1, 1);
    date2 = irxDate(2000, 1, 2);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2000, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 2, 28);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/360.0, TOL_PCT);

    date1 = irxDate(2001, 2, 28);
    date2 = irxDate(2001, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1.0/360.0, TOL_PCT);

    date1 = irxDate(2003, 11, 1);
    date2 = irxDate(2004, 5, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 182.0/360.0, TOL_PCT);

    date1 = irxDate(1999, 2, 1);
    date2 = irxDate(1999, 7, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 150.0/360.0, TOL_PCT);

    date1 = irxDate(1999, 7, 1);
    date2 = irxDate(2000, 7, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 366.0/360.0, TOL_PCT);

    date1 = irxDate(2002, 8, 15);
    date2 = irxDate(2003, 7, 15);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 334.0/360.0, TOL_PCT);

    date1 = irxDate(2003, 7, 15);
    date2 = irxDate(2004, 1, 15);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 184.0/360.0, TOL_PCT);

    date1 = irxDate(1999, 7, 30);
    date2 = irxDate(2000, 1, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 184.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 1, 30);
    date2 = irxDate(2000, 6, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 152.0/360.0, TOL_PCT);

    date1 = irxDate(1999, 11, 30);
    date2 = irxDate(2000, 4, 30);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 152.0/360.0, TOL_PCT);

    date1 = irxDate(2001, 3, 1);
    date2 = irxDate(2001, 2, 28);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -1.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 3, 1);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -1.0/360.0, TOL_PCT);

    date1 = irxDate(2001, 3, 1);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, -366.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2001, 2, 28);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 365.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2004, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1461.0/360.0, TOL_PCT);

    date1 = irxDate(1999, 2, 28);
    date2 = irxDate(2000, 2, 29);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 366.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2004, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1462.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 2, 29);
    date2 = irxDate(2003, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1096.0/360.0, TOL_PCT);

    date1 = irxDate(2000, 2, 28);
    date2 = irxDate(2004, 3, 1);
    rc = irxDayCountFraction(date1, date2, IRX_ACT_360, &dcf);
    BOOST_CHECK_EQUAL(rc, SUCCESS);
    BOOST_CHECK_CLOSE(dcf, 1463.0/360.0, TOL_PCT);
}

test_suite* init_unit_test_suite(int, char* []) {
    test_suite* tests = BOOST_TEST_SUITE("irx_dateutils tests");

    tests->add(BOOST_TEST_CASE(test_irxDayCountFractionActAct));
    tests->add(BOOST_TEST_CASE(test_irxDayCountFractionActActAFB));
    tests->add(BOOST_TEST_CASE(test_irxDayCountFractionAct360));

    return tests;
}
