//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZC3Stub.hpp
//
//   Description : Helper class for bootstrapping ALIB style zero curve
//
//   Author      : Richard Appleton
//
//   Date        : 2nd May 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_STUB_HPP
#define ZC3_STUB_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include <vector>


DRLIB_BEGIN_NAMESPACE

class DateTime;
class DayCountConvention;
class BadDayConvention;
class Holiday;


class MARKET_DLL ZC3Stub
{
    public:
        static int                 STUB_BASIS;
        static DayCountConvention* STUB_DCC;

        /*
         * Prepare a stub. We will always convert these into ACT/360
         * simple rates, taking account of currency basis.
         */
        ZC3Stub(const DateTime&       startDate,
             const DateTime&          endDate,
            double                    rate,
            const DayCountConvention& dcc,
            double                    basis,
            const DayCountConvention& basisDcc,
            const BadDayConvention&   badDayConv,
            const Holiday&            holidays);

        ZC3Stub(const ZC3Stub& other);
        ~ZC3Stub();
        ZC3Stub& operator=(const ZC3Stub& other);

        string toString() const;

        /*
         * Stubs are sorted by start date.
         */
        bool operator<(const ZC3Stub& other) const;

        const DateTime& getStartDate() const;
        const DateTime& getEndDate() const;
        double          getStubRate() const;

        bool isContiguousWith(const ZC3Stub& next) const;

    private:
        DateTime startDate; // adjusted for bad day convention
        DateTime endDate;   // adjusted for bad day convention
        double   stubRate;
};


typedef vector<ZC3Stub> ZC3StubArray;


DRLIB_END_NAMESPACE
#endif // ZC£_STUB_HPP
