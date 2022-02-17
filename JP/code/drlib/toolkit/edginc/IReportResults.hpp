//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RiskMgrInterface.hpp
//
//   Description : Defines interface to RiskMgr
//
//   Author      : Andrew J Swain
//
//   Date        : 19 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_IREPORTRESULTS_HPP
#define EDR_IREPORTRESULTS_HPP

#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

// FIXME: toy class to check integration of QSRM with SimDesk
class IReportResultsLite {
    public:
        virtual double * getRes() {return NULL;}
        virtual size_t   getLen() {return 0;}
        virtual ~IReportResultsLite() {}
};

/// Auxilary class for products that want to return results in reg.test framework
class IReportResults : public IReportResultsLite {
    public:
        virtual IObjectSP getResults() = 0;
};

DRLIB_END_NAMESPACE

#endif
