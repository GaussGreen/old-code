//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OutputRequestUtil.hpp
//
//   Description : Utility functions for output requests
//
//   Author      : Andrew J Swain
//
//   Date        : 22 August 2002
//
//
//----------------------------------------------------------------------------

#ifndef _OUTPUTREQUESTUTIL_HPP
#define _OUTPUTREQUESTUTIL_HPP

#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/BarrierLevel.hpp"
#include "edginc/ExpiryResult.hpp"

DRLIB_BEGIN_NAMESPACE


class RISKMGR_DLL OutputRequestUtil {
public:

    /** store results for the PAYMENT_DATES output request */
    static void recordPaymentDates(
        Control*             control,
        Results*             results,    
        const DateTimeArray* paydates);

    /** store results for the KNOWN_CASHFLOWS output request */
    static void recordKnownCashflows(
        Control*             control,
        Results*             results,    
        const string&        ccyName,
        const CashFlowArray* cfl);

    /** store results for the BARRIER_LEVEL output request */
    static void recordBarrierLevels(
        Control*                 control,
        Results*                 results,    
        const string&            assetName,
        const BarrierLevelArray* barriers);

    /** store ExpiryResultArrays, with the name passed as parameter */
    static void recordExpiryResultArray(
        Control*                 control,
        Results*                 results,    
        const ExpiryResultArray* erArray,
        const string&            outputName);

    /** store a named ExpiryResultArrays, with the name passed as parameter */
    static void recordNamedExpiryResultArray(
        Control*                 control,
        Results*                 results,    
        const ExpiryResultArray* erArray,
        const string&            outputName,
        const string&            name);

    /** Utility for collecting cashflows from different sources and ensuring
        that known cash flows really are known (eg if a cash flow on date A
        is known from one source but is only an estimate from somewhere else
        then it is not known - think Equity Swap with 4 legs) */
    class RISKMGR_DLL KnownCashFlows{
    public:
        //// Constructor - requires no arguments
        KnownCashFlows();
        /** Record a known cashflow on supplied date */
        void addKnownCashFlow(const string&   ccyISOCode,
                              const DateTime& date,
                              double          amount);
        /** Record a known cashflow */
        void addKnownCashFlow(const string&   ccyISOCode,
                              const CashFlow& cf);

        /** Record a unknown cashflow on the specified date*/
        void addUnknownCashFlowDate(const string&   ccyISOCode,
                                    const DateTime& date);

        /** Write known cashflows to Results object */
        void recordKnownCashFlows(Control*   control,
                                  Results*   results);
        ~KnownCashFlows();
    private:
        KnownCashFlows(const KnownCashFlows &rhs); // not implemented
        KnownCashFlows& operator=(const KnownCashFlows& rhs); //not implemented

        class Imp;
        auto_ptr<Imp>  my;
    };
};

DRLIB_END_NAMESPACE

#endif
