//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OutputRequestUtil.cpp
//
//   Description : Utility functions for output requests
//
//   Author      : Andrew J Swain
//
//   Date        : 22 August 2002
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include <algorithm>
#include <set>

DRLIB_BEGIN_NAMESPACE

/** store results for the PAYMENT_DATES output request */
void OutputRequestUtil::recordPaymentDates(Control*             control,
                                           Results*             results,    
                                           const DateTimeArray* paydates)
{
    static const string method = "OutputRequestUtil::recordPaymentDates";
    try {
        OutputRequest* request = 
            control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            if (paydates && paydates->size() == 0) {
                results->storeNotApplicable(request); 
            } 
            else {
                DateTimeListSP datelist(new DateTimeList(paydates));
                results->storeRequestResult(request, datelist); 
            }
        }       
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** store results for the KNOWN_CASHFLOWS output request */
void OutputRequestUtil::recordKnownCashflows(Control*             control,
                                             Results*             results,    
                                             const string&        ccyName,
                                             const CashFlowArray* cfl) 
{
    static const string method = "OutputRequestUtil::recordKnownCashflows";
    try {
        OutputRequest* request = 
            control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        // if we don't add anything here then 'NotApplicable' will be added
        // by 'RiskManager' at the end of the pricing
        if (request && cfl && cfl->size() != 0 ) {
            OutputNameSP name(new OutputName(ccyName));
            IObjectSP    cflows(new CashFlowList(cfl));
            results->storeRequestResult(request, cflows, name); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** store results for the BARRIER_LEVEL output request */
void OutputRequestUtil::recordBarrierLevels(Control*                 control,
                                            Results*                 results,    
                                            const string&            assetName,
                                            const BarrierLevelArray* barriers) 
{
    static const string method = "OutputRequestUtil::recordBarrierLevels";
    try {
        OutputRequest* request = 
            control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        // if we don't add anything here then 'NotApplicable' will be added
        // by 'RiskManager' at the end of the pricing
        if (request && barriers && !barriers->empty()) {
            OutputNameSP name(new OutputName(assetName));
            IObjectSP    barrs(copy(barriers));
            results->storeRequestResult(request, barrs, name); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** store ExpiryResultArrays, with the name passed as parameter */
void OutputRequestUtil::recordExpiryResultArray(
     Control*                 control,
     Results*                 results,    
     const ExpiryResultArray* erArray,
     const string&            outputName)
{
    static const string method = "OutputRequestUtil::recordExpiryResultArray";
    try {
        OutputRequest* request = control->requestsOutput(outputName);

        // if we don't add anything here then 'NotApplicable' will be added
        // by 'RiskManager' at the end of the pricing
        if (request && erArray && (erArray->size() != 0)) {
            IObjectSP    expResArray(copy(erArray));
            results->storeRequestResult(request, expResArray); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** store a named ExpiryResultArrays, with the name passed as parameter */
void OutputRequestUtil::recordNamedExpiryResultArray(
     Control*                 control,
     Results*                 results,    
     const ExpiryResultArray* erArray,
     const string&            outputName,
     const string&            name)
{
    static const string method = "OutputRequestUtil::recordNamedExpiryResultArray";
    try {
        OutputRequest* request = control->requestsOutput(outputName);

        // if we don't add anything here then 'NotApplicable' will be added
        // by 'RiskManager' at the end of the pricing
        if (request && erArray && (erArray->size() != 0)) {
            IObjectSP    expResArray(copy(erArray));
            OutputNameSP outName(new OutputName(name));
            results->storeRequestResult(request, expResArray, outName); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Ideally this struct lives in OutputRequestUtil::KnownCashFlows::Imp but
    assembler (solaris opt) can't cope with such long names */
/** Holds information per ccy (supports XCB's etc). Note that we 
    report in different curriencies because struck divs
    are paid using the fx rate on the day */
struct ORUCashFlowsForCcy{
    CashFlowArray   cfl;       // what get's paid when 
    set<DateTime>   estimates; // dates for which the amount is estimated
};

class OutputRequestUtil::KnownCashFlows::Imp {
public:
    typedef map<string, ORUCashFlowsForCcy> CashFlowsPerCcy;
    CashFlowsPerCcy  cflPerCcy; // payment info per ccy
};

OutputRequestUtil::KnownCashFlows::KnownCashFlows(): my(new Imp())
{}

/** Record a known cashFlow on supplied date */
void OutputRequestUtil::KnownCashFlows::addKnownCashFlow(
    const string&   ccyISOCode,
    const DateTime& date,
    double          amount){
    ORUCashFlowsForCcy& cfsPerCcy = my->cflPerCcy[ccyISOCode];
    cfsPerCcy.cfl.push_back(CashFlow(date, amount));
}

/** Record a known cashFlow */
void OutputRequestUtil::KnownCashFlows::addKnownCashFlow(
    const string&   ccyISOCode,
    const CashFlow& cf)
{
    ORUCashFlowsForCcy& cfsPerCcy = my->cflPerCcy[ccyISOCode];
    cfsPerCcy.cfl.push_back(cf);
}

/** Record a unknown cashFlow on the specified date*/
void OutputRequestUtil::KnownCashFlows::addUnknownCashFlowDate(
    const string&   ccyISOCode,
    const DateTime& date)
{
    ORUCashFlowsForCcy& cfsPerCcy = my->cflPerCcy[ccyISOCode];
    cfsPerCcy.estimates.insert(date);
}

/** Write known cashFlows to Results object */
void OutputRequestUtil::KnownCashFlows::recordKnownCashFlows(Control*   control,
                                                             Results*   results)
{
    /* loop through the different currencies */
    for (Imp::CashFlowsPerCcy::iterator iter = my->cflPerCcy.begin();
         iter != my->cflPerCcy.end(); ++iter){
        // for ease
        ORUCashFlowsForCcy& cashFlowsForCcy = iter->second;
        // set up array of known cash Flows
        CashFlowArray knownCFL;
        // then loop through all cash Flows
        for (int i = 0; i < cashFlowsForCcy.cfl.size(); i++){
            if (cashFlowsForCcy.estimates.find(cashFlowsForCcy.cfl[i].date)
                == cashFlowsForCcy.estimates.end()){
                // not an estimate => must be known
                knownCFL.push_back(cashFlowsForCcy.cfl[i]);
            }
        }
        // sort the cash Flow array
        sort(knownCFL.begin(),knownCFL.end(), CashFlow::lessThenForDates);
        // then aggregate
        CashFlow::aggregate(knownCFL);
        // then store
        recordKnownCashflows(control, results, iter->first, &knownCFL);
    }
}

OutputRequestUtil::KnownCashFlows::~KnownCashFlows(){}

DRLIB_END_NAMESPACE
