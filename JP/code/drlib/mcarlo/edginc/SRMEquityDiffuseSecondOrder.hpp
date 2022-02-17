//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEquityDiffuseSecondOrder.hpp
//
//   Description : A generator of paths using stochastic rates
//                 for Equity Assets and the 'second order method'
//                 for the SDE approximation
//
//   Date        : Aug 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMEQUITYDIFFUSESECONDORDER_HPP
#define EDR_SRMEQUITYDIFFUSESECONDORDER_HPP

DRLIB_BEGIN_NAMESPACE  

class SRMEquityDiffuseSecondOrder : public SRMEquityDiffuse
{
public:

    SRMEquityDiffuseSecondOrder(IQMCDiffusibleInterestRateSP domIR) : SRMEquityDiffuse(domIR) {}
    virtual ~SRMEquityDiffuseSecondOrder() {}

    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
    virtual void finalize(DateTimeArrayConstSP allDates);

};

DRLIB_END_NAMESPACE  

#endif
