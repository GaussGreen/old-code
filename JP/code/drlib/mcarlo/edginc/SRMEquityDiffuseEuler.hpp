//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEquityDiffuseEuler.hpp
//
//   Description : A generator of paths using stochastic rates
//                 for Equity Assets and the Euler method for the 
//                 SDE approximation
//
//   Date        : Aug 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMEQUITYDIFFUSEEULER_HPP
#define EDR_SRMEQUITYDIFFUSEEULER_HPP

DRLIB_BEGIN_NAMESPACE  

class SRMEquityDiffuseEuler : public SRMEquityDiffuse
{
public:

    SRMEquityDiffuseEuler(IQMCDiffusibleInterestRateSP domIR) : SRMEquityDiffuse(domIR) {}
    virtual ~SRMEquityDiffuseEuler() {}

    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
};

DRLIB_END_NAMESPACE  

#endif
