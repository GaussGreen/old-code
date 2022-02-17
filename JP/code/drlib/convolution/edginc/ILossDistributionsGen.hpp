//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 03-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ILOSSDISTRIBUTIONSGEN_HPP
#define QLIB_ILOSSDISTRIBUTIONSGEN_HPP

#include "edginc/ICreditLossGen.hpp"
#include "edginc/IDistribution1D.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Specialisation of ICreditLossGen responsible for creating (unconditional)
 * loss distributions.
 * */
class CONVOLUTION_DLL ILossDistributionsGen: public virtual ICreditLossGen {
public:

    virtual ~ILossDistributionsGen();

    // Remarks:
    // 1. We deliberately don't pass counterparty or recover notional
    //    information (methodology for "inner loss configs" is really not clear)
    // 2. No "control" and "results" inputs -> needed ?
    // 3. The output distributions array is a term structure of distributions
    //    indexed by the timeline
    virtual IDistribution1DArraySP createLossDistributions(
        const DateTimeArray& timeline) const = 0;
        
protected:
    ILossDistributionsGen();

private:    
    ILossDistributionsGen(const ILossDistributionsGen& rhs); // don't use
    ILossDistributionsGen& operator=(const ILossDistributionsGen& rhs); // don't use
};

DECLARE_REF_COUNT(ILossDistributionsGen);

DRLIB_END_NAMESPACE

#endif /*QLIB_ILOSSDISTRIBUTIONSGEN_HPP*/
