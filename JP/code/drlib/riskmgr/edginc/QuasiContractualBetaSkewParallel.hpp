//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : QuasiContractualBetaSkewParallel.hpp
//
//   Description : Concrete class of QuasiContractualBaseCorrelation that 
//                 computes a "BetaSkewParallel" sensitivity after the shift
//                 defined in QuasiContractualBaseCorrelation
//
//   Author      : Jose Hilera
//
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_QUASICONTRACTUALBETASKEWPARALLEL_HPP
#define QLIB_QUASICONTRACTUALBETASKEWPARALLEL_HPP

#include "edginc/QuasiContractualBaseCorrelation.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL QuasiContractualBetaSkewParallel : public QuasiContractualBaseCorrelation 
{
public:
     const static CClassConstSP TYPE;
     const static string        NAME;

    /** Returns the sensitivity to be computed in 
     * QuasiContractualBaseCorrelation after applying the shift */
    virtual SensitivitySP sensitivityToComputeAfterShift();
    
private:
    // For Reflection
    QuasiContractualBetaSkewParallel();
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);  
    
    // Fields
    double shiftSize;
};

DRLIB_END_NAMESPACE

#endif
