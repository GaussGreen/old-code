//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : QuasiContractualBetaSkewMatrix.hpp
//
//   Description : Concrete class of QuasiContractualBaseCorrelation that 
//                 computes a "BetaSkewMatrix" sensitivity after the shift
//                 defined in QuasiContractualBaseCorrelation
//
//   Author      : Jose Hilera
//
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_QUASICONTRACTUALBETASKEWMATRIX_HPP
#define QLIB_QUASICONTRACTUALBETASKEWMATRIX_HPP

#include "edginc/QuasiContractualBaseCorrelation.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL QuasiContractualBetaSkewMatrix : public QuasiContractualBaseCorrelation 
{
public:
     const static CClassConstSP TYPE;
     const static string        NAME;

    /** Returns the sensitivity to be computed in 
     * QuasiContractualBaseCorrelation after applying the shift */
    virtual SensitivitySP sensitivityToComputeAfterShift();
    
private:
    // For Reflection
    QuasiContractualBetaSkewMatrix();
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);  
    
    // Fields
    double shiftSize;

    /**
     * Optional field to determine if we want to do a "clever" tweak (tweak
     * only the relevant points) or if we want to tweak ALL the points (useful
     * as a debug tool to see for example the effect of strike interpolation)
     * */
    bool tweakAll;

    /** Optional field to determine the number of "neighbour" strikes to tweak
        at each side of the used strikes, when tweakAll is false */
    int numberOfNeighbours;
};

DRLIB_END_NAMESPACE

#endif
