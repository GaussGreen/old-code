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

#include "edginc/config.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/BetaSkewParallel.hpp"
#include "edginc/QuasiContractualBetaSkewParallel.hpp"

DRLIB_BEGIN_NAMESPACE

const string QuasiContractualBetaSkewParallel::NAME = 
    "QUASI_CONTRACTUAL_BETA_SKEW_PARALLEL";


QuasiContractualBetaSkewParallel::QuasiContractualBetaSkewParallel() : 
    QuasiContractualBaseCorrelation(TYPE, NAME), 
    shiftSize(BetaSkewParallel::DEFAULT_SHIFT)
{}


/** Returns the sensitivity to be computed in 
 * QuasiContractualBetaSkewParallel after applying the shift */
SensitivitySP QuasiContractualBetaSkewParallel::sensitivityToComputeAfterShift() {
    return SensitivitySP(new BetaSkewParallel(shiftSize));
}


void QuasiContractualBetaSkewParallel::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(QuasiContractualBetaSkewParallel, clazz);
    SUPERCLASS(QuasiContractualBaseCorrelation);
    EMPTY_SHELL_METHOD(defaultConstructor);
    
    FIELD(shiftSize, "Shift size");
    FIELD_MAKE_OPTIONAL(shiftSize);
}

IObject* QuasiContractualBetaSkewParallel::defaultConstructor() {
    return new QuasiContractualBetaSkewParallel();
}

CClassConstSP const QuasiContractualBetaSkewParallel::TYPE =
CClass::registerClassLoadMethod("QuasiContractualBetaSkewParallel",
                                typeid(QuasiContractualBetaSkewParallel),
                                QuasiContractualBetaSkewParallel::load);


bool QuasiContractualBetaSkewParallelLinkIn() {
    return QuasiContractualBetaSkewParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
