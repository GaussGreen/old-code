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

#include "edginc/config.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/BetaSkewMatrixTweak.hpp"
#include "edginc/QuasiContractualBetaSkewMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

const string QuasiContractualBetaSkewMatrix::NAME = 
    "QUASI_CONTRACTUAL_BETA_SKEW_MATRIX";


QuasiContractualBetaSkewMatrix::QuasiContractualBetaSkewMatrix() : 
    QuasiContractualBaseCorrelation(TYPE, NAME), 
    shiftSize(BetaSkewMatrixTweak::DEFAULT_SHIFT),
    tweakAll(DEFAULT_VALUE_FOR_TWEAKALL),
    numberOfNeighbours(DEFAULT_VALUE_FOR_NUM_NEIGHBOURS)
{}


/** Returns the sensitivity to be computed in 
 * QuasiContractualBetaSkewMatrix after applying the shift */
SensitivitySP QuasiContractualBetaSkewMatrix::sensitivityToComputeAfterShift() {
    return SensitivitySP(new BetaSkewMatrixTweak(shiftSize, 
                                                 tweakAll, 
                                                 numberOfNeighbours));
}


void QuasiContractualBetaSkewMatrix::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(QuasiContractualBetaSkewMatrix, clazz);
    SUPERCLASS(QuasiContractualBaseCorrelation);
    EMPTY_SHELL_METHOD(defaultConstructor);
    
    FIELD(shiftSize, "Shift size");
    FIELD_MAKE_OPTIONAL(shiftSize);

    FIELD(tweakAll,
                 "YES = tweaks ALL points, NO = tweaks only relevant points");
    FIELD_MAKE_OPTIONAL(tweakAll);
    FIELD(numberOfNeighbours,
                 "number of neighbour strikes to tweak at each side of "
                 "the used strikes, when tweakAll is false");
    FIELD_MAKE_OPTIONAL(numberOfNeighbours);
}

IObject* QuasiContractualBetaSkewMatrix::defaultConstructor() {
    return new QuasiContractualBetaSkewMatrix();
}

CClassConstSP const QuasiContractualBetaSkewMatrix::TYPE =
CClass::registerClassLoadMethod("QuasiContractualBetaSkewMatrix",
                                typeid(QuasiContractualBetaSkewMatrix),
                                QuasiContractualBetaSkewMatrix::load);


bool QuasiContractualBetaSkewMatrixLinkIn() {
    return QuasiContractualBetaSkewMatrix::TYPE != NULL;
}

DRLIB_END_NAMESPACE
