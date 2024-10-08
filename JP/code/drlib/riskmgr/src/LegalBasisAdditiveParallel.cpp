/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/GenericScalarOneSidedShift.hpp"
#include "edginc/LegalBasisAdditiveParallelTweak.hpp"

DRLIB_BEGIN_NAMESPACE


/**
 * Calculates derivative of instrument price w.r.t. shifts in the additive
 * coefficients of the Legal Basis adjustment done to the CDSParSpreadCurve
 *
 * The output generated by this sensitivity is:
 * <DL>
 * <DT>LEGAL_BASIS_ADDITIVE_PARALLEL
 * <DD>Change in price when the additive coefficients of the Legal Basis adjustment are shifted up by 1bp
 * </DL>
 *
 * The derivative is reported in basis point units, i.e. the value reported is 
 * 10000 times smaller than the actual partial derivative.  This is implemented 
 * by setting SENSITIVITY_UNIT to 1e-4: see GenericScalarOneSidedShift::SENSITIVITY_UNIT.
 */

typedef GenericScalarOneSidedShift<LegalBasisAdditiveParallelTweak>
    LegalBasisAdditiveParallel;

template<> CClassConstSP const LegalBasisAdditiveParallel::TYPE =
    CClass::registerClassLoadMethod("LegalBasisAdditiveParallel",
                                    typeid(LegalBasisAdditiveParallel),
                                    LegalBasisAdditiveParallel::load);

template<> const string LegalBasisAdditiveParallel::NAME = "LEGAL_BASIS_ADDITIVE_PARALLEL";

const double ONE_BASIS_POINT = 0.0001;
template<> const double LegalBasisAdditiveParallel::DEFAULT_SHIFT = (0.1*ONE_BASIS_POINT);  // Default = 0.1bp
template<> const double LegalBasisAdditiveParallel::SENSITIVITY_UNIT = ONE_BASIS_POINT;


/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the
 Windows exe.
 */
bool LegalBasisAdditiveParallelLinkIn() {
    return LegalBasisAdditiveParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
