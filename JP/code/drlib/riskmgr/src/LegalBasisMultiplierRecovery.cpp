/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/GenericScalarOneSidedShift.hpp"
#include "edginc/LegalBasisMultiplierRecoveryTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Calculates derivative of instrument price w.r.t. shifts in the multiplier
 * recovery coefficient of the Legal Basis adjustment done to the CDSParSpreadCurve
 *
 * The output generated by this sensitivity is:
 * <DL>
 * <DT>LEGAL_BASIS_MULTIPLIER_RECOVERY
 * <DD>Change in price when the multiplier recovery coefficient of the Legal Basis adjustment is shifted up by 1%
 * </DL>
 *
 * The derivative is reported in 1% units, i.e. the value reported is 
 * 100 times smaller than the actual partial derivative.  This is implemented 
 * by setting SENSITIVITY_UNIT to 0.01: see GenericScalarOneSidedShift::SENSITIVITY_UNIT.
 */

typedef GenericScalarOneSidedShift<LegalBasisMultiplierRecoveryTweak>
    LegalBasisMultiplierRecovery;

template<> CClassConstSP const LegalBasisMultiplierRecovery::TYPE =
    CClass::registerClassLoadMethod("LegalBasisMultiplierRecovery",
                                    typeid(LegalBasisMultiplierRecovery),
                                    LegalBasisMultiplierRecovery::load);

template<> const string LegalBasisMultiplierRecovery::NAME =
    "LEGAL_BASIS_MULTIPLIER_RECOVERY";

template<> const double LegalBasisMultiplierRecovery::DEFAULT_SHIFT = 0.01;    // Requirement: 1%
template<> const double LegalBasisMultiplierRecovery::SENSITIVITY_UNIT = 0.01; // Requirement: 1%


/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe.
 */
bool LegalBasisMultiplierRecoveryLinkIn() {
    return LegalBasisMultiplierRecovery::TYPE != NULL;
}

DRLIB_END_NAMESPACE
