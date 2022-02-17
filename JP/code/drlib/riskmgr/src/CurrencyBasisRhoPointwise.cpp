/**
 *
 */


#include "edginc/config.hpp"
#include "edginc/GenericVectorOneSidedShift.hpp"
#include "edginc/CurrencyBasisRhoPointwiseTweak.hpp"

DRLIB_BEGIN_NAMESPACE

typedef GenericVectorOneSidedShift<
    CurrencyBasisRhoPointwiseTweak> CurrencyBasisRhoPointwise;

template<> CClassConstSP const CurrencyBasisRhoPointwise::TYPE =
    CClass::registerClassLoadMethod("CurrencyBasisRhoPointwise",
                                    typeid(CurrencyBasisRhoPointwise),
                                    load);

template<> const string CurrencyBasisRhoPointwise::NAME = "CCY_BASIS_RHO_POINTWISE";

const double ONE_BASIS_POINT = 0.0001;
template<> const double CurrencyBasisRhoPointwise::SENSITIVITY_UNIT = ONE_BASIS_POINT;
template<> const double CurrencyBasisRhoPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;

/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe.
 */
bool CurrencyBasisRhoPointwiseLinkIn() {
    return CurrencyBasisRhoPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE
