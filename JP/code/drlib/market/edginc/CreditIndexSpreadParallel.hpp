/**
 * @file CreditIndexSpreadParallel.hpp
 */

#ifndef DRLIB_CreditIndexSpreadParallel_H
#define DRLIB_CreditIndexSpreadParallel_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

struct MARKET_DLL CreditIndexSpreadParallel: CObject {

    static CClassConstSP const TYPE;
    CreditIndexSpreadParallel(); ~CreditIndexSpreadParallel();

    typedef Void Qualifier;

    enum { discrete = 0 };
};

//declaration of template specialisation
template <>
RiskMappingMatrixConstSP RiskProperty<CreditIndexSpreadParallel>::riskMappingMatrix(
    IObjectConstSP world,
    OutputNameConstSP name) const;

DRLIB_END_NAMESPACE

#endif

