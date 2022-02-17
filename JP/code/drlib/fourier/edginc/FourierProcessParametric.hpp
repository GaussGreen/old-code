/**
 * @file FourierProcessParametric.hpp
 */

#ifndef QLIB_FourierProcessParametric_H
#define QLIB_FourierProcessParametric_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Fourier.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(FourierProcessParametric)

class FOURIER_DLL FourierProcessParametric: public FourierProcess {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    bool allowRiskMapping;

protected:

    FourierProcessParametric(CClassConstSP type);

public:

    ~FourierProcessParametric();

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this process
     *
     * Returns riskMappingAllowed or riskMappingDisallowed depending on
     * allowRiskMapping.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;
};

DRLIB_END_NAMESPACE

#endif
