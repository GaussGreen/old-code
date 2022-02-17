/**
 * @file MCPathConfigParametric.hpp
 */

#ifndef QLIB_MCPathConfigParametric_H
#define QLIB_MCPathConfigParametric_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/IModel.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(MCPathConfigParametric)
FORWARD_DECLARE(DependenceMaker)

/**
 * Base class for MCPathConfig's implementing pdfs controlled by parametric
 * vols: includes a flag for enabling/disabling RiskMapping
 */

class MCARLO_DLL MCPathConfigParametric: public MCPathConfig {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    bool allowRiskMapping;

public:

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingAllowed or riskMappingDisallowed, depending on the
     * allowRiskMapping flag.  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

protected:

    MCPathConfigParametric(CClassConstSP type);
    MCPathConfigParametric(CClassConstSP type, DependenceMakerSP dependenceMaker);

public:

    ~MCPathConfigParametric();
};

DRLIB_END_NAMESPACE

#endif
