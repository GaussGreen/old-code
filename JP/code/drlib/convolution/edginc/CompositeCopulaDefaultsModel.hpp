//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 27-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_COMPOSITECOPULADEFAULTSMODEL_HPP
#define QLIB_COMPOSITECOPULADEFAULTSMODEL_HPP

#include "edginc/IConditionalDefaultsModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * "Composite Copula" conditional defaults model ("CCM"):
 * 
 * CCM = Product of 3 copulae:
 * - C1: "dependence" copula where all names default together with a "catastrophic"
 *       recovery R_c
 * - C2: any "mixture" copula (eg: Credit Metrics or RFL), recovery used is a
 *       "non catastrophic" R_nc (internally calibrated)
 * - C3: "independence" copula where each name defaults independently with a
 *       "non catastrophic" R_nc (internally calibrated)
 * 
 * For a given sigle-name survival probability SP_i(T), we write:
 * SP_i(T) = SP1_i(T).SP2_i(T).SP3_i(T)
 * 
 * where
 *   SP1_i(T) is explicitely defined by a "floor curve" (par spreads)
 *   SP2_i(T) = (SP_i(T)/SP1_i(T))^(1-a)
 *   SP2_i(T) = (SP_i(T)/SP1_i(T))^a
 *   a        = "independence factor" (double in [0,1])
 * 
 * For a name i, define:
 * X1_i = Z1
 * X2_i = f(Z2, e2_i)
 * X3_i = e3_i
 * 
 * where
 *   Z          = (Z1, Z2) = 2D market factor
 *   e2_i, e3_i = iid N(0,1) variables independent of Z
 *   f          = "mixture" copula (eg: Credit Metrics or RFL)
 * 
 * Defaults happens when X1_i < c1_i or X2_i < c2_i or X3_i < c3_i
 * 
 * We use the following rule to define which recovery (R_c or R_nc) to use in
 * case of "multiple defaults" (eg: X1_i < c1_i AND X2_i < c2_i):
 * "If X1_i < c1_i, use R_c. Otherwise, use R_nc."
 * 
 * */
class CONVOLUTION_DLL CompositeCopulaDefaultsModel:
    public CObject,
    public virtual IConditionalDefaultsModel
{
public:

    /** Supported types for "mixture copula" C2 (see field 'mixtureCopulaType') */
    enum CopulaTypes
    {
        CREDIT_METRICS,
        RFL
    };

	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~CompositeCopulaDefaultsModel();

    /** [Implements IConditionalDefaultsModel] */
    virtual ICondLossDistributionsGenKeySP initialise(
        double defaultProba,
        double expectedLoss,
        double notional,
        CreditEngineParametersConstSP modelParameters,
        const DateTime& startDate,
        const DateTime& endDate) const;

    /** [Implements IConditionalDefaultsModel] */
    virtual double integrateCondFunction(
        const MFunctionND* condFunction,
        ICondLossDistributionsGenKeyArrayConstSP condKeys,
        const DateTime& time) const; 

    /** [Implements IConditionalDefaultsModel] */
    virtual bool isCondFunctionIntegrationTimeDependent() const;

    /** [Implements IConditionalDefaultsModel] */
    virtual int marketFactorDimension() const;

    /** [Implements IConditionalDefaultsModel] */
    virtual CClassConstSP engineParamsType() const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

private:
    CompositeCopulaDefaultsModel();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
    /** Market factor integrator - used for internal thresholds calibration */
    IIntegratorSP integrator;

    /** Market factor model - defines distribution of the market factor */
    IMarketFactorModelSP marketFactorModel;
    
    /** Define the type of the "mixture copula" C2 */
    // NB: this field could be of type IConditionalDefaultsModelSP, but then
    //     we would rely on the user to populate the IConditionalDefaultsModel with
    //     consistent marketFactorModel and integrator. Seems more natural to
    //     internally build a consistent IConditionalDefaultsModelSP mixtureCopula
    //     with the CCM marketFactorModel and integrator.
    CopulaTypes mixtureCopulaType;

    /** (Internally built) "mixture copula" model */
    IConditionalDefaultsModelSP mixtureCopula;
};

DECLARE(CompositeCopulaDefaultsModel);

DRLIB_END_NAMESPACE

#endif /*QLIB_COMPOSITECOPULADEFAULTSMODEL_HPP*/
