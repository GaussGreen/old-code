//----------------------------------------------------------------------------
//
//   Group       : QR Credit Hybrids
//
//   Description : Container class for Credit Metrics name-specific model params
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CM_ONLY_PARAMETERS_HPP
#define QLIB_CM_ONLY_PARAMETERS_HPP

#include "edginc/TweakableWith.hpp"
#include "edginc/MappingFunction.hpp"
#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/CCMBetaSens.hpp"
#include "edginc/CCMAbsoluteBetaTweak.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CmOnlyParameters : 
    public RationalisedCreditEngineParameters,
    virtual public CCMBetaSens::IShift,
    virtual public TweakableWith<CCMAbsoluteBetaTweak> 
{
public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~CmOnlyParameters();
    
    /** Public constructor */
    CmOnlyParameters(string name, double beta, CDoubleSP decretionBetaSP);

    /** basic validation */
    virtual void validatePop2Object();

    /** Returns the name of this object */
    virtual string getName() const;

    /** Sets beta */
    virtual void setBeta(const double& newBeta);

    /** Returns the beta  of this object */
    virtual const double getBeta() const;   

    /** Returns the decretion beta of this object */
    virtual const double getDecretionBeta() const;   

    /** CCMBetaSens::IShift implementation */
    virtual string sensName(CCMBetaSens* shift) const;

    /** CCMBetaSens::IShift implementation */
    virtual bool sensShift(CCMBetaSens* shift);

    /** TweakableWith<CCMAbsoluteBetaTweak> implementation */
    virtual string sensName(CCMAbsoluteBetaTweak* shift) const;

    /** TweakableWith<CCMAbsoluteBetaTweak> implementation */
    virtual bool sensShift(CCMAbsoluteBetaTweak* shift);

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** For reflection */
    CmOnlyParameters();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    double checkBeta(double newBeta) const;

    // ----------------
    // FIELDS
    // ----------------  
    /** Name of this market object */
    string name;

    double beta;

    double decretionBeta;
};


typedef smartPtr<CmOnlyParameters> CmOnlyParametersSP;
typedef smartConstPtr<CmOnlyParameters> CmOnlyParametersConstSP;

// Support for wrapper
typedef MarketWrapper<CmOnlyParameters> CmOnlyParametersWrapper;

#ifndef QLIB_CMNLYPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CmOnlyParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CmOnlyParameters>);
#endif

DRLIB_END_NAMESPACE

#endif //QLIB_CM_ONLY_PARAMETERS_HPP
