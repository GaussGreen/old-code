//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : Type of Credit Asset characterised by a single CDSParSpreads
//                 curve and some engine specific parameters.
//
//                 NOTE: This class is what used to be SingleCreditAsset, since:
//                  a) SingleCreditAsset was the only derived type,  
//                  b) Typically CreditAssets were dynamically cast into 
//                     SingleCreditAssets anyway, and 
//                  c) The CreditAsset abstraction does not make much sense.
//                 For backwards compatibility reasons, SingleCreditAsset still
//                 exists as a class and derives from this CreditAsset, but 
//                 provides no additional functionality.
//
//   Date        : March 2005
//
//----------------------------------------------------------------------------

#ifndef EDR_CREDIT_ASSET_HPP
#define EDR_CREDIT_ASSET_HPP

#include "edginc/CreditAsset.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/CreditEngineParameters.hpp"

DRLIB_BEGIN_NAMESPACE

class CreditAsset;
#ifndef QLIB_CREDITASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CreditAsset>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CreditAsset>);
#endif

/** A CreditAsset is capable of producing an expected loss along a timeline.
    For example, a CDO takes an array of these as its underlyings */
class MARKET_DLL CreditAsset: public MarketObject {
public:
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~CreditAsset();
    
    /** Pull out the CDSParSpreads and chain down to the engine parameters
        object */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the name of this object. This is the same as the
        CDSParSpreads object */
    virtual string getName() const;
    
    /** Returns the default date */
    virtual DateTime defaultDate() const;

    virtual double recoveryRate() const;

    /** If required, converts the engineParams from an old-style parameter
        into a new-style paramter. See CreditEngineParameters for details */
    void convertToNewParamStyle(const string nameOfPortfolio,
                                CDoubleSP beta, 
                                CDoubleSP decretionBeta);
    
    bool hasDefaulted() const;

    /** Returns the underlying spread curve */
    ICDSParSpreadsConstSP getParSpreadCurve() const;

    /** Returns the engine parameters for this asset (NB Model needs to
        have selected appropriate set during getMarket()). The parameter 
        can specify what type of engine parameters are required. An 
        exception is thrown if the relevant type is not found */
    CreditEngineParametersConstSP getEngineParams(
        CClassConstSP engineParamsType) const;

    CreditEngineParametersConstSP getEngineParams() const 
      { return engineParams->getEngineParams(); }

    /** set the engine parameters 
      * This function modifies the object and is used in the calibration
      * to override the engineparams by an initial guess 
      */
    void setEngineParams(CreditEngineParametersSP engineParams);


protected:
    CreditAsset(const CClassConstSP& clazz);

private:
    friend class StopCompilerWarning;
    CreditAsset();
    CreditAsset(const CreditAsset& rhs); // don't use
    CreditAsset& operator=(const CreditAsset& rhs); // don't use
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
// TODO should be private
public:
    /// fields ///
    ICDSParSpreadsWrapper         cdsParSpreads;///<! CDS par curve and recovery
    CreditEngineParametersWrapper engineParams;
};

typedef smartPtr<CreditAsset> CreditAssetSP;
typedef smartConstPtr<CreditAsset> CreditAssetConstSP;
//// support for wrapper class
typedef MarketWrapper<CreditAsset> CreditAssetWrapper;
typedef smartPtr<CreditAssetWrapper> CreditAssetWrapperSP;
// then define an array of CreditAssetWrapper. Note here is an array of
// smart pointers to market object wrappers (rather than an array of structures)
typedef array<CreditAssetWrapperSP,
              CreditAssetWrapper> CreditAssetWrapperArray;

DRLIB_END_NAMESPACE

#endif
