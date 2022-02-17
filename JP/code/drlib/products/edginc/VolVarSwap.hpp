//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolVarSwap.hpp
//
//   Description : Utility method for variance swaps
//
//   Author      : Andrew J Swain
//
//   Date        : 26 February 2004
//
//
//----------------------------------------------------------------------------

#ifndef VOLVARSWAP_HPP
#define VOLVARSWAP_HPP

#include "edginc/DateTime.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Asset.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Model.hpp"
#include "edginc/ModelLN.hpp"

#include "edginc/Generic1Factor.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/SensitiveStrikes.hpp"

#include "edginc/NumericalIntegrationLN.hpp"
#include "edginc/ImpliedIntegration.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/CreditSupport.hpp"

#include "edginc/AssetHistory.hpp"
#include "edginc/VarSwapBasis.hpp"
#include "edginc/EquityBase.hpp"

#include "edginc/VarSwapUtilities.hpp"
#include "edginc/VegaMatrixLite.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL ClosedFormIntegrateLN: public CModelLN{
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object();

    /** constructor */
    ClosedFormIntegrateLN();
    
    ClosedFormIntegrateLN(const string& volType); /* takes type of vol to use */

    ClosedFormIntegrateLN(double        stdDevInput,
                          double        stdDevNbUp,
                          double        stdDevNbDown,
                          int           nbStrikesPerLogStrike,
                          double        nbStrikesVega,
                          const string& integrationMethod,
                          const string& strikesForVegaMatrix,
                          double        absVolPrecision,
                          double        relVolPrecision,
                          bool          useBasis);

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        /** invoke the pricing */
        virtual void price(ClosedFormIntegrateLN* model,
                           Control*                control, 
                           CResults*               results) = 0;

        virtual ~IProduct() {};
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class ClosedFormIntegrateLNHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormIntegrateLN* model) const = 0;
    };

    /** calculate single price and store result in results */
    void Price(CInstrument*  instrument, 
               CControl*     control, 
               CResults*     results);
    
    /** computes input parameters for integration, ie lowStrike and highStrike and nb of steps */
    const void limits (const Asset*       asset,
                 const DateTime&    today,
                 const DateTime&    maturity,
                 double&            lowStrike,
                 double&            highStrike,
                 int&               nbSteps) const;    
    
    inline string integrationMethodGet(void) const {return integrationMethod;}

    double calcAbsPrecision(double& tenor) const;
    double calcRelPrecision(double& tenor) const;

    inline double minToleranceGet(void) const {return MIN_TOLERANCE;}

    /** for VEGA_MATRIX - chose some semi-arbitrary strikes to report against */
    DoubleArraySP sensitiveStrikes(OutputNameConstSP outputName,
                                   const Asset*     asset,
                                   const DateTime&  today,
                                   const DateTime&  maturity);

    void getMarket(const MarketData* market,
                   IInstrumentCollectionSP instruments);
    
    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    /** Tycho's original implementation */
    static const string DISCRETE_INTEGRATION;
    
    /** Returns default method by string */
    static string getDefaultIntegrationMethod();

    /** possibilities for vegamatrix */
    static const string VOLSURFACE_STRIKES;
    static const string EQUIDISTANT_STRIKES;
    static const string EQUIDISTANT_VEGA;

    //divMethodology used to compute prices with dividends
    static const string DIV_CONTINUOUS_APPROX;
    static const string DIV_BRUTE_FORCE;
    static const string DIV_DEFAULT;
 
    /** Returns default strikes for vegamatrix by string */
    static string getDefaultStrikesForVegaMatrix();

    bool getUseBasis() const;

    //returns divMethodology used to compute prices with dividends
    string getDivMethodology() const;

    /** Returns a series of dates and weights for the term structure of portfolios */
    
    //using log-forward coefficients (Manos initial implementation)
    void getPriceTSPortfolios(const CAsset*   asset,
                              const DateTime&       valueDate,
                              const DateTimeArray&  obsDates,
                              DateTimeArray&        datesTS,
                              DoubleArray&          weightsTS) const;

    //using log-pv coeffieints (Gad implementation to account for discrete dividends)
    void getPriceTSPortfolios(const YieldCurve*     discount,
                              const DateTime&       valueDate,
                              const DateTimeArray&  obsDates,
                              DateTimeArray&        datesTS,
                              DoubleArray&          weightsTS) const;
    
    /** Produces an integrator from model parameters and modifies integrationDomain
        depending on control. Original and resulting integrationDomain must be in
        relative strikes (fwd moneyness) */
    void getIntegrator(const Control*  control,
                       const CAsset*   asset,
                       const DateTime& today,
                       const DateTime& maturity,
                       VanillaContractsRecorderSP recorder,
                       Integrator1DSP& integrator,
                       Range&          integrationDomain) const;
    
private:
    /** "Default" */
    static const string DEFAULT_INTEGRATION_METHOD;
    static const string DEFAULT_VEGAMATRIX_STRIKES;

    static const double MIN_TOLERANCE;
    static const double MAX_TOLERANCE;
    static const double DEFAULT_ABSOLUTE_VOL_PRECISION;
    static const double DEFAULT_RELATIVE_VOL_PRECISION;

    static const string TERM_STRUCTURE_DEFAULT;
    static const string TERM_STRUCTURE_NONE;
    static const string TERM_STRUCTURE_ALL;
    static const string TERM_STRUCTURE_DEFAULT_TENOR;
    
    friend class ClosedFormIntegrateLNHelper;

    double  stdDevInput;            // stdDev used in order to compute lowStrike and highStrike 
    double  stdDevNbUp;             // nb of stdDevs to compute highStrike
    double  stdDevNbDown;           // nb of stdDev to compute lowStrike
    int     nbStrikesPerLogStrike;  // nb of steps to compute approximation
    double  nbStrikesVega;          // nb of strike for VegaMatrix
    string  integrationMethod;      // Tycho or IMSL Integrator (IntFuncInf, IntFuncGeneral, ClosedRomberg)
    string  strikesForVegaMatrix;   // volSurfaceStrikes, equidistantStrikes, equidistantVega
    double  absVolPrecision;        // abs precision in terms of vol, to be converted to abs precision in terms of integrator
    double  relVolPrecision;        // rel precision in terms of vol, to be converted to rel precision in terms of integrator

    OutputNameArrayConstSP names;
    DoubleArrayArray volSurfaceStrikes;

    bool    useBasis;           // Defaulted to true and must be overriden for XCB, Fund, CCY Prot and CCY Struck
    string  termStructMethod;   // NONE, ALL, BENCHMARKS

    string  divMethodology;  //divMethodology used to compute prices with dividends
    
    ClosedFormIntegrateLN(const ClosedFormIntegrateLN &rhs);
    ClosedFormIntegrateLN& operator=(const ClosedFormIntegrateLN& rhs);
};

typedef smartPtr<ClosedFormIntegrateLN> ClosedFormIntegrateLNSP;

///////////////////////////////////////////////////////////////////////////


/** Utility method for variance swaps */
class PRODUCTS_DLL VarianceSwapUtil {
public:
    static double futureVar(const CAsset*           asset,
                            const YieldCurve*       discount,
                            const DateTime&         valueDate,
                            const DateTime&         volDate,
                            const IModel*           model,
                            Control*                control);
    
    static double priceVarSwapSimple(double totalVar, 
        double strikeVol, 
        double notional, 
        double scaleFactor,
        bool dontScaleByStrike);
    
    static double futureFwdVar(
        const CAsset*        asset,
        const YieldCurve*    discount,
        const DateTime&      valueDate,
        const DateTime&      startDate,
        const DateTime&      maturity,
        const IModel*        model,
        Control*             control);
    
    static double priceFwdStartingVarSwap(const CAsset*        asset,
        const YieldCurve*    discount,
        double               notional,
        double               scaleFactor,
        bool                 dontScaleByStrike,
        bool                 noDivAdj,
        const DateTime&      valueDate,
        const DateTime&      startDate,
        const DateTime&      maturity,
        double               strikeVol,
        int                  observationsPerYear,
        int                  observationsInSwap,
        const IModel*        model,
        Control*             control); 


    static DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                             const CAsset*          asset,
                                             const DateTime&        valueDate,
                                             const DateTime&        maturityDate,
                                             const IModel*           model);

    static void validateBasis(const IModel*        model,
                              const MarketData*    market,
                              const CAssetWrapper& asset);
private:
    VarianceSwapUtil();
};


// all the shared code goes here on the shell class
class PRODUCTS_DLL VolVarShell: public Generic1Factor, 
                   virtual public LastSensDate,
                   virtual public ISensitiveStrikes {
public:
    static CClassConstSP const TYPE;

    class PRODUCTS_DLL DivAdjuster: public virtual DividendList::IDivAdjuster{
    public:
        // simple constructor
        DivAdjuster();

        /** Converts a continuous dividend to a $0 dividend in order for continuous divs
            to be ignored in div adjustment computation */
        virtual void adjustDividend(Dividend& div);
        
        virtual ~DivAdjuster();
    };

    /** instrument validation */
    virtual void Validate();

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    static DoubleArraySP getSensitiveStrikesHelper(OutputNameConstSP outputName,
                                                   const DateTime& valueDate,
                                                   const DateTime& maturityDate,
                                                   const CAsset*   asset,
                                                   const IModel*         model);

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    virtual DateTime endDate(const Sensitivity* sensControl) const;;

    // determine static hedge
    virtual DoubleArray staticHedge(
        CMarketDataSP      market,
        IModel*            model,
        const DoubleArray& strikes,
        const BoolArray&   isCall,
        double             loShares,
        double             hiShares);

protected:
    friend class VolVarShellHelper;
    friend class VarOptFP;
    friend class VarianceSwapCreditSupport;

    VolVarShell();
    VolVarShell(const VolVarShell& rhs);
    VolVarShell& operator=(const VolVarShell& rhs);

protected:
    VolVarShell(CClassConstSP clazz): Generic1Factor(clazz), 
                                      dontScaleByStrike(false),
                                      divAdjOnExDate(false),   // to be updated
                                      numPastReturns(0),
                                      numTotalReturns(0),
                                      isVanilla(false) {}

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;

    /** calculate historical vol */
    double historicalVol(int& numReturns) const;

    /** calculate historical vol, given a value date */
    double historicalVol(int&            numReturns,
                         const DateTime& valueDate,
                         double          spot) const; 

    bool isFwdStarting() const
    {
        return( samples[0].date > valueDate );
    }


    // fields common to vol & var swaps
    CashFlowArray samples;
    double        strikeVol;
    int           observationsPerYear;
    bool          subtractMeanVol;
    string        payoffType;
    bool          dontScaleByStrike;
    double        cap;
    bool          noDivAdj;
    bool          divAdjOnExDate;

    // Used by VanillaVarSwap only
    int           numPastReturns;     // number of returns up to next sample
    int           numTotalReturns;    // total (expected) number of returns over life of contract
    bool          isVanilla;

};


// VarCapCVModel for CV adjustment
class PRODUCTS_DLL VarCapCVModel: public CModel {
public:
    static CClassConstSP const TYPE;

    static const string DEFAULT;
    static const string CONTROL_VARIATE;
    static const string MULTIPLICATIVE_SCALING;
    static const string ADJUST_MEAN_VOL;

	class Product;

	class IIntoProduct;

    void validatePop2Object();

    virtual void Price(CInstrument* instrument,
                       CControl*    control,
                       CResults*    results);

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingDisallowed if either cap or swap does,
     * else riskMappingAllowed if either cap or swap does,
     * else riskMappingIrrelevant.
     *
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz);

    // for VarCapCVModel::IIntoProduct  
    static void loadIntoProduct(CClassSP& clazz);

    static IObject* defaultVarCapCVModel();

    // constructor 
    VarCapCVModel();
    VarCapCVModel(IModelSP  volModel,
                  IModelSP  fwdModel,
                  string    methodology);

    //registered fields
    IModelSP    volModel;      // FourierEngine
    IModelSP    fwdModel;      // ClosedFormIntegrateLN
    string      methodology;    // CONTROL_VARIATE or MULTIPLICATIVE_SCALING
};

typedef smartPtr<VarCapCVModel> VarCapCVModelSP;

// Should be implemented by the instrument - VarianceSwap
class PRODUCTS_DLL VarCapCVModel::IIntoProduct: virtual public CModel::IModelIntoProduct {
public:
    virtual Product* createProduct(const VarCapCVModel* model) const = 0;
};

// variance swap
class PRODUCTS_DLL VarianceSwap: public VolVarShell,
                    virtual public ISupportVegaMatrixLite,
                    virtual public CreditSupport::Interface,
                    virtual public NumericalIntegrationLN::IIntoProduct, 
                    virtual public ImpliedIntegration::IIntoProduct,
                    virtual public ClosedFormIntegrateLN::IIntoProduct,
                    virtual public FourierEngine::IIntoProduct,
                    virtual public IMCIntoProduct,
                    virtual public VarCapCVModel::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    /** Support VEGA_MATRIX_LITE */
    virtual void avoidVegaMatrixLite(const IModel* model);
    
    class VarSwapIntegrandWeight;
    DECLARE_REF_COUNT(VarSwapIntegrandWeight);

    void validatePop2Object();

    /** retrieve market data */
    virtual void GetMarket(const IModel*       model, 
							const CMarketDataSP market);

    // GetMarket method for VarCapCVModel
    void GetCVMarket(const VarCapCVModel*   cvModel, 
                     const CMarketDataSP    market); 

    CreditSupportSP createCreditSupport(CMarketDataSP market);
    
    /** Implementation of NumericalIntegrationLN::IntoProduct interface */
    virtual NumericalIntegrationLN::IProduct* createProduct(NumericalIntegrationLN* model) const ;

    /** Implementation of ImpliedIntegration::IntoProduct interface */
    virtual ImpliedIntegration::IProduct* createProduct(ImpliedIntegration* model) const ;

    /** Implementation of ClosedFormIntegrateLN::IntoProduct interface */
    virtual ClosedFormIntegrateLN::IProduct* createProduct(ClosedFormIntegrateLN* model) const ;

    /** Implementation of FourierEngine::IntoProduct interface */
    virtual FourierProduct* createProduct(const FourierEngine* model) const;

    /** Implementation of MonteCarlo::IntoProduct interface */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

    /** Implementation of VarCapCVModel::IntoProduct interface */
    virtual VarCapCVModel::Product* createProduct(const VarCapCVModel* model) const;

    double varFromCalls(double                           lowStrike, 
                        double                           highStrike, 
                        int                              nbSteps,
                        const DateTime&                  maturity,
                        double                           fwd,
                        string                           integrationMethod,
                        const bool                       allowNegFwdVar,
                        double                           absPrecisionAdjusted,
                        double                           relPrecisionAdjusted,
                        bool                             isDoingVegaMatrix,
                        VarSwapBasis::VarSwapBasisProcSP basisProc = VarSwapBasis::VarSwapBasisProcSP( ),
                        VanillaContractsRecorderSP       recorder = VanillaContractsRecorderSP(new VanillaContractsRecorder(1.0, false))) const;

    /** price given the future vol squared */
    void price(double varFuture, 
               IModel* model,
               Control* control, 
               Results* results, 
               VarSwapBasis::VarSwapBasisProcSP basisProc = VarSwapBasis::VarSwapBasisProcSP( ),
               VanillaContractsRecorderSP recorder = VanillaContractsRecorderSP(0)) const; 

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
        returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    static double futureVar(const CAsset*           asset,
                            const YieldCurve*       discount,
                            const DateTime&         valueDate,
                            const DateTime&         volDate,
                            const IModel*           model,
                            Control*                control); 

protected:
    friend class VarianceSwapHelper;
    friend class VarianceSwapNumerical;
    friend class VarOptFP;
    friend class VarSwapProduct_FE;
    friend class VarianceSwapClosedForm;
    friend class VarianceSwapMC;
    friend class VanillaVarSwap;
    friend class VanillaVSProd;
    friend class VarSwapCallsPutsAddin;
    friend class VarCapCVModel;
    friend class VarOptCVFP;
	friend class VarSwapHedgingSupport;
    friend class VarianceSwapLeastSquareFit;
    friend class VAsset;


    // Transient - only used in VarCapCVModel
    bool            useCV;
    CAssetWrapper   marketAsset;

    VarianceSwap():VolVarShell(TYPE) {
		payoffType = "FORWARD";
		cap = 2.5;
		noDivAdj = false;
	}

	VarianceSwap(CClassConstSP clazz): VolVarShell(clazz) {
		payoffType = "FORWARD";
		cap = 2.5;
		noDivAdj = false;
	}

    // Constructor used by VanillaVarSwap
    VarianceSwap(const DateTime& valueDate,
                 const DateTime& startDate, 
                 const bool fwdStarting, 
                 const double initialSpot, 
                 const string& ccyTreatment,
                 InstrumentSettlementSP instSettle, 
                 InstrumentSettlementSP premiumSettle, 
                 CAssetWrapper asset, 
                 YieldCurveWrapper discount,
                 CashFlowArray samples, 
                 const double strikeVol, 
                 const int observationsPerYear,
                 const bool subtractMeanVol, 
                 const string& payoffType, 
                 const bool dontScaleByStrike,
                 const double cap, 
                 const bool noDivAdj,
                 const bool divAdjOnExDate,
                 const bool isVanilla,
                 const int numPastReturns,
                 const int numTotalReturns); 

    VarianceSwap(const VarianceSwap& rhs);
    VarianceSwap& operator=(const VarianceSwap& rhs);

    static void load(CClassSP& clazz);

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;
};

typedef smartPtr<VarianceSwap> VarianceSwapSP;


// VanVSModel class - holds SPs to variance swap model (for variance swap)
// and fourier engine (for cap)
class PRODUCTS_DLL VanVSModel : public CModel {
public:
    static CClassConstSP const TYPE;

    static const string METHOD_NONE;
    
    virtual void validatePop2Object();
        
    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(VanVSModel* model,
                           Control*   control, 
                           Results*   results) const = 0;
        virtual ~IProduct(){};
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(const VanVSModel* model) const = 0;
    };

    // inherited from CModel 
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
     
    /** Method to set the domestic yield curve name in the mdf  */
    virtual void setDomesticYCName (string discountYieldCurveName) const;
    
    MarketObjectSP GetMarket(const MarketData*    market,
                             const string&        name,
                             const CClassConstSP& type) const;

    virtual void getMarket(const MarketData* market,
                           IInstrumentCollectionSP instruments);

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingDisallowed if either cap or swap does,
     * else riskMappingAllowed if either cap or swap does,
     * else riskMappingIrrelevant.
     *
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    // registration, invoked when class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(VanVSModel, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultVanVSModel);
        FIELD(varSwapModel,"Model to price Variance Swap");
        FIELD(capModel,"Model to price Cap");
        FIELD_MAKE_OPTIONAL(capModel);
        FIELD(methodology,"control variate method");
        FIELD_MAKE_OPTIONAL(methodology);
		FIELD(useCV,"whether use capCVModel to price Cap");
        FIELD_MAKE_TRANSIENT(useCV);
		FIELD(capCVModel,"Model to price Cap using control variate method");
        FIELD_MAKE_TRANSIENT(capCVModel);
    }

    // for VanVSModel::IIntoProduct  
    static void loadIntoProduct(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER_INTERFACE(VanVSModel::IIntoProduct, clazz);
        EXTENDS(Model::IModelIntoProduct);
    }

    static IObject* defaultVanVSModel(){
        return new VanVSModel();
    }
    // constructor 
    VanVSModel():
	CModel(TYPE), 
	capModel(0), 
	methodology(VanVSModel::METHOD_NONE), 
	useCV(false),
	capCVModel(0){}

    // constructor 
    VanVSModel(IModelSP varSwapModel);

    // registered fields
    IModelSP                varSwapModel;       // either ImpliedIntegration, NumericalIntegrationLN, ClosedFormIntegrateLN, 
                                                // MonteCarlo, FourierEngine
    IModelSP        capModel;                   // FourierEngine
    string          methodology;                // control variate method: "none" (no control variate), "control_variate", 
                                                // "mutilicative_scaling" and "adjust_mean_vol"
	// transient fields
	bool			useCV;						// whether to do control variate adjustment
	IModelSP		capCVModel;					// model to do control variate adjustment

protected:
    VanVSModel(CClassConstSP clazz): 
	CModel(clazz), 
	capModel(0),
	methodology(VanVSModel::METHOD_NONE), 
	useCV(false),
	capCVModel(0){}

private:

};
typedef smartPtr<VanVSModel> VanVSModelSP;
typedef smartConstPtr<VanVSModel> VanVSModelConstSP;

// Fourier payoff
class PRODUCTS_DLL VarOptFP: public FourierProduct, 
                public FourierProductIntegrator1D,
                public StFourierProductQuadVar,
                public FwdStFourierProductQuadVar {
private:
    const VarianceSwap* inst;
    DateTime			maturity;
    double				mult;
    double				years;
    double				totalYears;
    double				strike;
    double				volPast;
    double				pastWeight;
	bool				isVolSwap;
	bool				useMatytsinOnly;
	
public:
    // equivalent to InstIntoFourierProduct 
    VarOptFP(const VarianceSwap*	inst,
			 const bool				isVolSwap,
			 const bool				useMatytsinOnly,
             const DateTime&		matDate);

    virtual void validate(FourierProcess* process) const;

    virtual void price(const FourierEngine* model,
                       Control*             control, 
                       Results*             results);

    /** Constructs integrands for var option */
    Function1DDoubleArrayConstSP Integrand(const FourierEngine * model, 
                                           const Integrator1D*  integrator);

    const DateTime& getStartDate() const {return inst->samples[0].date;}
    
    /** Post process method for VolSwapFP */
    virtual void postResults(const FourierEngine* model,
                             const Integrator1D*  integrator,
                             const FourierProductIntegrator1D::IntegralArray& integrals,
                             CControl*            control,
                             CResults*            results);
};


DRLIB_END_NAMESPACE

#endif
