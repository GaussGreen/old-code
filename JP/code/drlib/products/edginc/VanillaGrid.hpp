//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VanillaGridGrid.hpp
//
//   Description : Grid of European vanillas
//
//   Author      : Regis Guichard
//
//   Date        : 27 May 02
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VANILLA_GRID_HPP
#define EDR_VANILLA_GRID_HPP
#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
//#include "edginc/SensitiveStrikes.hpp"
//#include "edginc/SampleList.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Results.hpp"
#include "edginc/Weightmatrix.hpp"

#include "edginc/FDModel.hpp"

DRLIB_BEGIN_NAMESPACE

class VanillaGrid;
typedef smartPtr<VanillaGrid> VanillaGridSP;

/** Grid of Vanilla's instrument. */
class PRODUCTS_DLL VanillaGrid: public CInstrument, 
                   public virtual CClosedFormLN::IIntoProduct, 
                   public virtual IMCIntoProduct,
                   public virtual FDModel::IIntoProduct,
                   public virtual FourierEngine::IIntoProduct,
        //#error "FIXME is this Theta::Shiftable or not??"
                   public virtual Theta::IShift,
//                 public virtual ISensitiveStrikes,
                   public virtual LastSensDate {
public:
    static CClassConstSP const TYPE;
    class LeastSquareFit;
    friend class LeastSquareFit;
    class LeastSquareSimple;
    friend class LeastSquareSimple;
    friend class VanillaGridProd;

    /** Container used to output the grid of prices or else */
    class PRODUCTS_DLL Output: public CObject{
    public:
        static CClassConstSP const TYPE;
       
        struct PRODUCTS_DLL ErrorCode{
            static const int IMP_VOL_ERROR;
            static const int IMP_VOL_NOT_REQUESTED;
            static const int PRICE_NOT_REQUESTED;
        };

        Output(int numCols, int numRows);

        int numRows() const;
        int numCols() const;
        
        void setValue(int iCol, int iRow, double value);
        double getValue(int iCol, int iRow) const;

        void setError(int iCol, int iRow, int code);
        int getError(int iCol, int iRow) const;

        bool isError(int iCol, int iRow) const;

    private:        
        Output();

        static void load(CClassSP& clazz);
        static IObject* defaultCtor();

        CDoubleMatrix values;
        IntArrayArray errors;
    };
    typedef smartPtr<Output> OutputSP;
    typedef smartConstPtr<Output> OutputConstSP;

    /**  LeastSquareFit objFunc class */
    class PRODUCTS_DLL LeastSquareFit: public Calibrator::ObjFuncLeastSquare,
                          public Calibrator::ObjFunc::IBootstrappable{
    public:
        friend class LeastSquareFitHelper;
        static CClassConstSP const TYPE;

        struct PRODUCTS_DLL FitType{
            enum {
                PRICE = 0,
                VEGA_NORMALIZED_PRICE,
                IMPLIED_VOL,
                NB_ENUMS
            };
            enum {
                defaultIndex = IMPLIED_VOL
            };
        };

        LeastSquareFit();

        LeastSquareFit(const IModel&       model,
                       const VanillaGrid& vanillaGrid,
                       const string&             fitType,
                       const string&             volType);

        virtual void validatePop2Object();

        /** Gives this class the chance to use the market */
        virtual void getMarket(MarketData* market);

        /** Gives this class a chance to do additional validation 
            after getMarket has been called */
        virtual void validate();

        /** Returns an IObect that contains all the IAdjustable objects 
            that the calibrator is to operate upon */
        virtual IObjectSP getAdjustableGroup();

        /** Returns the number of functions */
        virtual int getNbFuncs() const;

        /** Calculates the values of the objective functions */
        virtual void calcValue(CDoubleArray& funcvals) const;

        void getSensitivities(CDoubleMatrix& sensitivities,
                              CDoubleMatrix& valGuess,
                              DoubleArray shift,
                              int index) const;

        void calcInitialVals(CDoubleMatrix& valGuess,
                             int index) const;

        /** Makes additional validations for calibration with bootstrapping */
        virtual void validate(const Calibrator::InstanceID::IBootstrappableArray& ids) const;

        /** Gets the first maturity with non zero weight */
        const int getIdxFirstMat(int nberMat, int nberStrk) const;
        
        /** Updates the instrument for calibration with bootstrapping 
            before each calibration run */
        virtual void update(int idxSmileMat);

        /** Reset the instrument for calibration with bootstrapping 
            after each calibration run */
        virtual void reset();

        /** Given an array of maturities and an array of strikes, create a default
            weight matrix */
        static CDoubleMatrixSP createDefaultWeights(const DateTime&      baseDate,
                                                    double               spot,
                                                    const DateTimeArray& maturities,
                                                    const DoubleArray&   strikes);

        const static string VOLSURFACE;

        void useNormalizedWeights(string normalize);

        static const string NORMALIZE;
       

    private:
       
        /* Registered fields */
        IModelSP         model;
        VanillaGridSP    vanillaGrid;
        VanillaGridSP    vanillaGridOriginal;
        bool             useSameRandomNumbers;
        string           fitType;    
        OutputSP         mktVals;    // optional market/ref values (vols/prices)
        OutputSP         mktVegas;   // optional market vegas
        string           volType;

        /* Transient fields */
        CControlSP       control;
        int              nbFuncs;
        int              fitWhat;
        string           outReqName;
        string           normalizeWeights;
   
    };
    typedef smartPtr<LeastSquareFit> LeastSquareFitSP;

    /**  LeastSquareSimple objFunc class 
         Idea will be to act as a facade and delegate to LeastSquareFit where possible. */
    class PRODUCTS_DLL LeastSquareSimple: public Calibrator::ObjFuncLeastSquare,
                             public Calibrator::ObjFunc::IBootstrappable{
    public:
        friend class LeastSquareSimpleHelper;
        static CClassConstSP const TYPE;

        // XXX Shouldn't this be common with definitions in LeastSquareFit above?
        struct PRODUCTS_DLL FitType{
            enum {
                PRICE = 0,
                VEGA_NORMALIZED_PRICE,
                IMPLIED_VOL,
                NB_ENUMS
            };
            enum {
                defaultIndex = IMPLIED_VOL
            };
        };

        LeastSquareSimple();

        LeastSquareSimple(const IModel&                   model,
                          string                         fitType,
                          bool                           useSameRandomNumbers,
                          const CAssetWrapper&           asset,
                          const YieldCurveWrapper&       discount,
                          const InstrumentSettlementSP&  instSettle,
                          const WeightMatrixWrapper&     weightMatrix,
                          const string&                  volType);

        virtual void validatePop2Object();

        /** Gives this class the chance to use the market */
        virtual void getMarket(MarketData* market);

        /** Gives this class a chance to do additional validation 
            after getMarket has been called */
        virtual void validate();

        /** Returns an IObect that contains all the IAdjustable objects 
            that the calibrator is to operate upon */
        virtual IObjectSP getAdjustableGroup();

        /** Returns the number of functions */
        virtual int getNbFuncs() const;

        /** Calculates the values of the objective functions */
        virtual void calcValue(CDoubleArray& funcvals) const;

        void getSensitivities(CDoubleMatrix& sensitivities,
                              CDoubleMatrix& valGuess,
                              DoubleArray shift,
                              int index) const;

        void calcInitialVals(CDoubleMatrix& valGuess,
                             int index) const;

       
        /** Makes additional validations for calibration with bootstrapping */
        virtual void validate(const Calibrator::InstanceID::IBootstrappableArray& ids) const;

        /** Gets the first maturity with non zero weight */
        const int getIdxFirstMat(int nberMat, int nberStrk) const;
        
        /** Updates the instrument for calibration with bootstrapping 
            before each calibration run */
        virtual void update(int idxSmileMat);

        /** Reset the instrument for calibration with bootstrapping 
            after each calibration run */
        virtual void reset();

        IModelSP getModel();

        VanillaGridSP getVanillaGrid();
		WeightMatrixWrapper getWeightMatrix();
		CAssetWrapper getAsset();

        void useNormalizedWeights(string normalize);
        string getNormalizedWeights();
        
    private:
       
        /* Registered fields */
        IModelSP                       model;
        string                         fitType;    
        bool                           useSameRandomNumbers;
        CAssetWrapper                  asset;
        YieldCurveWrapper              discount;
        InstrumentSettlementSP         instSettle;
        WeightMatrixWrapper            weightMatrix;
        string                         volType;

        OutputSP               mktVals;    // optional market/ref values (vols/prices)
        OutputSP               mktVegas;   // optional market vegas

        /* Transient fields */
        VanillaGridSP    vanillaGrid;
        LeastSquareFitSP leastSquareFit;

        CControlSP       control;

        string normalizeWeightsLSS;
    };
    typedef smartPtr<LeastSquareSimple> LeastSquareSimpleSP;

#if 0
    static double priceSpread(const DateTime& valueDate,
                              const DateTime& startDate,
                              const DateTime& matDate,
                              bool isCall,
                              bool fwdStarting,
                              bool oneContract,
                              double notional,
                              double initialSpot,
                              double lowStrike,
                              double highStrike,
                              const InstrumentSettlement* instSettle,
                              const Asset* asset,
                              const YieldCurve* discount);
#endif

    /** instrument validation */
    virtual void Validate();

    /** input data validation */
    virtual void validatePop2Object();

    /** retrieve market data needed by VanillaGrid - just valueDate, asset and
        discount yield curve */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** Implementation of FD model interface */
    virtual FDProductSP createProduct(FDModel* model) const;

    /** Implementation of MonteCarlo::IntoProduct interface */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const;

    /** Implementation of FourierEngine::IntoProduct interface */
    virtual FourierProduct* createProduct(const FourierEngine* model) const;

#if 0
    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    virtual CSensControl* AlterControl(const IModel*       modelParams,
                                       const CSensControl* sensControl) const;
#endif

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    virtual DateTime endDate(const Sensitivity* sensControl) const;
    
    virtual DateTime getValueDate() const;

#if 0
    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;
#endif

    /** Make a simple started vanilla grid */
    VanillaGrid(const CAssetWrapper&          asset,
                const YieldCurveWrapper&      discount,
                const WeightMatrixWrapper&    weightMatrix,
                const InstrumentSettlementSP& instSettle = getDefaultSettlement(),
                const string&                 ccyTreatment = CAsset::CCY_TREATMENT_NONE);

    VanillaGrid(const DateTime&               valueDate,
                const CAssetWrapper&          asset,
                const YieldCurveWrapper&      discount,
                const WeightMatrixWrapper&     weightMatrix,
                const InstrumentSettlementSP& instSettle = getDefaultSettlement());

    VanillaGrid(bool                          fwdStarting,
                DateTime                      startDate,
                bool                          oneContract,
                double                        notional,
                double                        initialSpot,
                const CAssetWrapper&          asset,
                const YieldCurveWrapper&      discount,
                const WeightMatrixWrapper&     weightMatrix,
                const InstrumentSettlementSP& instSettle = getDefaultSettlement());

    // fully blown ctor
    VanillaGrid(bool                          fwdStarting,
                DateTime                      startDate,
                bool                          oneContract,
                double                        notional,
                double                        initialSpot,
                const DateTime&               valueDate,
                const CAssetWrapper&          asset,
                const string&                 ccyTreatment,
                const YieldCurveWrapper&      discount,
                const WeightMatrixWrapper&     weightMatrix,
                const InstrumentSettlementSP& instSettle,
                const InstrumentSettlementSP& premiumSettle);

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    /** The weight matrix */
    CDoubleMatrixConstSP getWeights() const;

    /** The notional */
    double getNotional() const;

private:
    friend class VanillaGridHelper;
    friend class VanillaGridClosedForm;
    friend class VanillaGridMC;
	friend class VanillaGridMCSV;
    friend class VanillaGridFP;
    friend class VanillaGridInstrumentCollection;

    VanillaGrid();
    VanillaGrid(CClassConstSP clazz);

    VanillaGrid(const VanillaGrid& rhs);
    VanillaGrid& operator=(const VanillaGrid& rhs);

    static void load(CClassSP& clazz);

    static InstrumentSettlementSP getDefaultSettlement();

    /* Registered fields */
    DateTime                startDate;
    bool                    isCall;
    bool                    fwdStarting;
    bool                    oneContract;
    double                  notional;
    double                  initialSpot;
    DateTime                valueDate;

    CAssetWrapper           asset;
    string                  ccyTreatment;

    YieldCurveWrapper       discount;
    InstrumentSettlementSP  instSettle;
    InstrumentSettlementSP  premiumSettle;
    WeightMatrixWrapper     weightMatrix;

    /*transient fields*/

    string                  strikeUnits;
    string                  instType;
    DateTimeArraySP         maturities;
    DoubleArraySP           strikes;
    CDoubleMatrixSP         weights;
    CDoubleMatrixSP         instStrikesUsed;
    IntArrayArraySP         instTypesUsed;

};

/** Fourier algorithm parameters for vanilla */
class PRODUCTS_DLL VanillaGridISAP: public CObject, 
                       public FourierEngine::ISAP {
public:
    static CClassConstSP const TYPE;
    friend class VanillaGridISAPHelper;
    friend class VanillaGridFP;
    
    void validatePop2Object(){
        static const string method = "VanillaGridISAP::validatePop2Object";    
        if(!Maths::isPositive(payoffToProcFreqBoundWeight) || !Maths::isPositive(1.0 - payoffToProcFreqBoundWeight)){
            throw ModelException(method,
                                 "payoffToProcFreqBoundWeight must be strictly between 0.0 and 1.0.");
        }
    }

protected:
    VanillaGridISAP():
    CObject(TYPE),
    useOneIntegral(true),
    payoffToProcFreqBoundWeight(0.9){}

private:
    bool useOneIntegral;
    double payoffToProcFreqBoundWeight;
};

typedef smartPtr<VanillaGridISAP> VanillaGridISAPSP;
typedef smartConstPtr<VanillaGridISAP> VanillaGridISAPConstSP;

DRLIB_END_NAMESPACE
#endif
