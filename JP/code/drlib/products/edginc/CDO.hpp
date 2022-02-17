//----------------------------------------------------------------------------
//
// Group       : Derivatives Research
//
// Filename    : CDO.hpp
//
// Description : Collateralized Debt Obligation Instrument. An instance of
//               a GeneralisedCDO where the underlying portfolio is a tranche
//
// Author      : Mark A Robson/Sebastien Hitier/Antoine Gregoire
//
// Date        : 28 Oct 2004
//
//----------------------------------------------------------------------------

#include "edginc/GeneralisedCDO.hpp"
#include "edginc/IForwardRatePricer.hpp"


DRLIB_BEGIN_NAMESPACE


/** Credit Tranche Instrument. */
class PRODUCTS_DLL CDO: 
	public GeneralisedCDO
{
public:
    //CDO needs access
    friend class MyConvolutionProduct;

    static CClassConstSP const TYPE;

	virtual ~CDO();

    /** Public constructor */
    CDO(double                       lowStrike,
        double                       highStrike,
        bool                         isLong,
        CBoolSP                      recoverNotional,
        smartPtr<CounterPartyCredit> cptyInfo,
        ICreditContingentLegSP       cLeg,
        ICreditFeeLegSP              fLeg,
        CDOPortfolioSP               portfolio,
        const string&                yieldCurveName,
        CIntSP                       triggerDelay = CIntSP(),
        CIntSP                       defaultToCalculationDelay = CIntSP(),
        DateTime                     lastTriggerDate = DateTime(),
        CDoubleSP                    temporaryLossAmount = CDoubleSP(),
        const string&                settlementBDC = "N",
        HolidayWrapper               settlementHols = HolidayWrapper());

    virtual void validatePop2Object();

    //// function used only in the TrancheIndexLeastSquareFit
    virtual void getMaturityAndStrikesPercent(DateTime& t, double &k1, double &k2) const;

    //// function to get the inner CDOPortfolio
    virtual CDOPortfolioSP getCDOPortfolio() const;


private:
    // usual methods
    CDO(const CDO& rhs); // don't use
    CDO& operator=(const CDO& rhs); // don't use

	// CDO();

	CDO(CClassConstSP clazz = TYPE);

    static IObject* defaultConstructor();
    /** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);

    //------
    // Fields
    //------
    double    lowStrike;        /** lowStrike and highStrike defines portion of
                                    portfolio that payoff depends upon */
    double    highStrike;

    /**
     * Essentially a CDO with maturity / spread / conventions instead of
     * fee / contingent legs
     * */
    class QuickPricer: public KComponent,	// for tree pricing use KComponent instead of CInstrument
        public virtual LastSensDate,
        public virtual Theta::Shift,
        public virtual ConvolutionEngine::IIntoProduct,
		public virtual FDModel::IIntoProduct // for tree pricing
    {
    public:
        static CClassConstSP const TYPE;

        /** Called immediately after object constructed */
        virtual void validatePop2Object();

        /** Retrieve market data from cache for all our components */
        virtual void GetMarket(
            const IModel* model,
            const CMarketDataSP market);

		/** Setup the component */ 
		// The derived function must also call KComponent::setup() and the setup() of its underlyings
		// In the KComponent framework this method does everything usually found in GetMarket except
		// for populating market wrappers (this is done in KComponent::GetMarket)
		virtual void setup(const IModel* model, const MarketData* market); 

        /** Called after market data has been retrieved */
        virtual void Validate();

        //// Required part of CInstrument
        virtual DateTime getValueDate() const;

        /** Returns the name of the instrument's discount currency. */
        virtual string discountYieldCurveName() const;

        /**
         * when to stop tweaking (need to change infrastructure to route through
         * model rather than instrument)
         * */
        virtual DateTime endDate(
            const Sensitivity* sensControl) const;

        //// Required part of Theta::Shift
        virtual bool sensShift(
            Theta* shift);

    private:
        /** Creates an instance of an ConvolutionProduct */
        virtual IGeneralisedConvolutionProduct* createProduct(
            ConvolutionEngineConstSP model) const;

		/** Creates as instance of SpreadLossTree product */
		virtual FDProductSP createProduct(FDModel * model) const;

        /** Invoked once at start up when this class is 'loaded' */
        static void load(CClassSP& clazz);
        static IObject* defaultConstructor();

        /** Private constructor (for reflection) */
        QuickPricer();

        //-------
        // Fields
        //-------

        DateTime                     today;          // Today [optional]
        DateTime                     tradeDate;      // date at which trade was started [optional]
        double                       lowStrike;      // Tranche low strike (absolute value, not a percentage) [mandatory]
        double                       highStrike;     // Tranche high strike (absolute value, not a percentage) [mandatory]
        bool                         isLong;         // Long or short the tranche deal [mandatory]
        CBoolSP                      recoverNotional;// Do we recover notional from the top [optional]
        smartPtr<CounterPartyCredit> cptyInfo;       // Counterparty [optional]
        CDOPortfolioSP               portfolio;      // Tranche portfolio [mandatory]
        ICreditLegConventionSP       legConventions; // Conventions that define fee / contingent legs [mandatory]
        ExpirySP                     expiry;         // Tranche expiry [mandatory]
        double                       spread;         // Spread (fee leg fixed coupon) [mandatory]
        double                       upfrontPayment; // Upfront fee payment [optional]
        double                       legNotional;    // Fee and contingent leg notional [optional]
        CDoubleSP                    feeLegNotionalOverride; // Override for the fee leg notional [optional]

        smartPtr<CDO>                cdo;            // CDO representation of this CDO::QuickPricer
                                                     // [transient but tweakable]

        // Default handling parameters
        /** Delay between default and eventDeterminationDate */
        CIntSP triggerDelay;

        /** Delay between credit event and calculation date */
        CIntSP defaultToCalculationDelay;

        /** Last date when a default occurred during the protection period
         * can be triggered */
        DateTime lastTriggerDate;

        /** Assumption for the loss amount between a credit event and the
         * related calculation date */
        CDoubleSP temporaryLossAmount;

        /** Bad day convention for the CDO */
        string settlementBDC;

        /** Holidays */
        HolidayWrapper settlementHols;
    };

    // QuickPricer needs to have access to some CDO private methods
    friend class QuickPricer;

public:
	using GeneralisedCDO::createProduct;
};

typedef smartPtr<CDO> CDOSP;


DRLIB_END_NAMESPACE
