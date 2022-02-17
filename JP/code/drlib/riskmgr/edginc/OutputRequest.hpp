
//

#ifndef EDG_OUTPUT_REQUEST_HPP
#define EDG_OUTPUT_REQUEST_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/OutputName.hpp"

DRLIB_BEGIN_NAMESPACE
class Results;
class IModel;
class CInstrument;
class OutputRequestHelper;

FORWARD_DECLARE(OutputRequest)

/** Identifies additional outputs which don't require tweaking - for example,
    indicative vol, fwd at maturity, effective strike etc. These will usually
    be added to the result set during the first pricing. */
class RISKMGR_DLL OutputRequest: public CObject{
public:
    friend class OutputRequestHelper;
    static CClassConstSP const TYPE;
  
    // definition of all output requests
    static const string LAST_INSTR_DATE;
    static const string FWD_AT_MAT;
    static const string IND_VOL;
    static const string AVG_VOL;
    static const string AVG_CORR;
    static const string AVG_VARIANCE;
    static const string AVG_FWD;
    static const string EFFECTIVE_STRIKE;
    static const string CLEAN_PRICE;

    /** for Constant Maturity CDS*/
    static const string CMCDS_ADJUSTED_RATES;
    static const string CMCDS_CAP_PV;
    static const string CMCDS_FLOOR_PV;
    static const string CMCDS_DISCOUNT_FACTORS;
    static const string CMCDS_SURV_PROBS;
    static const string CMCDS_ACCS_ON_DEFAULT;
    static const string CMCDS_FORWARD_ADJS;
    static const string CMCDS_MAT_ADJS;
  
    /**"Clean Price" according to CDS high-yield quoting convention. This is 1-F-A, where
       F is the total fee payed by the protection buyer, and A is the (positive) accrued 
       interest. F and A aree both per unit notional. Thus a par CDS on a payment date has
       clean price 100%.*/
    static const string CDS_CLEAN_PRICE;
    static const string DIRTY_PRICE;
    static const string QUOTED_PRICE;
    static const string OPTION_PRICE;
    static const string OPTION_VEGA;
    /**Simple delta of option, where this is meaningful (i.e. European BS model) */
    static const string OPTION_DELTA;
    /**Simple gamma of option, where this is meaningful (i.e. European BS model) */
    static const string OPTION_GAMMA;
    static const string ACCRUED_INTEREST;
    static const string DELAY_PRICE;
    static const string CALC_TIME;
    static const string PRICE_TIME;
    static const string PRICE_PCT_PAYOFF;
    static const string MAX_PAYOFF;
    static const string VALUE_STE;
    static const string ESW_LEG_PRICE;
    static const string CURRENT_SPREAD;
    static const string IND_CDS_PAR_SPREAD;
    static const string STRIKE;
    static const string ADJUSTED_STRIKE;
    /**Adjustment to forward. e.g. for spread struck CDS options on index underlyings, this is
       the adjustment needed because of the flat curve CDSW pricing convention. */
    static const string FORWARD_ADJUSTMENT;
    static const string EXPECTED_STRIKE;
    static const string EXPECTED_SCALE_FACTOR;
    static const string SCALE_FACTOR_UPTO_DATE;
    static const string PAYMENT_DATES;
    static const string KNOWN_CASHFLOWS;
    static const string LEGAL_TERMS_FV;
    static const string EVENTS;
    /**Time in A/365 units (uses DateTime::yearFrac) from value date to exercise time, 
       conditional that exercise actually does occur.*/
    static const string EXPECTED_TIME_TO_EXERCISE;
    /**Risk-neutral probability that exercise occurs.*/
    static const string EXERCISE_PROBABILITY;


    static const string IW_IND_VOL;
    static const string IW_STRIKE;
    static const string IW_INTEREST;
    static const string IW_PUT;
    static const string IW_CAPITAL;
    static const string IW_PVDIVS;
    static const string IW_FUTSTRIKE;


    static const string TICK_VALUE;
    static const string CLEAN_DEFAULT_SPREAD_CURVE;
    static const string CLEAN_PRICE_FOR_SCALING;
    /**Default probability implied by an 'E2C' style model.*/
    static const string DEFAULT_PROBABILITY;
    /**Default probability implied by a regular risky discount curve, i.e. par spread curve.*/
    static const string IMPLIED_DEFAULT_PROBABILITY;
    static const string IMPLIED_CDS_SPREAD;
    static const string HIST_DISC_RATE_TO_MAT;
    static const string HIST_SPOT_FX;
    static const string HIST_DELTA;
    static const string HIST_DIV_DIFF;
    static const string HIST_DIV_DIFF_ANSI_DATE;

    static const string DRO_HIST_DISC_RATE_TO_MAT;
    static const string DRO_HIST_SPOT_FX;
    static const string DRO_HIST_DELTA;
    static const string DRO_HIST_DIV_DIFF;

    static const string STRIKE_ADJUSTMENT_UPTO_DATE;
    static const string CURRENT_STRIKE_ADJUSTMENT;
    static const string COMPUTE_ESTIMATE;
    static const string COMPUTE_INDEX;
    static const string THEO_PRICE;
    static const string RECOVERY_VALUE;
    static const string BARRIER_LEVEL;
    static const string BARRIER_DISCONTINUITY_RISK;
    
    /** for convertibles */
    static const string PUT_STRIKE;
    static const string CALL_STRIKE;
    static const string YIELD_TO_MATURITY;
    static const string YIELD_TO_FIRST_PUT;
    static const string PARITY;
    static const string PARITY_CONVERSION_RATIO;
    static const string THETA_ACCRUED_INTEREST;
    static const string PUT_PROBABILITY;
    static const string PUT_STOCK_LEVEL;
    static const string NAKED_BOND_PRICE;
    static const string NAKED_BOND_PRICE2; // price calculated from non static spread. ie no simple disc cashflow
    static const string BOND_DURATION;
    static const string BOND_CONVEXITY;
    static const string ASSET_SWAP_SPREAD;
    static const string CONT_CONVERT_STATUS;  // status of contingent convertibility

    /** for options on convertibles */
    static const string OPTION_ON_CONVERTIBLE_BOND_FLOOR;
    static const string OPTION_ON_CONVERTIBLE_LOCKOUT_FEE;
    static const string OPTION_ON_CONVERTIBLE_STRIKE;
    static const string CONVERTIBLE_PRICE_DIRTY;
    static const string CONVERTIBLE_PRICE_CLEAN;
    static const string OPTION_ON_CONVERTIBLE_INTRINSIC_VALUE;

    /** for vol/var swap */
    static const string VOL_IN_FUTURE;
    static const string VOL_IN_PAST;
    static const string PAST_WEIGHT;
    static const string DISCOUNT_FACTOR;
    static const string TOTAL_VOL;
    static const string STRIKE_VOL;
    static const string STRIKE_REF;
    static const string EXPECTED_N;
    static const string IMPLIED_N;
    static const string CV_COEFF;
    static const string VARSWAP_FV;
    static const string VARCAP_FV;
    static const string VAR_SWAP_VOL_BASIS;
    static const string VAR_SWAP_CUTOFF;
    static const string VAR_OPTION_IND_VOL;
    static const string FSA_VALUE;
    static const string FSA_PRR;

    /** For Corridor Variance Swap */
    static const string EXPECTED_PCT_TIME_IN_RANGE;
    static const string CORRIDOR_VARIANCE_VALUE;
    static const string CORRIDOR_ACCRUAL_VALUE;
    
    /** for VIX Future */
    static const string PARVOL_SWAPVAR;
    static const string PARVOL_SWAPVAR_SV ;
    static const string SQT_SWAPVAR_SV ;
    static const string PARVOL_SWAPVAR_VIX_BASIS0;

    /** for Tracer */
    static const string TRACER_CONVERT_PRICE;
    static const string TRACER_DECS_PRICE;

    /** for SPI */
    static const string SPI_GAP_RISK;
    static const string SPI_GAP_RISK_PROFILE;
    static const string SPI_DYN_BASKET;
    static const string SPI_EQUITY_LEVELS;
    static const string SPI_ALLOCS;
    static const string SPI_BOND;
    static const string SPI_BOND_TODAY;
    static const string SPI_BOND_ALLOC;
    static const string SPI_EQUITY_NAV;
    static const string SPI_BOND_NAV;
    static const string SPI_EQUITY_NAV_PCT;
    static const string SPI_BOND_NAV_PCT;
    static const string SPI_UNBAL_EXPO; 
    static const string SPI_SUST_EXPO;
    static const string SPI_TARGET_EXPO; 
    static const string SPI_BOND_FLOOR;//historic one
    static const string SPI_BOND_FLOOR_LEVEL; // next one
    static const string SPI_LOCKED_IN_VALUE;
    static const string SPI_REPORT_DATE;
    static const string SPI_REPORT;
    static const string SPI_GAP_EVENT_STATS;
    static const string SPI_LOAN_COST_RATE;

    /** for Credit Default Swaption */
    static const string BASIS_POINT_VALUE;
    static const string FORWARD_CDS_SPREAD;
    static const string CDS_VALUE;

    /** TaxWrapper */
    static const string TAX_GROSS_UP;
    static const string TAX_NEXT_CASHFLOW;

    /** CorporateBond */
    static const string RISKY_STREAM_DETAILS;

    /** DDE */
    static const string DDE_PACKET;
    static const string DDE_EQ_VOL;
    static const string DDE_SPRD_BBONE;
    static const string DDE_SPRD_VOL;
    static const string DDE_EQ_SPRD_COR;
    static const string DDE_SPRD_FUNC;
    static const string DDE_LAST_GOOD_CALIB_DATE;
    static const string DDE_CALIB_SCHEDULE;
    static const string DDE_SPRD_DELTA;
    static const string DDE_VOL_DELTA;
    static const string DDE_VOL_SKEW;
    static const string DDE_ATMVOL_DELTA;
    static const string DDE_ATMVOL_SKEW;
    static const string DDE_IMPLIED_VOL;

    /** EDS */
    static const string EDS_UPFRONT;
    static const string EDS_RUNNING_FEE;
    static const string EDS_RISKY_DURATION;
    static const string EDS_DELTA_CURVE;
    static const string EDS_VEGA_CURVE;

	/** for InsuranceAnnuityGMB */
    static const string FEE_PRICE;

    /** for var swap hedging */
    static const string PL_AVERAGE;
    static const string PL_DISTRIBUTION;
    static const string PL_DELTA_HEDGE;
    static const string PATHS;
    static const string HEDGES;

    /** CorpAct (RiskArb) */
    static const string PRORATED_CASH_VALUE;

    /** for general debug results */
    static const string DBG;

    /** for things that physically deliver */
    static const string PHYSICAL_DELIVERY;

    /** for NAPOLI */
    static const string COUPON_DUE;
    static const string ACCRUAL_CALENDAR;

    /** for EGK **/
    static const string FIXED_LEG_VALUE;
    static const string LIBOR_LEG_VALUE;
    static const string LIBOR_FUNDING_TWEAK;
    
    /** for CDO */
    static const string TRANCHE_EXPECTED_LOSS_CURVE;
    static const string TRANCHE_EXPECTED_LOSS_CURVE_WITH_RECOVERED_NOTIONAL;
    static const string TRANCHE_CONTINGENT_LEG_PRICE;
    static const string TRANCHE_CONTINGENT_LEG_FV;
    static const string TRANCHE_FEE_LEG_PRICE;
    static const string TRANCHE_FEE_LEG_FV;
    static const string CPTY_CREDIT_CHARGE;
    static const string FV_MINUS_CCC;
    static const string VALUE_WITHOUT_CCC;
    static const string PAR_SPREAD_CURVE;
    static const string TRANCHE_RISKY_DURATION;
    static const string TRANCHE_IMPLIED_SPREAD;
    static const string TRANCHE_LOWER_BC_BETAS;
    static const string TRANCHE_UPPER_BC_BETAS;
    static const string TRANCHE_OUTSTANDING_NOTIONAL;

    /** for Credit TARN */
    static const string REDEEMING_PROBABILITY;

    /** for ABS CDO */
    static const string TRANCHE_ABS_BOTTOM_LOSS_CURVE;
    static const string TRANCHE_ABS_TOP_LOSS_CURVE;
    static const string TRANCHE_ABS_NOTIONAL_CURVE;
    
    /** for Bond Future */
    static const string BOND_FUTURE_Z_SPREAD;
    static const string PVBP;

    /** for CDS */
    static const string FEE_LEG_RISKLESS_FV;
    static const string CDS_RISKY_DURATION;

    /** for Rates */
    static const string UNDERLYING;
    static const string FWD_UNDERLYING;
    static const string FLOAT_LEG_VALUE;
    static const string COMPLEX_LEG_VALUE;

    static const string MOD_CORR_MATRIX_SQ_ERROR;
    static const string CORR_REALIZED;
    static const string CORR_IMPLIED;

    /** for EnergySwap */
    static const string FIXED_LEG_DETAILS;
    static const string FLOATING_LEG_DETAILS;

    /** for CreditIndexSwaps */
    static const string INDEX_BASIS;

    /** for RiskMapping */
    static const string RISK_MAPPED_SENS;

    /** for RiskMapping tests */
    static const string DEBUG_RELEVANT_RISK_MAPPING_MATRICES;

    /* for Conditional Loss Mapping */
	static const string CONDITIONAL_LOSS_SAMPLES;
	static const string CONDITIONAL_LOSS_MAP;
	static const string FLAT_CDO;
	static const string FLAT_CDO_TRANCHE_DECOMPOSITION;

    /** For KProxyGenerator (Radar representation) purposes; added by Kranthi K. Gade (07/25/2006) **/
    static const string RADAR_REP;
    static const string RADAR_DIAGNOSTICS;
    static const string RADAR_PACKET;
    static const string DURATION_WEIGHTED_AVG;

    /** for the horror that is FFX */
    static const string FFX_UFV;
    static const string FFX_PV;
    static const string FWD_FX_RATE;
    static const string FFX_FUTURE_VALUE;
    static const string USD_DISC_RATE;

    /** for SpreadLossTree */
    static const string CONDITIONAL_FWD;

    static const string CONTINGENT_LEG_FV;
    static const string FEE_LEG_FV;
    static const string ACCRUED_INT_FV;

    /** constructor */
    OutputRequest(const string& requestName);

    /** "constructor" returning smartPtr */
    static OutputRequestSP SP(const string& requestName);

    //// creates deep copy
    virtual IObject* clone() const;

    //// ensures that the output request is a valid one
    virtual void validatePop2Object();

    virtual void setHasFinished(bool hasFinished);

    virtual bool getHasFinished() const;

    //// what this request is known as eg IND_VOL
    virtual const string& getRequestName() const;

    //// where in the Results the output for this request is stored
    virtual const string& getPacketName() const;
    
    //// does the result from this output request need to be scaled
    virtual bool isScaleable() const;

    //// does the result from this output request need to be added across
    //// composite instruments
    virtual bool isAdditive() const;

    //// does the result from this output request need to be scaled
    ////  in the post processing
    virtual bool isScaleablePostProcess() const;

    /** scale the results in the Results Object for this output request
        by supplied factor. The singleInstStatistics
        indicates whether we're trying to caclulate the mean value of
        a result by pricing the same instrument multiple times */
    virtual void scaleResult(Results*     results,     // (M)
                             double       scaleFactor,
                             bool         singleInstStatistics) const;

    /** Modify results in the Results Object for this output request by
        adding all results in resultsToAdd as indicated by control */
    virtual void addResult(
        Results*           results,     // (M)
        const Results*     resultsToAdd,
        double             scaleFactor,
        bool               sameInstrument) const; /* are the 2 results for
                                                     the same instrument */

    /** scale the results in the Results Object for this output request
        by supplied factor.  */
    virtual void scalePostProcess(Results*     results,     // (M)
                                  double       scaleFactor) const;



    /** return list of all possible requests */
    static CStringArraySP allRequests();

    // once we're through pricing run some requests maye not have been done
    // they're either non-applicable (default behaviour) or need to be done
    // we call the calculator in the helper class
    void handleUnfulfilled(const IModel* model, const CInstrument* instrument,
                            Results* results);

protected:
    /** In case we want to derive from this class */
    OutputRequest(const CClassConstSP& clazz,
                  const string&        requestName);

private:
    string               requestName;
    // $unregistered
    bool                 hasFinished;    // indicates whether request has been 
                                         // calculated successfully
    const OutputRequestHelper* data;     // $unregistered

    /** default constructor for reflection */
    OutputRequest();

    // not implemented
    OutputRequest(const OutputRequest &rhs);
    OutputRequest& operator=(const OutputRequest& rhs);
};

#ifndef QLIB_OUTPUTREQUEST_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<OutputRequest>);
EXTERN_TEMPLATE(class RISKMGR_DLL array<OutputRequestSP _COMMA_ OutputRequest>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<OutputRequestArray>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<OutputRequestArray>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<OutputRequest>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<OutputRequestSP _COMMA_ OutputRequest>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<OutputRequestArray>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<OutputRequestArray>);
#endif

/** Normally OutputRequests are calculated by instruments/models. However,
    some (eg Legal terms) can be calculated generically. In such cases 
    create a class derived from this one. The initialisation code in 
    OutputRequest.cpp should then specify the Calculator to use. The
    default one (here) just returns N/A */
class RISKMGR_DLL OutputRequestCalculator {
public:
    virtual ~OutputRequestCalculator();

    /** Calculate OutputRequest for given request, for supplied model and
        instrument. Default implementation here stores N/A */
    virtual void calculate(OutputRequest* request, const IModel* model, 
                           const CInstrument* instrument, Results* results);

    OutputRequestCalculator() {};
private:
    // not implemented
    OutputRequestCalculator(const OutputRequestCalculator &rhs);
    OutputRequestCalculator& operator=(const OutputRequestCalculator& rhs);
};
typedef refCountPtr<OutputRequestCalculator> OutputRequestCalculatorSP;

DRLIB_END_NAMESPACE

#endif
