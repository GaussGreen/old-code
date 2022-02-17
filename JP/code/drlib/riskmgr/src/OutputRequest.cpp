

#include "edginc/config.hpp"
#define QLIB_OUTPUTREQUEST_CPP
#include "edginc/OutputRequest.hpp"
#include "edginc/Addin.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/EventsMgr.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

const string OutputRequest::CMCDS_ADJUSTED_RATES   = "CMCDS_ADJUSTED_RATES";
const string OutputRequest::CMCDS_CAP_PV           = "CMCDS_CAP_PV";
const string OutputRequest::CMCDS_FLOOR_PV         = "CMCDS_FLOOR_PV";
const string OutputRequest::CMCDS_DISCOUNT_FACTORS = "CMCDS_DISCOUNT_FACTORS";
const string OutputRequest::CMCDS_SURV_PROBS       = "CMCDS_SURV_PROBS";
const string OutputRequest::CMCDS_ACCS_ON_DEFAULT  = "CMCDS_ACCS_ON_DEFAULT";
const string OutputRequest::CMCDS_FORWARD_ADJS     = "CMCDS_FORWARD_ADJS";
const string OutputRequest::CMCDS_MAT_ADJS         = "CMCDS_MAT_ADJS";
const string OutputRequest::LAST_INSTR_DATE        = "LAST_INSTR_DATE";
const string OutputRequest::FWD_AT_MAT             = "FWD_AT_MAT";
const string OutputRequest::IND_VOL                = "IND_VOL";
const string OutputRequest::AVG_VOL                = "AVG_VOL";
const string OutputRequest::AVG_CORR               = "AVG_CORR";
const string OutputRequest::AVG_VARIANCE           = "AVG_VARIANCE";
const string OutputRequest::AVG_FWD                = "AVG_FWD";
const string OutputRequest::EFFECTIVE_STRIKE       = "EFFECTIVE_STRIKE";
const string OutputRequest::NAKED_BOND_PRICE       = "NAKED_BOND_PRICE";
const string OutputRequest::NAKED_BOND_PRICE2      = "NAKED_BOND_PRICE2";
const string OutputRequest::CDS_CLEAN_PRICE        = "CDS_CLEAN_PRICE";
const string OutputRequest::CLEAN_PRICE            = "CLEAN_PRICE";
const string OutputRequest::DIRTY_PRICE            = "DIRTY_PRICE";
const string OutputRequest::QUOTED_PRICE           = "QUOTED_PRICE";
const string OutputRequest::OPTION_PRICE           = "OPTION_PRICE";
const string OutputRequest::OPTION_VEGA            = "OPTION_VEGA";
const string OutputRequest::OPTION_GAMMA           = "OPTION_GAMMA";
const string OutputRequest::OPTION_DELTA           = "OPTION_DELTA";
const string OutputRequest::ACCRUED_INTEREST       = "ACCRUED_INTEREST";
const string OutputRequest::CONT_CONVERT_STATUS    = "CONT_CONVERT_STATUS";
const string OutputRequest::DELAY_PRICE            = "DELAY_PRICE";
const string OutputRequest::CALC_TIME              = "CALC_TIME";
const string OutputRequest::PRICE_TIME             = "PRICE_TIME";
const string OutputRequest::PRICE_PCT_PAYOFF       = "PRICE_PCT_PAYOFF";
const string OutputRequest::MAX_PAYOFF             = "MAX_PAYOFF";
const string OutputRequest::VALUE_STE              = "VALUE_STE";
const string OutputRequest::ESW_LEG_PRICE          = "ESW_LEG_PRICE";
const string OutputRequest::PARITY                 = "PARITY";
const string OutputRequest::PARITY_CONVERSION_RATIO= "PARITY_CONVERSION_RATIO";
const string OutputRequest::THETA_ACCRUED_INTEREST = "THETA_ACCRUED_INTEREST";
const string OutputRequest::PUT_PROBABILITY        = "PUT_PROBABILITY";
const string OutputRequest::PUT_STOCK_LEVEL        = "PUT_STOCK_LEVEL";
const string OutputRequest::CURRENT_SPREAD         = "CURRENT_SPREAD";
const string OutputRequest::IND_CDS_PAR_SPREAD     = "IND_CDS_PAR_SPREAD";
const string OutputRequest::PUT_STRIKE             = "PUT_STRIKE";
const string OutputRequest::CALL_STRIKE            = "CALL_STRIKE";
const string OutputRequest::STRIKE                 = "STRIKE";
const string OutputRequest::YIELD_TO_MATURITY      = "YIELD_TO_MATURITY";
const string OutputRequest::YIELD_TO_FIRST_PUT     = "YIELD_TO_FIRST_PUT";
const string OutputRequest::ADJUSTED_STRIKE        = "ADJUSTED_STRIKE";
const string OutputRequest::FORWARD_ADJUSTMENT     = "FORWARD_ADJUSTMENT";
const string OutputRequest::EXPECTED_STRIKE        = "EXPECTED_STRIKE";
const string OutputRequest::EXPECTED_SCALE_FACTOR  = "EXPECTED_SCALE_FACTOR";
const string OutputRequest::SCALE_FACTOR_UPTO_DATE = "SCALE_FACTOR_UPTO_DATE";
const string OutputRequest::PAYMENT_DATES          = "PAYMENT_DATES";
const string OutputRequest::KNOWN_CASHFLOWS        = "KNOWN_CASHFLOWS";
const string OutputRequest::BARRIER_LEVEL          = "BARRIER_LEVEL";
const string OutputRequest::LEGAL_TERMS_FV         = "LEGAL_TERMS_FV";
const string OutputRequest::EVENTS                 = "EVENTS";
const string OutputRequest::PVBP                   = "PVBP";
const string OutputRequest::EXPECTED_TIME_TO_EXERCISE = "EXPECTED_TIME_TO_EXERCISE";
const string OutputRequest::EXERCISE_PROBABILITY   = "EXERCISE_PROBABILITY";
const string OutputRequest::IW_INTEREST            = "IW_INTEREST";
const string OutputRequest::IW_PUT                 = "IW_PUT";
const string OutputRequest::IW_CAPITAL             = "IW_CAPITAL";
const string OutputRequest::IW_PVDIVS              = "IW_PVDIVS";
const string OutputRequest::IW_FUTSTRIKE           = "IW_FUTSTRIKE";


const string OutputRequest::BARRIER_DISCONTINUITY_RISK         =
"BARRIER_DISCONTINUITY_RISK";

const string OutputRequest::OPTION_ON_CONVERTIBLE_BOND_FLOOR   =
"OPTION_ON_CONVERTIBLE_BOND_FLOOR";
const string OutputRequest::OPTION_ON_CONVERTIBLE_LOCKOUT_FEE  =
"OPTION_ON_CONVERTIBLE_LOCKOUT_FEE";
const string OutputRequest::OPTION_ON_CONVERTIBLE_STRIKE       =
"OPTION_ON_CONVERTIBLE_STRIKE";
const string OutputRequest::CONVERTIBLE_PRICE_DIRTY            =
"CONVERTIBLE_PRICE_DIRTY";
const string OutputRequest::CONVERTIBLE_PRICE_CLEAN            =
"CONVERTIBLE_PRICE_CLEAN";
const string OutputRequest::OPTION_ON_CONVERTIBLE_INTRINSIC_VALUE       =
"OPTION_ON_CONVERTIBLE_INTRINSIC_VALUE";
const string OutputRequest::BOND_DURATION               = "BOND_DURATION";
const string OutputRequest::BOND_CONVEXITY              = "BOND_CONVEXITY";
const string OutputRequest::TICK_VALUE                  = "TICK_VALUE";
const string OutputRequest::CLEAN_DEFAULT_SPREAD_CURVE  =
 "CLEAN_DEFAULT_SPREAD_CURVE";
const string OutputRequest::VOL_IN_FUTURE               = "VOL_IN_FUTURE";
const string OutputRequest::VOL_IN_PAST                 = "VOL_IN_PAST";
const string OutputRequest::PAST_WEIGHT                 = "PAST_WEIGHT";
const string OutputRequest::DISCOUNT_FACTOR             = "DISCOUNT_FACTOR";
const string OutputRequest::TOTAL_VOL                   = "TOTAL_VOL";
const string OutputRequest::STRIKE_VOL                  = "STRIKE_VOL";
const string OutputRequest::STRIKE_REF                  = "STRIKE_REF";
const string OutputRequest::EXPECTED_N                  = "EXPECTED_N";
const string OutputRequest::IMPLIED_N                   = "IMPLIED_N";
const string OutputRequest::CV_COEFF                    = "CV_COEFF";
const string OutputRequest::VARSWAP_FV                  = "VARSWAP_FV";
const string OutputRequest::VARCAP_FV                   = "VARCAP_FV";
const string OutputRequest::VAR_SWAP_VOL_BASIS          = "VAR_SWAP_VOL_BASIS";
const string OutputRequest::VAR_SWAP_CUTOFF             = "VAR_SWAP_CUTOFF";
const string OutputRequest::VAR_OPTION_IND_VOL          = "VAR_OPTION_IND_VOL";
const string OutputRequest::EXPECTED_PCT_TIME_IN_RANGE  = "EXPECTED_PCT_TIME_IN_RANGE";
const string OutputRequest::CORRIDOR_VARIANCE_VALUE     = "CORRIDOR_VARIANCE_VALUE";
const string OutputRequest::CORRIDOR_ACCRUAL_VALUE      = "CORRIDOR_ACCRUAL_VALUE";
const string OutputRequest::IW_IND_VOL                  = "IW_IND_VOL";
const string OutputRequest::IW_STRIKE                   = "IW_STRIKE";


/** for VIX Future */
const string OutputRequest::PARVOL_SWAPVAR              = "PARVOL_SWAPVAR";
const string OutputRequest::PARVOL_SWAPVAR_SV           = "PARVOL_SWAPVAR_SV";
const string OutputRequest::SQT_SWAPVAR_SV              = "SQT_SWAPVAR_SV";
const string OutputRequest::PARVOL_SWAPVAR_VIX_BASIS0   = "PARVOL_SWAPVAR_VIX_BASIS0";

const string OutputRequest::CLEAN_PRICE_FOR_SCALING     = "CLEAN_PRICE_FOR_SCALING";
const string OutputRequest::DEFAULT_PROBABILITY         = "DEFAULT_PROBABILITY";
const string OutputRequest::IMPLIED_DEFAULT_PROBABILITY = "IMPLIED_DEFAULT_PROBABILITY";
const string OutputRequest::IMPLIED_CDS_SPREAD          = "IMPLIED_CDS_SPREAD";

const string OutputRequest::HIST_DISC_RATE_TO_MAT       = "HIST_DISC_RATE_TO_MAT";
const string OutputRequest::HIST_SPOT_FX                = "HIST_SPOT_FX";
const string OutputRequest::HIST_DELTA                  = "HIST_DELTA";
const string OutputRequest::HIST_DIV_DIFF               = "HIST_DIV_DIFF";
const string OutputRequest::HIST_DIV_DIFF_ANSI_DATE     = "HIST_DIV_DIFF_ANSI_DATE";

const string OutputRequest::DRO_HIST_DISC_RATE_TO_MAT       =
"DRO_HIST_DISC_RATE_TO_MAT";
const string OutputRequest::DRO_HIST_SPOT_FX                =
"DRO_HIST_SPOT_FX";
const string OutputRequest::DRO_HIST_DELTA                  =
"DRO_HIST_DELTA";
const string OutputRequest::DRO_HIST_DIV_DIFF               =
"DRO_HIST_DIV_DIFF";

const string OutputRequest::STRIKE_ADJUSTMENT_UPTO_DATE =
"STRIKE_ADJUSTMENT_UPTO_DATE";
const string OutputRequest::CURRENT_STRIKE_ADJUSTMENT   =
"CURRENT_STRIKE_ADJUSTMENT";
const string OutputRequest::COMPUTE_ESTIMATE     = "COMPUTE_ESTIMATE";
const string OutputRequest::COMPUTE_INDEX        = "COMPUTE_INDEX";
const string OutputRequest::THEO_PRICE           = "THEO_PRICE";
const string OutputRequest::TRACER_CONVERT_PRICE = "TRACER_CONVERT_PRICE";
const string OutputRequest::TRACER_DECS_PRICE    = "TRACER_DECS_PRICE";
const string OutputRequest::RECOVERY_VALUE       = "RECOVERY_VALUE";

const string OutputRequest::SPI_GAP_RISK              = "SPI_GAP_RISK";
const string OutputRequest::SPI_GAP_RISK_PROFILE      = "SPI_GAP_RISK_PROFILE";
const string OutputRequest::SPI_DYN_BASKET            = "SPI_DYN_BASKET";
const string OutputRequest::SPI_EQUITY_LEVELS         = "SPI_EQUITY_LEVELS" ;
const string OutputRequest::SPI_ALLOCS                = "SPI_ALLOCS";
const string OutputRequest::SPI_BOND                  = "SPI_BOND";
const string OutputRequest::SPI_BOND_TODAY            = "SPI_BOND_TODAY";
const string OutputRequest::SPI_BOND_ALLOC            = "SPI_BOND_ALLOC";
const string OutputRequest::SPI_EQUITY_NAV            = "SPI_EQUITY_NAV";
const string OutputRequest::SPI_BOND_NAV              = "SPI_BOND_NAV";
const string OutputRequest::SPI_EQUITY_NAV_PCT        = "SPI_EQUITY_NAV_PCT";
const string OutputRequest::SPI_BOND_NAV_PCT          = "SPI_BOND_NAV_PCT";
const string OutputRequest::SPI_UNBAL_EXPO            = "SPI_UNBAL_EXPO";
const string OutputRequest::SPI_SUST_EXPO             = "SPI_SUST_EXPO";
const string OutputRequest::SPI_TARGET_EXPO           = "SPI_TARGET_EXPO";
const string OutputRequest::SPI_BOND_FLOOR            = "SPI_BOND_FLOOR";
const string OutputRequest::SPI_BOND_FLOOR_LEVEL      = "SPI_BOND_FLOOR_LEVEL";
const string OutputRequest::SPI_LOCKED_IN_VALUE       = "SPI_LOCKED_IN_VALUE";
const string OutputRequest::SPI_REPORT_DATE           = "SPI_REPORT_DATE";
const string OutputRequest::SPI_REPORT                = "SPI_REPORT";
const string OutputRequest::SPI_GAP_EVENT_STATS       = "SPI_GAP_EVENT_STATS";
const string OutputRequest::SPI_LOAN_COST_RATE        = "SPI_LOAN_COST_RATE";

const string OutputRequest::TAX_GROSS_UP              = "TAX_GROSS_UP";
const string OutputRequest::TAX_NEXT_CASHFLOW         = "TAX_NEXT_CASHFLOW";
const string OutputRequest::ASSET_SWAP_SPREAD         = "ASSET_SWAP_SPREAD";

const string OutputRequest::BASIS_POINT_VALUE         = "BASIS_POINT_VALUE";
const string OutputRequest::CDS_VALUE                 = "CDS_VALUE";
const string OutputRequest::FORWARD_CDS_SPREAD        = "FORWARD_CDS_SPREAD";

const string OutputRequest::RISKY_STREAM_DETAILS      = "RISKY_STREAM_DETAILS";

const string OutputRequest::DDE_PACKET              = "DDE_PACKET";
const string OutputRequest::DDE_EQ_VOL              = "DDE_EQ_VOL";
const string OutputRequest::DDE_SPRD_BBONE          = "DDE_SPRD_BBONE";
const string OutputRequest::DDE_SPRD_VOL            = "DDE_SPRD_VOL";
const string OutputRequest::DDE_EQ_SPRD_COR         = "DDE_EQ_SPRD_COR";
const string OutputRequest::DDE_SPRD_FUNC           = "DDE_SPRD_FUNC";
const string OutputRequest::DDE_LAST_GOOD_CALIB_DATE = "DDE_LAST_GOOD_CALIB_DATE";
const string OutputRequest::DDE_CALIB_SCHEDULE      = "DDE_CALIB_SCHEDULE";
const string OutputRequest::DDE_SPRD_DELTA          = "DDE_SPRD_DELTA";
const string OutputRequest::DDE_VOL_DELTA           = "DDE_VOL_DELTA";
const string OutputRequest::DDE_VOL_SKEW            = "DDE_VOL_SKEW";
const string OutputRequest::DDE_ATMVOL_DELTA        = "DDE_ATMVOL_DELTA";
const string OutputRequest::DDE_ATMVOL_SKEW         = "DDE_ATMVOL_SKEW";
const string OutputRequest::DDE_IMPLIED_VOL         = "DDE_IMPLIED_VOL";

const string OutputRequest::EDS_UPFRONT        = "EDS_UPFRONT";
const string OutputRequest::EDS_RUNNING_FEE    = "EDS_RUNNING_FEE";
const string OutputRequest::EDS_RISKY_DURATION = "EDS_RISKY_DURATION";
const string OutputRequest::EDS_DELTA_CURVE    = "EDS_DELTA_CURVE";
const string OutputRequest::EDS_VEGA_CURVE     = "EDS_VEGA_CURVE";

const string OutputRequest::FEE_PRICE          = "FEE_PRICE";

const string OutputRequest::PL_AVERAGE          = "PL_AVERAGE";
const string OutputRequest::PL_DISTRIBUTION     = "PL_DISTRIBUTION";
const string OutputRequest::PL_DELTA_HEDGE      = "PL_DELTA_HEDGE";
const string OutputRequest::PATHS               = "PATHS";
const string OutputRequest::HEDGES              = "HEDGES";

// for CorpAct/RiskArb instrument
const string OutputRequest::PRORATED_CASH_VALUE = "PRORATED_CASH_VALUE";

const string OutputRequest::DBG = "DEBUG";
const string OutputRequest::PHYSICAL_DELIVERY = "PHYSICAL_DELIVERY";

const string OutputRequest::COUPON_DUE = "COUPON_DUE";
const string OutputRequest::ACCRUAL_CALENDAR = "ACCRUAL_CALENDAR";

const string OutputRequest::FIXED_LEG_VALUE = "FIXED_LEG_VALUE";
const string OutputRequest::LIBOR_LEG_VALUE = "LIBOR_LEG_VALUE";
const string OutputRequest::LIBOR_FUNDING_TWEAK = "LIBOR_FUNDING_TWEAK";

const string OutputRequest::TRANCHE_EXPECTED_LOSS_CURVE = "TRANCHE_EXPECTED_LOSS_CURVE";
const string OutputRequest::TRANCHE_EXPECTED_LOSS_CURVE_WITH_RECOVERED_NOTIONAL = "TRANCHE_EXPECTED_LOSS_CURVE_WITH_RECOVERED_NOTIONAL";
const string OutputRequest::TRANCHE_CONTINGENT_LEG_PRICE = "TRANCHE_CONTINGENT_LEG_PRICE";
const string OutputRequest::TRANCHE_CONTINGENT_LEG_FV = "TRANCHE_CONTINGENT_LEG_FV";
const string OutputRequest::TRANCHE_FEE_LEG_PRICE = "TRANCHE_FEE_LEG_PRICE";
const string OutputRequest::TRANCHE_FEE_LEG_FV = "TRANCHE_FEE_LEG_FV";
const string OutputRequest::CPTY_CREDIT_CHARGE = "CPTY_CREDIT_CHARGE";
const string OutputRequest::FV_MINUS_CCC = "FV_MINUS_CCC";
const string OutputRequest::VALUE_WITHOUT_CCC = "VALUE_WITHOUT_CCC";
const string OutputRequest::PAR_SPREAD_CURVE = "PAR_SPREAD_CURVE";
const string OutputRequest::TRANCHE_RISKY_DURATION = "TRANCHE_RISKY_DURATION";
const string OutputRequest::TRANCHE_IMPLIED_SPREAD = "TRANCHE_IMPLIED_SPREAD";
const string OutputRequest::TRANCHE_LOWER_BC_BETAS = "TRANCHE_LOWER_BC_BETAS";
const string OutputRequest::TRANCHE_UPPER_BC_BETAS = "TRANCHE_UPPER_BC_BETAS";
const string OutputRequest::TRANCHE_OUTSTANDING_NOTIONAL = "TRANCHE_OUTSTANDING_NOTIONAL";

const string OutputRequest::REDEEMING_PROBABILITY = "REDEEMING_PROBABILITY";

const string OutputRequest::TRANCHE_ABS_BOTTOM_LOSS_CURVE = "TRANCHE_ABS_BOTTOM_LOSS_CURVE";
const string OutputRequest::TRANCHE_ABS_TOP_LOSS_CURVE = "TRANCHE_ABS_TOP_LOSS_CURVE";
const string OutputRequest::TRANCHE_ABS_NOTIONAL_CURVE = "TRANCHE_ABS_NOTIONAL_CURVE";

const string OutputRequest::BOND_FUTURE_Z_SPREAD = "BOND_FUTURE_Z_SPREAD";
const string OutputRequest::FEE_LEG_RISKLESS_FV = "FEE_LEG_RISKLESS_FV";
const string OutputRequest::CDS_RISKY_DURATION = "CDS_RISKY_DURATION";
const string OutputRequest::FIXED_LEG_DETAILS = "FIXED_LEG_DETAILS";
const string OutputRequest::FLOATING_LEG_DETAILS = "FLOATING_LEG_DETAILS";
const string OutputRequest::FSA_VALUE = "FSA_VALUE";
const string OutputRequest::FSA_PRR = "FSA_PRR";

/** for Rates */
const string OutputRequest::UNDERLYING = "UNDERLYING";
const string OutputRequest::FWD_UNDERLYING = "FWD_UNDERLYING";
const string OutputRequest::FLOAT_LEG_VALUE = "FLOAT_LEG_VALUE";
const string OutputRequest::COMPLEX_LEG_VALUE = "COMPLEX_LEG_VALUE";

const string OutputRequest::MOD_CORR_MATRIX_SQ_ERROR = "MOD_CORR_MATRIX_SQ_ERROR";
const string OutputRequest::CORR_REALIZED = "CORR_REALIZED";
const string OutputRequest::CORR_IMPLIED = "CORR_IMPLIED";

const string OutputRequest::INDEX_BASIS = "INDEX_BASIS";

const string OutputRequest::RISK_MAPPED_SENS = "RISK_MAPPED_SENS";
const string OutputRequest::DEBUG_RELEVANT_RISK_MAPPING_MATRICES = "DEBUG_RELEVANT_RISK_MAPPING_MATRICES";

/*For Radar representation purposes */
const string OutputRequest::RADAR_REP = "RADAR_REP";
const string OutputRequest::RADAR_DIAGNOSTICS = "RADAR_DIAGNOSTICS";
const string OutputRequest::RADAR_PACKET = "RADAR";
const string OutputRequest::DURATION_WEIGHTED_AVG = "DURATION_WEIGHTED_AVG";

const string OutputRequest::FFX_UFV          = "FFX_UFV";
const string OutputRequest::FFX_PV           = "FFX_PV";
const string OutputRequest::FWD_FX_RATE      = "FWD_FX_RATE";
const string OutputRequest::USD_DISC_RATE    = "USD_DISC_RATE";
const string OutputRequest::FFX_FUTURE_VALUE = "FFX_FUTURE_VALUE";

const string OutputRequest::CONDITIONAL_LOSS_SAMPLES 
												= "CONDITIONAL_LOSS_SAMPLES";

const string OutputRequest::CONDITIONAL_LOSS_MAP     = "CONDITIONAL_LOSS_MAP";
const string OutputRequest::FLAT_CDO                 = "FLAT_CDO";

const string OutputRequest::FLAT_CDO_TRANCHE_DECOMPOSITION 
											= "FLAT_CDO_TRANCHE_DECOMPOSITION";

/** for SpreadLossTree */
const string OutputRequest::CONDITIONAL_FWD = "CONDITIONAL_FWD";


const string OutputRequest::CONTINGENT_LEG_FV = "CONTINGENT_LEG_FV";
const string OutputRequest::FEE_LEG_FV = "FEE_LEG_FV";
const string OutputRequest::ACCRUED_INT_FV = "ACCRUED_INT_FV";

typedef map<string, const OutputRequestHelper*> OutputRequestMap;

OutputRequestCalculator::~OutputRequestCalculator(){}

//// default implementation is just to return N/A
void OutputRequestCalculator::calculate(
    OutputRequest* request, const IModel* model,
    const CInstrument* instrument, Results* results) {
    results->storeNotApplicable(request);
}

class OutputRequestHelper{
public:
    /** Rules for combining results when the same instrument is priced
        multiple times and the results are aggregated  (eg
        splitting MC into blocks) */
    enum SameInstRule{
        NONE,        /* Just take 1st result and leave alone eg IND VOL */
        STATISTICAL, /* (eg a probability) we should scale and add
                        this result  */
        ADD          /* (eg price time) we should just add this result */

    };
    // **** fields ******
    // instances of this class are used to hold the details
    // of known OutputRequests
    // note const string references - this allows to define instances via
    // initialisation (this happens before main() is run). Cannot copy the
    // the strings since they are defined from strings in different files
    // and may have not been initialised before the code here is run
    const string& requestName;
    const string& packetName;
    const bool    scale;    // does this output scale (eg DELAY_PRICE)
    const bool    additive; // add for composites (eg ACCRUED_INTEREST)
    const bool    addNA;    // if additive false, does N/A + result = result?
    const bool    scalePostProcessing;  /* does this output scale for
                                           (converts) post-processing? */
    const SameInstRule sameInstRule; // combining results for the same inst

    OutputRequestCalculatorSP    calculator; // for handling 'top level' requests

    //// creates an instance of an OutputRequestHelper with everything
    //// specified and the default calculator
    OutputRequestHelper(
        const string& requestName,// eg FWD_AT_MAT
        // packetName: where in results it goes
        const string& packetName = Results::INSTRUMENT_PACKET,
        bool          scale    = false,    // scale
        bool          additive = false, // add for composites
        bool          scalePostProcessing = false,
        bool          addNA = false, // does N/A + result = result?
        SameInstRule  sameInstRule = NONE): // need to take mean?
        requestName(requestName), packetName(packetName),
        scale(scale), additive(additive), addNA(addNA),
        scalePostProcessing(scalePostProcessing),
        sameInstRule(sameInstRule),
        calculator(new OutputRequestCalculator()) {}

    //// creates an instance of an OutputRequestHelper with everything
    //// specified including a calculator
    OutputRequestHelper(
        OutputRequestCalculator* myCalculator,
        const string& requestName,// eg FWD_AT_MAT
        // packetName: where in results it goes
        const string& packetName = Results::INSTRUMENT_PACKET,
        bool          scale    = false,    // scale
        bool          additive = false, // add for composites
        bool          scalePostProcessing = false,
        bool          addNA = false, // does N/A + result = result?
        SameInstRule  sameInstRule = NONE): // need to take mean?
        requestName(requestName), packetName(packetName),
        scale(scale), additive(additive), addNA(addNA),
        scalePostProcessing(scalePostProcessing),
        sameInstRule(sameInstRule), calculator(myCalculator) {}

    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(OutputRequest, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultOutputRequest);
        FIELD(requestName, "Request's name");

        // build up the map of all known outputs;
        for (int i = 0; i < totalNumRequests; i++){
            allRequestsMap[allRequests[i].requestName] = &allRequests[i];
        }
    }

    static IObject* defaultOutputRequest(){
        return new OutputRequest();
    }

    // map of all known output requests
    static map<string, const OutputRequestHelper*> allRequestsMap;
    // array of all known output requests for easy initialisation
    static const OutputRequestHelper allRequests[];
    static const int totalNumRequests;
};

/** constructor */
OutputRequest::OutputRequest(const string&        requestName):
    CObject(TYPE), requestName(requestName), hasFinished(false), data(0) {
    validatePop2Object();
}

OutputRequestSP OutputRequest::SP(const string& requestName) {
    return OutputRequestSP(new OutputRequest(requestName));
}

/** In case we want to derive from this class */
OutputRequest::OutputRequest(const CClassConstSP& clazz,
                             const string&        requestName):
    CObject(clazz), requestName(requestName), hasFinished(false), data(0){
    validatePop2Object();
}

/** default constructor for reflection */
OutputRequest::OutputRequest():CObject(OutputRequest::TYPE),
                               hasFinished(false) {}

//// creates deep copy - need to override default to copy data pointer over
IObject* OutputRequest::clone() const{
    OutputRequest* copy = new OutputRequest();
    copy->requestName = requestName;
    copy->hasFinished = hasFinished;
    copy->data = data;
    return copy;
}

//// ensures that the output request is a valid one
void OutputRequest::validatePop2Object(){
    OutputRequestMap::const_iterator pos =
OutputRequestHelper::allRequestsMap.find(requestName);
    if (pos == OutputRequestHelper::allRequestsMap.end()) {
        throw ModelException("OutputRequest::validatePop2Object", requestName+
                             " is not a valid output request");
    }
    data = pos->second;
}

void OutputRequest::setHasFinished(bool hasFinished){
    this->hasFinished = hasFinished;
}

bool OutputRequest::getHasFinished() const{
    return hasFinished;
}

const string& OutputRequest::getRequestName() const{
    return requestName;
}

const string& OutputRequest::getPacketName() const{
    return data->packetName;
}

//// does the result from this output request need to be scaled
bool OutputRequest::isScaleable() const{
    return data->scale;
}

//// does the result from this output request need to be added across
//// composite instruments
bool OutputRequest::isAdditive() const{
    return data->additive;
}

//// does the result from this output request need to be scaled
////  in the post processing
bool OutputRequest::isScaleablePostProcess() const{
    return data->scalePostProcessing;
}

/** scale the results in the Results Object for this output request
    by supplied factor.  */
void OutputRequest::scaleResult(Results*     results,     // (M)
                                double       scaleFactor,
                                bool         singleInstStatistics) const{
    if ((!singleInstStatistics && data->scale) ||
        (singleInstStatistics &&
         data->sameInstRule == OutputRequestHelper::STATISTICAL)){
        if (data->packetName == requestName){
            // scale all results in packetName
            results->scale(requestName, scaleFactor);
        } else {
            // scale result requestName in packet packetName
            results->scale(data->packetName, requestName, scaleFactor);
        }
    }
}

/** scale the results in the Results Object for this output request
    by supplied factor.  */
void OutputRequest::scalePostProcess(Results*     results,     // (M)
                                     double       scaleFactor) const
{
    if (data->scalePostProcessing){
        if (data->packetName == requestName){
            // scale all results in packetName
            results->scale(requestName, scaleFactor);
        } else {
            // scale result requestName in packet packetName
            results->scale(data->packetName, requestName, scaleFactor);
        }
    }
}


/** Modify results in the Results Object for this output request by
    adding all results in resultsToAdd as indicated by control */
void OutputRequest::addResult(
    Results*           results,     // (M)
    const Results*     resultsToAdd,
    double             scaleFactor,
    bool               sameInstrument) const{ /* are the 2 results for
                                                 the same instrument */
    bool additive = data->additive;
    bool scale = data->scale;
    if (sameInstrument){
        if (data->sameInstRule == OutputRequestHelper::NONE){
            if (results->exists(this)){
                return; // nothing needs doing
            }
            scale = false; /* we are going to add this result in once (to
                              an empty result) so do not scale */
        } else {
            // eg consider prob of hitting barrier (scale = true) or
            // price time (scale = false)
           scale = data->sameInstRule == OutputRequestHelper::STATISTICAL;
        }
        // add the result in
        additive = true;
    }
    else if (!additive && data->addNA){
        // might be additive if results are missing or not applicable
        bool resultExists = results->exists(this);
        const NotApplicable* na = 0;
        // add the result if the result doesn't exist or if the result
        // is not applicable and the cause of it being not applicable is
        // not related to the output aggregation (otherwise you get, eg,
        // double + double -> N/A and then N/A + double -> double)
        if (!resultExists ||
            ((na = results->isNotApplicable(this)) ||
             (na = resultsToAdd->isNotApplicable(this))) &&
            na->getCause() != NotApplicable::aggregation){
            additive = true;
        }
    }
    if (additive){
        if (!scale) {
            // if we can't scale it, make sure we don't change it
            scaleFactor = 1.0;
        }
        if (data->packetName == requestName){
            // add all results in packetName
            results->add(requestName, resultsToAdd, scaleFactor);
        } else {
            // add result requestName in packet packetName
            results->add(data->packetName, requestName,
                         resultsToAdd, scaleFactor);
        }
    } else {
        results->storeNotApplicable(this, NotApplicable::aggregation);
    }
}

/** return list of all possible requests */
CStringArraySP OutputRequest::allRequests() {
    static const string method = "OutputRequest::allRequests";
    CStringArraySP allOfThem(new CStringArray(0));
    try {
        OutputRequestMap::const_iterator pos =
            OutputRequestHelper::allRequestsMap.begin();

        while (!(pos == OutputRequestHelper::allRequestsMap.end())) {
            allOfThem->push_back(pos->first);
            ++pos;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
    sort(allOfThem->begin(), allOfThem->end());
    return allOfThem;
}

void OutputRequest::handleUnfulfilled(const IModel* model,
                                        const CInstrument* instrument, Results* results) {
    try {
        data->calculator->calculate(this, model, instrument, results);
    } catch (exception& e) {
        results->storeRequestResult(this, IObjectSP(new Untweakable(e)));
    }
}

/** this is the list of all known output requests. Before adding another
    you have to decide where does the result live (the packet name),
    does the result scale (when multipled by 'analytics multiplier'), and
    does the result need to be combined across composite instruments? */
const OutputRequestHelper OutputRequestHelper::allRequests[] = {
    OutputRequestHelper(OutputRequest::LAST_INSTR_DATE,
                        Results::INSTRUMENT_PACKET, 
                        false, false, true),
    OutputRequestHelper(OutputRequest::FWD_AT_MAT,
                        OutputRequest::FWD_AT_MAT, // own packet
                        false, false, false),
    OutputRequestHelper(OutputRequest::IND_VOL,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* doesn't add in general */,
                        false /* doesn't scale pp*/,
                        true),/* but N/A + ind vol = ind vol */
    OutputRequestHelper(OutputRequest::AVG_VOL),
    OutputRequestHelper(OutputRequest::AVG_CORR),
    OutputRequestHelper(OutputRequest::AVG_VARIANCE),
    OutputRequestHelper(OutputRequest::AVG_FWD),
    OutputRequestHelper(OutputRequest::EFFECTIVE_STRIKE),
    OutputRequestHelper(OutputRequest::NAKED_BOND_PRICE,
                        Results::INSTRUMENT_PACKET, false, false, true),
    OutputRequestHelper(OutputRequest::NAKED_BOND_PRICE2,
                        Results::INSTRUMENT_PACKET, false, false, true),
    OutputRequestHelper(OutputRequest::CDS_CLEAN_PRICE,
                        Results::INSTRUMENT_PACKET, false, false, true),
    OutputRequestHelper(OutputRequest::CLEAN_PRICE,
                        Results::INSTRUMENT_PACKET, false, false, true),
    OutputRequestHelper(OutputRequest::DIRTY_PRICE,
                        Results::INSTRUMENT_PACKET, false, false, true),
    OutputRequestHelper(OutputRequest::QUOTED_PRICE,
                        Results::INSTRUMENT_PACKET, false, false, true),
    OutputRequestHelper(OutputRequest::OPTION_PRICE,
                        Results::INSTRUMENT_PACKET, false, false,
                        true), //// does this scale?
    OutputRequestHelper(OutputRequest::OPTION_VEGA,
                        Results::INSTRUMENT_PACKET, false, false,
                        true), //// does this scale?
	OutputRequestHelper(OutputRequest::OPTION_DELTA,
                        Results::INSTRUMENT_PACKET, false, false,
                        true), 
	OutputRequestHelper(OutputRequest::OPTION_GAMMA,
                        Results::INSTRUMENT_PACKET, false, false,
                        true),
    OutputRequestHelper(OutputRequest::TICK_VALUE),
    OutputRequestHelper(OutputRequest::ACCRUED_INTEREST,
                        Results::INSTRUMENT_PACKET,
                        true, true, true),
    OutputRequestHelper(OutputRequest::CONT_CONVERT_STATUS,
                        Results::INSTRUMENT_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DELAY_PRICE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* but not additive */,
                        true,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::PRICE_PCT_PAYOFF,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* but adds */, false),
    OutputRequestHelper(OutputRequest::MAX_PAYOFF,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* but adds */, false),
    OutputRequestHelper(OutputRequest::CALC_TIME,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale*/, true /* but adds */,
                        false,
                        false /* don't add NA */,
                        ADD /* add for same inst, but don't scale */),
    OutputRequestHelper(OutputRequest::PRICE_TIME,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale*/, true /* but adds */,
                        false,
                        false /* don't add NA */,
                        ADD /* add for same inst, but don't scale */),
    OutputRequestHelper(OutputRequest::COMPUTE_ESTIMATE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale*/, true /* but adds */,
                        false,
                        false /* don't add NA */,
                        ADD /* add for same inst, but don't scale */),
    OutputRequestHelper(OutputRequest::COMPUTE_INDEX,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale*/, false /* or adds */,
                        false), // MC takes care of aggregating this
    OutputRequestHelper(OutputRequest::VALUE_STE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false), // MC takes care of aggregating this
    OutputRequestHelper(OutputRequest::ESW_LEG_PRICE,
                        OutputRequest::ESW_LEG_PRICE,  // own packet
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::PARITY,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::PARITY_CONVERSION_RATIO,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::THETA_ACCRUED_INTEREST,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, true),
    OutputRequestHelper(OutputRequest::PUT_PROBABILITY,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::PUT_STOCK_LEVEL,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::CURRENT_SPREAD,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::IND_CDS_PAR_SPREAD,
                        OutputRequest::IND_CDS_PAR_SPREAD, /* own packet */
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::PUT_STRIKE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::CALL_STRIKE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::STRIKE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::YIELD_TO_MATURITY,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::YIELD_TO_FIRST_PUT,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::ADJUSTED_STRIKE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
	OutputRequestHelper(OutputRequest::FORWARD_ADJUSTMENT,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::EXPECTED_STRIKE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::EXPECTED_SCALE_FACTOR,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::SCALE_FACTOR_UPTO_DATE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::BARRIER_LEVEL,
                        OutputRequest::BARRIER_LEVEL, // own packet
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::BARRIER_DISCONTINUITY_RISK,
                        Results::INSTRUMENT_PACKET, // own packet
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::PAYMENT_DATES,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        true /* but "additive" */, false),
    OutputRequestHelper(OutputRequest::KNOWN_CASHFLOWS,
                        OutputRequest::KNOWN_CASHFLOWS, // own packet
                        true /* scales */, true /* and additive */, false),
    OutputRequestHelper(OutputRequest::OPTION_ON_CONVERTIBLE_BOND_FLOOR,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not "additive" */, true),
    OutputRequestHelper(OutputRequest::OPTION_ON_CONVERTIBLE_LOCKOUT_FEE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not "additive" */, true),
    OutputRequestHelper(OutputRequest::OPTION_ON_CONVERTIBLE_STRIKE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not "additive" */, true),
    OutputRequestHelper(OutputRequest::CONVERTIBLE_PRICE_DIRTY,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* but "additive" */, true),
    OutputRequestHelper(OutputRequest::CONVERTIBLE_PRICE_CLEAN,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* but "additive" */, true),
    OutputRequestHelper(OutputRequest::OPTION_ON_CONVERTIBLE_INTRINSIC_VALUE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not "additive" */, true),
    OutputRequestHelper(OutputRequest::BOND_DURATION,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, true),
    OutputRequestHelper(OutputRequest::BOND_CONVEXITY,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, true),
    OutputRequestHelper(OutputRequest::CLEAN_DEFAULT_SPREAD_CURVE,
                        OutputRequest::CLEAN_DEFAULT_SPREAD_CURVE, // own packet
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::VOL_IN_FUTURE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::VOL_IN_PAST,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::PAST_WEIGHT,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::DISCOUNT_FACTOR,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::TOTAL_VOL,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::STRIKE_VOL,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::STRIKE_REF,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::EXPECTED_N,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::IMPLIED_N, 
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::CV_COEFF,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::VARSWAP_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */,
                        true /* is additive */, false),
    OutputRequestHelper(OutputRequest::VARCAP_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */,
                        true /* is additive */, false),
    OutputRequestHelper(OutputRequest::VAR_SWAP_VOL_BASIS,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */,
                        false /* is additive */,
                        false /* doesn't scale pp*/,
                        true),/* but N/A + basis = basis */     // Used for composite VSWs
    OutputRequestHelper(OutputRequest::VAR_SWAP_CUTOFF,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */,
                        false /* is additive */,
                        false /* doesn't scale pp*/,
                        true),/* but N/A + cutoff = cutoff */   // Used for composite VSWs
    OutputRequestHelper(OutputRequest::VAR_OPTION_IND_VOL,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */,
                        false /* is additive */, false),
    OutputRequestHelper(OutputRequest::EXPECTED_PCT_TIME_IN_RANGE,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */,
                        false /* is additive */, false),
    OutputRequestHelper(OutputRequest::CORRIDOR_VARIANCE_VALUE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */,
                        true /* is additive */, 
                        true),
    OutputRequestHelper(OutputRequest::CORRIDOR_ACCRUAL_VALUE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */,
                        true /* is additive */, true
                        ),
    //VIX Future                                                  // Used for composite VIX Future
    OutputRequestHelper(OutputRequest::PARVOL_SWAPVAR,    
                            Results::INSTRUMENT_PACKET,
                            false /* doesn't scale */,
                            false /* not additive */, false),
    OutputRequestHelper(OutputRequest::PARVOL_SWAPVAR_SV,
                            Results::INSTRUMENT_PACKET,
                            false /* doesn't scale */,
                            false /* not additive */, false),
    OutputRequestHelper(OutputRequest::SQT_SWAPVAR_SV,
                            Results::INSTRUMENT_PACKET,
                            false /* doesn't scale */,
                            false /* not additive */, false),
    OutputRequestHelper(OutputRequest::PARVOL_SWAPVAR_VIX_BASIS0,
                            Results::INSTRUMENT_PACKET,
                            false /* doesn't scale */,
                            false /* not additive */, false),
    OutputRequestHelper(OutputRequest::HIST_DISC_RATE_TO_MAT,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::HIST_SPOT_FX,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::HIST_DELTA,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::HIST_DIV_DIFF,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::HIST_DIV_DIFF_ANSI_DATE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::DRO_HIST_DISC_RATE_TO_MAT,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::DRO_HIST_SPOT_FX,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::DRO_HIST_DELTA,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::DRO_HIST_DIV_DIFF,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::STRIKE_ADJUSTMENT_UPTO_DATE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::CURRENT_STRIKE_ADJUSTMENT,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::DEFAULT_PROBABILITY,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::IMPLIED_DEFAULT_PROBABILITY,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::IMPLIED_CDS_SPREAD,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::THEO_PRICE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */,
                        false /* but not additive */, true),
    OutputRequestHelper(OutputRequest::TRACER_CONVERT_PRICE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::TRACER_DECS_PRICE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    OutputRequestHelper(OutputRequest::RECOVERY_VALUE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, false),
    // ---------------- SPI stuff ------------------
    OutputRequestHelper(OutputRequest::SPI_GAP_RISK,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_GAP_RISK_PROFILE,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_DYN_BASKET,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_EQUITY_LEVELS,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_ALLOCS,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_BOND,
                        Results::INSTRUMENT_PACKET,
                        false  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_BOND_TODAY,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_BOND_ALLOC,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_EQUITY_NAV,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_BOND_NAV,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_EQUITY_NAV_PCT,
                        Results::INSTRUMENT_PACKET,
                        false  /* does NOT scale */,
                        false  /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_BOND_NAV_PCT,
                        Results::INSTRUMENT_PACKET,
                        false /* does NOT scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_UNBAL_EXPO,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_SUST_EXPO,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_TARGET_EXPO,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_BOND_FLOOR,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_BOND_FLOOR_LEVEL,
                        Results::INSTRUMENT_PACKET,
                        false  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_LOCKED_IN_VALUE,
                        Results::INSTRUMENT_PACKET,
                        true  /* does scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::SPI_REPORT_DATE,
                        Results::INSTRUMENT_PACKET,
                        false  /* does not scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        NONE  /* not statistical */),
    OutputRequestHelper(OutputRequest::SPI_REPORT,
                        Results::INSTRUMENT_PACKET,
                        false  /* does not scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        NONE  /* not statistical */),
    OutputRequestHelper(OutputRequest::SPI_GAP_EVENT_STATS,
                        Results::INSTRUMENT_PACKET,
                        false  /* does not scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        STATISTICAL  /* statistical sort of - needs some fiddling */),
    OutputRequestHelper(OutputRequest::SPI_LOAN_COST_RATE,
                        Results::INSTRUMENT_PACKET,
                        false  /* does not scale */,
                        false /* not additive */, false,
                        false /* don't add NA */,
                        NONE  /* not statistical */),
    // ------------- Tax Wrapper stuff ------------------
    OutputRequestHelper(OutputRequest::TAX_GROSS_UP,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */ , false /* additive */, false),
    OutputRequestHelper(OutputRequest::TAX_NEXT_CASHFLOW,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */ , true /* additive */, false),
    OutputRequestHelper(OutputRequest::ASSET_SWAP_SPREAD,
                        Results::INSTRUMENT_PACKET,
                        false /* doesnt scale */, true /* but adds */,
                        false /* doesn't scale pp*/),
    OutputRequestHelper(OutputRequest::BASIS_POINT_VALUE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesnt scale */, true /* but adds */,
                        false /* doesn't scale pp*/),
    OutputRequestHelper(OutputRequest::CDS_VALUE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesnt scale */, true /* but adds */,
                        false /* doesn't scale pp*/),
    OutputRequestHelper(OutputRequest::FORWARD_CDS_SPREAD,
                        Results::INSTRUMENT_PACKET,
                        false /* doesnt scale */, true /* but adds */,
                        false /* doesn't scale pp*/),
    OutputRequestHelper(OutputRequest::RISKY_STREAM_DETAILS,
                        Results::INSTRUMENT_PACKET,
                        false /* doesnt scale */, false /* but adds */,
                        false /* doesn't scale pp*/),
    OutputRequestHelper(OutputRequest::DDE_EQ_VOL,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_SPRD_BBONE,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_SPRD_VOL,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_EQ_SPRD_COR,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_SPRD_FUNC,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_LAST_GOOD_CALIB_DATE,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_CALIB_SCHEDULE,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_SPRD_DELTA,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_VOL_DELTA,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_VOL_SKEW,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_ATMVOL_DELTA,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_ATMVOL_SKEW,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::DDE_IMPLIED_VOL,
                        OutputRequest::DDE_PACKET,
                        false, false, false),
    OutputRequestHelper(OutputRequest::EDS_UPFRONT,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::EDS_RUNNING_FEE,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::EDS_RISKY_DURATION,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::EDS_DELTA_CURVE,
                        OutputRequest::EDS_DELTA_CURVE, // own packet
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::EDS_VEGA_CURVE,
                        OutputRequest::EDS_VEGA_CURVE, // own packet
                        false /* scales */, false /* and additive */, false),
	OutputRequestHelper(OutputRequest::FEE_PRICE,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::PL_AVERAGE,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::PL_DISTRIBUTION,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::PATHS,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::PL_DELTA_HEDGE,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::HEDGES,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::PRORATED_CASH_VALUE,
                        Results::INSTRUMENT_PACKET, false, false, true),
    OutputRequestHelper(OutputRequest::DBG,
                        OutputRequest::DBG, // own packet
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::PHYSICAL_DELIVERY,
                        OutputRequest::PHYSICAL_DELIVERY, // own packet
                        false /* doesn't scale */,
                        false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::COUPON_DUE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and additive */, false),
    OutputRequestHelper(OutputRequest::ACCRUAL_CALENDAR,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and additive */, false),
    OutputRequestHelper(OutputRequest::FIXED_LEG_VALUE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::LIBOR_LEG_VALUE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::LIBOR_FUNDING_TWEAK,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::TRANCHE_EXPECTED_LOSS_CURVE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_EXPECTED_LOSS_CURVE_WITH_RECOVERED_NOTIONAL,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_LOWER_BC_BETAS,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_UPPER_BC_BETAS,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_CONTINGENT_LEG_PRICE,
                        Results::TRANCHE_CONTINGENT_LEG_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_CONTINGENT_LEG_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_FEE_LEG_PRICE,
                        Results::TRANCHE_FEE_LEG_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_FEE_LEG_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_OUTSTANDING_NOTIONAL,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::REDEEMING_PROBABILITY,
                        Results::INSTRUMENT_PACKET,
                        false/* does not scale */, false /* and not additive */, false),
	OutputRequestHelper(OutputRequest::TRANCHE_ABS_BOTTOM_LOSS_CURVE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_ABS_TOP_LOSS_CURVE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_ABS_NOTIONAL_CURVE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::CPTY_CREDIT_CHARGE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* additive */, false),
    OutputRequestHelper(OutputRequest::FV_MINUS_CCC,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* additive */, false),
    OutputRequestHelper(OutputRequest::VALUE_WITHOUT_CCC,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* additive */, false),
    OutputRequestHelper(OutputRequest::PAR_SPREAD_CURVE,
                        OutputRequest::PAR_SPREAD_CURVE,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_IMPLIED_SPREAD,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::TRANCHE_RISKY_DURATION,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::BOND_FUTURE_Z_SPREAD,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */, false /* and not additive */, false),
    OutputRequestHelper(LegalTerms::createCalculator(), // note this is done at the top level
                        OutputRequest::LEGAL_TERMS_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* additive */,
                        false,
                        false /* don't add NA */,
                        STATISTICAL),
    OutputRequestHelper(OutputRequest::FEE_LEG_RISKLESS_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::CDS_RISKY_DURATION,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */,
                        false),
    OutputRequestHelper(OutputRequest::UNDERLYING,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::FWD_UNDERLYING,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::FLOAT_LEG_VALUE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::COMPLEX_LEG_VALUE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::MOD_CORR_MATRIX_SQ_ERROR,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */,
                        false),
    OutputRequestHelper(OutputRequest::CORR_REALIZED,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and not additive */,
                        false),
    OutputRequestHelper(OutputRequest::CORR_IMPLIED,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* and not additive */,
                        false),
    OutputRequestHelper(EventsMgr::createCalculator(), // note this is done at the top level
                        OutputRequest::EVENTS,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* additive */,
                        false,
                        false /* don't add NA */,
                        NONE),
    OutputRequestHelper(OutputRequest::FSA_VALUE,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* additive */,
                        false),
    OutputRequestHelper(OutputRequest::FSA_PRR,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* additive */,
                        false),
    OutputRequestHelper(OutputRequest::FIXED_LEG_DETAILS,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* additive */,
                        false),
    OutputRequestHelper(OutputRequest::FLOATING_LEG_DETAILS,
                        Results::INSTRUMENT_PACKET,
                        false /* scales */, false /* additive */,
                        false),
    OutputRequestHelper(OutputRequest::PVBP,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* additive */,
                        false),
    OutputRequestHelper(OutputRequest::INDEX_BASIS,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::EXPECTED_TIME_TO_EXERCISE,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::EXERCISE_PROBABILITY,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::CMCDS_ADJUSTED_RATES,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::CMCDS_FLOOR_PV,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::CMCDS_CAP_PV,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::CMCDS_DISCOUNT_FACTORS,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::CMCDS_SURV_PROBS,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::CMCDS_FORWARD_ADJS,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::CMCDS_MAT_ADJS,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::CMCDS_ACCS_ON_DEFAULT,
                        Results::INSTRUMENT_PACKET,
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::DURATION_WEIGHTED_AVG,
                        OutputRequest::DURATION_WEIGHTED_AVG, //own packet
                        false /* doesn't scale */,
                        false /* not additive */, 
                        false),
    OutputRequestHelper(OutputRequest::RISK_MAPPED_SENS),

    // for legacy FFX
    OutputRequestHelper(OutputRequest::FFX_UFV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::FFX_PV,
                        OutputRequest::FFX_PV, // own packet
                        true /* scales */, true /* and adds */,
                        false),
    OutputRequestHelper(OutputRequest::FWD_FX_RATE,
                        OutputRequest::FWD_FX_RATE, // own packet
                        false, false, false),
    OutputRequestHelper(OutputRequest::FFX_FUTURE_VALUE,
                        OutputRequest::FFX_FUTURE_VALUE, // own packet
                        true, true, false), // scales and adds
    OutputRequestHelper(OutputRequest::USD_DISC_RATE,
                        Results::INSTRUMENT_PACKET,
                        false, false, false),

    OutputRequestHelper(OutputRequest::DEBUG_RELEVANT_RISK_MAPPING_MATRICES),
    OutputRequestHelper(OutputRequest::IW_IND_VOL,
                        Results::INSTRUMENT_PACKET),
    OutputRequestHelper(OutputRequest::IW_STRIKE,
                        Results::INSTRUMENT_PACKET),
	OutputRequestHelper(OutputRequest::IW_INTEREST,
                        Results::INSTRUMENT_PACKET),
    OutputRequestHelper(OutputRequest::IW_PUT,
                        Results::INSTRUMENT_PACKET),
	OutputRequestHelper(OutputRequest::IW_CAPITAL,
                        Results::INSTRUMENT_PACKET),
    OutputRequestHelper(OutputRequest::DEBUG_RELEVANT_RISK_MAPPING_MATRICES),
	OutputRequestHelper(OutputRequest::IW_PVDIVS,
                        Results::INSTRUMENT_PACKET),
	OutputRequestHelper(OutputRequest::CONDITIONAL_LOSS_SAMPLES,
                        Results::INSTRUMENT_PACKET),
	OutputRequestHelper(OutputRequest::CONDITIONAL_LOSS_MAP,
                        Results::INSTRUMENT_PACKET),
	OutputRequestHelper(OutputRequest::FLAT_CDO,
                        Results::INSTRUMENT_PACKET),
	OutputRequestHelper(OutputRequest::FLAT_CDO_TRANCHE_DECOMPOSITION,
                        Results::INSTRUMENT_PACKET),
    OutputRequestHelper(OutputRequest::CONDITIONAL_FWD,
                        Results::INSTRUMENT_PACKET),
    OutputRequestHelper(OutputRequest::IW_FUTSTRIKE,
                        Results::INSTRUMENT_PACKET),
    /*Radar representation */
    OutputRequestHelper(OutputRequest::RADAR_REP,
                        OutputRequest::RADAR_PACKET, // own packet
                        false /* scales */, false /* and additive */, false),
    OutputRequestHelper(OutputRequest::RADAR_DIAGNOSTICS,
                        OutputRequest::RADAR_PACKET, // own packet
                        false /* scales */, false /* and additive */, false),
    /* Generic credit output requests */
    OutputRequestHelper(OutputRequest::CONTINGENT_LEG_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::FEE_LEG_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */, false),
    OutputRequestHelper(OutputRequest::ACCRUED_INT_FV,
                        Results::INSTRUMENT_PACKET,
                        true /* scales */, false /* and not additive */, false)

};

const int OutputRequestHelper::totalNumRequests =
sizeof(OutputRequestHelper::allRequests)/sizeof(OutputRequestHelper);

map<string, const OutputRequestHelper*> OutputRequestHelper::allRequestsMap;


CClassConstSP const OutputRequest::TYPE = CClass::registerClassLoadMethod(
    "OutputRequest", typeid(OutputRequest), OutputRequestHelper::load);

DEFINE_TEMPLATE_TYPE(OutputRequestArray);


// addin functions
class RequestListAddin: public CObject {
    static CClassConstSP const TYPE;

    /** the 'addin function' - list all sensitivities */
    static IObjectSP listRequest(RequestListAddin* params){
        static const string method = "RequestListAddin::listRequest";
        try {
            CStringArraySP request(OutputRequest::allRequests());
            return request;
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    RequestListAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(RequestListAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultRequestListAddin);

        Addin::registerClassObjectMethod("LIST_REQUEST",
                                         Addin::UTILITIES,
                                         "lists all requests",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)listRequest);
    }

    static IObject* defaultRequestListAddin(){
        return new RequestListAddin();
    }
};

CClassConstSP const RequestListAddin::TYPE = CClass::registerClassLoadMethod(
    "RequestListAddin", typeid(RequestListAddin), load);


DRLIB_END_NAMESPACE
