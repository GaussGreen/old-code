//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Description : General class that takes an ICreditLossConfig plus a fee
//               and contingent legs
//
// Date        : July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GeneralisedCDO.hpp"
#include "edginc/CCMPriceUtil.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/ICreditLossConfig.hpp"
#include "edginc/CreditTrancheLossConfig.hpp" // jlhp just temporarily
#include "edginc/CreditContingentLegBase.hpp"
#include "edginc/CDOTree.hpp"
#include "edginc/FtDProduct.hpp"
#include "edginc/CreditMetricsBaseCorrelation.hpp"   //
#include "edginc/BaseCorrelationOnlyParameters.hpp"  // Required for
#include "edginc/IndexWeights.hpp"                   // OutputRequest::
#include "edginc/CreditIndex.hpp"                    // DURATION_WEIGHTED_AVG
#include "edginc/NToDefaultLossConfig.hpp"
#include <map>
#include "edginc/ConditionalLossModel.hpp"
#include "edginc/ILossDistributionsGen.hpp"
#include "edginc/ArrayInstrumentCollection.hpp"
#include "edginc/Results.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/MCPathConfigCCM.hpp"
#include "edginc/FlatCDO2LossConfig.hpp"

DRLIB_BEGIN_NAMESPACE

typedef pair<ICDSParSpreadsWrapper, double> IndexEntry;
typedef vector<IndexEntry>                  IndexEntryArray;
typedef map<string, IndexEntryArray>        IndexInfoMap;


GeneralisedCDO::~GeneralisedCDO()
{}

/** Public constructor */
GeneralisedCDO::GeneralisedCDO(
        CClassConstSP                clazz,
        bool                         isLong,
        smartPtr<CounterPartyCredit> cptyInfo,
        ICreditContingentLegSP       cLeg,
        ICreditFeeLegSP              fLeg,
        CDOPortfolioSP               portfolio,
        CBoolSP                      recoverNotional,
        const string&                yieldCurveName,
        CIntSP                       triggerDelay,
        CIntSP                       defaultToCalculationDelay,
        DateTime                     lastTriggerDate,
        CDoubleSP                    temporaryLossAmount,
        const string&                settlementBDC,
        HolidayWrapper               settlementHols) :
    KComponent(yieldCurveName,"GeneralisedCDO",clazz),
    isLong(isLong),
    cptyInfo(cptyInfo),
    cLeg(cLeg),
    fLeg(fLeg),
    portfolio(portfolio),
    recoverNotional(recoverNotional),
    triggerDelay(triggerDelay),
    defaultToCalculationDelay(defaultToCalculationDelay),
    lastTriggerDate(lastTriggerDate),
    temporaryLossAmount(temporaryLossAmount),
    settlementBDC(settlementBDC),
    settlementHols(settlementHols),
	clProduct(false)
{
    // validatePop2Object() will be invoked from the constructor of
    // the concrete implementations of GeneralisedCDO that called
    // this method (eg, CDO)

}


void GeneralisedCDO::validatePop2Object() {
    static const string method("GeneralisedCDO::validatePop2Object");

    // Set default value if recoverNotional is not supplied.
    if (!recoverNotional) {
        recoverNotional.reset(CBool::create(false));
    }

    // If trigger delay is present, it must be non-negative
    if (!!triggerDelay && (triggerDelay->intValue() < 0)) {
        throw ModelException(method,
                             "The triggerDelay is negative (" +
                             Format::toString(triggerDelay->intValue()) +
                             "), and negative delays are not accepted.");
    }

    // If defaultToCalculationDelay is not present, set it to the triggerDelay
    // (which may or may not be present, too)
    if (!defaultToCalculationDelay) {
        defaultToCalculationDelay = triggerDelay;
    }
    else if (defaultToCalculationDelay->intValue() < 0) {
        // If defaultToCalculationDelay is present, it must be non-negative
        throw ModelException(method,
                             "The defaultToCalculationDelay is negative (" +
                             Format::toString(defaultToCalculationDelay->intValue()) +
                             "), and negative delays are not accepted.");
    }

    // If temporaryLossAmount is not present, assume 100%
    if (!temporaryLossAmount) {
        temporaryLossAmount.reset(CDouble::create(1.0));
    }
    else if ((temporaryLossAmount->doubleValue() < 0.0) ||
             (temporaryLossAmount->doubleValue() > 1.0))
    {
        throw ModelException(method,
                             "The temporaryLossAmount (" +
                             Format::toString(temporaryLossAmount->doubleValue()) +
                             ") is out of bounds [0, 1].");
    }
    settlementBadDayConv.reset(BadDayConventionFactory::make(settlementBDC));

	// populate settlement holidays
	if (settlementHols.isEmpty())
	{
		// default to no holidays
		settlementHols.setObject(MarketObjectSP(Holiday::noHolidays()));
	}

}


//// called after market data has been retrieved
void GeneralisedCDO::Validate() {
    static const string method("GeneralisedCDO::Validate");

    // lastTriggerDate is optional and may not be present. If so,
    // set it to the lastObservationDate
    if (lastTriggerDate.empty()) {
        lastTriggerDate = lastObservationDate();
    }

#if 0
    // MAR: removed this check since we do actually need to support this
    // There could still be lots of calculations to do
    if (port->nbName() > 0 && (port->nbName() == nbDefName))
        throw ModelException(method, "Portfolio contains only defaulted "
                             "names. Check this trade's booking.");
#endif
}


/** Returns the recover notional flag */
bool GeneralisedCDO::getRecoverNotional() const {
    return recoverNotional->boolValue();
}

//// Required part of CInstrument
DateTime GeneralisedCDO::getValueDate() const {
    return today;
}


//// Required part of Theta::Shift
bool GeneralisedCDO::sensShift(Theta* shift) {
    // alter immediate data
    today = shift->rollDate(today);
    return true; // then shift components
}

/** Retrieve market data from cache for all our components */
void GeneralisedCDO::GetMarket(const IModel* model, const CMarketDataSP market){

	static const string method = "GeneralisedCDO::GetMarket";

	if(FDModel::TYPE->isInstance(model))
	{
		// Call getMarket on underlying KComponent - needed for tree
		// for the moment only call this when the pricing model is of tree type
		// could remove this later if feel more confident with KComponent
		KComponent::GetMarket(model, market);

	}
	else
	{
		market->GetReferenceDate(today);

		const ConditionalLossModel* clModel =
			dynamic_cast<const ConditionalLossModel*>(model);

		if (clModel)
		{
			marketPointer = market;  // save for deferred call for second model

			GetMarket(clModel->getMCModel(), market);
			// or, can it be done right here?
			// GetMarket(clModel->getImpliedModel(), market);
		};

		portfolio->getMarket(model, market.get());
		if (cptyInfo.get()){
			cptyInfo->getMarket(model, market.get());
		}
		if (fLeg.get()){
			fLeg->getMarket(model, market.get());
		}

		if (discount.isEmpty())
		{
			// Note discount is now optional since it is member of KComponent
			throw ModelException(method, "discount curve is mandatory for "
                                 "GeneralisedCDO if not pricing on a tree");
		}
		discount.getData(model, market.get());

		settlementHols.getData(model, market);

	}

}


// Setup the component.
// The derived function must also call KComponent::setup() and the setup() of its underlyings
void GeneralisedCDO::setup(const IModel* model, const MarketData* market)
{
	static const string method = "GeneralisedCDO::setup";

	KComponent::setup(model, market);

	market->GetReferenceDate(today);

	portfolio->getMarket(model, market);
	if (cptyInfo.get()){
		cptyInfo->getMarket(model, market);
	}
	if (fLeg.get()){
		fLeg->getMarket(model, market);
	}


	// probably don't need to populate market wrappers here since done in KComponent::GetMarket
	discount.getData(model, market);

	settlementHols.getData(model, market);


}

//------------------------------------------
//  IBadDayAdjuster methods
//------------------------------------------
/** Returns "date" bad day adjusted using the bad day convention
 * and holidays in this object */
DateTime GeneralisedCDO::badDayAdjust(const DateTime& date) const {
    return settlementBadDayConv->adjust(date, settlementHols.get());
}

/** Add a number of business days to a date */
DateTime GeneralisedCDO::addBusinessDays(const DateTime& from, int busDays) const {
    return settlementHols->addBusinessDays(from, busDays);
}


/** when to stop tweaking (need to change infrastructure to route through
    model rather than instrument) */
DateTime GeneralisedCDO::endDate(const Sensitivity* sensControl) const {
    static const string method("GeneralisedCDO::endDate");

    // compute credit end date
    IModel* model = sensControl->getModel();
    if (ConvolutionEngine::TYPE->isInstance(model)) {
		const ConvolutionEngine* engine = DYNAMIC_CAST(ConvolutionEngine, model);
		return engine->endDate(this, sensControl);
	}
	else if (MonteCarlo::TYPE->isInstance(model)) {
		const MonteCarlo* mc = DYNAMIC_CAST(MonteCarlo, model);
		return mc->endDate(this, sensControl);
	}
	else if(SpreadLossTree::TYPE->isInstance(model)){
		return KComponent::endDate(sensControl);
	}
	else {
        throw ModelException(method, "End date only available for credit type "
                             "tweaks when a ConvolutionEngine, MonteCarlo or a SpreadLossTree"
                             "model has been specified");
    }
}


/** Returns the last 'scheduled' pay date (ie ignoring credit
    events).  Introduced for the output request IND_CDS_PAR_SPREAD
    - not clear if this is right date to use */
DateTime GeneralisedCDO::maturityDate() const {
    DateTime feeLegLastPayDate;
    DateTime ctgLegLastPayDate;

    if (fLeg.get()) {
        feeLegLastPayDate = fLeg->getLastPayDate();
    }
    if (cLeg.get()) {
        ctgLegLastPayDate =
            cLeg->lastPayDate(IBadDayAdjusterConstSP::attachToRef(this));
    }
    return feeLegLastPayDate.max(ctgLegLastPayDate);
}

/* return the max obs end date of the legs */
DateTime GeneralisedCDO::lastObservationDate() const {
    DateTime maxObsEndDate(today);

    if (cLeg.get()) {
        maxObsEndDate = maxObsEndDate.max(cLeg->lastObservationEndDate());
    }

    if (fLeg.get()) {
        DateTime maxObs = fLeg->getLastObservationDate();
        if (maxObsEndDate < maxObs){
            maxObsEndDate = maxObs;
        }
    }
    return maxObsEndDate;
}





CashFlowArraySP GeneralisedCDO::generateKnownCashflows(
    double               fLegOutstandingNotional,
    const CashFlowArray& fLegPastNotionalReductions,
    const CashFlowArray& cLegPastNotionalReductions,
    const BoolArray&     cLegPayPastNotionalReductions,
    IForwardRatePricerSP model) const
{
    static const string method("GeneralisedCDO::generateKnownCashflows");
    try {
        // initialise leg cashflows in case there arent any....
        // otherwise merge will cause a crash
        CashFlowArraySP fkcfl(new CashFlowArray(0));
        CashFlowArraySP ckcfl(new CashFlowArray(0));

        // fee leg known cashflows
        double initialNotional = portfolio->notional();

        if (fLeg.get()) {
            // get the cashflows - to change interface of generateKnownCashflows
            // method. Also take outstandingNotional
            DateTimeArray lossDates(
                CashFlow::dates(fLegPastNotionalReductions));
            DoubleArraySP lossAmounts(
                CashFlow::amounts(fLegPastNotionalReductions));
            double pastLoss = initialNotional - fLegOutstandingNotional;
            fkcfl = fLeg->generateKnownCashFlows(today,
                                                 initialNotional,
                                                 lossDates,
                                                 *lossAmounts,
                                                 pastLoss,
                                                 model);

            // apply the direction setting
            for (int i = 0; i < fkcfl->size(); i++) {
                (*fkcfl)[i].amount = - (*fkcfl)[i].amount;
            }
        }

        // contingent leg known cashflows
        if (cLeg.get()) {
            // get the cashflows
            ckcfl = cLeg->generateKnownCashFlows(today,
                                                 initialNotional,
                                                 cLegPastNotionalReductions,
                                                 cLegPayPastNotionalReductions,
                                                 IBadDayAdjusterConstSP::attachToRef(this));
        }

        // now merge the 2 lists together
        CashFlowArraySP kcfl(CashFlow::merge(fkcfl,ckcfl));
        // and combine coincidental dates
        CashFlow::aggregate(*kcfl);
        // and finish
        return kcfl;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


GeneralisedCDO::GeneralisedCDO(CClassConstSP clazz) :
	KComponent("","GeneralisedCDO",clazz),
    cptyInfo(0),
    cLeg(0),
    fLeg(0),
    settlementBDC("N"),
	clProduct(false)
{}

IObject* GeneralisedCDO::defaultConstructor() {
    return new GeneralisedCDO();
}

/** Invoked once at start up when this class is 'loaded' */
void GeneralisedCDO::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(GeneralisedCDO, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(ConvolutionEngine::IIntoProduct);
	IMPLEMENTS(FDModel::IIntoProduct);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(IBadDayAdjuster);
    IMPLEMENTS(IRebateCalculator);
    IMPLEMENTS(IProtectionProvider);
	IMPLEMENTS(IMCIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(today, "value date");
    FIELD_MAKE_OPTIONAL(today);
    FIELD(isLong, "Long or Short");
    FIELD(cptyInfo, "cpty charge info");
    FIELD_MAKE_OPTIONAL(cptyInfo);
    FIELD(cLeg, "contingent leg");
    FIELD_MAKE_OPTIONAL(cLeg);
    FIELD(fLeg, "Fee Leg");
    FIELD_MAKE_OPTIONAL(fLeg);
    FIELD(portfolio, "Generalised loss config");
	// Note that discount is now an optional field in KComponent
	FIELD(recoverNotional, "Recover notional on the fee leg? "
          "default false");
	FIELD_MAKE_OPTIONAL(recoverNotional);

    FIELD(triggerDelay, "Delay between market default and "
          "eventDeterminationDate");
    FIELD_MAKE_OPTIONAL(triggerDelay);
    FIELD(defaultToCalculationDelay, "Delay between market "
          "default and the associated calculation date");
    FIELD_MAKE_OPTIONAL(defaultToCalculationDelay);
    FIELD(lastTriggerDate, "Last date when a default occurred during the "
          "protection period can be triggered. Default: Protection end date.");
    FIELD_MAKE_OPTIONAL(lastTriggerDate);
    FIELD(temporaryLossAmount, "Assumption for the loss amount associated to a "
          "default, used to determine the fees in accrual periods between the "
          "credit event and the calculation date. Default: 100%");
    FIELD_MAKE_OPTIONAL(temporaryLossAmount);
    FIELD(settlementBDC, "Bad day convention. Used to adjuste the delays "
          "regarding defaults' settlements. Default: None");
    FIELD_MAKE_OPTIONAL(settlementBDC); // to be made mandatory ideally
    FIELD(settlementHols, "Holidays. Used to adjust the delays regarding "
          "defaults' settlements. Default: weekends only");
    FIELD_MAKE_OPTIONAL(settlementHols); // to be made mandatory ideally
    FIELD(settlementBadDayConv, "Bad day convention for settlements");
    FIELD_MAKE_TRANSIENT(settlementBadDayConv);

	FIELD(clProduct, "Is this a conditional Loss product. default: false");
    FIELD_MAKE_TRANSIENT(clProduct);
}

CClassConstSP const GeneralisedCDO::TYPE = CClass::registerClassLoadMethod(
    "GeneralisedCDO", typeid(GeneralisedCDO), load);


/*
 * A few remarks on tranche expected loss calls organization:
 *
 * To compute the tranche expected loss adjusted by the london floor method,
 * it is convenient to write the tranche expected loss calculation function
 * like a parameterised function call: F(param, K1, K2)
 * (where param is a parameter that does not change when strikes change).
 * The london floor methodology then consists in making up to 3 calls to F:
 * F(param, K1, Keq), F(param, Keq, Kse), F(param, Kse, K2)
 *
 * For the correlation skew methodology, it is convenient to introduce
 * a parameterized function G, which takes the beta vector as a
 * parameter and return the expected loss of an equity tranche with upper
 * attachment point K: G(param, beta, K)
 * (where param is a parameter that does not change when strike/beta change).
 * To price a [K1,K2] tranche, we then need 2 calls of G:
 * G(param, beta1, K1), G(param, beta2, K2)
 *
 * In terms of implementation of F or G, these functions can be implemented
 * either using the convolution algorithm or the fast loss calculation
 * algorithm.
 *
 * Optimization 1: for the london floor implementation using the convolution
 * algorithm, the 3 calls to F can be made based on the same loss density.
 * It is therefore much faster to calculate the loss density once, store
 * it in lambda, and then do the multiple calls to F.
 *
 * Optimization 2: tranche expected loss conditional on counterparty surviving
 * can be done much faster if it is done at the same time as the unconditional
 * one. Because of that, the function F and G return both unconditional and
 * cpty survival conditional expected loss (and both the conditional and
 * unconditional loss density are stored in lambda for the convolution algorithm
 * call).
 *
 * - Function and Structure Names used for London Floor
 *   -# convolution: convTrancheExpectedLoss,     ConvParam
 *   -# fast loss:   fastTrancheExpectedLoss,     FastParam
 * - Function and Structure Names used for beta skew: to be defined
 */



//------------------------------------------
//  IRebateCalculator method
//------------------------------------------
/** Compute the fee rebate payment due to the difference between the temporary
 * loss assumption and the real losses on the fee leg.
 * Note the payment will not be pv'd in the sense that, if the fees should
 * have been X dollars higher on the last fee payment date, the rebate will be
 * X dollars - the fact that the rebate will be paid several days/weeks after
 * the original fee payment date does not impact the rebate computation.
 * Actually this method is not at all aware of when the actual rebate payment
 * will take place */
double GeneralisedCDO::computeRebate(
    const CashFlowArrayConstSP assumedReductions,
    const CashFlowArrayConstSP realReductions,
    const double initialNotional,
    IForwardRatePricerSP model) const
{
    static const string method("GeneralisedCDO::computeRebate");
    try {
        double rebate = 0.0;

        if (!!fLeg) {
            // get the cashflows with the assumed reductions
            DateTimeArray assumedReductionDates(CashFlow::dates(*assumedReductions));
            DoubleArraySP assumedReductionAmounts(CashFlow::amounts(*assumedReductions));
            CashFlowArraySP assumedKcfl =
                fLeg->estimateKnownCashFlows(initialNotional,
                                             assumedReductionDates,
                                             *assumedReductionAmounts,
                                             model);

            // get the cashflows with the real reductions
            DateTimeArray realReductionDates(CashFlow::dates(*realReductions));
            DoubleArraySP realReductionAmounts(CashFlow::amounts(*realReductions));
            CashFlowArraySP realKcfl =
                fLeg->estimateKnownCashFlows(initialNotional,
                                             realReductionDates,
                                             *realReductionAmounts,
                                             model);


            // The rebate is the difference between the sum of the realKcfl
            // amounts minus the sum of the assumedKcfl amounts.
            // Note that the cashflows above have the oposite signs to what
            // they should, so perform the oposite subtraction
            int numKcfl = assumedKcfl->size();
            if (numKcfl != realKcfl->size()) {
                throw ModelException(method,
                                     "Internal error: The arrays of 'assumed' "
                                     "and 'real' losses have different lenghts: " +
                                     Format::toString(numKcfl) + " vs " +
                                     Format::toString(realKcfl->size()));
            }
            for (int i=0; i<numKcfl; ++i) {
                rebate += (*assumedKcfl)[i].amount - (*realKcfl)[i].amount;
            }
        }
        return rebate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns the max of maxDate and the last date from when a yield
    curve is used (eg for discounting or rate estimation).
    This method does not have to worry about what dates are used
    when the clean spread curves are built: Just used to
    control when to stop tweaking */
DateTime GeneralisedCDO::lastYCSensDate(const DateTime& maxDate) const {
    DateTime newMax(maxDate);
    if (cLeg.get()){
        newMax = cLeg->lastYCSensDate(maxDate, IBadDayAdjusterConstSP::attachToRef(this));
    }
    if (fLeg.get()){
        newMax = fLeg->lastYCSensDate(newMax);
    }
    return newMax;
}

typedef pair<double, double> doublePair;

typedef vector<doublePair> doublePairArray;

typedef vector<doublePairArray> doublePairArrayArray;

//static DateTime LARGEDATE = DateTime("01-Jul-3000", "EOD");
//static DateTime SMALLDATE = DateTime("01-Jul-1900", "EOD");

static const DateTime LARGEDATE = DateTime(10000000,0);
static const DateTime SMALLDATE = DateTime(100, 0);

static SimSeriesSP trivialSS(const DateTimeArray &dates)
{
	SimSeriesSP ss = SimSeriesSP(new SimSeries(1));

	ss->addDates(dates);
	return ss;
}

static AccrualPeriodArraySP convertDates(
	const DateTimeArray & dates)
{
	AccrualPeriodArraySP result = AccrualPeriodArraySP (
					new AccrualPeriodArray());
	int i;
	for (i = 0; i < dates.size()-1; i++)
	{
		result->push_back(AccrualPeriodSP(new AccrualPeriod(dates[i], dates[i+1])));
	};

	return result;
};

static void addCutOffs(list<DateTime> &timelineList,
					   const ICreditLossConfig *lossConfig)
{

	const CreditTrancheLossConfig *tranche;
	const CDOPortfolio *portfolio;
	const PortfolioName *name;

	if (tranche = dynamic_cast<const CreditTrancheLossConfig*>(lossConfig))
	{
		if (tranche->hasCutOffFlag()) {
			timelineList.push_back(tranche->getCutOffDate());
		};
		timelineList.push_back(tranche->getProtectionStartDate());

		addCutOffs(timelineList, tranche->getPortfolio().get());

	} else if (portfolio = dynamic_cast<const CDOPortfolio*>(lossConfig))
	{
		int i;
		for (i = 0; i < portfolio->numInnerLossConfigs(); i++)
		{
			addCutOffs(timelineList, portfolio->getInnerLossConfig(i).get());
		};
	} else if (name = dynamic_cast<const PortfolioName*>(lossConfig))
	{
		const DateTime& matCutOff = name->getNameMaturityCutOff();
		if (!matCutOff.empty())
			timelineList.push_back(matCutOff);
		timelineList.push_back(name->getProtectionStartDate());

	} else {

		throw ModelException("Unknown lossConfig type");
	};
};

static void
sortAndRemoveDuplicates(
	DateTimeArray& dateTimeArray, //output
	list<DateTime>& dateTimeList)
{
	//then sort the dates
	dateTimeList.sort();

	//remove the duplicates
	if (dateTimeList.size() > 1)
	{
		list<DateTime>::iterator iter = dateTimeList.begin();
		list<DateTime>::iterator lastIter = dateTimeList.begin();
		iter++;

		while (iter!= dateTimeList.end())
		{
			if (*iter == *lastIter) //duplicate
			{
				dateTimeList.erase(iter);
				iter = lastIter;
			}
			if (iter!= dateTimeList.end())
			{
				lastIter  = iter;
				iter++;
			}
		}
	}

	//now covert to dateTimeArray
	//what if it only has the value date inside it?
	dateTimeArray.clear();

	for (list<DateTime>::iterator iter = dateTimeList.begin(); iter!= dateTimeList.end(); ++iter)
		dateTimeArray.push_back(*iter);
}

class PRODUCTS_DLL GeneralisedCDO::ConditionalLossProduct :
	public virtual ConditionalLossModel::IProduct,
	public virtual VirtualDestructorBase
{
	const CInstrument *instrument;

	const ConditionalLossModel *model;

	CDOPortfolioSP regressor;

	DateTimeArraySP datesOfInterest;

public:

	virtual ~ConditionalLossProduct() {};

	ConditionalLossProduct(
		const GeneralisedCDO *inst,
		const ConditionalLossModel *mod
	) :
	instrument( inst ),
	model(mod)
	{};

public:

	virtual void price(Control* control, CResults* results);

	static CDOPortfolioSP getFlatLossConfig(const CreditTrancheLossConfig* lc);

}; // GeneralisedCDO::ConditionalLossProduct


class PRODUCTS_DLL GeneralisedCDO::ConditionalLossMCProduct :
	public virtual MCProductClient,
	public virtual VirtualDestructorBase
{
public:

	// need a CObject to use the Results infrastructure
	// number of rows = size of column corresponds to number of samples
	// number of columns = size of row corresponds to number of timepoints per sample
	class ResultType : public DoubleMatrix {

		DoubleMatrix values;

	public:

		static CClassConstSP const TYPE;

		ResultType(CClassConstSP clazz = TYPE) {};

		static IObject* defaultConstructor()
		{
			return new ResultType();
		};

		// DoubleMatrix takes #cols, #rows ...
		ResultType(
			const doublePairArrayArray &x,
			CClassConstSP clazz = TYPE
		) :
			DoubleMatrix( x[0].size(), x.size()),
			values( x[0].size(), x.size() )
		{
			int i, j;

			for (i = 0; i < numCols(); i++) // loop over columns = timepoints
				for (j = 0; j < numRows(); j++)   // loop over rows = samples
				{
					(*this)[i][j] = x[j][i].first;
					values[i][j] = x[j][i].second;
				}
		}

		// i = which column = which timepoint
		// j = which row = sample
		double X(int i, int j) const { return (*this)[i][j];};

		double Y(int i, int j) const { return values[i][j];};

		static void load(CClassSP& clazz);

		virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const
		{
			DoubleMatrix::outputWrite(linePrefix, prefix+"::Xvalues", stream);

			values.outputWrite(linePrefix, prefix+"::Yvalues", stream);
		};
	};

	DECLARE(ResultType);

	class Prices : public IMCPricesSimple,
		  		   public virtual VirtualDestructorBase
	{

		doublePairArrayArray losses;


	public:

		Prices(){};

		virtual ~Prices() {};

		virtual ResultTypeSP getResults() const {return ResultTypeSP(new ResultType(losses));};

		virtual void add(const doublePairArray &x);

		virtual void postResults(
			MonteCarlo*                mcarlo,  // the model
			IMCProduct*                 product,
			Control*                   control,
			IMCPrices&         prices,
			const string&              discountISOCode,
			Results*                   results) {};

		virtual void add(double x)
		{
			throw ModelException("add(double) not defined for ",
				"GeneralisedCDO::ConditionalLossMCProduct::Prices");
		};

		virtual void getResult(double *, double*) const {};

		virtual bool repriceForGreek(int) {return false;};

		virtual void reset() {};

		virtual IMCPrices* clone() const {return NULL;};

		virtual double lastPrice() const {return 0;};

		virtual int storagePerPath(IMCProduct*) const {return 2*sizeof(double);};

		virtual void configureCache(const IntArray&) {};

		virtual IMCPrices* emptyConstructor()const {return new Prices();};

		virtual double maxWithZero(double x) const {return (x>0)?x:0;};


	}; // GeneralisedCDO::ConditionalLossMCProduct::Prices

	DECLARE(Prices);

private:

	// these all refer to CreditTrancheLossConfigs
	const ICreditLossConfigSP instLC;
	const ICreditLossConfigSP regrLC;

	ICreditLossConfigIndexedSVGenConstSP instGen;
	ICreditLossConfigIndexedSVGenConstSP regrGen;

	ICreditLossConfigIndexedSVSP indexedInstSV;
	ICreditLossConfigIndexedSVSP indexedRegrSV;

	const GeneralisedCDO *instrument;

public:

	// this method gets SV's from the gen, and adds them to the collector
	// MC calls it
	virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

	// adds the values of both SV's to "prices", the array of results
	// somehow, the instance of prices knows to store the results (an array of doubles)
	// instead of averaging them
	virtual void payoff(
		const MCPathGenerator*  pathGen,
		IMCPrices& prices);

	virtual IMCPrices* createOrigPrices(
		int  nbIter,
        int  nbSubSamples,
        int  mode);

	virtual void recordExtraOutput(
	CControl*     control,
    Results*      results,
    const IMCPrices& prices) const;

	virtual ~ConditionalLossMCProduct() {};

	void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

	ConditionalLossMCProduct(
		const GeneralisedCDO *inst,
		const MonteCarlo *mod
	) :
	instrument(inst),
	MCProductClient(
		inst->getValueDate(),
		inst->discount.get(),
		trivialSS(*inst->auxTimeLine.get())
	),
	instLC( inst->portfolio ),
	regrLC( GeneralisedCDO::ConditionalLossProduct::getFlatLossConfig(
		dynamic_cast<const CreditTrancheLossConfig*>(inst->portfolio.get()) )
	)
	{
		instLC->getMarket(mod, inst->marketPointer.get());
		regrLC->getMarket(mod, inst->marketPointer.get());

		// create the SV gens
		CIntConstSP triggerDelay;
		CIntConstSP defaultToCalculationDelay;
		double temporaryLossAmount = 0;    // JCP
		DateTime lastTriggerDate;
		AccrualPeriodArrayConstSP accrualPeriods = convertDates(*inst->auxTimeLine.get());
		IBadDayAdjusterConstSP bda;
		IProtectionProviderConstSP protect;
		IRebateCalculatorConstSP rebateCalc;

		CreditTrancheLossConfig* tempInstLC =
			dynamic_cast<CreditTrancheLossConfig*>(instLC.get());

		if (!tempInstLC)
			throw ModelException("Not a Credit Tranche for CL method");

//	add tranchCutOff to the auxTimeline 		
		list<DateTime> dateTimeList; 
		for (int i=0; i<inst->auxTimeLine->size(); ++i)
			dateTimeList.push_back((*inst->auxTimeLine)[i]);
		
		addCutOffs(dateTimeList, inst->portfolio.get());
		
		DateTimeArraySP timeline(new DateTimeArray()); 
		sortAndRemoveDuplicates(
			*timeline, 
			dateTimeList);

// create Indexed SV Gen
		instGen = tempInstLC->createIndexedSVGen(
		     timeline,
		     triggerDelay,
		     defaultToCalculationDelay,
			 temporaryLossAmount,
		     lastTriggerDate,
		     accrualPeriods,
		     bda,
		     protect,
		     rebateCalc,
			 true // JCP recover notional
		);

		CDOPortfolio* tempRegrLC =
			dynamic_cast<CDOPortfolio*>(regrLC.get());

		if (!tempRegrLC)
			throw ModelException("Not a CDPPortfolio for CL method");

		regrGen = tempRegrLC->createIndexedSVGen(
		     timeline,
		     triggerDelay,
		     defaultToCalculationDelay,
			 temporaryLossAmount,
		     lastTriggerDate,
		     accrualPeriods,
		     bda,
		     protect,
		     rebateCalc,
			 true // JCP recover notional
		);

		indexedInstSV = instGen->createNewSV(NULL);
		indexedRegrSV = regrGen->createNewSV(NULL);
	};
};

void
GeneralisedCDO::ConditionalLossMCProduct::ResultType::load(CClassSP& clazz)
{
	REGISTER(ResultType, clazz);

	SUPERCLASS(DoubleMatrix);

	clazz->setPublic();

	EMPTY_SHELL_METHOD(defaultConstructor);

	FIELD(values, "Y values");

}

CClassConstSP const GeneralisedCDO::ConditionalLossMCProduct::ResultType::TYPE
		= CClass::registerClassLoadMethod(
			"GeneralisedCDO::ConditionalLossMCProduct::ResultType",
			typeid(GeneralisedCDO::ConditionalLossMCProduct::ResultType),
			load);


void
GeneralisedCDO::ConditionalLossMCProduct::Prices::add(
	const doublePairArray &x)
{
	losses.push_back(x);
}

class NameTableEntry {

public:

	string name;
	double notional;
	DateTime cutoffDate;
	DateTime startDate;
	const PortfolioName* pname;

	NameTableEntry(string nm):
		name(nm), notional(0),  cutoffDate(SMALLDATE), startDate(LARGEDATE), pname(NULL)
	{}
};

typedef  vector<NameTableEntry> NameTableArray;

int find(string name, NameTableArray &table)
{
	unsigned int i;

	for (i = 0; i < table.size(); i++)
	{
		if (table[i].name == name) return i;
	};

	table.push_back(NameTableEntry(name));

	return table.size()-1;
};

// insert/update name in flat representation
void handleName(const PortfolioName* pname,
				double amount,
				DateTime refDate,
				NameTableArray& table)
{
	int pos = find(pname->getName(), table);

	table[pos].notional += amount;

	if (!table[pos].cutoffDate.empty())
	{
		const DateTime& date = pname->getNameMaturityCutOff();	
		if (date.empty())
			table[pos].cutoffDate = date; 
		else 
			table[pos].cutoffDate = max(date, 
										table[pos].cutoffDate);
	}
	table[pos].startDate = min(table[pos].startDate,
							   pname->getProtectionStartDate());

	table[pos].pname = pname;
}

CDOPortfolioSP GeneralisedCDO::ConditionalLossProduct::getFlatLossConfig(
	const CreditTrancheLossConfig* lc)
{
	DateTime refDate = LARGEDATE; // ctlc.maturityDate();

	if (!lc)
		throw ModelException(
			"Need a CreditTrancheLossConfig to make a flat CDO");

	const CDOPortfolio* lc1 =
		dynamic_cast<const CDOPortfolio*>(lc->getPortfolio().get());

	if (!lc1)
		throw ModelException(
			"Odd: CreditTrancheLossConfig should contain a CDOPortfolio");

	NameTableArray names;

	int i;

	for (i = 0; i < lc1->numInnerLossConfigs(); i++)
	{
		const PortfolioName* pname;

		// components are either single names or CDO's
		if ( pname = dynamic_cast<const PortfolioName*>
							( lc1->getInnerLossConfig(i).get() ) )
		{
			handleName(
				pname,
				pname->getNameNotional(),
				refDate,
				names);

		} else if ( dynamic_cast<const CreditTrancheLossConfig*>
							( lc1->getInnerLossConfig(i).get() ) )
		{
			const CreditTrancheLossConfig* ctl =
				dynamic_cast<const CreditTrancheLossConfig*>
					( lc1->getInnerLossConfig(i).get() );

            double ntnl = ctl->trancheOutstandingNotional(false);

			// strike diff / total notional
			double factor = ntnl / ctl->portfolioNotional();

			const CDOPortfolio* lc2 =
				dynamic_cast<const CDOPortfolio*>(ctl->getPortfolio().get());

			if (!lc2)
				throw ModelException(
					string("Odd: CreditTrancheLossConfig ") +
					string("should contain a CDOPortfolio"));

			int j;
			for (j = 0; j < lc2->numInnerLossConfigs(); j++)
			{
				pname = dynamic_cast<const PortfolioName*>
							( lc2->getInnerLossConfig(j).get() );

				if (!pname)
					throw ModelException("Not a CDO^2");

				handleName( pname,
							pname->getNameNotional(),  // removed factor
							refDate,
							names);
			}
		} else {
			throw ModelException("Not a CDO^2");
		}
	}

	PortfolioNameArray pNames;

	for (i = 0; i < (int) names.size(); i++)
	{
		PortfolioName *copy =
			dynamic_cast<PortfolioName*>( names[i].pname->clone() );

		copy->AdjustNotionalAndDates(
			names[i].notional,
			names[i].startDate,
			names[i].cutoffDate);

		pNames.push_back(PortfolioNameSP(copy));

	};

	return CDOPortfolioSP(new CDOPortfolio(pNames));

};

IMCPrices*
GeneralisedCDO::ConditionalLossMCProduct::createOrigPrices(
		int  nbIter,
        int  nbSubSamples,
        int  mode)
{
	return new ConditionalLossMCProduct::Prices();
}

// returns a conditional loss product, like any other model
ConditionalLossModel::IProduct* GeneralisedCDO::createProduct(
	const ConditionalLossModel* model) const
{
	return new ConditionalLossProduct(this, model);
};


void
GeneralisedCDO::ConditionalLossMCProduct::pathGenUpdated(
	IStateVariableGen::IStateGen* newPathGen)
{
	static const string method = "ConditionalLossMCProduct::pathGenUpdated";
	try
	{
		if(newPathGen)
		{
			// fetching the credit loss config SV
			IStateVariableSP stateVariableInst =
				instGen->create(
					indexedInstSV,
					newPathGen);
			// this will return the object from the svDBase

			indexedInstSV =
				ICreditLossConfigIndexedSVSP(
					dynamic_cast<ICreditLossConfig::IIndexedSVGen::ISV*>
						(stateVariableInst.get()));

			// fetching the regressor SV
			IStateVariableSP stateVariableRegr =
				regrGen->create(
					indexedRegrSV,
					newPathGen);//will return the object from the svDBase

			indexedRegrSV =
				ICreditLossConfigIndexedSVSP(
					dynamic_cast<ICreditLossConfig::IIndexedSVGen::ISV*>
						(stateVariableRegr.get()));

		}
    }
	catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

void
GeneralisedCDO::ConditionalLossMCProduct::payoff(
		const MCPathGenerator*  pathGen,
		IMCPrices& prices)
{
	// prices is of type ConditionalLossMCProduct::Prices

	Prices *localPrices = dynamic_cast<Prices *> (&prices);

	const CreditLossConfigIndexedSVResultPath& resultsInst =
		indexedInstSV->getResults();

	const CreditLossConfigIndexedSVResultPath& resultsRegr =
		indexedRegrSV->getResults();

	doublePairArray losses;

	int i;
	double totalLossInst = 0;
	double totalLossRegr = 0;

	// over the conditional loss timeline
	DoubleArray lossInst(instrument->auxTimeLine->size());
	DoubleArray lossRegr(instrument->auxTimeLine->size());

	for (i = resultsInst.begin(); i <= resultsInst.end(); ++i)
	{
		if (-1 != resultsInst[i].index)
			lossInst[resultsInst[i].index] += resultsInst[i].lossChange;

		if (-1 != resultsRegr[i].index)
			lossRegr[resultsRegr[i].index] += resultsRegr[i].lossChange;
	};

	for(i = 1; i < instrument->auxTimeLine->size(); i++)
	{
		lossRegr[i] += lossRegr[i-1];

		lossInst[i] += lossInst[i-1];

		losses.push_back(doublePair(lossRegr[i], lossInst[i]));
	};

	localPrices->add(losses);
}

static
IntArraySP countNonTrivial(
	const GeneralisedCDO::ConditionalLossMCProduct::ResultType* lresults,
	IMFSamplingSP sampling)
{

	// determine min and max
	double minLoss, maxLoss;

	int i, last;

	// numbCols = # of time points, we use the last one
	last = lresults->numCols() - 1;

	for (i = 0; i < lresults->numRows(); i++)
	{
		if (i>0)
		{
			maxLoss = max(lresults->Y(last,i), maxLoss);
			minLoss = min(lresults->Y(last,i), minLoss);
		} else {
			minLoss = maxLoss = lresults->Y(last,i);
		};
	};

	DoubleArrayConstSP points  = sampling->getPoints();
	IntArrayConstSP numSamples = sampling->getCumNumberSamples();

	IntArraySP  result = IntArraySP(new IntArray(points->size()));

	int pathIdx = 0;
	int pointIndex = 0;
	for (pathIdx = 0; pathIdx < lresults->numRows(); pathIdx++)
	{
		while (pathIdx >= (*numSamples)[pointIndex])
		{
			pointIndex++;
			if (pointIndex >= numSamples->size())
			{
				// something's gone wrong
				throw ModelException("Drawing too many paths");
			};
		};

		if ( (lresults->Y(last, pathIdx) < maxLoss) && (lresults->Y(last, pathIdx) > minLoss) )
		{
			(*result)[pointIndex]++;
		}
	};

	return result;
}


static
FlatCDO2LossConfig::LinearLossProfileArraySP
processSamples(
	const GeneralisedCDO::ConditionalLossMCProduct::ResultType* lresults,
	int nBuckets)
{
	int i = 0;
	int j = 0;

	FlatCDO2LossConfig::LinearLossProfileArraySP slices =
		FlatCDO2LossConfig::LinearLossProfileArraySP(
				new FlatCDO2LossConfig::LinearLossProfileArray());

	// number of timepoints
	for (j = 0; j < lresults->numCols(); j++)
	{
		DoubleArray xs;
		DoubleArray ys;
		DoubleArray os;
		DoubleArray lowBounds;
		DoubleArray upBounds;

		DoubleArray avgx(nBuckets, 0);
		DoubleArray avgy(nBuckets, 0);
		DoubleArray avgy2(nBuckets, 0);

		DoubleArray bounds(nBuckets+1);

		IntArray count(nBuckets, 0);

		double maxLoss, minLoss;

		for (i = 0; i < lresults->numRows(); i++)
		{
			if (i>0)
			{
				maxLoss = max(lresults->X(j,i), maxLoss);
				minLoss = min(lresults->X(j,i), minLoss);
			} else {
				minLoss = maxLoss = lresults->X(j,i);
			};
		};

		double dx = (maxLoss - minLoss)/nBuckets;

		bounds[0] = minLoss;

		for (i = 1; i <= nBuckets; i++)
			bounds[i] = bounds[i-1] + dx;

		bounds[0] -= 1e6;
		bounds[nBuckets] += 1e6;
		// make limits wider to avoid comparison/roundoff errors

		// number of samples
		for (i = 0; i < lresults->numRows(); i++)
		{
			int  k = 0;

			// could be done faster with binsearch

			while ( bounds[k+1] < lresults->X(j,i) )
			{
				k++;

				// unlikely
				if (k >= nBuckets)
					break;
			};

			// unlikely
			if (k >= nBuckets)
					break;

			avgx[k] += (lresults->X(j,i) );

			avgy[k] += (lresults->Y(j,i) );

			avgy2[k] += (lresults->Y(j,i) )*(lresults->Y(j,i) );

			count[k]++;
		}

		for (i = 0; i < nBuckets; i++)
		{
			if (count[i] > 0)
			{
				avgx[i] /= count[i];
				avgy[i] /= count[i];
				avgy2[i] /= count[i];
				avgy2[i] = sqrt(avgy2[i] - avgy[i]*avgy[i]);
			};
		};

		// create profiles with the buckets with +ve number of hits
		for (i = 0; i < nBuckets; i++)
		{
			if (count[i] > 0)
			{
				xs.push_back(avgx[i]);
				ys.push_back(avgy[i]);
				os.push_back(avgy2[i]);

				if (i> 0)
					lowBounds.push_back( bounds[i] );
				else
					lowBounds.push_back( minLoss );

				if (i<nBuckets-1)
					upBounds.push_back( bounds[i+1] );
				else
					upBounds.push_back( maxLoss );
			}
		};

		slices->push_back(FlatCDO2LossConfig::LinearLossProfileSP(
			new FlatCDO2LossConfig::LinearLossProfile(xs, ys, lowBounds, upBounds, os)));
	};

	return slices;
}

void
GeneralisedCDO::ConditionalLossMCProduct::recordExtraOutput(
	CControl*     control,
    Results*      results,
    const IMCPrices& prices) const
{
	OutputNameConstSP outputName;

	const Prices* localPrices = dynamic_cast<const Prices*>(&prices);

	results->storeRequestResult(
		new OutputRequest(OutputRequest::CONDITIONAL_LOSS_SAMPLES),
        localPrices->getResults()
	);
};

void GeneralisedCDO::ConditionalLossProduct::price(
    Control*      control,
    CResults*     results)
{

	static const string method =
		"GeneralisedCDO::ConditionalLossProduct::price";

	Results localResults;

	// mcModel will create a ConditionalLossProduct
	// which knows about the two instruments

	// where MC calls pathUpdate and payoff, the results type has to be such
	// that it takes in pairs for add().
	// this is achieved by getting ????? (in MC) to use ConditionalLossPricing
	// probably this is what localControl needs to indicate to MonteCarlo

	// 1. create a conditional loss product
	// 2. call it's price method

	// the instrument member is const, so we need to do this

	try {

			Instrument *localInstrument = const_cast<Instrument *>(instrument);

		GeneralisedCDO* gcdo = dynamic_cast<GeneralisedCDO*> (localInstrument);

		const CreditTrancheLossConfig* config =
			dynamic_cast<const CreditTrancheLossConfig*>(gcdo->portfolio.get());

		if (!config)
			throw ModelException(
				"Expected a CreditTrancheLossConfig",
				"GeneralisedCDO::ConditionalLossProduct::price");

		CDOPortfolioSP flatPfolio = getFlatLossConfig(config);

		FlatCDO2LossConfig::LinearLossProfileArraySP slices;

		GeneralisedCDO::ConditionalLossMCProduct::ResultType* lresults;

		if (control->isPricing())
		{

			SensitivityArraySP sens =  SensitivityArraySP(new SensitivityArray());

			OutputRequestSP clRequest = OutputRequestSP(
				new OutputRequest(OutputRequest::CONDITIONAL_LOSS_SAMPLES));

			OutputRequestArraySP outReqs =
				OutputRequestArraySP(new OutputRequestArray());

			outReqs ->push_back(clRequest);

			Control localControl(
				sens,
				outReqs,
				false,
				"");

			gcdo->clProduct = true;

			gcdo->auxTimeLine = model->getTimeLine();

			IMCPathConfigSP ipc = model->getMCModel()->getPathConfig();

			MCPathConfigCCM* pcccm = dynamic_cast<MCPathConfigCCM*>(ipc.get());

			if (pcccm) {

				IMFSamplingSP isampling  = pcccm->getSampling();

				AdaptiveSampling* adSampling = dynamic_cast
					<AdaptiveSampling*>(isampling.get());

				GaussSampling* gaussSampling = dynamic_cast
					<GaussSampling*>(isampling.get());

				if (adSampling) {

					// generate a uniform sampling
					// run a small monte carlo
					/// process the results and config the sampling in the pathconfig

					IMFSamplingSP initSampling = adSampling->getInitialSampling();
					pcccm->setSampling(initSampling );

					Results samplingResults;

					model->getMCModel()->Price(gcdo, &localControl, &samplingResults);

					IObjectConstSP tupleList = samplingResults.retrieveRequestResult(
						OutputRequest::CONDITIONAL_LOSS_SAMPLES);

					const GeneralisedCDO::ConditionalLossMCProduct::ResultType* initResults =
						dynamic_cast<const
							GeneralisedCDO::ConditionalLossMCProduct::ResultType*>
								(tupleList.get());

					// isampling is the SP for adSampling
					IntArrayConstSP nonTrivial = countNonTrivial(initResults, initSampling);

					DoubleArrayConstSP initialPoints = initSampling->getPoints();

					adSampling->config(
						model->getMCModel()->getNbIters(true),
						initialPoints,
						nonTrivial);

					pcccm->setSampling(isampling);

				} else if (gaussSampling) {

					gaussSampling->config(model->getMCModel()->getNbIters(true));

				} else {

					throw ModelException("Not a known sampling scheme");
				}
			};

			// this will create the ConditionalLossMC product and run its
			// price method which stores the pairs of values
			model->getMCModel()->Price(gcdo, &localControl, &localResults);

			gcdo->clProduct = false;

			// now localResults has the lists of pairs (one list per maturity)

			IObjectConstSP tupleList = localResults.retrieveRequestResult(
				OutputRequest::CONDITIONAL_LOSS_SAMPLES);

			const GeneralisedCDO::ConditionalLossMCProduct::ResultType*
					lresultsTemp =
				dynamic_cast
					<const GeneralisedCDO::ConditionalLossMCProduct::ResultType*>
						(tupleList.get());

			lresults =
				const_cast<GeneralisedCDO::ConditionalLossMCProduct::ResultType*>
							(lresultsTemp);

			slices = processSamples(lresults, model->numberBuckets());

			gcdo->bucketsCL = FlatCDO2LossConfigSP(
				new FlatCDO2LossConfig(
						gcdo->today,
						flatPfolio,
						slices,
						model->getTimeLine()
					));

		} else {

			if (!gcdo->bucketsCL.get())
			{
				throw ModelException("Need to have run in pricing mode", method);
			};

		}; // end of if isPricing

		flatPfolio->getMarket(model->getImpliedModel(), gcdo->marketPointer.get());

		// Change the lossConfig in the local copy of the instrument

		ICreditLossConfigSP temp = gcdo->portfolio;

		gcdo->portfolio = gcdo->bucketsCL;

		// apply the alpha to the buckets

		CInstrumentArraySP insts = CInstrumentArraySP(new CInstrumentArray());

		insts->push_back(GeneralisedCDOSP::attachToRef(gcdo));

		ArrayInstrumentCollectionSP instruments =
			ArrayInstrumentCollectionSP(new ArrayInstrumentCollection(insts));

		CResultsArraySP resultsCollection = CResultsArraySP(new CResultsArray());

		resultsCollection->push_back(ResultsSP::attachToRef(results));

		// should this go here?
		gcdo->GetMarket(model->getImpliedModel(), gcdo->marketPointer);

		// price instrument using the effective curve, store price in results

		model->getImpliedModel()->PriceMulti(instruments,
				control,
				resultsCollection);

		gcdo->portfolio = temp;

		if (control->isPricing())
		{
			OutputRequest* request
					= control->requestsOutput(OutputRequest::CONDITIONAL_LOSS_MAP);

			if (request) {
				results->storeRequestResult(
					request,
					slices);
			}

			request = control->requestsOutput(OutputRequest::FLAT_CDO);

			if (request) {

				PortfolioName::NameViewArraySP names =
					PortfolioName::NameViewArraySP(new PortfolioName::NameViewArray());

				int k;
				for (k = 0; k < flatPfolio->numInnerLossConfigs(); k++)
				{
					const PortfolioName* name = dynamic_cast<const PortfolioName*>
						(flatPfolio->getInnerLossConfig(k).get());

					if (!name)
						throw ModelException("Not a flat CDO");

					names->push_back(name->getNameView());
				};

				results->storeRequestResult(
					request,
					names);
			}

			request = control->requestsOutput(OutputRequest::FLAT_CDO_TRANCHE_DECOMPOSITION);

			if (request) {
				results->storeRequestResult(
					request,
					gcdo->bucketsCL->getTrancheDecomposition());
			}

			request = control->requestsOutput(OutputRequest::CONDITIONAL_LOSS_SAMPLES);

			// localResults.retrieveRequestResult(
			//   	OutputRequest::CONDITIONAL_LOSS_SAMPLES);

			if (request) {
				results->storeRequestResult(
					request,
					GeneralisedCDO::ConditionalLossMCProduct::ResultTypeSP::attachToRef(lresults)
					//localResults.retrieveRequestResult(OutputRequest::CONDITIONAL_LOSS_SAMPLES)
					);
			}
		}
	} catch (exception& e) {

		throw ModelException(e, method);
	}
};

void GeneralisedCDO::ConditionalLossMCProduct::collectStateVars(
	IStateVariableCollectorSP svCollector) const
{
	static const string method = "ConditionalLossMCProduct::collectStateVars";
	try
	{
		instGen->collectStateVars(svCollector);
		regrGen->collectStateVars(svCollector);
    }
	catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

//------------------------------------------
//  IProtectionProvider method
//------------------------------------------
/* Checks if the input date is covered for protection */
bool GeneralisedCDO::isDateCoveredForProtection(const DateTime& date) const {
    return (!!cLeg && cLeg->isDateCoveredForProtection(date));
}


//--------------------------------------------
// GeneralisedCDO::MyConvolutionProduct::MyOutput methods
//--------------------------------------------
GeneralisedCDO::MyConvolutionProduct::MyOutput::MyOutput() :
    cLegPV(0.0),
    fLegPV(0.0),
    riskyDurationTotal(0.0),
    riskyNotionalsMean(0.0),
    risklessCFPV(0.0),
    cLegDebugUnitPrice(new DoubleArray()),
    cLegDebugUnitHistPrice(new DoubleArray()),
    fLegDebugUnitPrice(new DoubleArray()),
    fLegDebugUnitHistPrice(new DoubleArray())
{}

/** Returns the fair value (aka mark to market) given this
    Output object */
double GeneralisedCDO::MyConvolutionProduct::MyOutput::price() const {
    return (cLegPV - fLegPV);
}

/** Returns fee leg price for CCC view (positive cashflows only).
    Requires fLegDebugUnitPrice to be populated */
double GeneralisedCDO::MyConvolutionProduct::MyOutput::feeLegPriceForCCC(int longOrShort) {
    return priceCreditChargeView(*fLegDebugUnitPrice, -longOrShort);
}

/** Returns fee leg price for CCC view (positive cashflows only).
    Requires cLegDebugUnitPrice to be populated/ */
double GeneralisedCDO::MyConvolutionProduct::MyOutput::contingentLegPriceForCCC(int longOrShort) {
    return priceCreditChargeView(*cLegDebugUnitPrice, longOrShort);
}

/** Fudge for Kapital which wants fair value on [yield curve] spot
    date */
void GeneralisedCDO::MyConvolutionProduct::MyOutput::scaleByForwardFactor(
    double fwdFact)
{
    if (cLegDebugUnitPrice.get()){
        for (int i = 0; i < cLegDebugUnitPrice->size(); i++) {
            (*cLegDebugUnitPrice)[i]     *= fwdFact;
            (*cLegDebugUnitHistPrice)[i] *= fwdFact;
        }
    }
    if (fLegDebugUnitPrice.get()){
        for (int i = 0; i < fLegDebugUnitPrice->size(); i++) {
            (*fLegDebugUnitPrice)[i]     *= fwdFact;
            (*fLegDebugUnitHistPrice)[i] *= fwdFact;
        }
    }
    cLegPV *= fwdFact;
    fLegPV *= fwdFact;
    riskyDurationTotal *= fwdFact;
    risklessCFPV *= fwdFact;
}

/** Builds CDOPriceDetails object for fee leg */
CDOPriceDetails*
GeneralisedCDO::MyConvolutionProduct::MyOutput::makeFeeLegPriceDetails(
    int longOrShort) const
{
    return new CDOPriceDetails(
        fLegDebugUnitPrice,
        fLegDebugUnitHistPrice,
        priceCreditChargeView(*fLegDebugUnitPrice, -longOrShort));
}

/** Builds CDOPriceDetails object for contingent leg */
CDOPriceDetails*
GeneralisedCDO::MyConvolutionProduct::MyOutput::makeContingentLegPriceDetails(
    int longOrShort) const
{
    return new CDOPriceDetails(
        cLegDebugUnitPrice,
        cLegDebugUnitHistPrice,
        priceCreditChargeView(*cLegDebugUnitPrice, longOrShort));
}

//// just returns sum of unitPrice whose sign is the same as longOrShort
double GeneralisedCDO::MyConvolutionProduct::MyOutput::priceCreditChargeView(
    const DoubleArray& unitPrice,
    int                longOrShort)
{
    double priceCreditChargeView = 0.0;
    for (int i = 0; i < unitPrice.size(); i++) {
        if (Maths::sign(unitPrice[i]) == longOrShort){
            priceCreditChargeView += unitPrice[i];
        }
    }
    return priceCreditChargeView;
}


//---------------------------------
//GeneralisedCDO::MyConvolutionProduct methods
//---------------------------------
/** Construtor */
GeneralisedCDO::MyConvolutionProduct::MyConvolutionProduct(
    ConvolutionEngineConstSP convolutionEngine,
    const GeneralisedCDO*    cdo) :
        ConvolutionProduct(convolutionEngine, cdo->discount.getSP()),
        cdo(cdo),
        totalLongNotional(0.0),
        totalShortNotional(0.0),
        totalPastLongLoss(0.0),
        totalPastShortLoss(0.0),
        totalPastRecoveredNotional(0.0)
{
	static string method =
		"GeneralisedCDO::MyConvolutionProduct::MyConvolutionProduct";

    lossCfg = cdo->portfolio.get();
    CreditTrancheLossConfig* tranche =
        dynamic_cast<CreditTrancheLossConfig*>(lossCfg);

	FlatCDO2LossConfig* lossCfgTS =
			dynamic_cast<FlatCDO2LossConfig*>(lossCfg);

    if (tranche) {
        if (tranche->isCrossSub()) {
            throw ModelException(method, "Cannot price CrossSub tranches using "
                                 "this model");
        }

		portf.reset(dynamic_cast<const CDOPortfolio*>(tranche->getPortfolio().get()));
		if (!portf) {
			throw ModelException (method, "Convolution models can (currently) only "
                              "be applied to those tranches over CDOPortfolios, but"
                              "the tranche portfolio is not of type CDOPortfolio");
		}
		// Get the strikes from the tranche
		tranche->getTrancheStrikes(lowStrike, highStrike);
	}
    else if (lossCfgTS) {
		portf.reset(dynamic_cast<const CDOPortfolio*>
						(lossCfgTS->getPortfolio().get()));
		if (!portf) {
			throw ModelException (method, "Convolution models can (currently) only "
                              "be applied to those FlatCDO2s over CDOPortfolios, but"
                              "the FlatCDO2portfolio is not of type CDOPortfolio");
		}
		lowStrike = 0;
		highStrike = lossCfgTS->notional();
	}
    else {
        throw ModelException (method, "Convolution models can (currently) only "
                              "be applied to Tranches and FlatCDO2s ");
	}

    if (!portf)
        throw ModelException ("Internal error in GeneralisedCDO::"
                              "MyConvolutionProduct::MyConvolutionProduct 2");

    portf->getCompositionDetails(totalLongNotional,
                                 totalShortNotional,
                                 totalPastLongLoss,
                                 totalPastShortLoss,
                                 totalPastRecoveredNotional,
                                 totalPastLongRecoveredNotional,
                                 totalPastShortRecoveredNotional);

    /* Compute historic tranche notional reductions for the contingent
     * and fee legs */
    // JLH Maybe only need to do this if there is a contingent leg?

	CtgLegLossPerDefaultArraySP cLegPastTrancheLosses;

	if (tranche) {
        cLegPastTrancheLosses=
            lossCfg->historicContingentLegLosses(
                cdo->triggerDelay,
                cdo->defaultToCalculationDelay,
                cdo->lastTriggerDate,
                IBadDayAdjusterConstSP::attachToRef(cdo),
                cdo); // IProtectionProvider
    }
	else { // if (lossCfgTS)
		cLegPastTrancheLosses =
			lossCfgTS->historicContingentLegLosses(
				cdo->triggerDelay,
	            cdo->defaultToCalculationDelay,
		        cdo->lastTriggerDate,
			    IBadDayAdjusterConstSP::attachToRef(cdo),
				cdo); // IProtectionProvider
	}

    int numTrancheLosses = cLegPastTrancheLosses->size();
    cLegPastNotionalReductions.reset(new CashFlowArray(numTrancheLosses));
    cLegPayPastNotionalReductions.reset(new BoolArray(numTrancheLosses));
    double cumulativeAmount = 0.0;
    for (int i=0; i < numTrancheLosses; ++i) {
        // Here we need the cumulative losses on the contingent leg, so
        // add them up manually... while keeping track of which ones have
        // to be paid
        cumulativeAmount += (*cLegPastTrancheLosses)[i]->loss.amount;
        (*cLegPastNotionalReductions)[i] =
            CashFlow((*cLegPastTrancheLosses)[i]->loss.date, cumulativeAmount);
        (*cLegPayPastNotionalReductions)[i] =
            cdo->isDateCoveredForProtection(
				(*cLegPastTrancheLosses)[i]->defaultDate);
    }

    if (!(cdo->fLeg)) {
        // If there is no fee leg, do not bother computing the notional
        // reductions on the fee leg since a) they cannot be computed
        // properly (need the accrual periods) and b) they won't be used
        fLegPastReductions.reset(new FeeLegReductionPerDefaultArray(0));
    }
    else {
        // Obtain the fee leg accrual periods, since they are required to
        // determine when the losses/recovered amounts actually take effect
        // on the fee leg
        AccrualPeriodArrayConstSP accrualPeriods =
            cdo->fLeg->getAccrualPeriods();

        if (tranche) {
            fLegPastReductions = lossCfg->historicFeeLegReductions(
                cdo->triggerDelay,
                cdo->defaultToCalculationDelay,
                cdo->temporaryLossAmount->doubleValue(),
                cdo->lastTriggerDate,
                accrualPeriods,
			    IBadDayAdjusterConstSP::attachToRef(cdo),
                recoverNotional());
        }
		else if (lossCfgTS) {
            fLegPastReductions = lossCfgTS->historicFeeLegReductions(
				cdo->triggerDelay,
				cdo->defaultToCalculationDelay,
				cdo->temporaryLossAmount->doubleValue(),
				cdo->lastTriggerDate,
				accrualPeriods,
			    IBadDayAdjusterConstSP::attachToRef(cdo),
				recoverNotional());
		}
        else {
			throw ModelException("Not a tranche or a flatCDO2 loss config");
		}
    }
    // How to handle past decretion XXX
}


/** Returns the object that defines the losses for this product */
ICreditLossConfigConstSP
GeneralisedCDO::MyConvolutionProduct::getLossConfig() const {
    return cdo->getLossConfig();
}

/** Returns the last observation date for credit related data */
DateTime
GeneralisedCDO::MyConvolutionProduct::lastObservationDate() const {
    return cdo->lastObservationDate();
}

/** Returns the last 'scheduled' pay date (ie ignoring credit
    events).  Introduced for the output request IND_CDS_PAR_SPREAD
    - not clear if this is right date to use */
DateTime GeneralisedCDO::MyConvolutionProduct::lastPayDate() const {
    return cdo->maturityDate();
}

/** Returns the max of maxDate and the last date from when a yield
    curve is used (eg for discounting or rate estimation) by
    payoff method.  This method does not have to worry about what
    dates are used when the clean spread curves are built. Used to
    control when to stop tweaking */
DateTime GeneralisedCDO::MyConvolutionProduct::lastYCSensDate(
    const DateTime& maxDate) const
{
    return cdo->lastYCSensDate(maxDate);
}

/** Returns number of names contained within instrument */
int GeneralisedCDO::MyConvolutionProduct::numNames() const
{
    return lossCfg->numInnerLossConfigs();
}

/** Returns the representation of the underlying name corresponding to the
    specified index which must lie in range [0, numNames()-1] */
SingleCreditAssetConstSP
GeneralisedCDO::MyConvolutionProduct::nameAsset(int index) const
{
    return portf->nameAsset(index);
}

/** Returns the max of valueDate and protection start date
    for specified name (if any). The specified index must lie in
    range [0, numNames()-1]*/
const DateTime& GeneralisedCDO::MyConvolutionProduct::nameProtectionStartDate(
    int index) const
{
    return portf->getName(index)->getProtectionStartDate();
}

/** Returns the min of last date parameter and protection end date
    for specified name (if any). The specified index must lie in
    range [0, numNames()-1]. This if for name 'cutoff' */
const DateTime& GeneralisedCDO::MyConvolutionProduct::nameProtectionEndDate(
    int index,
    const DateTime& lastDate) const
{
    return portf->getName(index)->getProtectionEndDate(lastDate);
}

/** Returns the recovery for the name corresponding to the
    specified index which must lie in range [0, numNames()-1]. Note that
    the recovery can  be overridden at the trade level hence this method,
    which reflects any trade level overrides, should be used rather than
    the recovery off the par spread curve */
double GeneralisedCDO::MyConvolutionProduct::nameRecovery(int index) const {
    return portf->getName(index)->getNameRecovery();
}

/** Returns true if the name has defaulted where the name corresponds to the
    specified index which must lie in range [0, numNames()-1] */
bool GeneralisedCDO::MyConvolutionProduct::nameDefaulted(int index) const {
    return portf->getName(index)->hasDefaulted();
}

/** Returns the notional for the name corresponding to the
    specified index which must lie in range [0, numNames()-1] */
double GeneralisedCDO::MyConvolutionProduct::nameNotional(int index) const {
    return cdo->portfolio->getInnerLossConfig(index)->notional();
}

/** Returns the loss given default for the name corresponding to the
    specified index which must lie in range [0, numNames()-1] */
double GeneralisedCDO::MyConvolutionProduct::nameLossGivenDefault(
    int index) const
{
    return cdo->portfolio->getInnerLossConfig(index)->maxPossibleLoss();
}

/** Returns the 'beta' for the name corresponding to the specified
    index which must lie in range [0, numNames()-1] */
double GeneralisedCDO::MyConvolutionProduct::nameBeta(int index) const {
    return portf->getName(index)->getBeta();
}


/** Returns the 'loss decretion beta' */
double GeneralisedCDO::MyConvolutionProduct::lossDecretionBeta() const {
    return portf->lossDecretionBeta;
}


/** Returns the ranges for Loss Given Default for the name
    corresponding to the specified index which must lie in range
    [0, numNames()-1]. In the particular, the LGD of a name is given by
    nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
void GeneralisedCDO::MyConvolutionProduct::nameLGDRanges(
    int     index,       /* (I) */
    double& lgdNotional, /* (O) */
    double& lgdFloor,    /* (O) */
    double& lgdCap)      /* (O) */ const
{
    portf->getName(index)->
        nameLGDRanges(lgdNotional, lgdFloor, lgdCap);
}

/** Returns a CounterPartyCredit object containing data required for
    pricing CCC. Note may return null if no counter party information */
CounterPartyCreditConstSP
GeneralisedCDO::MyConvolutionProduct::getCounterParty() const
{
    return cdo->cptyInfo;
}

/** Returns the total notional of the portfolio */
double GeneralisedCDO::MyConvolutionProduct::portfolioNotional() const {
    return (totalLongNotional - totalShortNotional);
}

/** Returns the total long notional of the portfolio */
double GeneralisedCDO::MyConvolutionProduct::portfolioLongNotional() const {
    return totalLongNotional;
}

/** Returns the total short notional of the portfolio */
double GeneralisedCDO::MyConvolutionProduct::portfolioShortNotional() const {
    return totalShortNotional;
}

/** Returns the sum of the past losses of the portfolio which eat
    into the total notional returned by portfolioNotional() */
double GeneralisedCDO::MyConvolutionProduct::portfolioLoss() const {
    return (totalPastLongLoss - totalPastShortLoss);
}

/** Returns the sum of the long past losses of the portfolio which eat
    into the total notional returned by portfolioNotional() */
double GeneralisedCDO::MyConvolutionProduct::portfolioLongLoss() const {
    return totalPastLongLoss;
}

/** Returns the sum of the short past losses of the portfolio which eat
    into the total notional returned by portfolioNotional() */
double GeneralisedCDO::MyConvolutionProduct::portfolioShortLoss() const {
    return totalPastShortLoss;
}

/** Returns whether the instument requires notional to be
    recovered from the top of the portfolio NB may be later
    overridden by the model if it is deemed to be unnecessary */
bool GeneralisedCDO::MyConvolutionProduct::recoverNotional() const {
    return cdo->getRecoverNotional();
}

/** Returns the sum of the past recovered notional of the
    * portfolio which eat into the total notional returned by
    * portfolioNotional() */
double
GeneralisedCDO::MyConvolutionProduct::portfolioRecoveredNotional() const {
    return totalPastRecoveredNotional;
}

/** Returns the sum of the past long recovered notional of the
    * portfolio which eat into the total notional returned by
    * portfolioNotional() */
double
GeneralisedCDO::MyConvolutionProduct::portfolioLongRecoveredNotional() const
{
    return totalPastLongRecoveredNotional;
}

/** Returns the sum of the past short recovered notional of the
    * portfolio which eat into the total notional returned by
    * portfolioNotional() */
double
GeneralisedCDO::MyConvolutionProduct::portfolioShortRecoveredNotional() const
{
    return totalPastShortRecoveredNotional;
}

/** Returns the strikes of the tranche. These are relative to the values
    returned by portfolioRanges() */
void GeneralisedCDO::MyConvolutionProduct::trancheStrikes(
    double& lowStrk,
    double& highStrk) const
{
    lowStrk  = lowStrike;
    highStrk = highStrike;
}

/** map loss interpolation string to CCMPriceUtil enum */
CCMPriceUtil::ExpLossType
GeneralisedCDO::MyConvolutionProduct::mapLossInterpolationType(
    const string& lossInterpolation)
{
    if (lossInterpolation == EffectiveCurve::FLAT_FORWARD){
        return CCMPriceUtil::FLATFWD;
    }
    if (lossInterpolation == EffectiveCurve::LINEAR){
        return CCMPriceUtil::LINEAR;
    }
    throw ModelException("GeneralisedCDO::MyConvolutionProduct",
                         "Unknown interpolation type "+lossInterpolation);
}


/** Computes payoff for instrument given effectiveCurve */
ConvolutionProduct::OutputSP GeneralisedCDO::MyConvolutionProduct::payoff(
    EffectiveCurveSP     cLegEffectiveCurve,
    EffectiveCurveSP     fLegEffectiveCurve,
    double               forwardFactor,
    IForwardRatePricerSP model) const
{
    static const string routine("GeneralisedCDO::MyConvolutionProduct::payoff");
    MyOutputSP myOutput(new MyOutput());

    try {
        /* valDateCF: price CF after this date */
        const DateTime& valDateCF = cfCutOffDate();
        if (cdo->cLeg.get()) {
            /* calculate contingent leg price */
            double initialNotional = highStrike - lowStrike;
            double cLegOutstandingNotional = trancheOutstandingNotional(false);

            myOutput->cLegPV =
                cdo->cLeg->price(initialNotional,
                                 cLegOutstandingNotional,
                                 getToday(),
                                 valDateCF,
                                 cLegEffectiveCurve,
                                 *cLegPastNotionalReductions,
                                 *cLegPayPastNotionalReductions,
                                 true, // to do: only if needed
                                 *myOutput->cLegDebugUnitPrice,
                                 *myOutput->cLegDebugUnitHistPrice,
                                 IBadDayAdjusterConstSP::attachToRef(cdo.get()));
        }

        if (cdo->fLeg.get()) {
            /* To price the rebates the fee leg needs to be awarawaree of the
             * payment details (pay as you go, days delay / pay dates, etc).
             * Since they are implemented in the contingent leg AND it is too
             * late to move them into the instrument AND it has been estimated
             * that duplicating them in the fee leg would cause confussion
             * among the library users (specially the observation arrays), the
             * preferred approach is to get them from the contingent leg and
             * pass them to the fee leg. So do just that... */
            BoolArraySP     payAsYouGoArray;
            IntArraySP      numDelayDaysArray;
            DateTimeArraySP startDates;
            DateTimeArraySP endDates;
            DateTimeArraySP paymentDates;

            if (cdo->cLeg.get()) {
                cdo->cLeg->getPaymentInformation(payAsYouGoArray,
                                                 numDelayDaysArray,
                                                 startDates,
                                                 endDates,
                                                 paymentDates);
            } // else pass in null parameters

            const bool recNotional = recoverNotional();

            // calculate fee leg price
            double fLegOutstandingNotional =
                trancheOutstandingNotional(recNotional);

            // Compute the rebate payments
			CashFlowArraySP fLegRebates =
                lossCfg->historicRebatePayments(cdo.get(), // IRebateCalculator
                                                fLegPastReductions,
                                                model,
                                                recNotional);

            CashFlowArraySP fLegPastNotionalReductions =
                FeeLegReductionPerDefault::getReductions(fLegPastReductions,
                                                         true); // losses

            if (recNotional) {
                // get recovered notional cashflows
                CashFlowArraySP fLegPastRecoveredNotional =
                    FeeLegReductionPerDefault::getReductions(fLegPastReductions,
                                                             false); // rec notional

                // and merge them
                fLegPastNotionalReductions =
                    CashFlow::merge(fLegPastNotionalReductions,
                                    fLegPastRecoveredNotional);
                CashFlow::aggregate(*fLegPastNotionalReductions);
            }

            // Obtain the cumulative cashflows:
            CashFlowArraySP cumulativeFeeLegReductions =
                CashFlow::accumulateCashFlows(fLegPastNotionalReductions);

            myOutput->fLegPV =
                cdo->fLeg->price(getToday(),
                                 valDateCF,
                                 fLegEffectiveCurve,
                                 cdo->discount,
                                 lowStrike,
                                 highStrike,
                                 fLegOutstandingNotional,
                                 *cumulativeFeeLegReductions,
                                 myOutput->riskyDurationTotal,
                                 myOutput->riskyNotionalsMean,
                                 myOutput->risklessCFPV,
                                 true, // to do: only if needed
                                 *myOutput->fLegDebugUnitPrice,
                                 *myOutput->fLegDebugUnitHistPrice,
                                 fLegRebates,
                                 payAsYouGoArray,
                                 numDelayDaysArray,
                                 startDates,
                                 endDates,
                                 paymentDates,
                                 IBadDayAdjusterConstSP::attachToRef(cdo.get()),
                                 model);
        }
        // then scale appropriate results by forwardFactor
        myOutput->scaleByForwardFactor(forwardFactor);
        return myOutput;
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Compute the credit charge using supplied results and credit charge
    view type. To do: sort out rep of creditChargeViewType */
double GeneralisedCDO::MyConvolutionProduct::computeCreditCharge(
    const EffectiveCurveSP cptyCurve,
    const string&          creditChargeViewType,
    OutputSP               output,
    OutputSP               cccOutput) const
{
    MyOutput& myOutput = static_cast<MyOutput&>(*output);
    MyOutput& myCCCOutput = static_cast<MyOutput&>(*cccOutput);
    int longOrShort = cdo->isLong? 1 : -1;
    double priceUncond, priceCond;
    if (creditChargeViewType == ConvolutionProduct::CC_AGGRESSIVE) {
        priceUncond = myOutput.cLegPV - myOutput.fLegPV;
        priceCond = myCCCOutput.cLegPV - myCCCOutput.fLegPV;
    } else {
        priceUncond = myOutput.contingentLegPriceForCCC(longOrShort) -
            myOutput.feeLegPriceForCCC(longOrShort);
        priceCond =
            myCCCOutput.contingentLegPriceForCCC(longOrShort) -
            myCCCOutput.feeLegPriceForCCC(longOrShort);
    }
    /* determine the end of exposure to cpty */
    DateTime lastDay(cfCutOffDate());
    if (cdo->cLeg.get()) {
        lastDay = lastDay.max(cdo->cLeg->lastPayDate(
                      IBadDayAdjusterConstSP::attachToRef(cdo.get())));
    }

    if (cdo->fLeg.get()) {
        lastDay = lastDay.max(cdo->fLeg->getLastPayDate());
    }

    return cdo->cptyInfo->calculateCSACapPrice(cdo->isLong,
                                                myOutput.price(),
                                                priceCond,
                                                priceUncond,
                                                cptyCurve,
                                                cdo->today,
                                                lastDay);
}

/** Populate CResults object with price and any output requests etc. Note
    cccOutput may be null if CCC not being computed */
void GeneralisedCDO::MyConvolutionProduct::storeExtraResults(
    CResults*  results,
    CControl*  control,
    OutputSP   output,
    OutputSP   cccOutput,
    IForwardRatePricerSP model,   // may be 0
    const IConvolutionModel* convolutionModel) const
{
    MyOutput& myOutput = static_cast<MyOutput&>(*output);
    MyOutput* myCCCOutput = static_cast<MyOutput*>(cccOutput.get());
    int longOrShort = cdo->isLong? 1 : -1;
    OutputRequest* request; // reused across different OutputRequests

    // TRANCHE_IMPLIED_SPREAD
    request=control->requestsOutput(OutputRequest::TRANCHE_IMPLIED_SPREAD);
    if (request) {
        if (myOutput.riskyDurationTotal != 0.0){
            double modifiedContingentLegPrice =
                myOutput.cLegPV - myOutput.risklessCFPV;
            results->storeRequestResult(request,
                                        modifiedContingentLegPrice/
                                        myOutput.riskyDurationTotal);
        } else {
            results->storeNotApplicable(request);
        }
    }
    // TRANCHE_RISKY_DURATION
    request=control->requestsOutput(OutputRequest::TRANCHE_RISKY_DURATION);
    if (request) {
        double outstdNtl = trancheOutstandingNotional(recoverNotional());
        double initialNtl = trancheNotional();
        if (myOutput.riskyNotionalsMean != 0.0 && outstdNtl != 0.0 && initialNtl != 0.0){
            double outstdNtlRatio = outstdNtl / initialNtl;

            results->storeRequestResult(request,
                                        myOutput.riskyDurationTotal/
                                        (outstdNtlRatio*myOutput.riskyNotionalsMean));
        } else {
            results->storeNotApplicable(request);
        }
    }
    // TRANCHE_OUTSTANDING_NOTIONAL
    request=control->requestsOutput(OutputRequest::TRANCHE_OUTSTANDING_NOTIONAL);
    if (request) {
        results->storeRequestResult(request, trancheOutstandingNotional(recoverNotional()));
    }
    // TRANCHE_CONTINGENT_LEG_PRICE
    request = control->requestsOutput(OutputRequest::
                                        TRANCHE_CONTINGENT_LEG_PRICE);
    if (request) {
        CDOPriceDetailsSP price(
            myOutput.makeContingentLegPriceDetails(longOrShort));
        CDOPriceDetailsSP priceCond(
            myCCCOutput?
            myCCCOutput->makeContingentLegPriceDetails(longOrShort):
            new CDOPriceDetails());
        CDOLegOutputSP cLegOutput(new CDOLegOutput(price, priceCond));
        results->storeRequestResult(request, cLegOutput);
    }
    // TRANCHE_FEE_LEG_PRICE
    request = control->requestsOutput(OutputRequest::TRANCHE_FEE_LEG_PRICE);
    if (request) {
        CDOPriceDetailsSP price(
            myOutput.makeFeeLegPriceDetails(longOrShort));
        CDOPriceDetailsSP priceCond(
            myCCCOutput?
            myCCCOutput->makeFeeLegPriceDetails(longOrShort):
            new CDOPriceDetails());
        CDOLegOutputSP fLegOutput(new CDOLegOutput(price, priceCond));
        results->storeRequestResult(request, fLegOutput);
    }
    // TRANCHE_CONTINGENT_LEG_FV
    request =
        control->requestsOutput(OutputRequest::TRANCHE_CONTINGENT_LEG_FV);
    if (request) {
        results->storeRequestResult(request, myOutput.cLegPV);
    }
    // TRANCHE_FEE_LEG_FV
    request = control->requestsOutput(OutputRequest::TRANCHE_FEE_LEG_FV);
    if (request) {
        results->storeRequestResult(request, myOutput.fLegPV);
    }

    // KNOWN_CASHFLOWS
    CashFlowArraySP cfls;
    request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
    OutputRequest* request2 =
        control->requestsOutput(OutputRequest::PAYMENT_DATES);
    if (request || request2) {
        double fLegOutstandingNotional = trancheOutstandingNotional(true);
        CashFlowArraySP fLegPastNotionalReductions =
            FeeLegReductionPerDefault::getReductions(fLegPastReductions,
                                                         true); // losses
        if (recoverNotional()) {
            // get recovered notional cashflows
            CashFlowArraySP fLegPastRecoveredNotional =
                FeeLegReductionPerDefault::getReductions(fLegPastReductions,
                                                         false); // rec notional

            // and merge them
            fLegPastNotionalReductions =
                CashFlow::merge(fLegPastNotionalReductions,
                                fLegPastRecoveredNotional);
            CashFlow::aggregate(*fLegPastNotionalReductions);
        }
        // Obtain the cumulative cashflows:
        CashFlowArraySP cumulativeFeeLegReductions =
            CashFlow::accumulateCashFlows(fLegPastNotionalReductions);

        CashFlowArraySP cfls(cdo->generateKnownCashflows(
            fLegOutstandingNotional,
            *cumulativeFeeLegReductions,
            *cLegPastNotionalReductions,
            *cLegPayPastNotionalReductions,
            model));

        if (request){
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    cdo->discount->getCcy(),
                                                    cfls.get());
        }
        // PAYMENT_DATES
        if (request2) {
            DateTimeArray paydates = CashFlow::dates(*cfls);
            OutputRequestUtil::recordPaymentDates(control,
                                                  results,
                                                  &paydates);
        }
    }

    // DURATION_WEIGHTED_AVG
    request = control->requestsOutput(OutputRequest::DURATION_WEIGHTED_AVG);
    if (request) {
        // This output request only makes sense if model = base correlation and
        // the loss config is a tranche.
        // This is REALLY hacky!
        const CreditMetricsBaseCorrelation* bc =
            dynamic_cast<const CreditMetricsBaseCorrelation*>(convolutionModel);

        CreditTrancheLossConfig* tranche =
            dynamic_cast<CreditTrancheLossConfig*>(lossCfg);
        if (!bc || !tranche) {
            results->storeNotApplicable(request);
        }
        else {
            try {
                computeAndStoreDurWeightAvg(bc, tranche, results, request);
            }
            catch (exception& ) {
                results->storeNotApplicable(request);
            }
        }
    }
}


void GeneralisedCDO::MyConvolutionProduct::computeAndStoreDurWeightAvg(
    const CreditMetricsBaseCorrelation* bc,
    const CreditTrancheLossConfig* tranche,
    CResults* results,
    OutputRequest* request) const
{
    const map<string, double>& dwas = bc->getDurationWeightedAverage(
        CreditTrancheLossConfigConstSP::attachToRef(tranche),
        lastObservationDate(),
        cdo->discount.getSP());

    map<string, double>::const_iterator dwasIter = dwas.begin();
    for (; dwasIter != dwas.end(); dwasIter++) {
        string name = (*dwasIter).first;
        OutputNameConstSP myOutputName(new OutputName(name));
        IObjectSP result(CDouble::create((*dwasIter).second));

        results->storeRequestResult(request, result, myOutputName);
    }

    return;
}


/** Returns the object that defines the losses for this product */
ICreditLossConfigConstSP GeneralisedCDO::getLossConfig() const {
    return portfolio;
}


/** Creates an instance of an ConvolutionProduct */
IGeneralisedConvolutionProduct* GeneralisedCDO::createProduct(
    ConvolutionEngineConstSP model) const
{
    // Convolution Engines can currently be used to price tranches and NtDs.
    // Here need to differentiate both cases and return the right type of
    // product
    if (NToDefaultLossConfig::TYPE->isInstance(portfolio)) {
        return new FtDProduct(model, this);
    }
    // else assume tranche and return a ConvolutionProduct instance
    return new MyConvolutionProduct(model, this);
}


//------------------------------------------------------------------------------
//                             FDModel::IIntoProduct Methods
//------------------------------------------------------------------------------
FDProductSP GeneralisedCDO::createProduct(FDModel * model) const
{
	return FDProductSP(new CDOTree(GeneralisedCDOConstSP::attachToRef(this),
                                   model));
}


//------------------------------------------------------------------------------
//                             CDOLegOutput
//------------------------------------------------------------------------------

//// Simple constructor - just takes references to CDOPriceDetailsSP
CDOLegOutput::CDOLegOutput(CDOPriceDetailsSP price,
                           CDOPriceDetailsSP priceCond)
: CObject(TYPE),
  price(price),
  priceCond(priceCond)
{}

CDOLegOutput::CDOLegOutput()
: CObject(TYPE)
{}

IObject* CDOLegOutput::defaultConstructor() {
    return new CDOLegOutput();
}

/** Invoked once at start up when this class is 'loaded' */
void CDOLegOutput::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CDOLegOutput, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(price,     "price");
    FIELD(priceCond, "price cond ND");
}

CClassConstSP const CDOLegOutput::TYPE = CClass::registerClassLoadMethod(
    "CDOLegOutput", typeid(CDOLegOutput), load);


//------------------------------------------------------------------------------
//                             CDOPriceDetailsSP
//------------------------------------------------------------------------------

CDOPriceDetails::CDOPriceDetails() :
    CObject(TYPE),
    price(0.0),
    histPrice(0.0),
    priceCreditChargeView(0.0)
{}

    /** Creates CDOPriceDetails from supplied inputs. Note that the DoubleArrays
        must be of the same length. */
CDOPriceDetails::CDOPriceDetails(
    DoubleArraySP unitPrice,  /* price for each leg unit */
    DoubleArraySP unitHistPrice, /*price for each leg unit due to historical default*/
    double        priceCreditChargeView) :
    CObject(TYPE),
    price(0.0),
    histPrice(0.0),
    priceCreditChargeView(priceCreditChargeView),
    unitPrice(unitPrice),
    unitHistPrice(unitHistPrice)
{
    int nbPrices = unitPrice->size();
    for (int i = 0; i < nbPrices; i++) {
        price     += (*unitPrice)[i];
        histPrice += (*unitHistPrice)[i];
    }
}

IObject* CDOPriceDetails::defaultConstructor() {
    return new CDOPriceDetails();
}

/** Invoked once at start up when this class is 'loaded' */
void CDOPriceDetails::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CDOPriceDetails, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(price,                 "total leg price");
    FIELD(histPrice,             "total leg price due to historical default");
    FIELD(priceCreditChargeView, "total price of CCC view (positive cf only)");
    FIELD       (unitPrice,             "price for each leg unit");
    FIELD       (unitHistPrice,         "price for each leg unit due to historical default");

    FIELD_MAKE_OPTIONAL(unitPrice);
    FIELD_MAKE_OPTIONAL(unitHistPrice);
}

CClassConstSP const CDOPriceDetails::TYPE = CClass::registerClassLoadMethod(
    "CDOPriceDetails", typeid(CDOPriceDetails), load);



/* For class loading (avoid having header file) */
bool GeneralisedCDOLoad() {
    return (GeneralisedCDO::TYPE != 0);
}

void
GeneralisedCDO::generateMCTimelines(
	DateTimeArraySP& modifiedFeeLegObservationDates,
	DateTimeArraySP& productTimeline,
	DateTimeArraySP& discFactorTimeline,
	IntArraySP& dateToDiscFactorIndex
	) const
{
	static const string routine = "GeneralisedCDO::generateMCTimeline";
	try
	{
//check for whether either of FeeLeg or ContingentLeg are present
		if  ( (this->fLeg.get() == 0) && (this->cLeg.get() == 0) )
			throw ModelException(
				routine,
				"Cannot generate MC timeline as both fee and contingent leg are missing");

		list<DateTime> discFactorTimelineList;
		list<DateTime> timelineList;

//add valueDate to the lists
		DateTime valueDate =this->getValueDate();
		timelineList.push_back(valueDate);
		discFactorTimelineList.push_back(valueDate);

//FIRST add all fee leg timeline dates that occur after today
		if (this->fLeg.get()) //fee leg exists
		{
//get the cashflow dates - to update the discount factor timeline
			IForwardRatePricerSP model; //later this would be the interest rate model
			AbstractCashFlowArrayConstSP abstractCashFlowArray =
				this->fLeg->getCashFlows(model);

			if (abstractCashFlowArray.get())
			{
				for (int i=0; i< abstractCashFlowArray->size(); ++i)
					discFactorTimelineList.push_back((*abstractCashFlowArray)[i]->getPayDate());
			}

//get only the risky payment dates and accrual periods
            DateTimeArraySP feePayDates = fLeg->getRiskyCashFlowDates();
            AccrualPeriodArrayConstSP accruals = fLeg->getRiskyAccrualPeriods();
			modifiedFeeLegObservationDates = fLeg->getRiskyObservationDates();
			CouponNotionalTypesArraySP couponNotionalTypes = fLeg->getRiskyCouponNotionalTypes();

			int numFees = feePayDates->size();
            //check accruals are the same size
            if (	(accruals->size() != numFees) ||
					(modifiedFeeLegObservationDates->size() != numFees )  ||
					(couponNotionalTypes->size() != numFees ) )
            {
                throw ModelException(
                    routine,
                    "Internal error : mismatch in fee dates");
            }

            for (int i=0; i< numFees; ++i)
			{
				if ( (*feePayDates)[i] >= valueDate)
				{
					if ( (*couponNotionalTypes)[i] == AVERAGE )
					{
						//when both start and end date are both in the future
						if ( ((*accruals)[i])->startDate() >= valueDate)
						{
							DateTime midDate = 	((*accruals)[i])->startDate().midPoint(((*accruals)[i])->endDate());
							timelineList.push_back(midDate);
							(*modifiedFeeLegObservationDates)[i] = midDate;
						}
						else
						{
						//when start date in the past and end date in the future
							timelineList.push_back((*accruals)[i]->startDate());
							if ( (*accruals)[i]->endDate() > valueDate)
							{
								DateTime midDate = 	valueDate.midPoint((*accruals)[i]->endDate());
								timelineList.push_back(midDate);
								(*modifiedFeeLegObservationDates)[i] = midDate;
							}
							else //when both start and end date are in the past
							{
								(*modifiedFeeLegObservationDates)[i] = DateTime();//empty date
							}

						}
					}
					else
						timelineList.push_back((*modifiedFeeLegObservationDates)[i]);
				}
			}
		}

//SECOND, add the past default dates
		CtgLegLossPerDefaultArraySP cLegPastTrancheLosses =
			this->portfolio->historicContingentLegLosses(
				this->triggerDelay,
				this->defaultToCalculationDelay,
				this->lastTriggerDate,
				IBadDayAdjusterConstSP::attachToRef(this),
				this); // IProtectionProvider

		if (cLegPastTrancheLosses.get())
		{
			for (unsigned int i=0; i< cLegPastTrancheLosses->size() ; ++i)
			{
				DateTime date = (*cLegPastTrancheLosses)[i]->defaultDate;
				timelineList.push_back(date);
			}
		}

//THIRD, then add all contingent leg timeline dates that occur after today
		BoolArraySP     lpayAsYouGoArray;
		IntArraySP      lnumDelayDaysArray;
		DateTimeArraySP lstartDateArray;
		DateTimeArraySP lendDateArray;
		DateTimeArraySP lpaymentDateArray;

		if (this->cLeg.get()) //contingent leg exists
		{
			this->cLeg->getPaymentInformation(
				lpayAsYouGoArray,
				lnumDelayDaysArray,
				lstartDateArray,
				lendDateArray,
				lpaymentDateArray);

			//do a safety check
			if ( (lstartDateArray->size() != lendDateArray->size()) ||
				(lendDateArray->size() != lpaymentDateArray->size()) )
				throw ModelException(routine, "Contingent leg: StartDates, "
                                     "EndDates, PayDates are of different sizes");

			for (int i=0; i< lendDateArray->size(); ++i)
			{
				DateTime& date = (*lendDateArray)[i];

				if (date > valueDate)
					timelineList.push_back(date);

				date = (*lstartDateArray)[i];

				if (date > valueDate)
					timelineList.push_back(date);

				date = (*lpaymentDateArray)[i];

				if (date > valueDate)
					discFactorTimelineList.push_back(date);
			}
		}

		addCutOffs(timelineList, portfolio.get());

		sortAndRemoveDuplicates(
			*productTimeline,
			timelineList);

		if (this->cLeg.get()) //contingent leg exists
		{
//now if payAsYouGo, we need to have a payDate (in the discount factor timeline) for each productTimeline period
			int lastIndex = 0;
			for (int i=0; i< productTimeline->size(); ++i)
			{
//search whether that timeline date belongs to a contingent leg period
				for (int k = lastIndex; k < lstartDateArray->size(); ++k)
				{
					const DateTime& date = (*productTimeline)[i]; //stochastic defaults occur only at the product timeline
					if (  ( (i==0)   &&  ((*productTimeline)[0] == (*lstartDateArray)[k]))  ||
						(( date >  (*lstartDateArray)[k] ) && ( date <= (*lendDateArray)[k])) )
					{
						lastIndex = k;//belongs to the kth contingentleg index

						DateTime payDate;
						if ((*lpayAsYouGoArray)[k])
							payDate = date.rollDate((*lnumDelayDaysArray)[k]);//we apply the delay
						else //it gets paid at the next payment period
							payDate = (*lpaymentDateArray)[k];

						discFactorTimelineList.push_back(payDate);
						k = lstartDateArray->size();//to exit the loop
					}
				}
			}
		}

		sortAndRemoveDuplicates(
			*discFactorTimeline,
			discFactorTimelineList);

//the following is done as if the trade ends today the the MC framework would conclude that the trade has no future and
//and price object will not get updated.
		if (productTimeline->back() == valueDate)
			productTimeline->push_back(valueDate.rollDate(1));

//builds a map between a date and the Index on the timeline
		int maxDateValue = discFactorTimeline->back().getDate();
		dateToDiscFactorIndex->resize(maxDateValue+1, -1);

		for (int j=0; j <discFactorTimeline->size(); ++j)
		{
			int dateValue = (*discFactorTimeline)[j].getDate();
			(*dateToDiscFactorIndex)[dateValue] = j;
		}
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}
}

/** utility method that sorts a list, removes duplicates
	and converts it into an arraySP
*/
IMCProduct*
GeneralisedCDO::createProduct(const MonteCarlo* model) const
{
	static const string method = "GeneralisedCDO::createProduct";

	if (clProduct)
		return new ConditionalLossMCProduct(this, model);
	try
	{
		DateTimeArraySP discFactorTimeline = DateTimeArraySP(new DateTimeArray());
		IntArraySP dateToDiscFactorIndex = IntArraySP(new IntArray());
		DateTimeArraySP productTimeline = DateTimeArraySP(new DateTimeArray());
		DateTimeArraySP modifiedFeeLegObservationDates = DateTimeArraySP(new DateTimeArray());

		this->generateMCTimelines(
			modifiedFeeLegObservationDates,
			productTimeline,
			discFactorTimeline,
			dateToDiscFactorIndex);

		// this simSeries will be used to initialise the MCProductClient object
		// JCP do we need number of Inners? 1 is just as good
		// as this really does nothing

		int numberOfInners = this->portfolio->numInnerLossConfigs();
		SimSeriesSP sims(new SimSeries(1));

		sims->addDates(clProduct?*auxTimeLine:*productTimeline);

		return new GeneralisedCDO::GeneralisedCDOMC(
			GeneralisedCDOConstSP(this),
			sims,
			model->getNbIters(true),
			model->getNbSubSamples(),
			clProduct?auxTimeLine:productTimeline,
			discFactorTimeline,
			dateToDiscFactorIndex,
			modifiedFeeLegObservationDates);
    }
	catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

// ########################################################################

// ##  Methods of GeneralisedCDOMC::GeneralisedCDOMC follow
GeneralisedCDO::GeneralisedCDOMC::GeneralisedCDOMC(
	const GeneralisedCDOConstSP& cdo,
	const SimSeriesConstSP&   simSeries,
	int numberOfIterations,
	int numberOfSubSamples,
	const DateTimeArrayConstSP& productTimeline,
	const DateTimeArrayConstSP& discFactorTimeline,
	const IntArrayConstSP& dateToDiscFactorIndex,
	const DateTimeArrayConstSP& modifiedFeeLegObservationDates
	): //needed for the MCProductClient
	MCProductClient(
		cdo->getValueDate(),
		cdo->discount.get(),
		simSeries),
		cdo(cdo),
		productTimeline(productTimeline), //must pass another timeline for the discFactor
		modifiedFeeLegObservationDates(modifiedFeeLegObservationDates),
		trancheFeeLegRFPrice(numberOfIterations, numberOfSubSamples),
		trancheFeeLegPrice(numberOfIterations, numberOfSubSamples),
		unitCouponTrancheFeeLegPrice(numberOfIterations, numberOfSubSamples),
		trancheContingentLegPrice(numberOfIterations, numberOfSubSamples)
{
	static const string method = "GeneralisedCDO::GeneralisedCDOMC constructor";

	try
	{
		int timelineSize = productTimeline->size();

		ICreditLossConfigConstSP creditLossConfig = cdo->portfolio;
		AccrualPeriodArrayConstSP accrualPeriods;

		if (this->cdo->fLeg.get())
		{
			accrualPeriods = cdo->fLeg->getAccrualPeriods();

			creditFeeLegSVGen = this->cdo->fLeg->createSVGen(
															cdo->getValueDate(),
															*modifiedFeeLegObservationDates,
															*productTimeline,
															*dateToDiscFactorIndex,
															creditLossConfig->notional());
			creditFeeLegSV = creditFeeLegSVGen->createNewSV(creditFeeLegSV, 0);
		}

		// check whether it can support MC
		const ICreditLossConfigIndexedSVGenMC* creditLossConfigIndexedSVGenMC =
			dynamic_cast<const ICreditLossConfigIndexedSVGenMC*>
				(creditLossConfig.get());

		if (creditLossConfigIndexedSVGenMC == 0)
			throw ModelException(
						method,
						"Embedded loss config cannot support MC");

		// create the state variable gen and the state variable for
		// the CreditLossConfig
		indexedSVGen = creditLossConfigIndexedSVGenMC->createIndexedSVGen(
			productTimeline,
			cdo->triggerDelay,
			cdo->defaultToCalculationDelay,
			cdo->temporaryLossAmount->doubleValue(),
			cdo->lastTriggerDate,
			accrualPeriods,
			IBadDayAdjusterConstSP(cdo),
			IProtectionProviderConstSP(cdo),
			IRebateCalculatorConstSP(cdo),
			cdo->getRecoverNotional());
			//true); // JCP recover notional);

		// We create the SV here. This is because the product will have the
		// entire hierarchy of SVs. And the engine willonly update the lowest
		// level of this SV.
		indexedSV = indexedSVGen->createNewSV(NULL);
		double lossConfigNotional = indexedSV->getNotional();

//create the state variable gen for the DiscountFactor
		discFactorSVGen = SVGenDiscFactorConstSP(
							new SVGenDiscFactor(
								cdo->today,
								cdo->discount.getSP(),
								*discFactorTimeline));

//aggregate contingent leg information that is to used for pricing
		if (this->cdo->cLeg.get())
		{
			CreditContingentLegBase* contLeg = 0;
			try
			{
				contLeg =  dynamic_cast<CreditContingentLegBase*>(this->cdo->cLeg.get());
			}
			catch(exception& e1)
			{
				throw ModelException(
					e1,
					method,
					"Unable to downcast to cdoContingentLeg base");
			}
			this->cdo->cLeg->getPaymentInformation(
				this->lpayAsYouGoArray,
				this->lnumDelayDaysArray,
				this->lstartDateArray,
				this->lendDateArray,
				this->lpaymentDateArray);

//fetch the notionals
			contNotionals = contLeg->notionals();
			if (this->lpayAsYouGoArray->size() != this->contNotionals.size())
				throw ModelException(
					method,
					"getPaymentInformation(...) and notionals(..) returns arrays of different sizes");

//compute the discount factors
			cLegDiscountFactors.resize(timelineSize);
			for (int i=0; i< timelineSize; ++i)
				cLegDiscountFactors[i] = cdo->discount->pv((*productTimeline)[i]);

			contLegDiscFactorIndex.resize(timelineSize, -1);

//INITIALISING STEPS TO EASE CONTINGENT LEG CALCULATIONS
//we also assume that there is no overlap in contingent leg periods and they are in ascending order of start Dates

//calculate the payLoss values and contingent leg scaling factors
			contingentLegScalingFactors.resize(timelineSize, 0);
			payLoss.resize(timelineSize, false);
			int lastIndex = 0;

//Note that we have dates on the productTimeline that are before today
			for (int i=0; i< productTimeline->size(); ++i)
			{
//search whether that timeline date belongs to a contingent leg period
				for (int k = lastIndex; k < lstartDateArray->size(); ++k)
				{
					const DateTime& date = (*productTimeline)[i]; //stochastic defaults occur only at the product timeline
					if (  ( (i==0)   &&  ((*productTimeline)[0] == (*lstartDateArray)[k]))  ||
						(( date >  (*lstartDateArray)[k] ) && ( date <= (*lendDateArray)[k])) )
					{
						lastIndex = k;//belongs to the kth contingentleg index
						payLoss[i] = true;
						contingentLegScalingFactors[i] =  contNotionals[k]/lossConfigNotional;

						DateTime payDate;
						if ((*lpayAsYouGoArray)[k])
							payDate = date.rollDate((*lnumDelayDaysArray)[k]);//we apply the delay to the default date
						else //it gets paid at the next payment period
							payDate = (*lpaymentDateArray)[k];

//now search for the timeline date that is nearest to the payDate
						contLegDiscFactorIndex[i] = payDate.findNearest(*discFactorTimeline);
						k = lstartDateArray->size();//to exit the loop
					}
				}
			}
		}
	}
	catch (exception& e)
	{
        throw ModelException(e, method);
    }
}

void
GeneralisedCDO::GeneralisedCDOMC::collectStateVars(IStateVariableCollectorSP svCollector) const
{
	static const string method = "GeneralisedCDOMC::collectStateVars";
	try
	{
		this->indexedSVGen->collectStateVars(svCollector);
		svCollector->append(discFactorSVGen.get());
    }
	catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

void
GeneralisedCDO::GeneralisedCDOMC::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen)
{
	static const string method = "GeneralisedCDOMC::pathGenUpdated";
	try
	{
		if(newPathGen)
		{
//fetching the credit loss config SV
			IStateVariableSP stateVariable =
				indexedSVGen->create(
					indexedSV,
					newPathGen);//will return the object from the svDBase
			indexedSV = ICreditLossConfigIndexedSVSP(dynamic_cast<ICreditLossConfig::IIndexedSVGen::ISV*>(stateVariable.get()));

//fetching the discount factor SV
			stateVariable =
				discFactorSVGen->create(
					discFactorSV,
					newPathGen);//will return the object from the svDBase
			discFactorSV = SVDiscFactorSP(dynamic_cast<SVDiscFactor*>(stateVariable.get()));
		}
    }
	catch (exception& e)
	{
		throw ModelException(e, method);
	}
}

void
GeneralisedCDO::GeneralisedCDOMC::payoff(
	const IPathGenerator*  pathGen,
	IMCPrices& prices)
{
	static const string method = "GeneralisedCDOMC::payoff";
	try
	{
		const CreditLossConfigIndexedSVResultPath& results = indexedSV->getResults();

		double pv = 0;
		double contingentLegPV = 0;

//check whether discountFactor state variable got populated
		if (discFactorSV.get() == 0)
			throw ModelException(
				method,
				"Discount factor state variable not populated");

//computing the contingentLeg
		if (this->cdo->cLeg.get())
		{
			int i;
			for (i = results.begin(); i <= results.end(); ++i)
			{
				int productIndex = results[i].index;
				if (productIndex != -1) //credit event happened
				{
					if (payLoss[productIndex]) //put a check for whether contingent leg exists
					{
						int discFactorIndex = contLegDiscFactorIndex[productIndex];
						contingentLegPV += discFactorSV->getDF(discFactorIndex) *
											contingentLegScalingFactors[productIndex] *
											results[i].lossChange;
					}
				}
			}
		}

//computing the fee leg
		double feeLegPV = 0;
		double riskFreePV = 0;
		double unitCouponPV = 0;

		if (creditFeeLegSV.get())
			feeLegPV = creditFeeLegSV->pv(
											riskFreePV,
											unitCouponPV,
											true,
											cdo->today,
											results,
											*discFactorSV
											);

//aggregate the pvs of the two legs
		pv += contingentLegPV;
		trancheContingentLegPrice.add(contingentLegPV);

		pv -= feeLegPV;
		trancheFeeLegPrice.add(feeLegPV);
		trancheFeeLegRFPrice.add(riskFreePV);
		unitCouponTrancheFeeLegPrice.add(unitCouponPV);

		prices.add(pv);
	}
	catch (exception& e)
	{
		throw ModelException(e, method);

	}
}

/** stores extra output to be computed for the   */
void
GeneralisedCDO::GeneralisedCDOMC::recordExtraOutput(
		Control* control,
		Results* results,
		const IMCPrices& prices) const
{
	double stdErr;

//store contingent leg price
	double contingentLeg;
	trancheContingentLegPrice.getResult(&contingentLeg, &stdErr);
	OutputRequest* request =
		control->requestsOutput(OutputRequest::TRANCHE_CONTINGENT_LEG_PRICE);
    if (request)
		results->storeRequestResult(request, contingentLeg);


//store fee leg price
	double feeLeg;
	trancheFeeLegPrice.getResult(&feeLeg, &stdErr);
	request =
		control->requestsOutput(OutputRequest::TRANCHE_FEE_LEG_PRICE);
    if (request)
		results->storeRequestResult(request, feeLeg);

//calculate and store tranche implied spread
	double riskFreeCF;
	double unitCouponPV;
	request =
		control->requestsOutput(OutputRequest::TRANCHE_IMPLIED_SPREAD);
	if (request)
	{
		trancheFeeLegRFPrice.getResult(&riskFreeCF, &stdErr);
		unitCouponTrancheFeeLegPrice.getResult(&unitCouponPV, &stdErr);

		if (unitCouponPV != 0)
		{
			double impliedSpread = (contingentLeg - riskFreeCF) /unitCouponPV;
			results->storeRequestResult(request, impliedSpread);
		}
	}
}


DateTime
GeneralisedCDO::GeneralisedCDOMC::endDate(const Sensitivity*  sensControl) const
{
	return productTimeline->back();
}

DRLIB_END_NAMESPACE
