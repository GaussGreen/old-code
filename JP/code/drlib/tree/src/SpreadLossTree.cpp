//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : SpreadLossTree.cpp
//
//   Description : Tree for bivariate spread and portfolio loss process
//
//   Author      : Matthias Arnsdorf
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/TimeLine.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/SpreadLossTree.hpp"
#include "edginc/ConstAuxSpreadProcess.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"
#include <algorithm>
#include "edginc/Timer.hpp"

DRLIB_BEGIN_NAMESPACE

// flag for counting number of generator matrices created
//#define COUNT_GENERATORS
#ifdef COUNT_GENERATORS
static int GEN_COUNTER =0;
#endif

const double SpreadLossTree::DEFAULT_MIN_TRANS_PROB = 1e-12;
const int SpreadLossTree::DEFAULT_MAX_ORDER = 300;
const int SpreadLossTree::DEFAULT_PADE_ORDER = 4;
const double SpreadLossTree::DT_ACCURACY = 1e6;
const int SpreadLossTree::MAX_CACHE_SIZE = 100;

/** default constructor */
IObject* SpreadLossTree::defaultConstructor(void) { return new SpreadLossTree(); }

/** constructor */
SpreadLossTree::SpreadLossTree() : 
FDModel(TYPE),
spreadProcess(new ConstAuxSpreadProcess(0.01)),
discYC(0),
condDefTransProbs(),
minTransProb(DEFAULT_MIN_TRANS_PROB),
maxExpansionOrder(DEFAULT_MAX_ORDER),
padeOrder(DEFAULT_PADE_ORDER),
usePade(true),
maxLossForGenerator(-1), 
bucketStartLoss(1),
bucketSize(1),
useEffCurves(false),
maxStrike(0),
maxDefRangeSize(-1)
{
	// intialise variables
	elHighStrikes.resize(0);
	elLowStrikes.resize(0);
	elWithRRHighStrikes.resize(0);
	elWithRRLowStrikes.resize(0);
	trancheEL.resize(0);
	trancheELwithRR.resize(0);
	coeffs.resize(0);

	

}

/** destructor */
SpreadLossTree::~SpreadLossTree()
{	
	// clear cache
	contagionCache.clear();
	condDefProbCache.clear();
	numDefBrancheCache.clear();
};

void SpreadLossTree::load(CClassSP& clazz)
{
        REGISTER(SpreadLossTree, clazz);
        SUPERCLASS(FDModel);
		EMPTY_SHELL_METHOD(defaultConstructor);
        clazz->setPublic(); // make visible to EAS/spreadsheet

		FIELD(numNames, "Number of credits in the portfolio");

		FIELD(stepsPerYear, "Number of timeline steps per year");

		FIELD(maxNumDefs, "Maximum number of defaults per tree time period");
		FIELD(minTransProb, "Minimum transition probabiltiy [default = "+
			Format::toString(DEFAULT_MIN_TRANS_PROB)+ "]");

		

		FIELD(usePade, "TRUE: use Pade approximation. FALSE: use exapnsion. [default = TRUE]");
		FIELD_MAKE_OPTIONAL(usePade);

		FIELD(padeOrder, 
			"Order for Pade approximation, [default = " +
			Format::toString(DEFAULT_PADE_ORDER)+ "]");
		FIELD_MAKE_OPTIONAL(padeOrder);

		FIELD(maxLossForGenerator,
			"Max loss level to use for generator (to use with Pade in calibration). [default: high strike of current instrument]");
		FIELD_MAKE_OPTIONAL(maxLossForGenerator);

		FIELD(maxExpansionOrder, 
			"Maximum order for expansion of conditional default transition probabilities, [default = " +
			Format::toString(DEFAULT_MAX_ORDER)+ "]");
		FIELD_MAKE_OPTIONAL(maxExpansionOrder);

		FIELD_MAKE_OPTIONAL(minTransProb);

		FIELD(recoveryModel, "Model for the recovery rate");
		
		FIELD(spreadProcess, "Auxiliary spread process [default = none]");
		FIELD_MAKE_OPTIONAL(spreadProcess); 

        FIELD(curveToDiscount, "Discount curve");

		FIELD(localContagion, "Calibrated local contagion surface");

        FIELD(bucketStartLoss,"Loss level at which we start bucketing. [default = 1]");
        FIELD_MAKE_OPTIONAL(bucketStartLoss);

        FIELD(bucketSize, "Number of defaults per bucket. [default = 1]");
        FIELD_MAKE_OPTIONAL(bucketSize);

		FIELD(useEffCurves, "True if using effective curve pricing [default = false]");
		FIELD_MAKE_OPTIONAL(useEffCurves);

		// Transient fields
		FIELD(valueDate,"");
		FIELD_MAKE_TRANSIENT(valueDate);

		FIELD(discYC,"");
		FIELD_MAKE_TRANSIENT(discYC);
		
		FIELD(elLowStrikes,"");
		FIELD_MAKE_TRANSIENT(elLowStrikes);
		
		FIELD(elHighStrikes,"");
		FIELD_MAKE_TRANSIENT(elHighStrikes);

		FIELD(elWithRRLowStrikes,"");
		FIELD_MAKE_TRANSIENT(elWithRRLowStrikes);
		
		FIELD(elWithRRHighStrikes,"");
		FIELD_MAKE_TRANSIENT(elWithRRHighStrikes);

		FIELD(coeffs,"");
		FIELD_MAKE_TRANSIENT(coeffs);

        FIELD(defaults,"");
        FIELD_MAKE_TRANSIENT(defaults);

		FIELD(totalLosses,"");
		FIELD_MAKE_TRANSIENT(totalLosses);
		FIELD(totalRecRates,"");
		FIELD_MAKE_TRANSIENT(totalRecRates);
		
		FIELD(contagionDates,"");
		FIELD_MAKE_TRANSIENT(contagionDates);

		FIELD(spreadDates,"");
		FIELD_MAKE_TRANSIENT(spreadDates);

		FIELD(maxDefsForCache,"");
		FIELD_MAKE_TRANSIENT(maxDefsForCache);

		FIELD(maxStrike,"");
		FIELD_MAKE_TRANSIENT(maxStrike);

        FIELD(maxDefIdxForCache,"");
        FIELD_MAKE_TRANSIENT(maxDefIdxForCache);

        FIELD(bucketStart,"");
        FIELD_MAKE_TRANSIENT(bucketStart);

        FIELD(spreadRangeBot,"");
        FIELD_MAKE_TRANSIENT(spreadRangeBot);

        FIELD(spreadRangeTop,"");
        FIELD_MAKE_TRANSIENT(spreadRangeTop);

        FIELD(defaultRangeBot,"");
        FIELD_MAKE_TRANSIENT(defaultRangeBot);

        FIELD(defaultRangeTop,"");
        FIELD_MAKE_TRANSIENT(defaultRangeTop);


}

CClassConstSP const SpreadLossTree::TYPE = CClass::registerClassLoadMethod(
    "SpreadLossTree", typeid(SpreadLossTree), load);



/************************************************************************/
/* check inputs and do initialisation of variables that depend on
	model only (not instruments which drive the timline)                */
/************************************************************************/
void SpreadLossTree::validatePop2Object()
{
	//TODO get rid of maxNumDefs as input
	if(true)
	{
		// allow all names to default in single time period
		maxNumDefs = numNames;
	}
	else
	{
		// bound checks on maxNumDefs
		if(maxNumDefs < 1) 
		{
			throw ModelException("maxNumDefs ("+Format::toString(maxNumDefs)+") < 1");
		}
		if(maxNumDefs > numNames)
		{
			throw ModelException("maxNumDefs ("+Format::toString(maxNumDefs)+") > numNames ("
				+Format::toString(numNames)+ ")");
		}

	}
	

	// check that time interpolation for localContagion is flat (since need this for generator caching) */
	if(localContagion->getMatInterpType() != SkewSurface::FLAT_1D_INTERPOLATION)
	{
		throw ModelException("matInterpolationType in the localContagion surface needs to be FLAT");
	}

    if(bucketStartLoss > 1 || bucketStartLoss < 0)
    {
        throw ModelException("bucketStartLoss ("+Format::toString(bucketStartLoss) +
            ") is out of bounds [0, 1]");
    }

    // floor bucket size
    bucketSize = Maths::max(1, bucketSize); 
  	
	// now can setup product independent model data
	setupModel();

}


/************************************************************************/
/* setup product independent model parameters							*/
/************************************************************************/
void SpreadLossTree::setupModel()
{
	static const string method = "SpreadLossTree::setupModel";
	try
	{
		
		// intialise recovery model
		recoveryModel->setNumNames(numNames);

		// get recovery arrays
		totalLosses	= recoveryModel->getTotalLosses();
		totalRecRates =  recoveryModel->getTotalRecoveries();
		//recRates = recoveryModel->getRecoveries();

		// set up contagion dates 
		contagionDates = *localContagion->getMaturities();

		// remove duplicates and sort (note that maturities in skew surface are per quote not per unique maturity)
		DateTime::doSortUniq(contagionDates);

		// set up maxDefsForCache
		if(maxLossForGenerator > 1)
		{
			throw ModelException("maxLossForGenerator (" +
				Format::toString(maxLossForGenerator) +") > 1");
		}

		//set maxDefsForCache
		if(maxLossForGenerator < 0.0)
		{
			maxDefsForCache = -1;
		}
		else
		{
			// find first default level giving loss >= maxLossForGenerator
			maxDefsForCache = 0;
			while(maxDefsForCache < numNames-1 && (*totalLosses)[maxDefsForCache] < maxLossForGenerator)
				maxDefsForCache++;
		}

        // find number of defaults corresponding to bucketStartLoss
        bucketStart = 0;
        while(bucketStart < numNames -1 && (*totalLosses)[bucketStart] < bucketStartLoss )
            bucketStart++;

        // set up the slice of underlying default values
        defaults.resize(bucketStart+1);
        // populate the possible default levels
        // defaults levels start at 0
        // default levels move in steps of 1 up to bucket start
        // after that they move in steps of 'bucketSize'
        int i;
        for(i  = 0; i <= bucketStart; i++) // note that bucketStart <= numNames
        {
            defaults[i] = i;
        }

        // now start the buckets
        for(i = bucketStart+1 ; defaults[i-1] + bucketSize < numNames ; i++)
        {
            defaults.push_back(defaults[i-1] + bucketSize);
        }
        // final point if needed
        if(defaults.back() < numNames)
            defaults.push_back(numNames);

        
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}

}


/** Invoked after instrument has got its market data. */
void SpreadLossTree::getMarket(const MarketData* market, IInstrumentCollectionSP instruments)
{
	static const string method = "SpreadLossTree::getMarket";
	try 
	{
		// intialise cache
		contagionCache.clear();
		condDefProbCache.clear();
		numDefBrancheCache.clear();

		// populate value date
		valueDate = market->GetReferenceDate();
		
		// populate discount curve
		curveToDiscount.getData(this, market);

		// discount yield curve
		discYC = curveToDiscount.getSP();	

		// pass market data to spreadProcess
		spreadProcess->getMarketData(market, stepsPerYear);

		// get the dates at which the spread slice values change
		spreadDates = spreadProcess->getSpreadDates();
		

	}
	catch (exception& e) {
		throw ModelException(e, method);
	}
}




/** retrieving market data */
void SpreadLossTree::retrieveFactor(){
	static const string method = "SpreadLossTree::retrieveFactor";
	try{
		
	}
	catch (exception& e){
		throw ModelException(e, method);
	}
}


/************************************************************************/
/*   create range for 2d spread - default slice                         */
/************************************************************************/
TreeSliceGeneral::RangeSP SpreadLossTree::createRange2D(
	int bot1, int top1,
	int bot2, int top2
	)
{
	TreeSliceGeneral::RangeSP outRange = TreeSliceGeneral::Range::create(
		1, // count
		bot1,	top1, 
		bot2,	top2 
		);

	
	return outRange;
}


// FDModel INTERFACE //////////////////////////////////////////////

/** creates new 2-d slice */
TreeSliceSP SpreadLossTree::createSlice(
    const string & curveToDEV,
    const string & factorName1,
    const string & factorName2,
    const string & factorName3 ) const
{
	int dimBits = 3; // TODO: check this is correct for 2-d slice
	return TreeSliceGeneral::create(*range, dimBits);
}


/** creates new  generalSlice */
TreeSliceGeneralSP SpreadLossTree::createGeneralSlice() const
{
	return DYNAMIC_POINTER_CAST<TreeSliceGeneral>(createSlice());
}



/** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
//TODO get rid of this once have moved to TreeSliceGeneral
double SpreadLossTree::getPrice0( const TreeSlice & price ) const
{	
	const TreeSliceGeneral & p = static_cast< const TreeSliceGeneral & >( price );
	return p(0,0);
}


/** return slice of total losses for each node for maximum loss range*/
TreeSliceSP SpreadLossTree::getLossSlice()
{
	TreeSliceGeneralSP lossSlice = createGeneralSlice();

	int d,s;
    QLIB_VERIFY(maxDefRangeSize >=0,"maxDefRangeSize is not initialised");
	// loop through defaults
	for(d = 0; d < maxDefRangeSize ;d++)
	{
		// loop through spreads
		for(s = spreadRangeBot[currStep]; s <= spreadRangeTop[currStep]; s++)
		{
			(*lossSlice)(s,d) = (*totalLosses)[defaults[d]];
		}
	}

	return lossSlice;

}

/** return slice of total recovery for each node */
TreeSliceSP  SpreadLossTree::getRecoverySlice()
{
	TreeSliceGeneralSP rrSlice = createGeneralSlice();

	int d,s;
	// loop through defaults
    QLIB_VERIFY(maxDefRangeSize >=0,"maxDefRangeSize is not initialised");
	for(d = 0; d < maxDefRangeSize ;d++)
	{
		//loop through spreads
		for(s = spreadRangeBot[currStep]; s <= spreadRangeTop[currStep]; s++)
		{
			(*rrSlice)(s,d) = (*totalRecRates)[defaults[d]];
		}
	}

	return rrSlice;

}

/** Get slice representing zero/discounter between useDate and matDate 
*  If "slice" is null it creates a new slice (no ownership kept), 
*  otherwise it uses the slices passed to return the zero */
void SpreadLossTree::getZero(
					 DateTime useDate, 
					 DateTime matDate, 
					 string curveName, 
					 TreeSliceSP &slice) 
{
	static const string method = "SpreadLossTree::getZero";
	
	if(curveName != discYC->getName())
	{
		throw ModelException("Only support discounting by " + discYC->getName(), method);
	}

	double df = discYC->pv(useDate, matDate);

	if(!slice)
	{
		slice = createSlice();
	}
	*slice = df;

					 
}

/************************************************************************/
/*   collect model initialization data, set up timeline                 */
// override FDModel::initModel() so that can have own timeline def                                                                   */
/************************************************************************/
void SpreadLossTree::initModel()
{

	static const string method = "SpreadLossTree::initModel()";
	try
	{
		// check cache size -------------------------------------------------------------

		if(condDefProbCache.size() > MAX_CACHE_SIZE)
		{
			condDefProbCache.clear();
			numDefBrancheCache.clear();
			contagionCache.clear();
		}
        // The two checks belwo are for sanity ,they should not occur
		if(numDefBrancheCache.size() > MAX_CACHE_SIZE)
		{
			throw ModelException("numDefBrancheCache size exceeds max cache size");	
		}
		if(contagionCache.size() > MAX_CACHE_SIZE)
		{
			throw ModelException("contagionCache size exceeds max cache size");		
		}

	

		// reset strikes
		elLowStrikes.resize(0);
		elHighStrikes.resize(0);
		elWithRRLowStrikes.resize(0);
		elWithRRHighStrikes.resize(0);

		// set up model timeline --------------------------------------------------------

		// set timeMetric in FDModel (needed to build timeline) */
		timeMetric = TimeMetricSP(new TimeMetric(1.0, Holiday::noHolidays()));

		// set number of steps per year in FDModel - this is what will drive the timeline construction
		stepsPerYearFD = stepsPerYear;
		
        // add valueDate to crit dates in case not already included
        addCritDate(valueDate);

		// intialise segments in FDModels
		// timeline dates are allocated within the segment such that the total number
		// of steps per year is stepsPerYearFD
		// note however that there is a minimum of 3 dates per segment
		initSegments( 
            critDates,
			IntArray(0)		// density 
			);
		

		// call implementation on FDModel to setup timeline 
		FDModel::initModel();

		// setup spread Process ////////////////////////////////////////////////////////////
		spreadProcess->setupModel(timeLine);


		//** number of dates in timeLine */
		int numDates = getLastStep() + 1;

		// currStep holds the index corresponding to the current time step
		// Initialise to final step
		currStep = getLastStep();

		// set up ranges /////////////////////////////////////////////////////////////////////

		/// SPREAD RANGE -----------------------------------------

		spreadRangeTop = spreadProcess->getTopLimits();
		spreadRangeBot = spreadProcess->getBotLimits();

		if(spreadRangeTop.size() != numDates)
		{
			throw ModelException("spreadRangeTop has incorrect size");
		}
		if(spreadRangeBot.size() != numDates)
		{
			throw ModelException("spreadRangeBot has incorrect size");
		}
		

		// set lowest range at final step
		minSpreadRangeBot = spreadRangeBot.back();

		// maximumum number of spread nodes
		maxSpreadSize = spreadRangeTop.back() - minSpreadRangeBot +1; //(*spreadProcess->getRange())->size[0];

		// number of branches in spread tree
		numSpreadBranches = spreadProcess->getNumBranches();


		/// DEFAULTS RANGE ---------------------------------------------
		if(maxStrike == 0) 
		{
			throw ModelException("maxStrike has not been set");
		}

		/** number of defaults that correspond to first loss >= maxStrike */
		int defBound= 0;
		while(defBound < numNames && (*totalLosses)[defBound] < maxStrike) defBound++;

		// can now clear maxStrike for next pricing call
		maxStrike = 0;

		// floor maxDefsForCache by defBound 
		maxDefsForCache = Maths::max(maxDefsForCache, defBound);

        
        // corresponding index such that defaults[maxDefIdxForCache] >= maxDefsForCache;
        if(maxDefsForCache <= bucketStart)
        {
            maxDefIdxForCache = maxDefsForCache;
        }
        else
        {
            for(maxDefIdxForCache = bucketStart; 
                maxDefIdxForCache < numNames && defaults[maxDefIdxForCache] < maxDefsForCache; 
                maxDefIdxForCache++);
        }

		
        // maximum top index for defaults range
        int finalRangeIndex;
        if(defBound <= bucketStart)
        {
            finalRangeIndex = defBound;
        }
        else
        {
            // find first index such that defaults[index] >= defBound;
            for(finalRangeIndex = bucketStart; 
                finalRangeIndex < numNames && defaults[finalRangeIndex] < defBound; 
                finalRangeIndex++);
        }

		maxDefRangeSize = finalRangeIndex + 1;


		// 2-dim SPREAD / DEFAULTS  RANGE -----------------------------------------

		/** 2-d range for full tree slices */
		// This sets the number of nodes in tree for each dimension and each time slice
		// bot/top refer to the bottom and top indicies in the treeSlice.
		range = createRange2D(
			spreadRangeBot.back(), 
			spreadRangeTop.back(), 
			0,
			finalRangeIndex);


	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

/************************************************************************/
/*  finalize model initialization										*/
// called after product initialisation
/************************************************************************/
void SpreadLossTree::finaliseModel(CControl*    control)
{
	static const string method ="SpreadLossTree::finaliseModel";

	int i,t;

    //don't want to muck around with FD models
    //but this limits the type of cashflows supportable
    //by the tree
    ClosedFormForwardRatePricerSP fwdModel =
        ClosedFormForwardRatePricerSP(
            new ClosedFormForwardRatePricer());


	
	int numDates = currStep +1;

	// set up temp slice which can be used to store temporary slice values
	tempSlice = createGeneralSlice();



	// if using effective curve pricing -----------------------------------------------------
	// set up expected loss containers if contingent or fee legs have been inserted
	// allocates maximum amount of memory that is needed at start
	if(useEffCurves)
	{
		if(elLowStrikes.size() == 0)
		{
			throw ModelException("No strikes have been added for effective curve pricing");
		}
		if(elLowStrikes.size() != elHighStrikes.size())
		{
			throw ModelException("low and high strikes have different size", method);
		}
		if(elWithRRLowStrikes.size() != elWithRRHighStrikes.size())
		{
			throw ModelException("low and high strikes with RR have different size", method);
		}

		trancheEL.resize(elLowStrikes.size());

		for(i = 0; i < (int) trancheEL.size(); i++)
		{
			trancheEL[i].resize(numDates);
			for(t = 0 ;t < numDates; t++)
			{
				// for each tranche and each time point hold tree slice
				// TODO make this only dependent on default states
				trancheEL[i][t] = createGeneralSlice(); 
			}
		}	



		trancheELwithRR.resize(elWithRRLowStrikes.size());

		trancheELwithRR.resize(elWithRRLowStrikes.size());
		for(i = 0; i < (int) trancheELwithRR.size(); i++)
		{
			trancheELwithRR[i].resize(numDates);
			for(t = 0 ;t < numDates; t++)
			{
				// for each tranche and each time point hold tree slice
				// TODO make this only dependent on default states
				trancheELwithRR[i][t] = createGeneralSlice(); 
			}
		}
	}
	
	
	// set up coeff matrix to max size (used for expansion )
	if(!usePade)
		coeffs.resize((maxNumDefs+2)*maxExpansionOrder);

	// contagion date projection gives for each timeline point the index in contagionDates
	// of the contagionDate that is just greater than the current timeLine date
	contagionDateProjection = DateTime::getCeilingProjection(getDates(), contagionDates);

	// spread date projection gives for each timeline point the index in spreadDates
	// of the spreadDate that is just greater than the current timeLine date
	spreadDateProjection = DateTime::getCeilingProjection(getDates(), spreadDates);


    // CONDITIONAL QUANTITITES ---------------------------------------------

    // clear cache
    stateProbCache.clear();
    marginalDefProbCache.clear();

    

    // calc and store conditional trans probs --------------------------------
    // intialise lastKey 
    lastKey = -1;

    // define bottom range for defaults (needed when calculating condDefTransProbs)
    // bottom of default range is always 0
    defaultRangeBot.resize(numDates,0);
    
    // backwards time loop store keys in keyArray;
    IntArray keyArray(numDates);
    for(t = getLastStep(); t >=0 ; t--)
    {
        spreadProcess->updateSpread(t, FDProduct::BWD);
        keyArray[t] = updateCondProbs(t, FDProduct::BWD);
    }
    
    // do tree clipping in default dimension -------------------------------
   
    defaultRangeTop.resize(numDates);
    defaultRangeTop[0] = 0;
    int s, maxBranches;
    for(t = 1 ; t <= getLastStep() ; t++)
    {
        int prevTop = defaultRangeTop[t-1];

        // check if we are already at maximum level
        if(prevTop < maxDefRangeSize)
        {
            // need to get maxNumDefBranches from last step
            IntArrayArraySP numDefBranches = numDefBrancheCache[keyArray[t-1]];

            maxBranches = 0;
            for(s = 0; s < numDefBranches->size(); s++ )
            {
                if((*(*numDefBranches)[s])[prevTop] > maxBranches)
                    maxBranches = (*(*numDefBranches)[s])[prevTop];
            }

        }
        defaultRangeTop[t] = Maths::min(prevTop + maxBranches, maxDefRangeSize-1);
    }
    

}


/************************************************************************/
/*  Solver Loop: Update tree by sweeping timeline                       */                                        
/************************************************************************/
void SpreadLossTree::roll()
{
	int t; // counter
	int T = getLastStep();
	string prodType, prodOutputName;

	try {

        
		const FDProductArray & products = getProducts();  //for easy reference
		FDProduct::UpdateType type = FDProduct::BWD_T;

		//Timer timer; // TESTING

		// loop backwards through timeline
		for (t = T; t >= 0; t--) 
		{

		// expand slices to the correct dimension for dev
			// TODO do we need this????
			for(size_t k=0; k < products.size(); k++) 
            {
				if (products[k]->isCalcOff())
					continue;

				const vector< TreeSliceSP > & slices = products[k]->getSlicesToDEV();
				for (int j=0; j<(int)slices.size(); j++) {
					//sliceExpand(*slices[j], slices[j]->getCurveToDEV());
					sliceExpand(*slices[j], discYC->getName());
				}
			}

			// update engine prob's and internal variables
			update(t, type);

			// for each product
			for(size_t k=0; k < products.size(); k++) 
            {
				if (products[k]->isCalcOff())
					continue;

				// to assist debugging print out type and name if available of current product
				prodType = typeid((*products[k])).name();
				prodOutputName = products[k]->getOutputName();

				// product with no slices still calls preCalc() and update()
				products[k]->preCalc(t);

				// DEV slices
				const vector< TreeSliceSP > & slices = products[k]->getSlicesToDEV();
				for (int j=0; j<(int)slices.size(); j++) {
					//sliceDev(*slices[j], slices[j]->getCurveToDEV());
					sliceDev(*slices[j], discYC->getName());
				}
			

				// now call update on product
				products[k]->update(t, type);
			}
			updateFinalize(t, type);
			type = FDProduct::BWD;
		} // end time loop

		//TESTING
		//double timeVal = timer.calcTime();
		//ErrorHandler::writeMsg("Roll time: " +Format::toString(timeVal) );
		//cout << "Roll time: " << timeVal << endl;

#ifdef COUNT_GENERATORS
        ErrorHandler::writeMsg("Num generators: " +Format::toString(GEN_COUNTER) );
        cout << "Num generators: " << GEN_COUNTER << endl;
#endif

        
	}
	catch (exception& e) {
		throw ModelException(e, "SpreadLossTree::Solver::roll failed for product type \"" + 
			prodType + "\", product outputName \"" +
			prodOutputName + "\" on timestep " +
			Format::toString(t) + " of " +
			Format::toString(T)+ " (model is "+getClass()->getName()+")");
	}

}

/************************************************************************/
/*  Update tree variables at the "top" of the time loop                 */
/************************************************************************/
void SpreadLossTree::update(int step, FDProduct::UpdateType type)
{
	static const string method = "SpreadLossTree::update";
	try
	{
		// set time index in slice to current step
		currStep = step;

		// update spread process
		spreadProcess->update(step, type); 

		// update conditional default transition probs
        //TODO improve this to use keyArray
        updateCondProbs(currStep, type);
        
		if(useEffCurves)
		{
            if(type != FDProduct::FWD) // only support BWD induction for now
            {
			    // update expected losses if using effective curve pricing
			    updateExpLoss();
		    }
        }

	}
	catch (exception& e){
		throw ModelException(e, method);
	}

}

/************************************************************************/
/*  updates the conditional transition probabilities and returns the key */
/* assumes that spreads have been updated                                */
/************************************************************************/
int SpreadLossTree::updateCondProbs(int timeIdx, FDProduct::UpdateType type)
{
    static const string method = "SpreadLossTree::updateCondProbs";
    try
    {
        if(timeIdx == getLastStep()) return -1; // no probs for last step
        
        int i,j;
        // check if condDefProbs are in cache, update otherwise

        /** year frac between current time step and next time step on timeline */
        double dt = timeLine->TradeYrFrac[timeIdx +1];

        /** key for condDefProbCache */
        int key = ((int)(dt*DT_ACCURACY)) * 
            (contagionDateProjection[timeIdx]+1)*
            (spreadDateProjection[timeIdx] + 7); // TODO check this is unique


        if(key < 0)
        {
            throw ModelException("key <0");
        }
        // check if key has changed from last update */
        bool haveNewKey = (key != lastKey);
        lastKey = key;

        // if key has not changed no more needs to be done. Otherwise need to check cache
        if(haveNewKey)
        {
            /** contagion values at contagionDate */
            DoubleArraySP contagionVals(new DoubleArray(0));
            // get relevant contagion strikes and value
            DoubleArray dummyStrikes(0); // not used
            localContagion->getStrikesAndSkews(
                dummyStrikes, 
                *contagionVals, 
                contagionDates[ contagionDateProjection[timeIdx] ]
                );

                // append number of possible defaults (2nd and 3d dim of condDefTransPorbs)
                contagionVals->push_back((double) maxDefIdxForCache);

                /** flag to indicate if we need to recalculate transition probs */
                // need to recalc if either
                // 1) there is no entry in cache for the key
                // 2) the current contagion vals are different to the ones corresponding to the cache entry
                bool recalcProbs = true;

                /** flag to indicate whether we need new memory allocation */
                // need new memory if either
                // 1) there is no entry in cache for the key.
                // 2) The size of the stored memory is less than what is needed now.
                //	  This can be if maxDefsForCache has changed between pricings.
                bool needNewMemory = true;

                // set the flags
                if(	condDefProbCache.find(key) != condDefProbCache.end() )
                {
                    if( *contagionCache[key] == *contagionVals	)
                    {
                        recalcProbs = false;
                        needNewMemory = false;
                    }
                    else
                    {
                        // check if default dim is current maxDefsForCache+1
                        // this check is needed in case maxLossForGenerator was not set or set wrongly
                        if(contagionCache[key]->back() == maxDefIdxForCache)
                        {
                            needNewMemory = false;
                        }
                    }
                }


                if(needNewMemory)
                {
                    // need to allocate new memory 
                    // (*condDefTransProbs)[s][n][d1] is prob of d1+n defaults at t+dt given 

                    // size of spread dimension at timeIdx
                    int spreadDimSize = spreadRangeTop[timeIdx] - spreadRangeBot[timeIdx] +1;
                    // note that this only works if we are moving bwds in tree since then the
                    // spread tree is contracting. This means that the first time we allocate memory
                    // will also have the maximum spread dim for this spread level
                    QLIB_VERIFY(type != FDProduct::FWD, "Cannot update probs while moving fwds in tree");

                    // spread level s and default level d1 at t
                    condDefTransProbs = DoubleArrayArrayArraySP(new DoubleArrayArrayArray(spreadDimSize));
                    for(i = 0; i < spreadDimSize; i++)
                    {
                        (*condDefTransProbs)[i].resize(maxDefIdxForCache+1);

                        for(j = 0 ;j< maxDefIdxForCache +1; j++)
                            (*condDefTransProbs)[i][j].resize(maxDefIdxForCache + 1);
                    }
                    condDefProbCache[key] = condDefTransProbs;

                    // set up number of default branches
                    numDefaultBranches = IntArrayArraySP(new IntArrayArray(spreadDimSize));
                    for(i=0; i< spreadDimSize; i++)
                    {
                        (*numDefaultBranches)[i] = IntArraySP(new IntArray(maxDefIdxForCache+1));
                    }

                    numDefBrancheCache[key] = numDefaultBranches;
                }
                else // can just use memory already allocated to cache
                {
                    condDefTransProbs = condDefProbCache[key];
                    numDefaultBranches = numDefBrancheCache[key];

                }

                if(recalcProbs)
                {
                    // do new calculation of condDefTransProbs
                    if(usePade)
                    {
                        updateCondDefTransProb_Pade(timeIdx, condDefTransProbs, numDefaultBranches);
                    }
                    else
                    {
                        updateCondDefTransProb_expansion(timeIdx, condDefTransProbs, numDefaultBranches);
                    }

                    // update contagion val cache
                    contagionCache[key] = contagionVals;

                }

        }

#ifdef DEBUG
        // sanity check that condDefTransProbs have been set
        if(!condDefTransProbs)
        {
            throw ModelException("condDefTransProbs are NULL");
        }
#endif

       return key;
       
    }
    catch (exception & e)
    {
        throw ModelException(e, method);
    }
}


/************************************************************************/
/*  update expected loss curves                                         */
/************************************************************************/
void SpreadLossTree::updateExpLoss()
{
	static const string method = "SpreadLossTree:updateExpLoss";

	
	// EV previous expected losses 
	int tIdx, i;
	for(tIdx = currStep+1; tIdx <= getLastStep(); tIdx++)
	{
		//loop through tranches with no RR write down
		for(i = 0; i < (int) trancheEL.size() ; i++)
		{
			sliceEv(*trancheEL[i][tIdx]);
		}
		//loop through tranches with RR write down
		for(i = 0; i < (int) trancheELwithRR.size() ; i++)
		{
			sliceEv(*trancheELwithRR[i][tIdx]);
		}
	}

	int s,d;
	double trancheLoss;
	double nodeLossVal;
	
	// add expected loss for current slice
	//loop through tranches with no RR write down
	for(i = 0; i < (int) trancheEL.size() ; i++)
	{

		/** pointer to TreeSlice values */
		//double * pTrancheELSlice =  getValPtr(*trancheEL[i][currStep]);
		/** expected loss slice for ease */
		TreeSliceGeneral& trancheELSlice = *trancheEL[i][currStep];

		// loop through defaults at currStep
		for(d = defaultRangeBot[currStep]; d <= defaultRangeTop[currStep] ;d++)
		{
			nodeLossVal = (*totalLosses)[defaults[d]];

			// calc tranche loss
			trancheLoss =  
				Maths::min(nodeLossVal, elHighStrikes[i]) 
				- Maths::min(nodeLossVal, elLowStrikes[i]);

			// loop through spreads and allocate trancheLoss
			for(s = spreadRangeBot[currStep]; s <= spreadRangeTop[currStep]; s++)
			{
				trancheELSlice(s,d) = trancheLoss;
			}
		} // end default loop
	} // end tranche loop

	double nodeRRVal;
	//loop through tranches with  RR write down
	for(i = 0; i < (int) trancheELwithRR.size() ; i++)
	{
		/** pointer to TreeSlice values */
		//double * pTrancheELSlice =  getValPtr(*trancheELwithRR[i][currStep]);

		TreeSliceGeneral& trancheELSlice = *trancheELwithRR[i][currStep];

		// loop through defaults
		for(d = defaultRangeBot[currStep]; d <= defaultRangeTop[currStep] ;d++)
		{
			nodeLossVal = (*totalLosses)[defaults[d]];
			nodeRRVal	= (*totalRecRates)[defaults[d]];

			// calc tranche loss
			trancheLoss =  
				Maths::min(nodeLossVal, elWithRRHighStrikes[i]) 
				- Maths::min(nodeLossVal, elWithRRLowStrikes[i])
				+ Maths::min(nodeRRVal, 1 - elWithRRLowStrikes[i]) 
				- Maths::min(nodeRRVal, 1 - elWithRRHighStrikes[i]);

			// loop through spreads and allocate trancheLoss
			for(s = spreadRangeBot[currStep]; s <= spreadRangeTop[currStep]; s++)
			{
				trancheELSlice(s,d) = trancheLoss;
			}
		} // end default loop
	} // end tranche loop

}



/************************************************************************/
/*  update contingent legs		                                        */
/************************************************************************/
void SpreadLossTree::updateContLeg(TreeSlice& valueSlice,
								   double lowStrike,
								   double highStrike)
{

	static const string method ="SpreadLossTree::updateContLeg";
	try
	{
		// need to ev the difference in losses between currStep and nextStep
		// assume that default happens in the middle of the period
		// for now we assume that we have protection over the entire model timeline
		// TODO change this

		// nothing to do at the back of the tree 
		if (currStep == getLastStep()) return;

		/** size of spread offset in condDefProbs */
		int spreadOffset = -(condDefTransProbs->size()-1)/2;

#ifdef DUBUG
        // check that spreadOffset is correct
        QLIB_VERIFY(spreadRangeBot[currStep] - spreadOffset >= 0, "spreadRangeBot out of bounds");
        QLIB_VERIFY(spreadRangeTop[currStep] - spreadOffset < condDefTransProbs.size(), 
            "spreadRangeTop out of bounds");
#endif

		// Do DEV on existing fee leg slice-------------------------


		// get discount factor to middle of period
		DateTime currDate = timeLine->StepDates[currStep];
		double df = discYC->pv(currDate, 
			currDate.midPoint(timeLine->StepDates[currStep+1]));
	
		// counters
		int l;	
		int sIdx, dIdx;


		// cast to general slice
		TreeSliceGeneral& contLegSlice = dynamic_cast<TreeSliceGeneral&>(valueSlice);
		// TODO check pointer


		/** tranche loss at currStep */
		double currTrancheLoss;

		/** tranche loss at next step */
		double nextTrancheLoss;
		
		/** expected value */
		double expVal = 0;

		// for ease
		double trancheSize = highStrike - lowStrike;
		
		double nodeLossVal;

		// loop through slice ------------------------------------------------
		
		int maxDefBranches;
		double condProb;
		double lastProb;
		// loop through  defaults at currStep
		for(dIdx = defaultRangeBot[currStep]; dIdx <= defaultRangeTop[currStep]; dIdx++)
		{	

			// calc tranche loss at currStep
			nodeLossVal = (*totalLosses)[defaults[dIdx]];
			
			currTrancheLoss =  
				Maths::min(nodeLossVal, highStrike) 
				- Maths::min(nodeLossVal, lowStrike);
			
			
			// loop through  spread at currStep
			for(sIdx = spreadRangeBot[currStep]; sIdx <= spreadRangeTop[currStep]; sIdx++)
			{
				maxDefBranches = Maths::min(
					(*(*numDefaultBranches)[sIdx-spreadOffset])[dIdx], defaultRangeTop[currStep+1] - dIdx + 1);

				expVal = 0;
				
				// loop through number of possible defaults
				lastProb = 1;
				for(l = 0; l < maxDefBranches;l++)
				{
					// calc tranche loss at next step
					// TODO could calculate this upfront 
					nodeLossVal = (*totalLosses)[defaults[dIdx + l]];
		
					nextTrancheLoss =  
						Maths::min(nodeLossVal, highStrike) 
						- Maths::min(nodeLossVal, lowStrike);
		
					if(l < maxDefBranches -1)
					{
						condProb = (*condDefTransProbs)[sIdx-spreadOffset][l][dIdx];
						lastProb -= condProb;
					}
					else
					{
						condProb = lastProb;
					}
					expVal += 
						condProb * 
						(nextTrancheLoss - currTrancheLoss);


				} // end outer default loop

				// update the slice value
				contLegSlice(sIdx,dIdx) += df*expVal / trancheSize;

			} // end inner spread loop
		} // end inner default loop

	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}

}


/** Implement any updating logic that must occur at the "bottom" of the time loop after
all product updates */
void SpreadLossTree::updateFinalize(int t, FDProduct::UpdateType type)
{
	//TODO
}

/** Expand the dimension of the slice so that it can be DEV'd. */
void SpreadLossTree::sliceExpand(TreeSlice& slice, const string& curveName) const
{
	try {
		TreeSliceGeneral& sliceGeneral = dynamic_cast<TreeSliceGeneral&>(slice);
		if(sliceGeneral.getDim() != 2)
		{
			int dimBits = 3; // TODO : check this
			sliceGeneral.expand(dimBits);
		}
	
	}
	catch (exception& e) {
		throw ModelException(e, "SpreadLossTree::sliceExpand");
	}
}



/** same as sliceDEV with df = 1 */
void SpreadLossTree::sliceEv(TreeSlice& slice, FDProduct::UpdateType type) const
{
	sliceDev(slice, "NONE", type);
}


/************************************************************************/
// Calculates discounted expected value of slice.                           
// Can be used in fwd or bkwd direction. In Bckwd direction we obtain      
// the standard DEV from                                                   
// t+1 to t. In fwd case we obtain the values of the slice at t+1 given       
// the slice values at t, ie the the value at node (sIdx, dIdx) at t is  
// given by the weighted average of all values at t with branches that    
// end at (sIDx, dIdx). The weighting is df^{-1} (prob of branch).                                                                     */
/************************************************************************/ 
void SpreadLossTree::sliceDev(TreeSlice& slice, 
                              const string& curveName, 
                              FDProduct::UpdateType type) const
{
	static const string method = "SpreadLossTree::sliceDev";
	try
	{
	    /** are we using fwd induction */
        bool isFwd = (type == FDProduct::FWD);

        // cast slice to GeneralTreeSlice so that can access values
        TreeSliceGeneral & sliceGeneral = dynamic_cast< TreeSliceGeneral & >(slice);

		// nothing to do at the back of the tree if are using bwd induction
		if (currStep == getLastStep() ) return;
		
        // if using fwd induction then only need to intialise at front of tree
        if(currStep == 0 && isFwd)
        {
            sliceGeneral(0,0) = 1.0;    
        }
		
        /** size of spread offset in condDefProbs */
		int spreadOffset = -(condDefTransProbs->size()-1)/2;
#ifdef DUBUG
        // check that spreadOffset is correct
        QLIB_VERIFY(spreadRangeBot[currStep] - spreadOffset >= 0, "spreadRangeBot out of bounds");
        QLIB_VERIFY(spreadRangeTop[currStep] - spreadOffset < condDefTransProbs.size(), 
            "spreadRangeTop out of bounds");
#endif

		// get discount factor
		double df;
		if(curveName == "NONE")
		{
			df = 1.0;
		}
		else if(curveName != discYC->getName())
		{
			throw ModelException("Only support discounting by " + discYC->getName());
		}
		else
		{
			df = discYC->pv(getDate(currStep), getDate(currStep+1));
		}
        QLIB_VERIFY(df != 0.0, "df = 0" );
		
		// counters for spread and default dimensions
		int sIdx, dIdx;	


        /** lastStep = currStep + 1 for bwd induction and currStep for fwd induction */
        int lastStep = isFwd ? currStep  : currStep + 1;
		
        
		// copy slice values at lastStep---------------------------------------
		// loop through second dim: loss
		for(dIdx = defaultRangeBot[lastStep]; dIdx <= defaultRangeTop[lastStep]; dIdx++)
		{
			// loop through first dim: spread 
			for(sIdx = spreadRangeBot[lastStep]; sIdx <= spreadRangeTop[lastStep]; sIdx++)
			{
				(*tempSlice)(sIdx, dIdx) = sliceGeneral(sIdx, dIdx); 
			} 
		}

        // in fwd case need to now initialise sliceGeneral to 0
        if(isFwd)
        {
            sliceGeneral = 0;
            sliceExpand(sliceGeneral, discYC->getName());
        }
        /** conditional default transition prob */
        double condProb;
        

        /** expected value */
        double evValue;

		// get spread process parameters ----------------------------
		DoubleArrayArrayArraySP spreadTransProbs = spreadProcess->getTransitionProbs();
		int maxContagionDefs = (*spreadTransProbs)[0].size()-1;
		
		IntArrayArraySP spreadOuterBotLimits = spreadProcess->getOuterRangeBot();

		// DEV the slice -------------------------------------------
        // we loop through the slice at timeStep where the branches originate
    
		// loop through second dim: loss at currStep
		for(dIdx = defaultRangeBot[currStep]; dIdx <= defaultRangeTop[currStep]; dIdx++)
		{
			// loop through first dim: spread at currStep
			for(sIdx = spreadRangeBot[currStep]; sIdx <= spreadRangeTop[currStep]; sIdx++)
			{
				// outer loop over values at nextStep ---------------------
                evValue = 0;
				
                /** maximum number of default branches at this node */
				int maxDefBranches = Maths::min(
					(*(*numDefaultBranches)[sIdx-spreadOffset])[dIdx], defaultRangeTop[currStep+1] - dIdx + 1);
				
                /** probability on last branch determined by prob conservation */
				double lastProb = 1.0;
                
				// loop through number of possible defaults
                int numDefIdx;
				for(numDefIdx = 0; numDefIdx < maxDefBranches; numDefIdx++)
				{
                    /** number of defaults corresponding to numDefIdx */
                    int numDefs = defaults[dIdx + numDefIdx] - defaults[dIdx];

                    /** bottom of outer spread range */
					int outerSpreadBot = 
						(*(*spreadOuterBotLimits)[sIdx - minSpreadRangeBot])[Maths::min(numDefs,maxContagionDefs)];
					
                    /** top of outer spread range */
                    int outerSpreadTop = outerSpreadBot + numSpreadBranches - 1;
					
                    // calc condProb
					if(numDefIdx < maxDefBranches -1)
					{
						condProb = (*condDefTransProbs)[sIdx-spreadOffset][numDefIdx][dIdx];
                        condProb = Maths::collar(condProb,1.0,0.0);
						lastProb -= condProb;
                        
#ifdef DEBUG
                        if(lastProb < 0)
                        {
                            ErrorHandler::writeMsg("lastProb ("+Format::toString(lastProb)+") < 0 at: t:"
                                +Format::toString(currStep) + " dIdx: "+Format::toString(dIdx)+
                                " sIdx:"+Format::toString(sIdx));
                        }
#endif
					}
					else
					{
                        condProb = Maths::collar(lastProb,1.0,0.0);
					}
                    

					
                    // outer spread loop
                    int sIdxOuter;
					for(sIdxOuter = outerSpreadBot; sIdxOuter <= outerSpreadTop; sIdxOuter++) 
					{
                        double spreadProb = 
                            (*spreadTransProbs)[sIdx-minSpreadRangeBot][Maths::min(numDefs,maxContagionDefs)][sIdxOuter - outerSpreadBot];
                
						if(isFwd)
                        {
                            //update the slice
                            sliceGeneral(sIdxOuter, (dIdx + numDefIdx)) += 
                                condProb * spreadProb * (*tempSlice)( sIdx, dIdx) / df;
                        }
                        else
                        {
                            // update the expected value
                            evValue += condProb * spreadProb * (*tempSlice)( sIdxOuter, (dIdx + numDefIdx) );
                        }
                        
							
					} // end outer spread loop	sIdOuter

				} // end outer default loop: numDefIdx

				// update slice value if doing bwd induction
                if(!isFwd)
                {
				    sliceGeneral(sIdx,dIdx) = df*evValue;
                }
				
			} 
		}

	}
	catch (exception& e){
		throw ModelException(e, method);
	}

}



/************************************************************************/
/**   update the condDefTransProbs vector of slices                     */
//	
// Prob P_n[dIdx, sIdx] of n defaults given spread (sIdx) and default level (dIdx)
// is given by expansion of order x
// P_n[dIdx, sIdx] = \sum_{j=0}^{\infty} P^j_n[dIdx, sIdx]
//
// coefficents are defined recursively:
//
// p^0_0[dIdx, sIdx] = 1
//
// P^j_0[dIdx, sIdx] = P^{j-1}_0 ( -f[dIdx]S[sIdx]dt ) / j
//
// P^0_n[dIdx, sIdx] = P^0_{n-1}[dIdx, sIdx] f[dIdx + (n-1)] S[sIdx] dt / n
//
// P^j_n[dIdx, sIdx] = ( f[dIdx+n-1] P^j_{n-1} - f[dIdx+n] P^{j-1}_n ) S[sIdx] dt /(n+j)
//
// here f[] is the locFacArray and S[sIdx] is the spread level at sIdx
/************************************************************************/
void SpreadLossTree::updateCondDefTransProb_expansion(
    int tIdx,
    DoubleArrayArrayArraySP probsToUpdate,
	IntArrayArraySP defBranchesToUpdate) 
{

	static const string method = "SpreadLossTree::updateCondDefTransProb_expansion";
	try
	{
#ifdef DEBUG
		if(tIdx >= getLastStep())
		{
			throw ModelException("tIdx >= getLastStep()");
		}
#endif

#ifdef COUNT_GENERATORS
		GEN_COUNTER++;
#endif
		int dIdx, sIdx; // default and spread counter
		int j,n ;

		/** size of spread offset in condDefProbs */
        // Note that this is dependent on time step
		int spreadOffset = -(probsToUpdate->size()-1)/2;
        int sIdxBot = spreadOffset;
        int sIdxTop = -spreadOffset;

		/** time step at next step */
		int nextStep  = tIdx+1;

		/** year frac between current time step and next time step on timeline */
		double dt = timeLine->TradeYrFrac[nextStep];


		// loop through default states and create array of local factors
		// local factors are those that only depend on default levels and not spread
		DoubleArray locFacArray(maxDefIdxForCache+1);

		// need local factor for each default node at nextStep
		for (dIdx = 0; dIdx <= maxDefIdxForCache; dIdx++)
		{
			locFacArray[dIdx] = localFactor(defaults[dIdx], tIdx);
		}

		// local variables used in loop below;
		
		/** spread at spreadIdx */
		double s;	
		
		/** conditional probability */
		double condProb;

		/** order of expansion of previous probability */
		double prevOrder;


		// first coefficient
		coeffs[0] = 1.0;

		// holder for current coefficient
		double coeff;

       
		
		for(dIdx = 0; dIdx <= maxDefIdxForCache; dIdx++)
		{
            /** maximum number of default probabilities we need to calc at this node */
			// need to cap total number of defaults by numNames
			// note that numProbs >= 1 since always have prob of no default
			int numProbs = (maxDefIdxForCache - dIdx) + 1 ;

			if(numProbs == 1)
			{
				// in this case have single branch so probability along branch must be 1
				for( sIdx = sIdxBot; sIdx <= sIdxTop; sIdx++)
				{
					(*probsToUpdate)[sIdx - spreadOffset][0][dIdx] = 1.0;
					(*(*defBranchesToUpdate)[sIdx-spreadOffset])[dIdx] = numProbs;
				}

				return; // nothing else to do
			}


			// now loop through spreads
            for( sIdx = sIdxBot; sIdx <= sIdxTop; sIdx++)
            {
				// spread at spreadIdx 
				s = spreadProcess->spread(sIdx);		

				prevOrder = -1;
				
                /** last brache prob calculated by conservation of probs */
                double lastProb = 1.0;
				// loop through default branches ----------------
				for(n = 0; ( n < numProbs-1 && lastProb > minTransProb ); n++)
				{
					// intialise j = 0 state
					if(n>0) // note n=0, j=0 state is already initialised to 1
					{
						coeffs[n] = 
							coeffs[ n-1 ] * 
							locFacArray[dIdx + n-1] * s *dt / ((double)n) ;
					}
					condProb = coeffs[n];
				
					// loop over higher orders
					for(j = 1; j <= maxExpansionOrder ; j++)
					{
                        // note that fac1 is 0 if n=0 since prevOrder will be -1
						double fac1 = (j <= prevOrder) ?
							coeffs[j*numProbs + (n-1)] * locFacArray[ dIdx+(n-1) ] /((double)(n+j)) :
							0; // assume that higher order coefficients are 0

						// calculate next order in expansion
                        double fac2 = coeffs[(j-1)*numProbs + n ] * locFacArray[ dIdx+n ] /((double)(n+j));
						
                        coeff = ( fac1 - fac2 ) * s * dt;

						coeffs[j*numProbs + n] = coeff;
						// update probability
						condProb += coeff;
								
						if ((j > prevOrder - n) && ( fabs(coeff) <= minTransProb )){
							break;
						}
		
					}
					
					prevOrder = j-1;

                   
#ifdef DEBUG   /// SOME CONSISTENCY CHECKS ----------------------------------------------

					if(coeff > minTransProb)
					{
						//cout << "coeff not converged: " << coeff 
						//	<< " n: " << n << " sIdx: " << sIdx << " dIdx: " << dIdx 
						//	<< " t: " << tIdx << endl;

						ErrorHandler::writeMsg("coeff (" +
							Format::toString(coeff) +
							") not converged at n:" +
							Format::toString(n) + " sIdx:" +
							Format::toString(sIdx) +" dIdx:" +
                            Format::toString(dIdx)+
                            ", t: "+ Format::toString(tIdx)
							);

					}

					if(condProb < - minTransProb || condProb > 1+minTransProb  )
					{
						//cout << "condProb out of bounds: " << condProb 
						//	<< " n: " << n << " sIdx: " << sIdx << " dIdx: " << dIdx 
						//	<< " t: " << tIdx << endl;
						
						ErrorHandler::writeMsg("conDefTransProb (" +
							Format::toString(condProb) +
							") out of bounds at n:" +
							Format::toString(n) + " sIdx:" +
							Format::toString(sIdx) +" dIdx:" +
                            Format::toString(dIdx)+
                            ", t: "+ Format::toString(tIdx)
							);
						

					}
#endif  // end consistency check -----------------------------

					// cap and floor in case of some (small) numerical violations of prob bounds 
					condProb = Maths::max(0.0, Maths::min( lastProb, condProb));


					// assign probability
					(*probsToUpdate)[sIdx - spreadOffset][n][dIdx] = condProb;	
                    lastProb -=condProb;
                    lastProb = Maths::max(0.0, lastProb); // for small numerical errors
			

				}  // end outer default loop (n)
				// final prob given by conservation if lastProb >0
                if(lastProb>0)
                {
				    (*probsToUpdate)[sIdx-spreadOffset][n][dIdx] = Maths::max(0.0, Maths::min( 1.0, lastProb));
				    (*(*defBranchesToUpdate)[sIdx-spreadOffset])[dIdx] = n+1;
                }
                else
                {
                    (*(*defBranchesToUpdate)[sIdx-spreadOffset])[dIdx] = n;
                }

			}	// end spread loop (sIdx)
		} // end default state loop (dIdx)
	}
	catch (exception& e)
	{
		throw ModelException(e, method);		
	}

}



/************************************************************************/
/**   update the condDefTransProbs vector of slices                     */
//	
// uses Pade approximation
/************************************************************************/
void SpreadLossTree::updateCondDefTransProb_Pade(
    int tIdx,
    DoubleArrayArrayArraySP probsToUpdate,
	IntArrayArraySP defBranchesToUpdate) 
{

	static const string method = "SpreadLossTree::updateCondDefTransProb_Pade";
	try
	{
#ifdef DEBUG
		if(tIdx >= getLastStep())
		{
			throw ModelException("tIdx >= getLastStep()");
		}
#endif
#ifdef COUNT_GENERATORS
		GEN_COUNTER++;
#endif
		int dIdx, sIdx; // default and spread counter
		int n;

		/** size of spreadOffset in condDefProbs */
		int spreadOffset = -(probsToUpdate->size()-1)/2;
        int sIdxBot = spreadOffset;
        int sIdxTop = -spreadOffset;

		/** time step at next step */
		int nextStep  = tIdx+1;

		/** year frac between current time step and next time step on timeline */
		double dt = timeLine->TradeYrFrac[nextStep];


		/** size of generator */
		int genSize = maxDefIdxForCache + 1;

		// allocation of the generator 
        // (nothing will be done if it is already of the right size)
		generator.resizeOnly(genSize, genSize);
        // set all values to 0.0
        generator.fill(0.0);

		// need local factor for each default node at nextStep
		// for dIdx = defaultRange top have absorbing state hence entire row is 0;
        double maxAbsValueOfEigenvalue = 0.0;
		for (dIdx = 0; dIdx < genSize-1; dIdx++)
		{
			const double locFac = localFactor(defaults[dIdx], tIdx);
			// generator[col][row]
			generator[dIdx][dIdx] = -locFac; 
			generator[dIdx+1][dIdx] = locFac;
            if (maxAbsValueOfEigenvalue < fabs(locFac))
            {
                maxAbsValueOfEigenvalue = fabs(locFac);
            }
		}
	    // Compute Spectral norm // YV: this is an approximation ( a ceiling )
        const double _spectralNormOfGenerator = maxAbsValueOfEigenvalue;

		// loop through spreads at tIdx and do Pade exponentiation
		for( sIdx = sIdxBot; sIdx <= sIdxTop; sIdx++)
		{
			
			// spread at spreadIdx 
			const double s = spreadProcess->spread(sIdx);		

			// do pade
			generator.expPade(  dt*s, 
                                _spectralNormOfGenerator, 
                                padeOrder, 
                                &expGen,
                                &aux4exp1,
                                &aux4exp11,
                                &aux4exp2,
                                &aux4exp21,
                                &aux4exp3);
			
	
			for (dIdx = defaultRangeBot[tIdx]; dIdx <= maxDefIdxForCache; dIdx++)
			{

                /** maximum number of default probabilities we need to calc at this node */
                int maxNumProbs = (maxDefIdxForCache - dIdx) +1 ;
				double lastProb = 1.0;
				// loop through inner probs
				for(n = 0; ( n < maxNumProbs-1 && lastProb > minTransProb ); n++)
				{
					double condProb = expGen[dIdx+n][dIdx];
				
					// assign probability
					(*probsToUpdate)[sIdx - spreadOffset][n][dIdx] = condProb;	
					lastProb -= condProb;
                    lastProb = Maths::max(0.0,lastProb);

#ifdef DEBUG   /// SOME CONSISTENCY CHECKS ----------------------------------------------

					if(condProb < - minTransProb || condProb > 1+minTransProb  )
					{
						cout << "condProb out of bounds: " << condProb 
							<< " n: " << n << " sIdx: " << sIdx << " dIdx: " << dIdx 
							<< " t: " << tIdx << endl;

						ErrorHandler::writeMsg("condDefTransProb (" +
							Format::toString(condProb) +
							") out of bounds at n:" +
							Format::toString(n) + " sIdx:" +
							Format::toString(sIdx) +" dIdx:" +
							Format::toString(dIdx)
							);
					}
#endif  // end consistency check -----------------------------

				}  // end outer default loop (n)

				// final prob given by conservation
				(*probsToUpdate)[sIdx-spreadOffset][n][dIdx] = Maths::max(0.0, Maths::min( 1.0, lastProb));
				
				// update number of default branches relevant for calculations
				// note that if using pade the top of the default range is constant at defBound
				// hence defBranchesToUpdate only needs to be updated  when there is a new generator
				// (or in future when we are pricing new instrument with same model)
				// TODO make this contingent on defBound
				(*(*defBranchesToUpdate)[sIdx-spreadOffset])[dIdx] = n+1;
					
			} // end default loop


		}	// end spread loop (sIdx)
		
	}
	catch (exception& e)
	{
		throw ModelException(e, method);		
	}

}

/*****************************************************************************************/
/** factor adjustment to next-to-default intensity that depends only on the default levels*/
// This function needs to be calibrated
/*****************************************************************************************/
double SpreadLossTree::localFactor(
				   int numDefs,
				   int timeIndex
				   ) const
{

	static const string method  = "SpreadLossTree::localFactor";
	
	/** localContagion value */
	double contagion = 1.0;

    double totalLoss = (*totalLosses)[numDefs];

	// if localContagion is provided then get interpolated values
	if(localContagion.get() != NULL)
	{

		DateTime date = getDate(timeIndex);
		
		contagion = localContagion->getSkew( totalLoss, date, SkewSurface::NORMAL );
	}
	
    return contagion * (numNames - numDefs);
}

/************************************************************************/
/*  get the state probabilities for a given time                       */
/************************************************************************/
TreeSliceGeneralSP SpreadLossTree::getStateProbSlice(int step)
{
    static const string method = "SpreadLossTree::getStateProbSlice";
    try
    {
        //check cache
        if(stateProbCache.find(step) != stateProbCache.end())
        {
            return stateProbCache[step];
        }
        else // not in cache so calculate
        {
            /** slice containing state probs */
            TreeSliceGeneralSP  stateProbs = createGeneralSlice();
            // intialise
            *stateProbs = 0;

            // fwd loop through timeline up to date to update state probs
            int t;
            for(t = 0 ; t < step ; t++ )
            {
                // update transition probabilities
                update(t, FDProduct::FWD);

                // update state probs from time t to t+1
                sliceExpand(*stateProbs, discYC->getName());

                // note that this updates the probs from t to t+1
                sliceEv(*stateProbs, FDProduct::FWD);

                
#ifdef DEBUG 
                //check normalisation
                int sIdx,dIdx;
                double probSum =0;
                for(dIdx = defaultRangeBot[t+1] ; dIdx <= defaultRangeTop[t+1]; dIdx++)
                {
                    for(sIdx = spreadRangeBot[t+1]; sIdx <= spreadRangeTop[t+1] ; sIdx++)
                    {
                        probSum += (*stateProbs)(sIdx, dIdx);
                    }
                }
                if(probSum < 1-minTransProb || probSum > 1+minTransProb)
                {
                    ErrorHandler::writeMsg("probSum (" +
                        Format::toString(probSum) +
                        ") out of bounds at sIdx:" +
                        Format::toString(sIdx) +" dIdx:" +
                        Format::toString(dIdx)+
                        ", t: "+ Format::toString(t+1)
                        );
                }
#endif
            }

            // add to cache
            stateProbCache[step] = stateProbs;
            return stateProbs;

        }
        
    }
    catch (exception & e)
    {
       throw ModelException(e, method);
    }
}

/************************************************************************/
/*  get the marginal default prob for a given time                      */
/* P(d) = sum_s P(s ^ d)                                                */
/************************************************************************/
DoubleArraySP SpreadLossTree::getMarginalDefaultProbs(int step)
{
    static const string method = "SpreadLossTree::getMarginalDefaultProbs";
    try
    {
        //check cache
        if(marginalDefProbCache.find(step) != marginalDefProbCache.end())
        {
            return marginalDefProbCache[step];
        }
        else // not in cache so calculate
        {
           TreeSliceGeneral & stateProbs = *getStateProbSlice(step);

           //output intialised to 0
           DoubleArraySP marginals(new DoubleArray(defaults.size(),0));
           //loop through defaults
           int d;
           for(d = 0; d <= defaultRangeTop[step]; d++)
           {
                //loop through spread states at step
               int sIdx;
               for(sIdx = spreadRangeBot[step] ; sIdx <= spreadRangeTop[step]; sIdx++ )
               {
                    (*marginals)[d] += stateProbs(sIdx, d);
               }
           }
            

            // add to cache
            marginalDefProbCache[step] = marginals;
            return marginals;

        }

    }
    catch (exception & e)
    {
        throw ModelException(e, method);
    }
}

/************************************************************************/
/*  get values of slice conditional on loss states                      */
/* for given quantity x(s(t),d(t)) dependent on default level d(t)      */
/* and spread level s(t) and time t we have:                            */
/*                                                                      */
/* P(x | d(t)) = \sum_s(t) P(x | s(t), d(t) ) * P(s(t) ^ d(t)) / P(d(t))*/
/************************************************************************/
DoubleArraySP SpreadLossTree::getValConditionalOnLoss(const TreeSlice & valueSlice,
                                                        int step)
                                      
{
    static const string method = "SpreadLossTree::getValConditionalOnLoss";
    try
    {
        // cast valSlice so that can access values
        const TreeSliceGeneral & valSliceGeneral = dynamic_cast< const TreeSliceGeneral & >(valueSlice);
        //get state probs
        TreeSliceGeneral & stateProbs = *getStateProbSlice(step);

        //get marginals
        DoubleArraySP marginals = getMarginalDefaultProbs(step);

        // output initialised to 0
        DoubleArraySP condVals(new DoubleArray(defaults.size(),0));

       //loop through defaults at step
        int dIdx, sIdx;
        for(dIdx = 0; dIdx <= defaultRangeTop[step] ; dIdx++)
        {
            if( (*marginals)[dIdx] < DEFAULT_MIN_TRANS_PROB)
            {
                (*condVals)[dIdx] = 0;
            }
            else
            {
                //loop through spreads
                for(sIdx = spreadRangeBot[step]; sIdx <=spreadRangeTop[step] ; sIdx++)
                {
                    double stateProb = stateProbs(sIdx, dIdx);
                    (*condVals)[dIdx] += valSliceGeneral(sIdx, dIdx) * stateProb;
                }
                (*condVals)[dIdx] /= (*marginals)[dIdx];

            }
            
            
        }

        return condVals;
    }
    catch (exception & e)
    {
        throw ModelException(e, method);
    }

}

                                      
                                      
//=====================================================================================
//
// insert strikes 
//
//====================================================================================

/** adds strikes to either contLeg...Strike or feeLeg...Strikes if not there already.
returns the index in the array at which the strikes are located */
int SpreadLossTree::addStrikes(
							   double lowStrike,
							   double highStrike,
							   bool recoverNotional
							   )
{
	static const string method ="SpreadLossTree::addStrikes";

	if(lowStrike < 0 || lowStrike > 1)
	{
		throw ModelException("lowStrike out of bounds");
	}
	if(highStrike < 0 || highStrike > 1)
	{
		throw ModelException("highStrike out of bounds");
	}
	if(lowStrike >= highStrike)
	{
		throw ModelException("lowStrike >= highStrike");
	}


	// reference to DoubleArray that we want to add low strikes to
	DoubleArray & lowStrikes = recoverNotional ? elWithRRLowStrikes :elLowStrikes ;

	// reference to DoubleArray that we want to add high strikes to
	DoubleArray & highStrikes = recoverNotional ?  elWithRRHighStrikes : elHighStrikes ;

	// first check if they are included already
	int index = findPair(lowStrike, highStrike, lowStrikes, highStrikes);

	if(index < 0)
	{
		// if arrive here it means that pair of strikes is new
		lowStrikes.push_back(lowStrike);
		highStrikes.push_back(highStrike);

		// return final index
		return lowStrikes.size() - 1;

	}
	else
	{
		// pair is already included
		return index;
	}


}

/** add max strike */
void SpreadLossTree::addMaxStrike(double strike)
{
	if(strike > maxStrike) maxStrike = strike;
}

/** searches for pair of values in pair of arrays
returns index of first pair in the arrays if found otherwise returns -1
*/
int SpreadLossTree::findPair(double val1, double val2, const DoubleArray & array1, const DoubleArray & array2) const
{
	static const string method = "SpreadLossTree::findPair";

	int n = array1.size();
	if(n != array2.size())
	{
		throw ModelException("array1.size() != array2.size()", method);
	}

	for(int i = 0;i< n ; i++)
	{
		if(val1 == array1[i] && val2 == array2[i] )
		{
			return i;
		}
	}

	// if we are here then no matching pair was found
	return -1;

}

//===========================================================================
//
// access to contingent and fee leg values when using effective curve pricing
//
//===========================================================================


/************************************************************************/
/* Get a reference to a specific contingentleg's price at a				*/
/* given time step.														*/
/* to be used with effective curve pricing only							*/
/************************************************************************/
void SpreadLossTree::getCDOContingentLegValue(
	TreeSlice& valueSlice,
	const TrancheContingentLeg& contLeg, 
	int step)
{
	static const string method = "SpreadLossTree::getCDOContingentLegValue";
	try
	{
		// for ease
		string name = contLeg.name;

		// step need to be equal to currStep to make sure that effective curves have been calculated
		if(step != currStep)
		{
			throw ModelException("step != currStep");
		}

		if(!useEffCurves)
		{
			throw ModelException("This method is only for effective curve pricing");
		}

		// if we are here it means that we are using effective curve pricing --------------------

		// variables for ease
		double lowStrike = contLeg.lowStrike;
		double highStrike= contLeg.highStrike;

		/** size of the tranche as fraction of 1 */
		double trancheSize = highStrike - lowStrike;

		// first get index for tranche
		int trancheIdx = findPair(lowStrike, highStrike, elLowStrikes, elHighStrikes);
		
		// check that tranche strikes have been inserted. This shouldn't happen if previous check was ok
		if(trancheIdx < 0)
		{
			throw ModelException(
				"Strikes for TrancheContingentLeg with name "+name + " have not been inserted");
		}
		
		
		/** dates for expected losses */
		DateTimeArray elDates(getLastStep() - step + 1);
		int t;
		for(t = step; t <= getLastStep(); t++)
		{
			elDates[t - step] = getDate(t);
		}

		/** contingent leg effective curve */
		EffectiveCurveSP contLegEffCurve;
		
		
		// effective curve for the relevant tranche need to be extracted from expected loss Curves
		/** effective curve */
		DoubleArray effCurve(0);
		effCurve.resize(getLastStep() - step + 1);


		// get bot and top slice ranges at step = currStep
		/** spread bottom index */
		int sBot = spreadRangeBot[currStep];
		/** spread top index */
		int sTop = spreadRangeTop[currStep];
		/** default bottom index */
		int dBot = defaultRangeBot[currStep];
		/** default top index */
		int dTop = defaultRangeTop[currStep];

		// cast to general slice
		TreeSliceGeneral& contLegSlice = dynamic_cast<TreeSliceGeneral&>(valueSlice);

		// now loop through slice values
		int s,d;
		for(s = sBot; s <= sTop; s++)
		{
			for(d = dBot; d <= dTop; d++)
			{
				// loop through timeline from step to end
				// assume her that entire memory is initialised at the start
				for(t = step ; t <= getLastStep(); t++)
				{
					// value pointer
					//double * pTrancheEL = getValPtr(*trancheEL[trancheIdx][t]);
					TreeSliceGeneral& trancheELslice = *trancheEL[trancheIdx][t];

					effCurve[t - step] = Maths::max(0.0, 1.0 - trancheELslice(s,d)/trancheSize );
				}

				// construct contingent leg effective curve */
				contLegEffCurve = EffectiveCurveSP(new EffectiveCurve(
					elDates[0],							// we want to value to this date
					curveToDiscount.getSP(),
					elDates,
					effCurve,
					EffectiveCurve::FLAT_FORWARD		// interpolation type
					));

				// dummy variables
				CashFlowArray pastTrancheLosses(0);
				BoolArray  payPastTrancheLosses(0);
				DoubleArray debugUnitPrice(0);
				DoubleArray debugUnitHistPrice(0);

				/** value of contingent leg */
				// TODO check the inputs
				double contLegVal = contLeg.contingentLeg->price(
					contLeg.highStrike - contLeg.lowStrike,					// initial notional
					effCurve[0]*trancheSize,		// outstanding notional at start
					elDates[0],												// want to value to this date
					elDates[0],												// valDateCF ??
					contLegEffCurve,										// effective curve
					pastTrancheLosses,										//	pastTrancheLosses
					payPastTrancheLosses,									//  payPastTrancheLosses
					false,													// compute debug prices
					debugUnitPrice,											// debug unit price
					debugUnitHistPrice,										// debug unit hist price
					contLeg.bda     										// Bad day adjuster
					);

				// assign value to relevant contLegSlice
				contLegSlice(s,d) = contLegVal;
				
				
			}
		}

	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
	
}


/************************************************************************/
/* Get a reference to a specific  fee leg's price at a					*/
/* given time step.	Uses effective curve pricing!						*/
/************************************************************************/
void SpreadLossTree::getCreditFeeLegValue(
	TreeSlice& valueSlice,	// (O) slice to fill
	const TrancheFeeLeg& feeLeg, 
	int step,
	IForwardRatePricerSP model) //required for fee leg
{
	static const string method = "SpreadLossTree::getCreditFeeLegValue";
	try
	{
		// for ease
		string name = feeLeg.name;

		// step need to be equal to currStep to make sure that effective curves have been calculated
		if(step != currStep)
		{
			throw ModelException("step != currStep");
		}

		
		if(!useEffCurves)
		{
			throw ModelException("This method is only used with effective curve pricing");
		}

		// If we are here it means we are pricing using effective curves ---------------------

		// variables for ease
		double lowStrike = feeLeg.lowStrike;
		double highStrike= feeLeg.highStrike;

		/** size of the tranche as fraction of 1 */
		double trancheSize = highStrike - lowStrike;

		// first get index for tranche
		int trancheIdx = findPair(
			lowStrike, highStrike, 
			feeLeg.recoverNotional ? elWithRRLowStrikes : elLowStrikes, 
			feeLeg.recoverNotional ? elWithRRHighStrikes : elHighStrikes
			);
		
		// check that feeLeg has been inserted
		if (trancheIdx < 0) 
		{
			throw ModelException(
				"Strikes for TrancheFeeLeg with name "+ feeLeg.name + " have not been inserted");
		}


		/** dates for expected losses */
		DateTimeArray elDates(0);
		elDates.resize(getLastStep() - step + 1);
		int t;
		for(t = step; t <= getLastStep(); t++)
		{
			elDates[t - step] = getDate(t);
		}

		/** contingent leg effective curve */
		EffectiveCurveSP feeLegEffCurve;

		
		// effective curve for the relevant tranche need to be extracted from expected loss Curves
		/** effective curve */
		DoubleArray effCurve(0);
		effCurve.resize(getLastStep() - step + 1);

		

		// get bot and top slice ranges at step = currStep
		/** spread bottom index */
		int sBot = spreadRangeBot[currStep];
		/** spread top index */
		int sTop = spreadRangeTop[currStep];
		/** default bottom index */
		int dBot = defaultRangeBot[currStep];
		/** default top index */
		int dTop = defaultRangeTop[currStep];


		TreeSliceGeneral& feeLegSlice = dynamic_cast<TreeSliceGeneral&>(valueSlice);
		// TODO check pointer

		// now loop through slice values
		int s,d;
		for(s = sBot; s <= sTop; s++)
		{
			for(d = dBot; d <= dTop; d++)
			{
				// loop through timeline from step to end
				// assume her that entire memory is initialised at the start
				for(t = step ; t <= getLastStep(); t++)
				{
					TreeSliceGeneral& trancheELslice = feeLeg.recoverNotional ?
						*trancheELwithRR[trancheIdx][t] :
						*trancheEL[trancheIdx][t];


					effCurve[t - step] = Maths::max(0.0, 1.0 - trancheELslice(s,d)/trancheSize );
				}
					

				// construct contingent leg effective curve */
				feeLegEffCurve = EffectiveCurveSP(new EffectiveCurve(
					elDates[0],							// we want to value to this date
					curveToDiscount.getSP(),
					elDates,
					effCurve,
					EffectiveCurve::FLAT_FORWARD		// interpolation type
					));

				
				// output variables NOT USED
				double riskyDurationTotal;								// (O)
				double riskyNotionalsMean;								// (O)
				double	risklessCFPV;

				// dummy variables
				CashFlowArray pastTrancheLosses(0);
				DoubleArray debugUnitPrice(0);
				DoubleArray debugUnitHistPrice(0);
				CashFlowArraySP rebatePayments;

				BoolArraySP     payAsYouGoArray;
				IntArraySP      numDelayDaysArray;
				DateTimeArraySP startDates;
				DateTimeArraySP endDates;
				DateTimeArraySP paymentDates;

				/** value of fee leg */
				double feeLegVal = feeLeg.feeLeg->price(
					elDates[0],												// today
					elDates[0],												// valDateCF
					feeLegEffCurve,											// effective curve
					curveToDiscount,										// YC wrapper
					feeLeg.lowStrike,
					feeLeg.highStrike,
					effCurve[0]*trancheSize,								// outstanding notional
					pastTrancheLosses,										// pastTrancheLosses
					riskyDurationTotal,										// (O)
					riskyNotionalsMean,										// (O)
					risklessCFPV,											// (O)
					false,													// compute extra?
					debugUnitPrice,											// debugUNitPrice
					debugUnitHistPrice,										// debugUnitHistPrice
					rebatePayments,											// rebate payments
					payAsYouGoArray,
					numDelayDaysArray,
					startDates,
					endDates,
					paymentDates,
					feeLeg.bda,  											//IBadDayAdjuster
					model); 
				

				// assign value to relevant feeLegSlice
				feeLegSlice(s,d) = feeLegVal;
				
				
			}
		}
		
	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

//==========================================================================================
//
// Primitive Instruments
//
//==========================================================================================

/************************************************************************/
/* TrancheContigentLeg                                              */
/************************************************************************/
/**public constructor */
TrancheContingentLeg::TrancheContingentLeg(
	ICreditContingentLegSP contingentLeg,		// underlying contingent leg
	double lowStrike,							// %
	double highStrike,							// %
	const IBadDayAdjusterConstSP  bda
	) :  
CObject(TYPE), 
contingentLeg(contingentLeg),
lowStrike(lowStrike),
highStrike(highStrike),
bda(bda)
{
	static int uniqueId=0;
	name = "TrancheContingentLeg-"+Format::toString(++uniqueId);

}

/** returns critical dates */
DateTimeArraySP TrancheContingentLeg::getCritDates() const {
	// TODO get protection start and end form ICDOContingentLeg
	DateTimeArraySP critDates(0);
	return critDates;
}

FDProductSP TrancheContingentLeg::createProduct( FDModel * model ) const {
	return FDProductSP(new TrancheContingentLegProd(TrancheContingentLegConstSP(this), model));
}

void TrancheContingentLeg::load(CClassSP& clazz)
{
	REGISTER(TrancheContingentLeg, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(FDModel::IIntoProduct);
	EMPTY_SHELL_METHOD(defaultConstructor);

	FIELD(contingentLeg, "Contingent leg");

	FIELD(lowStrike,"");
	FIELD(highStrike,"");

	FIELD(bda, "bad day adjuster");

	FIELD(name,""); 
	FIELD_MAKE_TRANSIENT(name);
	
}

CClassConstSP const TrancheContingentLeg::TYPE = CClass::registerClassLoadMethod(
	"TrancheContingentLeg", typeid(TrancheContingentLeg), load);

/** private constructor */
TrancheContingentLeg::TrancheContingentLeg() :
CObject(TYPE)
{}
/************************************************************************/
/* TrancheContingentLegProd                                             */
/************************************************************************/
TrancheContingentLegProd::TrancheContingentLegProd(const TrancheContingentLegConstSP &inst, FDModel* model)
: FDProduct(model), inst(inst)
{
	tree = dynamic_cast<SpreadLossTree*>(model);
	if (!tree) throw ModelException("TrancheContingentLegProd", "Model must derive from SpreadLossTree");

    // get start date
    startDate = inst->contingentLeg->firstObservationStartDate();
}

void TrancheContingentLegProd::init(Control*) const {
	try {

		// add critical dates to model timeline
		tree->addCritDate(inst->contingentLeg->firstObservationStartDate());
		tree->addCritDate(inst->contingentLeg->lastObservationEndDate());

		// add max strike to limit tree size
		tree->addMaxStrike(inst->highStrike);

		
	}
	catch (exception& e) {
		throw ModelException(e, "TrancheContingentLegProd::init, Error for instrument type " +
			inst->getClass()->getName() + " and name " +
			inst->getName());
	}
}

void TrancheContingentLegProd::initProd(){
	// create a single slice to store the tree contingent leg values
	mainSlice = tree->createSlice();
	mainSlice->name = "TrancheContingentLeg";
	// initialise
	*mainSlice = 0;

	if(tree->useEffCurves)
	{
		tree->addStrikes(inst->lowStrike, inst->highStrike, false);
	}
	else
	{
		// add slices to DEV
		startDEV(mainSlice); 
	}

}

void  TrancheContingentLegProd::update(int & step, FDProduct::UpdateType)
{

	if(!tree->useEffCurves)
	{
        // update between protection dates
        const DateTime & protStart = inst->contingentLeg->firstObservationStartDate();
        const DateTime & protEnd = inst->contingentLeg->lastObservationEndDate();
        if(tree->getDate(step).isGreaterOrEqual(protStart) && 
            tree->getDate(step).isLess(protEnd))
        {
		    tree->updateContLeg(*mainSlice, inst->lowStrike, inst->highStrike);
        }
        
	}

}

const TreeSlice & TrancheContingentLegProd::getValue(int step) const
{
	if(tree->useEffCurves)
	{
		tree->getCDOContingentLegValue(
			*mainSlice,
			*inst, 
			step);
	}
	
	return *mainSlice;
	
}



/************************************************************************/
/* TrancheFeeLeg	                                                    */
/************************************************************************/
/**public constructor */
TrancheFeeLeg::TrancheFeeLeg(
	ICreditFeeLegSP feeLeg,				// underlying fee leg
	double lowStrike,					// %
	double highStrike,					// %	
	bool recoverNotional,
	const IBadDayAdjusterConstSP bda
	) :  
CObject(TYPE), 
feeLeg(feeLeg),
lowStrike(lowStrike),
highStrike(highStrike),
recoverNotional(recoverNotional),
bda(bda)
{
	static int uniqueId=0;
	name = "TrancheFeeLeg-"+Format::toString(++uniqueId);

}

/** returns critical dates */
DateTimeArraySP TrancheFeeLeg::getCritDates() const {
	// TODO get cpn dates from ICreditFeeLeg
	DateTimeArraySP critDates(0);
	return critDates;
}

FDProductSP TrancheFeeLeg::createProduct( FDModel * model ) const {
	return FDProductSP(new TrancheFeeLegProd(TrancheFeeLegConstSP(this), model));
}

void TrancheFeeLeg::load(CClassSP& clazz)
{
	REGISTER(TrancheFeeLeg, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(FDModel::IIntoProduct);
	EMPTY_SHELL_METHOD(defaultConstructor);

	FIELD(feeLeg, "Fee leg");

	FIELD(lowStrike,"");
	FIELD(highStrike,"");
	FIELD(recoverNotional, "Is the notional written down by the recoverd amount?");

	FIELD(bda, "bad day adjuster");

	FIELD(name,"");
	FIELD_MAKE_TRANSIENT(name);
}

CClassConstSP const TrancheFeeLeg::TYPE = CClass::registerClassLoadMethod(
	"TrancheFeeLeg", typeid(TrancheFeeLeg), load);

/** private constructor */
TrancheFeeLeg::TrancheFeeLeg() :
CObject(TYPE)
{}

/************************************************************************/
/* TrancheFeeLegProd                                             */
/************************************************************************/
TrancheFeeLegProd::TrancheFeeLegProd(const TrancheFeeLegConstSP &inst, FDModel* model)
: FDProduct(model), inst(inst)
{
	tree = dynamic_cast<SpreadLossTree*>(model);
	if (!tree) throw ModelException("TrancheFeeLegProd", "Model must derive from SpreadLossTree");

    // start date is first accrual date 
    startDate = inst->feeLeg->getAccrualPeriods()->front()->startDate();
}

void TrancheFeeLegProd::init(Control*) const 
{
	try {

		tree->addMaxStrike(inst->highStrike);		

		// add flow ates to model timeline
		tree->addCritDates(*inst->feeLeg->getCashFlowDates());
	
	}
	catch (exception& e) {
		throw ModelException(e, "TrancheFeeLegProd::init, Error for instrument type " +
			inst->getClass()->getName() + " and name " +
			inst->getName());
	}
}

void TrancheFeeLegProd::initProd()
{
	// create a single slice to store the tree fee leg values
	mainSlice = tree->createSlice();
	mainSlice->name = "TrancheFeeLeg";
	// intialise
	*mainSlice = 0;


	if(tree->useEffCurves)
	{
		tree->addStrikes(inst->lowStrike, inst->highStrike,
			inst->recoverNotional		
			);	
	}
	else
	{
		// add slices to DEV
		startDEV(mainSlice); 

		ClosedFormForwardRatePricerSP fwdModel =
			ClosedFormForwardRatePricerSP(
			new ClosedFormForwardRatePricer());

		// Risky -----------------
		CashFlowArraySP flows = inst->feeLeg->getRiskyCashFlows(fwdModel);
		DateTimeArray riskyDates = CashFlow::dates(*flows);

		riskyFlowAmounts.resize(riskyDates.size());
		int i;
		for(i = 0 ; i < riskyDates.size(); i++)
		{
			riskyFlowAmounts[i] = (*flows)[i].amount;
		}


		// riskless --------------
		flows = inst->feeLeg->getRisklessCashFlows(fwdModel);
		DateTimeArray risklessDates = CashFlow::dates(*flows);

		risklessFlowAmounts.resize(risklessDates.size());
		for(i = 0 ; i < risklessDates.size(); i++)
		{
			risklessFlowAmounts[i] = (*flows)[i].amount;
		}

		// model timeline for ease
		DateTimeArray  modelDates = tree->getDates();

		// define projections of cashflow dates onto model timeline
		riskyDateIdx = DateTime::getProjection(modelDates, riskyDates);
		risklessDateIdx = DateTime::getProjection(modelDates, risklessDates);



		// set up slice with outstanding fraction ----------------------------------------
		// need to calculate only once since default slice values stay the same

		outstandingFraction = tree->createSlice();
		outstandingFraction->name = "outstandingFraction";

		// variables for ease
		double lowStrike = inst->lowStrike;
		double highStrike = inst->highStrike;
		double trancheSize = highStrike - lowStrike;


		TreeSliceSP losses = tree->getLossSlice();

		*outstandingFraction = 1.0
				- ( smin( *losses, highStrike ) - smin( *losses, lowStrike ) ) / trancheSize;

		if(inst->recoverNotional)
		{
			TreeSliceSP rr = tree->getRecoverySlice();
			// if reducing by recovered notional as well need to add term
			*outstandingFraction -= 
				( smin( *rr, 1 - lowStrike) - smin( *rr, 1 - highStrike) ) / trancheSize;
		}
	}


}

void TrancheFeeLegProd::update(int & step , FDProduct::UpdateType)
{
	static const string method ="TrancheFeeLegProd::update";
	try
	{

		if(!tree->useEffCurves)
		{	
			//  add any cashflows that happen at step ------------------

			//check if currStep corresponds to feeLeg date
			/** index in riskyFlows for current date */ 
			int riskyIdx = riskyDateIdx[step];
			/** index in risklessFlows for current date */ 
			int risklessIdx = risklessDateIdx[step];

			if(riskyIdx == -1 && risklessIdx == -1) return; // not on flow date so nothing else to do


			// variables for ease
			double riskyAmount = riskyIdx >= 0 ? 
				riskyFlowAmounts[riskyIdx] :
				0;
			double risklessAmount = risklessIdx >= 0 ?
				risklessFlowAmounts[risklessIdx] :
				0;

			*mainSlice += (*outstandingFraction)*riskyAmount + risklessAmount;

			//TODO accrual

		}
	

	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
}

const TreeSlice & TrancheFeeLegProd::getValue(int step) const
{
	if(tree->useEffCurves)
	{
		//don't want to muck around with FD models
		//but this limits the type of cashflows supportable
		//by the tree
		ClosedFormForwardRatePricerSP fwdModel =
			ClosedFormForwardRatePricerSP(
			new ClosedFormForwardRatePricer());


		tree->getCreditFeeLegValue(
			*mainSlice,
			*inst, 
			step,
			fwdModel);

	}
	
	return *mainSlice;
	
	
}



DRLIB_END_NAMESPACE
