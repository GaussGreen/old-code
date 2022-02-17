//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : OUAuxSpreadProcess.hpp
//
//   Description : OU auxiliary spread process. Instance of IAuxSpreadProcess
//
//   Author      : Matthias Arnsdorf
//
//   Date        : September 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/OUAuxSpreadProcess.hpp"

DRLIB_BEGIN_NAMESPACE

/** default value for maximum number of defaults */
const int OUAuxSpreadProcess::DEFAULT_MAX_DEFAULTS = 1;

/** default value for number of branches for tree (assumed fixed) */
const int OUAuxSpreadProcess::DEFAULT_NUM_BRANCHES = 3;

/** default val for minSpreadVal */
const int OUAuxSpreadProcess::DEFAULT_MIN_SPREAD_VAL = 0;

/** default val for q */
const double OUAuxSpreadProcess::DEFAULT_Q = 0.0;

/** scale factor for dR: dR^2 = SCALE * avegrageDT * vol^2  */
// 3 is value given by HULL
const double OUAuxSpreadProcess::SCALE = 3.0;

/** if |prevDT - dt| < DT_ACCURACY they will be assumed same */
const double OUAuxSpreadProcess::DT_ACCURACY = 1e-10;


/** private constructor */
OUAuxSpreadProcess::OUAuxSpreadProcess() : 
CObject(TYPE), 
fwdDates(0),
cleanSpreads(),
maxDefaults(DEFAULT_MAX_DEFAULTS),
timeline(),
fwds(0),
dR(0),
contagionJump(0),
q(DEFAULT_Q),
avgDT(0),
topLimits(0),
botLimits(0),
prevDT(-1)
{
}

/** Destructor */
OUAuxSpreadProcess::~OUAuxSpreadProcess()
{}


void OUAuxSpreadProcess::load (CClassSP& clazz) {
	clazz->setPublic(); // make visible to EAS/spreadsheet
	REGISTER(OUAuxSpreadProcess, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(IAuxSpreadProcess);
	EMPTY_SHELL_METHOD(defaultOUAuxSpreadProcess);

	FIELD(indexCurve, "Mean reversion level as term structure of index spreads");

	FIELD(vol, "log volatility (Black vol)");

	FIELD(meanReversion, "speed of mean reversion");

    FIELD(contagionJump, "Jump as % of spread level when have 1 or more defaults. [default = 0]");
    FIELD_MAKE_OPTIONAL(contagionJump);

    FIELD(maxDefaults, "Above this number all defaults will cause the same contagion spread jump. [default = 1]");
    FIELD_MAKE_OPTIONAL(maxDefaults);

    FIELD(q, "q-mapping parameter. Log-normal: q = 0, Normal: q = 0, Sub-nomral q = -1. [default = "+
		Format::toString(DEFAULT_Q)+"]");
	FIELD_MAKE_OPTIONAL(q);

	FIELD(minSpreadVal,"mimumum value for spreads. [default = " + Format::toString(DEFAULT_MIN_SPREAD_VAL) + "]");
	FIELD_MAKE_OPTIONAL(minSpreadVal);

	FIELD(fwdDates,"");
	FIELD_MAKE_TRANSIENT(fwdDates);

	FIELD(cleanSpreads,"");
	FIELD_MAKE_TRANSIENT(cleanSpreads);

	FIELD(qMap,"");
	FIELD_MAKE_TRANSIENT(qMap);

	FIELD(transitionProbs,"");
	FIELD_MAKE_TRANSIENT(transitionProbs);

	FIELD(avgDT,"");
	FIELD_MAKE_TRANSIENT(avgDT);

	FIELD(dR,"");
	FIELD_MAKE_TRANSIENT(dR);

    FIELD(scaledJump,"");
    FIELD_MAKE_OPTIONAL(scaledJump);


}

IObject* OUAuxSpreadProcess::defaultOUAuxSpreadProcess() {
	return new OUAuxSpreadProcess();
}

CClassConstSP const OUAuxSpreadProcess::TYPE = 
CClass::registerClassLoadMethod("OUAuxSpreadProcess", 
								typeid(OUAuxSpreadProcess), 
								load);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
* script to force the linker to include this file */
bool OUAuxSpreadProcessLoad() 
{
	return (OUAuxSpreadProcess::TYPE != 0);
}



/************************************************************************/
/*   validate                                                           */
/************************************************************************/
void OUAuxSpreadProcess::validatePop2Object()
{
	static const string method = "OUAuxSpreadProcess::validatePop2Object";
	try
	{
		// check that numBranches is odd
		if( (int)((((double)DEFAULT_NUM_BRANCHES) + 1.0) / 2.0) !=  (DEFAULT_NUM_BRANCHES + 1) / 2 )
		{
			throw ModelException("Number of tree branches needs to be odd");
		}

		// chcek vol
		if(vol < 0)
		{
			throw ModelException("vol (" + Format::toString(vol) +") needs to be positive");
		}
		// check mean reversion
		if(meanReversion < 0)
		{
			throw ModelException("meanReversion (" + Format::toString(meanReversion) +
				") needs to be positive");
		}

       

		// check minSpreadVal
		if(minSpreadVal < 0)
		{
			throw ModelException("minSpreadVal (" + Format::toString(minSpreadVal) +
				") needs to be positive");
		}

		

        // maxDefaults
        if(maxDefaults< 0)
        {
            throw ModelException("maxDefaults ("+Format::toString(maxDefaults)+") needs to be >= 0");
        }
        // set max defaults to 0 if have no contagion
        if(contagionJump == 0)
        {
            maxDefaults = 0;
        }
        
        
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}

/************************************************************************/
/* chance to pass market data to the spread process                     */
/************************************************************************/
void OUAuxSpreadProcess::getMarketData(const MarketData* market,
									   int stepsPerYear)
{
	// Dummy model 
	CModelSP  dummy(new NonPricingModel());
	indexCurve.getData(dummy.get(), market);


	// intialise fwds -------------------------------------------------------------
	DefaultRatesSP defRates = indexCurve.getSP()->defaultRates();
	CashFlowArraySP cleanSpreadFlows = defRates->getCleanSpreadCurve();

	fwdDates = CashFlow::dates(*cleanSpreadFlows);
	cleanSpreads = CashFlow::amounts(*cleanSpreadFlows);

	// calc tree spacing
	avgDT = 1.0 / (double)stepsPerYear;

	dR = vol * sqrt( SCALE * avgDT );

    scaledJump = dR > 0 ? contagionJump/dR : 0 ;
	
}

/************************************************************************/
/*  initialise the model. Called prior to any requests                  */
/************************************************************************/
void OUAuxSpreadProcess::setupModel(
									TimeLineSP timeLine		// timeline to use
									)
{
	try
	{
		timeline = timeLine;
		
		// check that getMarketData has been called
		if (!cleanSpreads)
		{
			throw ModelException("indexCurve has not been initialised. Has getMarketData() been called?");
		}
		

		// assume that the clean spread are FLAT_FORWARD (this is correct currently: 19.9.06)

		/** fwdIdx[t] gives index of first date int spreadDates >= stepDates[t] */ 
		vector<int> fwdIdx = DateTime::getCeilingProjection(timeLine->StepDates, fwdDates);

		fwds.resize(fwdIdx.size());
		int t = 0;
		for(; t < (int)fwdIdx.size(); t++)
		{
			fwds[t] = (*cleanSpreads)[fwdIdx[t]];
		}


		// now calculate tree limits
		calcLimits();


		// initialise range -----------------------------------------------------------
		
		
		// range is centred around 0
		/** top of tree slice range for final slice */
		int rangeTop = topLimits.back();
		/** bottom of tree slice range for final slice */
		int rangeBot = botLimits.back();

		minRangeBot = rangeBot;

		/** number of tree nodes at final tree slice */
		int finalNumNodes = rangeTop - rangeBot + 1;

		range = TreeSliceGeneral::Range::create(
			1, // count
			rangeBot,	rangeTop
			);
		
		//set up  spread values 
		spreads.resize(finalNumNodes);

		// set up qMap ----------------------------------------
		/** Ito adjustment so that mean revert to fwd level */
		//double adj = (1-q) * vol * vol / (2 * meanReversion);
        double adj = 0; //TEST
		
		qMap.resize(finalNumNodes);
		int sIdx;
		for(sIdx = rangeBot; sIdx <= rangeTop; sIdx++)
		{
			if(q==1.0)
			{
                qMap[sOffset(sIdx)] = 1.0 + sIdx*dR;		
			}
			else
			{
				qMap[sOffset(sIdx)] = 1.0 + ( exp( (1-q) * (sIdx*dR  - adj ) ) - 1.0 ) / (1-q) ;
			}
			
		}

		// initialise shifts to 0
		shifts.resize(finalNumNodes);
		int i,j;
		for(i = 0; i < finalNumNodes; i++)
		{
			shifts[i].resize(maxDefaults + 1);
			for(j = 0; j <= maxDefaults; j++)
			{
				shifts[i][j] = 0;
			}
		}

		outerBot = IntArrayArraySP(new IntArrayArray(finalNumNodes));
		
		for(sIdx = rangeBot; sIdx <= rangeTop; sIdx++)
		{
			(*outerBot)[sOffset(sIdx)] = IntArraySP( new IntArray(maxDefaults + 1));
			for(j = 0; j <= maxDefaults; j++)
			{
				(*(*outerBot)[sOffset(sIdx)])[j] = sIdx + shifts[sOffset(sIdx)][j] - (DEFAULT_NUM_BRANCHES - 1)/2;
			}
		}


		// allocate memory for the transition probs -----------------------------------
		transitionProbs = DoubleArrayArrayArraySP(new DoubleArrayArrayArray(finalNumNodes));
		int bIdx;
		for(i = 0; i < finalNumNodes; i++)
		{
			(*transitionProbs)[i].resize(maxDefaults + 1);
			for(bIdx = 0; bIdx < maxDefaults + 1; bIdx++)
			{
				(*transitionProbs)[i][bIdx].resize(DEFAULT_NUM_BRANCHES);
			}
		}

		// initialise prevDT
		prevDT = -1;
	}
	catch (exception & e)
	{
		throw ModelException(e, "OUAuxSpreadProcess::setupModel");
	}
}


/************************************************************************/
/*  fwd loop to calc tree limits                                       */
/*  defines topLimits array												*/
/* botLimits are assumed to be equal and opposite						*/
/************************************************************************/
void OUAuxSpreadProcess::calcLimits()
{
	try
	{
		int numSteps = timeline->NumOfStep;

		topLimits.resize(numSteps+1);
		botLimits.resize(numSteps+1);
		scaledVar.resize(numSteps);

		topLimits[0] = 0;
		botLimits[0] = 0;
		
		double dt;

		/** first value of j at which we can switch geometry */
		int jMin;

		int t = 1;
		for(; t <= numSteps; t++)
		{
			dt = timeline->TradeTime[t] - timeline->TradeTime[t-1];
			scaledVar[t-1] = dt / (avgDT * SCALE); // = vol^2dt/dR^2
			
			if(meanReversion > 0)
			{
				jMin = 1 + (int) ( (1.0 - sqrt(1.0 - scaledVar[t-1])) / (meanReversion * dt) );

				topLimits[t] = topLimits[t-1] >= jMin ?
					topLimits[t-1] :
				topLimits[t-1] + 1;

				botLimits[t] = -topLimits[t];

			}
			else
			{
				topLimits[t] = topLimits[t-1] + 1;

				botLimits[t] = -topLimits[t];
			}
			

		}


	}
	catch(exception & e)
	{
		throw ModelException(e, "OUAuxSpreadProcess::calcLimits");
	}

}

/************************************************************************/
/*  updates tree variables at each time step                            */
/************************************************************************/
void OUAuxSpreadProcess::update(int timeStep, FDProduct::UpdateType type)
{
	try
	{
#ifdef DEBUG
		// check that setupModel has been called
		if(!timeline)
		{
			throw ModelException("setupModel has not been called");
		}
		if(timeStep > timeline->NumOfStep)
		{
			throw ModelException("illegal timeStep");
		}
#endif

		// update the current spread values
		updateSpread(timeStep, type);

		if(timeStep == timeline->NumOfStep)
			return; // nothing else to do

		// if we are here we are before last time step --------------------------------
		
		/** time step */
        double dt = timeline->TradeYrFrac[timeStep+1];
   
        // value at t of scaled variance
		double scaledVarVal = scaledVar[timeStep];

		/** mean reversion* dt */
		double mrDt = meanReversion * dt;
       

        // calculate the center nodes
        calcShifts(timeStep, dt);


		// calculate transition probs
		
		/** number of defaults */
		int nd;
		for(nd = 0; nd <= maxDefaults; nd++)
		{
			// first update the outer edges
			int sIdx = topLimits[timeStep];

			/** mean of R over time step divided by dR */ 
			double scaledMean = - mrDt * sIdx + scaledJump*nd;

			calcTransProbs(
				scaledMean, scaledVarVal, shifts[sOffset(sIdx)][nd], (*transitionProbs)[sOffset(sIdx)][nd]
				);

			sIdx = botLimits[timeStep];

			scaledMean = - mrDt * sIdx + scaledJump*nd;
			calcTransProbs(
				scaledMean, scaledVarVal, shifts[sOffset(sIdx)][nd], (*transitionProbs)[sOffset(sIdx)][nd]
				);

			// only need to update the remaining probs if dt has changed or if nd > 0
			if(!Maths::areEqualWithinTol(dt, prevDT, DT_ACCURACY) || nd > 0)
			{
				for(sIdx = botLimits[timeStep]+1; sIdx <= topLimits[timeStep]-1 ; sIdx++ )
				{
					scaledMean = -mrDt * sIdx + scaledJump*nd;
					calcTransProbs(
						scaledMean, scaledVarVal, shifts[sOffset(sIdx)][nd], (*transitionProbs)[sOffset(sIdx)][nd]
						);
				}

			}
			
		}
        prevDT = dt;

	}
	catch(exception & e)
	{
		throw ModelException(e, "OUAuxSpreadProcess::update");
	}
}

/************************************************************************/
/*   calculate centreNode array  and bot limits                         */
/************************************************************************/
void OUAuxSpreadProcess::calcShifts(int timeStep, double dt)
{
	static const string method = "OUAuxSpreadProcess::calcShifts";


#ifdef DEBUG
	if(timeStep >= timeline->NumOfStep)
	{
		throw ModelException("illegal timeStep", "OUAuxSpreadProcess::calcCentreNodes");
	}
#endif

	

	int bot = botLimits[timeStep];
	int top = topLimits[timeStep];
    int nextBot = botLimits[timeStep+1];
    int nextTop = topLimits[timeStep+1];

    double mrDt = meanReversion*dt;

    /** offset between center and edge nodes */
    int brancheOffset = (DEFAULT_NUM_BRANCHES - 1)/2;


	// note that shifts are initialised to 0
	
    int numDefs;

	// do tree edges first. These are same regardless of the number of defaults
    // since we have absorbing state at the edges
    // This assumes brancheOffset = 1. To change this need to update more than the edges
    QLIB_VERIFY(brancheOffset==1, "brancheOffset needs to be 1 here.");

    for(numDefs = 0; numDefs <= maxDefaults; numDefs++)
    {
	    if(top < topLimits[timeStep +1]) // standard branching
	    {
		    shifts[sOffset(top)][numDefs] = 0;
		    shifts[sOffset(bot)][numDefs] = 0;

		    (*(*outerBot)[sOffset(top)])[numDefs] = top - brancheOffset;
		    (*(*outerBot)[sOffset(bot)])[numDefs] = bot - brancheOffset;
	    }
	    else 
	    {
		    shifts[sOffset(top)][numDefs] = - 1; // branch down
		    (*(*outerBot)[sOffset(top)])[numDefs] = top - 1 - brancheOffset;
    		
		    shifts[sOffset(bot)][numDefs] = + 1; // branch up
		    (*(*outerBot)[sOffset(bot)])[numDefs] = bot + 1 - brancheOffset;
	    }
    }

    // loop through inside of slice
    // for numDefs = 0 shift is always 0 here and doesn't need to be updated
    // unless we are doing fwd induction in which case need to set the next to
    // edge node shifts to 0
    
    if(timeStep > 0) // note in this case we have at least 3 branches
    {
        shifts[sOffset(top-1)][0] = 0;
        shifts[sOffset(bot+1)][0] = 0;

        (*(*outerBot)[sOffset(top-1)])[0] = (top-1) - brancheOffset;
        (*(*outerBot)[sOffset(bot+1)])[0] = (bot+1) - brancheOffset;
        
    }

    //TEST
    if(timeStep == timeline->NumOfStep-1)
    {
        for(int sIdx = bot+1; sIdx <= top-1; sIdx++)
        {
            shifts[sOffset(sIdx)][0] = 0;
            (*(*outerBot)[sOffset(sIdx)])[0] = sIdx - brancheOffset;
        }
    }

    
    // now consider numDefs >0
    // note that we need nextBot+brancheOffset <= sIdx + shift <= nextTop-brancheOffset 
    // so that branches are contained within tree limits
    // only need to update if dt has changed. Note that prev dt = -1 at back of tree
   
    for(int sIdx = bot+1; sIdx <= top-1; sIdx++)
    {
        for(numDefs = 1; numDefs <= maxDefaults; numDefs++)
        {
            double scaledMean = - mrDt * sIdx + scaledJump*numDefs;
            int shift = Maths::round(scaledMean);

            shifts[sOffset(sIdx)][numDefs] = 
                Maths::collar(shift, nextTop-brancheOffset-sIdx, nextBot+brancheOffset-sIdx);

            (*(*outerBot)[sOffset(sIdx)])[numDefs] = 
                sIdx + shifts[sOffset(sIdx)][numDefs] - brancheOffset;

        }
    }
   

}

/************************************************************************/
/*  calculate transition probabilites                                   */
/************************************************************************/
void OUAuxSpreadProcess::calcTransProbs(
					double scaledMean,			// mean over period divided by dR
					double scaledVarVal,			// dt/avgDT/SCALE
					int shift,					// shift at current node
					DoubleArray & transProbs	// array of transition probs for each branch	
					)
{
	static const string method = "OUAuxSpreadProcess::calcTransProbs";

#ifdef DEBUG
	if(transProbs.size() != 3)
	{
		throw ModelException("transProbs has wrong size", method);
	}
#endif

    // adjust mean by shift
    double adjMean = scaledMean - shift;
	double adjMeanSq = adjMean*adjMean;
	

    // calc probs
    // need to cap and floor in case have jump outside of tree bounds
	transProbs[0] = (scaledVarVal - adjMean + adjMeanSq) / 2;	    // bottom branch
    transProbs[0] = Maths::collar(transProbs[0], 1.0, 0.0);

    transProbs[2] = (scaledVarVal + adjMean + adjMeanSq) / 2;		    // top branch
    transProbs[2] = Maths::collar(transProbs[2], 1-transProbs[0], 0.0);

	transProbs[1] = 1 - transProbs[0] - transProbs[2];					// middle branch
	transProbs[1] = Maths::collar(transProbs[1], 1-transProbs[0], 0.0); // just in case there is some weird numerical stuff happening
		


#ifdef DEBUG
	if(transProbs[0] < 0 || transProbs[0] > 1)
	{
		throw ModelException("transProbs[0] = "+
			Format::toString(transProbs[0]) +" out of bounds");
	}
	if(transProbs[1] < 0 || transProbs[1] > 1)
	{
		throw ModelException("transProbs[1] = "+
			Format::toString(transProbs[1]) +" out of bounds");
	}
	if(transProbs[2] < 0 || transProbs[2] > 1)
	{
		throw ModelException("transProbs[2] = "+
			Format::toString(transProbs[2]) +" out of bounds");
	}
#endif
}

/************************************************************************/
/*  updates the allows spread values at the time step                   */
/************************************************************************/
void OUAuxSpreadProcess::updateSpread(int timeStep, FDProduct::UpdateType type)
{
	try
	{
        bool isFwd = (type == FDProduct::FWD);
        /** lastStep = timeStep -1 for FWD and timeStep+1 for BWD */
        int lastStep = isFwd ? timeStep - 1 : timeStep +1;
		
        // for ease
        int finalStep = timeline->NumOfStep;
        // check is update is needed
        if(isFwd)
        {
             if(timeStep > 0 && fwds[timeStep] == fwds[lastStep]) return;
        }
        else
        {
            if(timeStep < finalStep && fwds[timeStep] == fwds[lastStep]) return;
        }
        
		
		// loop over tree nodes and assign spreads
		// spreads are floored at minSpreadVal
        // limits that we loop over depends on update direction
        // if going backwards only need to update cuurent limits since tree is contracting
        // if going fwds need update max range since tree is expanding
        int bot = isFwd ? botLimits[finalStep] : botLimits[timeStep];
        int top = isFwd ? topLimits[finalStep] : topLimits[timeStep];
		
        /** current value of fwd */
        double fwd = fwds[timeStep];
        int rIdx;
		for( rIdx = bot; rIdx <= top; rIdx++ )
		{
			spreads[sOffset(rIdx)] = Maths::max( minSpreadVal, fwd * qMap[sOffset(rIdx)] ); 
		}
	

	}
	catch(exception & e)
	{
		throw ModelException(e, "OUAuxSpreadProcess::updateSpread");
	}

}

DRLIB_END_NAMESPACE
