//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : SplineTrancheQuoteInterpolator.cpp
//
//   Description : Interpolation of CDOQuotes using mulit-spline approach
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/CCMPriceUtil.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/FlatExpLossPrior.hpp"
#include "edginc/GridFactory.hpp"
#include "edginc/SplineTrancheQuoteInterpolator.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

const double SplineTrancheQuoteInterpolator::DEFAULT_SMOOTHING = 0.99;

void SplineTrancheQuoteInterpolator::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SplineTrancheQuoteInterpolator, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ITrancheQuoteInterpolator);
    EMPTY_SHELL_METHOD(defaultSplineTrancheQuoteInterpolator);

	FIELD(prior, "Expected Loss prior. [Default: FlatExpLossPrior(0)]");
	FIELD_MAKE_OPTIONAL(prior);

	FIELD(smoothingWeight, "smoothing weight for spline in [0,1[. [Default = "+
		Format::toString(DEFAULT_SMOOTHING)+"]");

	FIELD_MAKE_OPTIONAL(smoothingWeight);

	FIELD(interpInterval,"Time interpolation discretisation interval (>= 1Q) [default = 1Y]");
	FIELD_MAKE_OPTIONAL(interpInterval)

	FIELD(interpDates ,"Override dates at which to interpolate [default: quote dates].");
	FIELD_MAKE_OPTIONAL(interpDates);
	
	FIELD(calibrateToIndex,"Calibrate to index swap spread (TRUE) or super senior quote (FALSE). [default = TRUE]");
	FIELD_MAKE_OPTIONAL(calibrateToIndex);

	FIELD(calibrateToMid, "Calibrate to Mid (TRUE) or to bid/ask quotes (FALSE). [default = TRUE]");
	FIELD_MAKE_OPTIONAL(calibrateToMid);

	FIELD(minSepTime,"Minimum (absolute) increase in expected loss of tranches between time points. [default = 10e-8]");
	FIELD_MAKE_OPTIONAL(minSepTime);

	FIELD(minSepStrike,"Minimum (absolute) differnce in expected loss of consecutive tranches (unit size). [default = 10e-8]");
	FIELD_MAKE_OPTIONAL(minSepStrike);

	FIELD(lastWeights, "Prior weights for loss points of tranches at last date in timeline. [default = 1.0]");
	FIELD_MAKE_OPTIONAL(lastWeights);
	

}

IObject* SplineTrancheQuoteInterpolator::defaultSplineTrancheQuoteInterpolator() {
    return new SplineTrancheQuoteInterpolator();
}

CClassConstSP const SplineTrancheQuoteInterpolator::TYPE = 
    CClass::registerClassLoadMethod("SplineTrancheQuoteInterpolator", 
                                    typeid(SplineTrancheQuoteInterpolator), 
                                    load);

/** private constructor */
SplineTrancheQuoteInterpolator::SplineTrancheQuoteInterpolator() : 
CObject(TYPE),
prior(new FlatExpLossPrior(0.0)),
smoothingWeight(DEFAULT_SMOOTHING),
interpInterval("1Y"),
interpDates(0), 
calibrateToIndex(true),
calibrateToMid(true),
minSepTime(0.0000001),
minSepStrike(0.0000001),
lastWeights(0)
{}

/** public constructor */
SplineTrancheQuoteInterpolator::SplineTrancheQuoteInterpolator(
		const DateTimeArray & interpDates,
		ITrancheQuoteInterpolatorSP prior,
		double smoothingWeight,
		bool calibrateToIndex,
		bool calibrateToMid,
		double minSepTime,
		double minSepStrike) : 
CObject(TYPE),
prior(prior),
smoothingWeight(smoothingWeight),
interpInterval("1Y"),
interpDates(new DateTimeArray(interpDates)),
calibrateToIndex(calibrateToIndex),
minSepTime(minSepTime),
minSepStrike(minSepStrike),
lastWeights(0)
{
	validatePop2Object();
}

/** Destructor */
SplineTrancheQuoteInterpolator::~SplineTrancheQuoteInterpolator()
{}

/** validate */
void SplineTrancheQuoteInterpolator::validatePop2Object()
{
	static const string method = "SplineTrancheQuoteInterpolator::validatePop2Object";
	try
	{
		// check that have at least 1 intepDate if included
		if(interpDates.get() !=NULL)
		{
			if(interpDates->size() < 1) throw ModelException("interpDates need to include at least one date");
		

			// check that interpDate are strictly increasing
			DateTime::ensureStrictlyIncreasing(*interpDates,
												   string("interpDates are not strictly increasing"),
												   true);
		}

		// check smoothingWeight
		if(smoothingWeight < 0 || smoothingWeight >=1) 
		{
			throw ModelException("smoothingWeight ("+Format::toString(smoothingWeight)+") is out of bounds [0,1[");
		}

		
		// check lastWeigts
		if(lastWeights.get() == NULL)
		{
			lastWeights = DoubleArraySP(new DoubleArray(1));
			(*lastWeights)[0] = 1.0;
		}


	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
/** create splineDates, i.e, node points of splines
	These include: valueDate, all quotes dates, user defined interpDates if given
	and schedule determined by interpInterval if available
 */
///////////////////////////////////////////////////////////////////////////////////////
DateTimeArraySP SplineTrancheQuoteInterpolator::createSplineDates(const CDOQuotes & quotes,
																  const ICDSParSpreads & index
																  ) const
{
	static const string method ="SplineTrancheQuoteInterpolator::createSplineDates";
	try 
	{

		int i;

		/** Collection of dates to merge to produce spline dates */
		vector<DateTimeArray> datesToMerge;

		/** value date */
		DateTime valueDate = quotes.getValueDate();
		
		// include valueDate in DatesToMerge
        DateTimeArray valDateArray(1);
        valDateArray[0] = valueDate;
        datesToMerge.push_back(valDateArray);

		
		/** tranche quotes object */
		CDOTrancheQuotesArrayConstSP trancheQuotes = quotes.getTrancheQuotes();

		// loop over all quote dates
		for(i =0;i< trancheQuotes->size();i++)
		{
			// add quote dates  to datesToMerge
			datesToMerge.push_back(*(*trancheQuotes)[i]->getMaturityDates());
		}
		
		
		/** index quote dates */ 
		int numIndexQuotes = index.getParSpreadsExpiries()->size(); 
		DateTimeArray indexQuoteDates(numIndexQuotes);
		
		for(i = 0;i<numIndexQuotes;i++)
		{
			indexQuoteDates[i] = (*index.getParSpreadsExpiries())[i]->toDate(valueDate);
		}

        // add to datesToMerge
        datesToMerge.push_back(indexQuoteDates);


		DateTimeArraySP intervalArray(0);
		if(interpDates.get() != NULL) // have explicit set of interpDates 
		{
			// include interpDates
			datesToMerge.push_back(*interpDates);
		}
		else //  Make interpInterval schedule if do not have explicit interpDates input
		{
			// schedule is from valueDate to last quote date with interpInterval
			// stub is at the beginnnig. This ensures that date are on IMM dates
			
			// find last quote date
			DateTime lastQuoteDate= valueDate;
			
			// first check through all tranche quote dates
			for(i =0;i< trancheQuotes->size();i++)
			{
				DateTime trancheDate = (*trancheQuotes)[i]->getMaturityDates()->back();
				
				if(trancheDate > lastQuoteDate) lastQuoteDate = trancheDate;
			}


			// now check through index swap dates
			DateTime indexDate = indexQuoteDates.back();
			if(indexDate > lastQuoteDate) lastQuoteDate = indexDate;


			// start date for schedule should be at least 1M away from start
			MaturityPeriod startInterval(1,"M");
			DateTime startDate = startInterval.toDate(valueDate);

			// decompose interpInterval into count and interval by first creating MaturityPeriod form string interpInterval
			int count;
			string interval;

			MaturityPeriod period(interpInterval);
			period.decompose(count, interval);
			
			/** interval schedule */
			intervalArray = DateTimeArraySP(
				SwapTool::dateArray(
					startDate,
					lastQuoteDate,
					count,
					interval,
					false				// stub at end
					)
				);
        

			// now replace first date in interval array (startDate) with value date 
			(*intervalArray)[0] = valueDate;

			// add to dates to merge
			datesToMerge.push_back(*intervalArray);  
		}
		
		
      
        // merge dates ---------------------------------------------------------
        return DateTimeArraySP(new DateTimeArray(DateTime::merge(datesToMerge)));

	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
	
}

//////////////////////////////////////////////////////////////////////////////////////
/**
	return surface of expected losses that reprice the quotes
	*/
//////////////////////////////////////////////////////////////////////////////////////
ExpectedLossSurfaceSP SplineTrancheQuoteInterpolator::getELSurface(
	   const CDOQuotes & quotes,
	   const ICDSParSpreads & index	
		) const
{
	static const string method ="SplineTrancheQuoteInterpolator::getELSurface";
	try 
	{

		int i,ti; 
		

		/** Dates at which spline values are defines */
		DateTimeArraySP splineDates = createSplineDates(quotes, index);
		

		// variables for ease ----------------------------------------------

		/** value date for quotes */
		DateTime valueDate = quotes.getValueDate();	

		/** number of interpolation dates */
		int N = splineDates->size();   // number of timeGrid points

		/** year fracs corresponding to splineDates */
 		DoubleArray timeGrid(N);
		for(ti = 0;ti<N;ti++)
		{
			timeGrid[ti] = valueDate.yearFrac((*splineDates)[ti]);
		}
		
		/** tranche quotes object */
		CDOTrancheQuotesArrayConstSP trancheQuotes = quotes.getTrancheQuotes();

		// create strikes for tranches that we are calibrating -------------------------------------
		
		int numTranches = trancheQuotes->size();
		/** low strikes ordered from low to high */
		DoubleArray lowStrike(numTranches);

		/** corresponding high strikes */
		DoubleArray highStrike(numTranches);
		
		// populate
		for(i = 0;i< numTranches;i++)
		{
			lowStrike[i] = (*trancheQuotes)[i]->getLowStrike();
			highStrike[i] = (*trancheQuotes)[i]->getHighStrike();
			
			if(highStrike[i] - lowStrike[i] <= 0 || highStrike[i] - lowStrike[i] > 1)
			{
				throw ModelException("trancheSize["+Format::toString(i)+
					"] ("+Format::toString(highStrike[i]-lowStrike[i])+") is out of bounds ]0,1]");
			}
		}
		
		/** array containing 0 and 1 */
		DoubleArray boundaryStrikes(2);
		boundaryStrikes[0] = 0.0;
		boundaryStrikes[1] = 1.0;

		/** strikes of tranches we are calibrating */
		// create strikes array by merging low, high and boundary strikes 
		DoubleArraySP strikes = GridFactory::mergeNoDuplicate(
			*GridFactory::mergeNoDuplicate(lowStrike,highStrike),
			boundaryStrikes
			);

		int numStrikes = strikes->size();
		/** low strikes of tranches that we are interpolating */
		DoubleArray lowS(numStrikes-1);
		/** high strikes of tranches that we are interpolating */
		DoubleArray highS(numStrikes-1);

		for(i = 0;i< numStrikes-1;i++)
		{
			lowS[i] = (*strikes)[i];
			highS[i] = (*strikes)[i+1];
		}

	
		// solve for spline curves
		ConstrainedCubicSplineSP spline = getSpline(
											quotes,
											index,
											lowS,
											highS,
											splineDates,
											timeGrid
											);
	
			
		return populateELSurface(lowS,highS,spline, splineDates, index.getRecovery());
	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
}
///////////////////////////////////////////////////////////////////////////////////////
/**
 get solved spline object 

 spline curves are normalised expected loss curves (i.e. for tranches of unit notional)
 we sove for one curve per {lowS,higS} pair

  each individual spline curve is defined at the timeGrid points
 
*/
/////////////////////////////////////////////////////////////////////////////////////
ConstrainedCubicSplineSP SplineTrancheQuoteInterpolator::getSpline( 
		  const CDOQuotes & quotes,
		  const ICDSParSpreads & index,
		  const DoubleArray & lowS,
		  const DoubleArray & highS,
		  DateTimeArraySP	splineDates,
		  const DoubleArray & timeGrid
		  ) const
{
	static const string method = "SplineTrancheQuoteInterpolator::getSpline";
	try
	{
		int i;
		/** number of strikes = number of tranches to interpolate */
		int numStrikes = lowS.size();

		/** Array of prior loss points */
		DoubleArray priorLosses(0);
		/** Array of prior weights */
		DoubleArray priorWeights(0);
		// call routine to define priorLosses and priorWeights array given the expectedLossPrior
		setPrior(
			quotes,
			index,
			lowS,
			highS,
			splineDates,
			priorLosses,
			priorWeights);

		/** smoothing parameter used in spline */
		double smoothing = smoothingWeight;
		
		
		//set up multi - spline class for numStrike number of curves --------------------------------------------
		ConstrainedCubicSplineSP spline(new ConstrainedCubicSpline(
			timeGrid,
			priorLosses,			 
			priorWeights,
			smoothing,
			numStrikes
			)
			);   

		// set zero constraints, i.e. EL(t = 0) = 0  for all curves-------------------------------------------------
		// assumes that first element in timeGrid is the valueDate which is the case by construction of splineDates

		DoubleArray targets(numStrikes, 0.0);
		IntArray targetIndexInGrid(numStrikes, 0);
		IntArray splineIndex(numStrikes);
		
		for(i = 0;i<numStrikes;i++)
		{
			splineIndex[i] = i;
		}
		spline->setTargets(targets, targetIndexInGrid, splineIndex);


	

		// set monotonicity constraints ---------------------------------------------------
		/** minimum separation between neighbouring points in time dimension */
		spline->setIncreasing(minSepTime); 
		

		// set remaining constraints -----------------------------------------------
		/** array to hold equality constraint values */
		DoubleArray eqVals(0);
		/** array to hold equality constraint functions */
		DoubleArray eqFuncs(0);

		/** array to hold inequality constraint values */
		DoubleArray ineqVals(0);
		/** array to hold inequality constraint functions */
		DoubleArray ineqFuncs(0);

		/** bid spread constraints */
		DoubleArray bidVals(0);
		DoubleArray bidFuncs(0);

		/** ask spread constraints */
		DoubleArray askVals(0);
		DoubleArray askFuncs(0);

		//get par spread constraints ------------------------------------------------------------
		
		// if we are calibrating to mid then we have an equality constraint
		if(calibrateToMid)
		{	
			// populate the spread constraints
			getSpreadConstraints(
				quotes,
				index,
				lowS,
				highS,
				splineDates,
				timeGrid,
				"MID",
				eqVals, 
				eqFuncs);
		}
		else // have two inequality constraints
		{
			
			// bid constraints
			getSpreadConstraints(
				quotes,
				index,
				lowS,
				highS,
				splineDates,
				timeGrid,
				"BID",
				bidVals, 
				bidFuncs);
			
			// add to inequality constraints
			ineqVals.insert(ineqVals.begin(), bidVals.begin(),bidVals.end());
			ineqFuncs.insert(ineqFuncs.begin(),bidFuncs.begin(),bidFuncs.end());

			// ask constraints
			getSpreadConstraints(
				quotes,
				index,
				lowS,
				highS,
				splineDates,
				timeGrid,
				"ASK",
				askVals, 
				askFuncs);
			
			// add to inequality constraints
			ineqVals.insert(ineqVals.begin(),askVals.begin(),askVals.end());
			ineqFuncs.insert(ineqFuncs.begin(),askFuncs.begin(),askFuncs.end());
		
		}
		
		
		// get concavity constraints ---------------------------------------------------
		DoubleArray concConstraintVals(0);
		DoubleArray concConstraintFuncs(0);

		getConcavityConstraints(
			lowS,
			highS,
			minSepStrike,			
			index.getRecovery(),
			splineDates->size(),
			concConstraintVals, concConstraintFuncs
			);
		// add to inequality constraints
		ineqVals.insert(ineqVals.begin(),concConstraintVals.begin(),concConstraintVals.end());
		ineqFuncs.insert(ineqFuncs.begin(),concConstraintFuncs.begin(),concConstraintFuncs.end());

		
		// get bound constraint: all (scaled) EL <= 1.0.
		// Only need to implement this for most junior tranche at final maturity because of the conc and increasing constratints. 
		DoubleArray boundVals(0);
		DoubleArray boundFuncs(0);

		getBoundConstraint(numStrikes, splineDates->size(), boundVals, boundFuncs);
		
		// add to inequality constraints
		ineqVals.insert(ineqVals.begin(),boundVals.begin(),boundVals.end());
		ineqFuncs.insert(ineqFuncs.begin(),boundFuncs.begin(),boundFuncs.end());

			
		// add remaining  constraints and solve -------------------------------------------------------
		spline->solve(
			eqFuncs,
			eqVals,
			ineqFuncs,
			ineqVals
			);
		
		return spline;
		

	}
	catch(IMSLException& eIMSL)
	{
		if(eIMSL.getErrorCode() == IMSL_SYSTEM_INCONSISTENT)
		{
			eIMSL.addMsg("Tranche Quote interpolation failed. Input quotes might be inconsistent.");
			eIMSL.addMsg("Try a lower recovery rate or finer time discretisation.");
			throw IMSLException(eIMSL);
		}
		else if(eIMSL.getErrorCode() == IMSL_NO_MORE_PROGRESS)
		{
			eIMSL.addMsg("Tranche Quote interpolation failed. Input quotes might be inconsistent.");
			eIMSL.addMsg("Try a different prior.");
			throw IMSLException(eIMSL);
		}
		else
		{
			throw ModelException("imsl_d_quadratic_prog failed with unrecognized error code ("+
				Format::toString(eIMSL.getErrorCode())+")");
		}
	}
	catch (exception & e)
	{
		throw ModelException(e,method);
	}
}
////////////////////////////////////////////////////////////////////////////////////
/*
   calculate output expected loss surface

  last curve is for strike 1-RR if have tranche above 1-RR or else last highStrike
  */
////////////////////////////////////////////////////////////////////////////////////
ExpectedLossSurfaceSP SplineTrancheQuoteInterpolator::populateELSurface(
		const DoubleArray & lowS,				// low strikes
		const DoubleArray & highS,				// high strikes
		ConstrainedCubicSplineSP  spline,	// spline that has been solved
		DateTimeArraySP	splineDates,				// dates at which spline is defined	
		double RR								// recovery rate
		) const
{
	static const string method = "SplineTrancheQuoteInterpolator::populateELSurface";
	try
	{
		int si,ti;


		// need contiguous strikes starting at 0 to calculate base curves
		if(lowS[0] != 0.0) throw ModelException("First tranche has to start at 0");
		for(si=1;si<lowS.size();si++)
		{
			if(lowS[si] != highS[si-1]) throw ModelException("HighStrike in tranche " +
				Format::toString(si) +
				"is not equal to lowStrike in tranche "+
				Format::toString(si-1)
				);
		}

		/** number of spline dates */
		int numTimes = splineDates->size();
		/** number of strikes */
		int numStrikes = 0;
		

		// base strikes for which we can calculate expected losses -----------------------------------
		DoubleArraySP strikes = DoubleArraySP(new DoubleArray(0));
		
		//find copy highStrikes STRICTLY before 1-RR
		int index = 0;
		while (index < highS.size() && highS[index] < 1-RR)
		{
			// copy strike STRICTLY before 1-RR
			strikes->push_back(highS[index]);
			index++;
		}
		
		
		// check if have tranche above or equal to  1-RR. in this case index < highS.size()
		if(index < highS.size())
		{
			// in this case we add 1-RR curve 
			strikes->push_back(1-RR);
			// also add 100% curve
			if(strikes->back() < 1.0) strikes->push_back(1.0);
		}
				
		numStrikes = strikes->size();

		// expected losses of base tranches ------------------------------------------------------
		DoubleArrayArraySP baseLosses = DoubleArrayArraySP(new DoubleArrayArray(numStrikes));

		DoubleArray normalisedTrancheEL(numTimes);

		// populate curves up to index
		for(si = 0;si<index;si++)
		{
			(*baseLosses)[si] = DoubleArray(numTimes);
			for(ti = 0;ti<numTimes;ti++)
			{
				// retrieve losses from spline
				spline->getValueArray(si, normalisedTrancheEL);

				(*baseLosses)[si][ti] = normalisedTrancheEL[ti]*(highS[si]-lowS[si]);

				if(si>0) (*baseLosses)[si][ti] += (*baseLosses)[si-1][ti];
			}
		}

		

		//if have tranches above 1-RR then add 1-RR and 100% tranche (these are the same)
		if(index < highS.size())
		{
			/** curves for strikes >= 1-RR are equal to finalCurve */
			DoubleArray finalCurve(numTimes);
			for(ti = 0;ti<numTimes;ti++)
			{
				// retrieve losses from spline
				spline->getValueArray(index, normalisedTrancheEL);
				finalCurve[ti] = normalisedTrancheEL[ti]*(highS[index]-lowS[index]);

				if(index>0) finalCurve[ti] += (*baseLosses)[index-1][ti];
			}

			// populate final two baseLosses
			for(si = index;si<numStrikes;si++)
			{
				(*baseLosses)[si] = finalCurve;	
			}
		}
		
		return ExpectedLossSurfaceSP(new ExpectedLossSurface(splineDates, strikes, baseLosses,RR));


	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
}


///////////////////////////////////////////////////////////////////////////////////////
/**
define prior array according to prior Type
*/
/////////////////////////////////////////////////////////////////////////////////////
void SplineTrancheQuoteInterpolator::setPrior( 
											  const CDOQuotes & quotes,
											  const ICDSParSpreads & index,
											  const DoubleArray & lowS,
											  const DoubleArray & highS,
											  DateTimeArraySP splineDates,		
											  DoubleArray & priorValues,
											  DoubleArray & priorWeights
											  ) const
{
	static const string method = "SplineTrancheQuoteInterpolator::setPrior";
	try
	{
		int i,j;
		// variables for ease ----------------------------------------------
		
		int N = splineDates->size();

		/** number of strikes */
		int numStrikes = lowS.size();

		/** number of elements in prior array */
		int length = N*numStrikes;

		priorValues.resize(length);

		priorWeights.resize(length);
	
		ExpectedLossSurfaceSP priorLossSurface = prior->getELSurface(
			quotes,
			index
			);
		
		double high;
		double low;
		DateTime date;
		//loop through strikes
		for(j = 0;j<numStrikes;j++)
		{
			high = highS[j];
			low = lowS[j];
			// loop through times and copy to prior
			for(i = 0; i< N;i++)
			{
				date = (*splineDates)[i];
				priorValues[j*N + i] = (priorLossSurface->getBaseEL(high,date) - priorLossSurface->getBaseEL(low,date))/(high-low);
				
				if(i == N-1) // special case for last date. Use lastWeights.
				{
					double lw = 1.0; // last weight
					if(j<lastWeights->size())
					{
						lw = (*lastWeights)[j];
					}
					priorWeights[j*N + i] = lw*(1.0 - smoothingWeight);
				}
				else
				{
					priorWeights[j*N + i] = 1.0 - smoothingWeight;
				}
				
			}
		}

	

	} catch (exception & e)
	{
		throw ModelException(e,method);
	}
}



///////////////////////////////////////////////////////////////////////////////////////////
/**
  get concavity constraints
*/
//////////////////////////////////////////////////////////////////////////////////////////
void SplineTrancheQuoteInterpolator::getConcavityConstraints(
		const DoubleArray & lowS,				// low strikes
		const DoubleArray & highS,				// high strikes
		double minSep,							// absolute minimum separation between tranche EL's (unit notional tranches)
		double RR,								// recovery rate
		int N,									// num splineDates
		DoubleArray & constraintVals,			// (O) constraint values
		DoubleArray & constraintFuncs			// (O) constraint functionals
		) const		
{
	static const string method ="SplineTrancheQuoteInterpolator::getConcavityConstraints";
	try
	{
		int i,k;
		int row;
		
		// variables for ease -------------------------------------------------	

		int numStrikes = lowS.size();
		
		/** number of elements in single spline vector */
		int singleSplineLength = 2*N-2;
		/** number of elements in multi-spline vector = number of columns in constraintFunc matrix */
		int numCols = singleSplineLength*numStrikes;
		
		/** number of rows in constraint matrix corresponding to flattend spreadContraintFuncs 
		This is the number of constraints: N-1 constraints for all curves except most senior one */
		int numRows = (numStrikes-1)*(N-1);

	
		constraintVals.resize(numRows);
		constraintFuncs.resize(numRows*numCols);

		// initialise to 0
		for(i = 0;i<constraintFuncs.size();i++) constraintFuncs[i] = 0.0;

		// loop through tranches and set constraints ---------------------------
		// last tranche is treated seperately if last highStrike > 1-RR

		
		// For tranche [k,l] with l>1-R we have tighter constraint
		// we have: unscaled EL[k,l] = unscaled EL[k,1-R]
		// This because we assume constant recovery and hence there can be no losses above 1-R
		// hence for the scaled losses we have: EL[k,1-R] = EL[k,l]*(l-k)/((1-R)-k)
		// We have condition that all other tranche EL > EL[k,1-R]	

		/** adjustment factor for last tranche if it has upper strike > 1-R  */
		double adjustFac = 1;
		
		double trancheSize;
		for(k = 0;k<numStrikes-1;k++)
		{
			// check that tranches are not overlapping
			if(lowS[k+1] < highS[k]) throw ModelException("Cannot have overlapping tranches");

			// set sdjust factor if needed. Check if next tranche has detach point > 1-RR
			if(highS[k+1] > 1-RR)
			{
				if(lowS[k+1] < 1-RR) // otherwise do nothing. Loss is 0 anyway
				{
					trancheSize = highS[k+1] - lowS[k+1];
					adjustFac = trancheSize/(1-RR-lowS[k+1]);
				}
			}
			//loop through timeline points starting at 2nd point
			for(i= 1;i<N;i++)
			{
				// row for constraint
				row = k*(N-1) + (i-1);
				//set constraint vals
				constraintVals[row] = minSep;
				
				// constraint x_{i}(k) - x_{i}(k+1) >= minSep
				constraintFuncs[row*numCols + k*singleSplineLength + i] = 1.0;
				constraintFuncs[row*numCols +(k+1)*singleSplineLength + i] = -adjustFac; // adjust tranche loss if needed
			}
		}
		
		
	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}

/////////////////////////////////////////////////////////////////////////////////////////
/**
  get Bound constraint
*/
/////////////////////////////////////////////////////////////////////////////////////////
void SplineTrancheQuoteInterpolator::getBoundConstraint(
			 int numStrikes,				// number of tranches
			 int N,							// num splineDates
			 DoubleArray & constraintVals,	// constraint values
			 DoubleArray & constraintFuncs	// constraint functionals
			 ) const
{
	static const string method ="SplineTrancheQuoteInterpolator:getBoundConstraint";
	try
	{	
		int i;
		// variables for ease -------------------------------------------------	
		
		/** number of elements in single spline vector */
		int singleSplineLength = 2*N-2;
		/** number of elements in multi-spline vector = number of columns in constraintFunc matrix */
		int numCols = singleSplineLength*numStrikes;
		
		/** number of rows in constraint matrix corresponding to flattend spreadContraintFuncs 
		Have 1 constraint for most equity tranche and highest maturity */
		int numRows = 1;

	
		constraintVals.resize(numRows);
		constraintFuncs.resize(numRows*numCols);

		// initialise to 0
		for(i = 0;i<constraintFuncs.size();i++) constraintFuncs[i] = 0.0;

		// set constraint ---------------------------
		// -x_{i}(k) >= -1.0
		
			
		constraintVals[0] = -1.0;
			
		constraintFuncs[N-1] = -1.0; // last maturity for first tranche
		
		
		
	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 populate the linear spread constraint and constraint values 

 Have one linear constraint per quote to calibrate.

 A linear constraint is defined by a constraint functional c(i) and constraint value v such that:

				\sum_{i=0}^{(2N-2)*numSplines} c(i)x(i) = v

  Here x[i] is the multi-spline vector. This consists of numSplines cppies of a single spline vector.

  a single spline vector if defined N values and N-2 second derivatives consecutively,
  where N is the number of timeGrid points

  
  */
////////////////////////////////////////////////////////////////////////////////////////////////////////
void SplineTrancheQuoteInterpolator::getSpreadConstraints(
		const CDOQuotes &  quotes,					// tranche quotes
		const ICDSParSpreads & index,				// index swap spreads
		const DoubleArray & lowS,					// low strikes		
		const DoubleArray & highS,					// high strikes 
		DateTimeArraySP	splineDates,				// splineDates
		const DoubleArray & timeGrid,				// interpolation time grid
		const string & quoteType,					// MID, BID or ASK
		DoubleArray & spreadConstraintVals,			// (O) constraint values	 [v]
		DoubleArray & spreadConstraintFuncs			// (O) constraint functionals [c(i)]
		) const
{
	static const string method ="SplineTrancheQuoteInterpolator::getSpreadConstraints";
	try
	{
		int k,i;
		
		// variables for ease -------------------------------------------------

		/** recovery rate */
		double RR = index.getRecovery();

		/** tranche quotes object */
		CDOTrancheQuotesArrayConstSP trancheQuotes = quotes.getTrancheQuotes();
	
		
		// dicount curve set up ----------------------------------------
		/** value date */
		DateTime valueDate = quotes.getValueDate();

		/** yield curve */
		YieldCurveConstSP yieldCurveSP = quotes.getDiscount();

		int numZeros= yieldCurveSP->zeroDates().size();
		/** times (Act.365 year fracs) at which have zero rates */
		DoubleArray yearFracsDf(numZeros);
		/** discount factors with associated with yearFracsDf */
		DoubleArray dfCurve(numZeros);

		// populate
		for(i = 0;i<numZeros;i++) 
		{
			dfCurve[i] = yieldCurveSP->pv(yieldCurveSP->zeroDates()[i]);
			yearFracsDf[i] = valueDate.yearFrac(yieldCurveSP->zeroDates()[i]);
		}

			


		// count number of quotes to calibrate ----------------------------------------------------
		
		/** number of input tranche curves which we calibrate to */
		int numTrancheTargetCurves = trancheQuotes->size();
		if(calibrateToIndex && trancheQuotes->back()->getHighStrike() == 1.0) 
		{
			// in this case the last curve is a super senior curve. Since we are calibrating to index it is ignored
			numTrancheTargetCurves -= 1;
		}

		/** total number of quotes that we are calibrating to */
		int totalNumQuotes = 0;
		
		// first add standard quotes
		for(i = 0;i< numTrancheTargetCurves;i++)
		{
			totalNumQuotes += (*trancheQuotes)[i]->getMaturityDates()->size();
		}


		// when calibrating to index swap spread need to add index quotes
		/** quote dates for index swap quotes (used later) */
		DateTimeArray indexQuoteDates(0);
		/** number of index quotes */
		int numIndexQuotes = 0; 
		
		if(calibrateToIndex) 
		{
			numIndexQuotes = index.getParSpreadsExpiries()->size(); 
			totalNumQuotes += numIndexQuotes; 
			
			indexQuoteDates.resize(numIndexQuotes);
			for(i = 0;i<numIndexQuotes;i++)
			{
				indexQuoteDates[i] = (*index.getParSpreadsExpiries())[i]->toDate(valueDate);
			}
		}
		
	
		
		// set dimensions of output arrays and initialise -------------------------------------------------------
	
		int N = timeGrid.size();
		/** number of elements in single spline vector */
		int singleSplineLength = 2*N-2;
		
		
		/** number of tranches to interpolate */
		int numInterpTranches = highS.size();
		/** number of elements in multi-spline vector = number of columns in constraint matrix */
		int numCols = singleSplineLength*numInterpTranches;
		
		/** number of rows in constraint matrix corresponding to flattend spreadContraintFuncs 
		This is the number of constraints */
		int numRows = totalNumQuotes;
	

		// resize output arrays 
		spreadConstraintVals.resize(numRows);
		spreadConstraintFuncs.resize(numRows*numCols);

		// initialise to 0;
		for(i = 0;i<spreadConstraintFuncs.size();i++) spreadConstraintFuncs[i] = 0;


		// define auxilliary variables needed later -------------------------------------------	
		
		/** start index for given curve in row of constraintFuncs */
		// for each tranche have numQuotes*numCols constraints
		// each constraint is numCols long
		int firstConstraintCol;


		/** spread to calibrate to */
		double targetSpread;
		/** upfront to calibrate to */
		double targetUpfront;
		/** date for quote we are calibrating to */
		DateTime quoteDate;

		/** row in constraint matrix we are at. This just counts the number of constraints */
		int row = 0;
	
		// loop through quoted tranches and add constraints --------------------------------------------------
		// we have one constraint for each quote for each tranche
		for(k = 0;k<numTrancheTargetCurves;k++)
		{
			
			/** number of quotes in this tranche */
			int numQuotes = (*trancheQuotes)[k]->getMaturityDates()->size();
			
			// loop through quote for each tranche
			for(i = 0; i< numQuotes;i++)
			{
				// position in contraintFunc vector to populate constraint
				firstConstraintCol = row*numCols;
				
				// quote data
				quoteDate = (*(*trancheQuotes)[k]->getMaturityDates())[i];
				
				if(quoteType=="MID")
				{
					targetSpread = (*trancheQuotes)[k]->getSpread(quoteDate);
					targetUpfront = (*trancheQuotes)[k]->getUpfront(quoteDate);
				}
				else if(quoteType=="ASK")
				{
					if((*trancheQuotes)[k]->hasBidAsk())
					{
						targetSpread = (*trancheQuotes)[k]->getAskSpread(quoteDate);
						targetUpfront = (*trancheQuotes)[k]->getAskUpfront(quoteDate);
					}
					else
					{
						throw ModelException("calibrateToMid is FALSE but there are no bid/ask spreads for tranche ("+Format::toString(k)+")");
					}
				}
				else if(quoteType=="BID") 
				{
					if((*trancheQuotes)[k]->hasBidAsk())
					{
						targetSpread = (*trancheQuotes)[k]->getBidSpread(quoteDate);
						targetUpfront = (*trancheQuotes)[k]->getBidUpfront(quoteDate);
					}
					else
					{
						throw ModelException("calibrateToMid is FALSE but there are no bid/ask spreads for tranche ("+Format::toString(k)+")");
					}
				}
				else throw ModelException("quoteType "+quoteType+" is not supported");
				
	
				/** genertate a generic fee leg with unit coupon and 0 upfront for quote instrument */
				// Dynamic_cast is temporary until proper PV method has been included in CreditFeeLeg
				ICreditFeeLegSP feeLeg = quotes.generateFeeLegOverride(1,0,quoteDate);
                CreditFeeLegWithPVSP indexFeeLegSP = 
					CreditFeeLegWithPVSP::dynamicCast(feeLeg);

				// Now set the constraint in the relevant position in spreadConstraintFuncs
				// Also return the constraint value
				spreadConstraintVals[row] = linearSpreadConstraint(
					targetSpread,		
					targetUpfront,
					quoteDate,
					(*trancheQuotes)[k]->getLowStrike(),
					(*trancheQuotes)[k]->getHighStrike(),
					lowS,
					highS,
					splineDates,
					timeGrid, 
					yieldCurveSP,
					yearFracsDf,
					dfCurve,
					indexFeeLegSP,
					RR,
					quoteType,
					firstConstraintCol,
					spreadConstraintFuncs
					);
				
				// update row to next constraint
				row += 1;
			} // end loop over quote dates
		} // end loop over tranches
			
		 // Add index swap constraints if calibrating to index
		if(calibrateToIndex)
		{
			double lowStrikeTarget = 0.0; // low strike for index  = [0,1] tranche
			double highStrikeTarget = 1.0; // high strike for index = [0,1] tranche

			// loop through quotes used to calibrate ss tranche
			for(i = 0; i < numIndexQuotes;i++)
			{				
				// position in contraintFunc vector to populate constraint
				firstConstraintCol = row*numCols; 
				
				// quote details
				quoteDate = indexQuoteDates[i];
				targetSpread = (*index.getParSpreads())[i];
				targetUpfront = 0;
	

				/** genertate a generic fee leg with unit coupon and 0 upfront for quote instrument */
				// Dynamic_cast is temporary until proper PV method has been included in CreditFeeLeg
                ICreditFeeLegSP feeLeg = quotes.generateFeeLegOverride(1,0,quoteDate);
				CreditFeeLegWithPVSP indexFeeLegSP = 
					CreditFeeLegWithPVSP::dynamicCast(feeLeg);

				// Now set the constraint in the relevant position in spreadConstraintFuncs
				// Also return the constraint value
				spreadConstraintVals[row] = linearSpreadConstraint(
					targetSpread,		
					targetUpfront,
					quoteDate,
					lowStrikeTarget,
					highStrikeTarget,
					lowS,
					highS,
					splineDates,
					timeGrid, 
					yieldCurveSP,
					yearFracsDf,
					dfCurve,
					indexFeeLegSP,
					RR,
					quoteType,
					firstConstraintCol,
					spreadConstraintFuncs
					);

				// update row to next constraint
				row +=1;
			} // end loop over quotes
		} // end index swap calibration

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  linearSpreadConstraint 

  fills in (2N-2)*numStrikes values in constraintArray

  Note on linearSpreadConstraint:

  The constraintArray is a flattened matrix where the number of rows = number of constraints
  linearSpreadConstraint populates one of these `rows'. This row corresponds to a range in constraintrray
  starting at startIndex.

  Let the index i correspond to the elements in constraintArray that correspond to the `row' (constraint) of interest
  then the constraintVal and constraintArray are defined such that:

  constraintVal = \sum_i x[i]*constraintArray[i] 

  here x[i] are the spline values. A single spline is parameterized by N values and N-2 second derivatives consecutively
  Since we are dealing with mulit-spline, x[i] corresponds to consecutive spline parametrisations for all tranches.
  Hence x[i] is an array of length numCols = numInterpTranches*(2N-2).
  This hence also the number of elements of constraintArray that are populated. In other words we populate the indicies
  in the range i \in [startIndex, startIndex + numCols [.

  linearSpreadConstraint imposes the constraint that the tranche curves we are solving for satisfy the
  given tranche spread and upfront at given maturityIndex.
		

  Note that here expected loss always means expected loss of the contingent leg, which for the SS tranche is not
  the same as the expected loss on the fee leg. The two can be related if we assume constant rcovery.

  We find the constraintArray in the following way:
  We note first that the tranche both fee and contingent legs are linear in x[i].
  (Note that the x[i] represent the contingent leg expected losses / tranche notional size for each tranche)

  
  i.e we have the constraint contLeg(x[i]) - spread*feeLeg(x[i]) = upfront.

  which because of linearity can be written:

  \sum_i contLegFactor[i]*x[i] -spread*(feeLegFactor[i]*x[i]+ constFactor) = upfront. 

  (note that the contingentLeg has no constFac since contLeg = 0 if all expected losses are 0)

  Hence we deduce: constraintVal[i] = contLegFactor[i] - spread*feeLegFactor[i]
  and constraintVal = upfront + spread*constFactor

  To determine the various factors we value the contingent and feelegs with specail choices of x[i].

  To obtain the constFactor we choose x[i] = 0 for all 0. We then have feeLeg(x[.]) = constFactor.

  We now choose  x[i] = 1 if i = j and 0 otherwise.
  we then have contLeg(x[.]) = contLegFactor[j] and feeLeg(x[.]) = feeLegFactor[j] + constFactor.

  As a last point we note that in practice the contLegs and feeLegs are defined as functions of the effectiveCurve
  
	The relation between the x[i] and the effective curve is given by the function effCurveValue().

 */
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double SplineTrancheQuoteInterpolator::linearSpreadConstraint(
		double spread,							// target spread
		double upfront,							// target upfront
		const DateTime & quoteDate,			    // date corersponding to quotes
		double lowStrikeTarget,					// low strike corresponding to quotes
		double highStrikeTarget,				// high strike corresponding to quotes
		const DoubleArray & lowS,				// low strikes
		const DoubleArray & highS,				// high strikes
		DateTimeArraySP splineDates,			// spline Dates
		const DoubleArray & timeGrid,			// time grid for interpolation
		 YieldCurveConstSP & yieldCurveSP,		// yield curve
		 const DoubleArray & yearFracsDf,		// Act/365F year fraction corresponding to zeroDates
		 const DoubleArray & dfCurve,			// array of discount factors defined at yearFracsDf times
		 CreditFeeLegWithPVSP indexFeeLegSP,	// index fee leg
		 double RR,								// recovery
		const string & quoteType,				// MID, BID or ASK
		int startIndex,							// first index in constraintArray to populate. This should be start of constraintFunc row
		DoubleArray & constraintArray			// (O) constraint array to populate
		) const
{
	static const string method ="SplineTrancheQuoteInterpolator::linearSpreadConstraint";

	try
	{
        // Require a model for calculating fees
        // so just provide the default closed form
        ClosedFormForwardRatePricerSP cfPricer =
            ClosedFormForwardRatePricerSP(
                new ClosedFormForwardRatePricer());

		int i,k;
		int N = timeGrid.size();
		int singleSplineLength = 2*N-2;
		int numStrikes = lowS.size();
		
		DateTime valueDate = splineDates->front();


		/** index of quoteDate in splineDates */
		int maturityIndex = quoteDate.find(*splineDates);

		/** factor to take into account whether we are using mid, bid or ask
		* have MID: \sum constraintFunc[i]x[i] = constraintVal
		* BID:		\sum constraintFunc[i]x[i] >= constraintVal
		* ASK:		\sum -constraintFunc[i]x[i] >= - constraintVal */
		double quoteFac = 1.0;
		if(quoteType == "ASK") quoteFac = -1.0;


		//validation -----------------------------------------------------------
		if(startIndex <0) throw ModelException("startIndex<0");
		if(startIndex > constraintArray.size() - (2*N-2)) throw ModelException("startIndex too high");
		if(maturityIndex < 0 || maturityIndex > N-1) 
		{
			throw ModelException("maturityIndex "+ Format::toString(maturityIndex) + 
				" outside [0, "+Format::toString(N)+"]");
		}
	
					
		
		/** effective curve effCurve(t) is defined as expected notional (= 1-exploss) at t divided by the tranche size 
			if calibrateToIndex = true then it is the effective curve of the [0,100] tranche 
			otherwise it is the effective curve of the super senior tranche */
		DoubleArray effCurve(N);	
		

		// initialise (all epected losses are 0)
		for(i = 0; i< N;i++)  effCurve[i] = 1.0; 

		/** fee leg when all expected losses = 0 */
		double constFactor = indexFeeLegSP->pvRisky(
			(*splineDates), 
			  effCurve, 
			  valueDate,
			  yieldCurveSP,
              cfPricer);
	


		/** effective curve value of contingent leg */
		double contLegValue;
		/** effective curve value of fee leg */
		double feeLegValue;
		
		/** linear factor for contingent leg */
		double contLegFactor; 
		/** linear factor for fee leg */
		double feeLegFactor;

		// loop through tranches -----------------------------------------------------------
		for(k = 0; k< numStrikes; k++) 
		{
			
			// get effective curve values of contingent and fee legs (this is indep of maturity)
			effCurveValue(
				lowStrikeTarget,
				highStrikeTarget,
				lowS[k],
				highS[k],
				RR,
				contLegValue,
				feeLegValue
				);

			/** protection start date T+1 for now until have proper contingent leg object */
			double protectionStart = valueDate.yearFrac(valueDate.rollDate(1));
			
			// loop through maturities ----------------------------------------------------
			for(i = 0; i < maturityIndex+1 ;i++)
			{

				// cont leg factor ---------------------------------------------
				
				effCurve[i] = contLegValue;
				
				if(effCurve[i] == 1.0)
				{
					// in this case there are no losses and the contingent leg value is 0
					contLegFactor = 0;
				}
				else
				{
					contLegFactor = CCMPriceUtil::contingentPrice(
												protectionStart,			// protection start
												timeGrid[maturityIndex],	// protection end
												0,							// delay
												yearFracsDf, dfCurve,
												timeGrid, effCurve,
												CCMPriceUtil::LINEAR);
				}
				
			
				// fee leg factor  ----------------------------
			
				effCurve[i] = feeLegValue;

				if(effCurve[i] == 1.0)
				{
					// in this case there are no losses and fee leg = constFactor
					feeLegFactor = constFactor;
				}
				else
				{
					feeLegFactor = indexFeeLegSP->pvRisky(
						(*splineDates), 
						  effCurve, 
						  valueDate,
						  yieldCurveSP,
                          cfPricer);
				}
				
				// need to take off the constFactor
				feeLegFactor -= constFactor;

				//reset effCurve
				effCurve[i] = 1.0;


				// set coefficient in constraintArray ---------------
				constraintArray[startIndex + k*singleSplineLength + i] = (contLegFactor - spread*feeLegFactor)*quoteFac;
			
			} // end loop over maturities

		} // end loop over tranches
		
					

		// constraintVal ---------------------------------------------------------------------------------
		// this is given by upfront + spread*feeleg(1 effCurve)
		double constraintVal = (upfront + (spread)*constFactor)*quoteFac;

		return constraintVal;

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/**
function that produces effective curve value for one time t(it doesn't depend on which time point)
of a chosen tranche with strikes: lowStrikeTarget, highStrikeTarget.

  Gives the effective curve values for both the contingent and fee legs.

We assume that the contingent leg effective curve values of all tranche are 1 at t except for the (test) tranche given by:
lowStrikeTest and highStrikeTest

For the test tranche the contingent leg effective curve value is assumed to be 0.
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////
void SplineTrancheQuoteInterpolator::effCurveValue(
						double lowStrikeTarget,		/** low strike of tranche for which we want effective curve */
						double highStrikeTarget,	/** high strike of tranche for which we want effective curve */
						double lowStrikeTest,		/** low strike of tranche which has 0 effective curve val */
						double highStrikeTest,		/** high strike of tranche that has 0 effective curve val */
						double RR,					/** recovery rate */
						double & contLegValue,		/** (O) contingent leg eff curve value */
						double & feeLegValue		/** (O) fee leg eff curve value */
						) const
{

	static const string method = "SplineTrancheQuoteInterpolator::effCurveValue";
	try
	{
		/** target tranche size */
		double trancheSizeTarget = highStrikeTarget - lowStrikeTarget;
		if(trancheSizeTarget <= 0) throw ModelException("trancheSizeTarget <= 0");

		/** test tranche size */
		double trancheSizeTest = highStrikeTest - lowStrikeTest;
		if(trancheSizeTest <= 0) throw ModelException("trancheSizeTest <= 0");

		

		// we are assuming that lowS and highS describe consecutive tranches 
		// These include lowStrikeTarget, hihgStrikeTarget and lowStrikeTest, highStrikeTest

		// contingent leg ---------------------------------------------------------------
		if(lowStrikeTest >= highStrikeTarget || highStrikeTest <= lowStrikeTarget)
		{
			// in this case the test tranche is outside the target tranche.
			// hence the target tranche does not depend on the test tranche
			contLegValue = 1.0;
		}
		else
		{
			// the target tranche is a sum of consecutive tranche, one of which is the test tranche.
			// since all tranches expcept the test tranche are assumed to have value 1 we have:
			if(trancheSizeTest > trancheSizeTarget) throw ModelException("trancheSizeTest () > trancheSizeTarget ()");
			contLegValue = 1 - trancheSizeTest/trancheSizeTarget;
		}

		// now consider fee leg -------------------------------------------------
		if(highStrikeTarget <= 1-RR)
		{
			// in this case fee leg and cont leg effective curves are the same
			feeLegValue = contLegValue;
		}
		else
		{

			// demand that highStrikeTarget = 1
			if(highStrikeTarget < 1.0) throw ModelException("Only high strike > 1-RR allowed is 1");
			// demand that lowStrikeTarget <= 1-RR
			if(lowStrikeTarget > 1-RR) throw ModelException("If highStrikeTarget = 1 need corresponding low strike to be < 1-RR");


			if(RR<=0) throw ModelException("RR<=0"); // Note that if RR=0 this should have been caught in previous case since highStrikeTarget <=1
			// convert contLegValue to feeLegValue using functionin ExpectedLossSurface

			// need baseELLow = EL[0,A], A = (1-RR)/RR*(1-lowStrikeTarget) >= 1-RR since lowStrikeTarget <= 1-RR
			// and baseELHigh = EL[0,B] , B = (1-RR)/RR*(1-highStrikeTarget) = 0 since highStrikeTarget = 1
			
			// since A >= 1-RR EL[0,A] =  EL[0,1] = EL[testTranche]  = trancheSizeTest since we are assuming effCurve = 0 for this tranche 
			// this is since all tranche except tranche size test are assumed to have a effCurve value of 1
			// whereas the testTranche has an effCurve value of 0
			double baseELLow = trancheSizeTest;
			
			// since B = 0 => baseELHigh = 0
			double baseELHigh = 0;
				

			feeLegValue = ExpectedLossSurface::feeLegEffCurve(
				(1-contLegValue)*trancheSizeTarget,  // unscaled expected loss of target tranche
				trancheSizeTarget,
				baseELLow,
				baseELHigh,
				RR);
		}
													

	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}





	
/////////////////////////////////////////////////////////////////////////////////////////////////////


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool SplineTrancheQuoteInterpolatorLoad() {
    return (SplineTrancheQuoteInterpolator::TYPE != 0);
}


DRLIB_END_NAMESPACE


