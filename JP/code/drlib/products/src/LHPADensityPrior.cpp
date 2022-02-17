//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : LHPADensityPrior.cpp
//
//   Description : Probability density prior using large homogeneous pool approximation
//					Currently supports Gaussian Copula
//
//   Author      : Matthias Arnsdorf
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/LHPADensityPrior.hpp"
#include "edginc/GridFactory.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/mathlib.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE


/** maximum value for any exponent  */
const double LHPADensityPrior::maxExp = 500;

/** high bound for density */
const double LHPADensityPrior::huge = exp(maxExp);

/** default numIntervals */
const int LHPADensityPrior::def_numIntervals = 200;

/** default sigma */
const double LHPADensityPrior::def_sigma = 6;

void LHPADensityPrior::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(LHPADensityPrior, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IDensityPrior);
    EMPTY_SHELL_METHOD(defaultLHPADensityPrior);

	//FIELD(lowBound, "Lower bound of range on which density is defined [default = 0]");
	//FIELD_MAKE_OPTIONAL(lowBound);
	
	FIELD(highBound, "Upper bound of range on which density is defined [default = 1]");
	FIELD_MAKE_OPTIONAL(highBound);

	FIELD(numIntervals, "Number of deiscretiation intervals between lowBound and highBound [default = "+
		Format::toString(def_numIntervals) + "]");
	FIELD_MAKE_OPTIONAL(numIntervals);

	FIELD(sigma, "Number of standard deviations for which we sample cummulative normal. [default ="+
		Format::toString(def_sigma) + "]");
	FIELD_MAKE_OPTIONAL(sigma);

	FIELD(betaDates, "Dates for betas");
	
	FIELD(betas, "Term structure of betas");

	FIELD(floor, "Floor for the density [default = 0]");
	FIELD_MAKE_OPTIONAL(floor);

	

	FIELD(index, "Index swap quotes");
	
	FIELD(market, "");

	FIELD(grid,"");
	FIELD_MAKE_TRANSIENT(grid);

	FIELD(gridDate,"");
	FIELD_MAKE_TRANSIENT(gridDate);

	
	FIELD(cummDensity,"");
	FIELD_MAKE_TRANSIENT(cummDensity);


}

IObject* LHPADensityPrior::defaultLHPADensityPrior() {
    return new LHPADensityPrior();
}

CClassConstSP const LHPADensityPrior::TYPE = 
    CClass::registerClassLoadMethod("LHPADensityPrior", 
                                    typeid(LHPADensityPrior), 
                                    load);


/** private constructor */
LHPADensityPrior::LHPADensityPrior() : 
CObject(TYPE), 
lowBound(0), 
highBound(1), 
numIntervals(def_numIntervals),
sigma(def_sigma),
floor(0),
grid(0),
gridDate(0),
cummDensity(0)
{
}

/** public constructor */
LHPADensityPrior::LHPADensityPrior(double lowBound, 
								   double highBound, 
								   int numIntervals,
								   double sigma,
								   double floor) : CObject(TYPE), 
lowBound(lowBound), 
highBound(highBound), 
numIntervals(numIntervals),
floor(floor),
sigma(sigma),
grid(0),
gridDate(0),
cummDensity(0)
{
	validatePop2Object();
}

/** Destructor */
LHPADensityPrior::~LHPADensityPrior()
{}

/** validate */
void LHPADensityPrior::validatePop2Object()
{

	

	// Dummy model 
	smartPtr<NonPricingModel> dummy(new NonPricingModel());
	
	// get market data -----------------------------------------
	index.getData(dummy.get(), market.get());

	if(lowBound != 0.0)
	{
		throw ModelException("lowBound needs to be 0");
	}
	if(lowBound >= highBound)
	{
		throw ModelException("lowBound (" + Format::toString(lowBound) +
			") >= highBound (" + Format::toString(highBound)+ ")");
	}
	if (numIntervals<1) 
	{ 
		throw ModelException("numIntervals ("+Format::toString(numIntervals)+
			") needs to be > 0");
	}

	if(betas->size() != betaDates->size())
	{
		throw ModelException("betas and betaDates have different length");
	}

	DateTime::ensureStrictlyIncreasing(*betaDates,
		"BetaDates not strictly increasing",
		true
		);

	if(sigma<=0)
	{
		throw ModelException("sigma needs to be strictly positive");
	}

	
}

/** return grid points between which density is linear. */
DoubleArraySP LHPADensityPrior::getGridPoints(
	double start, 
	double end,
	const DateTime & date)
{
	
	static const string method ="LHPADensityPrior::getGridPoints";
	try
	{
		if(start > end) throw ModelException("Need start <= end");
		
		
		// calc grid if not already done so
		if(gridDate.get() == NULL || *gridDate != date)
		{
			calcGridPoints(date);
		}


		/** output array */
		DoubleArraySP out = DoubleArraySP(new DoubleArray(0));

		//loop through grid
		int i;
		for(i=0; i < grid->size(); i++)
		{
			// grid point
			double x = (*grid)[i];
			// check if in range [start,end]
			if(x >=start && x <=end)
			{
				out->push_back(x);
			}
		}
	
		
		return out;
	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
}

/** calculate the grid points at date */
void LHPADensityPrior::calcGridPoints(const DateTime & date)
{
	static const string method ="LHPADensityPrior::calcGridPoints";
	try
	{

		gridDate = DateTimeSP(new DateTime(date));

		// first get beta at date needed for grid calculation
		// do linear interpolation or flat extrapolation of beta
		double beta;
		int lowI = date.findLower(*betaDates);

		if(lowI < 0) 
		{
			beta = (*betas)[0];
		}
		else if(lowI >= betas->size()-1)
		{
			beta = betas->back();
		}
		else
		{
			beta = (*betas)[lowI];
			beta += ((*betas)[lowI+1] - (*betas)[lowI])/(*betaDates)[lowI].yearFrac((*betaDates)[lowI+1]);
		}

		// cap beta
		beta = Maths::min(beta, 0.999);

		
		// are assuming that lowBound is 0;
		if(lowBound != 0.0) throw ModelException("lowBound needs to be 0");
		if(highBound == 0.0) throw ModelException("highBound needs to be > 0");



		/** default prob */
		double prob = 1 - index.getSP()->survivalProb(date);

		// Ninv of prob
		double nInvP = N1InverseBetter(prob);

		double y;
		double phi;

		int i;

#if 0 // this mehtod is currently not used
		// build grid by using inverse of LHPA cumm density
		// calc cumulative density at density grid 
		
		
		// construct homogeneous array of cumm density values. 
		cummDensity = GridFactory::homogeneous(0,1.0,numIntervals+1);

		// initialise grid
		grid = DoubleArraySP(new DoubleArray(cummDensity->size()));
		for(i = 0; i< grid->size(); i++)
		{
			/** cumm density Point */
			phi = (*cummDensity)[i];

			// check boundary cases
			if(phi == 0.0)
			{
				(*grid)[i] = lowBound; // grid starts at 0;
			}
			else if(phi==1.0)
			{
				(*grid)[i] = 1.0*(highBound-lowBound);
			}
			else
			{
				// phi = N((sqrt(1-beta^2)*NInv(x/(hB-lB) - nInv(p))/beta)
				// invert this to get x
				y = N1InverseBetter(phi);	

				(*grid)[i] = (highBound-lowBound)*N1((beta*y + nInvP)/sqrt(1-beta*beta));
			}
		}
#endif

		// method using homogeneous grid between given standard deviation around centre
		
		/** lower range bound */
		double low = N1((-sigma*beta + nInvP)/sqrt(1-beta*beta))*(highBound-lowBound);
		/** upper range bound */
		double high = N1((sigma*beta + nInvP)/sqrt(1-beta*beta))*(highBound-lowBound);

		// homogeneous grid
		grid = GridFactory::homogeneous(low,high,numIntervals+1);

		// add lowBound and highBound points
		//note that lowBound forced to 0 hence 0< lowBound and high < highBound
		grid->push_back(highBound);
		grid->insert(grid->begin(),lowBound);

		// calc cummDensity at gridPoints
		cummDensity = DoubleArraySP(new DoubleArray(grid->size()));
		double x;
		for(i = 0;i< cummDensity->size();i++)
		{
			x = (*grid)[i];
			y = N1InverseBetter(x/(highBound-lowBound));
			(*cummDensity)[i] =  N1((sqrt(1-beta*beta)*y - nInvP)/beta);
		}

	} catch (exception & e)
	{
		throw ModelException(e, method);
	}
}


/** return density at gridPoints */
DoubleArraySP LHPADensityPrior::getDensity(const DoubleArray & gridPoints, 
											  const DateTime & date			
											  ) const
{
	static const string method ="LHPADensityPrior::getDensity";
	try
	{
		// check if gridPoints have been calculated and corresponding to date
		if(gridDate.get() == NULL || *gridDate != date)
		{
			throw ModelException("gridPoints have not been calculated");
		}


		// construct piecewise constant density that coincides with cummDensity
		/** density function */
		DoubleArray density;
		density.resize(grid->size());

		if(grid->size()<2) throw ModelException("grid->size() <2");

		/** step size in grid */
		double gridStep = (*grid)[1] - (*grid)[0];
		if(gridStep == 0.0) throw ModelException("gridStep=0");
		
		/** value of integrated density over segment */
		double intDense = ((*cummDensity)[1] - (*cummDensity)[0]);

		// first segment   
		// approximate density as constant in first period
		density[0] = intDense/gridStep;
		
		int i;
		for(i = 1; i < density.size();i++)
		{
			gridStep = (*grid)[i] - (*grid)[i-1];
			if(gridStep == 0.0) throw ModelException("gridStep=0");
			intDense = (*cummDensity)[i] - (*cummDensity)[i-1];

			// chose density so that if assume piecewise linear between points will give
			// same integrated value
			// (density[i]+density[i-1])/2*gridStep = intDense
			density[i] = 2*intDense/gridStep - density[i-1];

			// floor at floor
			density[i] = Maths::max(floor, density[i]);
			
		}


		/** output density array */
		// these are linearly interpolated density points at input gridPoints
		DoubleArraySP out(new DoubleArray(gridPoints.size()));
		
		
		if(grid->size() <2) throw ModelException("need at least 2 grid elements");
		
		double x;
		int lowI, highI;
		// populate out array 
		for(i = 0 ; i < out->size(); i++) 
		{
			x = gridPoints[i];
			if(x< lowBound || x > highBound)
			{
				throw ModelException("gridPoint at index "+ 
					Format::toString(i)+ " (" + 
					Format::toString(x) +") is not in range [" +
					Format::toString(lowBound) +", " +
					Format::toString(highBound) +"]");
			}

			// find lowest index highI such that grid[highI] => gridPoints[i]
			highI = distance(grid->begin(), lower_bound(grid->begin(),grid->end(), x) );
			
			if(highI >= grid->size() || highI < 0)
			{
				throw ModelException("highI is out of bounds"); // this shouldn' t happen
			}

			lowI = highI-1;
			if(lowI < 0) lowI = 0;

			double outTemp = LinearInterp(
				x,
				(*grid)[lowI],
				(*grid)[highI], // can cope with grid points being the same
				density[lowI],
				density[highI]
				);	

			(*out)[i] = outTemp;
		}

		return out;
	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
}


/** set range on which density is defined */
void LHPADensityPrior::setRange(double lowB, double highB)
{

	static const string method ="LHPADensityPrior::setRange";
	try
	{
		lowBound = lowB;
		highBound = highB;

		if(lowBound != 0.0) throw ModelException("lowBound needs to be 0");

		if(lowBound >= highBound) 
		{
			throw ModelException("lowBound ("+ Format::toString(lowBound) +
				") >= highBound (" + Format::toString(highBound) + ")");
		}

		// calc new grid
		grid = GridFactory::homogeneous(lowBound,highBound,numIntervals+1);
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}



/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool LHPADensityPriorLoad() {
    return (LHPADensityPrior::TYPE != 0);
}


DRLIB_END_NAMESPACE
