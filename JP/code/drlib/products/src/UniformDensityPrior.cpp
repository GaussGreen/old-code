//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : UniformDensityPrior.cpp
//
//   Description : Uniform probability density prior
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/UniformDensityPrior.hpp"

DRLIB_BEGIN_NAMESPACE



/** private constructor */
UniformDensityPrior::UniformDensityPrior() : CObject(TYPE), lowBound(0), highBound(1)
{}

/** public constructor */
UniformDensityPrior::UniformDensityPrior(double lowBound, double highBound) : CObject(TYPE), 
lowBound(lowBound), highBound(highBound)
{
	validatePop2Object();
}

/** Destructor */
UniformDensityPrior::~UniformDensityPrior()
{}

/** validate */
void UniformDensityPrior::validatePop2Object()
{
	if(lowBound >= highBound)
	{
		throw ModelException("lowBound (" + Format::toString(lowBound) +
			") >= highBound (" + Format::toString(highBound)+ ")");
	}
}

/** return grid points between which density is linear.
	For uniform distribution these are just points {start,end} */
DoubleArraySP UniformDensityPrior::getGridPoints(
	double start, 
	double end,
	const DateTime & date // ignored here
	) 
{
	static const string method ="UniformDensityPrior::getGridPoints";
	try
	{
		if(start > end) throw ModelException("Need start <= end");

		/** output array */
		DoubleArraySP out(new DoubleArray(2));
		
		// floor at 0 and cap at 1
		(*out)[0] = Maths::max(lowBound,Maths::min(start,highBound));
		(*out)[1] = Maths::max(lowBound,Maths::min(end,highBound));
		
		return out;
	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
}
/** return uniform (constant) density at gridPoints */
DoubleArraySP UniformDensityPrior::getDensity(const DoubleArray & gridPoints, 
											  const DateTime & date			// Not needed for uniform
											  ) const
{
	static const string method ="UniformDensityPrior::getDensity";
	try
	{
		
		double norm = highBound - lowBound;

		/** output density array */
		DoubleArraySP out(new DoubleArray(gridPoints.size()));

		// populate out array 
		for(int i = 0;i< out->size();i++)
		{
			// check gridPoints
			if(gridPoints[i] < lowBound || gridPoints[i] > highBound)
			{
				throw ModelException("gridPoint at index "+ 
					Format::toString(i)+ " (" + 
					Format::toString(gridPoints[i]) +") is not in range [" +
					Format::toString(lowBound) +", " +
					Format::toString(highBound) +"]");
			}
			
			
			(*out)[i] = 1.0/norm;			
		}

		return out;
	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}
}


/** set range on which density is defined */
void UniformDensityPrior::setRange(double lowB, double highB)
{

	static const string method ="UniformDensityPrior::setRange";
	try
	{
		lowBound = lowB;
		highBound = highB;

		if(lowBound >= highBound) 
		{
			throw ModelException("lowBound ("+ Format::toString(lowBound) +
				") >= highBound (" + Format::toString(highBound) + ")");
		}
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}


void UniformDensityPrior::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(UniformDensityPrior, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IDensityPrior);
    EMPTY_SHELL_METHOD(defaultUniformDensityPrior);

	FIELD(lowBound, "Lower bound of range on which density is defined [default = 0]");
	FIELD_MAKE_OPTIONAL(lowBound);
	
	FIELD(highBound, "Upper bound of range on which density is defined [default = 1]");
	FIELD_MAKE_OPTIONAL(highBound);

}

IObject* UniformDensityPrior::defaultUniformDensityPrior() {
    return new UniformDensityPrior();
}

CClassConstSP const UniformDensityPrior::TYPE = 
    CClass::registerClassLoadMethod("UniformDensityPrior", 
                                    typeid(UniformDensityPrior), 
                                    load);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool UniformDensityPriorLoad() {
    return (UniformDensityPrior::TYPE != 0);
}


DRLIB_END_NAMESPACE
