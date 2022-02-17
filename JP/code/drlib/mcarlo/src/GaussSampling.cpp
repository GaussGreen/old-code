//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids DR
//
//   Filename    : AdaptiveSampling.cpp
//
//   Description : 
//
//   Date        : Oct 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GaussSampling.hpp"
#include "edginc/imsl.h"

DRLIB_BEGIN_NAMESPACE

IObject* GaussSampling::defaultGaussSampling()
{
	return new GaussSampling();
};

CClassConstSP const GaussSampling::TYPE =
	CClass::registerClassLoadMethod(
		"GaussSampling", 
		typeid(GaussSampling), 
		load);

void
GaussSampling::load(CClassSP& clazz)
{
    clazz->setPublic();
    
	REGISTER(GaussSampling, clazz);

	SUPERCLASS(CObject);

	IMPLEMENTS(IMFSampling);
	
	FIELD(numPoints, "Number of grid points for market factor");

	FIELD(minX, "Low end of sampling grid");

	FIELD(maxX, "High end of sampling grid");

	FIELD(bias, "Shift in the probability density function");

	EMPTY_SHELL_METHOD(defaultGaussSampling);
}

void 
GaussSampling::validatePop2Object()
{
// Here we check that numPoints > 0, minX < maxX (== if numPoints == 1)

// when we know how many samples there are to compute, we can find which value 
// of bias gives us that many samples considering the roundoff
	
};

bool GaussSamplingLoad()
{
	return GaussSampling::TYPE != NULL;
};

void GaussSampling::config(int nSamples)
{
	int i, total = 0;

	allPoints = DoubleArraySP(new DoubleArray(numPoints));
	cumSamplesPerPoint = IntArraySP(new IntArray(numPoints));

	double dx = (numPoints>1) ? 
					(maxX-minX)/(numPoints-1) 
					: 0;

	double norm = imsl_d_normal_cdf(maxX+bias);

	for (i = 0; i < numPoints; i++)
	{
		(*allPoints)[i] = minX + i*dx;

		total =
			floor(nSamples*imsl_d_normal_cdf((*allPoints)[i]+bias)/norm);

		// store the cummulative number of paths in numSamples
		(*cumSamplesPerPoint)[i] = total;
	}
	(*cumSamplesPerPoint)[numPoints-1] = nSamples;
};

DRLIB_END_NAMESPACE