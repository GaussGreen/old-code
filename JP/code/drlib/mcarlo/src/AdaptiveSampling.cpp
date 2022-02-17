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
#include "edginc/AdaptiveSampling.hpp"
#include "edginc/imsl.h"

DRLIB_BEGIN_NAMESPACE



IObject* AdaptiveSampling::defaultAdaptiveSampling()
{
	return new AdaptiveSampling();
};

CClassConstSP const AdaptiveSampling::TYPE =
	CClass::registerClassLoadMethod(
		"AdaptiveSampling", 
		typeid(AdaptiveSampling), 
		load);

void
AdaptiveSampling::load(CClassSP& clazz)
{
    clazz->setPublic();
    
	REGISTER(AdaptiveSampling, clazz);

	SUPERCLASS(CObject);

	IMPLEMENTS(IMFSampling);

	FIELD(minX, "lower bound for sampling, default is -7.5");
	FIELD_MAKE_OPTIONAL(minX);

	FIELD(maxX, "upper bound for sampling, default is 7.5");
	FIELD_MAKE_OPTIONAL(maxX);

	FIELD(numberInitialPoints, "Number of uniform points in first sample, "
							   "default is 100");
	FIELD_MAKE_OPTIONAL(numberInitialPoints);
    
	FIELD(samplesPerPoint, "Number of samples per point in first sample, "
						   "default is 5");
	FIELD_MAKE_OPTIONAL(samplesPerPoint);

	FIELD(stretchFactor, "Factor for min and max values of Gaussian scheme, "
						 "default is 0.25");
	FIELD_MAKE_OPTIONAL(stretchFactor);

	FIELD(outsideRatio, "% of non-Gaussian points, default is 0.2");
	FIELD_MAKE_OPTIONAL(outsideRatio);

	FIELD(lowGaussian, "Low end of Gaussian sampling");
	FIELD_MAKE_TRANSIENT(lowGaussian); 
	
	FIELD(hiGaussian, "Hi end of Gaussian sampling");
	FIELD_MAKE_TRANSIENT(hiGaussian);  
	
	FIELD(lowEnd, "Actual low end of interval sampled");
	FIELD_MAKE_TRANSIENT(lowEnd);     

	FIELD(highEnd, "Actual high end of interval sampled");
	FIELD_MAKE_TRANSIENT(highEnd);     

	FIELD(allPoints, "Set of sampled points");
	FIELD_MAKE_TRANSIENT(allPoints);  

	FIELD(cumSamplesPerPoint, "Cum number of samples per point");
	FIELD_MAKE_TRANSIENT(cumSamplesPerPoint); 

	EMPTY_SHELL_METHOD(defaultAdaptiveSampling);
}

void 
AdaptiveSampling::validatePop2Object()
{
// Here we check that numPoints > 0, minX < maxX (== if numPoints == 1)

// when we know how many samples there are to compute, we can find which value 
// of bias gives us that many samples considering the roundoff
	
};

IMFSamplingSP AdaptiveSampling::getInitialSampling() const
{
	DoubleArraySP points = DoubleArraySP(new DoubleArray(numberInitialPoints));

	IntArraySP samples = IntArraySP(new IntArray(numberInitialPoints));

	double dx = (maxX-minX)/(numberInitialPoints-1);

	int i;
	for (i = 0; i < numberInitialPoints; i++)
	{
		(*points)[i] = minX+i*dx;
		(*samples)[i] = samplesPerPoint;
	};

	return IMFSamplingSP(new ArbSampling(points, samples, false));
	
};


bool AdaptiveSamplingLoad()
{
	return AdaptiveSampling::TYPE != NULL;
};

void 	
AdaptiveSampling::config(
		int					nSamples,
		DoubleArrayConstSP	points,
		IntArrayConstSP		numNonTrivial)
{
	int nPoints = points->size();

	int first = 0;
	while ((*numNonTrivial)[first] == 0)
	{
		first ++;
		if (first == numNonTrivial->size())
		{
			first = 0;
			break;
		}
	};

	int last = nPoints-1;
	while ((*numNonTrivial)[last] == 0)
	{
		last --;
		if (last == -1)
		{
			last = nPoints-1;
			break;
		}
	};
	
	lowGaussian = (*points)[first] 
					- stretchFactor * ((*points)[last] - (*points)[first]);
	hiGaussian  = (*points)[last] 
					+ stretchFactor * ((*points)[last] - (*points)[first]);
	
	// we floor the gaussian start at the specified min
	lowGaussian = max(minX, lowGaussian);
	lowEnd = min(minX, lowGaussian);

	// and cap the gaussian end at the specified max
	hiGaussian = min(maxX, hiGaussian);
	highEnd = max(maxX, hiGaussian);
	
	int lowPoints, hiPoints;

	if ( (lowEnd == lowGaussian) && (highEnd == hiGaussian))
	{
		lowPoints = hiPoints = 0;
	} else {

		lowPoints = nSamples*outsideRatio*(lowGaussian -lowEnd) 
						/ ( (lowGaussian -lowEnd) + (highEnd - hiGaussian) );

		hiPoints = nSamples*outsideRatio*(lowGaussian -lowEnd) 
						/ ( (lowGaussian -lowEnd) + (highEnd - hiGaussian) );
	}

	int gaussPoints = nSamples - lowPoints - hiPoints;

	allPoints = DoubleArraySP(new DoubleArray(nPoints));
	cumSamplesPerPoint  = IntArraySP(new IntArray(nPoints));

	int i;
	int nLow = 0;
	int nHigh = 0;

	double dx = (highEnd - lowEnd)/(nPoints-1);

	for (i = 0; i < nPoints; i++)
	{
		(*allPoints)[i] = lowEnd + i*dx;

		if ((*allPoints)[i] < lowGaussian) {
			nLow++;
		};

		if ((*allPoints)[i] > hiGaussian) {
			nHigh++;
		}
	}

	int nGaussian = nPoints - nLow - nHigh;

	for (i = 0; i < nLow; i++)
	{
		(*cumSamplesPerPoint)[i] = (1.*lowPoints)/(nLow-i);
		
		lowPoints -= (*cumSamplesPerPoint)[i];

		if (i)
			(*cumSamplesPerPoint)[i] += (*cumSamplesPerPoint)[i-1];
	};

	double norm = imsl_d_normal_cdf((*allPoints)[nLow+nGaussian-1]);

	for (i = nLow; i < nPoints-nHigh; i ++)
	{
		double total = floor(gaussPoints*imsl_d_normal_cdf( (*allPoints)[i])) 
						/ norm;

		if (nLow>0)
		{
			(*cumSamplesPerPoint)[i] = total + (*cumSamplesPerPoint)[nLow-1];
		} else
			(*cumSamplesPerPoint)[i] = total;
	}

	if (nLow>0)
		(*cumSamplesPerPoint)[nPoints-nHigh-1] = gaussPoints + (*cumSamplesPerPoint)[nLow-1];
	else
		(*cumSamplesPerPoint)[nPoints-nHigh-1] = gaussPoints;

	for (i = nPoints-nHigh; i < nPoints; i ++)
	{
		(*cumSamplesPerPoint)[i] = (1.*hiPoints)/(nPoints-i);
		
		hiPoints -= (*cumSamplesPerPoint)[i];

		if (i)
			(*cumSamplesPerPoint)[i] += (*cumSamplesPerPoint)[i-1];
	}
};


DRLIB_END_NAMESPACE