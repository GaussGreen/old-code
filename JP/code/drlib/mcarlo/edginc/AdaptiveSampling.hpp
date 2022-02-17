//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Filename    : GaussSampling.hpp
//
// Description : 
//
// Date        : Oct 2006
//
//----------------------------------------------------------------------------

#ifndef ADAPTIVESAMPLING_HPP
#define ADAPTIVESAMPLING_HPP

#include "edginc/ArbSampling.hpp"
#include "edginc/Array.hpp"
#include "edginc/IMFSampling.hpp"

DRLIB_BEGIN_NAMESPACE

class StratMCRandomNoCache;
// Gaussian strat sampling
class MCARLO_DLL AdaptiveSampling : public CObject,
										public virtual IMFSampling
{
	// first 4 parameters define the initial grid (uniform) and how it's 
	// sampled. 

	double minX;
	
	double maxX;

	int numberInitialPoints;

	int samplesPerPoint;

	double stretchFactor; // normally 25%
	double outsideRatio;  // normally 20%

	// we know these once the scheme is configured
	double lowGaussian;
	double hiGaussian;
	double lowEnd;
	double highEnd;

	DoubleArraySP allPoints;
	IntArraySP    cumSamplesPerPoint;


public:

	friend class StratMCRandomNoCache;

	static CClassConstSP const TYPE;

	AdaptiveSampling(
		int nPoints = 100,
		double minX = -7.5,
		double maxX = 7.5,
		int samplesPerPoint = 5,
		double stretchFactor = 0.25,
		double outsideRatio = 0.2,
		CClassConstSP clazz = TYPE):
	CObject(clazz),
	numberInitialPoints(nPoints),
	minX(minX),
	maxX(maxX),
	samplesPerPoint(samplesPerPoint),
	stretchFactor(stretchFactor),
	outsideRatio(outsideRatio){};

	virtual ~AdaptiveSampling() {};

	virtual void validatePop2Object(); 

	static void load(CClassSP& clazz);

	static IObject* defaultAdaptiveSampling();

	// the adaptive sampling takes in the results of a preview simulation :
	// set of market values and the number of interesting results for each
	// for now, we are only using the highest and lowest market values that 
	// yield interesting results
	void config(
		int				    nSamples,
		DoubleArrayConstSP	points,
		IntArrayConstSP		numNonTrivial);

	IMFSamplingSP getInitialSampling() const;

	virtual DoubleArrayConstSP getPoints() const {return allPoints;};

	virtual IntArrayConstSP getCumNumberSamples() const {return cumSamplesPerPoint;};
	
	virtual int getNumSamples() const {
		return (*cumSamplesPerPoint)[cumSamplesPerPoint->size()-1];
	};

};

DECLARE(AdaptiveSampling);

DRLIB_END_NAMESPACE

#endif