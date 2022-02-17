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

#ifndef GAUSSSAMPLING_HPP
#define GAUSSSAMPLING_HPP

// #include "edginc/ArbSampling.hpp"  JCP
#include "edginc/Array.hpp"
#include "edginc/IMFSampling.hpp"

DRLIB_BEGIN_NAMESPACE


// Gaussian strat sampling
class MCARLO_DLL GaussSampling : public CObject,
							    public virtual IMFSampling
{
	int numPoints;
	int numSamples;
	double minX;
	double maxX;
	double bias;

	DoubleArraySP allPoints;
	IntArraySP    cumSamplesPerPoint;

public:

	friend class StratMCRandomNoCache;

	static CClassConstSP const TYPE;

	GaussSampling(
		int nPoints = 100,
		double minX = -7.5,
		double maxX = 7.5,
		double bias = 0.0,
		CClassConstSP clazz = TYPE):
	CObject(clazz),
	numPoints(nPoints),
	minX(minX),
	maxX(maxX),
	bias(bias) {};

	virtual ~GaussSampling() {};

	virtual void validatePop2Object(); 

	static void load(CClassSP& clazz);

	static IObject* defaultGaussSampling();

	void config(int nSamples);

	virtual DoubleArrayConstSP getPoints() const {return allPoints;};

	virtual IntArrayConstSP getCumNumberSamples() const {return cumSamplesPerPoint;};

	virtual int getNumSamples() const {
		return (*cumSamplesPerPoint)[cumSamplesPerPoint->size()-1];
	};

};

DECLARE(GaussSampling);

DRLIB_END_NAMESPACE

#endif