//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Filename    : IMFSampling.hpp
//
// Description : 
//
// Date        : Oct 2006
//
//----------------------------------------------------------------------------

#ifndef ARBSAMPLING_HPP
#define ARBSAMPLING_HPP

#include "edginc/ArbSampling.hpp"
#include "edginc/Array.hpp"
#include "edginc/IMFSampling.hpp"

DRLIB_BEGIN_NAMESPACE

class MCARLO_DLL ArbSampling : public CObject,
	public virtual IMFSampling
{

	DoubleArraySP points;
	IntArraySP    samples;

public:
	
	static CClassConstSP const TYPE;

	ArbSampling(CClassConstSP clazz = TYPE):
	CObject(clazz) {};

	ArbSampling(
		DoubleArrayConstSP inPoints,
		IntArrayConstSP    inSamples,
		bool               isCum,
		CClassConstSP clazz = TYPE):
		CObject(clazz),
		points(DoubleArraySP(new DoubleArray(*(inPoints.get())))),
		samples(IntArraySP(new IntArray(*(inSamples.get()))))
		{
			if (isCum) return;

			int i;
			for (i = 1; i < samples->size(); i++)
			{
				(*samples)[i] += (*samples)[i-1];
			}
		
		};

	virtual void validatePop2Object(); 

	static void load(CClassSP& clazz);

	static IObject* defaultArbSampling();

	virtual DoubleArrayConstSP getPoints() const {return points;};

	virtual IntArrayConstSP getCumNumberSamples() const {return samples;};

	virtual int getNumSamples() const {return (*samples)[samples->size()-1];};

};

DECLARE(ArbSampling);

DRLIB_END_NAMESPACE

#endif