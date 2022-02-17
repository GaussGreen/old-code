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

#ifndef IMFSAMPLING_HPP
#define IMFSAMPLING_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class IMFSampling : public virtual IObject {

public:

	static CClassConstSP const TYPE;

	virtual DoubleArrayConstSP getPoints() const = 0;

	virtual IntArrayConstSP getCumNumberSamples() const = 0;

	virtual int getNumSamples() const = 0;

};

DECLARE(IMFSampling);

DRLIB_END_NAMESPACE

#endif