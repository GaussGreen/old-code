//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : GridFactory.hpp
//
//   Description : create various grids useful for for integration 
//
//   Author      : Matthias Arnsdorf
//
//   Date        : October 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/AtomicArray.hpp"



#ifndef QLIB_GRID_FACTORY
#define QLIB_GRID_FACTORY


DRLIB_BEGIN_NAMESPACE


class UTIL_DLL GridFactory {
public:
	
	/** exponential grid */
	static DoubleArraySP exponential(double start, double end, int numPoints);
	/** homogeneous grid */
	static DoubleArraySP homogeneous(double start, double end, int numPoints);

	/** grid of type exp(-x^2) from -\infty to 0*/
	static DoubleArraySP semiGaussian(double start, double end, int numPoints);

	/** merges two strictly increasing doubleArrays and removes duplicates */
	static DoubleArraySP mergeNoDuplicate(const DoubleArray & array1, const DoubleArray & array2);

private:

	static const double TINY;

};

DRLIB_END_NAMESPACE

#endif
