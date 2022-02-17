//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : GridFactory.cpp
//
//   Description : create various Grids useful for for integration 
//
//   Author      : Matthias Arnsdorf
//
//   Date        : October 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GridFactory.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE


const double GridFactory::TINY = 1E-10;

/*************************************************************************************
  creates exponentially spaced  Grid of doubles 
  map uniform [0,1] onto [start, end] via y[i] = start + (exp(i*dx*b)-1)
  b = ln(end-start+1), dx = 1/(numPoints-1)
  y[0] = start, y[1] = end
	numPoints includes start and end points
*************************************************************************************/
DoubleArraySP GridFactory::exponential(double start,
											double end, 
											int numPoints
											) 
{
	static const string method = "GridFactory::exponential";
	try{

		if(numPoints<2) throw ModelException("numPoints < 2");
		if(start >= end) throw ModelException("start >= end");

		double dx = 1.0/((double) numPoints-1.0);
		
		double b  = log(end-start+1);

		DoubleArraySP gridPoints(new DoubleArray(numPoints));

		int i; // counter
		(*gridPoints)[0] =start; 
		for(i = 1; i< numPoints-1;i++) (*gridPoints)[i] = start+exp(b*dx*i)-1;
		(*gridPoints)[numPoints-1] = end; // to avoid numerical issues and guarantee last element is end

		return gridPoints;

	} catch(exception& e){
		throw ModelException(e,method);
	}
	

}

DoubleArraySP GridFactory::homogeneous(double start, double end, int numPoints)
{
	static const string method = "GridFactory::homogeneous";
	try
	{

		if(numPoints<2) throw ModelException("numPoints < 2");
		if(start >= end) throw ModelException("start >= end");

		double dx = (end-start)/((double) numPoints-1.0);

		DoubleArraySP gridPoints(new DoubleArray(numPoints));

		int i; // counter
		(*gridPoints)[0] =start; 
		for(i = 1; i< numPoints-1;i++) (*gridPoints)[i] = start+i*dx;
		(*gridPoints)[numPoints-1] = end; // to avoid numerical issues and guarantee last element is end

		return gridPoints;

	} catch(exception& e){
		throw ModelException(e,method);
	}
}


DoubleArraySP GridFactory::semiGaussian(double start, double end, int numPoints)
{
	static const string method = "GridFactory::semiGaussian";
	try
	{

		if(numPoints<2) throw ModelException("numPoints < 2");
		if(start >= end) throw ModelException("start >= end");

		double A = (end-start-TINY);

		double x_start = -sqrt(-log(TINY/A)); // -\infty
		if(x_start >=0) throw ModelException("x_start should be negative");

		double dx = -x_start/((double) numPoints-1.0);

		DoubleArraySP gridPoints(new DoubleArray(numPoints));

		int i; // counter
		(*gridPoints)[0] =start; 
		for(i = 1; i< numPoints-1;i++) (*gridPoints)[i] = start-TINY + A*exp(-Maths::square(x_start+ i*dx));
		(*gridPoints)[numPoints-1] = end; // to avoid numerical issues and guarantee last element is end

		return gridPoints;

	} catch(exception& e){
		throw ModelException(e,method);
	}
}


/***********************************************************************

  merges two INCREASING arrays and removes duplicates

  ********************************************************************/
DoubleArraySP GridFactory::mergeNoDuplicate(const DoubleArray & array1, const DoubleArray & array2)
{
	static const string method = "GridFactory::mergeNoDuplicate";
	try
	{
		int i;
		DoubleArraySP out = DoubleArraySP(new DoubleArray(0));

		/** index for array2 */
		int index2 = 0;
		
		// loop through array1
		for(i = 0; i< array1.size();i++)
		{	
			if(i>0 && array1[i] < array1[i-1]) throw ModelException("array1 is not increasing");
			// loop through array2
			while(index2 < array2.size() && array2[index2] <= array1[i])
			{
				if(index2>0 && array2[index2] < array2[index2-1]) throw ModelException("array2 is not increasing");
				// add  elements of array2 to out if not duplicate
				if(out->size()<1 || array2[index2] != out->back())
				{
					out->push_back(array2[index2]);
				}
				index2++;
			}
			// check if array1[i] has been included. If not need to add to out
			if(out->size()<1 || array1[i] != out->back())
			{
				out->push_back(array1[i]);
			}
		}

		// loop through any remainder of array2
		for(i = index2; i < array2.size();i++)
		{
			if(i>0 && array2[i] < array2[i-1]) throw ModelException("array2 is not increasing");
				// add  elements of array2 to out if not duplicate
			// add  elements of array2 to out if not duplicate
			if(out->size()<1 || array2[i] != out->back())
			{
				out->push_back(array2[i]);
			}

		}
		return out;

	} catch(exception& e)
	{
		throw ModelException(e,method);
	}
}

DRLIB_END_NAMESPACE
