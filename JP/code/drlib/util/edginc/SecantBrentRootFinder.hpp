//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : SecantBrentRootFinder.hpp
//
//   Description : Port of ALIB root finder.
//
//   Author      : Richard Appleton
//
//   Date        : 10th June 2005
//
//----------------------------------------------------------------------------

#ifndef SECANT_BRENT_ROOT_FINDER_HPP
#define SECANT_BRENT_ROOT_FINDER_HPP

#include "edginc/config.hpp"
#include "edginc/RootFinder.hpp"


DRLIB_BEGIN_NAMESPACE


/*
 * ALIB converges in some instances where ZBrent does not. Therefore ALIB
 * root finding code has been ported as well.
 */
class UTIL_DLL SecantBrentRootFinder
{
public:
    SecantBrentRootFinder(double tolerance = 1E-10, int maxIterations = 50);

    double solve(
        const Func1D::NoDeriv& fn, 
        const double           guess, 
        const double           lowerBound, 
        const double           upperBound,
        double                 initialXStep = 0.0001);

private:
    double secantMethod(
        const Func1D::NoDeriv& fn, 
		double*                xPoints,        // (I/O) Array of x points
		double*                yPoints,        // (I/O) Array of y points
		double                 boundLo,        // (I) Lower bound
		double                 boundHi,        // (I) Upper bound
		bool&                  foundIt,        // (O) If solution was found
		bool&                  bracketed);     // (O) if root was bracketed

	double brentMethod(
        const Func1D::NoDeriv& fn, 
        double                 *xPoints,       // (I) Array of x values
        double                 *yPoints);      // (I) Array of y values

    double* xPoints;             // array of x values
    double* yPoints;             // array of y values
    double  tolerance;           // y accuracy tolerance
    int     maxIterations;
};




DRLIB_END_NAMESPACE

#endif // SECANT_BRENT_ROOT_FINDER_HPP
