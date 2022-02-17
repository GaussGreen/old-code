//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FindRoot.hpp - from CMLib FindRoot.h
//
//   Description : Base stuff for finding roots
//
//   Author      : CMLib
//
//   Date        : Dec 17, 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FINDROOT_HPP
#define EDR_FINDROOT_HPP

#include "float.h"

// in debug mode, keep track of points to display 
// why solver did not converge
#ifdef _DEBUG
#define EDR_STORE_SOLVER_POINTS
#endif

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL FindRoot{
public:
    static const double LOWER_BOUND; // = -DBL_MAX;
    static const double UPPER_BOUND; // = DBL_MAX;
    
    static const double X_TOLERANCE; //  = 1e-08;
    static const double F_TOLERANCE; //  = 1e-12;

    static const int    MAX_NUM_OF_ITERATIONS; //  = 50;
    static const int    FIXED_NUM_OF_ITERATIONS; //  = 20;

    static const double INITIAL_X_STEP_SIZE; //  = 1E-8;
    static const double INITIAL_F_DERIVATIVE; //  = 0;

    /** This is what the root solvers return */
    struct Point {
        double x;
        double y;
        
        Point();
        Point(double x, double y);
    };

private:
    FindRoot(); // do not create instances of this class
};

DRLIB_END_NAMESPACE

#endif
