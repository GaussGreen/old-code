//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Maths.hpp
//
//   Description : Various useful maths type functions etc
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MATHS_HPP
#define EDG_MATHS_HPP

#include <float.h>
#include <math.h>
#include "edginc/Format.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/Complex.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

/** Various useful maths type functions etc */

class TOOLKIT_DLL Maths {
public:
    static const double PI;
    static const double ROOT_TWO_PI;

    /** replacement for MAX macro */
    static double max(double a, double b);

    /** replacement for MAX macro */
    static int max(int a, int b);

    /** replacement for MIN macro */
    static double min(double a, double b);

    /** replacement for MIN macro */
    static int min(int a, int b);

    static double collar(double amt, 
				  	     double cap, 
					     double floor);

    static int collar(int amt, 
				  	  int cap, 
					  int floor);

    //constrains the return value y to be 0 >= y >= u-l
    static double shiftedCollar(double x, double l, double u);

    //constrains the return value y to be 0 >= y >= u-l
    static int shiftedCollar(int x, int l, int u);


    /** Payoff function that generalises the usual collar function to
        negative strikes. k2>k1 whatever the signs. */
    static double creditCollar(double x,    // (I) argument
                               double k1,   // (I) lower strike
                               double k2);  // (I) upper strike


    /** rounds double to nearest integer. X.5 is rounded up */
    static int round(double x); 

    /** replacement for SIGN macro */
    inline static int sign(int a){
        return ((a >= 0)? 1: -1);
    }

    inline static double sign(double a){
        return ((a >= 0.0)? 1.0: -1.0);
    }

    inline static double sign(double a, double b){
        return ((b) >= 0.0 ? fabs(a) : -fabs(a));
    }

    /** replacement for SQUARE macro */
    static double square(double a);

    inline static int square(int a){
        return (a * a);
    }

    inline static Complex square(const Complex& a){
        return (a * a);
    }

    /** replacement for ARE_ALMOST_EQUAL macro. Returns true if 'a' and 'b'
        are within DBL_EPSILON of each other */
    static bool equals(double a, double b);

    /** replacement for IS_ALMOST_ZERO macro. Returns true if 'a'
        is within DBL_EPSILON of 0.0  */
    static bool isZero(double a);

    inline static bool isZero(const Complex& a){
        return Maths::isZero(Complex::absSquare(a));
    }

    /** Is it negative? nicer way of saying  x < - DBL_EPSILON */
    static bool isNegative(double a);

    /** Is it positive? nicer way of saying  x > DBL_EPSILON */
    static bool isPositive(double a);

    /** replacement for ARE_EQUAL_TO_TOL. Returns true if a = b within
        tolerance */
    inline static bool areEqualWithinTol(double a, double b, double tol){
        double diff = a - b;
        return ((diff <= tol) && (diff >= -tol));
    }

    /** Throws exception if check fails */
    inline static void checkNonNegative(double a, string name) {
        if (isNegative(a)){
            throw ModelException(name + " ("+
                                 Format::toString(a)+") is negative");
        }
    }

    inline static void checkNonNegative(int a, string name) {
        if (a<0){
            throw ModelException(name + " ("+
                                 Format::toString(a)+") is negative");
        }
    }
    
    inline static void checkNonNegative(const CDoubleArray& v, string name) {
        for (int i = 0; i < v.size(); ++i){
            if (isNegative(v[i])){
                throw ModelException(name + "[" + Format::toString(i+1) + "]"
                                     + " ("+ Format::toString(v[i])+") is negative");
            }
        }
    }
    
    /** Throws exception if check fails */
    inline static void checkPositive(double a, string name) {
        if (!isPositive(a)){
            throw ModelException(name + " ("+
                                 Format::toString(a)+") is non-positive");
        }
    }

    inline static void checkPositive(int a, string name) {
        if (a <= 0){
            throw ModelException(name + " ("+
                                 Format::toString(a)+") is non-positive");
        }
    }

    inline static void checkPositive(const CDoubleArray& v, string name) {
        for (int i = 0; i < v.size(); ++i){
            if (!isPositive(v[i])){
                throw ModelException(name + "[" + Format::toString(i+1) + "]"
                                     + " ("+ Format::toString(v[i])+") is non-positive");
            }
        }
    }

    ////Determines whether given double-precision floating point value is finite
    static bool finite(double x);

    class TOOLKIT_DLL SpecialFunc{
    public:
        static double gamma(double x);
    };

    /** calculates the product of 2 arrays of same sizes */
    static DoubleArray ArrayProduct(const DoubleArray& array1, const DoubleArray& array2);

    /** calculates the scalar product of 2 arrays of same sizes */
    static double ArrayScalarProduct(const DoubleArray& array1, const DoubleArray& array2);
    
    /** returns the maximum element of a double array 
        not very efficient algorithm */
    static double ArrayMax(const DoubleArray& array1);
};

#if !defined(DEBUG) || defined(QLIB_MATHS_CPP)
OPT_INLINE bool Maths::isZero(double a){
    return (a <= DBL_EPSILON && a >= -DBL_EPSILON);
}
/** replacement for MAX macro */
OPT_INLINE  double Maths::max(double a, double b){
    return ((a > b)? a: b);
}
/** replacement for MIN macro */
OPT_INLINE  double Maths::min(double a, double b){
    return ((a < b)? a: b);
}
OPT_INLINE bool Maths::isPositive(double a){
    return (a > DBL_EPSILON);
}
OPT_INLINE bool Maths::isNegative(double a){
    return (a < - DBL_EPSILON);
}
OPT_INLINE double Maths::square(double a){
    return (a * a);
}
OPT_INLINE int Maths::max(int a, int b){
    return ((a > b)? a: b);
}
OPT_INLINE int Maths::min(int a, int b){
    return ((a < b)? a: b);
}
OPT_INLINE bool Maths::equals(double a, double b){
    double diff = a - b;
    return (isZero(diff));
}
OPT_INLINE double Maths::collar(double amt, 
						  	    double cap, 
							    double floor){
    return (Maths::max(floor,Maths::min(amt,cap)));
}
OPT_INLINE int Maths::collar(int amt, 
						  	 int cap, 
							 int floor){
    return (Maths::max(floor,Maths::min(amt,cap)));
}

OPT_INLINE double Maths::shiftedCollar(double x, double l, double u) {
    return collar(x,u,l) - l;
}

//constrains the return value y to be 0 >= y >= u-l
OPT_INLINE int Maths::shiftedCollar(int x, int l, int u) {
    return collar(x,u,l) - l;
}

/** rounds double to nearest integer. X.5 is rounded up */
OPT_INLINE int Maths::round(double x) {
    return (int) ( x + 0.5 );
}

/** Payoff function that generalises the usual collar function to
    negative strikes. k2>k1 whatever the signs. */
OPT_INLINE double Maths::creditCollar(double x,    // (I) argument
                                      double k1,   // (I) lower strike
                                      double k2)   // (I) upper strike
{  
    return Maths::shiftedCollar(x, k1, k2) + 
           Maths::min(0.0, k1) -
           Maths::min(0.0, k2);
}

#endif

DRLIB_END_NAMESPACE

#endif

