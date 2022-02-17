//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Maths.cpp
//
//   Description : Various useful maths type functions etc
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_MATHS_CPP
#include "edginc/Maths.hpp"
#include "edginc/imslerror.hpp"
#if defined(sun) || defined(__CYGWIN__)
// for finite() function - but not on Linux
#include <ieeefp.h>
#endif
#include <math.h>

DRLIB_BEGIN_NAMESPACE

const double Maths::PI          = 3.14159265358979;
const double Maths::ROOT_TWO_PI = 2.506628274631;

double Maths::SpecialFunc::gamma(double x){
    try{
        double rtn = imsl_d_gamma(x);
        IMSLError::throwExceptionIfError();
        return rtn;
    }
    catch(exception& e){
        throw ModelException::addTextToException(e,
                                                 "Maths::SpecialFunc::gamma: "
                                                 "Failed with input value "
                                                 + Format::toString(x));
    }
}

/** returns the product of 2 arrays of same size */
DoubleArray Maths::ArrayProduct(const DoubleArray& array1,
                                const DoubleArray& array2) {
    static const string method = "Maths::ArrayProduct";
    try {
        // validation: same number of elements and not empty arrays
        if (array1.size() != array2.size()) {
            throw ModelException(method, "the first array (size " + Format::toString(array1.size()) + 
                                 ") and the second array (size " + Format::toString(array2.size()) +
                                 ") must have the same size");
        }
        if (array1.size() == 0) {
            throw ModelException(method, "the two arrays must contain some elements whereas they are empty");
        }
       
        // result
        DoubleArray result(array1.size());
        
        for(int i=0; i<array1.size(); i++){
            result[i] = array1[i] * array2[i];
        }
        
        return result;
    }
    catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** returns the scalar product of 2 arrays of same size */
double Maths::ArrayScalarProduct(const DoubleArray& array1,
                                const DoubleArray& array2) {
    static const string method = "Maths::ArrayScalarProduct";
    try {
        DoubleArray product = Maths::ArrayProduct(array1,array2);
       
        // result
        double result(0.0);
        
        for(int i=0; i<product.size(); i++){
            result += product[i];
        }
        
        return result;
    }
    catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** returns the maximum of an non-empty array */
double Maths::ArrayMax(const DoubleArray& array1) {
    static const string method = "Maths::ArrayMax";
    try {
        // validation: not empty array
        if (array1.size() == 0) {
            throw ModelException(method, "the array must contain some elements whereas it is empty");
        }
       
        // result 
        double result = array1[0];
        
        for(int i=1; i<array1.size(); i++){
            result = Maths::max(array1[i], result);
        }
        
        return result;
    }
    catch(exception& e) {
        throw ModelException(e, method);
    }
}


////Determines whether given double-precision floating point value is finite
bool Maths::finite(double x){
#if defined(UNIX) || defined(__CYGWIN__)
    return ::finite(x) != 0;
#elif defined(WIN32)
    return _finite(x) != 0;
#else
#error finite not defined on this platform
#endif
}

DRLIB_END_NAMESPACE
