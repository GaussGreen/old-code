//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : UtilFuncs.hpp
//
//   Description : some useful simple ustilities
//
//   Author      : Ning Shen
//
//
//----------------------------------------------------------------------------
#ifndef EDG_UTILFUNCS_HPP
#define EDG_UTILFUNCS_HPP

#include <math.h>
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/**------------------------------------------------------------------------------
*   Description  :	Search for the neighbouring items in an array.
*
*   Parameters   :	target		target value for the search
*					valueList   Array to search, must be in ascending or descending order
*					start		Index of first element of array to search
*					end	      	Index of last element to search
*					mode		Search mode
*								mode = -1,  position of the largest item smaller than or equal to target
*								mode = 0, position where target is found in the list
*								mode = 1,  position of the smallest item larger than or equal to target
*
*   Returns     :  	>=0 	position of the found item.
*					<0		errors
*------------------------------------------------------------------------------*/
UTIL_DLL int Neighbour(double target, const double* valueList, int start, int end,  int mode);
#if defined(__GNUC__) && ( __GNUC__ >= 3)
    // added to allow input array - where iterators aren't typedefd to pointers
int Neighbour(double target, vector<double>::const_iterator valueList,
              int start, int end,  int mode);
#endif
UTIL_DLL int Neighbour(const DateTime& target, const DateTimeArray& valueList, int start, int end,  int mode);

/** Smax smoothing method */
UTIL_DLL double SMax(double x, double y, double smoothOver);

/*****  Smooth_Step  ********************************************************/
/*
*       Julia's smooth step function (cf. corresponding memo).
*/
UTIL_DLL double SmoothValue (
                        double	UpValue,       /* (I) Up value              */
                        double  UpLevel,
                        double	DownValue,     /* (I) Down value            */
                        double  DownLevel,
                        double	SmoothLevel);   /* (I) Index value           */

/** Impoved linear interpolation with 4 points. size must be 4.
return: <0 failure; 0 - no improvement, 1 - improved linear interp */
inline int LinearInterpImprove(const double* x, const double* y, int size, double s, double& result, bool payoffPositive)
{   // the result, if successful, is bound within triangle formed by
    // (x[1], y[1]), (x[2], y[2]) and intersection of lines 1->2 and 3->4.
    const double minResolution = 1.0e-16;
    int status = 0;

    if (size != 4)
        status = -1; // failure, only works for 4 pts
    else if (s < x[1] || x[2] < s)
        status = -2; // failure, must be within 2 inner pts
    else if (fabs(x[1] - x[2]) < minResolution) 
        result = (y[1]+y[2])/2.0; // just average
    else if (fabs(x[1] - x[0]) < minResolution || fabs(x[2] - x[3]) < minResolution)
        result = y[1]+ (y[2]-y[1])*(s-x[1])/(x[2]-x[1]); // use simple linear here
    else
    {
        double  a0 = (y[1] - y[0])/(x[1] - x[0]);
        double  a1 = (y[2] - y[1])/(x[2] - x[1]);
        double  a2 = (y[3] - y[2])/(x[3] - x[2]);
        result = a1*(s - x[1]) + y[1]; // linear interp
        if ((a0>a1 && a1>a2) || (a0<a1 && a1<a2))
        {
            double cross = (a0*x[0] - a2*x[2] + y[2] - y[0])/(a0-a2);
            double tmp;
            if (s < cross)
                tmp = (result + a0*(s - x[0]) + y[0])/2.0;
            else
                tmp = (result + a2*(s - x[2]) + y[2])/2.0;

            if (payoffPositive && tmp < 0.0)
                status = 1; // not improved but return result with no change.
            else
            {
                result = tmp;
                status = 1; // improved
            }
        }
    }
    return status;
}


/** returns max[cp*(s-k), 0] for option payoff or cp*(s-k) for fwd payoff */
inline double GetIntrinsic(double s, double k, bool isCall, bool isCallOrPut)
{
    double result;

    result = (isCall? s-k : k-s);
    if (isCallOrPut && result <0.0)
        result = 0.0;
    
    return result;
}

DRLIB_END_NAMESPACE
#endif
