//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : UtilFuncs.cpp
//
//   Description : some useful simple ustilities
//
//   Author      : Ning Shen
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"

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
#if defined(__GNUC__) && ( __GNUC__ >= 3)
    // added to allow input array - where iterators aren't typedefd to pointers
int Neighbour(double target, vector<double>::const_iterator valueList,
              int start, int end,  int mode){
    return Neighbour(target, &*valueList, start, end, mode);
}
#endif

int Neighbour(double target, const double* valueList, int start, int end,  int mode)
{
    int mid, low, high;
    bool ascend;	// true=ascending order, false=descending order
    int result = -9;

    if (end <start)
        //	throw ModelException("Locate", " end < start.");
        return -1;
    if (valueList == 0)
        // throw ModelException("Locate", " invalid valueList for searching.");
        return -2;

    ascend = (valueList[end] > valueList[start]);
    low = start-1;
    high = end+1;

    while (high - low > 1)
    {
        mid = ((high + low) >> 1);
        if ((target > valueList[mid]) == ascend)
            low = mid;
        else
            high = mid;
    }

    if (low == start-1 || high == end+1)	// only one point can possibly be used
    {
        mid = (low == start-1 ? start : end);
        if ((mode == 0 && target == valueList[mid])
            ||(mode == -1 && target >= valueList[mid])
            ||(mode == 1 && target <= valueList[mid]))
            result = mid;
        return result;
    }

    switch (mode)
    {
    case 0:
        if (target == valueList[low])
            result = low;
        else if (target == valueList[high])
            result = high;
        break;
    case -1:
        result = (ascend ? low : high);
        break;
    case 1:
        result = (ascend ? high : low);
        break;
    default:
        result = -3;
        // throw ModelException("Locate", " invalid mode request.");
    }
    return result;
}

/**------------------------------------------------------------------------------
*   Description  :	Search for the neighbouring items in an array.
*
*   Same as above, but for TimeDate Array.
*------------------------------------------------------------------------------*/
int Neighbour(const DateTime& target, const DateTimeArray& valueList, int start, int end,  int mode)
{
    int mid, low, high;
    bool ascend;	// true=ascending order, false=descending order
    int result = -9;

    if (end <start)
        //	throw ModelException("Locate", " end < start.");
        return -1;
    if (valueList.size()<1)
        // throw ModelException("Locate", " invalid valueList for searching.");
        return -2;

    ascend = (valueList[end] > valueList[start]);
    low = start-1;
    high = end+1;

    while (high - low > 1)
    {
        mid = ((high + low) >> 1);
        if ((target > valueList[mid]) == ascend)
            low = mid;
        else
            high = mid;
    }

    if (low == start-1 || high == end+1)	// only one point can possibly be used
    {
        mid = (low == start-1 ? start : end);
        if ((mode == 0 && target == valueList[mid])
            ||(mode == -1 && target >= valueList[mid])
            ||(mode == 1 && target <= valueList[mid]))
            result = mid;
        return result;
    }

    switch (mode)
    {
    case 0:
        if (target == valueList[low])
            result = low;
        else if (target == valueList[high])
            result = high;
        break;
    case -1:
        if (target == valueList[low])
            result = low;
        else if (target == valueList[high])
            result = high;
        else
            result = (ascend ? low : high);
        break;
    case 1:
        if (target == valueList[high])
            result = high;
        else if (target == valueList[low])
            result = low;
        else
            result = (ascend ? high : low);
        break;
    default:
        result = -3;
        // throw ModelException("Locate", " invalid mode request.");
    }
    return result;
}

/****************************************************************
 * DRSAbs function. Smooth Absolute Value
 ****************************************************************/
double SmoothAbs(double  x)                       /* (I) */
{
    double  smoothAbs    = 0;

    if (x < 0)
        x = -x;
    
    if (x < 0.5)
    {
        smoothAbs  = - x * x + 3;
        smoothAbs *=   x * x;
        smoothAbs +=   0.875;
        smoothAbs /=     3;
    }
    else if (x < 1)
    {
        smoothAbs  = x - 4;
        smoothAbs *= x;
        smoothAbs += 6;
        smoothAbs *= x;
        smoothAbs -= 1;
        smoothAbs *= x;
        smoothAbs += 1;
        smoothAbs /= 3;
    }
    else    /* x >= 1*/
    {
        smoothAbs = x;
    }
    return (smoothAbs);
}
/****************************************************************
 * DRSMax function. Smooth Maximum Funtion
 ****************************************************************/
double SMax(double     x,                       /* (I) */
            double     y,                       /* (I) */
            double     smoothOver)              /* (I) */
{
    const double minDiff = 1e-10;
    double  smax    = 0;
    if (fabs(smoothOver) < minDiff)
        smax = (x>y? x:y);
    else
    {
        smax  = x + y + smoothOver * SmoothAbs((x - y)/smoothOver);
        smax /= 2.;
    }

    if (smax < 0.0)
        smax = 0.0; // protecting from rounding
    return (smax);
}

/*****  Smooth_Step  ********************************************************/
/*
*       Julia's smooth step function (cf. corresponding memo).
*/
double SmoothValue (
                        double	UpValue,       /* (I) Up value              */
                        double  UpLevel,
                        double	DownValue,     /* (I) Down value            */
                        double  DownLevel,
                        double	SmoothLevel)   /* (I) Index value           */
{

    double  x, y;
    double result;          
    double range;

    range = (UpLevel - DownLevel) / 2.0;
                     
    x = SmoothLevel - (UpLevel + DownLevel) / 2.0;
    if (range==0.) {
        if (x < 0.) { 
            return DownValue;
        }
        return UpValue;
    } 

    x /= range;
    
    if (x < -1.)                    /* Step function = 0 below step */
    { 
        result = DownValue;
    }
    else if (x > 1.)                /* Step function = 1 above step */
    { 
        result = UpValue;
    }
    else                            /* Step function: polynomial in between */
    {	
        y = 0.5 + x / 16. * (15. + x * x * (-10. + 3. * x * x));
        
        result = (UpValue - DownValue) * y + DownValue;
        
    }  /* if then else */	

    return result;

}  /* Smooth_Step */

DRLIB_END_NAMESPACE
