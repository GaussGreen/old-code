//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TermStructure.cpp
//
//   Description : term structure data class
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/TermStructure.hpp"

DRLIB_BEGIN_NAMESPACE

// copy constructor
CTermStructure::CTermStructure(const CTermStructure &source)
{
    *this = source;
}

// assignement operator
CTermStructure& CTermStructure::operator =(const CTermStructure &source)
{
    RefDate = source.RefDate;
    Dates = source.Dates;
    X = source.X;
    Y = source.Y;
    IsTermStructure = source.IsTermStructure;
    return *this;
}

int CTermStructure::size() const 
{
    ASSERT(X.size() == Y.size());
    ASSERT(Dates.size() == X.size() || !IsTermStructure);
    return X.size();
}

DateTime CTermStructure::GetDate(int i) const
{
    if (!IsTermStructure)
        throw ModelException("CTermStructure::GetDate","this is a simple X-Y data with no dates");

    return Dates[i];
}
double CTermStructure::GetX(int i) const
{
    return X[i];
}
double CTermStructure::GetY(int i) const
{
    return Y[i];
}

/**======================================================================
*	populate term structure data
*=======================================================================*/
void CTermStructure::Populate(const DateTime&                  ref, 
                              const int                        size, 
                              vector<DateTime>::const_iterator date,
                              vector<double>::const_iterator   yrs, 
                              vector<double>::const_iterator   val)
{
    if (size <= 0)
    {// data invalid
        throw ModelException("CTermStructure::Populate", "No valid data");
    }
    IsTermStructure = true;
    RefDate = ref;
    Dates.resize(size);
    X.resize(size);
    Y.resize(size);

    for (int i=0; i<size; i++)
    {
        Dates[i] = date[i];
        X[i] = yrs[i];
        Y[i] = val[i];
    }
}

/**======================================================================
*	Name:       :	CTermStructure::Shift
*	Description	:	pt			0=parallel shift
*								>0 shift the point
*					shiftType	1=normal displacement shift
*								2=shift sqrt(Values), used for VolSq
*
*	Return		:	none
*=======================================================================*/
void CTermStructure::Shift(double amount, int pt, int shiftType)
{
    int num = X.size();
    if (num == 0 || pt <0 || pt > num)
    {
        throw ModelException("CTermStructure::Shift", "invalid point or no data.");
    }
    if (shiftType != 1 && shiftType != 2)
    {
        throw ModelException("CTermStructure::Shift", "unknown shift type.");
    }

    if (pt >0) // point shift
    {
        if (shiftType ==1)
            Y[pt-1] += amount;
        else
        {
            Y[pt-1] = sqrt(Y[pt-1]) + amount;
            Y[pt-1] *= Y[pt-1];
        }
    }
    else // parallel
    {
        if (shiftType ==1)
            for (int i=0; i< num; i++)
                Y[i] += amount;
        else
        {
            for (int i=0; i< num; i++)
            {
                Y[i] = sqrt(Y[i]) + amount;
                Y[i] *= Y[i];
            }
        }
    }
}

/** ======================================================================
*	Description	:	generic interpoation
*======================================================================*/
double CTermStructure::Interpolate(double x, InterpType interp)
{
    int num = X.size();
    if (num <1)
    {
        throw ModelException("CTermStructure::Interpolate", "no data.");
        return -1; // never reach here
    }
    else if (num==1)
        return Y[0];

    if(interp == FLAT_FWD)
        return InterpFlatFwd(0.0, x);
    else if (interp == STAIRS)
        return InterpStairs(x);
    else if (interp == CUBIC_SPLINE)
        return InterpCubicSpline(x);
    else if (interp == LINEAR)
        return InterpLinear(x);
    else
        return InterpPoly(x, interp);
}

/*======================================================================
*	Description	:	flat y value till next x value. 
*                   If array has only one value it is returned.
*                   If the requested x value is smaller than that of first array element
*                      thefirst element value returned !!
*======================================================================*/
double CTermStructure::InterpStairs(double t) const
{
    double result;
    int end;

    int num = X.size();
    if (num < 1)
        throw ModelException("CTermStructure::InterpStepReset", "no data.");

    else if (num == 1) 
        return Y[0];

    end = Neighbour(t, &X[0], 0, num-1,  1);
    if (end <0)
        result = Y[num-1];// beyond end point, flat value used
    else
    {
        if (X[end] > t)
            end = Maths::max(0,end-1); // the first item will be returned even if it is after t.
        result = Y[end];
    }
    return result;
}

/*======================================================================
*	Description	:	flat fwd (piece wise linear) interpolation. used for zero rate and VolSquare..
*					interpolate between available points and flat outside.
*=======================================================================*/
double CTermStructure::InterpFlatFwd(double t1, double t2) const
{
    double fwd, v1, v2;
    int i, start1=0, start2=0;

    int num = X.size();
    if (num < 1)
        throw ModelException("CTermStructure::InterpFlatFwd", "no data.");

    else if (num == 1) 
        return Y[0];
    if (t2 < t1 || t1 < 0.0)
        throw ModelException("CTermStructure::InterpFlatFwd", "t2 > t1 and t1 >0 required.");

    if (t2 <= X[0])
        return Y[0];
    if (t1 >= X[num-1]) 
        return Y[num-1];
    if (t2 == t1) // shift t2 by one day if t2 == t1 to help interpolation
        t2 += 1.0/365.0;
    for (i=0; i<num; i++)
    {
        if (X[i] < t1)
            start1 = i;
        if (X[i] < t2)
            start2 = i;
        else
            break;
    }
    if (X[start1] < 0.0)
        throw ModelException("CTermStructure::InterpFlatFwd", "negative time found.");


    fwd = (Y[start1+1]*X[start1+1] - Y[start1]*X[start1])/(X[start1+1]-X[start1]);
    if (t1 <= X[0])
        v1 = Y[0];
    else
        v1 =  (Y[start1]*X[start1] + fwd*(t1-X[start1]))/t1;
    if (start1 == start2)
        v2 =  (Y[start1]*X[start1] + fwd*(t2-X[start1]))/t2;
    else if (start2 == num-1)
        v2 = Y[num-1];
    else
    {
        fwd = (Y[start2+1]*X[start2+1] - Y[start2]*X[start2])/(X[start2+1]-X[start2]);
        v2 =  (Y[start2]*X[start2] + fwd*(t2-X[start2]))/t2;
    }

    return (v2*t2 - v1*t1)/(t2-t1);
}

/*======================================================================
*	Description	:	linear interpolations, linear between available points and flat outside.
*=======================================================================*/
double CTermStructure::InterpLinear(double x) const
{
    double result;
    int start, end;

    int num = X.size();
    if (num < 1)
    {
        throw ModelException("CTermStructure::InterpLinear", "no data.");
        return -1; // never reach here
    }
    else if (num == 1) 
        return Y[0];
    start = Neighbour(x, &X[0], 0, num-1,  -1);
    end = Neighbour(x, &X[0], 0, num-1,  1);

    if (start == end)
        return Y[start];

    if (start <0)
        return Y[0];
    if (end <0)
        return Y[num-1];

    result = Y[start] + (Y[end]-Y[start])*(x - X[start])/(X[end]-X[start]);
    return result;

}

/*======================================================================
*	Description	:	polynomial interpolations.
*=======================================================================*/
double CTermStructure::InterpPoly(double x, int order)
{
    double dy, y;
    int start, end;

    int num = X.size();
    if (num <= order)
        throw ModelException("CTermStructure::InterpPoly", "not enough data points for interpolation.");
    return -1; // never reach here

    start = Neighbour(x, &X[0], 0, num-1,  -1);
    end = Neighbour(x, &X[0], 0, num-1,  1);

    if (start <0 && end <0)
        throw ModelException("CTermStructure::InterpPoly", "no valid data.");

    if (start == end)
        return Y[start];

    start -= order/2;
    end += order/2;

    if (start <0 ) start = 0;
    if (end - start > order )
        end -= 1;
    if (end >= num) end = num -1;

    if (end - start < order)
        throw ModelException("CTermStructure::InterpPoly",
                             "cannot interpolate at this order due to insufficient data points.");

    polint(&*X.begin()-1+start, &*Y.begin()-1+start, 
           end - start + 1, x, &y, &dy);

    return y;
}

/**======================================================================
*	cubic spline for one point. using NR routine
*=======================================================================*/
double CTermStructure::InterpCubicSpline(double x)
{
    const double yp1 = 1e30, ypn = 1e30; // use naural boundaries, see NR
    double y;

    int num = X.size();
    double *y2 = new double[num];

    spline(&*X.begin()-1, &*Y.begin()-1, num, yp1, ypn, y2-1);
    splint(&*X.begin()-1, &*Y.begin()-1, y2-1, num, x, &y);

    delete [] y2;

    return y;
}

/**======================================================================
*	cubic spline for an array. using NR routine.
*=======================================================================*/
void CTermStructure::InterpCubicSpline(const int size, const double* x, double* result)//array form
{
    const double yp1 = 1e30, ypn = 1e30; // use naural boundaries, see NR
    double y;

    int num = X.size();
    double *y2 = new double[num];

    spline(&*X.begin()-1, &*Y.begin()-1, num, yp1, ypn, y2-1);

    for (int i=0; i<size; i++)
    {
        splint(&*X.begin()-1, &*Y.begin()-1, y2-1, num, x[i], &y);
        result[i] = y;
    }

    delete [] y2;
}

DRLIB_END_NAMESPACE
