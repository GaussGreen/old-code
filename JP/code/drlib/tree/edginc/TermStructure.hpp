//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TermStructure.hpp
//
//   Description : term structure data class
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#ifndef EDG_TERMSTRUCTURE_HPP
#define EDG_TERMSTRUCTURE_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"

// leave here for now
typedef enum {CUBIC_SPLINE=-2, STAIRS=-1, FLAT_FWD=0, LINEAR=1, QUADRATIC=2, POLY3=3, POLY4=4, POLY5=5} InterpType;

DRLIB_BEGIN_NAMESPACE

/////////////////////////////////////////////////////////
// CTermStructure  simple term structure data class
/////////////////////////////////////////////////////////
class TREE_DLL CTermStructure
{
public:
    CTermStructure(){};
	CTermStructure(const CTermStructure &source);
	CTermStructure& operator =(const CTermStructure &source);

    int size() const;

    /** poppulate term structure data */
    void Populate(const DateTime&                  ref, 
                  const int                        size, 
                  vector<DateTime>::const_iterator date,
                  vector<double>::const_iterator   yrs, 
                  vector<double>::const_iterator   val);

    /** interpolation method */
	double  Interpolate(double x, InterpType interp);

    void    Shift(double amount, int pt, int shiftType=1);
	double  InterpFlatFwd(double t1, double t2) const; // flat forward, for zero rate and vol^2
	double  InterpLinear(double x) const;
	double  InterpStairs(double x) const; // for exercise strikes or conversion prices
	double  InterpPoly(double x, int order);
	double  InterpCubicSpline(double x);
	void    InterpCubicSpline(const int size, const double* x, double* result); //array form

    // accesss functions
	DateTime    GetDate(int i) const;
	double      GetX(int i) const;
	double      GetY(int i) const;

private:
    bool                IsTermStructure; // if the data is a term structure (having Dates)
    DateTime            RefDate; // refererence date for XValue=0 when used for term structure data
	vector<DateTime>    Dates;
	vector<double>      X; // measured from RefDate if used for term structure
	vector<double>      Y;
};

DRLIB_END_NAMESPACE
#endif
