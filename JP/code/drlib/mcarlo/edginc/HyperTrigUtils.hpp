//----------------------------------------------------------------------------
//
//   Group       : xAsset SRM
//
//   Filename    : HyperTrigUtils.hpp 
//
//   Description : Tools for manipulating hyper trig local vol
//

//
//----------------------------------------------------------------------------


#ifndef HYPER_TRIG_UTILS_HPP
#define HYPER_TRIG_UTILS_HPP
#include <map>
#include <vector>
#include "edginc/Maths.hpp"
DRLIB_BEGIN_NAMESPACE
   
class QuadraticPoly  
{
public:

    double m_deg0, m_deg1, m_deg2;
    double evaluate(double x) const
    {
        //Note: two multiplications rather than three with usual factorisation.
        return (m_deg0 +  x * (m_deg1 + m_deg2 * x));
    }
};

//The quadratic polys are truncated Taylor series
//centred on the x_i: polynomial indeterminate 'is' (x - x_i)
class FunctionInterpolant
{
public:
    double          m_centre;
    QuadraticPoly   m_taylor_series;
};


//typedef std::pair<double, HyperTrigUtils::QuadraticPoly> KInverseInterpolant;   
/** Models the state of the hyperbolic local vol 1 + a tanh cx + b(1 - sech cx),
    and functions involving the local vol */
class MCARLO_DLL HyperLocalVolState 
{
    //Size of table for inverting 'K'
    int m_inverseTableSize;
public:

    typedef std::vector<FunctionInterpolant> InterpolantMapType;

    HyperLocalVolState(double a, double b, double c, int table_size = 20);

    //store the a,b,c smile parameters which this inverse map refers to
    double m_a, m_b, m_c;

    InterpolantMapType m_quadraticInterpolants;

    /** Inverse of KFuncExp for smile params given by idx smileParamIdx in member KFuncInverses */
    void evaluateInverseKExp(double Kt, double *val) const;

    /** Calculates the general smile mapping function K(x)
    Taken from: Hyb3_Kfunc */
    void evaluateK(double  M, double  logM, double  *val) const;

    /** Related to KFunc by KFuncExp(Z) = KFunc(e^Z) */
    void evaluateKExp(double Z, double  *val) const;

    /** 1 + a1 * tanh(a3*u) + a2 * (1 - sech(a3*u)),
        and its derivative */
    void evaluateLocalVol(double u, double *s, double *sdash) const;

private:

    /** Uses interval bisection to locate the interpolant with centre closest to Kt */
    int findClosestInterpolantIdx(double Kt) const;

    /** Create interpolating polynomial and add to the vector.
    x is the image of y under KfuncExp. */
    void addInterpolant(double x, double y, double sm, double smDash); 
};

DRLIB_END_NAMESPACE

#endif
