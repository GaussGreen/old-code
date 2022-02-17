//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Function.cpp
//
//   Description : 
//
//   Date        : 22 May 02
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Function.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Format.hpp"
#include "edginc/FunctionOperations.hpp"

DRLIB_BEGIN_NAMESPACE

MFunctionND::MFunctionND(int N,      // nb of vars
                         int M):     // nb of funcs
N(N), M(M), 
ranges(*InfiniteRange::createInfiniteRangeArray(N)){} // infinite intervals of definition by default

MFunctionND::MFunctionND(int N,      // nb of vars
                         int M,      // nb of funcs
                         const RangeArray&/*&*/ intervals):   // intervals of definition
N(N), M(M), 
ranges(intervals){}

const RangeArray&/*&*/ MFunctionND::getIntervals() const{
    return ranges;
}

void MFunctionND::setIntervals(RangeArray& _ranges) { ranges = _ranges; }

int MFunctionND::getNbVars() const{
    return N;
}

void MFunctionND::setNbVars(int nbVars) { N=nbVars; setIntervals(*InfiniteRange::createInfiniteRangeArray(N)); }


/**
 * Split this integrand into a vector of "1D" integrands
 * (eg: if integrand is a function f:R^N->R^M, to1DIntegrands() will
 * return M integrands that will be functions R^N->R).
 * 
 * WARNING: This is a generic implementation where each R^N->R integrand
 * still makes a call to original f:R^N->R^M (and then projects result
 * R^M to R), meaning that performance is generally very poor.
 * This method should be overriden on a case by case basis.
 * */
IIntegrandArrayConstSP MFunctionND::to1DIntegrands() const
{
    IIntegrandArraySP result(new IIntegrandArray(M));
    for (int i = 0; i < M; ++i) {
        (*result)[i] = FunctionOperations::project(*this, i);
    }
    return result;
}

int MFunctionND::getNbFuncs() const{
    return M;
}

void MFunctionND::setNbFuncs(int nbFuncs) { M=nbFuncs; }

FunctionNDDouble::FunctionNDDouble(int N):
    MFunctionND(N,      // nb of vars
                1){}    // nb of funcs

FunctionNDDouble::FunctionNDDouble(int N, const RangeArray& intervals):
    MFunctionND(N,           // nb of vars
                1,           // nb of funcs
                intervals){} // intervals of definition

void FunctionNDDouble::operator()(const CDoubleArray&  x,
                                  CDoubleArray&        f) const{
    f[0] = operator()(x);
}

Function1DDouble::Function1DDouble():
    FunctionNDDouble(1){}

Function1DDouble::Function1DDouble(const Range& interval):
    FunctionNDDouble(1,      // nb of vars
                     RangeArray(1, RangeSP( new Range(interval)))){} // intervals of definition


const Range& Function1DDouble::getInterval() const{
    return *(FunctionNDDouble::getIntervals()[0]);
}


double Function1DDouble::operator()(const CDoubleArray&  x) const{
    return operator()(x[0]);
}


Function1DComplex::Function1DComplex():
MFunctionND(1,      // nb of vars
            2){}    // nb of funcs

Function1DComplex::Function1DComplex(const Range& interval):
MFunctionND(1,      // nb of vars
            2,      // nb of funcs
            RangeArray(1, RangeSP(new Range(interval)))){}   // intervals of definition

const Range& Function1DComplex::getInterval() const{
    return *(MFunctionND::getIntervals()[0]);
}


void Function1DComplex::operator()(const CDoubleArray&  x,
                                   CDoubleArray&        f) const{
    Complex rtn(operator()(x[0]));
    f[0] = rtn.real();
    f[1] = rtn.imag();
}

FunctionND::FunctionND(int N):
MFunctionND(N, N){}

FunctionND::FunctionND(int N,
                       const RangeArray& intervals):
MFunctionND(N, N, intervals){}

int FunctionND::getN() const{
    return getNbVars();
}

FunctionNDWithJacobian::FunctionNDWithJacobian(int N):
FunctionND(N){}

FunctionNDWithJacobian::FunctionNDWithJacobian(
   int N,
   const RangeArray& intervals):
FunctionND(N, intervals){}

const FunctionNDWithJacobian* FunctionNDWithJacobian::me = 0;

void FunctionNDWithJacobian::vecfunc(int n, double inx[], double outfvec[]){
    // copy input x over
    int i = 0, inr = 1;
    for (; i < n; ++i, ++inr){
        (*me->x)[i] = inx[inr];
    }
    // compute fvec
    (*me)(*me->x, *me->fvec);
    // copy fvec over
    for (i = 0, inr = 1; i < n; ++i, ++inr){
        outfvec[inr] = (*me->fvec)[i];
    }
}

void FunctionNDWithJacobian::jac(const DoubleArray&  inx,
                                 const DoubleArray&  inf,
                                 DoubleMatrix&       fjac) const{
    // turn fjac into numerical recipe style (ie [1..n][1..n]) matrix
    int n = getN();
    vector<double*> df(n);
    int i = 0;
    for (; i < n; ++i){
        df[i] = fjac[i]-1;
    }
    // allocate memory once for x and fvec vectors
    me = this;
    if (!me->x){
        me->x = DoubleArraySP(new DoubleArray(n));
        me->fvec = DoubleArraySP(new DoubleArray(n));
    }
    // compute jacobian using numerical recipes fdjac
    double* x = const_cast<double*>(&inx[0]-1);
    double* f = const_cast<double*>(&inf[0]-1);
    fdjac(n, x, f, &df[0]-1, vecfunc);
    // done
    me = 0;
}

void FunctionNDWithJacobian::operator()(const DoubleArray&  x,
                                        DoubleArray&        f,
                                        DoubleMatrix&       fjac) const{
    // calculate f
    (*this)(x, f);
    // and fjac
    jac(x, f, fjac);
}

DRLIB_END_NAMESPACE
