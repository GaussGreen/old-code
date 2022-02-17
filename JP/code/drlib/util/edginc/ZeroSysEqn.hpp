//----------------------------------------------------------------------------
//
//   Group       : Core Analytics Team
//
//   Filename    : ZeroSysEqn.hpp
//
//   Description : Wrapper around IMSL_ZERO_SYS_EQN
//
//   Author      : Richard Appleton
//
//   Date        : 19th June 2006
//
//----------------------------------------------------------------------------
#include "edginc/AtomicArray.hpp"
#include "edginc/IMSLException.hpp"


#ifndef EDR_ZERO_SYS_EQN_HPP
#define EDR_ZERO_SYS_EQN_HPP

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL ZeroSysEqnFunc
{
public:
    virtual void operator()(int n, double* x, double* errors) = 0;
};

class UTIL_DLL ZeroSysEqnSolver
{
public:
    ZeroSysEqnSolver(double tolerance, int maxIterations);
    virtual ~ZeroSysEqnSolver();
    DoubleArraySP solve(ZeroSysEqnFunc& f, DoubleArray& guess); // throws IMSLException

private:
    ZeroSysEqnSolver();

    double tolerance;
    int    maxIterations;
};


DRLIB_END_NAMESPACE

#endif  // ZERO_SYS_EQN
