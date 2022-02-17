//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RootFinderND.hpp
//
//   Date        : 17 June 03
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Range.hpp"
#include "edginc/Function.hpp"
#include "edginc/OneToOneMapping.hpp"
#include "edginc/DECLARE.hpp"

#ifndef EDR_ROOT_FINDER_ND_HPP
#define EDR_ROOT_FINDER_ND_HPP

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL RootFinderND{
public:
    virtual void solve(const FunctionND& func,
                       DoubleArray&      x) const = 0;   // (I/O)
    virtual ~RootFinderND(){}
};

class UTIL_DLL RootFinderNDWithJacobian{
public:
    virtual void solve(const FunctionNDWithJacobian& func,
                       DoubleArray&                  x) const = 0;   // (I/O)
    virtual ~RootFinderNDWithJacobian(){}
};

/** Wrapper around numerical recipes mnewt routine */
class UTIL_DLL RootFinderNDNewton: public RootFinderNDWithJacobian{
public:

    RootFinderNDNewton(int    ntrial,
                       double tolx,
                       double tolf);

    virtual void solve(const FunctionNDWithJacobian& func,
                       DoubleArray&                  x) const;   // (I/O)

    struct Info{
        int iter;
        double errx;    // absolute diff sum
        double errf;    // absolute diff sum
    };
    DECLARE_REF_COUNT(Info);

    InfoConstSP getInfo() const;

private:
    static void usrfun(double *x,int n,double *fvec,double **fjac);

    int    ntrial;
    double tolx;
    double tolf;

    // transient
    mutable InfoSP          info;
    mutable OneToOneMappingConstArray mappings;
    mutable const FunctionNDWithJacobian* func;
    mutable DoubleArray     x;
    mutable DoubleArray     fvec;
    mutable CDoubleMatrixSP fjac;
    // pointer to myself for use in static method
    static const RootFinderNDNewton* me;
};

/** Wrapper around numerical recipes newt routine */
class UTIL_DLL RootFinderNDNewtonSafe: public RootFinderNDWithJacobian{
public:

    RootFinderNDNewtonSafe(int    ntrial,
                           double tolx,
                           double tolf,
                           double stpmx);

    virtual void solve(const FunctionNDWithJacobian& func,
                       DoubleArray&                  x) const;   // (I/O)

    struct Info{
        int    iter;
        double errx;    // max relative diff
        double errf;    // max absolute diff
    };
    DECLARE_REF_COUNT(Info);

    InfoConstSP getInfo() const;

private:
    static void vecfunc(int, double [], double []);
    static void jac(int, double [], double [], double **,
                    void (*)(int, double [], double []));

    int    ntrial;
    double tolx;
    double tolf;
    double stpmx;

    // transient
    mutable InfoSP          info;
    mutable OneToOneMappingConstArray mappings;
    mutable const FunctionNDWithJacobian* func;
    mutable DoubleArray     x;
    mutable DoubleArray     fvec;
    mutable CDoubleMatrixSP fjac;
    // pointer to myself for use in static method
    static const RootFinderNDNewtonSafe* me;
};

DRLIB_END_NAMESPACE

#endif




