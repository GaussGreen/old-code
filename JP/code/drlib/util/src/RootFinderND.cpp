//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RootFinderND.cpp
//
//   Date        : 17 June 03
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RootFinderND.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

// NEWTON-RAPHSON
RootFinderNDNewton::RootFinderNDNewton(
    int    ntrial,
    double tolx,
    double tolf):
ntrial(ntrial),
tolx(tolx),
tolf(tolf){}

RootFinderNDNewton::InfoConstSP RootFinderNDNewton::getInfo() const{
    return info;
}

void RootFinderNDNewton::solve(const FunctionNDWithJacobian& func,
                               DoubleArray&                  inx) const{   // (I/O)
    static const string method("RootFinderNDNewton::solve");
    try{
        int n = inx.size();
        if (n != func.getN()){
            throw ModelException(method,
                                 "func's size ("
                                 + Format::toString(func.getN())
                                 + ") does not match x's ("
                                 + Format::toString(n)
                                 + ")");
        }
        info = InfoSP(new Info());
        me = this;
        // create mappings, map intial vars
        // anc calculate tolerance in mapped vars
        me->mappings.resize(n);
        const RangeArray& intervals = func.getIntervals();
        double maxDeriv = 0.0;
        int i = 0;
        for (; i < n; ++i){
            me->mappings[i] = OneToOneMapping::create(*(intervals[i]));
            inx[i] = me->mappings[i]->inverse(inx[i]);
            double deriv = me->mappings[i]->derivative(inx[i]);
            maxDeriv = Maths::max(maxDeriv, fabs(deriv));
        }
        double scaledtolx = tolx / maxDeriv;
        // initialize transients
        me->func = &func;
        me->x.resize(n);
        me->fvec.resize(n);
        me->fjac = CDoubleMatrixSP(new DoubleMatrix(n, n));
        // solve
        try{
            mnewt(ntrial, 
                  &inx[0]-1, 
                  n, 
                  scaledtolx,
                  tolf,
                  &info->iter,
                  &info->errx,
                  &info->errf,
                  usrfun);
        }
        catch(exception&){
            // map vars back
            for (i = 0; i < n; ++i){
                inx[i] = (*me->mappings[i])(inx[i]);
            }
            me = 0;
            throw;
        }
        // map vars back
        for (i = 0; i < n; ++i){
            inx[i] = (*me->mappings[i])(inx[i]);
        }
        me = 0;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

const RootFinderNDNewton* RootFinderNDNewton::me = 0;

void RootFinderNDNewton::usrfun(double *inx,int n,double *outfvec,double **outfjac){
    // map vars
    int i = 0, inr = 1;
    for (; i < n; ++i, ++inr){
        me->x[i] = (*me->mappings[i])(inx[inr]);
    }
    // call func
    (*me->func)(me->x,
                me->fvec,
                *me->fjac);
    // output and adjust jacobian by the derivatives of the mappings
    int j = 0, jnr = 1;
    for (; j < n; ++j, ++jnr){
        double deriv = me->mappings[j]->derivative(inx[jnr]);
        outfvec[jnr] = me->fvec[j];
        for (i = 0, inr = 1; i < n; ++i, ++inr){
            outfjac[inr][jnr] = (*me->fjac)[i][j] * deriv;
        }
    }
}

// GLOBALLY CONVERGENT NEWTON-RAPHSON
RootFinderNDNewtonSafe::RootFinderNDNewtonSafe(
    int    ntrial,
    double tolx,
    double tolf,
    double stpmx):
ntrial(ntrial),
tolx(tolx),
tolf(tolf),
stpmx(stpmx){}

RootFinderNDNewtonSafe::InfoConstSP RootFinderNDNewtonSafe::getInfo() const{
    return info;
}

void RootFinderNDNewtonSafe::solve(const FunctionNDWithJacobian& func,
                                   DoubleArray&                  inx) const{   // (I/O)
    static const string method("RootFinderNDNewtonSafe::solve");
    try{
        int n = inx.size();
        if (n != func.getN()){
            throw ModelException(method,
                                 "func's size ("
                                 + Format::toString(func.getN())
                                 + ") does not match x's ("
                                 + Format::toString(n)
                                 + ")");
        }
        info = InfoSP(new Info());
        me = this;
        // create mappings and map intial vars
        me->mappings.resize(n);
        const RangeArray& intervals = func.getIntervals();
        double maxDerivTimesRatio = 0.0;
        int i = 0;
        for (; i < n; ++i){
            me->mappings[i] = OneToOneMapping::create(*(intervals[i]));
            double oldx = inx[i];
            inx[i] = me->mappings[i]->inverse(inx[i]);
            double deriv = me->mappings[i]->derivative(inx[i]);
            maxDerivTimesRatio 
                = Maths::max(maxDerivTimesRatio, 
                             fabs(deriv) 
                             * fabs(inx[i]) / Maths::max(1.0, fabs(oldx)));
        }
        double scaledtolx = tolx / maxDerivTimesRatio;
        // initialize remaining transients
        me->func = &func;
        me->x.resize(n);
        me->fvec.resize(n);
        me->fjac = CDoubleMatrixSP(new DoubleMatrix(n, n));
        // solve
        int check;
        try{
            newt(&inx[0]-1, 
                 n,
                 &info->iter,
                 &info->errf,
                 &info->errx,
                 &check,
                 vecfunc,
                 jac,
                 ntrial, 
                 tolf,
                 0.01 * tolf, // TOLMIN
                 scaledtolx,
                 stpmx);
        }
        catch(exception&){
            // map vars back
            for (i = 0; i < n; ++i){
                inx[i] = (*me->mappings[i])(inx[i]);
            }
            me = 0;
            throw;
        }
        // map vars back
        for (i = 0; i < n; ++i){
            inx[i] = (*me->mappings[i])(inx[i]);
        }
        me = 0;
        if (check){ // failure
            throw ModelException(method, 
                                 "newt has converged to a local minimum.\n"
                                 "Try restarting from a different initial guess");
        }
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

const RootFinderNDNewtonSafe* RootFinderNDNewtonSafe::me = 0;

void RootFinderNDNewtonSafe::vecfunc(int n, double inx[], double outfvec[]){
    // map vars
    int i = 0, inr = 1;
    for (; i < n; ++i, ++inr){
        me->x[i] = (*me->mappings[i])(inx[inr]);
    }
    // call func
    (*me->func)(me->x,
                me->fvec);
    // output
    for (i = 0, inr = 1; i < n; ++i, ++inr){
        outfvec[inr] = me->fvec[i];
    }
}

void RootFinderNDNewtonSafe::jac(int n, double inx[], double infvec[], double **outfjac,
                                 void (*notused)(int, double [], double [])){
    // map vars
    int i = 0, inr = 1;
    for (; i < n; ++i, ++inr){
        me->x[i] = (*me->mappings[i])(inx[inr]);
        me->fvec[i] = infvec[inr];
    }
    // call func
    me->func->jac(me->x,
                  me->fvec,
                  *me->fjac);
    // output and adjust jacobian by the derivatives of the mappings
    int j = 0, jnr = 1;
    for (; j < n; ++j, ++jnr){
        double deriv = me->mappings[j]->derivative(inx[jnr]);
        for (i = 0, inr = 1; i < n; ++i, ++inr){
            outfjac[inr][jnr] = (*me->fjac)[i][j] * deriv;
        }
    }
}

DRLIB_END_NAMESPACE
