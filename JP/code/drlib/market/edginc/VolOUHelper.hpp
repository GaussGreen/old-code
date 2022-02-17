//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolOUHelper.hpp
//
//   Description : Base class for OU-type vol helpers
//
//   Date        : 28 April 03
//
//
//----------------------------------------------------------------------------

#ifndef VOL_OU_HELPER_HPP
#define VOL_OU_HELPER_HPP

#include "edginc/Object.hpp"
#include "edginc/Complex.hpp"

DRLIB_BEGIN_NAMESPACE

/** Base class for OU-type vol helpers */
class MARKET_DLL VolOUHelper: public CObject{
public:
    static CClassConstSP const TYPE;

    /** Computes comp1 and comp2 in 
            E_0 exp(u1 * log(S_t / F_t) + u2 * sigma^2(t) + u3 \int_0^t sigma^2(s)ds) 
                = exp(comp1 + u1 * log(S_0 / F_0) + comp2 * sigma^2(0) */
    void calcJointCumulantComponents(double         tau,
                                     const Complex& u1,     // log spot
                                     const Complex& u2,     // instantaneous variance
                                     const Complex& u3,     // integrated variance
                                     Complex&       comp1,
                                     Complex&       comp2);

    /** Computes comp1 and comp2 in 
            E_0 exp(u1 * log(S_t / F_t) + u2 * sigma^2(t) + u3 \int_0^t sigma^2(s)ds) 
                = exp(comp1 + u1 * log(S_0 / F_0) + comp2 * sigma^2(0) 
        where S is the weighted/scaled log spot - cf superposition */
    void calcJointCumulantComponents(double         tau,
                                     const Complex& u1,     // log spot
                                     const Complex& u2,     // instantaneous variance
                                     const Complex& u3,     // integrated variance
                                     double         weight,
                                     Complex&       comp1,
                                     Complex&       comp2);

protected:
    VolOUHelper(const CClassConstSP& clazz,
                double               lambda,
                double               rho);

    VolOUHelper(const CClassConstSP& clazz);

    /** Computes lambda * \int_0^t k(f(s)) where
        k(x) = nu * x / (alpha - x), alpha > x and    
        f(s) = c1 + c2 exp(-lambda(t-s)) */
    virtual Complex integral(double         tau,
                             const Complex& c1,
                             const Complex& c2) = 0;

    /** Returns log E_0 exp(u Z(1)) = nu * x / (alpha - x) for alpha > x*/
    virtual Complex cumulantBDLP(const Complex& u) = 0;

    // fields
    double  lambda;
    double  rho;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif
