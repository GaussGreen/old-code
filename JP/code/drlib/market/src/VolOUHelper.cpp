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

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/VolOUHelper.hpp"


DRLIB_BEGIN_NAMESPACE

VolOUHelper::VolOUHelper(const CClassConstSP& clazz,
                         double               lambda,
                         double               rho):
CObject(clazz),
lambda(lambda),
rho(rho){}

VolOUHelper::VolOUHelper(const CClassConstSP& clazz):
CObject(clazz),
lambda(0.0),
rho(0.0){}

/** Computes comp1 and comp2 in 
        E_0 exp(u1 * log(S_t / F_t) + u2 * sigma^2(t) + u3 \int_0^t sigma^2(s)ds) 
            = exp(comp1 + u1 * log(S_0 / F_0) + comp2 * sigma^2(0) */
void VolOUHelper::calcJointCumulantComponents(double         tau,
                                              const Complex& u1,     // log spot
                                              const Complex& u2,     // instantaneous variance
                                              const Complex& u3,     // integrated variance
                                              Complex&       comp1,
                                              Complex&       comp2){
    static const string method = "VolOUHelper::calcJointCumulantComponents";    
    try {
        Complex c = (0.5 * (u1 * u1 - u1) + u3) / lambda;
        Complex c1 = u1 * rho + c;
        Complex c2 = u2 - c;

        double effTime = lambda * tau;
        double expMEffTime = exp(-effTime);

        // call integral routine
        comp1 = - u1 * effTime * cumulantBDLP(rho) + integral(tau, c1, c2);

        comp2 = c * (1.0 - expMEffTime) + u2 * expMEffTime;
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** Computes comp1 and comp2 in 
        E_0 exp(u1 * log(S_t / F_t) + u2 * sigma^2(t) + u3 \int_0^t sigma^2(s)ds) 
            = exp(comp1 + u1 * log(S_0 / F_0) + comp2 * sigma^2(0) */
void VolOUHelper::calcJointCumulantComponents(double         tau,
                                              const Complex& u1,     // log spot
                                              const Complex& u2,     // instantaneous variance
                                              const Complex& u3,     // integrated variance
                                              double         weight,
                                              Complex&       comp1,
                                              Complex&       comp2){
    static const string method = "VolOUHelper::calcJointCumulantComponents";    
    try {
        Complex z1 = weight * u1;
        Complex z3 = u3 - (0.5 * weight * (weight - 1.0)) * u1;
        Complex c = (0.5 * (z1 * z1 - z1) + z3) / lambda;
        Complex c1 = z1 * rho + c;
        Complex c2 = u2 - c;

        double effTime = lambda * tau;
        double expMEffTime = exp(-effTime);

        // call integral routine
        comp1 = - u1 * effTime * cumulantBDLP(rho * weight) + integral(tau, c1, c2);

        comp2 = c * (1.0 - expMEffTime) + u2 * expMEffTime;
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

/** Invoked when Class is 'loaded' */
void VolOUHelper::load(CClassSP& clazz){
    REGISTER(VolOUHelper, clazz);
    SUPERCLASS(CObject);
    FIELD(lambda, "lambda");
    FIELD(rho, "rho");
}

CClassConstSP const VolOUHelper::TYPE =
CClass::registerClassLoadMethod("VolOUHelper", typeid(VolOUHelper), load);

DRLIB_END_NAMESPACE
