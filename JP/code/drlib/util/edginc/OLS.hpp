//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OLS.hpp
//
//   Description : Ordinary Least Squares class
//
//----------------------------------------------------------------------------

#ifndef OLS_HPP
#define OLS_HPP

#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL OLS: public CObject {
public:
    static CClassConstSP const TYPE;
    
    /** default constructor */
    OLS();
    
    /** full C++ constructor */
    OLS(const CDoubleMatrixSP& x, const DoubleArraySP& y, bool addIntercept = true);
    
    /** destructor */
    ~OLS() {};

    // validate and populate the OLS object
    virtual void validatePop2Object();
    
    // take x (explanatory variables) and y (explained variable) as DoubleMatrixSP and DoubleArraySP and
    // fill the fields x and y of the OLS object
    // calculate the regression coefficients   
    void refitData(const CDoubleMatrixSP& xx, const DoubleArraySP& yy);

    // return the fitted value at point x
    const double getFittedValues(const DoubleArraySP& xObs) const;

    // return the fitted values for a set of observations x
    const DoubleArray getFittedValues(const CDoubleMatrixSP& xObs) const;

    // return a DoubleArray containing the regression coefficients
    const DoubleArray& getCoefficients() const;

    // return a DoubleArray containing the regression coefficients
    static IObjectSP getAddinCoefficients(OLS* params);

    // return the number of observations
    const int getNbObservations() const;

    // return the number of regressors
    const int getNbRegressors() const;

private:
    // fields
    CDoubleMatrixSP x;  // matrix containing the values of the explanatory variables
    DoubleArraySP y;    // vector containing the values of the explained variable

    // Optional fiels
    bool addIntercept;  // Set to TRUE to append an intercept regressor
    
    // transient fields
    DoubleArray xInternal; // vector <- matrix containing the values of the explanatory variables
    DoubleArray coeffs; // regression coefficients
    
    int nbObservations; // number of observations
    int nbRegressors; // number of explanatory variables

    /** Default constructor */
    static IObject* defaultOLS();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

typedef smartConstPtr<OLS> OLSConstSP;
typedef smartPtr<OLS> OLSSP;
#ifndef QLIB_OLS_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<OLS>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<OLS>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<OLS>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<OLS>);
#endif

DRLIB_END_NAMESPACE

#endif
