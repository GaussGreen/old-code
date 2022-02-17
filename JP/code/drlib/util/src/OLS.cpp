//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OLS.cpp
//
//   Description : Ordinary Least Squares
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_OLS_CPP
#include "edginc/OLS.hpp"
#include "edginc/Addin.hpp"
#include "edginc/imslerror.hpp"

DRLIB_BEGIN_NAMESPACE

/** default constructor */
OLS::OLS(): CObject(TYPE), addIntercept(true) {}

/** full C++ constructor */
OLS::OLS(const CDoubleMatrixSP& x, const DoubleArraySP& y, bool addIntercept) : 
CObject(TYPE), x(x), y(y), addIntercept(addIntercept) {
    validatePop2Object();
}

// validate and populate the OLS object
void OLS::validatePop2Object() {
    static const string routine("OLS::validatePop2Object");
    try {
        // check if the inputs x and y have the right size
        if (this->y->size() != this->x->numRows()) {
            throw ModelException(routine, "the number of observations in x and y are not the same");
        }
        
        // number of explanatory variables and number of observations
        nbObservations = y->size(); // number of observations
        nbRegressors = x->numCols(); // number of explanatory variables

        // vector <- matrix containing the values of the explanatory variables
        xInternal = DoubleArray(nbObservations * nbRegressors);
    
        int jRegressors;
        for (int iObservations=0 ; iObservations<nbObservations ; iObservations++) {
            for (jRegressors=0 ; jRegressors<nbRegressors ; jRegressors++) {
                xInternal[iObservations * nbRegressors + jRegressors] = (*x)[jRegressors][iObservations];
            }
        }

        // calculate regression coefficients using isml regression function
        int nbCoeffs = nbRegressors;
        if(addIntercept) {
            nbCoeffs ++;    // number of coeffs = nbRegressors + 1 (for intercept)
        }

        coeffs = DoubleArray(nbCoeffs);

        double* coefficients;
        if(addIntercept) {
            // IMSL will ask for intercept by default
            coefficients = imsl_d_regression(nbObservations, nbRegressors, 
                &xInternal[0], &(*y)[0], 0);
        } else {
            // Do not request the intercept from IMSL
            coefficients = imsl_d_regression(nbObservations, nbRegressors, 
                &xInternal[0], &(*y)[0], IMSL_NO_INTERCEPT, 0);
        }
        IMSLError::throwExceptionIfError();
        // Interestingly, some IMSL warnings (not errors) return a NULL
        // pointer. Throw an exception in that case
        if(!coefficients) {
            throw ModelException("IMSL OLS regression failed");
        }
        

        for (jRegressors=0 ; jRegressors < nbCoeffs; jRegressors++) {
            coeffs[jRegressors] = coefficients[jRegressors];
        }

        delete[] coefficients;
    
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

// take x (explanatory variables) and y (explained variable) as DoubleMatrixSP and DoubleArraySP and
// fill the fields x and y of the OLS object
// calculate the regression coefficients                           
void OLS::refitData(const CDoubleMatrixSP& x, const DoubleArraySP& y) {
    
    static const string method("OLS::refitData");

    // check if the inputs x and y have the right size
    if (this->y->size() != y->size()) {
        throw ModelException(method, "the number of observations in y is incorrect");
    }
    if (this->x->numRows() != x->numRows()) {
        throw ModelException(method, "the number of observations in x is incorrect");
    }
    if (this->x->numCols() != x->numCols()) {
        throw ModelException(method, "the number of regressors in x is incorrect");
    }
    
    // fill the fields x and y of the OLS object
    this->x = x;
    this->y = y;
    
    // construct the vector xInternal from the matrix x 
    int jRegressors;
    for (int iObservations=0 ; iObservations<nbObservations ; iObservations++) {
        for (jRegressors=0 ; jRegressors<nbRegressors ; jRegressors++) {
            xInternal[iObservations * nbRegressors + jRegressors] = (*x)[jRegressors][iObservations];
        }
    }
    
    // calculate regression coefficients using isml regression function
    int nbCoeffs = nbRegressors + 1; // number of coeffs = nbRegressors + 1 (for intercept)
    
    double* coefficients = imsl_d_regression (nbObservations, nbRegressors, &xInternal[0], &(*y)[0], 0);
    
    for (jRegressors=0 ; jRegressors < nbCoeffs; jRegressors++) {
        coeffs[jRegressors] = coefficients[jRegressors];
    }

    delete[] coefficients;
}

// return the fitted value at point x
const double OLS::getFittedValues(const DoubleArraySP& xObs) const {
    
    static const string method("OLS::getfittedValues");

    // check if the size of xObs is equal to the number of explanatory variables
    if (xObs->size() != nbRegressors) {
        throw ModelException(method, "the number of regressors in xObs is incorrect");
    }

    // calculate the fitted value corresponding to xObs
    double fittedValue = coeffs[0];

    for (int jRegressors=0 ; jRegressors < nbRegressors; jRegressors++) {
        fittedValue += coeffs[jRegressors + 1] * (*xObs)[jRegressors];
    }

    return fittedValue;
}

// return the fitted values for a set of observations x
const DoubleArray OLS::getFittedValues(const CDoubleMatrixSP& xObs) const {
static const string method("OLS::getfittedValues");

    // check if the size of xObs is equal to the number of explanatory variables
    if (xObs->numCols() != nbRegressors) {
        throw ModelException(method, "the number of regressors in xObs is incorrect");
    }

    // number of observations
    int nbObs = xObs->numRows();

    // calculate the fitted value corresponding to xObs
    DoubleArray fittedValues(nbObs);

    for (int iObs=0 ; iObs < nbObs; iObs++) {
        fittedValues[iObs] = coeffs[0];
        for (int jRegressors=0 ; jRegressors < nbRegressors; jRegressors++) {
            fittedValues[iObs] += coeffs[jRegressors + 1] * (*xObs)[jRegressors][iObs];
        }
    }
    return fittedValues;
}

// return the number of observations
const int OLS::getNbObservations() const {
    return nbObservations;
}

// return the number of regressors
const int OLS::getNbRegressors() const {
    return nbRegressors;
}

// return a DoubleArray containing the regression coefficients
const DoubleArray& OLS::getCoefficients() const {
    return coeffs;
}

// return a IObjectSP containing the regression coefficients
IObjectSP OLS::getAddinCoefficients(OLS* params){
    static const string routine = "OLS::getCoefficients";
    try {
        DoubleArraySP coeffs(new DoubleArray(params->getCoefficients()));
        return coeffs;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Default constructor */
IObject* OLS::defaultOLS(){
    return new OLS();
}

// invoked when Class is 'loaded'
void OLS::load(CClassSP& clazz) {
    REGISTER(OLS, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultOLS);
    // mandatory fields
    FIELD(x, "Explanatory variables data matrix");
    FIELD(y, "Explained variable data vector");
    // optional fields
    FIELD(addIntercept, "Add intercept to regressors");
    FIELD_MAKE_OPTIONAL(addIntercept);
    // transient fields
    FIELD(xInternal, "Vector form of matrix x");
    FIELD_MAKE_TRANSIENT(xInternal);
    FIELD(coeffs, "Regression coefficients");
    FIELD_MAKE_TRANSIENT(coeffs);
    FIELD(nbObservations, "Number of observations");
    FIELD_MAKE_TRANSIENT(nbObservations);
    FIELD(nbRegressors, "Number of explanatory variables");
    FIELD_MAKE_TRANSIENT(nbRegressors);
    
    clazz->setPublic(); // make visible to EAS/spreadsheet

    // register addin methos getCoefficients
    Addin::registerClassObjectMethod("OLS_COEFFICIENTS",
                                     Addin::UTILITIES,
                                     "Regression coefficients",
                                     TYPE,
                                     true,
                                     Addin::expandSimple,
                                     (Addin::ObjMethod*)getAddinCoefficients);

}

CClassConstSP const OLS::TYPE = CClass::registerClassLoadMethod(
    "OLS", typeid(OLS), OLS::load);

//////////////////////////////////////////////////////////////////////////////

// ADDIN methods for OLS regression
class OLSAddin: public CObject{
    static CClassConstSP const TYPE;

    OLSSP olsregression;
    CDoubleMatrixSP xObs;
 
    static IObjectSP getFittedValues(OLSAddin* params){
        static const string routine = "OLSAddin::getFittedValues";
        try {
            DoubleArraySP fittedValues(new DoubleArray(params->olsregression->getFittedValues(params->xObs)));
            return fittedValues;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }
    
    // for reflection
    OLSAddin(): CObject(TYPE) {}

    // invoked when Class is 'loaded'
    static void load(CClassSP& clazz){
        REGISTER(OLSAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultOLSAddin);
        FIELD(olsregression, "OLS regression");
        FIELD(xObs, "Explanatory variables");
        
        // register addin methos getFittedValues
        Addin::registerClassObjectMethod("OLS_FITTEDVALUES",
                                         Addin::UTILITIES,
                                         "Fitted values",
                                         TYPE,
                                         true,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)getFittedValues);
    }

    static IObject* defaultOLSAddin(){
        return new OLSAddin();
    }
};

CClassConstSP const OLSAddin::TYPE = CClass::registerClassLoadMethod(
    "OLSAddin", typeid(OLSAddin), load);

DRLIB_END_NAMESPACE
