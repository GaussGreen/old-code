/*
Filename: RadarRepCreator.hpp
Description:  This class is instantiated inside 'KRadarRepGenerator' to create proxies (which 
are objects which can be used to estimate approximately the value of a complex product
using some 'underlying fitting variable' products.  See RegressionRadarRep.cpp for an implementation
of a radar.
Author : Kranthi K. Gade/Vladimir A. Grebinskiy
Date: July 14, 2006.
*/

#include "edginc/config.hpp"
#include "edginc/RadarRepCreator.hpp"
#include "edginc/IRadarRep.hpp"
#include "edginc/regression.h"
#include "edginc/IFittingVarTransform.hpp"
#include "edginc/IFunctionBasis.hpp"

DRLIB_BEGIN_NAMESPACE
IRadarRepSP RadarRepCreator::getRadarRep() {
    if (!regressionDone) {
        RegressionCoeff m_coef;
        
        SVDRegression svdReg;
        svdReg.DoTheRegression(basisValues, dealValues, m_coef);
        
        myRadarRep = IRadarRepSP(new RadarRep(m_coef, m_f, m_transform));
        regressionDone = true;
    }
    return myRadarRep;
}

    /*This function takes in a vector of fitting variables, applies fitting variable transformation
    to them before evaluating the basis functions at those points -- the resulting vector
    is added to 'basisValues' vector. 
    */
void RadarRepCreator::addSample(FittingArray fittingArr, double val) {       
    TransformedArray tarr = (*m_transform)(fittingArr);
       
    BasisValArray barr = (*m_f)(tarr);
       
    basisValues.push_back(vector<double>(barr.begin(), barr.end()));
    dealValues.push_back(val);
}

void RadarRepCreator::clearSamples(void) {
    dealValues.clear();
    basisValues.clear();
    regressionDone = false;
}

DRLIB_END_NAMESPACE
