/*
Filename: RadarRepCreator.hpp
Description:  This class is instantiated inside 'KRadarRepGenerator' to create proxies (which 
are objects which can be used to estimate approximately the value of a complex product
using some 'underlying fitting variable' products.  See RegressionRadarRep.cpp for an implementation
of a radar.
Author : Kranthi K. Gade/Vladimir A. Grebinskiy
Date: July 14, 2006.
*/

#ifndef _RADAR_REP_CREATOR_HPP
#define _RADAR_REP_CREATOR_HPP

#include "edginc/IRadarRepCreator.hpp"
#include "edginc/IRadarRep.hpp"
#include "edginc/IFunctionBasis.hpp"
#include "edginc/IFittingVarTransform.hpp"
#include "edginc/regression.h"

DRLIB_BEGIN_NAMESPACE

class  RADAR_DLL RadarRepCreator : public virtual IRadarRepCreator {
public:
    RadarRepCreator(
        IFunctionBasisSP f_cont, 
        IFittingVarTransformSP transform) :
            m_transform(transform), 
            m_f(f_cont), regressionDone(false) 
    {
        
    }

    //virtual IRadarRepSP getRadarRep() const {
    virtual IRadarRepSP getRadarRep();

    /*This function takes in a vector of fitting variables, applies fitting variable transformation
    to them before evaluating the basis functions at those points -- the resulting vector
    is added to 'basisValues' vector. 
    */
    virtual void addSample(FittingArray fittingArr, double val);

    virtual void clearSamples(void);
   
private:
    IFittingVarTransformSP m_transform;
    IFunctionBasisSP m_f;

    IRadarRepSP myRadarRep;
    bool regressionDone;

    vector<vector<double> > basisValues;
    vector<double> dealValues;
};

DECLARE_REF_COUNT(RadarRepCreator);

DRLIB_END_NAMESPACE

#endif
