//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestLNStrike.cpp
//
//   Description : Abstract LN vol request for a strike-centric world
//
//   Author      : Andrew J Swain
//
//   Date        : 21 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestLNStrike.hpp"
#include "edginc/Maths.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

VolRequestLNStrike::~VolRequestLNStrike(){}

VolRequestLNStrike::VolRequestLNStrike(const CClassConstSP& clazz): 
    CVolRequestLN(clazz){}

/** Returns the sensitive strike */
void VolRequestLNStrike::getSensitiveStrike(double        spot,
                                            DoubleArraySP strikes) const {
    // first get the raw data
    sensitiveStrikes(spot, strikes);

    // for efficiency
    if (strikes->size() > 1){
        // now sort it and remove duplicates
        sort(strikes->begin(), strikes->end());
        // strip duplicates
        vector<double>::iterator iter(strikes->begin());
        for (++iter; iter < strikes->end(); /* inc in loop body */){
            if(Maths::equals(*iter, *(iter-1))) {
                iter = strikes->erase(iter);
            } else {
                ++iter;
            }
        }
    }
}

static void VolRequestLNStrikeLoad(CClassSP& clazz){
    REGISTER(VolRequestLNStrike, clazz);
    SUPERCLASS(CVolRequestLN);
}

CClassConstSP const VolRequestLNStrike::TYPE = CClass::registerClassLoadMethod(
    "VolRequestLNStrike", typeid(VolRequestLNStrike), VolRequestLNStrikeLoad);

// initialise type for array of VolRequestLNStrike
DEFINE_TEMPLATE_TYPE(VolRequestLNStrikeArray);

DRLIB_END_NAMESPACE

