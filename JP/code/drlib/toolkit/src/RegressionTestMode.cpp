//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RegressionTestMode.cpp
//
//   Description : sets up environment used by regression tester
//
//   Author      : Mark A Robson
//
//   Date        : 26 Oct 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RegressionTestMode.hpp"

DRLIB_BEGIN_NAMESPACE

/** turns it on */
void RegressionTestMode::on(){
    inRegressionTester = true;
}

/** turns it off */
void RegressionTestMode::off(){
    inRegressionTester = false;
}

/** is it on or off ? */
bool RegressionTestMode::isOn(){
    return inRegressionTester;
}

bool RegressionTestMode::inRegressionTester = false; // default value

DRLIB_END_NAMESPACE




