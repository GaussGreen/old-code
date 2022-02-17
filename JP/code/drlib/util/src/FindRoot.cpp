//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FindRoot.cpp - from CMLib FindRoot.h
//
//   Description : Base stuff for finding roots
//
//   Author      : CMLib
//
//   Date        : Dec 17, 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FindRoot.hpp"

DRLIB_BEGIN_NAMESPACE

const double FindRoot::LOWER_BOUND = -DBL_MAX;
const double FindRoot::UPPER_BOUND = DBL_MAX;
    
const double FindRoot::X_TOLERANCE = 1e-08;
const double FindRoot::F_TOLERANCE = 1e-12;

const int    FindRoot::MAX_NUM_OF_ITERATIONS = 50;
const int    FindRoot::FIXED_NUM_OF_ITERATIONS = 20;

const double FindRoot::INITIAL_X_STEP_SIZE = 1E-8;
const double FindRoot::INITIAL_F_DERIVATIVE = 0;

FindRoot::Point::Point() {}
FindRoot::Point::Point(double x, double y): x(x), y(y){}

DRLIB_END_NAMESPACE


