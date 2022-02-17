//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TwoSidedDeriv.cpp
//
//   Description : Interface for sensitivities that are calculated via a two
//                 sided tweaking algorithm
//
//   Author      : Mark A Robson
//
//   Date        : 11 July 2002
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/TwoSidedDeriv.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for sensitivities that are calculated via a two sided
    tweaking algorithm (eg Delta, FXCrossGamma, CrossGamma) */

ITwoSidedDeriv::~ITwoSidedDeriv(){}
ITwoSidedDeriv::ITwoSidedDeriv() {}

DRLIB_END_NAMESPACE
