//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotCrossDerivative.hpp
//
//   Description : Interface for sensitivities that involve tweaking two
//   (or more?) pieces of market data (eg FXCrossGamma, CrossGamma) in order 
//   to tweak a cross derivative. Used by 'Quick X Gamma' in Monte Carlo
//
//   Author      : Mark A Robson
//
//   Date        : 8 November 2002
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/SpotCrossDerivative.hpp"

DRLIB_BEGIN_NAMESPACE

ISpotCrossDerivative::~ISpotCrossDerivative(){}
CClassConstSP const ISpotCrossDerivative::TYPE =
CClass::registerInterfaceLoadMethod(
    "ISpotCrossDerivative", typeid(ISpotCrossDerivative), 0);


DRLIB_END_NAMESPACE
