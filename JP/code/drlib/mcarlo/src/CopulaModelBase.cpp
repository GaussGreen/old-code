//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CopulaModel.cpp
//
//   Description : Declares the CopulaModelBase, data and functions common to
//                 all copulas in CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CopulaModelBase.hpp"
#include "edginc/DefaultRates.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/imsl.h"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

// clears the array first, for now

void CopulaModelBase::addNames(const CreditAssetWrapperArray& xnames)
{
	names = xnames;
};

DRLIB_END_NAMESPACE

