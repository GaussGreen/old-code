//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CopulaModel.hpp
//
//   Description : Declares the ICopulaModelBase interface, a set of functions
//				   common to all copulas in CCM
//
//   Date        : Oct 2006
//
//----------------------------------------------------------------------------

#ifndef ICOPULAMODELBASE_HPP
#define ICOPULAMODELBASE_HPP

#include "edginc/MCRandom.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class ICopulaModelBase 
{

public:

	virtual const CreditAssetWrapperArray& getNames() const = 0;

	virtual const CleanSpreadCurveArray& getCurves() const = 0;

	virtual const DateTime& getValueDate() const = 0;

	virtual int getNumAssets() const = 0;

	virtual void addNames(const CreditAssetWrapperArray& xnames) = 0;

	virtual bool nameHasDefaulted(int i) const = 0;

};

DECLARE(ICopulaModelBase);

DRLIB_END_NAMESPACE

#endif
