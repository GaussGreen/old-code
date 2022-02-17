//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CopulaModel.hpp
//
//   Description : Declares the CopulaModelBase, data and functions common to
//                 all copulas in CCM
//
//   Date        : Feb 2006
//
//----------------------------------------------------------------------------

#ifndef COPULAMODELBASE_HPP
#define COPULAMODELBASE_HPP

#include "edginc/DoubleMatrix.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/BaseSimulation.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/ICopulaModelBase.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 1e-16
#define ZERO_EVENT 10

// An interface class for all copula models: that model crdEvTime dependency 
// among assets

typedef enum{DEPENDENT, INDEPENDENT, RFL, CCM} COPULA_MODEL;

class CopulaModelBase : public virtual VirtualDestructorBase,
					    public virtual ICopulaModelBase
{

public:

	DateTime valueDate;
	CreditAssetWrapperArray names; 

	CleanSpreadCurveArraySP defaultRatesArray;
	array<DoubleArray> recoveries;

	BoolArray hasDefaulted;
	DateTimeArray defaultDates;

public:

	CopulaModelBase(
		DateTime  valDate,
		CleanSpreadCurveArraySP defaultRates
	) :
	valueDate(valDate),
	defaultRatesArray(defaultRates)
	{};

	CopulaModelBase():
		defaultRatesArray(CleanSpreadCurveArraySP(new CleanSpreadCurveArray(0)))
	{};

	virtual const CreditAssetWrapperArray& getNames() const 
		{return names;};

	virtual const CleanSpreadCurveArray& getCurves() const 
		{return *(defaultRatesArray.get());};

	virtual const DateTime& getValueDate() const 
		{return valueDate;};

	virtual int getNumAssets() const {return names.size();};

	virtual void addNames(const CreditAssetWrapperArray& xnames);

	static CleanSpreadCurveSP cleanSpreadCurveFromDF(
		DateTime valueDate,
		DateTimeArray dates,
		DoubleArray   dfs);

	virtual bool nameHasDefaulted(int i) const {return hasDefaulted[i];};

};

DECLARE(CopulaModelBase);

DRLIB_END_NAMESPACE

#endif
