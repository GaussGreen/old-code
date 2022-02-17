//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyCapFloor.hpp
//
//   Description : Energy cap instrument
//
//   Author      : Sean Chen
//
//   Date        : Dec 23,2005
//
//----------------------------------------------------------------------------

#ifndef ENERGY_CAPFLOOR_HPP
#define ENERGY_CAPFLOOR_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormEnergy.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/FirstSensDate.hpp"
#include "edginc/EnergySwap.hpp"
#include "edginc/EnergyImpliedVolSurface.hpp"

DRLIB_BEGIN_NAMESPACE

/** Vanilla instrument */
class PRODUCTS_DLL EnergyCapFloor: public CInstrument,
                public ClosedFormEnergy::IIntoProduct, 
                public LastSensDate,
                public FirstSensDate
{
    friend class EnergyCapFloorHelper;
    friend class EnergyCapFloorClosedForm;

public:
    static CClassConstSP const TYPE;

    /** instrument validation */
    virtual void Validate();

    /** input data validation */
    virtual void validatePop2Object();

    /** retrieve market data needed by Vanilla - just valueDate, asset and
        discount yield curve */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);

    pair<double, double> price(bool isSmileOff) const;

    string discountYieldCurveName() const;

    /** Implementation of ClosedFormEnergy::IntoProduct interface */
    virtual ClosedFormEnergy::IProduct* createProduct(ClosedFormEnergy* model) const;

    virtual DateTime getValueDate() const;
 
    virtual DateTime endDate(const Sensitivity* sensControl) const;
    virtual DateTime beginDate(const SensControl* sensControl) const;
    void addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const;

	virtual double oneDCalculateBasketVol(
	                             const vector<int>& weights,
	                             const vector<double>& futures,
	                             const vector<double>& vols,
	                             const vector<double>& t) const;
	virtual double twoDCalculateBasketVol(
	                             const vector<int>& weights,
	                             const vector<double>& futures,
	                             const vector<double>& vols,
	                             const vector<double>& t) const;
	virtual double genericCalculateBasketVol(
	                             const vector<int>& weights,
	                             const vector<double>& futures,
	                             const vector<double>& vols,
	                             const vector<double>& t) const;

    virtual ~EnergyCapFloor();

private:

    EnergyCapFloor();
    static void load(CClassSP& clazz);

protected:

    EnergyCapFloor(CClassConstSP clazz);
    
    EnergySwapSP energySwap;

    EnergyImpliedVolSurfaceWrapper    volData;

    double                  strike;
    bool                    isCap;
    bool                    isBuy;
    double                  notional;

    // optional below.
    DateTime                valueDate; // $unregistered
    DateTime                expiryDate;  // default calced from unerlyer $unregistered
    DateTime                paymentDate; // default calced from unerlyer $unregistered
    
};

typedef smartPtr<EnergyCapFloor> EnergyCapFloorSP;

DRLIB_END_NAMESPACE
#endif
