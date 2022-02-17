//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyVanilla.hpp
//
//   Description : Energy Vanilla instrument
//
//   Author      : Sean Chen
//
//   Date        : July 29,2005
//
//----------------------------------------------------------------------------

#ifndef ENERGY_VANILLA_HPP
#define ENERGY_VANILLA_HPP

#include "edginc/Class.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/ClosedFormEnergy.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/FirstSensDate.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/EnergyContractLabel.hpp"
#include "edginc/EnergyImpliedVolSurface.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

/** Vanilla instrument */
class PRODUCTS_DLL EnergyVanilla: public CInstrument,
                public ClosedFormEnergy::IIntoProduct, 
                public LastSensDate,
				public FirstSensDate
{
    friend class EnergyVanillaHelper;
    friend class EnergyVanillaClosedForm;

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

    /** Implementation of ClosedFormEnergy::IntoProduct interface */
    virtual ClosedFormEnergy::IProduct* createProduct(ClosedFormEnergy* model) const;

    virtual DateTime getValueDate() const;
 
    virtual DateTime endDate(const Sensitivity* sensControl) const;
	virtual DateTime beginDate(const SensControl* sensControl) const;
	string discountYieldCurveName() const;
	void addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const;

	virtual ~EnergyVanilla();

private:

    EnergyVanilla();
    static void load(CClassSP& clazz);

protected:

    EnergyVanilla(CClassConstSP clazz);
    
    EnergyUnderlyerWrapper  energyUnderlyer;
    EnergyFuturesCurveWrapper  futuresCurve;
    string     contractLabel;
    EnergyImpliedVolSurfaceWrapper    volData;
	YieldCurveWrapper       yieldCurve;

    double                  strike;
    bool                    isCall;
    bool                    isEuropean;
    bool                    isBuy;
    double                  notional;

    // optional below.
    DateTime                valueDate;
    DateTime                expiryDate;  // default calced from unerlyer
    DateTime                paymentDate; // default calced from unerlyer
//    InstrumentSettlementSP  instSettle;  // default being physical

    // transient
    DateTime                contractDate;
    
};

typedef smartPtr<EnergyVanilla> EnergyVanillaSP;

DRLIB_END_NAMESPACE
#endif
