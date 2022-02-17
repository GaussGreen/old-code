//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergySwapCoupon.hpp
//
//   Description : Energy Swap Coupon
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#ifndef _ENERGYSWAPCOUPON_H_
#define _ENERGYSWAPCOUPON_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/EnergyFixing.hpp"
#include "edginc/EnergyDWMYPeriod.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL EnergySwapCoupon : public CObject
{

public:

	static CClassConstSP const TYPE;

	friend class EnergySwapCouponHelper;

    EnergySwapCoupon(CClassConstSP clazz,
		const EnergyUnderlyerSP& energyUnderlyer,
        const DateTime& couponStartDate,
        const DateTime& couponEndDate,
        const DateTime& notionalStartDate,
        const DateTime& notionalEndDate,
        const DateTime& cashflowDate,
		const EnergyDWMYPeriod& averagingPeriod, 
        const EnergyDWMYPeriod& averagingFrequency,
        int avgDaysB4End);

    //EnergySwapCoupon(const EnergySwapCoupon& v);

    //EnergySwapCoupon&  operator=(const EnergySwapCoupon& v);

    virtual ~EnergySwapCoupon();

	virtual IObject* clone() const;

    const EnergyUnderlyerSP& getEnergyUnderlyer() const;

    DateTime getCashflowDate() const;

    DateTime getNotionalStartDate() const;

    DateTime getNotionalEndDate() const;

    
    DateTime getCouponStartDate() const;

    DateTime getCouponEndDate() const;

    int getNumFixings() const;

    int getNumTotalCoupons() const;

	virtual void calculateNumFixings() {}

	double getCashflowAmount(double couponRate,
                             double notional,
                             const string& notionalType) const;

	virtual double getRate(const DateTime today,
				   const EnergyFuturesCurveSP& energyFuturesCurve,
		           const string& currency, 
				   double pastAverage, 
				   bool pastAverageInclToday ) { return 0.0;}


    void setCouponStartDate(const DateTime& startDate);
    void setCouponEndDate(const DateTime& endDate);
    void setNumFixings(int numFixings);    
    void setNumTotalCoupons(int numTotalCoupons);

	static void load(CClassSP& classToLoad){
           // empty for now
    }

	virtual DateTime getFixingDate(int index) const { return DateTime();}
	virtual double getFixingRate(int index) const { return 0.0;}
	virtual string getFixingLabel(int index) const { return "";}

	virtual int getNumPastFixings() const { return 0;}
	virtual double getActualPastAverage() const { return 0.0; }

protected:

	EnergySwapCoupon();
	EnergySwapCoupon(CClassConstSP clazz);

	EnergyUnderlyerSP energyUnderlyer;
    DateTime couponStartDate;
    DateTime couponEndDate;
    DateTime notionalStartDate;
    DateTime notionalEndDate;
    DateTime cashflowDate;

	int numFixings;
	int numTotalCoupons;

	EnergyDWMYPeriod averagingPeriod;
    EnergyDWMYPeriod averagingFrequency;
    int avgDaysB4End;

private:

};
    
typedef smartPtr<EnergySwapCoupon> EnergySwapCouponSP;
typedef smartConstPtr<EnergySwapCoupon> EnergySwapCouponConstSP;

//**************************************************************************************
// FIXED LEG
//**************************************************************************************
    
class PRODUCTS_DLL EnergySwapCouponFixed : public EnergySwapCoupon
{

public:

	static CClassConstSP const TYPE;

	friend class EnergySwapCouponFixedHelper;

    EnergySwapCouponFixed(
        const EnergyUnderlyerSP& energyUnderlyer,
        const DateTime& couponStart,
        const DateTime& couponEnd,
        const DateTime& notionalStart, 
        const DateTime& notionalEnd, 
        const DateTime& cashflowDate,
        const EnergyDWMYPeriod& averagingPeriod,
        const EnergyDWMYPeriod& averagingFrequency,
        int avgDaysB4End);
    
    virtual ~EnergySwapCouponFixed();
    
    EnergySwapCouponFixed(const EnergySwapCouponFixed& v);
    EnergySwapCouponFixed& operator=(const EnergySwapCouponFixed& v);
    
    void calculateNumFixings();
 
private:

	EnergySwapCouponFixed();

};
    
typedef smartPtr<EnergySwapCouponFixed> EnergySwapCouponFixedSP;
typedef smartConstPtr<EnergySwapCouponFixed> EnergySwapCouponFixedConstSP;

//***********************************************************
// FLOATING
//***********************************************************

class PRODUCTS_DLL EnergySwapCouponFloating : public EnergySwapCoupon
{

public:

	static CClassConstSP const TYPE;

	friend class EnergySwapCouponFloatingHelper;

    EnergySwapCouponFloating(
        const EnergyUnderlyerSP& energyUnderlyer,
        const DateTime& couponStart, 
        const DateTime& couponEnd,
        const DateTime& notionalStart,
        const DateTime& notionalEnd,
        const DateTime& cashflowDate,
        const EnergyDWMYPeriod& averagingPeriod, 
        const EnergyDWMYPeriod& averagingFrequency,
        int   avgDaysB4End,
		const EnergyFixing& fixing);
    
    virtual ~EnergySwapCouponFloating();
    

	IObject* clone() const;

    // stores fixing dates and resolves contract number to the appropriate contract label
    // contract label is then stored in the fixing, which is more generic than storing fixing number
    // only used for the generic schedule
    void createEnergyFixings(
        const DateTime& today,
        const vector<DateTime>& fixingDates,
        const vector<int>& contractNumbers);    
    
    void calculateFixingDates(const EnergyFixing& fixing);
    
    EnergyFixing createFixing(const EnergyFixing& fixing, 
                              const DateTime& fixingDate) const;

	double getRate(const DateTime today,
				   const EnergyFuturesCurveSP& energyFuturesCurve,
		           const string& currency, 
				   double pastAverage, 
			       bool pastAverageInclToday );

	double determineFixingRate(const DateTime& today,
                               const EnergyFixing& fixing,
                               const EnergyFuturesCurveSP& futuresCurve,
                               const double pastAverage,
	                           bool pastAverageInclToday,
							   bool& isPastAverage) const;

    void calculateNumFixings();

	DateTime getFixingDate(int index) const;
	double getFixingRate(int index) const;
	string getFixingLabel(int index) const;

	int getNumPastFixings() const;
	double getActualPastAverage() const;
    
private:

	EnergySwapCouponFloating();

	vector<EnergyFixing> fixings; // $unregistered

	int numPastFixings; // $unregistered

	double actualPastAverage; // $unregistered

};
    
typedef smartPtr<EnergySwapCouponFloating> EnergySwapCouponFloatingSP;
typedef smartConstPtr<EnergySwapCouponFloating> EnergySwapCouponFloatingConstSP;


DRLIB_END_NAMESPACE

#endif
