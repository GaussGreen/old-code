//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyStreamSchedule.hpp
//
//   Description : Energy Swap stream schedule
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergyStreamSchedule_H_
#define _EnergyStreamSchedule_H_

#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyDWMYPeriod.hpp"
#include "edginc/EnergySwapCoupon.hpp"
#include "edginc/EnergyContractLabel.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL EnergyStreamSchedule : public CObject
{

public:

	static CClassConstSP const TYPE;
	

	friend class EnergyStreamScheduleHelper;

    virtual ~EnergyStreamSchedule();

    EnergyUnderlyerSP getEnergyUnderlyer() const;

    EnergyStreamSchedule(const EnergyUnderlyerSP&   energyUnderlyer,
                         const DateTime&            today,
                         const DateTime&            startDate,
                         const DateTime&            endDate,
                         int                        settleDays,
                         const EnergyDWMYPeriod&    paymentFrequency,
	                     int                        avgDaysB4End,
	                     const EnergyDWMYPeriod&    averagingPeriod,
                         const EnergyDWMYPeriod&    averagingFrequency,
                         int                        rollDay,
                         const string&              notionalSetting);


	// a version for children
	EnergyStreamSchedule(CClassConstSP clazz,
                         const EnergyUnderlyerSP&   energyUnderlyer,
	                     const DateTime&            today,
                         const DateTime&            startDate,
                         const DateTime&            endDate,
	                     int                        settleDays,
	                     const EnergyDWMYPeriod&    paymentFrequency,  
	                     int                        avgDaysB4End, // 0, 1, ...days from END
	                     const EnergyDWMYPeriod&    averagingPeriod,
	                     const EnergyDWMYPeriod&    averagingFrequency,
                         int                        rollDay,
	                     const string&              notionalSetting);


	virtual IObject* clone() const;

	virtual EnergySwapCoupon* createCouponPayment(const DateTime& today,
	                                  const EnergyUnderlyerSP& energyUnderlyer,
	                                  const DateTime& couponStart, 
	                                  const DateTime& couponEnd, 
	                                  const DateTime& notionalStart,
	                                  const DateTime& notionalEnd,
	                                  const DateTime& cashflowDate,
	                                  const EnergyDWMYPeriod& averagingPeriod,
	                                  const EnergyDWMYPeriod& averagingFrequency,
									  int   avgDaysB4End){return 0;};

	virtual double getPV(const EnergyFuturesCurveSP& energyFuturesCurve,
		                                       const YieldCurveSP& yieldCurve, 
			                                   double rate,
		                                       double notionalAmount,
											   const string& notionalType,
											   const string& dealDerection){return 0.0;}

	virtual double  getPV(const EnergyFuturesCurveSP& energyFuturesCurve,
		                                           const YieldCurveSP& yieldCurve, 
								                   double notionalAmount, 
												   const string& notionalType,
								                   const string& dealDerection,
								                   double pastAverage, 
												   bool pastAverageInclToday){return 0.0;}


    static void load(CClassSP& classToLoad){
           // empty for now
    }

	virtual void getScheduleDetails(
		StringArray& thePaymentDates,
		StringArray& theCouponStartDates,
		StringArray& theCouponEndDates,
		StringArray& theNotionalStartDates,
		StringArray& theNotionalEndDates,
		IntArray& theNumFixings,
		StringArray& theFixingDates,
		StringArray& theFixingLabels,
		DoubleArray& theFixingRates ) const ;

	void buildStreamSchedule();

	vector<EnergySwapCouponSP> getCoupons() const;

protected:
	
	EnergyStreamSchedule(CClassConstSP clazz);
	EnergyStreamSchedule();

    EnergyUnderlyerSP          energyUnderlyer;
    DateTime                   today;
    DateTime                   startDate;
	DateTime                   endDate;
    int                        settleDays;
    EnergyDWMYPeriod           paymentFrequency;
    int                        avgDaysB4End;
  	EnergyDWMYPeriod           averagingPeriod;
    EnergyDWMYPeriod           averagingFrequency;
    int                        rollDay; 
    string                     notionalSetting;

	// derived members
	DateTime                   rollDate;
	DateTime                   alignDate;
	int                        hasStub;
	DateTime                   currentPeriodStart;

	vector<EnergySwapCouponSP>   coupons; // $unregistered

private:


    void buildMonthlySchedule();

    DateTime findRollDate();

};

typedef MarketWrapper<EnergyStreamSchedule> EnergyStreamScheduleWrapper;
typedef smartPtr<EnergyStreamSchedule> EnergyStreamScheduleSP;
typedef smartConstPtr<EnergyStreamSchedule> EnergyStreamScheduleConstSP;

//****************************************************************************
// FIXED
//****************************************************************************

class PRODUCTS_DLL EnergyStreamScheduleFixed : public EnergyStreamSchedule
{

public:

	static CClassConstSP const TYPE;

	friend class EnergyStreamScheduleFixedHelper;


    EnergyStreamScheduleFixed(const EnergyUnderlyerSP&   energyUnderlyer,
                              const DateTime&            today,
                              const DateTime&            startDate,
                              const DateTime&            endDate,
                              int                        settleDays,
                              const EnergyDWMYPeriod&    paymentFrequency,
                              int                        avgStartDaysB4End,
  	                          const EnergyDWMYPeriod&    averagingPeriod,
                              const EnergyDWMYPeriod&    averagingFrequency,
                              int                        rollDay,
                              const string&              notionalSetting);

    EnergySwapCoupon* createCouponPayment(
	                     const DateTime& /*today*/,
	                     const EnergyUnderlyerSP& energyUnderlyer,
	                     const DateTime& couponStart, 
	                     const DateTime& couponEnd, 
	                     const DateTime& notionalStart,
	                     const DateTime& notionalEnd,
	                     const DateTime& cashflowDate,
	                     const EnergyDWMYPeriod& averagingPeriod,
	                     const EnergyDWMYPeriod& averagingFrequency,
	                     int   avgDaysB4End );

    double getPV(const EnergyFuturesCurveSP& energyFuturesCurve,
		                      const YieldCurveSP& yieldCurve, 
			                  double fixedRate,
		                      double notionalAmount, 
							  const string& notionalType,
							  const string& dealDerection);
	
    EnergyStreamScheduleFixed();

    ~EnergyStreamScheduleFixed();

	void getScheduleDetails(
		StringArray& thePaymentDates,
		StringArray& theCouponStartDates,
		StringArray& theCouponEndDates,
		StringArray& theNotionalStartDates,
		StringArray& theNotionalEndDates,
		IntArray& theNumFixings,
		StringArray& theFixingDates,
		StringArray& theFixingLabels,
		DoubleArray& theFixingRates ) const ;

private:

	double rate; // only for FIXED_LEG_DETAILS $unregistered
};

typedef MarketWrapper<EnergyStreamScheduleFixed> EnergyStreamScheduleFixedWrapper;
typedef smartPtr<EnergyStreamScheduleFixed> EnergyStreamScheduleFixedSP;
typedef smartConstPtr<EnergyStreamScheduleFixed> EnergyStreamScheduleFixedConstSP;

//***************************************************************************
// FLOATING STREAM .......
//***************************************************************************

class PRODUCTS_DLL EnergyStreamScheduleFloating : public EnergyStreamSchedule
{

public:

	static CClassConstSP const TYPE;

	friend class EnergyStreamScheduleFloatingHelper;

	
    EnergyStreamScheduleFloating(const EnergyUnderlyerSP&   energyUnderlyer,
                                 const DateTime&            today,
                                 const DateTime&            startDate,
                                 const DateTime&            endDate,
                                 int                        settleDays,
                                 const EnergyDWMYPeriod&    paymentFrequency,
                                 int                        avgDaysB4End,
                                 const EnergyDWMYPeriod&    averagingPeriod,
                                 const EnergyDWMYPeriod&    averagingFrequency,
                                 int                        rollDay,
                                 int                        nearbyRel,
                                 const EnergyContractLabel& nearbyAbsLabel,
                                 const string&              notionalSetting);

 
    EnergySwapCoupon* createCouponPayment(const DateTime& today,
	                                  const EnergyUnderlyerSP& energyUnderlyer,
	                                  const DateTime& couponStart, 
	                                  const DateTime& couponEnd, 
	                                  const DateTime& notionalStart,
	                                  const DateTime& notionalEnd,
	                                  const DateTime& cashflowDate,
	                                  const EnergyDWMYPeriod& averagingPeriod,
	                                  const EnergyDWMYPeriod& averagingFrequency,
	                                  int   avgDaysB4End);

    
	IObject* clone() const;

	double getPV(const EnergyFuturesCurveSP& energyFuturesCurve,
		                                           const YieldCurveSP& yieldCurve, 
								                   double notionalAmount, 
												   const string& notionalType,
								                   const string& dealDerection,
								                   double pastAverage, 
								                   bool pastAverageInclToday);

    EnergyStreamScheduleFloating();

    virtual ~EnergyStreamScheduleFloating();

	void getScheduleDetails(
		StringArray& thePaymentDates,
		StringArray& theCouponStartDates,
		StringArray& theCouponEndDates,
		StringArray& theNotionalStartDates,
		StringArray& theNotionalEndDates,
		IntArray& theNumFixings,
		StringArray& theFixingDates,
		StringArray& theFixingLabels,
		DoubleArray& theFixingRates ) const ;

private:

	int             nearbyRel;
    EnergyContractLabel    nearbyAbsLabel;

};


typedef MarketWrapper<EnergyStreamScheduleFloating> EnergyStreamScheduleFloatingWrapper;
typedef smartPtr<EnergyStreamScheduleFloating> EnergyStreamScheduleFloatingSP;
typedef smartConstPtr<EnergyStreamScheduleFloating> EnergyStreamScheduleFloatingConstSP;

DRLIB_END_NAMESPACE

#endif


