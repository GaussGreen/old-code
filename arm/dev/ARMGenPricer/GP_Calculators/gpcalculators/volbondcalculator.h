#ifndef _INGPCALCULATORS_VOLBONDCALCULATOR_H
#define _INGPCALCULATORS_VOLBONDCALCULATOR_H

#include "gencalculator.h"
#include "gpnumlib/typedef.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_SamplerBase;		
class ARM_ImpSampler;
class ARM_PathScheme;

class ARM_VolBondCalculator : public ARM_GenCalculator
{

public :

	enum ColAlias
	{
		ResetDate = 0,
		StartDate,
		FirstResetDate,
		FirstStartDate,
		PayDate,
		Nominal,
		IT,
		ZcPay,
		CMSRate,
		MinRate,
		MaxRate,
		Coupon	
	};
/*
	enum ProductToPriceCols
	{
		CouponToPrice,
	};
*/
	enum mdmKeysAlias
    {
        ZcCurveKey =  0,
		VolCurveKey,
        MrsKey,
		VolOfVolKey,
		CorrelationKey,
		VolMrsKey,

        NbKeys
    };

	enum payoffType
	{
		TypeI = 0,
		TypeII,
		LookBack
	};



	ARM_VolBondCalculator
		(
			// market data manager
			const ARM_MarketData_ManagerRep&	MktDataManager,
			const ARM_StringVector&				MdmKeys,

			// deal parameters
			const double &						Nominal,

			const ARM_Date	&					StartDate,
			const ARM_Date	&					EndDate,

			const long &						PayFreq,
			const long &						ResetFreq,

			const long &						DayCount,
			const string &						Tenor,

			const long &						IntRule,
			const long &						StubRule,

			const double &						ResetGap,

			const string &						PayCalendar,
			const string &						ResetCalendar,

			// modelparams: Runge-Kutta
			const long &	 					SolverType, 
			const ARM_GP_Vector &				SolverParams, 

			// modelparams: MonteCarlo
			const long &						NbSteps,
			const long &						StepsPerYear,
			const long &		 				BucketSize,
			const vector<string> &				RandGenerator,

			const int &							payofftype,

			// output
			const std::deque<bool>&				ProductsToPrice		
		);

	ARM_VolBondCalculator(const ARM_VolBondCalculator&);
	
	ARM_VolBondCalculator& operator=(const ARM_VolBondCalculator&);
	
	virtual ~ARM_VolBondCalculator();

	virtual ARM_Object*	Clone() const;

	virtual void		CreateAndSetModel();
	virtual void		CreateAndSetCalibration();
	virtual void		Calibrate();
    virtual double		Price();
    virtual void		CheckData(); /// Check internal data consistency
	virtual void		CheckCalibrationData(); /// Check calibration data consistency
	virtual void		CheckMktData(); // Check market data consistency
	virtual ARM_RowInfo	MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	
	virtual ARM_RowInfo ColumnNames() const ;
	virtual ARM_DateStripCombiner DatesStructure() const ;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const ;
	virtual void ComputePricingData() const;
    virtual void UpdateModel() ;
    virtual void UpdateCalibration(bool isUpdateStrike=true) ;

	 virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  ;
	 virtual ARM_DateStripPtr GetRefFundDateStrip() const ;
	 virtual ARM_DateStripPtr GetRefDateStrip() const ;
	 virtual ARM_GP_Vector GetvCpnNominal() const  ;
	 virtual ARM_GP_Vector GetvFundNominal() const ;
	 virtual ARM_GP_Vector GetvFundSpread() const ;

	 virtual ARM_ZeroCurve* GetDomesticZeroCurve() const  ;
	 virtual ARM_ZeroCurve* GetForeignZeroCurve() const  ;
	 virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  ;
	 virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve() const  ;
	 virtual ARM_Forex* GetForex() const  ;

private:
	void							SetColNamesTable();
	size_t							Index(const char* str) const;
	virtual void					InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;	

	vector<string>					ColNamesTable;

	ARM_Currency					itsCcy;

	bool							itsHasBeenPriced;

	// deal parameters
	double							itsNominal;

	ARM_Date						itsStartDate;
	ARM_Date						itsEndDate;

	long							itsPayFreq;
	long							itsResetFreq;

	long							itsDayCount;
	string							itsTenor;

	long							itsIntRule;
	long							itsStubRule;

	double							itsResetGap;

	string							itsPayCalendar;
	string							itsResetCalendar;

	int								itsPayoffType;

	// modelparams: Runge-Kutta
	long							itsSolverType; 
	ARM_GP_Vector					itsSolverParams; 

	// modelparams: MonteCarlo
	long							itsNbSteps;
	long							itsStepsPerYear;
	long							itsBucketSize;
	vector<string>					itsRandGenerator;

	// output
	std::deque<bool>				itsProductsToPrice;

	// Results
	double							itsPrice;

	double							itsFundingPrice;
	double							itsGlobalCapPrice;
	double							itsCouponPrice;
	double							itsSwapPrice;
	double							itsProductPrice;
	double							itsSumCoupPrice;
	double							itsDFPrice;
	double							itsFundRatePrice;
	double							itsIndexPrice;

	mutable ARM_CountedPtr<ARM_DateStripCombiner>		itsEventScheduleForwardStart;
};

CC_END_NAMESPACE()

# endif

