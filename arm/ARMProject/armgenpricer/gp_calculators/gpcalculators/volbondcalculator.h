/*!
 *
 * Copyright (c) NATIXIS July 2007 Paris
 *
 *	\file gencsocalculator.cpp
 *
 *  \brief file for the virtual calculator for Vol Bond Calculator
 *
 *	\author  Frédéric Degraeve
 *	\version 1.1
 *	\date July 2007
 */

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

	// ZcCurveKey should be equal to ZcCurveKey2
	// VolCurveKey should be equal to VolCurveKey2
	// etc...

	enum mdmKeysAlias
    {
        ZcCurveKey =  0,

		VolCurveKey,
		VolOfVolKey,
		CorrelationKey,
        MrsKey,
		VolMrsKey,

        NbKeys
    };

	enum mdmKeysCalibrationAlias
    {
        ZcCurveKey2 =  0,

		VolCurveKey2,
		VolOfVolKey2,
		CorrelationKey2,
        MrsKey2,
		VolMrsKey2,

		BSModelKey,
		VolCalibKey,
		VovCalibKey,
		CorrelCalibKey,

        NbKeysCalib
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
			const double &						Strike,
			const double &						Leverage,

			const ARM_Date	&					StartDate,
			const ARM_Date	&					EndDate,

			const long &						PayFreq,
			const long &						ResetFreq,

			const long &						PayDayCount,
			const long &						ResetDayCount,

			const string &						Tenor,

			const long &						IntRule,
			const long &						StubRule,

			const double &						ResetGap,

			const string &						PayCalendar,
			const string &						ResetCalendar,

			// modelparams: Runge-Kutta
			const long &	 					SolverType, 
			const std::vector<double> &				SolverParams, 

			// modelparams: MonteCarlo
			const long &						NbSteps,
			const long &						StepsPerYear,
			const long &		 				BucketSize,
			const vector<string> &				RandGenerator,

			const int &							payofftype,
			const string&						underlying,

			const bool&							CalibrationFlag,

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
	 virtual std::vector<double> GetvCpnNominal() const  ;
	 virtual std::vector<double> GetvFundNominal() const ;
	 virtual std::vector<double> GetvFundSpread() const ;

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
	double							itsStrike;
	double							itsLeverage;

	ARM_Date						itsStartDate;
	ARM_Date						itsEndDate;

	long							itsPayFreq;
	long							itsResetFreq;

	long							itsPayDayCount;
	long							itsResetDayCount;

	string							itsTenor;

	long							itsIntRule;
	long							itsStubRule;

	double							itsResetGap;

	string							itsPayCalendar;
	string							itsResetCalendar;

	int								itsPayoffType;
	string							itsUnderlying;

	// modelparams: Runge-Kutta
	long							itsSolverType; 
	std::vector<double>					itsSolverParams; 

	// modelparams: MonteCarlo
	long							itsNbSteps;
	long							itsStepsPerYear;
	long							itsBucketSize;
	vector<string>					itsRandGenerator;

	bool							itsCalibrationFlag;

	// output
	std::deque<bool>				itsProductsToPrice;

	// Results
	double							itsPrice;

	mutable ARM_CountedPtr<ARM_DateStripCombiner>		itsEventScheduleForwardStart;
};

CC_END_NAMESPACE()

# endif

