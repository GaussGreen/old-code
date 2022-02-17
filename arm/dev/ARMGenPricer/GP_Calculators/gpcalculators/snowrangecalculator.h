/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file SnowRangeCalculator.h
 *
 *  \brief
 *
 *	\author  P. LAM
 *	\version 1.0
 *	\date mar 2006
 */

#ifndef _INGPCALCULATORS_SNOWRANGECALCULATOR_H
#define _INGPCALCULATORS_SNOWRANGECALCULATOR_H

#include "gencalculator.h"
#include "gpbase/countedptr.h"

/// kernel
#include "util/refvalue.h"
#include "util/exercise.h"
#include "ccy/currency.h"
#include "inst/portfolio.h"
#include "mod/bssmiled.h"

/// STL
CC_USING_NS(std,pair)

/// forward declaration
class ARM_Portfolio;
class ARM_BSSmiledModel;

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStrip;
class ARM_MCMethod;

class ARM_SnowRangeCalculator : public ARM_GenCalculator
{
protected:

	static const string SnowRangeColNamesTable [];
	static const int SnowRangeProductToPriceColumns [];

public:
	
	enum SnowRangeColAlias
	{
		ResetDate = 0,
		StartDate,
		EndDate,
		DealEndDate,
		PayDate,
		Notional,
		Spread,
		Strike,
		Ratchet,
		CashFlow, 
		FundingBasis,
		CouponBasis,
		Numeraire,
		Index,
		Libor,
		DF,
		Funding,
		RatioADV,
		RatioARR,
		NbPeriod,
		FCorr,
		CouponLev,
		Coupon,
		StdCoupon,
		PlainCoupon,
		Swap,
		CallOption,
		CallSwap,
	};

	enum SnowRangeProductToPriceAlias
	{
		CallSwapPrice,
		CallOptionPrice,
		SwapPrice,
		PlainCouponPrice,
		StdCouponPrice,
		CouponPrice,
		FundingPrice,
		LiborPrice,
		IndexPrice,
		DFPrice,
		RatioADVPrice,
		RatioARRPrice,
		FCorrPrice,
		CouponLevPrice,

		NbProductsToPrice,
	};

	enum mdmKeysAlias
    {
		YcKey =  0,
		CfModelKey,
		OswModelKey,
		MrsKey,
		HumpKey,

		NbKeys
    };

	//Constructors
	ARM_SnowRangeCalculator(const ARM_MarketData_ManagerRep& mktDataManager );

	ARM_SnowRangeCalculator(const ARM_SnowRangeCalculator&);
	
	ARM_SnowRangeCalculator& operator=(const ARM_SnowRangeCalculator&);

	ARM_SnowRangeCalculator(  ARM_Currency& ccy,
							  ARM_Date&	startDate,
							  ARM_Date&	endDate,
							  int payRec,
							  ARM_ReferenceValue& notional,
							  string fundingIndexTerm,
							  int fundingDayCount,
							  string couponIndexTerm,
							  int couponDayCount,
							  int resetFreq,
							  int payFreq,
							  int resetTiming,
							  int payTiming,
							  int resetGap,
							  string resetCal,
							  string payCal,
							  int adjRule,
							  int intRule,
							  ARM_ReferenceValue& spread,
							  ARM_ReferenceValue& strike,
							  ARM_ReferenceValue& ratchet,
							  ARM_ReferenceValue& cashFlow,
							  ARM_ReferenceValue& fixedRate,
							  ARM_ReferenceValue& leverage,
							  ARM_Vector* snowRangeParams,
							  ARM_Vector* calibParams,
							  string modelName,
							  ARM_Vector* modelParams,
							  int nbSteps,
							  string generatorType,
							  string inversionMethod,
							  bool antithetic,
							  int samplerType,
							  ARM_StringVector& mdmKeys,
							  const ARM_MarketData_ManagerRep& mktDataManager,
							  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice);

	//Destructor
	~ARM_SnowRangeCalculator();

	virtual ARM_Object*			Clone() const;
	
	virtual ARM_CstManagerPtr	CreateCstManager();

	ARM_StringVector			PricedColumnNames() const;

	//Outputs	
	virtual ARM_RowInfo				ColumnNames() const;
	virtual ARM_RowInfo				MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	virtual ARM_DateStripCombiner	DatesStructure() const;
	virtual ARM_DateStripCombiner	CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	/// pricing function
	virtual void					ComputePricingData() const;
	virtual void					CreateAndSetModel();
	virtual void					CreateAndSetCalibration();
	void							CreateCalibMethodForSBGM();
    virtual void					UpdateModel();
    virtual void					UpdateCalibration(bool isUpdateStrike=true);
	virtual void					Calibrate();
    virtual double					Price();
    virtual void					CheckData(); /// Check internal data consistency
	virtual void					CheckMktData(); // Check mkt data consistency
	/// Initialisation to 0 of all columns of the deal description that could be priced
    virtual void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	// Some utilities
	virtual string					ExportShortName() const { return "LCSRG";}
	string							DealDesDataToString(const string& indent = "", const string& nextIndent = "") const;
	virtual string					GeneralDataToString(const string& indent = "", const string& nextIndent = "") const;
	virtual void					View(char* id = NULL, FILE* ficOut = NULL) const;
	void							GenerateProductDescription(ARM_Portfolio* portfolio);

	void CleanUp();
    void CopyNoCleanUp(const ARM_SnowRangeCalculator& rhs);

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  { return ARM_DateStripPtr(NULL);};
	inline virtual ARM_DateStripPtr GetRefFundDateStrip() const { return ARM_DateStripPtr(NULL);};
	inline virtual ARM_DateStripPtr GetRefDateStrip() const { return ARM_DateStripPtr(NULL);};

	/// Vector
	inline virtual ARM_GP_Vector GetvCpnNominal() const  { return ARM_GP_Vector(0);};
	inline virtual ARM_GP_Vector GetvFundNominal() const { return ARM_GP_Vector(0);};
	inline virtual ARM_GP_Vector GetvFundSpread() const { return ARM_GP_Vector(0);};

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve() const  { return NULL;};
	inline virtual ARM_Forex* GetForex() const  { return NULL;};
	
protected:

	// General Data
    ARM_Currency				itsCcy;
	ARM_Date					itsStartDate;           
	ARM_Date					itsEndDate;
	ARM_Date					itsSchedEndDate;
	int							itsPayRec;
	ARM_ReferenceValue			itsNotional;

	// Funding Data
	string						itsFundingIndexTerm;
	int							itsFundingDayCount;
	
	// Coupon Data
	string						itsCouponIndexTerm;
	int							itsCouponDayCount;
	int							itsResetFreq;
	int							itsPayFreq;
	int							itsResetGap;
	int							itsPayGap;
	int							itsResetTiming;
	int							itsPayTiming;
	string						itsResetCal;
	string						itsPayCal;
	int							itsAdjRule;
	int							itsIntRule;
	ARM_ReferenceValue			itsSpread;
	ARM_ReferenceValue			itsFixedRate;
	ARM_ReferenceValue			itsLeverage;

	// Range Data
	ARM_ReferenceValue			itsStrike;
	ARM_ReferenceValue			itsRatchet;
	ARM_ReferenceValue			itsCashFlow;
	ARM_Vector*					itsSnowRangeParams;

	// Model Params
	string						itsModelName;
	ARM_Vector*					itsModelParams;
	
	// Calibration Params
	ARM_Vector*					itsCalibParams;
	int							itsNbSteps;
	string						itsGeneratorType;
	string						itsInversionMethod;
	bool						itsAntithetic;
	int							itsSamplerType;

	ARM_NumMethodPtr			itsMCMethod;

	// Pricing Data
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/				itsProductsToPrice;

	bool						itsHasBeenPriced;

	double						itsCallSwapPrice;
	double						itsCallOptionPrice;
	double						itsSwapPrice;
	double						itsPlainCouponPrice;
	double						itsStdCouponPrice;
	double						itsCouponPrice;
	double						itsFundingPrice;
	double						itsLiborPrice;
	double						itsIndexPrice;
	double						itsDFPrice;
	double						itsRatioADVPrice;
	double						itsRatioARRPrice;
	double						itsFCorrPrice;
	double						itsCouponLevPrice;

	// Date Strip
	mutable ARM_DateStripPtr	itsProductDateStrip;
	mutable ARM_DateStripPtr	itsCalibDateStrip;
};

CC_END_NAMESPACE()

#endif
