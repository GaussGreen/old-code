/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file GlobalCapcalculator.h
 *
 *  \brief
 *
 *	\author  P. LAM & H. ESSEFI
 *	\version 1.0
 *	\date feb 2006
 */

#ifndef _INGPCALCULATORS_GLOBALCAPCALCULATOR_H
#define _INGPCALCULATORS_GLOBALCAPCALCULATOR_H

#include "gencalculator.h"
#include "gpbase/countedptr.h"

#include "volcube.h"

/// kernel
#include "util/refvalue.h"
#include "util/exercise.h"
#include "ccy/currency.h"
#include "inst/portfolio.h"
#include "mod/bssmiled.h"
#include "inst/globalcap.h"

/// STL
CC_USING_NS(std,pair)

/// forward declaration
class ARM_Portfolio;
class ARM_BSSmiledModel;

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStrip;
class ARM_MCMethod;

class ARM_GlobalCapCalculator : public ARM_GenCalculator
{
protected:

	static const string GlobalCapColNamesTable [];
	static const int GlobalCapProductToPriceColumns [];

public:
	
	enum GlobalCapColAlias
	{
		ResetDate = 0,
		StartDate,
		EndDate,
		PayDate,
		Notional,
		Spread,
		Fixed,
		CapLev,
		FundLev, 
		Strike,
		Basis,
		Index,
		DF,
		FundRate,
		CouponRate,
		Funding,
		SumCoup,
		Coupon,
		Swap,
		GlobalCap,
		Product,
	};

	enum GlobalCapProductToPriceAlias
	{
		ProductPrice,
		FundingPrice,
		CouponPrice,
		SwapPrice,
		GlobalCapPrice,
		IndexPrice,
		FundRatePrice,
		SumCoupPrice,
		DFPrice,

		NbProductsToPrice,
	};

	enum mdmKeysAlias
    {
        YcKey =  0,
		CfModelKey,
        MrsKey,
        
        NbKeys
    };

	//Constructors
	ARM_GlobalCapCalculator(const ARM_MarketData_ManagerRep& mktDataManager );

	ARM_GlobalCapCalculator(const ARM_GlobalCapCalculator&);
	
	ARM_GlobalCapCalculator& operator=(const ARM_GlobalCapCalculator&);

	ARM_GlobalCapCalculator(  ARM_Currency& ccy,
							  ARM_Date&	startDate,
							  ARM_Date&	endDate,
							  int payRec,
							  ARM_ReferenceValue& notional,
							  ARM_INDEX_TYPE indexTerm,
							  int dayCount,
							  int payFreq,
							  int resetGap,
							  int payGap,
							  int resetTiming,
							  int payTiming,
							  string resetCal,
							  string payCal,
							  int adjRule,
							  int intRule,
							  ARM_ReferenceValue& fundLev,
							  ARM_ReferenceValue& capLev,
							  ARM_ReferenceValue& capFixed,
							  ARM_ReferenceValue& capStrike,
							  ARM_ReferenceValue& capSpread,
							  ARM_Vector* globalCapParams,
							  ARM_ReferenceValue* pastFixings,
							  int nbSteps,
							  vector<string> randGenerator,
							  int samplerType,
							  ARM_Vector* calibParams,
							  ARM_StringVector& mdmKeys,
							  const ARM_MarketData_ManagerRep& mktDataManager,
							  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice);

	ARM_GlobalCapCalculator(  ARM_GlobalCap* globalCap,
							  ARM_ReferenceValue& fundLev,
							  ARM_ReferenceValue& capLev,
							  int nbSteps,
							  vector<string> randGenerator,
							  int samplerType,
							  ARM_Vector* calibParams,
							  ARM_StringVector& mdmKeys,
							  const ARM_MarketData_ManagerRep& mktDataManager,
							  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice);
	
	ARM_GlobalCapCalculator(  ARM_Date& asOfDate,
							  ARM_GlobalCap* globalCap,
							  ARM_ReferenceValue& fundLev,
							  ARM_ReferenceValue& capLev);

	//Destructor
	~ARM_GlobalCapCalculator();

	virtual ARM_Object*			Clone() const;
	
	virtual ARM_CstManagerPtr	CreateCstManager();

	ARM_StringVector			PricedColumnNames() const;

	//Accessors
	inline  const ARM_Currency&		GetCcy() const {return itsCcy;};
    inline  const ARM_Date&			GetStartDate() const { return itsStartDate;}
	inline  const ARM_Date&			GetEndDate() const { return itsEndDate;}   
	inline  const int				GetPayRec() const {return itsPayRec;}

	inline  int						GetDayCount() const {return itsDayCount;}
	inline  int						GetFreq() const {return itsPayFreq;}
	inline  int						GetResetGap() const {return itsResetGap;}
	inline  int						GetPayGap() const {return itsPayGap;}
	inline  int						GetPayTiming() const {return itsPayTiming;}
	inline  string					GetPayCal() const {return itsPayCal;}
	inline  int						GetAdjRule() const {return itsAdjRule;}
	inline	int						GetIntRule() const {return itsIntRule;}

	//Get Step-up vectors
	virtual inline ARM_ReferenceValue&		GetNotional(){return itsNotional;}
	virtual inline ARM_ReferenceValue&		GetFundLev(){return itsFundLev;}
	virtual inline ARM_ReferenceValue&		GetCapStrike(){return itsCapStrike;}
	virtual inline ARM_ReferenceValue&		GetCapSpread(){return itsCapSpread;}
	virtual inline ARM_ReferenceValue&		GetCapFixed(){return itsCapFixed;}
	virtual inline ARM_ReferenceValue&		GetCapLev(){return itsCapLev;}

	//Mutators
	void SetCcy(const ARM_Currency& ccy){itsCcy=ccy;}
	void SetStartDate(ARM_Date& startDate){itsStartDate = startDate;}
	void SetEndDate(ARM_Date& endDate){itsEndDate = endDate;}
	void SetPayRec(int payRec){itsPayRec = payRec;}
	void SetNotional(ARM_ReferenceValue& notional){itsNotional = notional;}
	
	void SetIndexType(ARM_INDEX_TYPE indexType) {itsIndexType = indexType;}
	void SetDayCount(int dayCount) {itsDayCount = dayCount;}
	void SetPayFreq(int payFreq) {itsPayFreq = payFreq;}
	void SetResetGap(int resetGap) {itsResetGap=resetGap;}
	void SetPayGap(int payGap) {itsPayGap=payGap;}
	void SetResetTiming(int resetTiming) {itsResetTiming = resetTiming;}
	void SetPayTiming(int payTiming) {itsPayTiming = payTiming;}
	void SetResetCal(string resetCal) {itsResetCal = resetCal;}
	void SetPayCal(string payCal) {itsPayCal = payCal;}
	void SetAdjRule(int adjRule) {itsAdjRule = adjRule;}
	void SetIntRule(int intRule) {itsIntRule = intRule;}
	void SetFundLev(ARM_ReferenceValue& fundLev) {itsFundLev = fundLev;}

	void SetCapStrike(ARM_ReferenceValue& capStrike){itsCapStrike = capStrike;}
	void SetCapSpread(ARM_ReferenceValue& capSpread){itsCapSpread = capSpread;}
	void SetCapFixed(ARM_ReferenceValue& capFixed){itsCapFixed = capFixed;}
	void SetCapLev(ARM_ReferenceValue& capLev){itsCapLev = capLev;}
 
	void SetCalibParams(ARM_Vector* calibParams){itsCalibParams = (ARM_Vector*) calibParams->Clone();}
	void SetGlobalCapParams(ARM_Vector* globalCapParams){itsGlobalCapParams = (ARM_Vector*) globalCapParams->Clone();}
	void SetPastFixings(ARM_ReferenceValue* pastFixings)
	{
		if (pastFixings)
		{
			itsPastFixings = (ARM_ReferenceValue*) pastFixings->Clone();
		}
		else
		{
			itsPastFixings = NULL;
		}
	}


	void InitGlobalCapFromSummit( ARM_ZeroCurve* zc, 
								  ARM_VolCurve* capVol,
								  ARM_VolCurve* rhoCap,
								  ARM_VolCurve* nuCap,
								  ARM_VolCurve* betaCap,
								  double mrs,
								  ARM_Vector* calibParams,
								  int nbSteps,
								  vector<string> randGenerator,
								  int samplerType,
								  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice);

	//Outputs	
	virtual ARM_RowInfo				ColumnNames() const;
	virtual ARM_RowInfo				MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	virtual ARM_DateStripCombiner	DatesStructure() const;
	virtual ARM_DateStripCombiner	CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	/// pricing function
	virtual void					ComputePricingData() const;
	virtual void					CreateAndSetModel();
	virtual void					CreateAndSetCalibration();
    virtual void					UpdateModel();
    virtual void					UpdateCalibration(bool isUpdateStrike=true);
	virtual void					Calibrate();
    virtual double					Price();
    virtual void					CheckData(); /// Check internal data consistency
	virtual void					CheckMktData(); // Check mkt data consistency
	/// Initialisation to 0 of all columns of the deal description that could be priced
    virtual void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	// Some utilities
	virtual string					ExportShortName() const { return "LCGCF";}
	string							DealDesDataToString(const string& indent = "", const string& nextIndent = "") const;
	virtual string					GeneralDataToString(const string& indent = "", const string& nextIndent = "") const;
	virtual void					View(char* id = NULL, FILE* ficOut = NULL) const;
	void							GenerateProductDescription(ARM_Portfolio* portfolio);

	void CleanUp();
    void CopyNoCleanUp(const ARM_GlobalCapCalculator& rhs);

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
	int							itsPayRec;
	ARM_ReferenceValue			itsNotional;

	// Funding Leg / Swap Var leg
	ARM_INDEX_TYPE				itsIndexType;
	int							itsDayCount;
	int							itsPayFreq;
	int							itsResetGap;
	int							itsPayGap;
	int							itsResetTiming;
	int							itsPayTiming;
	string						itsResetCal;
	string						itsPayCal;
	int							itsAdjRule;
	int							itsIntRule;
	ARM_ReferenceValue			itsFundLev;

	// Cap Data
	ARM_ReferenceValue			itsCapStrike;
	ARM_ReferenceValue			itsCapSpread;
	ARM_ReferenceValue			itsCapFixed;
	ARM_ReferenceValue			itsCapLev;

	// Model Params
	int							itsNbSteps;
	vector<string>				itsRandGenerator;
	int							itsSamplerType;

	// Calib Params
	ARM_Vector*					itsCalibParams;

	ARM_NumMethodPtr			itsMCMethod;

	ARM_Vector*					itsGlobalCapParams;
	ARM_ReferenceValue*			itsPastFixings;

	// Pricing Data
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/				itsProductsToPrice;

	bool						itsHasBeenPriced;

	double						itsFundingPrice;
	double						itsGlobalCapPrice;
	double						itsCouponPrice;
	double						itsSwapPrice;
	double						itsProductPrice;
	double						itsSumCoupPrice;
	double						itsDFPrice;
	double						itsFundRatePrice;
	double						itsIndexPrice;

	// Date Strip
	mutable ARM_DateStripPtr	itsExerDateStrip;
};

CC_END_NAMESPACE()

#endif
