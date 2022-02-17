/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gencalculator.h
 *
 *  \base class for Caption calculators
 *	\author  y khlif
 *	\version 1.0
 *	\date November 2004
 */

#ifndef _INGPCALCULATORS_CAPTIONCALCULATOR_H
#define _INGPCALCULATORS_CAPTIONCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gencalculator.h"

#include "gpbase/countedptr.h"

/// kernel
#include <ccy/currency.h>

/// STL
CC_USING_NS(std,pair)

class ARM_Portfolio;
class ARM_Swaption;
class ARM_Swap;
class ARM_CapFloor;

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStrip;
class ARM_SFRM;

class ARM_CaptionCalculator : public ARM_GenCalculator
{
	static const string CaptionColNamesTable [];

public:
    enum CaptionColAlias
    {
        EventDate=0,
		PayDate,    
		FWDStartDate,
		FWDEndDate,
		IT,
		Strike,
		CpnSpread,
		Notional,
		FundStartDate,
		FundEndDate,
		FundSpread,	
		FundNotional,
		NotionalExchange,
		CapIndex,
		CpnNPVIntermediate,
		CpnNPV,
		FundFragmentNPV,
		FundNPV,
		Caplet,
		CapIntermediate,
		Cap,
		FeesPayDate,
		Fees,
		TotalNPV,
		ExerciseFees,
		BermudaProfile,
		Caption,
		CaptionStrike,
		NumeraireDate,
		ExerciseCondition,
		ExerciseConditionOfIndex,
		ProbaOfExercise
    };

	enum mdmKeysAlias
    {
        YcKey=0,        
        CfModelKey,
        OswModelKey,
		BasisKey,
        FundingKey,        
        ForexKey,
		MrsKey,
		BetaKey,
		CorrelKey,

        NbKeys
    };

	enum productToPriceAlias
	{
		CaptionPrice,
		CapPrice,		
		CouponLegPrice,
		FundingLegPrice,
		CaptionStrikes,
		ExerciseProbas,

		NbProductsToPrice
	};

	enum CalibrationMode
	{
		SWOPT_At_Exer,				
		AllSWOPTION,
		NOCalib,
		Calib
	};

private:
    
public:

	/// constructor, copy constructor, assignment constructor, destructor
	ARM_CaptionCalculator( const ARM_Date& startDate,
			const ARM_Date& endDate,
            const string& cpnIdxTerm,
			int payRec,
			int CF,
			const ARM_Curve& couponProfile,
			const string& FundIdxTerm,
			int NotifDays,
			int NonCall,
			const ARM_Curve& exerciseProfile,		
			const ARM_Curve& notionalProfile,
			int cpnDayCount,
			int cpnResetTiming,
			const string& cpnResetCal,
			const string& cpnPayCal,
			const ARM_Curve& cpnSpreadProfile,
			int fundDayCount,
			const string& fundResetCal,
			const string& fundPayCal,
			const ARM_Curve& fundSpreadProfile,
			int factorNb,
			int SFRMVolType,
			CalibrationMode SWOPTCalibMode,
			CalibrationMode BETACalibFlag,
			std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
			const ARM_MarketData_ManagerRep& mktDataManager,
			const ARM_StringVector& keys);
		
	ARM_CaptionCalculator( const ARM_Date& asOfDate,
			const ARM_Date& startDate,
			const ARM_Date& endDate,
			const string& cpnIdxTerm,
			int payRec,
			int CF,
			const ARM_Curve& couponProfile,
			const string& FundIdxTerm,
			int NotifDays,
			int NonCall,
			const ARM_Curve& exerciseProfile,		
			const ARM_Curve& notionalProfile,
			int cpnDayCount,
			int cpnResetTiming,
			const string& cpnResetCal,
			const string& cpnPayCal,
			const ARM_Curve& cpnSpreadProfile,
			int fundDayCount,
			const string& fundResetCal,
			const string& fundPayCal,
			const ARM_Curve& fundSpreadProfile);

	ARM_CaptionCalculator(const ARM_CaptionCalculator&);
	ARM_CaptionCalculator& operator=(const ARM_CaptionCalculator&);
	~ARM_CaptionCalculator();

// FIXMEFRED: mig.vc8 (23/05/2007 11:38:36): missing return type
	void InitCaptionForSummit(ARM_StringVector& mdmKeys,
						 ARM_MarketData_ManagerRep* mktDataManager,
						 int factorNb,
						 int SFRMVolType,
						 CalibrationMode SWOPTCalibMode,
						 CalibrationMode BETACalibFlag,
						 std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice);

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
	virtual void View(char* id = NULL, FILE* ficOut = NULL) const;

	/// accessors
	inline int	 GetNbSteps() const { return itsNbSteps; }
	inline void SetNbSteps ( int NbSteps ) { itsNbSteps = NbSteps; }

	inline int GetSFRMVolType() const {return itsSFRMVolType;}
    inline void SetSFRMVolType(int voltype) {itsSFRMVolType=voltype;}

	inline CalibrationMode GetSwoptCalibMode() const {return itsSwoptCalibMode;}
    inline void SetSwoptCalibMode(CalibrationMode mode) {itsSwoptCalibMode=mode;}

	inline CalibrationMode GetBetaCalibFlag() const {return itsBetaCalibFlag;}
    inline void SetBetaCalibFlag(CalibrationMode mode) {itsBetaCalibFlag=mode;}

	void SetCaptionToPrice(bool toPrice);
    bool IsCaptionToPrice() const;

	void SetCapToPrice(bool toPrice);
    bool IsCapToPrice() const;

	void SetCouponLegToPrice(bool toPrice);
    bool IsCouponLegToPrice() const;

	void SetFundingLegToPrice(bool toPrice);
    bool IsFundingLegToPrice() const;

	void SetCaptionStrikesToPrice(bool toPrice);
    bool IsCaptionStrikesToPrice() const;

	void SetExerciseProbasToPrice(bool toPrice);
    bool IsExerciseProbasToPrice() const;

	/// Get & Set MRS param for SUMMIT interface (deal dependent data)
    const ARM_ModelParam& GetMRS() const;
    void SetMRS(ARM_ModelParam* mrsParam);

	/// Get & Set Beta param for SUMMIT interface (deal dependent data)
    const ARM_ModelParam& GetBeta() const;
    void SetBeta(ARM_ModelParam* mrsParam);

	/// Get & Set Correlation param for SUMMIT interface (deal dependent data)
    const ARM_ModelParam& GetCorrel() const;
    void SetCorrel(ARM_ModelParam* correlationParam);

	ARM_CalibMethod* GetCFCalibMethod() const;
    const ARM_StdPortfolioPtr GetCFPortfolio() const;
    void SetCFPortfolio(const ARM_StdPortfolio& port);

	ARM_CalibMethod* GetSWOPTCalibMethod() const;
    const ARM_StdPortfolioPtr GetSWOPTPortfolio() const;

	/// ----------- virtual functions
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	virtual void CreateAndSetModel();
	virtual void UpdateModel();
	virtual void CreateAndSetCalibration();
	virtual void UpdateCalibration(bool isUpdateStrike=true);
    virtual double Price();
	virtual void Calibrate();
	virtual void ComputePricingData() const;
    virtual ARM_GenSecurity* GetSubGenSecurity(int indexType) const {return NULL;};
	
	/// not private because used in the interface
    void SetSWOPTPortfolio(const ARM_StdPortfolio& port);

private:
	/// Caption financial datas
	ARM_Date                    itsStartDate;           // start date
	ARM_Date                    itsEndDate;             // end date
	string						itsCpnIdxTerm;
	int                         itsPayRec;
	int                         itsCapOrFloor;
	ARM_Curve					itsCoupon;  
	string						itsFundIdxTerm;
	int							itsNotifDays;
	int							itsNonCall;
	ARM_Curve					itsExerStyle;
	ARM_Curve					itsNotional;
	int							itsCpnDayCount;
	int							itsCpnResetTiming;
	string                      itsCpnResetCal;         
    string                      itsCpnPayCal;           
	ARM_Curve					itsCpnSpread;
	int							itsFundDayCount;
	string                      itsFundResetCal;         
    string                      itsFundPayCal;
	ARM_Curve					itsFundSpread;
	int                         itsPayRecFund;
	int							itsFactorNb;
	int							itsCpnFreq;
	int							itsFundFreq;
	ARM_INDEX_TYPE				itsCpnIndexType;
	ARM_Curve					itsStrike;
			           
	// MC Nb iterations
	int itsNbSteps;

	/// Flag to specify product to price
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ itsProductsToPrice;
	bool itsHasBeenPriced;

    /// Flags to enable auto-calibrations
	int itsSFRMVolType;
	CalibrationMode itsSwoptCalibMode;		
	CalibrationMode itsBetaCalibFlag;		
    
	/// For Summit support : funding and basis currency (not known from the MdM
    /// because may not be already built)
    ARM_Currency   itsFundingCcy;
    ARM_Currency   itsBasisCcy;

	/// To memorise the line in the deal description corresponding to the first call
    /// Mutable because may be updated in const functions
    CC_IS_MUTABLE size_t itsFirstCallIdx;

	int itsNbdDeedResets;

	/// Names of model used in basis case (chosen among MdM keys)
    mdmKeysAlias itsCpnModelKey;
    mdmKeysAlias itsFundingModelKey;
    mdmKeysAlias itsBasisRefModelKey;

	//ARM_DateStripPtr principalDateStrip = datesStructure.GetDateStrip(Principal_SCHED);
	ARM_DateStripPtr itsCpnDateStrip;
	ARM_DateStripPtr itsFundDateStrip;
	ARM_DateStripPtr itsCalibPorfolioDateStrip;
	ARM_DateStripPtr itsPrincipalDateStrip;

	ARM_StdPortfolioPtr itsSwaptionPF;
	ARM_StdPortfolioPtr itsCapFloorPF;

	// Pricing Data
	double itsCaptionPrice;
	ARM_VectorPtr itsCapPrice;
	ARM_VectorPtr itsCouponLegPrice;
	ARM_VectorPtr itsFundingLegPrice;
	ARM_VectorPtr itsCaptionStrikes;
	double itsExerciseProbas;
	

	/// Set the keys or name alias of model (usefull in multi-currency context)
	void SetModelKeys();

	/// Specialised version for datas consistency
    virtual void CheckData();
	virtual void CheckMktData() {};

	/// Utilities
    ARM_INDEX_TYPE GetCpnIndexType();
	ARM_INDEX_TYPE GetFundIndexType();

	// Return the reference model
	ARM_SFRM* GetRefModel();
	
	void GenerateExerciseSyle(const ARM_Curve& exerciseProfile);
	void GenerateCapStrikes(void);
	void GenerateEquivalentFundSpreads(void);
	std::vector<double>& GenerateEquivalentSwoptStrikes(void);
	void ComputeEquivalentSwoptStrikes(std::vector<double>& equivStrikes);

	void CreateConvexityModel(bool withSmile);

	double GetSubPrice(int startRowIdx,int endRowIdx,CaptionColAlias columnName,
                       const string& evalDateStr,const ARM_DealDescription& dealDesc);
    	
	void CreateCaptionCapletFloorletAndDigital();

	void CreateDiagonalSwaption();
	
	void ComputeCaptionCapletFloorletPrices(
		bool isFreezeWeights, 
		bool isInitParam,
		bool isArrearCpn,
		CalibrationMode capFloorCalibMode,
		CalibrationMode digitalCalibMode, 
		bool computeSwaptionPrice);

	void ComputeDiagonalSwaptionPrice(
		bool isFreezeWeights,
		bool isArrearCpn);

	/// calibration part
	void CreateEmptyCalibration();
	void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  { return ARM_DateStripPtr(NULL);};
	inline virtual ARM_DateStripPtr GetRefFundDateStrip() const { return ARM_DateStripPtr(NULL);};
	inline virtual ARM_DateStripPtr GetRefDateStrip() const { return ARM_DateStripPtr(NULL);};

	/// Vector
	inline virtual std::vector<double> GetvCpnNominal() const  { return std::vector<double>(0);};
	inline virtual std::vector<double> GetvFundNominal() const { return std::vector<double>(0);};
	inline virtual std::vector<double> GetvFundSpread() const { return std::vector<double>(0);};

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve() const  { return NULL;};
	inline virtual ARM_Forex* GetForex() const  { return NULL;};
};

CC_END_NAMESPACE()

#endif