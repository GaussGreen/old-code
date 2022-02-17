/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file callablesnowballcalculator.h
 *
 *  \base class for Callable snowball calculators
 *	\author  jp riaudel
 *	\version 1.0
 *	\date March 2005
 */

#ifndef _INGPCALCULATORS_CALLABALESNOWBALLCALCULATOR_H
#define _INGPCALCULATORS_CALLABALESNOWBALLCALCULATOR_H


#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/env.h"
#include "gencalculator.h"

#include "gpbase/countedptr.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/gpmatrix.h"

#include "gpinfra/curvemodelparam.h"

#include "gpmodels/typedef.h"



/// STL
CC_USING_NS(std,pair)

class ARM_CapFloor;
class ARM_Swap;
class ARM_SumOpt;



CC_BEGIN_NAMESPACE( ARM )

class ARM_SFRM;
class ARM_ModelParamsSFRM;


class ARM_CallableSnowBallCalculator : public ARM_GenCalculator 
{
private:
	static const string CallableSnowBallColNamesTable [];
	static const int Product1ToPriceColumns[];
	static const int Product2ToPriceColumns[];
	static const int AnaProductToPriceColumns[];

public:

    enum CallableSnowBallColAlias
    {
        EventDate=0,
		Index,
		SumOptStartDate,
		EndDate,
        CpnFwdStartDate,
        CpnFwdEndDate,
        CpnPayDate,
        CpnIT,
        FundingStartDate,
        FundingEndDate,
        NextCpnFwdStartDate,
        NextCpnFwdEndDate,
        NextCpnPayDate,
        NextCpnIT,
        PrevCpnFwdStartDate,
        PrevCpnFwdEndDate,
		ProductEndDate,
        DFPay,
		DFNum,
		DFNextPay,
		Nominal,
		NextNominal,
        Fees,
		Const,
		LevPrevCpn,
		LevNewOpt,
		StrikeOpt,
		MinCpn,
		MaxCpn,
		NextConst,
		NextLevPrevCpn,
		NextLevNewOpt,
		NextStrikeOpt,
		NextMinCpn,
		NextMaxCpn,
		PrevConst,
		PrevLevPrevCpn,
		PrevLevNewOpt,
		PrevStrikeOpt,
		PrevMinCpn,
		PrevMaxCpn,
		CpnIndex,
		Coupon,
		CouponSimple,
		CouponOption,
		StrikeSum,
		PaidCoupon,
		LBCoupon,
		LBCouponPaid,
		UBCoupon,
		UBCouponPaid,
		CouponOpt,
		CouponOptPaid,
		LBCouponAvg,
		LBCouponAvgPaid,
		AnaLBCouponAvg,
		AnaLBCoupon,
		CouponOptAvg,
		CouponOptAvgPaid,
		PrevCpnIndex,
		PrevCoupon,
		PrevCouponSimple,
		PrevCouponOpt,
		NextCpnIndex,
		NextSimpleCoupon,
		NextCpnIndexWithoutAdj,
		NextSimpleCouponWithoutAdj,
		StrikeCap,
		NextCouponCap,
		StrikeFloor,
		NextCouponFloor,
		NextCouponOpt,
		NextCoupon,
		FundingCoeff,
		FundingMargin,
		Funding,
		FundingPaid,
		CF,
		Option,
		Bermuda,
		CF1,
		CF2,
		AMCIndex,
		Option2,
		Bermuda2
    };

    enum mdmKeysAlias
    {
        YcKey=0,
        OswModelKey,
        CfModelKey,
        MrsKey,
		BetaKey,
		CorrelKey,
		HumpKey,
		BetaCorrelKey,
		ReCorrelKey,
        FundingKey,
		YcBasisDomKey,
        YcBasisFundKey,
		ForexKey,
			
		NbKeys,
    };

	enum productToPriceAlias
	{
		SnowBallSJPrice,
		FundingPrice,		
		SnowBallCallablePrice,
		NoOptSBPrice,
		StackySBPrice,
		CouponOptPrice,
		LBCouponAvgPrice,
		CouponOptAvgPrice,

		NbProductsToPrice
	};

	enum anaProductToPriceAlias
	{
		AnaLBCouponAvgPrice,
		AnaLBCouponPrice,
		NbAnaProductsToPrice
	};

	enum CalibMode
	{
		CAP,
		SWAPTION,
		SUMOPT,
		NOCALIB,
		CAPSWAPTION,

		NbCalibMode
	};

	enum ControlVariableMode
	{
		NO,
		SB,
		CSB,

		NbControlVariableMode
	};

	enum TriggerMode
	{
		TCoupon,
		TApproxCoupon
	};

	/// constructor, copy constructor, assignment constructor, destructor
	ARM_CallableSnowBallCalculator(const ARM_Date& startDate,
								   const ARM_Date& endDate,
								   int cpnFreq,
								   int cpnDaycount,
								   int cpnIndexTiming,
								   const string& cpnIdxTerm,
								   int cpnIndexDaycount,
								   const string& cpnResetCal,
								   const string& cpnPayCal,
								   int cpnIntRule,
								   int cpnIndexResetLag,
								   int payRec,
								   int CF,
								   const ARM_Curve& notionalProfile,
								   const ARM_Curve& fundnotionalProfile,
								   const ARM_Curve& constProfile,
								   const ARM_Curve& lPrevCpnProfile,
								   const ARM_Curve& lNewOptProfile,
								   const ARM_Curve& strikeOptProfile,
								   const ARM_Curve& minCpnProfile,
								   const ARM_Curve& maxCpnProfile,
								   int fundFreq,
								   int fundDaycount,
								   const ARM_Curve& fundingCoeffProfile,
								   const ARM_Curve& fundMarginProfile,
								   int NotifDays,
								   int NonCall,
								   const string& exerciseCal,
								   bool callSBOrStacky,
								   const ARM_Curve& feesProfile,
								   int NbPathBounding,
								   int NbPathPricing,
								   int NbMaxBucket,
								   int FixSteps,
								   const string& USMethod,
								   CalibMode calibMode,
								   bool betaCalib,
								   bool fixBoundary,
								   bool fixBeta,
								   ControlVariableMode conntrolVariableFlag,
								   bool fixControlVariable,
								   const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
								   const ARM_MarketData_ManagerRep& mktDataManager,
								   const ARM_StringVector& mdmKeys,
								   const string& genType1 = "Sobol",
								   const string& genType2 = "MRGK5",
								   int firstNbTimes = 0,
								   int firstNbDims = 0,
								   const string& PathScheme = "Incremental",
								   const string& PathOrder = "BucketOrder",
								   bool	antithetic = true,
								   ARM_ModelType modelType = ARM_PricingModelType::SFRM2F,
								   TriggerMode triggerMode = TCoupon,
								   int calibSwoptFreq=-1,
								   const string& regressors = "",
								   const vector< double >& hkData = vector< double >(0));

	ARM_CallableSnowBallCalculator(ARM_Currency& ccy,
								   const ARM_Date& startDate,
								   const ARM_Date& endDate,
								   int cpnFreq,
								   int cpnDaycount,
								   int cpnIndexTiming,
								   const string& cpnIdxTerm,
								   int cpnIndexDaycount,
								   const string& cpnResetCal,
								   const string& cpnPayCal,
								   int cpnIntRule,
								   int cpnIndexResetLag,
								   int payRec,
								   int CF,
								   const ARM_Curve& notionalProfile,
								   const ARM_Curve& constProfile,
								   const ARM_Curve& lPrevCpnProfile,
								   const ARM_Curve& lNewOptProfile,
								   const ARM_Curve& strikeOptProfile,
								   const ARM_Curve& minCpnProfile,
								   const ARM_Curve& maxCpnProfile,
								   int fundFreq,
								   int fundDaycount,
								   const ARM_Curve& fundingCoeffProfile,
								   const ARM_Curve& fundMarginProfile,
								   int NotifDays,
								   int NonCall,
								   const string& exerciseCal,
								   bool callSBOrStacky,
								   const ARM_Curve& feesProfile,
								   int NbPathBounding,
								   int NbPathPricing,
								   int NbMaxBucket,
								   int FixSteps,
								   const string& USMethod,
								   CalibMode calibMode,
								   bool betaCalib,
								   bool fixBoundary,
								   bool fixBeta,
								   ControlVariableMode conntrolVariableFlag,
								   bool fixControlVariable,
								   const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
								   const string& genType1 = "Sobol",
								   const string& genType2 = "MRGK5",
								   int firstNbTimes = 0,
								   int firstNbDims = 0,
								   const string& PathScheme = "Incremental",
								   const string& PathOrder = "BucketOrder",
								   bool	antithetic = true,
								   ARM_ModelType modelType = ARM_PricingModelType::SFRM2F,
								   TriggerMode triggerMode = TCoupon,
								   int calibSwoptFreq=-1,
								   const string& regressors = "",
								   const vector< double >& hkData = vector< double >(0));

	ARM_CallableSnowBallCalculator(const ARM_CallableSnowBallCalculator& rhs);
	ASSIGN_OPERATOR(ARM_CallableSnowBallCalculator)
	~ARM_CallableSnowBallCalculator();

	void InitCSBFromSummit(ARM_ZeroCurve* zcCpn,
						   ARM_VolCurve* swoptVC,
						   ARM_VolCurve* capVC,
						   ARM_VolLInterpol* capRho = NULL,
						   ARM_VolLInterpol* capNu = NULL,
						   ARM_VolLInterpol* capBeta = NULL,
						   ARM_VolLInterpol* swaptRho = NULL,
						   ARM_VolLInterpol* swaptNu = NULL,
						   ARM_VolLInterpol* swaptBeta = NULL,
						   double hump = -10000.0,
						   double betaCorrel = -10000.0,
						   double reCorrel = -10000.0,
						   int SABRSigmaOrAlpha = 1);


	/// Standard ARM support
	virtual ARM_Object* Clone() const;
	virtual void View(char* id = NULL, FILE* ficOut = NULL) const;

	/// ----------- virtual functions
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	// This function seems to be use less but I have some compiler crash and its the only way 
	// I found to fix it.
	// I think the MiddleRows function is really too long
	void MiddleRowsExt(int eventIdx, const ARM_DateStripCombiner& datesStructure, vector< string >& rowDescVec, vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;
	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	
	virtual void CreateAndSetModel();
	ARM_PricingModel* BuildModel_SFRM2F(ARM_ZeroCurve* pCurve);
	ARM_PricingModel* BuildModel_HK(ARM_ZeroCurve* pCurve);
	ARM_PricingModel* BuildModel_SBGM(ARM_ZeroCurve* pCurve);

	virtual void UpdateModel();
	void UpdateModel_SFRM2F();
	void UpdateModel_SBGM();
	void UpdateModel_HK();
		
	virtual void CreateAndSetCalibration();
	void CreateAndSetCalibration_SFRM2F();
	void CreateAndSetCalibration_HK();
	void CreateAndSetCalibration_SBGM();
	void CreateAndSetDensityCalibration(ARM_StdPortfolioPtr pfPtr);

	virtual void Calibrate();
	void Calibrate_SFRM2F();
	void Calibrate_SBGM();
	void Calibrate_HK();

	virtual void UpdateCalibration(bool isUpdateStrike=true);
	void UpdateCalibration_SFRM2F(bool isUpdateStrike=true);
	void UpdateCalibration_SBGM(bool isUpdateStrike=true);
	void UpdateCalibration_HK(bool isUpdateStrike=true);


	void CheckCalibMode();
    virtual double Price();
	
	virtual void ComputePricingData() const;
    virtual ARM_GenSecurity* GetSubGenSecurity(int indexType) const {return NULL;};

	void FilterCoupon(const ARM_GP_Vector& coupon, const ARM_GP_Vector& cpnIndex, double cpnMin, double cpnMax, ARM_GP_Vector& filterCoupon, ARM_GP_Matrix& filterCpnIndex );

	ARM_DateStrip* GetCalibSchedule() const;
	
	ARM_VanillaSecurityDensity* GetMarketDensity(
			double julianRD, double julianSD, double julianED,
			ARM_VolCurve* pAlpha, ARM_VolCurve* pBeta, ARM_VolCurve* pRho, ARM_VolCurve* pNu, 
			int sabrFlag, ARM_ZeroCurve* pCurve) const;


private:

	ARM_Date                    itsStartDate;           // start date
	ARM_Date                    itsEndDate;             // end date

    int                         itsPayRec;

    int                         itsCpnDayCount;         // coupon day count
    int                         itsCpnFreq;             // reset & pay frequency
    int                         itsCpnIndexTiming;      // fixing method (advance or arrears)
    string                      itsCpnIdxTerm;			// index term (only MM index)
    int                         itsCpnIndexDayCount;    // index day count (only MM index)
    string                      itsCpnResetCal;         // reset calendar
    string                      itsCpnPayCal;           // payment calendar
    int                         itsCpnResetLag;         // cpn reset lag
    int                         itsCpnIntRule;          // cpn intRule (adj or unadj)

	ARM_Curve					itsNotionalProfile;
	ARM_Curve					itsConstProfile;
	ARM_Curve					itsLPrevCpnProfile;
	ARM_Curve					itsLNewOptProfile;
	ARM_Curve					itsStrikeOptProfile;
	ARM_Curve					itsMinCpnProfile;
	ARM_Curve					itsMaxCpnProfile;

	int							itsCapOrFloor;			// 1 pour snowBall cap, -1 pour snowBall floor

	ARM_Curve					itsFundCoeffProfile;	// funding coeff curve
	ARM_Curve					itsFundNotionalProfile;	// funding coeff curve
    ARM_Curve					itsInitFundMarginProfile;   // funding spread curve
	ARM_Curve					itsFundMarginProfile;   // funding spread curve
    int                         itsFundFreq;            // funding frequency
    int                         itsFundDayCount;        // funding day count

    int                         itsNoticeDays;          // notification lag in days prior to reverse index reset
    string                      itsExerCal;             // notification calendar
	bool						itsCallSBOrStacky;		// Flag to call the SB or the Stacky

    ARM_Curve					itsFeesProfile;             // exercise fees curve

	int							itsNbNonCall;			// nb de periodes sans call

	ARM_GP_Vector				itsSBRank2Pos;
	ARM_GP_Vector				itsSBPos2Rank;
	ARM_GP_Vector				itsExerRank2Pos;
	ARM_GP_Vector				itsExerPos2Rank;

	/// Precomuted Vector
	std::deque<bool>			itsIsCpnFlow;
	ARM_GP_Vector				itsResetDatesVec;
	ARM_GP_Vector				itsStartDatesVec;
	ARM_GP_Vector				itsEndDatesVec;
	ARM_GP_Vector				itsPeriodsVec;
	ARM_GP_Vector				itsPaymentDatesVec;
	ARM_GP_Vector				itsConstVec;
	ARM_GP_Vector				itsLevPrevVec;
	ARM_GP_Vector				itsLevNewVec;
	ARM_GP_Matrix				itsSumOptionCoeffs;
	ARM_GP_Vector				itsStrikeOptVec;
	ARM_GP_Vector				itsCpnMinVec;
	ARM_GP_Vector				itsCpnMaxVec;
	ARM_GP_Vector				itsStrikeSum;
	ARM_GP_Vector				itsFeesVec;
	ARM_GP_Vector				itsvCpnNominal;

	ARM_GP_Vector				itsvInitialFundNominal;
	ARM_GP_Vector               itsvInitialFundSpread;


	int itsLastReset;

	ARM_DateStripCombiner		itsDateStructure;
	int							itsNbCoupons;

	/// MC Nb iterations
	int							itsNbPathPricing;
	int							itsNbPathBounding;
	int							itsNbMaxBucket;
	int							itsFixSteps;
	string						itsUSMethod;
	CalibMode					itsCalibMode;
	bool						itsBetaCalib;
	bool						itsFixBoundary;
	bool						itsFixBeta;
	ControlVariableMode			itsControlVariableFlag;
	bool						itsFixControlVariable;
	TriggerMode					itsTriggerMode;
	int							itsCalibSwoptFreq;
	string						itsRegressors;

	ARM_PricingModelPtr			itsSumOptModel;
	ARM_CalibMethodPtr			itsSumOptCalibMethod;

	string						itsGenType1;
	string						itsGenType2;
	int							itsFirstNbTimes;
	int							itsFirstNbDims;
	string						itsPathScheme;
	string						itsPathOrder;
	bool						itsAntithetic;

	ARM_ModelType				itsModelType;
	vector< double >			itsHkData;

	ARM_GP_Vector				itsFixBetaTimes;
	ARM_GP_Vector				itsFixBetaValues;

	ARM_DateStripPtr itsFundDateStrip;
	ARM_DateStripPtr itsStructDateStrip;

	/// Flag to specify product to price
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ itsProductsToPrice;
	bool itsHasBeenPriced;

	/// To memorise the line in the deal description corresponding to the first call
    /// Mutable because may be updated in const functions
    CC_IS_MUTABLE size_t itsFirstCallIdx;

	// Pricing Data
	double itsSnowBallSJPrice;
	double itsSnowBallSJStdDev;
	ARM_GP_VectorPtr itsSnowBallSJRA;
	double itsFundingPrice;
	double itsFundingStdDev;
	ARM_GP_VectorPtr itsFundingRA;
	double itsSnowBallCallablePrice;
	double itsSnowBallCallableStdDev;
	ARM_GP_VectorPtr itsSnowBallCallableRA;
	double itsNoOptSBPrice;
	double itsNoOptSBStdDev;
	ARM_GP_VectorPtr itsNoOptSBRA;
	double itsStackySBPrice;
	double itsStackySBStdDev;
	ARM_GP_VectorPtr itsStackySBRA;
	double itsCouponOptPrice;
	double itsCouponOptStdDev;
	ARM_GP_VectorPtr itsCouponOptRA;
	double itsLBCouponAvgPrice;
	double itsLBCouponAvgStdDev;
	ARM_GP_VectorPtr itsLBCouponAvgRA;
	double itsCouponOptAvgPrice;
	double itsCouponOptAvgStdDev;
	ARM_GP_VectorPtr itsCouponOptAvgRA;
	double itsAnaLBCouponAvgPrice;
	double itsAnaLBCouponPrice;
	ARM_GP_VectorPtr itsAnaLBCouponRA;

	ARM_VectorPtr itsCVBetas;
	/// Utilities
    ARM_INDEX_TYPE GetIndexType();

	// Return the reference model
	ARM_SFRM* GetRefModel();

	/// Specialised version for datas consistency
    virtual void CheckData();
	virtual void CheckMktData() {};

	ARM_DateStripCombiner InitDateStrip();    

	ARM_StdPortfolioPtr CreateCapPortfolio(bool atmOrStacky);
	ARM_StdPortfolioPtr CreateSwoptPortfolio();
	ARM_StdPortfolioPtr CreateSumOptionPortfolio();
	ARM_StdPortfolioPtr CreateBetaCapPortfolio();

	void UpdateBetaCapNominal(ARM_CalibMethod* calibMethod);
	void ComputePortPrices(bool isInitParam, CalibMode calibMode, ARM_CalibMethod* calibMethod );

	void FixBetaModelParam();

	ARM_CalibMethod* GetVolCalibMethod() const;
	ARM_CalibMethod* GetBetaCalibMethod() const;

	/// calibration part
	ARM_CalibMethod* CreateEmptyCalibration( CalibMode calibMode, bool betaCalib );
	void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;
	ARM_ModelParamsSFRM* CreateSFRMModelParams( CalibMode calibMode );

	void CreateSumOptionModel();
	void UpdateSumOptionStrike(double Strike);
	double ComputeSumOptionPrice(ARM_SumOpt* sumOpt);

	void ComputeSumOptCoeffs();
	ARM_StringVector CreateColumnNames();
	ARM_CstManagerPtr CreateCstManager();

	const int* ProductToPriceColumns() const;

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const	{ return  itsFundDateStrip;   };
	inline virtual ARM_DateStripPtr GetRefFundDateStrip()	const			{ return  itsFundDateStrip;     };
	inline virtual ARM_DateStripPtr GetRefDateStrip() const			{ return  itsStructDateStrip;     };

	inline virtual ARM_GP_Vector GetvCpnNominal() const				{ return  itsvCpnNominal;				};
	inline virtual ARM_GP_Vector GetvFundNominal() const			{ return  itsvInitialFundNominal;		};
	inline virtual ARM_GP_Vector GetvFundSpread() const				{ return  itsvInitialFundSpread;		};

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const			{ return dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));    };
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve()	const			{ return dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));     };
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));    };
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve()	const	{ return dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisFundKey]));     };
	inline virtual ARM_Forex* GetForex() const							{ return dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));     };

};

CC_END_NAMESPACE()

#endif