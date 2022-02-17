/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *      \file tarncalculator.h
 *
 *  \brief
 *
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPCALCULATORS_TARNCALCULATOE_H
#define _INGPCALCULATORS_TARNCALCULATOE_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gencalculator.h"
#include "gpbase/assignop.h"

#include "gpbase/countedptr.h"

#include "gpmodels/typedef.h"

/// kernel
#include <ccy/currency.h>

/// STL
CC_USING_NS(std,pair)

class ARM_Portfolio;
class ARM_Swaption;
class ARM_Swap;
class ARM_CapFloor;

/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string YC_BASIS_KEY_NAME      = "YC_BASIS_";
const string OSWMODEL_KEY_NAME      = "OSWMOD_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string MRS_KEY_NAME           = "MRS_";
const string BETA_KEY_NAME          = "BETA_";
const string CORREL_KEY_NAME        = "CORREL_";
const string FOREX_KEY_NAME         = "FOREX_";
const string HUMP_KEY_NAME          = "HUMP_";
const string BETACORREL_KEY_NAME    = "BETACORREL_";


CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStrip;
class ARM_SFRM;

class ARM_TARNCalculator : public ARM_GenCalculator
{
protected:
	static const string TARNColNamesTable [];
	static const int ProductToPriceColumns[];


public:
	enum productToPriceAlias
	{
		TARNPrice,
		SwapPrice,
		LifeTimeCapPrice,
		LifeTimeFloorPrice,
		DigitalFundingPrice,
		DigitalPrice,
		FundingPrice,
		ExerciseStrikesPrice,
		ExerciseProbasPrice,
		ExerciseTimeAverage,
		LastLifeTimeCapPrice,
		NbProductsToPrice
	};


    enum TARNColAlias
    {
        EventDate=0,
		NextEventDate,
        StartDate,
        EndDate,
        PayDate,
        IT,
        IndexStartDate,
		IndexEndDate,
		ITStdSwap,
        Strike,
        Leverage,
		CouponMin,
		CouponMax,
		LevPrev,
		LifeTimeCapTarget,
		LifeTimeFloorTarget,
        Nominal,
		Fees,
        FundingSpread,
        FundingStartDate,
        FundingEndDate,
        FundingAnnuity,
        FundingVarFlow,
        FundingSpreadFlow,
		Funding,
		FundingPaid,
		DFPay,
		CpnIndex,
		CouponWithoutIT,
		Coupon,
		PaidNb,
		CouponCF,
		PaidCoupon,
		SumCoupons,
		DigitalUp,
		DigitalDown,
		IsAlived,
		HasTriggered,
		EstimExerciseStrike,
		StdFixFlow,
		StdFixRFFlow,
		Swap,
		TargetCap,
		TargetFloor,
		ControlVariable,
		ControlVariable2,
		ExerciseStrike,
		ExerciseProba,
		ExerciseTime,
		LifeTimeCap,
		LastLifeTimeCap,
		LifeTimeFloor,
		DigitalFunding,
		Digital,
        TARN
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
        BasisKey,
		FundBasisKey,
        ForexKey,
		
        NbKeys
    };

	enum CalibrationMode
	{
		NOCalib,
		RFStrike,
		ExerStrike,
		ATMStrike,
		BlackShift,
		GaussShift
	};

private:
    
public:

	/// constructor, copy constructor, assignment constructor, destructor
	ARM_TARNCalculator(const ARM_Date& asOfDate,
		const ARM_Date& startDate,
        const ARM_Date& endDate,
        const ARM_Curve& strike,
        int payRec,
        int cpnDayCount,
        int cpnFreq,
        int cpnTiming,
        const string& cpnIndexTerm,
        int cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
        int cpnResetGap,
		int intRule,
        const ARM_Curve& leverage,
		const ARM_Curve& couponMin,
		const ARM_Curve& couponMax,
		double lifeTimeCap,
		bool globalCapFlag,
        double lifeTimeFloor,
        const ARM_Curve& fundSpread,
        int fundFreq,
        int fundDayCount,
        const ARM_Curve& nominal,
		const ARM_Curve& fees,
		const std::vector<double>& nbIterations,
		CalibrationMode capFloorCalibMode,
		CalibrationMode digitalCalibMode,
		bool oswCalibFlag,
		bool controlVariableFlag,
		bool digitalSmoothing,
		bool smiledFRMRescalling,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
		const ARM_Currency& cpnCcy,
		const ARM_Currency& fundingCcy,
		const ARM_Currency& basisCcy,
	    const ARM_Currency& domesticCcy,
		const ARM_Currency& foreignCcy,
		const ARM_Curve& fundnominal,
		const ARM_StringVector *mdmKeys = NULL,
		const ARM_MarketData_ManagerRep *mktDataManager = NULL,
		//Used for summit specific deals
		bool isCustomResetFlag = false,
		std::vector<double>& customResetDates = std::vector<double>(0));

		ARM_TARNCalculator(const ARM_Date& asOfDate,
		const ARM_Date& startDate,
        const ARM_Date& endDate,
        const ARM_Curve& strike,
        int payRec,
        int cpnDayCount,
        int cpnFreq,
        int cpnTiming,
        const string& cpnIndexTerm,
        int cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
        int cpnResetGap,
		int intRule,
        const ARM_Curve& leverage,
		const ARM_Curve& couponMin,
		const ARM_Curve& couponMax,
		double lifeTimeCap,
		bool globalCapFlag,
        double lifeTimeFloor,
        const ARM_Curve& fundSpread,
        int fundFreq,
        int fundDayCount,
        const ARM_Curve& nominal,
		const ARM_Curve& fees,
		const std::vector<double>& nbIterations,
		ARM_ModelType modelType,
		CalibrationMode capFloorCalibMode,
		CalibrationMode digitalCalibMode,
		bool digitalSmoothing,
		bool oswCalibFlag,
		bool controlVariableFlag,
		bool smiledFRMRescalling,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		bool antithetic,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
		const ARM_Currency& cpnCcy,
		const ARM_Currency& fundingCcy,
		const ARM_Currency& basisCcy,
	    const ARM_Currency& domesticCcy,
		const ARM_Currency& foreignCcy,
		const ARM_Curve& fundnominal,
		const ARM_StringVector *mdmKeys = NULL,
		const ARM_MarketData_ManagerRep *mktDataManager = NULL,
		//Used for summit specific deals
		bool isCustomResetFlag = false,
		std::vector<double>& customResetDates = std::vector<double>(0));

	/// constructor for SnowBall
	ARM_TARNCalculator(const ARM_Date& startDate,
        const ARM_Date& endDate,
        const ARM_Curve& strike,
        int payRec,
        int cpnDayCount,
        int cpnFreq,
        int cpnTiming,
        const string& cpnIndexTerm,
        int cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
        int cpnResetGap,
		int intRule,
		int stubRule,
        const ARM_Curve& leverage,
		const ARM_Curve& couponMin,
		const ARM_Curve& couponMax,
		const ARM_Curve& levPrev,
		double lifeTimeCap,
		bool globalCapFlag,
        double lifeTimeFloor,
        const ARM_Curve& fundSpread,
        int fundFreq,
        int fundDayCount,
        const ARM_Curve& nominal,
		const ARM_Curve& fees,
		const ARM_Curve& fundNominal,
		const std::vector<double>& nbIterations,
		ARM_ModelType modelType,
		CalibrationMode capFloorCalibMode,
		CalibrationMode digitalCalibMode,
		bool oswCalibFlag,
		bool controlVariableFlag,
		bool digitalSmoothing,
		bool smiledFRMRescalling,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		bool antithetic,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
        const ARM_MarketData_ManagerRep& mktDataManager);

	ARM_TARNCalculator(const ARM_TARNCalculator&);
	ASSIGN_OPERATOR(ARM_TARNCalculator)
	~ARM_TARNCalculator();

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
	virtual void View(char* id = NULL, FILE* ficOut = NULL) const;

	void InitTARNForSummit(ARM_ZeroCurve    *zcCpn, 
		                   ARM_VolCurve     *swoptVC, 
						   ARM_VolCurve     *capVC, 
						   ARM_VolLInterpol *capRo       = NULL, 
						   ARM_VolLInterpol *capNu       = NULL,
		                   ARM_VolLInterpol *capBeta     = NULL,
						   ARM_VolLInterpol *swoptRo     = NULL, 
						   ARM_VolLInterpol *swoptNu     = NULL,
		                   ARM_VolLInterpol *swoptBeta   = NULL,
		         		   ARM_ZeroCurve    *zcFund      = NULL, 
		                   ARM_ZeroCurve    *zcCpnBasis  = NULL, 
		                   ARM_ZeroCurve    *zcFundBasis = NULL, 
						   double fxSpot				 = 0.0,
						   int hedgeUpdate				 = 0,
                           ARM_ModelType modelType       = ARM_PricingModelType::SFRM2F,
                           double betaCorrel             = -10000.0,
                           double hump                   = -10000.0,
                           int SABRSigmaOrAlpha = 1);

	std::vector<double> GetNbIterations() const { return itsNbIterations; }
	void SetNbIterations( int nbIterations ) { itsNbIterations[0] = nbIterations; itsNbIterations[1] = nbIterations; }

	CalibrationMode GetCapFloorCalibMode() const {return itsCapFloorCalibMode;}
    void SetCapFloorCalibMode(CalibrationMode mode) {itsCapFloorCalibMode=mode;}

	CalibrationMode GetDigitalCalibMode() const {return itsDigitalCalibMode;}
    void SetDigitalCalibMode(CalibrationMode mode) {itsDigitalCalibMode=mode;}

	bool GetOSWCalibFlag() const {return itsOSWCalibFlag;}
    void SetOSWCalibFlag(bool flag) {itsOSWCalibFlag=flag;}

    void SetTARNToPrice(bool toPrice);
    bool IsTARNToPrice() const;

	void SetSwapToPrice(bool toPrice);
    bool IsSwapToPrice() const;

	void SetLifeTimeCapToPrice(bool toPrice);
    bool IsLifeTimeCapToPrice() const;

	void SetLifeTimeFloorToPrice(bool toPrice);
    bool IsLifeTimeFloorToPrice() const;

	void SetDigitalFundingToPrice(bool toPrice);
    bool IsDigitalFundingToPrice() const;

	void SetDigitalToPrice(bool toPrice);
    bool IsDigitalToPrice() const;

	void SetFundingToPrice(bool toPrice);
    bool IsFundingToPrice() const;

	void SetExerciseStrikesToPrice(bool toPrice);
    bool IsExerciseStrikesToPrice() const;

	void SetExercisesProbasToPrice(bool toPrice);
    bool IsExercisesProbasToPrice() const;

	void SetExercisesTimeAverageToPrice(bool toPrice);
    bool IsExercisesTimeAverageToPrice() const;

	void SetHasToPrice(bool toReset);
	bool IsResetPricing() const;

	void SetProductToPrice(vector<string>& flags);
	
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	ARM_StringVector CreateColumnNames();

	virtual string CreateExStrikeDescription(bool isFirstExer, const string& cpnModelName) const;
	virtual double GetC0() const {return 0.0;};
	virtual double GetCoefReverse() const {return 1.0;};

	const string& GetCpnIndexTerm() const {return itsCpnIndexTerm;};
	void SetCpnIndexTerm(const string& idxTerm){itsCpnIndexTerm=idxTerm;};

	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	void SetCustomResetFlag(bool flag){itsCustomResetFlag = flag;}
	const bool GetCustomResetFlag() const{return itsCustomResetFlag;}
	void SetCustomResetDates(std::vector<double> customResetDates){itsCustomResetDates = customResetDates;}
	std::vector<double> const GetCustomResetDates()const{return itsCustomResetDates;}
	
	// virtual for SNOWBALL
	virtual string CreateCpnDescription(bool isFirstExer) const;
	virtual string CreateCFCpnDdescription( bool isFirstExer, const string& cpnModelName) const;
	virtual void SetCVSwapMiddleText( ARM_DealDescriptionPtr& dealDesc, size_t knockOutIndex ) const;
	
	/// Get & Set MDM param for SUMMIT interface (deal dependent data)
    const ARM_ModelParam& GetMRS() const;
    void SetMRS(ARM_ModelParam* mrsParam);
    const ARM_ModelParam& GetBeta() const;
    void SetBeta(ARM_ModelParam* mrsParam);
    const ARM_ModelParam& GetCorrel() const;
    void SetCorrel(ARM_ModelParam* correlationParam);

    // For SBGM
    void SetBetaCorrel(ARM_ModelParam* betaCorrelParam);
    void SetHump(ARM_ModelParam* humpParam);

	virtual void CreateAndSetModel();
	ARM_PricingModel* GetModel_SFRM2F(ARM_ZeroCurve* pCurve);
	ARM_PricingModel* GetModel_HK(ARM_ZeroCurve* pCurve);
	ARM_PricingModel* GetModel_SBGM(ARM_ZeroCurve* pCurve);
	
	// Calibration
	ARM_CalibMethod* GetCFCalibMethod(bool calibSwaption) const;
    const ARM_StdPortfolioPtr GetCFPortfolio() const;
    void SetCFPortfolio(const ARM_StdPortfolio& port);

	ARM_CalibMethod* GetOSWCalibMethod() const;
    const ARM_StdPortfolioPtr GetOSWPortfolio() const;
    void SetOSWPortfolio(const ARM_StdPortfolio& port);

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

	bool	FindIdxInVector( const std::vector<double>& timeVector, double time, double daysNb, size_t& idx) const;
	double  FindValueInVector( const std::vector<double>& timeVector, double time, double daysNb ) const;

    virtual double Price();
	virtual void ComputePricingData() const;
    virtual ARM_GenSecurity* GetSubGenSecurity(int indexType) const {return NULL;};

	/// Specialised version for datas consistency
    virtual void CheckData();
	virtual void CheckMktData();

	/// Set the keys or name alias of model (usefull in multi-currency context)
	void SetModelKeys();

    void ARM_TARNCalculator::CleanAfterPricing(void);

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const	{ return itsFundDateStrip;};
	inline virtual ARM_DateStripPtr GetRefFundDateStrip() const			{ return itsFundDateStrip;};
	inline virtual ARM_DateStripPtr GetRefDateStrip() const			{ return itsStructDateStrip;};

	/// Vector
	inline virtual std::vector<double> GetvCpnNominal() const  { return itvCpnsNominal;};
	inline virtual std::vector<double> GetvFundNominal() const { return itsvFundNominal;};
	inline virtual std::vector<double> GetvFundSpread() const { return itsvInitialFundSpread;};

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const			{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));    };
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve()	const			{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));     };
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));    };
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve()	const	{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundBasisKey]));     };
	inline virtual ARM_Forex* GetForex() const							{ return dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));     };


private:

	/// TARN financial datas
	ARM_Date                    itsStartDate;           // start date
	ARM_Date                    itsEndDate;             // end date

    ARM_Curve					itsNominal;             // nominal curve
	ARM_Curve					itsFees;				// fees curve
    ARM_Curve					itsStrike;              // fix coupon rate or reverse coupon strike cst...
    int                         itsPayRec;

    ARM_Date                    itsFixEndDate;          // fix coupon leg end date
    int                         itsFixDayCount;         // fix coupon day count

    int                         itsCpnDayCount;         // reverse coupon day count
    int                         itsCpnFreq;             // reset & pay frequency (fix or reverse coupon)
    int                         itsCpnTiming;           // reverse fixing method (advance or arrears)
    string                      itsCpnIndexTerm;        // reverse index term (only MM index)
    int                         itsCpnIndexDayCount;    // reverse index day count (only MM index)
    string                      itsCpnResetCal;         // reverse reset calendar
    string                      itsCpnPayCal;           // reverse payment calendar
    int                         itsCpnResetGap;         // reverse index reset gap
	int							itsIntRule;				// reverse floater int rule
	int							itsStubRule;            // ability to have shortstart .... 	
    ARM_Curve					itsLeverage;            // reverse index leverage curve
	ARM_Curve					itsCpnMin;				// coupon min
	ARM_Curve					itsCpnMax;				// coupon max
	ARM_Curve					itsLevPrev;				// leverage coeff curve for the previous coupon
	double						itsLifeTimeCapTarget;   // TARN life time cap target
	bool						itsGlobalCapFlag;		// global cap flag
    double						itsLifeTimeFloorTarget; // TARN life time floor target

    ARM_Curve					itsFundSpread;          // funding spread curve
    int                         itsFundFreq;            // funding frequency
    int                         itsFundDayCount;        // funding day count

    ARM_Curve					itsFundNominal;         // funding nominal curve

	//Summit specific deals
	bool						itsCustomResetFlag;     //true in case of custom reset dates
	std::vector<double>				itsCustomResetDates;   
	
	ARM_DateStripPtr itsFundDateStrip;
	ARM_DateStripPtr itsStructDateStrip;

	std::vector<double>   itvCpnsNominal;
	std::vector<double>	itsvFundNominal;
	std::vector<double>   itsvInitialFundSpread;

	// MC Parameters
	string itsGenType1;
	string itsGenType2;
	int itsFirstNbTimes;
	int itsFirstNbDims;
	string itsPathScheme;
	string itsPathOrder;
	bool itsAntithetic;

	// MC Nb iterations
	std::vector<double>				itsNbIterations;
	ARM_ModelType				itsModelType;

protected:
	/// Flag to specify product to price
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ itsProductsToPrice;

private:
	bool itsHasBeenPriced;

    /// Flags to enable auto-calibrations
	CalibrationMode itsCapFloorCalibMode;		/// Volatility bootstrapping
	CalibrationMode itsDigitalCalibMode;		/// Beta bootsrapping
    bool itsOSWCalibFlag;						/// Mean reversion optimization
	bool itsControlVariableFlag;				/// Control Variable flag
	bool itsSmiledFRMRescallingFlag;			/// Smiled FRM Rescalling
	ARM_VectorPtr itsCVRefPrices;
	bool itsDigitalSmoothingFlag;				/// DigitalSmoothing

	// min, max or ref values for calibration (mean rev, beta and correl)
	double minMeanRev;
	double maxMeanRev;
	double minBeta;
	double maxBeta;
	double correl;

	/// To memorise the line in the deal description corresponding to the first call
    /// Mutable because may be updated in const functions
    CC_IS_MUTABLE size_t itsFirstEventIdx;
	CC_IS_MUTABLE size_t itsFirstRFIdx;

protected:
	/// Names of model used in basis case (chosen among MdM keys)
    mdmKeysAlias itsCpnModelKey;
    mdmKeysAlias itsFundingModelKey;
    mdmKeysAlias itsBasisRefModelKey;

private:
	 /// Vector of underlying diagonal swaps (used for equivalent strike & vol computations)
    vector< pair<ARM_Swap*,int> >         itsStdSwaps;
	vector< pair<ARM_CapFloor*,int> >	  itsCaplets;

	// Pricing Data
	double itsTARNPrice;
	double itsTARNStdDev;
	ARM_GP_VectorPtr itsTARNRA;
	double itsSwapPrice;
	double itsSwapStdDev;
	ARM_GP_VectorPtr itsSwapRA;
	double itsLifeTimeCapPrice;
	double itsLifeTimeCapStdDev;
	ARM_GP_VectorPtr itsLifeTimeCapRA;
	double itsLastLifeTimeCapPrice;
	double itsLastLifeTimeCapStdDev;
	ARM_GP_VectorPtr itsLastLifeTimeCapRA;
	double itsLifeTimeFloorPrice;
	double itsLifeTimeFloorStdDev;
	ARM_GP_VectorPtr itsLifeTimeFloorRA;
	double itsDigitalFundingPrice;
	double itsDigitalFundingStdDev;
	ARM_GP_VectorPtr itsDigitalFundingRA;
	double itsDigitalPrice;
	double itsDigitalStdDev;
	ARM_GP_VectorPtr itsDigitalRA;
	double itsFundingPrice;
	double itsFundingStdDev;
	ARM_GP_VectorPtr itsFundingRA;
	ARM_VectorPtr itsExerciseStrikes;
	ARM_VectorPtr itsExerciseProbas;
	double itsExerciseTimeAverage;

	/// Utilities
    ARM_INDEX_TYPE GetIndexType();

	// Return the reference model
	ARM_SFRM* GetRefModel();

	void CreateConvexityModel(bool withSmile);
	double GetSubPrice(int startRowIdx,int endRowIdx,TARNColAlias columnName,
                       const string& evalDateStr,const ARM_DealDescription& dealDesc) const;
    void ComputeEquivalentDatas(const ARM_DealDescription& dealDesc,
                                const pair<ARM_Swap*,int>& stdSwap,
                                ARM_Swaption* swaption, const string& evalDateStr,
                                std::vector<double>& equivDatas);
	
	ARM_StdPortfolioPtrVector CreateTARNCapletFloorletAndDigital();
	ARM_StdPortfolioPtr CreateDiagonalSwaption();
	ARM_StdPortfolioPtr CreateSumOptionPortfolio();
	void ComputeTARNCapletFloorletPrices(
		bool isFreezeWeights, 
		bool isInitParam,
		bool isArrearCpn,
		CalibrationMode capFloorCalibMode,
		CalibrationMode digitalCalibMode, 
		bool computeSwaptionPrice);

	void ComputeDiagonalSwaptionPrice(
		bool isFreezeWeights,
		bool isArrearCpn);

	void PrepareMeanRevForSwaptionCalib();
	void ComputeMCExerciseStrikes();
	void ComputeControlVariable(std::vector<double>& cvprices);
	
	/// calibration part
	void CreateEmptyCalibration(CalibrationMode capFloorCalibMode, CalibrationMode digitalCalibMode,bool swaptionCalibFlag);
	void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	ARM_DateStrip* GetCalibSchedule(std::vector<double>&& pStrike) const;
	ARM_VanillaSecurityDensity* GetMarketDensity(
			double julianRD, double julianSD, double julianED,
			ARM_VolCurve* pAlpha, ARM_VolCurve* pBeta, ARM_VolCurve* pRho, ARM_VolCurve* pNu, 
			int sabrFlag, ARM_ZeroCurve* pCurve, double strike) const;
};

CC_END_NAMESPACE()

#endif
