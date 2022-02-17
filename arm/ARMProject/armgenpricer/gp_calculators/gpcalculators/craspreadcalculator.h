/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file CRASpreadCalculator.cpp
 *
 *  \brief
 *  Calculator to valuate Callable Range Accrual on spread
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date August 2005
 */

#ifndef _INGPCALCULATORS_CRASPREADCALCULATOR_H
#define _INGPCALCULATORS_CRASPREADCALCULATOR_H

#include "cralocalcalculator.h"

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gencalculator.h"
#include "gpbase/countedptr.h"
#include "gpbase/gpmatrix.h"
#include "gpmodels/marketirmodel.h"
#include "gpmodels/typedef.h"
#include "gpinfra/curvemodelparam.h"

//#include "volcube.h"

/// STL
CC_USING_NS(std,pair)

class ARM_SpreadOption;
class ARM_Swaption;
class ARM_ZeroCurve;
class ARM_CorridorDblCondition;

CC_BEGIN_NAMESPACE( ARM )
class ARM_CorrelManager;
struct ARM_VanillaSpreadOptionArg;
struct ARM_VanillaSwaptionArg;

class ARM_CRASpreadCalculator : public ARM_CRALocalCalculator
{
public:

	enum mdmKeysAlias
    {
        YcKey =  0,
		OswModelKey,
		CfModelKey,
		MrsKey,
		VolRatioKey,
		MrsSpreadKey,
		CorrelKey,
		CorrelRateSpreadKey,
		YcFundKey,
		YcBasis,
		YcFundBasis,
		ForexKey,
        
        NbKeys
    };

	enum CalibrationType
	{
		DIAG_CALIBRATION = 0,
		BASKET_CALIBRATION,
		BASKET_CALIBRATION_SIMPLIFIED,
		DIAG_SPREAD_LONG,
		DIAG_SPREAD_SHORT,
		DIAG_SPREAD_INDEX,
		SHORT_LONG_SPREAD
	};

	enum CalibStrikeType
	{
		ATM = 0,
		EQUIVALENT,
		ZERO,
		CAP,
		FLOOR,
		FRONTIER
	};

	//Constructor 
	ARM_CRASpreadCalculator( 
		const ARM_Currency& ccy,
		const ARM_Date&	startDate,
		const ARM_Date&	endDate,
		int payRec,
		const ARM_ReferenceValue&  notional,
		int callFreq,
		int callNotice,
		const string& callCal,
		const ARM_ReferenceValue&  callFees,
		int fundFreq,
		const ARM_ReferenceValue& fundSpread,
		int fundDayCount,
		int cpnDayCount,
		int cpnPayFreq,
		const string& cpnResetCal,
		const string& cpnPayCal,
		int payIndex,
		int payIndexResetTiming,
		const ARM_ReferenceValue& boostedFixRate,
		const ARM_ReferenceValue& payIndexMult,
		const ARM_ReferenceValue&  ,
		const ARM_ReferenceValue& cpnBarrierUp,
		int cpnResetFreq,
		int cpnResetTiming,
		int cpnPayGap,
		int refIndex1,
		double refCoeff1,
		int refIndex2,
		double refCoeff2,
		int refIndex3 = K_FIXED,
		const ARM_ReferenceValue& barrierDown3 = ARM_ReferenceValue(),
		const ARM_ReferenceValue& barrierUp3 = ARM_ReferenceValue(),
		bool triple = false,
		bool vms = false,
		const ARM_ReferenceValue& boostedFixRate2 = ARM_ReferenceValue(),
		const ARM_ReferenceValue& cpnBarrierDown2 = ARM_ReferenceValue(),
		const ARM_ReferenceValue& cpnBarrierUp2 = ARM_ReferenceValue(),
		const ARM_ReferenceValue& boostedFixRate3 = ARM_ReferenceValue(),
		const ARM_ReferenceValue& cpnBarrierDown3 = ARM_ReferenceValue(),
		const ARM_ReferenceValue& cpnBarrierUp3 = ARM_ReferenceValue(),
		const ARM_ReferenceValue& refTenor1 = ARM_ReferenceValue(),
		ARM_ModelType		modelType = ARM_PricingModelType::HWM1F,
		CalibrationType	calibrationType = DIAG_CALIBRATION,
		vector<CalibStrikeType>	calibStrikeType = vector<CalibStrikeType>(),
		ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod = ARM_MarketIRModel::MONEYNESS,
		const ARM_StringVector& mdmKeys = ARM_StringVector(),
		const ARM_MarketData_ManagerRep& mktDataManager = ARM_MarketData_ManagerRep(ARM_Date()),
		const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice = std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/(),
		const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& localCalibFlags = std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/(5,false),
		const std::vector<double>& optimResetData = std::vector<double>(0),
		bool isExerciseProbas=false,
		size_t exerciseProbaOffset=0,
		ARM_Currency* fundccy = NULL,
		const ARM_ReferenceValue& fundNominal= ARM_ReferenceValue());

	ARM_CRASpreadCalculator(ARM_OptionPortfolio* optionPortfolio,
							ARM_StringVector& mdmKeys,
							ARM_ModelType modelType,
							CalibrationType	calibrationType,
							vector<CalibStrikeType>	calibStrikeType,
							ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod,
							ARM_ReferenceValue& payIndexMult,
							const ARM_MarketData_ManagerRep& mktDataManager,
							const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
							const std::vector<double>& optimResetData = std::vector<double>(0));

	ARM_CRASpreadCalculator(const ARM_Date& asOfDate,
							ARM_OptionPortfolio* optionPortfolio,
							ARM_ReferenceValue& payIndexMult,
							const std::vector<double>& optimResetData = std::vector<double>(0));

	ARM_CRASpreadCalculator(const ARM_Date& asOfDate,
							ARM_Swaption* swaption,
							ARM_ReferenceValue& payIndexMult,
							const std::vector<double>& optimResetData = std::vector<double>(0));

// FIXMEFRED: mig.vc8 (25/05/2007 14:53:08):missing return type
	void InitCRASpreadFromSummit(ARM_ZeroCurve* zc,
							ARM_VolCurve* capVol,
							ARM_VolCurve* rhoCap,
							ARM_VolCurve* nuCap,
							ARM_VolCurve* betaCap,
							ARM_VolCurve* swoptVol,
							ARM_VolCurve* rhoSwopt,
							ARM_VolCurve* nuSwopt,
							ARM_VolCurve* betaSwopt,
							int SABRSigmaOrAlpha,
							ARM_VolCurve* convAdjustVolCap,
							ARM_VolCurve* convAdjustVolSwopt,
							long convAdjustType,
							ARM_CorrelManager* correlCorr,
							ARM_VolCurve* correlDiagCap,
							ARM_VolCurve* correlDiagSwopt,
							double mrs,
							double volRatio,
							double mrsSpread,
							double correl,
							ARM_ModelType		modelType,
							CalibrationType	calibrationType,
							vector<CalibStrikeType>	calibStrikeType,
							ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod,
							const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
							const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& localCalibFlags);

	ARM_CRASpreadCalculator(const ARM_CRASpreadCalculator&);
	
	ARM_CRASpreadCalculator& operator = (const ARM_CRASpreadCalculator&);
	
	~ARM_CRASpreadCalculator();

	virtual ARM_Object*		Clone() const;

	// Create the date structure
	ARM_DateStripCombiner DatesStructure() const;

	// Create the constant manager
	virtual ARM_CstManagerPtr CreateCstManager();

	// Return the index type for the swaption
	ARM_INDEX_TYPE GetIndexType();

    // Overloaded  methods
	virtual void UpdateCalibration(bool isUpdateStrike=true);
	virtual void CreateAndSetCalibration();
	virtual void CreateAndSetModel();
    virtual void UpdateModel();
	virtual void Calibrate();
	virtual ARM_StringVector PricedColumnNames() const;
	virtual ARM_RowInfo	MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;

	void SetPayIndexMult(ARM_ReferenceValue& payIndexMult){itsPayIndexMult = payIndexMult;};
	virtual inline ARM_ReferenceValue&		GetPayIndexMult(){return itsPayIndexMult;}

	//Accessors / Mutators for reference index 2 .
	void SetRefIndexType2(int refIndexType){itsRefIndexType2 = refIndexType;};
	void SetRefTerm2(string& refTerm){itsRefTerm2 = refTerm;};
	void SetRefDayCount2(int refDayCount){itsRefDayCount2 = refDayCount;};
	void SetRefCoeff2(double refCoeff1){itsRefCoeff2 = refCoeff1;};

	inline int						GetRefIndexType2() const {return itsRefIndexType2;}
	inline string					GetRefTerm2() const {return itsRefTerm2;}
	inline  int						GetRefDayCount2() const {return itsRefDayCount2;}
	inline double					GetRefCoeff2() const {return itsRefCoeff2;};


	ARM_DateStripPtr GetExerDateStrip() const { return itsExerDateStrip; }
	ARM_DateStripPtr GetFundDateStrip() const { return itsFundDateStrip; }
	ARM_DateStripPtr GetStructDateStrip() const { return itsStructDateStrip; }

	ARM_DateStripPtr GetExerDateStripUnadj() const { return itsExerDateStripUnadj; }
	ARM_DateStripPtr GetDateStrip1() const { return itsDateStrip1; }
	ARM_DateStripPtr GetDateStrip2() const { return itsDateStrip2; }
	ARM_DateStripPtr GetDateStrip3() const { return itsDateStrip3; }
	ARM_GP_VectorPtr GetIndexVector() const { return itsIndexVector; }

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const			{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));			};
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve()	const			{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcFundKey]));     };
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasis]));		};
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve()	const	{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcFundBasis]));   };
	inline virtual ARM_Forex* GetForex() const							{ return NULL;}//dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));			};

	// Computes attributes vectors from object ARM_Curves;
	void ComputeProductVectorsFromCurves() ;

	// Create the Corridor Spread Option portfolio
	void CreateEmptyCalibration();
	void CreateEmptyCalibrationFor2F();
	void CreateCSOPortfolio();
	void CreateCSOPortfolio_FakeForDateStripOnly();
	void CreateCSOPortfolio_ForVMS();
	void CreateCSOPortfolio_FakeForDateStripOnly_ForVMS();

	void SetDiagStrikeToFrontier(ARM_StdPortfolioPtr& swaptionPF);

	// Compute the corridor spread option prices
	void ComputeCSOPrices();

	// Accessor to the corridor spread option portfolio
	ARM_StdPortfolioPtr GetCSOPF()	const { return itsCSOPF; };
	ARM_StdPortfolioPtr GetCSOPF2()	const { return itsCSOPF2; };
	ARM_StdPortfolioPtr GetCSOPF3()	const { return itsCSOPF3; };

	// Create the swaption porfolio
	ARM_StdPortfolioPtr CreateSwaptionPortfolio();

	// Compute the swaption portfolio prices
	void ComputeSwaptionPrices();
	void ComputeCalibPortfolioPrices();
	
	// Accessor to the swaption portfolio
	ARM_StdPortfolioPtr GetSwoptPF() const;


	// Create the CMS spread option / CMS caplet porfolio (underlying double corridor conditions)
	ARM_StdPortfolioPtr CreateSOPortfolio(bool isSOCond=true,bool isDouble=true,bool isLong=true);
	ARM_StdPortfolioPtr GetSOPF(bool isSOCond=true) const;
	ARM_StdPortfolioPtr GetCMSLONGPortfolio()		const ;
	ARM_StdPortfolioPtr GetCMSSHORTPortfolio()		const ;


	void ComputeSOPrices(ARM_StdPortfolioPtr& pf,CalibStrikeType cst=ATM,bool isSOBarrier=true);
	void ComputeCRA2NormalSensitivities();
	void ComputeCRA2Probabilities(ARM_CorridorDblCondition* dblCorridor,
								  std::vector<double>& probaIn,std::vector<double>& probaInFLT);
	double ComputeResidualUnderlying(int exerIdx, double& floatLegPv, double& exoLegPv);

	ARM_Date	GetNextNoticeDate() const;
	double		GetFirstProba() const;
	void		ComputeUnderlying(double& mktPrice,double& modelPrice) const;
	
	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  { return itsFundDateStrip;};
	inline virtual ARM_DateStripPtr GetRefFundDateStrip() const { return itsFundDateStrip;};
	inline virtual ARM_DateStripPtr GetRefDateStrip() const { return itsFundDateStrip;};

	/// Vector
	inline virtual std::vector<double> GetvCpnNominal() const  { return itsvFundNominal;};
	inline virtual std::vector<double> GetvFundNominal() const { return itsvRealFundNominal;};
	inline virtual std::vector<double> GetvFundSpread() const { return itsvRealFundSpread;};

	// checks if the payments are floating
	bool IsVariableCorridor() const;
	bool IsVariableCpn(size_t cpnIndex) const;
	bool IsDbleCorridor() const { return itsRefIndexType3 != K_FIXED; }
	bool IsTripleRange() const  { return itsIsTripleRange; }
	bool IsVms() const  { return itsIsVms; }

	bool SetSwoptCalib(bool swoptCalib) {itsSwoptCalib=swoptCalib; return !swoptCalib || swoptCalib && GetSwoptPF() != ARM_StdPortfolioPtr(NULL); }
	void SetCorrelUnSqueezer(bool correlUnSqueezer) { itsCorrelUnSqueezer=correlUnSqueezer; }
	void SetVolUnSqueezer(bool volUnSqueezer) { itsVolUnSqueezer=volUnSqueezer; }	
	void SetPVAdjuster(bool PVAdjuster) { itsPVAdjuster=PVAdjuster; }	

    // In the case of CRA or Swaption from CCSO
	double InterpolMeanRevParam(ARM_OptionPortfolio* optionPortfolio, 
                                ARM_VolLInterpol* meanRevParam, 
                                ARM_Date asOfDate);

    double InterpolMeanRevParam(ARM_Swaption* bermuda, 
                                ARM_VolLInterpol* meanRevParam, 
                                ARM_Date asOfDate);

	/// --------------------------------------------------
	/// ---- variable notional swaption calibration 
	/// --------------------------------------------------
	void ComputeCSONormalSensitivities();
	ARM_Swaption* CreateVarNotionalSwaptionAtExer(int exerIdx);
	ARM_VanillaSpreadOptionArg* CreateVanillaArgCSO(ARM_SpreadOption* spreadOption);


	/// Static function to convert a general swaption to a variable notional one
	static ARM_Swaption* ConvertToVarNotionalSwaption(ARM_ZeroCurve* curve, const ARM_VanillaSwaptionArg& oswVanillaArg);


	// Some utilities
	virtual string					ExportShortName() const { return "LCCSO"; }
	string							toString(const string& indent, const string& nextIndent) const;
	virtual void					View(char* id = NULL, FILE* ficOut = NULL) const;

private:
	/// The Corridor Spread Option portfolio
	ARM_StdPortfolioPtr 	itsCSOPF;
	ARM_StdPortfolioPtr 	itsCSOPF2;
	ARM_StdPortfolioPtr 	itsCSOPF3;
	ARM_ModelType 				itsModelType;
	bool 					itsSwoptCalib;
	CalibrationType			itsCalibType;
	vector<CalibStrikeType>	itsCalibStrikeType;
	ARM_MarketIRModel::VnsPricingMethod  itsVnsPricingMethod; /// ATM or MONEYNESS

	double				itsCoeff1;
	double				itsCoeff2;
	int					itsRefIndex1;
	int					itsRefIndex2;
	int					itsPayIndex;
	ARM_ReferenceValue	itsPayIndexMult;

	int					itsRefIndexType2;        // Libor, Cms
	string		        itsRefTerm2;             // 6M, 5Y
	int					itsRefDayCount2;         // ACT360, ACT365
	double				itsRefCoeff2;            // Used for CRA and Local CRA.

	/// Double corridor stuff
	int					itsRefIndex3;
	int					itsRefIndexType3;        // Libor, Cms
	string		        itsRefTerm3;             // 6M, 5Y
	ARM_ReferenceValue	itsRateBarrierDown;
	ARM_ReferenceValue	itsRateBarrierUp;
	bool				itsCorrelUnSqueezer;
	bool				itsPVAdjuster;

	bool				itsVolUnSqueezer;
	bool				itsBootstrapOptimizer;
	bool				itsSwitchOldCalibration;
			
	/// Triple Range stuff
	bool				itsIsTripleRange;
	ARM_ReferenceValue	itsBoostedFixRate2;
	ARM_ReferenceValue	itsBoostedFixRate3;
	ARM_ReferenceValue	itsCpnBarrierDown2;
	ARM_ReferenceValue	itsCpnBarrierDown3;
	ARM_ReferenceValue	itsCpnBarrierUp2;
	ARM_ReferenceValue	itsCpnBarrierUp3;

	/// VMS stuff
	ARM_ReferenceValue	itsTenor;
	bool				itsIsVms;
	
	// Exercise
	std::vector<double>		itsvExerFees;		// may not be used...
	ARM_BoolVector		itsvIsExerDate;
	size_t				itsExerSize;
	
	// Funding
	double				itsFundLeverage;
	ARM_IntVector		itsvFundIndex;
	std::vector<double>		itsvFundSpread;
	std::vector<double>		itsvFundNominal;
	std::vector<double>		itsvRealFundSpread;
	std::vector<double>		itsvRealFundNominal;
	size_t				itsFundSize;

	// Coupon
	ARM_IntVector		itsvCpnIndex;
	std::vector<double>		itsvCpnNominal;
	size_t itsCpnSize;

	/// --------------------------------------------------
	/// ---- variable notional swaption calibration 
	/// --------------------------------------------------
	
	// vanilla arg used to store CMS date strips
	ARM_VanillaArgPtr itsVanillaArgCSO;

	// sensi of cpns w.r.t long & short cms rate
	std::vector<double> itsFixedCpnLongSensi;
	std::vector<double> itsFixedCpnShortSensi;
	std::vector<double> itsVarCpnLongSensi;
	std::vector<double> itsVarCpnShortSensi;
	std::vector<double> itsPayIndexFwd;/// size = itsCpnSize...
	
	// cpn value (no discount)
	std::vector<double> itsFixedCpnValue;
	std::vector<double> itsVarCpnValue;

	std::vector<double> itsFixedCpnLongSensi2;
	std::vector<double> itsFixedCpnShortSensi2;
	std::vector<double> itsVarCpnLongSensi2;
	std::vector<double> itsVarCpnShortSensi2;
	std::vector<double> itsPayIndexFwd2;
	std::vector<double> itsFixedCpnValue2;
	std::vector<double> itsVarCpnValue2;
	std::vector<double> itsFixedValues2;

	std::vector<double> itsFixedCpnLongSensi3;
	std::vector<double> itsFixedCpnShortSensi3;
	std::vector<double> itsVarCpnLongSensi3;
	std::vector<double> itsVarCpnShortSensi3;
	std::vector<double> itsPayIndexFwd3;
	std::vector<double> itsFixedCpnValue3;
	std::vector<double> itsVarCpnValue3;
	std::vector<double> itsFixedValues3;

	mutable ARM_DateStripPtr itsFundDateStrip;
	mutable ARM_DateStripPtr itsStructDateStrip;
	mutable ARM_DateStripPtr itsExerDateStripUnadj;

	mutable ARM_DateStripPtr itsDateStrip1;
	mutable ARM_DateStripPtr itsDateStrip2;
	mutable ARM_DateStripPtr itsDateStrip3;

	ARM_GP_VectorPtr itsIndexVector;
	void Init(void);

	std::vector<double> itsOptimResetData;
	void SkipResetDates(const std::vector<double>& resets,const std::vector<double>& payments,
						double asOfDate, double skipSize, double limit,
						size_t &resetIdx, ARM_IntVector& resetsToErase) const;

	void OptimiseResetDates(ARM_SpreadOption* sec=NULL,ARM_DateStripPtr& sched=ARM_DateStripPtr(NULL)) const;

	void UpdateLegSchedules(const ARM_IntVector& newToOldReset,
							const ARM_IntVector& newToOldFlowEndDate,
							ARM_SwapLeg* leg=NULL,ARM_DateStripPtr& sched=ARM_DateStripPtr(NULL)) const;

	void AgregateSameResetDates(ARM_SpreadOption* sec);

};

CC_END_NAMESPACE()

#endif
