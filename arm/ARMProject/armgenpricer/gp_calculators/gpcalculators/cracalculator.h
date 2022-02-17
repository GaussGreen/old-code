/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file CRAcalculator.h
 *
 *  \brief
 *
 *	\author  H. BAKHTRI & M. ABDELMOUMNI & P. LAM
 *	\version 1.0
 *	\date June 2005
 */

#ifndef _INGPCALCULATORS_CRACALCULATOR_H
#define _INGPCALCULATORS_CRACALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gencalculator.h"
#include "gpbase/countedptr.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/datestripcombiner.h"

/// kernel
#include "util/refvalue.h"
#include "util/exercise.h"
#include "ccy/currency.h"
//#include "inst/optionportfolio.h"
//#include "mod/bssmiled.h"
//#include "mod/armfrmmodel.h"



/// STL
CC_USING_NS(std,pair)

/// forward declaration
class ARM_Portfolio;
class ARM_Swap;
class ARM_Swaption;
class ARM_VolCurve;
class ARM_ZeroCurve;
class ARM_VolLInterpol;
class ARM_BSModel;
class ARM_BSSmiledModel;
class ARM_CapFloor;

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStrip;
class ARM_SFRM;
class ARM_OptionPortfolio;
class ARM_FRMModel;
class ARM_IRIndex;
class ARM_CRACalculator : public ARM_GenCalculator
{
protected:

	static const string CRAColNamesTable [];
	static const int CRAProductToPriceColumns [];

public:

	enum PayIndexType
    {
        Fixed = K_FIXED,
		Libor = K_LIBOR,
		Cms   = K_CMS,
    };

	enum IndexLeg
	{
		Fund=0,
		Boosted,
		RefIndex,
	};

	enum mdmKeysAlias
    {
        YcKey =  0,
		OswModelKey,
		CfModelKey,
        MrsKey,
		BetaKey,
		YcFundKey,
		YcBasis,
		YcFundBasis,
		ForexKey,
        
        NbKeys
    };
	
	enum CRAColAlias
	{
		ResetDate = 0,
		StartDate,
		NextStartDate,
		MaturityDate,
		Fees,
		Funding,
		Corridor,
		ExoticSwap,
		ExoticSwap1,
		Option, 
		Bermuda,
		FundingLeg,
		CorridorLeg,
		ExoticSwap2,
		Option2, 
		Bermuda2,
		FwdCorridor,
		FwdFunding,
	};

	enum CRAProductToPriceAlias
	{
		BermudaPrice,
		CorridorPrice,
		FundingPrice,
		FwdCorridorPrice,
		FwdFundingPrice,
		Bermuda2Price,
		ExoticSwapPrice,

		NbProductsToPrice,
	};

	ARM_CRACalculator(const ARM_MarketData_ManagerRep& mktDataManager );

	//Constructor 
	ARM_CRACalculator(ARM_Currency& ccy,
					  ARM_Date&	startDate,
					  ARM_Date&	endDate,
					  int payRec,
					  ARM_ReferenceValue&  notional,
					  int callFreq,
					  int callNotice,
					  string callCal,
					  ARM_ReferenceValue&  callFees,
					  int fundFreq,
					  ARM_ReferenceValue& fundSpread,
					  int fundDayCount,
					  ARM_ReferenceValue& cpnSpread,
					  int cpnPayFreq,
					  string cpnResetCal,
					  string cpnPayCal,
					  int boostedIndexType,
					  ARM_ReferenceValue& boostedFixRate,
					  string boostedVarTerm,
					  int boostedResetGap,
					  int boostedResetTiming,
					  int boostedDayCount,
					  int boostedAdjRule,
					  int boostedIntRule,
					  ARM_ReferenceValue& cpnBarrierDown,
					  ARM_ReferenceValue& cpnBarrierUp,
					  int cpnResetFreq,
					  int cpnResetTiming,
					  int cpnResetGap,
					  int refIndexType,
					  string refTerm,
					  int refDayCount,
					  double refCoeff,
					  double meanRevMin,
                      double meanRevMax,
                      double betaMin,
                      double betaMax,
                      ARM_Vector* calibSecPFParams,
                      int nbSteps,
                      int flagToGenerateOSWATM,
					  ARM_StringVector& mdmKeys,
					  const ARM_MarketData_ManagerRep& mktDataManager,
				      const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
					  int reCalibMrs = -1,
					  int reCalibBeta = -1,
					  bool isStdCalib = true,
					  ARM_Currency* fundccy = NULL,
					  const ARM_ReferenceValue& fundNominal= ARM_ReferenceValue());

	// Simple Contsructor just to transfer datas usefull for dealdescription
	ARM_CRACalculator(const ARM_Currency& ccy,
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
					  const ARM_ReferenceValue& cpnSpread,
					  int cpnPayFreq,
					  const string& cpnResetCal,
					  const string& cpnPayCal,
					  int boostedIndexType,
					  const ARM_ReferenceValue& boostedFixRate,
					  const string& boostedVarTerm,
					  int boostedResetGap,
					  int boostedResetTiming,
					  int boostedDayCount,
					  int boostedAdjRule,
					  int boostedIntRule,
					  const ARM_ReferenceValue& cpnBarrierDown,
					  const ARM_ReferenceValue& cpnBarrierUp,
					  int cpnResetFreq,
					  int cpnResetTiming,
					  int cpnResetGap,
					  int refIndexType1,
					  const string& refTerm1,
					  int refDayCount1,
					  double refCoeff1,
					  bool isPortfolioNoticeDays,
					  const ARM_StringVector& mdmKeys,
					  const ARM_MarketData_ManagerRep& mktDataManager,
				      const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
					  int reCalibMrs = -1,
					  int reCalibBeta = -1,
					  bool isStdCalib = true,
					  ARM_Currency* fundccy = NULL,
					  const ARM_ReferenceValue& fundNominal= ARM_ReferenceValue());

	//Summit  CRA Constructor 
	ARM_CRACalculator(ARM_OptionPortfolio* optPortfolio,
					  double meanRevMin,
                      double meanRevMax,
                      double betaMin,
                      double betaMax,
                      ARM_Vector* calibSecPFParams,
                      int nbSteps,
                      int flagToGenerateOSWATM,
					  ARM_StringVector& mdmKeys,
					  const ARM_MarketData_ManagerRep& mktDataManager,
				      const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
					  int reCalibMrs = -1,
					  int reCalibBeta = -1,
					  bool isStdCalib = true);

	ARM_CRACalculator(const ARM_Date& asOfDate,
					  ARM_Security* security);

	ARM_CRACalculator(const ARM_CRACalculator&);
	
	ARM_CRACalculator& operator=(const ARM_CRACalculator&);
	
	void Set(ARM_OptionPortfolio* optPortfolio,
			 double meanRevMin,
			 double meanRevMax,
			 double betaMin,
			 double betaMax,
			 ARM_Vector* calibSecPFParams,
			 int nbSteps,
			 int flagToGenerateOSWATM,
			 ARM_StringVector& mdmKeys,
			 const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
			 int reCalibMrs,
			 int reCalibBeta,
			 int volType,
			 bool isStdCalib);

	void Set( ARM_Currency& ccy,
			  ARM_Date&	startDate,
			  ARM_Date&	endDate,
			  int payRec,
			  ARM_ReferenceValue&  notional,
			  int callFreq,
			  int callNotice,
			  string callCal,
			  ARM_ReferenceValue&  callFees,
			  int fundFreq,
			  ARM_ReferenceValue& fundSpread,
			  int fundDayCount,
			  ARM_ReferenceValue& cpnSpread,
			  int cpnPayFreq,
			  string cpnResetCal,
			  string cpnPayCal,
			  int boostedIndexType,
			  ARM_ReferenceValue& boostedFixRate,
			  string boostedVarTerm,
			  int boostedResetGap,
			  int boostedResetTiming,
			  int boostedDayCount,
			  int boostedAdjRule,
			  int boostedIntRule,
			  ARM_ReferenceValue& cpnBarrierDown,
			  ARM_ReferenceValue& cpnBarrierUp,
			  int cpnResetFreq,
			  int cpnResetTiming,
			  int cpnResetGap,
			  int refIndexType1,
			  string refTerm1,
			  int refDayCount1,
			  double refCoeff1,
			  double meanRevMin,
			  double meanRevMax,
			  double betaMin,
			  double betaMax,
			  int nbSteps,
			  ARM_Vector* calibSecPFParams,
			  int flagToGenerateOSWATM,
			  ARM_StringVector& mdmKeys,
			  const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
			  int reCalibMrs,
			  int reCalibBeta,
			  int volType,
			  bool isStdCalib);

	~ARM_CRACalculator();

	virtual ARM_Object*		Clone() const;
	
	virtual ARM_CstManagerPtr CreateCstManager();

	//ACCESSORS
	//virtual bool					GetIsSpreadIndex() const {return itsIsSpreadIndex;}
	virtual bool					GetIsPortfolioNoticeDays() const {return itsIsPortfolioNoticeDays;}

	ARM_Swaption*					GetSwaption(){return itsSwaption;};
	ARM_OptionPortfolio*			GetOptionPortfolio(){return itsOptionPortfolio;};
	inline	const int				GetNbSteps(){return itsNbSteps;}; 
	inline  const int				GetToCalibrateBeta(){return itsToCalibrateBeta;};
	inline  const int				GetToCalibrateMeanReversion(){return itsToCalibrateMeanReversion;};

	inline  const ARM_Currency&		GetCcy() const {return itsCcy;};
    inline  const ARM_Date&			GetStartDate() const { return itsStartDate;}
	inline  const ARM_Date&			GetEndDate() const { return itsEndDate;}   
	inline  const int				GetPayRec() const {return itsPayRec;}

	inline  int						GetCallFreq() const {return itsCallFreq;}
	inline  int						GetCallNotice() const {return itsCallNotice;}
	inline const char*				GetCallCal() const {return itsCallCal.c_str();}

	ARM_INDEX_TYPE					GetIndexType(int type);	

	inline  int						GetFundFreq() const {return itsFundFreq;}
	inline  int						GetFundDayCount() const {return itsFundDayCount;}

	inline  int						GetBoostedIndexType() const {return itsBoostedIndexType;}
	inline  int						GetCpnPayFreq() const {return itsCpnPayFreq;}
	inline  int						GetBoostedDayCount() const {return itsBoostedDayCount;}
	inline  int						GetBoostedIntRule() const {return itsBoostedIntRule;} 
	inline string					GetBoostedVarIndexTerm() const {return itsBoostedVarTerm;}
	inline int						GetBoostedResetTiming() const {return itsBoostedResetTiming;}
	inline int						GetCpnResetFreq() const {return itsCpnResetFreq;}
	inline int						GetCpnResetTiming() const {return itsCpnResetTiming;}
	inline int						GetCpnResetGap() const {return itsCpnResetGap;}
	inline int						GetRefIndexType1() const {return itsRefIndexType1;}
	inline string					GetRefTerm1() const {return itsRefTerm1;}
	inline  int						GetRefDayCount1() const {return itsRefDayCount1;}
	inline double					GetRefCoeff1() const {return itsRefCoeff1;};
	inline string					GetCpnResetCal() const {return itsCpnResetCal;}
	inline string					GetCpnPayCal() const {return itsCpnPayCal;}

	inline bool						GetSfrmStdCalib() const {return itsSFRMStdCalib;};

	//Get Step-up vectors
	virtual inline ARM_ReferenceValue&		GetCallFees(){return itsCallFees;}
	virtual inline ARM_ReferenceValue&		GetNotional(){return itsNotional;}
	virtual inline ARM_ReferenceValue&		GetFundSpread(){return itsFundSpread;}
	virtual inline ARM_ReferenceValue&		GetCpnSpread(){return itsCpnSpread;}
	virtual inline ARM_ReferenceValue&		GetBoostedFixRate(){return itsBoostedFixRate;}
	virtual inline ARM_ReferenceValue&		GetCpnBarrierDown(){return itsCpnBarrierDown;}
	virtual inline ARM_ReferenceValue&		GetCpnBarrierUp(){return itsCpnBarrierUp;}

	//MUTATORS
	void SetSwaption(ARM_Swaption* swaption);

	void SetOptionPortfolio(ARM_OptionPortfolio* optPortfolio)
	{
		/*if (itsOptionPortfolio)
		{
		   delete itsOptionPortfolio;
		   itsOptionPortfolio = NULL;
		}
		if (optPortfolio)
		   itsOptionPortfolio = (ARM_OptionPortfolio *) optPortfolio->Clone();*/
	};
	void SetCcy(const ARM_Currency& ccy){itsCcy=ccy;};
	void SetStartDate(ARM_Date& startDate){itsStartDate = startDate;};
	void SetEndDate(ARM_Date& endDate){itsEndDate = endDate;};
	void SetPayRec(int payRec){itsPayRec = payRec;};
	void SetNotional(ARM_ReferenceValue& notional){itsNotional = notional;};
	void SetCallFreq(int callFreq){itsCallFreq = callFreq;};
	void SetCallNotice(int callNotice){itsCallNotice = callNotice;};
	void SetCallCal(string& callCal){itsCallCal = callCal;};
	void SetCallFees(ARM_ReferenceValue& callFees){itsCallFees = callFees;};
	void SetFundFreq(int fundFreq){itsFundFreq = fundFreq;};
	void SetFundSpread(ARM_ReferenceValue& fundSpread){itsFundSpread = fundSpread;};
	void SetFundDayCount(int fundDayCount){itsFundDayCount = fundDayCount;};
	void SetCpnSpread(ARM_ReferenceValue& cpnSpread){itsCpnSpread = cpnSpread;};
	void SetCpnPayFreq(int cpnPayFreq){itsCpnPayFreq = cpnPayFreq;};
	void SetCpnResetCal(string& cpnResetCal){itsCpnResetCal = cpnResetCal;};
	void SetCpnPayCal(string& cpnPayCal){itsCpnPayCal = cpnPayCal;};
	void SetBoostedIndexType(int boostedIndexType){itsBoostedIndexType = boostedIndexType;};
	void SetBoostedFixRate(ARM_ReferenceValue& boostedFixRate){itsBoostedFixRate = boostedFixRate;};
	void SetBoostedVarTerm(string& boostedVarTerm){itsBoostedVarTerm = boostedVarTerm;};
	void SetBoostedResetGap(int boostedResetGap){itsBoostedResetGap = boostedResetGap = boostedResetGap;};
	void SetBoostedResetTiming(int boostedResetTiming){itsBoostedResetTiming = boostedResetTiming;};
	void SetBoostedDayCount(int boostedDayCount){itsBoostedDayCount = boostedDayCount;};
	void SetBoostedAdjRule(int boostedAdjRule){itsBoostedAdjRule = boostedAdjRule;};
	void SetBoostedIntRule(int boostedIntRule){itsBoostedIntRule = boostedIntRule;};
	void SetCpnBarrierDown(ARM_ReferenceValue& cpnBarrierDown){itsCpnBarrierDown = cpnBarrierDown;};
	void SetCpnBarrierUp(ARM_ReferenceValue& cpnBarrierUp){itsCpnBarrierUp = cpnBarrierUp;};
	void SetCpnResetFreq(int cpnResetFreq){itsCpnResetFreq = cpnResetFreq;};
	void SetCpnResetTiming(int cpnResetTiming){itsCpnResetTiming = cpnResetTiming;};
	void SetCpnResetGap(int cpnResetGap){itsCpnResetGap = cpnResetGap;};
	void SetRefIndexType1(int refIndexType){itsRefIndexType1 = refIndexType;};
	void SetRefTerm1(string& refTerm){itsRefTerm1 = refTerm;};
	void SetRefDayCount1(int refDayCount){itsRefDayCount1 = refDayCount;};
	void SetRefCoeff1(double refCoeff1){itsRefCoeff1 = refCoeff1;};
	void SetMeanRevMin(double meanRevMin){itsMeanRevMin = meanRevMin;};
	void SetMeanRevMax(double meanRevMax){itsMeanRevMax = meanRevMax;};
	void SetBetaMin(double betaMin){itsBetaMin = betaMin;};
	void SetBetaMax(double betaMax){itsBetaMax = betaMax;};
	void SetCalibSecPFParams(ARM_Vector* calibSecPFParams){itsCalibSecPFParams = calibSecPFParams;};
	void SetNbSteps(int nbSteps){itsNbSteps = nbSteps;};
	void SetFlagToGenerateOSWATM(int flagToGenerateOSWATM){itsFlagToGenerateOSWATM = flagToGenerateOSWATM;};
	void SetProductsToPrice(const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice){itsProductsToPrice = productsToPrice;};
	void SetVolType(int volType){itsVolType = volType;};
	void SetIsPortfolioNoticeDays(bool flag){itsIsPortfolioNoticeDays = flag;};
	void SetSfrmStdCalib(bool stdCalib){itsSFRMStdCalib = stdCalib;};

	void ResetHasBeenPriced() {itsHasBeenPriced=false;}

	//Some pertinent checks
	void							CheckCRAInputs();
	virtual void					CheckData();
	virtual void					CheckMktData();

	//Deal Description
	virtual ARM_StringVector		PricedColumnNames() const;
	virtual ARM_RowInfo				ColumnNames() const;
	virtual ARM_RowInfo				MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	void							InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;
	virtual ARM_DateStripCombiner	DatesStructure() const;
	virtual ARM_DateStripCombiner	CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	//Calibration & Model
	ARM_OptionPortfolio* GenerateOptionPortfolio(void);
	void GenerateProductDescription(ARM_OptionPortfolio* port);
	void GenerateProductDescription(ARM_Swaption* swaption);

	void CreateCalibMethods(void);
	ARM_CalibMethod* GetOSWCalibMethod(void) const;
	const ARM_StdPortfolioPtr GetOSWPortfolio(void) const;

	virtual ARM_PricingModel* GetSFRMModel(void) const;

    // Overloaded  methods
	virtual void UpdateCalibration(bool isUpdateStrike=true);
	virtual void CreateAndSetModel();
	virtual void CreateAndSetCalibration();
    virtual void UpdateModel();
	virtual void Calibrate();

	// Utilitiecs
	ARM_IRIndex GetIndex();
    void Initialize();
	void Pre_InitialiseFRMModel();
	virtual std::vector<double>& ComputeExerciseDatesAndSetFees() const;

	ARM_DateStripPtr GetExerDateStrip() const { return itsExerDateStrip; }
	void SetExerDateStrip(const ARM_DateStripPtr& dateStrip) const { itsExerDateStrip = dateStrip; }

	inline ARM_StdPortfolioPtr GetSigmaPF(){return itsSigmaPf;};
	inline ARM_StdPortfolioPtr GetMrsPF(){return itsBetaPf;};
	inline ARM_StdPortfolioPtr GetBetaPF(){return itsBetaPf;};

	const ARM_ModelParam& GetMRS() const;
    void  SetMRS(ARM_ModelParam* mrsParam);
    const ARM_ModelParam& GetBeta() const;
    void  SetBeta(ARM_ModelParam* mrsParam);
	
	int GetReCalibMrs(void){return itsReCalibMrs;};
	void SetReCalibMrs(int reCalibFlag){itsReCalibMrs = reCalibFlag;};
	int GetReCalibBeta(void){return itsReCalibBeta;};
	void SetReCalibBeta(int reCalibFlag){itsReCalibBeta = reCalibFlag;};
	void SetReCalibFlags(int reCalibMrs, int reCalibBeta);

	void Bump(ARM_ZeroCurve* ZeroCurve,
			  ARM_VolCurve* CapVol,
			  ARM_VolCurve* RhoCap,
			  ARM_VolCurve* NuCap,
			  ARM_VolCurve* SwptVol,
			  ARM_VolCurve* betaCap,
			  ARM_VolCurve* RhoSwopt,
			  ARM_VolCurve* NuSwopt,
			  ARM_VolCurve* BetaSwopt);

	//Outputs
	virtual void					ComputePricingData() const;
	virtual double					Price();	
	
	// Some utilities
	virtual string					ExportShortName() const { return "LCCRA";}
	string							toString(const string& indent = "", const string& nextIndent = "") const;
	virtual string					GeneralDataToString(const string& indent = "", const string& nextIndent = "") const;
	string							DealDesDataToString(const string& indent = "", const string& nextIndent = "") const;
	virtual void					View(char* id = NULL, FILE* ficOut = NULL) const;

	void CleanUp();
    void CopyNoCleanUp(const ARM_CRACalculator& rhs);

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  { return ARM_DateStripPtr(NULL);};
	inline virtual ARM_DateStripPtr GetRefFundDateStrip() const { return ARM_DateStripPtr(NULL);};
	inline virtual ARM_DateStripPtr GetRefDateStrip() const { return ARM_DateStripPtr(NULL);};

	/// Vector
	inline virtual std::vector<double> GetvCpnNominal() const  { return std::vector<double>(0);};
	inline virtual std::vector<double> GetvFundNominal() const { return std::vector<double>(0);};
	inline virtual std::vector<double> GetvFundSpread() const { return std::vector<double>(0);};

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const			{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));			};
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve()	const			{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcFundKey]));     };
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasis]));		};
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve()	const	{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcFundBasis]));   };
	inline virtual ARM_Forex* GetForex() const							{ return dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));			};

	
protected:

	ARM_Swaption*				itsSwaption;
	// PORTFOLIOS
	ARM_OptionPortfolio*		itsOptionPortfolio;

    ARM_StdPortfolioPtr			itsSigmaPf;
	ARM_StdPortfolioPtr			itsMrsPf;
	ARM_StdPortfolioPtr			itsBetaPf;

	// CALIBRATION PARAMETERS
    ARM_Vector*					itsCalibSecPFParams;

    // flags that depends on Min and Max values :
	int							itsToCalibrateMeanReversion;
    int							itsToCalibrateBeta;
	// flags set by user :
	int							itsReCalibMrs; 
	int							itsReCalibBeta; 

    ARM_FRMModel*				itsFRMModel;

    ARM_Vector*					itsSigmaCurve;
    ARM_Vector*					itsBetaOrShift;

    int							itsVolType; // = K_DIAG;

    int							itsFlagToGenerateOSWATM; // Init to NO (0)

    double						itsMeanRevMin;
    double						itsMeanRevMax;
    double						itsBetaMin;
    double						itsBetaMax;

    int							itsNbSteps;
	double						itsMeanRev;

	// GENERAL DATAS
    ARM_Currency				itsCcy;
	ARM_Date					itsStartDate;           
	ARM_Date					itsEndDate;    
	int							itsPayRec;
	ARM_ReferenceValue			itsNotional;
	ARM_ReferenceValue			itsRealFundNotional;

	// CALL DATAS
	int							itsCallFreq;
	int							itsCallNotice; // Positive value
	string						itsCallCal;
	CC_IS_MUTABLE ARM_ReferenceValue			itsCallFees;   // NoticeDates, Associated Fees 

	// FUNDING LEG / Swap Var leg
	int							itsFundFreq;
	ARM_ReferenceValue			itsFundSpread;
	int							itsFundDayCount; //ACT360,ACT365
	int							itsFundAdjRule;  //F, MF...
	int							itsFundIntRule;  //ADJ/UNADJ
	int							itsFundStubRule; //SS, LS...
	
	// COUPON LEG  / Swap fix leg
	ARM_ReferenceValue			itsCpnSpread;
	int							itsCpnPayFreq; 
	string						itsCpnResetCal; //Used in case of variable payment: Fixing calendar
	string                      itsCpnPayCal;	//Payment calendar
	int							itsCpnAdjRule;  //F, MF...
	int							itsCpnIntRule;  //ADJ/UNADJ
	int							itsCpnStubRule; //SS, LS...

	// Boosted index type
	int							itsBoostedIndexType; // Fixed, Libor, Cms.
	ARM_ReferenceValue			itsBoostedFixRate;	 // Used in case of fixed boosted payment. 5%
	string                      itsBoostedVarTerm;   // Used in case of variable boosted payment: 6M, 5Y
	int							itsBoostedResetGap;
	int							itsBoostedResetTiming;   
	int							itsBoostedDayCount;      // A360, ACT365: used to compute the coverage (IT).
	int							itsBoostedAdjRule;       // F, MF, P, PF
	int							itsBoostedIntRule;       // ADJ/UNADJ

	// Corridor index
	ARM_ReferenceValue			itsCpnBarrierDown;
	ARM_ReferenceValue			itsCpnBarrierUp;
	int							itsCpnResetFreq; 
	int							itsCpnResetTiming;      // ADV/ARR
	int							itsCpnResetGap;

	int							itsRefIndexType1;       // Libor, Cms
	string                      itsRefTerm1;            // 6M, 5Y
	int							itsRefDayCount1;        // ACT360, ACT365
	double						itsRefCoeff1;            //Used for CRA and Local CRA.

	bool						itsIsPortfolioNoticeDays;

	// Pricing Data
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/				itsProductsToPrice;

	bool						itsHasBeenPriced;

	double						itsBermuda1Price;
	double						itsBermuda2Price;
	double						itsFundingPrice;
	double						itsCorridorPrice;
	double						itsExoticSwapPrice;
	double						itsFwdCorridorPrice;
	double						itsFwdFundingPrice;

	std::vector<double>				itsCorridorLegPrices;

	// Date Strip
	mutable ARM_DateStripPtr	itsExerDateStrip;

	bool						itsSFRMStdCalib;
};

CC_END_NAMESPACE()

#endif
