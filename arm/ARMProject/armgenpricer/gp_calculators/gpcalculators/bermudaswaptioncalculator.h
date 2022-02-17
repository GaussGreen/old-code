/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file bermudaswaptioncalculator.cpp
 *
 *  \brief
 *
 *	\author  H. BAKHTRI & A. TRIKI
 *	\version 1.0
 *	\date June 2005
 */

#ifndef _INGPCALCULATORS_BERMUDASWAPTIONCALCULATOR_H
#define _INGPCALCULATORS_BERMUDASWAPTIONCALCULATOR_H

#include "gencalculator.h"

//gpbase
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/env.h"
#include "gpbase/countedptr.h"
#include "gpbase/gpmatrix.h"

//gpmodels
#include "gpmodels/MarketIRModel.h"
#include "gpmodels/typedef.h"

//kernel
#include <util/refvalue.h>
#include <util/exercise.h>
#include <ccy/currency.h>

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
class ARM_CapFloor;

CC_BEGIN_NAMESPACE( ARM )

class ARM_DateStrip;
class ARM_SFRM;

class ARM_BermudaSwaptionCalculator : public ARM_GenCalculator
{
protected:

	static const string ControlVariateColNamesOutputs [];

	//GenSecType = 1
	static const string BermudaSwaptionColNamesTable1 [];
	static const string ControlVariateColNamesTable1 [];
	
	//GenSecType = 2
	static const string BermudaSwaptionColNamesTable2 [];
	static const string BermudaSwaptionProbaColNamesTable2 [] ;
	static const string BermudaSwaptionProbaPricedColNamesTable2 [];


	//GenSecType = 3
	static const string BermudaSwaptionColNamesTable3 [];
	static const string ControlVariateColNamesTable3 [];

	//GenSecType = 4
	static const string BermudaSwaptionColNamesTable4 [];


public:

	//GenSecType = 1
    enum BermudaSwaptionColAlias1
    {
        ResetDate1=0,
		StartDate1,
		PayDate1,
		SwapStartDate1,
		NextSwapStartDate1,
		LastSwapStartDate1,
		Strike1,
		Fees1,
		Notional1,
		StdFixSwaplet1,
		StdVarSwaplet1,
		StdSwaplet1,
		StdSwap1,
		Option1,
		Bermuda1,
    };

	//GenSecType = 2
	enum BermudaSwaptionColAlias2
    {
        ResetDate2 = 0,
		StartDate2,
		EndDate2,
		SwapStartDate2,
		NextSwapStartDate2,
		LastSwapStartDate2,
		Strike2,
		Fees2,
		Notional2,
		StdFixSwaplet2,
		StdVarSwaplet2,
		StdSwaplet2,
		StdSwap2,
		Option2,
		Bermuda2,
		SwapRate2,
		Frontier2,
		DealEndDate2,
		ExerciseCondition2,
    };

	//GenSecType = 3
	enum BermudaSwaptionColAlias3
    {
       ResetDate3 = 0,
	   StartDate3,
	   EndDate3,
	   SwapStartDate3,
	   Strike3,
	   Fees3,
	   StdSwaplet3,
	   ReverseStdSwaplet3,
	   OptionAnnulation3,
	   Bermuda3,
    };

	//GenSecType = 4
	enum BermudaSwaptionColAlias4
    {
		ResetDate4 = 0,
		StartDate4,
		PastStartDate4,
		SwapStartDate4,
		EndDate4,
		Strike4,
		Fees4,
		DfStart4,
		CvgStart4,
		DfEnd4,
		CvgEnd4,
		Annuity4,
		Flows4,
		CashFlow4,
		Bermuda4,
		SwapRate4,
		Frontier4,
		DealEndDate4,
    };

	enum mdmKeysAlias
    {
        YcKey =  0,
		OswModelKey,
		CfModelKey,
       	NormalModel,
		BetaKey,
		fundingKey,
		YcBasisDomKey,
        YcBasisFundKey,
		ForexKey,
        
        NbKeys
    };

	enum numMethodType
	{
		MonteCarlo = 0,
		Tree,
		PDE,
	};

	enum portfolioMode
	{
		Summit = 0,
		Manual,
		Faster,
		NewMode,
	};

	enum calibStrikes
	{
		Calib = 0,
		ATM,
		Moneyness,
	};

	enum calibMrsMode
	{
		Column = 0,
		Correl,
		Surdiag,
		AntiDiag,
	};

	//Constructor 
	ARM_BermudaSwaptionCalculator(ARM_Currency& ccy,
								  ARM_Date&	startDate,
								  ARM_Date&	endDate,
								  ARM_ReferenceValue&  nominal,
								  ARM_ReferenceValue&  strike,
								  int payRec,
								  int stubRule,
								  int callFreq,
								  int callNotice,
								  string callCal,
								  ARM_Date firstCallDate,
								  ARM_Date lastCallDate,
								  ARM_ReferenceValue& fees,
								  int fixFreq,
								  int fixBasis,
								  string fixPayCal,
								  int fixAdjRule,
								  int fixRule,
								  int fixPayGap,
								  bool isZc,
								  int varFreq,
								  int varBasis,
								  string varResetCal,
								  string varPayCal,
								  string varIndexTerm,
								  int varAdjRule,
								  int varRule,
								  int varResetGap,
								  int varPayGap,
								  ARM_ReferenceValue& varSpread,
								  int genSecType,
								  vector<int>* controlVariate,
								  vector<double>* controlPrices,
								  ARM_StringVector& mdmKeys,
								  const ARM_MarketData_ManagerRep& mktDataManager,
								  int modelType,
								  vector<ARM_ReferenceValue*>* modelParams,
								  bool mrsCalibFlag,
								  bool betaCalibFlag,
								  int numMethodType,
								  int amcIter,
								  int mcIter,
								  int maxBucketSize,
								  string genType1,
								  string genType2,
								  string pathOrder,
								  string pathScheme,
								  int firstNbDims,
								  int treeSteps,
								  vector<int>& portfolioMode,
								  bool fixBoundaryFlag = false,
								  bool approxMarginFlag = false,
								  bool freezeBetasFlag = false,
								  bool calculateProbaFlag = false);
				
	//Summit Constructor 
	ARM_BermudaSwaptionCalculator(const ARM_Date&	asOfDate,
								  ARM_Currency& ccy,
								  ARM_Date&	startDate,
								  ARM_Date&	endDate,
								  ARM_ReferenceValue&  nominal,
								  ARM_ReferenceValue&  strike,
								  int payRec,
								  int stubRule,
								  int callFreq,
								  int callNotice,
								  string callCal,
								  ARM_Date firstCallDate,
								  ARM_Date lastCallDate,
								  ARM_ReferenceValue& fees,
								  int fixFreq,
								  int fixBasis,
								  string fixPayCal,
								  int fixAdjRule,
								  int fixRule,
								  int fixPayGap,
								  bool isZc,
								  int varFreq,
								  int varBasis,
								  string varResetCal,
								  string varPayCal,
								  string varIndexTerm,
								  int varAdjRule,
								  int varRule,
								  int varResetGap,
								  int varPayGap,
								  ARM_ReferenceValue& varSpread,
								  int genSecType,
								  bool isPastZc = false,
								  ARM_Date& pastStartDate = ARM_Date());

	ARM_BermudaSwaptionCalculator(const ARM_Date& asOfDate,
								  ARM_Swaption* swaption,
								  int fixFreqIfZC = K_DEF_FREQ);

	void InitBermudaSwaptionForSummit(  ARM_StringVector& mdmKeys,
									    ARM_MarketData_ManagerRep* mktDataManager,
									    vector<int>* controlVariates,
									    vector<double>* controlPrices,
									    vector<ARM_ReferenceValue*>* modelParams,
										bool mrsCalibFlag,
										bool atmCalibFlag,
										int numMethodType,
										int amcIter,
										int	mcIter,
										int	maxBucketSize,
										string genType,
										int	treeSteps,
										vector<int> portfolioMode,
										bool boundaryFlag,
										bool approxMarginFlag,
										bool freezeBetasFlag,
										int	modelType,
										bool calculateProbaFlag);

	void InitBermudaSwaptionForSummit(  vector<int>* controlVariates,
										vector<double>* controlPrices,
									    vector<ARM_ReferenceValue*>* modelParams,
										bool mrsCalibFlag,
										bool atmDiagonalFlag,
										int numMethodType,
										int amcIter,
										int mcIter,
										int maxBucketSize,
										string genType,
										int treeSteps,
										vector<int> portfolioMode,
										bool boundaryFlag,
										bool approxMarginFlag,
										bool freezeBetasFlag,
										int modelType,
										bool calculateProbaFlag,
										ARM_ZeroCurve*		zcCpn, 
										ARM_VolCurve*		swoptVC, 
										ARM_VolCurve*		capVC, 
										ARM_VolLInterpol*	capRo = NULL, 
										ARM_VolLInterpol*	capNu = NULL,
										ARM_VolLInterpol*	capBeta	= NULL,
										ARM_VolLInterpol*	swoptRo	= NULL, 
										ARM_VolLInterpol*	swoptNu	= NULL,
										ARM_VolLInterpol*	swoptBeta = NULL,
										ARM_MarketIRModel*	normalModel = NULL,
										int 				hedgeUpdate	= 0,
										int					SABRSigmaOrAlpha = 1);

	ARM_BermudaSwaptionCalculator(const ARM_BermudaSwaptionCalculator&);	
	ASSIGN_OPERATOR(ARM_BermudaSwaptionCalculator)	
	~ARM_BermudaSwaptionCalculator();

	ARM_ExerciseStyle* CreateExerciseStyle() const;

	ARM_ReferenceValue*		ConvertVariableMargin(ARM_ReferenceValue& marginCurve, int initFreq, int newFreq, int varBasis, int fixBasis );
	ARM_ReferenceValue*		ConvertVariableMarginToStrike(ARM_ReferenceValue& marginCurve );

	ARM_ReferenceValue*		ComputeMeanNotionalCurve(ARM_ReferenceValue& notionalCurve, int newFreq  );

	virtual ARM_Object*		Clone() const;
	
	//ACCESSORS / MUTATORS
    //Generic Security Type
	const int				GetGenSecType() const;
    void					SetGenSecType(int genSecType);
	void					UpdateGenSecType(int genSecType);

    //General Datas
	inline  const ARM_Date& GetStartDate() const { return itsStartDate;}
	inline  const ARM_Date& GetEndDate() const { return itsEndDate;}   
	const ARM_Currency&		GetCcy() const {return itsCcy;};
	void					SetCcy(const ARM_Currency& ccy){itsCcy=ccy;};
	inline  const int		GetPayRec() const {return itsPayRec;}
	
	//Call Datas
	inline  const ARM_Date& GetFirstCallDate() const { return itsFirstCallDate;}
	inline  const ARM_Date& GetLastCallDate() const { return itsLastCallDate;} 
    inline  const bool		GetFreezeBetasFlag() const { return itsFreezeBetasFlag;}  	
	inline ARM_VectorPtr	GetBetasCoeff () {return itsBetasCoeff;}

	//Underlying Fix Leg Datas
	inline  const int		GetFixFreq() const { return itsFixFreq;}
	inline  void			SetFixFreq(int freq) { itsFixFreq = freq;}
	inline  const int		GetFixBasis() const { return itsFixBasis;}
	inline  const string&	GetFixResetCal() const {return itsVarResetCal;}
	inline  const string&	GetFixPayCal() const {return itsFixPayCal;}
	inline  const int		GetFixAdjRule() const { return itsFixAdjRule;}
	inline  const int		GetFixRule() const { return itsFixRule;}  // ADJ , UNADJ	
	inline  void			SetZc(bool value){itsZC=value;}
	inline  const bool		GetZCFlag() const { return itsZC;}      

	//Underlying Var Leg Datas
	inline  const int		GetCallFreq() const { return itsCallFreq;}
	inline  void			SetCallFreq (int callFreq) {itsCallFreq = callFreq;}
	inline  const int		GetVarFreq() const { return itsVarFreq;}
	inline  void			SetVarFreq (int newFreq) {itsVarFreq = newFreq;}
	inline  const int		GetVarBasis() const { return itsVarBasis;}
	inline  const string&	GetVarResetCal() const {return itsVarResetCal;}
	inline  const string&	GetVarPayCal() const {return itsVarPayCal;}
	const   string&			GetVarIndexTerm() const {return itsVarIndexTerm;};
	void					SetVarIndexTerm(const string& idxTerm){itsVarIndexTerm=idxTerm;};
	ARM_INDEX_TYPE			GetIndexType();	
	inline  const int		GetVarAdjRule() const { return itsVarAdjRule;}
	inline  const int		GetVarRule() const { return itsVarRule;}
	inline  const int		GetVarResetGap() const { return itsVarResetGap;}
	inline  const int		GetStubRule() const { return itsStubRule;}
	
	//Some pertinent checks
	void					CheckSwaptionInputs();
	virtual void			CheckData();
	virtual void			CheckMktData();

	//Calib Parameters
	vector<ARM_ReferenceValue*>* GetModelParams() const;
	void					SetModelParams(vector<ARM_ReferenceValue*>* modelParams);
	double					GetDefaultMrs() const;
	inline void				SetInitMrs(double initMrs);
	inline double			GetInitVol(double date = 0.0){return (*itsModelParams)[0]->Interpolate(date);}
	inline double			GetInitMrs(double date = 0.0){return (*itsModelParams)[1]->Interpolate(date);}
	inline double			GetInitBeta(double date = 0.0){return (*itsModelParams)[2]->Interpolate(date);}
	inline double			GetInitTheta(double date = 0.0){return (*itsModelParams)[3]->Interpolate(date);}
	inline double			GetInitSkew(double date = 0.0){return (*itsModelParams)[4]->Interpolate(date);}
	inline double			GetInitMrsL(double date = 0.0){return (*itsModelParams)[5]->Interpolate(date);}
	inline double			GetInitMrsU(double date = 0.0){return (*itsModelParams)[6]->Interpolate(date);}
	inline bool				GetMrsCalibFlag() const{return itsMrsCalibFlag;}
	inline void				SetMrsCalibFlag(bool flag){itsMrsCalibFlag = flag;}
	inline void				SetATMDiagonalFlag(bool flag){itsATMDiagonalFlag = flag;}
	inline bool				GetATMDiagonalFlag() const{ return itsATMDiagonalFlag;}
	inline int				GetNumMethodType() const {return itsNumMethodType;}
	inline void				SetNumMethodType(int numMethodType) {itsNumMethodType = numMethodType;}
	inline int				GetAmcIter() const {return itsAmcIter;}
	inline void				SetAmcIter(int iter) {itsAmcIter = iter;}
	inline void				SetMcIter(int iter) {itsMcIter = iter;}
	inline int				GetMcIter() const {return itsMcIter;}
	inline int				GetMaxBucketSize() const {return itsMaxBucketSize;}
	inline void				SetMaxBucketSize(int max) {itsMaxBucketSize = max;}
	inline int				GetTreeSteps() const {return itsTreeSteps;}
	inline void				SetTreeSteps(int steps) {itsTreeSteps = steps;}
	inline const int		GetOSWPortfolioMode() const { return itsPortfolioMode[0];}
	inline void				SetPortfolioMode(vector<int> portfolioMode){ itsPortfolioMode = portfolioMode;}
	
	inline const int		GetMrsPortfolioMode() const {return itsPortfolioMode[0];}
	inline const int		GetMrsCalibFreq() const {return itsPortfolioMode[1];}
	inline const int		GetMrsCalibStrikes() const {return itsPortfolioMode[2];}
	inline const int		GetMrsCalibMode() const { return itsPortfolioMode[3];}

	inline const bool		GetFixBoundaryFlag() const { return itsFixBoundaryFlag;}      
	inline void				SetFixBoundaryFlag(bool flag){ itsFixBoundaryFlag = flag;}      
	inline const bool		GetApproxMarginFlag() const { return itsApproxMarginFlag;}      
	inline void				SetApproxMarginFlag(bool approx) {itsApproxMarginFlag = approx;}  
	inline void				SetFreezeBetasFlag(bool flag) {itsFreezeBetasFlag = flag;}  
	inline const int		GetAtmMrsCalibType() const {return itsAtmMrsCalibType;}
	inline const int		GetVolatilityLag() const { return itsPortfolioMode[4];}

	inline bool				GetCalculateProbaFlag() { return itsCalculateProbaFlag;}
	inline void				SetCalculateProbaFlag(bool flag) { itsCalculateProbaFlag = flag;}
    inline const int		GetFactorNb() const { return ((itsModelType==ARM_PricingModelType::SFRM2F)?2:1);}

	void					SetPreCalibMethod(ARM_CalibMethodPtr calibMethod) {itsPreCalibMethod = calibMethod;}
	ARM_CalibMethodPtr	    GetPreCalibMethod() { return itsPreCalibMethod ;}

	ARM_DateStripPtr		GetCallDateStrip() { return itsCallDateStrip;}

	int						AdjustStubRule(ARM_Date& firstDate, ARM_Date& lastDate);
	inline void				SetStubRule(int stub){itsStubRule = stub;}

	//Deal Description
	ARM_StringVector		PricedColumnNames() const;
	ARM_StringVector		ControlVariateColumnNames() const;
	virtual ARM_RowInfo		ColumnNames() const;
	virtual ARM_RowInfo		MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	void					InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;
	virtual ARM_DateStripCombiner	DatesStructure() const;
	virtual ARM_DateStripCombiner	CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	//Bermuda Swaption constant manager: Strike, Notional, Margin.
	ARM_CstManagerPtr		CreateCstManager();

	//To move
	ARM_Curve*				RefValueToARM_Curve(ARM_ReferenceValue* refVal, ARM_Date& date, ARM_Interpolator<double,double>* interpolator);

	//Get Step-up vectors
	inline ARM_ReferenceValue& GetStrike(){return itsStrike;}
	inline ARM_ReferenceValue& GetNominal(){return itsNominal;}
	inline ARM_ReferenceValue& GetSpread(){return itsVarSpread;}
	inline ARM_ReferenceValue& GetFees(){return itsFees;}
	inline void				   SetVarSpread(ARM_ReferenceValue& varSpread){itsVarSpread = varSpread;}

	//Right rate to propagate
	void SetEquivalentVarIndexTerm();
	
	//Calibration part
	ARM_ReferenceValue* CreateEquivalentNotional(double strike, int FixFreq, ARM_ReferenceValue* initNotionalCurve);
	void CreateEquivalentSwaption(ARM_ExerciseStyle* exerciseDates);
	
	void CreateEmptyCalibration();
	ARM_CalibMethod* GetOSWCalibMethod() const;
    const ARM_StdPortfolioPtr GetOSWPortfolio() const;
	const ARM_StdPortfolioPtr GetStdPortfolio() const;
    void SetOSWPortfolio(const ARM_StdPortfolio& port);
	const ARM_StdPortfolioPtr GetSTMPortfolio() const;
    void SetSTMPortfolio(const ARM_StdPortfolio& port);
	ARM_StdPortfolioPtr CreateOSWPortfolio();
	ARM_StdPortfolioPtr CreateStdOSWPortfolio();
	ARM_StdPortfolioPtr CreateSTMPorftolio(ARM_StdPortfolioPtr port);

	inline void			 SetExoSwaption(ARM_Swaption* exoSwaption){itsExoSwaption = exoSwaption;}
	inline ARM_Swaption* GetExoSwaption(){return itsExoSwaption;}
	inline ARM_Swap*	 GetUnderSwap(){return itsUnderSwap;}
	inline void			 SetUnderSwap(ARM_Swap* underSwap){itsUnderSwap = underSwap;}
//	inline int		 GetNbExercise() { return ((itsExoSwaption->GetExercisingDates())->size());}
	
	ARM_Swaption* CreateForwardSwaption(ARM_Swaption* swaption);
	void GenerateProductDescription();

    const ARM_ModelParam& GetBeta() const;
    void SetBeta(ARM_ModelParam* mrsParam);

	inline string   GetGenType1() const {return itsGenType1;}
	inline string   GetGenType2() const {return itsGenType2;}
	inline string   GetPathOrder() const {return itsPathOrder;}
	inline string   GetPathScheme() const {return itsPathScheme;}
	inline int		GetFirstNbTimes() const {return itsFirstNbTimes;}
	inline void		SetGenType(string type) {itsGenType1 = type;}
	inline int		GetModelType() const {return itsModelType;}
	inline void		SetModelType(int type) {itsModelType = type;}

	//Calibration & Model
	virtual void UpdateCalibration(bool isUpdateStrike=true);
	virtual void CreateAndSetModel();
	virtual void CreateAndSetCalibration();
    virtual void UpdateModel();
	virtual void Calibrate();
	virtual void CreateAndSetCalibration_Basket();
	virtual void CreateAndSetCalibration_Frontier();
	ARM_StdPortfolioPtr SwaptionPortfolio_Frontier();
	
	//Mean Reversion Finder
	double MeanReversionRoot(double targetPrice, double fOrXTolerance, int maxIter);

	//Outputs
	virtual void				ComputePricingData() const;

	virtual double				Price();
	inline double				GetBermudaSwaptionPrice(){return itsBermudaSwaptionPrice;}
	inline void					SetBermudaSwaptionPrice(double price){itsBermudaSwaptionPrice = price;}
	inline double				GetBermudaSwaptionStdDev(){return itsBermudaSwaptionStdDev;}
	inline void					SetBermudaSwaptionStdDev(double stdDev){itsBermudaSwaptionStdDev = stdDev;}

	inline 	double				GetControlVariatePrice(int pi){return (*itsControlVariatePrices)[pi];} 
	inline 	void				SetControlVariatePrice(int i, double pricei){(*itsControlVariatePrices)[i] = pricei;} 
	inline const int			GetNbControlVariate(){return (*itsControlVariate).size();}
	inline const vector<int>*	GetControlVariate() const { return itsControlVariate;}
	inline void					SetControlVariate(vector<int>* ctrlVar){ delete itsControlVariate;
																		 itsControlVariate = ctrlVar;}
	inline void					SetControlVariatePrices(vector<double>* ctrlVar){ delete itsControlVariatePrices;
																				  itsControlVariatePrices = ctrlVar;}

	void						ComputeControlVariate(std::vector<double>& cvprices);
	std::vector<double>				GetCVPrices();
    double						GetCVPrice(int cvi);

	//Some utilities
	virtual string				ExportShortName() const { return "LCBSA";}
	string						toString(const string& indent, const string& nextIndent) const;
	virtual void				View(char* id = NULL, FILE* ficOut = NULL) const;

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



private:

	//Bermuda Swaption general datas
	ARM_Currency				itsCcy;
	ARM_Date					itsPastStartDate;	//Used for ZC startted deals.
	ARM_Date                    itsStartDate;           
	ARM_Date                    itsEndDate;            
	ARM_ReferenceValue			itsNominal;			//Notional Step-Up
    ARM_ReferenceValue			itsStrike;			//Strike Step-Up         
  	int                         itsPayRec;
	int                         itsStubRule;		//Underlying stub rule.
	
	bool						itsZC;				//GenSecType=4
	bool						itsPastZC;			//Used for started bermuda ZC.

	//Call datas
	int							itsCallFreq;
	int							itsCallNotice;		//Positive value
	string                      itsCallCal;
	ARM_Date					itsFirstCallDate;
	ARM_Date					itsLastCallDate;
	CC_IS_MUTABLE ARM_ReferenceValue	itsFees;			//NoticeDate, Associated Fees 

	//Underlying fix leg
	int                         itsFixFreq;			//FixPayFreq
	int                         itsFixBasis;		//A360, A365, ...    
    string						itsFixPayCal;		//Payment calendar
	int							itsFixAdjRule;		//Following, modified following
	int							itsFixRule;			//ADJ, UNADJ
	int							itsFixPayGap;		//Payment days gap
	
   	//Underlying var leg
	int                         itsVarFreq;			//VarFixFreq = VarPayFreq 
	int							itsVarBasis;		//A360, A365 ...
	string						itsVarResetCal;		//Fixing calendar
	string                      itsVarPayCal;		//Payment calendar
	string                      itsVarIndexTerm;	//Fixed Index
 	int							itsVarAdjRule;		//Following, modified following...
	int							itsVarRule;			//ADJ or UNADJ
	int                         itsVarResetGap;		//Fixing Days Gap
	int							itsVarPayGap;		//Payment Days Gap
    ARM_ReferenceValue			itsVarSpread;		//Spread Step-Up

  	mdmKeysAlias				itsModelKey;

	int							itsGenSecType;

	bool						itsHasBeenPriced;

	//Model Parameters
	vector<ARM_ReferenceValue*>*	itsModelParams;		//5) InitVol, InitMrs, InitBeta, InitTheta, InitSkew, InitMrsL, InitMrsU

	bool						itsMrsCalibFlag;	//5) CalibMRS (Y/N question).

	bool						itsATMDiagonalFlag; //5) CalibBETA (Y/N question).

	//Numeric Flags
	int							itsNumMethodType;
	int							itsAmcIter;
	int							itsMcIter;
	int							itsMaxBucketSize;
	int							itsTreeSteps;
	bool						itsFixBoundaryFlag;
	string						itsGenType1;
	string						itsGenType2;
	string						itsPathOrder;
	string						itsPathScheme;
	int							itsFirstNbTimes;
	int							itsModelType;

	//ExoSwaption to construct CalibPortfolio as Summit
	ARM_Swaption*				itsExoSwaption;
	ARM_Swap*					itsUnderSwap;

	//Portfolio Mode
	vector<int>					itsPortfolioMode;	//3) Mode(NEWMODE/SUMMIT), CalibFreq, CalibStrikes, MRSCalibMode

	int							itsAtmMrsCalibType;

	bool						itsApproxMarginFlag;
	bool						itsFreezeBetasFlag;
	bool						itsCalculateProbaFlag;
	
	//// For Probabilities Calculation
	ARM_GP_Matrix				itsProbaMatrix;
	size_t						itsNbCalculatedProbabilities;
	ARM_VectorPtr				itsCumProba;

	//Pricing Data
	double						itsDefaultMrs;		// Default MRS used in the Construction of the Calculator
	double						itsBermudaSwaptionPrice;
	double						itsBermudaSwaptionStdDev;
	ARM_GP_VectorPtr			itsBermudaSwaptionRA;
	vector<double>*				itsControlVariatePrices;
	vector<int>*				itsControlVariate;
	ARM_VectorPtr				itsBetasCoeff;
	
	ARM_CalibMethodPtr			itsPreCalibMethod;	// For Mean Reversion Calibration using a std portfolio

	CC_IS_MUTABLE ARM_DateStripPtr	itsCallDateStrip;
	ARM_DateStripPtr			itsFundDateStrip;
	ARM_DateStripPtr			itsStructDateStrip;
};

CC_END_NAMESPACE()

#endif
