/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file prcscalculator.h
 *
 *  \brief file for the Power Reverse Dual Currencies Calculator
 *	\author  JM.Prie & E.Ezzine
 *	\version 1.0
 *	\date February 2006
 */


#ifndef _INGPCALCULATORS_PRDCALCULATOR_H
#define _INGPCALCULATORS_PRDCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gencalculator.h"

// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"


#include "typedef.h"
CC_USING_NS(std,pair)


/// Kernel forward declaration
class ARM_PowerReverse;
class ARM_DFBSModel;
class ARM_VolLInterpol;
class ARM_StdPortfolio;
class ARM_Option;

CC_BEGIN_NAMESPACE( ARM )

/// GP forward declaration
class ARM_DateStrip;
class ARM_MarketHybridModel;
struct ARM_VanillaIrFxSwaption;
class ARM_BS_Model;


///-----------------------------------------------------------------------------
/// \class ARM_PRDCalculator
/// \brief
///  Class that implements a Power Reverse Dual Currencies calculator
///-----------------------------------------------------------------------------
class ARM_PRDCalculator : public ARM_GenCalculator 
{
protected: 
    static const string PRDColNamesTable [];
    static const string PRDProfileNamesTable [];

public:

    enum PRCSColAlias
    {
        EventDate=0,
        FxResetDate,        
        StartDate,
        FundingEndDate,
		FundingPayDate,
        EndDate,
		CpnPayDate,
		CpnIT,
        CpnLeverage,
		FxLowStrike,
        FxHighStrike,
        FundingLeg,		
		FxLowCall,
		FxHighCall,
        FixedLeg,
        FxLowCallStrip,
        FxHighCallStrip,
        RedemptionResetDate,
        RedemptionPayDate,
        RedemptionStrike,
		Nominal,
        Redemption,
        CpnLeg,
        PRDSwap,
		FundingFlow,
		FixFlow,
        FxStrip,
		PRDCFirstSwap,
        PRDCFirstEuropean,
    };

    enum mdmKeysAlias
    {
        YcDomKey=0,
        YcForKey,
		YcFundKey,
        ForexKey,
		FunForexKey,
        YcBasisDomKey,
        YcBasisForKey,
		YcBasisFundKey,
        OswDomModelKey,
        OswForModelKey,
        FxModelKey,
        CorrelMatrixKey,
        MrsDomKey,
        MrsForKey,
        QFxKey,
        QDomKey,
        QForKey,
        LocalFxModelKey,
        FlooredFxModelKey,
        CappedFxModelKey,
        RedemptionFxModelKey,
		MarketIrModelKey,
	BSFxVol,

        NbKeys
    };
    
    enum PRDCProfileType
    {
        FxDateStripProfile = 0,
		CpnDateStripProfile,
		FundingDateStripProfile,
		NominalProfile,
		FxLowStrikeProfile,
        FxHighStrikeProfile,
        FxLeverageProfile,
        MarginProfile,
		MarginTimesNotional,
		FundLvgeTimesNotional,
		CpnMinTimesNotional,
    };
        
	/// constructor for deal from term sheet
	ARM_PRDCalculator(const ARM_Date& asofDate,
		const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& DomCcy,
		const ARM_Currency& ForCcy,
		const ARM_Currency& FundCcy,
		int cpnDayCount,
		int cpnFreq,
		int FxResetGap,
		int stubRule,
		int resetTiming,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& domesticCpn,
		const ARM_Curve& foreignCpn,
		const ARM_Curve& MinCpn,
		const ARM_Curve& MaxCpn,
		const ARM_Curve& initialFx,
		int fundFreq,
		int fundDayCount,
		int compFreq,
		int compType,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,		
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees,
		int redemptionGap,
		double redemptionStrike,
		const ARM_RedemptionType& redemptionType    = ARM_PRCSRedemptionType::standard,		
		const ARM_StringVector& productsToPrice		= ARM_StringVector(1,PRDColNamesTable[PRDCFirstEuropean]),
		bool fxLocalModelFlag						= true,
		bool basisIRCalibFlag					    = true,
		const ARM_Curve& fundlevrage                = ARM_FlatCurve(1));

	void Init(const ARM_MarketData_ManagerRep& mktDataManager,
		const ARM_GP_Vector& schedulerDatas = ARM_GP_Vector(0),
		const ARM_GP_Vector& truncatorDatas = ARM_GP_Vector(0),
		bool markovianDriftSamplerFlag = true,
		ARM_PRDCCalibType calibType = ARM_PRCSCalibTypes::ATMCalib,
		const ARM_GP_Vector& calibDatas = ARM_GP_Vector(0));

	/// Constructor for trade id from Summit
	ARM_PRDCalculator(ARM_PowerReverse* powRev,
		ARM_DFBSModel* model,
		const ARM_ObjectVector& otherMktDatas,
		const ARM_GP_Vector& schedulerDatas = ARM_GP_Vector(0),
		const ARM_GP_Vector& truncatorDatas = ARM_GP_Vector(0),
		const ARM_StringVector& columnsToPrice=ARM_StringVector(1,PRDColNamesTable[PRDCFirstEuropean]),
		bool markovianDriftSamplerFlag = true,
		bool fxLocalModelFlag = false,
		ARM_PRDCCalibType calibType = ARM_PRCSCalibTypes::ATMCalib, 
		const ARM_GP_Vector& fxATSCalibDatas = ARM_GP_Vector(0),
		bool basisIRCalibFlag = true,
		ARM_BasisType basisType = ARM_PRCSBasisType::flowByflow);


	///copy constructor, assignment constructor, destructor
	ARM_PRDCalculator( const ARM_PRDCalculator& rhs );
	ASSIGN_OPERATOR(ARM_PRDCalculator)
	~ARM_PRDCalculator();

    /// Portfolio accessors
    const ARM_StdPortfolioPtr GetOSWPortfolio(mdmKeysAlias oswModelKey) const;
    const ARM_StdPortfolioPtr GetFxPortfolio() const;
	void SetOSWPortfolio(const ARM_StdPortfolio& pf, mdmKeysAlias oswModelKey);
    void SetFxPortfolio(const ARM_StdPortfolio& pf);

    const ARM_StdPortfolioPtr GetExtraFxPortfolio() const;
    const ARM_StdPortfolio* const GetFlooredFxPortfolio() const { return itsFlooredFxOptionPF; }
    const ARM_StdPortfolio* const GetCappedFxPortfolio() const { return itsCappedFxOptionPF; }
    const ARM_StdPortfolio* const GetRedemptionFxPortfolio() const { return itsRedemptionFxOptionPF; }

	inline ARM_DFBSModel* GetAnalyticalModel()  const { return itsAnalyticalModel;}
	inline ARM_PowerReverse* GetPowerReverseSwap()  const  {  return itsPowerReverseSwap;}
	
    /// CalibMethod accessors
    ARM_CalibMethod* GetFxCalib() const;
    void SetFxCalib(ARM_CalibMethod* fxCalib) const;
	const ARM_CalibMethod* GetExtraFxCalib() const;

    /// Initialisation to 0 of all columns of the deal description that could be priced
    void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	/// Core of the calculator
	ARM_BS_Model* CreateConvAdjstModel();
	virtual void CreateAndSetModel();
	virtual void CreateAndSetCalibration();
	virtual void UpdateModel();
	virtual void UpdateCalibration(bool isUpdateStrike=true);
    virtual void Calibrate();
    virtual double Price();
    virtual void CheckData();
    virtual void CheckMktData();
    virtual void ComputePricingData() const;

    /// Forward FX vol computation of PRDC FX option strip
    ARM_VolLInterpol* ComputeATMFxOptionVols();

	/// to calculate fx Volatility at the strike/barrier
	virtual void UpdateFxOption(size_t eventIdx,ARM_Option*& fxOption) {ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + 
		" : GetStrikeAtATSFx is unimplemented for PRDC");}

	/// Underlying computing
	ARM_Vector* ComputeAll();

    /// PRDC deal description creation functions
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;
	virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

    /// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_PRDCalculator(*this); }

	virtual string ExportShortName() const { return "LPRDC";}

	/// Accessors
	inline ARM_IntVector GetvCpnIndex() const { return itsvCpnIndex;}
	inline size_t GetNbNoCall() const { return itsNbNoCall;}

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

	/// PRDC financial datas
	ARM_Date			itsStartDate;           // start date
	ARM_Date			itsFixedEndDate;        // fixed end date
	ARM_Date			itsEndDate;             // end date
	ARM_Date            itsRefDate;             // the first endDate
						
	ARM_Currency		itsDomesticCcy;         // domestic Currency
	ARM_Currency		itsForeignCcy;			// foreign Currency
	ARM_Currency		itsFundingCcy;			// funding Currency
						
	int					itsCpnDaycount;         // reverse dual coupon day count
    int					itsCpnFreq;             // reset & pay frequency 
	int					itsFxResetGap;
    string				itsCpnResetCal;         // reverse dual reset calendar
    string				itsCpnPayCal;           // reverse dual payment calendar
	int					itsStubRule;            // ability to have shortstart .... 	
	int                 itsResetTiming;
						
	int					itsFundFreq;
	int					itsFundDaycount;
	int                 itsCompFreq;
	int                 itsCompType;
						
	int					itsExerciseFreq;
	int					itsNoticeGap;
	int					itsPayRec;						
	size_t				itsNbNoCall;
    ARM_Curve			itsFees;

	ARM_Curve			itsCpnnominalCv;
	ARM_Curve			itsDomesticCpnCv;
	ARM_Curve			itsForeignCpnCv;
	ARM_Curve			itsMinCpnCv;
	ARM_Curve			itsMaxCpnCv;
	ARM_Curve			itsInitialFxCv;
	ARM_Curve			itsFundnominalCv;
	ARM_Curve			itsFundSpreadCv;
	ARM_Curve			itsFundlevrageCv;
						
	int					itsNbFixFlows;   // to store the number of  foxed period
	/// Redemption features
	ARM_RedemptionType	itsRedemptionType;
	double				itsRedemptionStrike;
	int					itsRedemptionGap;
	ARM_Date            itsResetRedemptionDate; 

	/// Used only for underlying computing via Kernel
	ARM_GP_Vector   itsDomesticCpn;
	ARM_GP_Vector   itsForeignCpn;
	ARM_GP_Vector   itsMinCpn;
	ARM_GP_Vector   itsMaxCpn;
	ARM_GP_Vector   itsInitialFx;

	//// --------------------------------------------------
	/// ---- Vector conversions --------------------------
	/// --------------------------------------------------
	/// Convert ARM_Curves into relevant vectors
	/// this will be much easier for computations
	ARM_GP_Vector	itsvFundSpread;
	ARM_GP_Vector   itsvFundLeverage;
	ARM_GP_Vector   itsvInitialFundSpread; /// to store the initial curve for hedging
	ARM_GP_Vector	itsvFundNominal;
	ARM_GP_Vector   itsvInitialFundNominal;
	size_t			itsFundSize;
	ARM_IntVector	itsvFundIndex;
	ARM_BasisType   itsBasisType;
	ARM_Curve       itsfundMargin;
	
	ARM_GP_Vector	itsvCpnNominal;
	ARM_GP_Vector	itsvLeverage;
	ARM_GP_Vector   itsvFixCpn;
	ARM_GP_Vector	itsvLowStrikeCpn;
    ARM_GP_Vector	itsvHighStrikeCpn;
	ARM_BoolVector	itsvCpnIsCapped;	
	ARM_BoolVector	itsvCpnIsFloored;	
	ARM_BoolVector	itsvCpnIsFixed;	    
	size_t			itsCpnSize;
	ARM_IntVector	itsvCpnIndex;
	
	ARM_BoolVector	itsvIsExerDate;
	size_t			itsExerSize;

	/// Dates Strip
	ARM_DateStripPtr itsFundDateStripFrom;
	ARM_DateStripPtr itsFundDateStrip;
	ARM_DateStripPtr itsStructDateStrip;
	ARM_DateStripPtr itsForexDateStrip;
	ARM_DateStripPtr itsExerciseDateStrip;

    /// Scheduler and truncator parameter default 
	/// values may be changed from outside
    ARM_GP_Vector       itsSchedulerDatas;
    ARM_GP_Vector       itsTruncatorDatas;

    /// Flags for simpler models
    bool itsDomModelHWFlag;         /// Domestic IR model is degenerated in pure H&W version
    bool itsForModelHWFlag;         /// Foreign IR model is degenerated in pure H&W version

    /// Column names to get multi-prices
    ARM_StringVector itsColumnsToPrice;

    /// Flag to switch between MarkovianDriftSampler and DriftedMeanRevertingSampler
    bool itsMarkovianDriftSamplerFlag;

    /// Portfolios of Fx option for floored & capped coupons
    ARM_StdPortfolio*   itsFlooredFxOptionPF;
    ARM_StdPortfolio*   itsCappedFxOptionPF;
    ARM_StdPortfolio*   itsRedemptionFxOptionPF;

    /// Flag to allow local model for FX calibration & pricing
    bool itsFxLocalModelFlag;

    /// Datas to manage calibration
    ARM_PRDCCalibType itsCalibType;
    ARM_GP_Vector itsCalibDatas;

    /// Flag to use or not basis curve in swaption calibration
    bool itsBasisIRCalibFlag;

    /// Additional calib method using ATM Fx options of the PRD leg
    ARM_CalibMethod*   itsATMDoubleCalibMethod;

    /// ATM foward FX vol computed from boostrapped spot Fx and Zc vols
    ARM_VolLInterpol* itsATMFwdFxVol;

	/// Market hybrid model & portfolio for new FX calibration
	ARM_MarketHybridModel*  itsMktHybridModel;
    ARM_CalibMethod*		itsHybridBasketCalibMethod;

	/// Analytical model to price the underlying
	ARM_DFBSModel* itsAnalyticalModel;
	ARM_PowerReverse* itsPowerReverseSwap;
	ARM_SwapLeg* itsFundingLeg;

	ARM_GP_Vector     itsFxResetTimeBefore;

	/// to keep old interface by using bps
	bool itsTradeFromDataBase;

	/// to avoid de re-calculate the option price
	bool itsHasBeenComputed;

	/// Convert ARM_Curves into relevant vectors
	void ComputeProductVectorsFromCurves();

	/// Create in internal a generic BSModel to price multicurrencies swap
	ARM_DFBSModel* CreateAnalyticalModel(mdmKeysAlias key=FxModelKey,bool useMktDataMgrCorrels=false);
	/// Create the underlyind power reverse swap
	void CreateUnderlying(const ARM_Date& fisDate);
	void CreateFundingLeg();

	/// Interface with existing SUMMIT/Kernel structures
    void FillMarketDataManager(const ARM_ObjectVector& mktDatas);
	ARM_CstManagerPtr ComputeCstManager();

	/// Get schedule and market data
	void ComputeBoosterDatas();
	void CreateDataFromUnderlying();

	/// Convert le basis from funding Ccy to domestic Ccy
	void ComputeDomesticBasis();


    /// Calibration product creation and target pricing
    pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr > CreateDiagonalSwaption();
    ARM_CalibMethod* CreateIRCalibration(const ARM_StdPortfolioPtr& diagonalSwaptionPF,
										ARM_Currency* ccy, 
										int modelIdx);

    double ComputeFxEquivStrike(size_t eventIdx);
    virtual ARM_StdPortfolioPtr CreateFxOption();
    ARM_CalibMethod* CreateFxCalibration(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx);

	ARM_StdPortfolioPtr CreateHybridBasketOption();
	void ComputeHybridBasketPrices();
	void ComputeHybridBasketCorrelations(double evalTime,
		const vector< ARM_VanillaIrFxSwaption* >& hybridBaskets);

	string HybridBasketDump();

    void CreateLocalFxOption();

	/// compute the prices for calibMethod
    void ComputeIROptionPrices(ARM_CalibMethod* calibMethod,
								mdmKeysAlias oswModelIdx, 
								bool isModelHW,
								bool isFreezeWeights, 
								bool isInitVolParam);
    void ComputeFxOptionPrices(ARM_CalibMethod* calibMethod, 
								bool isFreezeWeights, 
								bool isInitVolParam);
    void ComputeLocalFxOptionPrices();

    /// Specific calibration for local Fx model
    void CalibrateLocalFxModel();

    /// Calibration errors checking
    void CheckCalibErrors();
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

