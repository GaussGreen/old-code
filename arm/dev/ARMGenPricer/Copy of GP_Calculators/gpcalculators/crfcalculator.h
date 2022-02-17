/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file crfcalculator.h
 *
 *  \brief file for the CRF Calculator
 *	\author  JM Prie
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPCALCULATORS_CRFCALCULATOR_H
#define _INGPCALCULATORS_CRFCALCULATOR_H

#include "gpbase/port.h"
#include "gencalculator.h"
#include "gpbase/assignop.h"
#include "typedef.h"

#include "gpmodels/typedef.h"
#include "gpbase/countedptr.h"
#include "gpmodels/argconvdefault.h"

/// kernel
#include <util/refvalue.h>


CC_USING_NS(std,pair)

/// forward declaration
class ARM_Portfolio;
class ARM_Swap;
class ARM_Swaption;
class ARM_VolCurve;
class ARM_ZeroCurve;
class ARM_VolLInterpol;
class ARM_BSModel;

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_DateStrip;
class ARM_CRFCalculator;
class ARM_CurveModelParam;


///-----------------------------------------------------------------------------
/// \class CRFStrikeAdjFct
/// \brief
///  Functor class to implement the computation of vanilla price error given by
///  the model of a CRF calculator and the market price
///-----------------------------------------------------------------------------
class CRFStrikeAdjFct
{
private:
    size_t              itsProdIdx;
    ARM_CRFCalculator*  itsCRFCalculator;

public:
    CRFStrikeAdjFct(size_t prodIdx=0, ARM_CRFCalculator* crfCalculator=NULL);
    virtual ~CRFStrikeAdjFct();

    void SetProdIdx(size_t prodIdx) {itsProdIdx=prodIdx;}
    double operator () ( double x ) const;
};



///-----------------------------------------------------------------------------
/// \class ARM_CRFCalculator
/// \brief
///  Class that implements a Callable Reverse Floater calculator
///-----------------------------------------------------------------------------
class ARM_CRFCalculator : public ARM_GenCalculator 
{
    static const string CRFColNamesTable [];

public:

    enum CRFColAlias
    {
        EventDate=0,
        StartDate,
        EndDate,
        PayDate,
        IT,
        IndexStartDate,
        Strike,
        CpnMin,
        CpnMax,
        Leverage,
        AdjCapSpread,
        AdjFloorSpread,
        Nominal,
        FundingSpread,
        FundingNominal,
        FundingStartDate,
        FundingEndDate,
        FundingAnnuity,
        FundingNominalExchange,
        FundingVarFlow,
        FundingSpreadFlow,
        Funding,
        DFPay,
        CpnIndex,
        CpnIndexFlow,
        StrikeFlow,
        StdVarFlow,
        StdFixFlow,
        Caplet,
        Floorlet,
        StdFlow,
        RFFlow,
        StdSwaplet,
        StdSwap,
        RFSwaplet,
        RFSwap,
        Fee,
        StdBermuda,
		FinalDate,
		SwapRate,		
        RFBermuda,
		Frontier,
        StdBermudaPrice,
        RFBermudaPrice
    };

    enum mdmKeysAlias
    {
        YcKey=0,
        OswModelKey,
        CfModelKey,
        MrsKey,
        FundingKey,
        BasisKey,
        ForexKey,

        NbKeys,
    };

    enum modelsAlias
    {
        myRefModel=0,         // the stochastic IR model
        myIrMarginModel,      // the forward margin model to generated the 2nd IR model
        myBasisMarginModel,   // the forward margin model for basis swap
        myForexModel,         // the forex model

        NbModels
    };

	/// constructor, copy constructor, assignment constructor, destructor
	ARM_CRFCalculator(const ARM_Date& startDate,
                const ARM_Date& endDate,
                const ARM_ReferenceValue& strike,
                int payRec,
                const ARM_Date& fixEndDate,
                int fixDayCount,
                int cpnDayCount,
                int cpnFreq,
                int cpnTiming,
                const string& cpnIndexTerm,
                int cpnIndexDayCount,
                const string& cpnResetCal,
                const string& cpnPayCal,
				int stubRule,
                int cpnResetGap,
                const ARM_ReferenceValue& leverage,
                const ARM_ReferenceValue& cpnMin,
                const ARM_ReferenceValue& cpnMax,
                const ARM_ReferenceValue& fundSpread,
                int fundFreq,
                int fundDayCount,
                const ARM_ReferenceValue& nominal,
                int exerGap,
                int nbNonCall,
                const ARM_ReferenceValue& exerFee,
                ARM_SigmaCalibType oswCalibType,
                ARM_MRSCalibType mrscalibType,
				ARM_MRSStrikeCalibType mrsStrikeType,
                bool capCalibFlag,
                bool floorCalibFlag,        
				const ARM_MarketData_ManagerRep& mktDataManager,
                const ARM_StringVector& mdmKeys,
                const ARM_ReferenceValue& fundnominal = ARM_ReferenceValue(),
                ARM_ModelType modelType = ARM_PricingModelType::HWM1F,
				bool skewAutoCalFlag = true,
				long isFrontier = 0);

    /// SUMMIT version
	ARM_CRFCalculator(const ARM_Date& asOfDate,
                const ARM_Date& startDate,
                const ARM_Date& endDate,
                const ARM_ReferenceValue& strike,
                int payRec,
                const ARM_Date& fixEndDate,
                int fixDayCount,
                int cpnDayCount,
                int cpnFreq,
                int cpnTiming,
                const string& cpnIndexTerm,
                int cpnIndexDayCount,
                const string& cpnResetCal,
                const string& cpnPayCal,
				int stubRule,
                int cpnResetGap,
                const ARM_ReferenceValue& leverage,
                const ARM_ReferenceValue& cpnMin,
                const ARM_ReferenceValue& cpnMax,
                const ARM_ReferenceValue& fundSpread,
                int fundFreq,
                int fundDayCount,
                const ARM_ReferenceValue& nominal,
                const ARM_ReferenceValue& exerFee,
                const ARM_Currency& cpnCcy,
                const ARM_Currency& fundingCcy,
                const ARM_Currency& basisCcy,
                const ARM_Currency& domesticCcy=ARM_Currency(),
                const ARM_Currency& foreignCcy=ARM_Currency(),
                const ARM_ReferenceValue& fundnominal = ARM_ReferenceValue(),
				ARM_ModelType modelType = ARM_PricingModelType::HWM1F,
                bool SkewReCalFlag = false);

    /// SUMMIT version with customized schedule
	ARM_CRFCalculator(const ARM_DateStrip& fixLegSched,
					  const ARM_DateStrip& RFLegSched,
					  const ARM_DateStrip& ExerSched,
					  const ARM_Date& asOfDate,
                const ARM_ReferenceValue& strike,
                int payRec,
                const ARM_Date& fixEndDate,
                int fixDayCount,
                int cpnDayCount,
                int cpnFreq,
                int cpnTiming,
                const string& cpnIndexTerm,
                int cpnIndexDayCount,
                const string& cpnResetCal,
                const string& cpnPayCal,
				int stubRule,
                int cpnResetGap,
                const ARM_ReferenceValue& leverage,
                const ARM_ReferenceValue& cpnMin,
                const ARM_ReferenceValue& cpnMax,
                const ARM_ReferenceValue& fundSpread,
                int fundFreq,
                int fundDayCount,
                const ARM_ReferenceValue& nominal,
                const ARM_ReferenceValue& exerFee,
                const ARM_Currency& cpnCcy,
                const ARM_Currency& fundingCcy,
                const ARM_Currency& basisCcy,
                const ARM_Currency& domesticCcy=ARM_Currency(),
                const ARM_Currency& foreignCcy=ARM_Currency(),
                const ARM_ReferenceValue& fundnominal = ARM_ReferenceValue(),
				ARM_ModelType modelType = ARM_PricingModelType::HWM1F);

	ARM_CRFCalculator( const ARM_CRFCalculator& rhs );
	ASSIGN_OPERATOR(ARM_CRFCalculator)
	~ARM_CRFCalculator();

    /// Get & Set MRS param for SUMMIT interface (deal dependent data)
    const ARM_CurveModelParam& GetMRS() const;
    void SetMRS(ARM_ModelParam* mrsParam);
	void SetMRS(ARM_ReferenceValue* mrsValues);
	void SetMRS(double mrsValue);
	double InterpolMeanRevParam(ARM_ReferenceValue* mrsParam);

	void InitCRFForSummit(ARM_ZeroCurve* zcCpn, 
		                  ARM_VolCurve* swoptVC, 
						  ARM_VolCurve* capVC, 
						  ARM_VolLInterpol* capRho = NULL, 
						  ARM_VolLInterpol* capNu = NULL,
						  ARM_VolLInterpol* capBeta = NULL,
						  ARM_VolLInterpol* swaptRho = NULL,
						  ARM_VolLInterpol* swaptNu = NULL,
						  ARM_VolLInterpol* swaptBeta = NULL,
		                  ARM_ZeroCurve* zcFund   = NULL, 
		                  ARM_ZeroCurve* zcCpnBasis  = NULL, 
		                  ARM_ZeroCurve* zcFundBasis  = NULL, 
						  double forex = 0.0,
                          int UPDATE = 0,
                          long SkewReCalFlag   = -1,
                          int SABRSigmaOrAlpha = 1);

    /// For SUMMIT because forex isn't known at CRF construction time
    void UpdateWithForex(double spotValue);

    ARM_CalibMethod* GetOSWCalibMethod() const;
    const ARM_StdPortfolioPtr GetOSWPortfolio() const;
    void SetOSWPortfolio(const ARM_StdPortfolio& port);
	
    ARM_CalibMethod* GetSTMCalibMethod() const;
    const ARM_StdPortfolioPtr GetSTMPortfolio() const;
    void SetSTMPortfolio(const ARM_StdPortfolio& port);

    const vector< ARM_VanillaArg* >& GetVanillaArgVect(void) const {return itsVanillaArgVect;}
    const ARM_StdPortfolio* GetCFPortfolio(void) const;
    void SetCFPortfolio(const ARM_StdPortfolio& port);

	const ARM_StdPortfolio* GetSkewPortfolio(void) const;

	const ARM_StdPortfolio* GetOSWSkewPortfolio(void) const;
	void SetOSWSkewPortfolio(const ARM_StdPortfolio& portfolio);
	
	/// To pre-Initialize and pre-Calibrate  QGM Skew
	void PreCalibrateSkew();

    inline void SetModelType(const ARM_ModelType modelType)				{	itsModelType = modelType;		}	
	void SetModelType(const string& modelTypeStr)						{	itsModelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(modelTypeStr);	}
    inline ARM_ModelType GetModelType() const							{	return itsModelType;			}
	inline void SetMRSStrikeType(const ARM_MRSStrikeCalibType mrsType)	{	itsMRSStrikeType = mrsType;		}
	inline ARM_MRSStrikeCalibType GetMRSStrikeCalibType() const			{	return itsMRSStrikeType;		}
    inline ARM_SigmaCalibType GetOSWCalibFlag() const					{	return itsOSWCalibFlag;			}
    inline void SetOSWCalibFlag(ARM_SigmaCalibType flag)				{	itsOSWCalibFlag=flag;			}
    inline ARM_MRSCalibType GetSTMCalibFlag() const						{	return itsMRSCalibType;			}
	inline void SetSTMCalibFlag(ARM_MRSCalibType mrsType) 				{	itsMRSCalibType = mrsType;		}
    inline bool GetCapCalibFlag() const									{	return itsCapCalibFlag;			}
    inline void SetCapCalibFlag(bool flag)								{	itsCapCalibFlag=flag;			}
    inline bool GetFloorCalibFlag() const								{	return itsFloorCalibFlag;		}
    inline void SetFloorCalibFlag(bool flag)							{	itsFloorCalibFlag=flag;			}
	inline bool GetSkewCalibFlag() const								{	return itsFloorCalibFlag;		}
    inline void SetSkewCalibFlag(bool skewFlag)							{	itsSkewCalFlag=skewFlag;		}
	/// C'est borrin ......
	inline void SetSkewReCalibFlag(bool skewFlag)						{	itsSkewReCalFlag=skewFlag;		}
	inline void SetNbIterFrontier (long isfrontier)                     { itsIsFrontier = isfrontier;       }
	/// Get the B&S models from the market data manager					
	ARM_BSModel* GetOSWBSModel() const;
	void SetCalibParam(ARM_ModelParam* skewParam);
	void SetSTMAndUpdateCalibFlag(ARM_MRSCalibType type);				

    void SetCRFToPrice();
    bool IsCRFToPrice() const;
    void SetCapToPrice();
    bool IsCapToPrice() const;
    void SetFloorToPrice();
    bool IsFloorToPrice() const;
    void SetFundingToPrice();
    bool IsFundingToPrice() const;
    void SetStdLegToPrice();
    bool IsStdLegToPrice() const;
    void SetRFLegToPrice();
    bool IsRFLegToPrice() const;
    void SetStdSwapToPrice();
    bool IsStdSwapToPrice() const;
    void SetRFSwapToPrice();
    bool IsRFSwapToPrice() const;
    void SetBermudaToPrice();
    bool IsBermudaToPrice() const;

	/// Core of the calculator
	virtual void CreateAndSetModel();
	virtual void UpdateModel();
	virtual void CreateAndSetCalibration();
	void CreateAndSetCalibration_Frontier();
	virtual void UpdateCalibration(bool isUpdateStrike=true);
    virtual double Price();
	virtual void ComputePricingData() const;

    /// To get a subGenSecurity
    virtual ARM_GenSecurity* GetSubGenSecurity(int indexType) const ;

    /// CRF deal description creation functions
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;
	virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;
    /// Auto-calibration of strike spreads
    virtual void Calibrate();

    /// Standard ARM support
	virtual ARM_Object* Clone() const;
	virtual void View(char* id = NULL, FILE* ficOut = NULL) const;

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

private:

    /// CRF financial datas
	ARM_Date                    itsStartDate;           // start date
	ARM_Date                    itsEndDate;             // end date

    ARM_ReferenceValue          itsNominal;             // nominal curve
    ARM_ReferenceValue          itsFundNominal;         // funding nominal curve
    ARM_ReferenceValue          itsStrike;              // fix coupon rate or reverse coupon strike cst...
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
	int                         itsStubRule;            // ability to have shortstart .... 
    int                         itsCpnResetGap;         // reverse index reset gap
    ARM_ReferenceValue          itsLeverage;            // reverse index leverage curve
    ARM_ReferenceValue          itsCpnMin;              // reverse coupon floor level curve
    ARM_ReferenceValue          itsCpnMax;              // reverse coupon cap level curve

    ARM_ReferenceValue          itsFundSpread;          // funding spread curve
    int                         itsFundFreq;            // funding frequency
    int                         itsFundDayCount;        // funding day count

    int                         itsExerGap;             // notification lag in days prior to reverse index reset
    string                      itsExerCal;             // notification calendar
	CC_IS_MUTABLE int			itsCurrentTenor;

    /// These attributes are mutable because the SUMMIT constructor inputs
    /// its notification dates. The CRF schedule builder may lead to
    /// additional non callable reset dates and it then will update this curve and compute non call periods
    CC_IS_MUTABLE int                   itsNbNonCall;           // number of non call period (w.r.t. coupon frequency)
    CC_IS_MUTABLE ARM_ReferenceValue    itsExerFee;             // exercise fees curve

    /// Vector of underlying diagonal swaps (used for equivalent strike & vol computations)
    vector< pair<ARM_Swap*,int> >         itsStdSwaps;
    vector< pair<ARM_Swap*,int> >         itsFundStdSwaps;
    vector< pair<ARM_Swap*,int> >         itsIndexStdSwaps;
    ARM_IntVector                         itsOtherStdSwapsIdx;
	ARM_IntVector                         itsKeepOSWIdx;
	ARM_GP_Vector                         itsvStrike;

    /// Vector of implied caplet & floorlet
    vector< ARM_VanillaArg* >   itsVanillaArgVect;  // for GP pricing
    ARM_StdPortfolio*           itsCapFloorPF;      // for ARM target price
	ARM_StdPortfolio*           itsOWSSkewPF;      // for ARM target price

    /// Flag to specify product to price
    int itsProductToPrice;

    /// Flags to enable auto-calibrations
    CC_IS_MUTABLE ARM_SigmaCalibType itsOSWCalibFlag;	/// Volatility bootstrapping
    ARM_MRSCalibType itsMRSCalibType;					/// MRS optimisation
    bool itsCapCalibFlag;								/// Caplet strike spread fitting
    bool itsFloorCalibFlag;								/// Floorlet strike spread fitting
	bool itsSkewCalFlag;								/// Calibrate the Skew
	bool itsSkewReCalFlag;                              /// C'est borrin mais....
	long itsIsFrontier;                                  /// C'est borrin mais....
  

    /// Names of model used in basis case (chosen among MdM keys)
    mdmKeysAlias itsCpnModelKey;
    mdmKeysAlias itsFundingModelKey;
    mdmKeysAlias itsBasisRefModelKey;


    /// To memorise the line in the deal description corresponding to the first call
    /// Mutable because may be updated in const functions
    CC_IS_MUTABLE size_t itsFirstCallIdx;


    /// Convexity model for equivalent strike computation for
    /// diagonal swaption autocalibration (only used in in-arrears case !)
    /// It will a B&S model but not at the moment because the GP doesn't support B&S !
	ARM_PricingModelPtr itsConvexityModel;

    /// Type of model pricing //HWM or QGM
    ARM_ModelType itsModelType;

	/// Type of MRS Calibration
	ARM_MRSStrikeCalibType itsMRSStrikeType;


    /// Set the keys or name alias of model (usefull in multi-currency context)
    void SetModelKeys();


    /// Specialised version for datas consistency
    virtual void CheckData();
	virtual void CheckMktData();


    /// To handle the exercise schedule given by SUMMIT
    void ReplaceExerDates(ARM_DateStrip* schedule,const ARM_DateStrip* cpnSchedule) const;
    void ReplaceExerDatesWithExerGap(ARM_DateStrip* schedule) const;
 
    /// Initialisation to 0 of all columns of the deal description that could be priced
    virtual void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;


    /// Utilities
    ARM_INDEX_TYPE GetIndexType();

    /// Calibration product creation
    int GetNextValidDateRow(const ARM_DealDescription& dealDescr,CRFColAlias columnName,
                            int prevRowIdx,int step,int stopIdx);
    ARM_StdPortfolioPtr CreateDiagonalSwaption();
    ARM_StdPortfolioPtr CreateShortTermVanilla(const ARM_StdPortfolioPtr& diagonalSwaptionPF);
	ARM_StdPortfolio* CreateDiagFWDVanilla(const ARM_StdPortfolio& portfolio);
    void CreateImpliedCapletFloorlet();
	double CreateEquivalentStrike(int startRowIdx, int endRowIdx);

    /// Calibration target prices for hedge purpose
    double GetSubPrice(int startRowIdx,int endRowIdx,CRFColAlias columnName,
                       const string& evalDateStr,
                       const ARM_DealDescription& dealDesc,
                       const string& payModelName, 
                       ARM_PricingModelPtr modelpricing) const;

	void ComputeEquivalentDatas(ARM_Swap* swap,
								ARM_Swaption* swaption,
								ARM_GP_Vector& equivDatas);
    void ComputeEquivalentDatas(const ARM_DealDescription& dealDesc,
                                const string& payModelName,
                                const pair<ARM_Swap*,int>& stdSwap,
                                ARM_Swaption* swaption, 
								const string& evalDateStr,
                                ARM_GP_Vector& equivDatas);

	void ComputeEquivalentStdDatas(
								const ARM_DealDescription& dealDesc,
								const string& payModelName,
								const pair<ARM_Swap*,int>& stdSwap,
								ARM_Swaption* swaption, 
								const string& evalDateStr, 
								ARM_GP_Vector& equivDatas);

    void ComputeEquivalentDatas(const ARM_DealDescription& dealDesc,
                                const string& payModelName,
                                const pair<ARM_Swap*,int>& stdSwap,
                                const pair<ARM_Swap*,int>& fundStdSwap,
                                const pair<ARM_Swap*,int>& indexStdSwap,
                                ARM_Swaption* swaption, 
								const string& evalDateStr,
                                ARM_GP_Vector& equivDatas);

	void ReselectSwaption(ARM_Swaption* swaption,double swapRate,double vol,bool isStd,int idx,ARM_IntVector& removedOSWIdx);

    void ComputeDiagonalSwaptionPrice(bool isFreezeWeights, bool isInitParam, bool isUpdateStrike);
	ARM_CurveModelParam* CheckPortfolioOSW( ARM_IntVector& removedOSWIdx );
    void ComputeShortTermVanillaPrice(bool isFreezeWeights, bool isInitParam, bool isUpdateStrike = true);
	void ComputeDiagFWDVanillaPrice(bool isFreezeWeights, bool isInitParam);
    void ComputeImpliedCapletFloorletPrices(bool isFreezeWeights);

	double CreateFloorletAndDelta(int startRowIdx,int endRowIdx,const string& evalDateStr);
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

