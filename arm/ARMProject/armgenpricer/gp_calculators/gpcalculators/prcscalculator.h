/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file prdccalculator.h
 *
 *  \brief file for the Power Reverse Dual Currencies Calculator
 *	\author  J-M Prié
 *	\version 1.0
 *	\date February 2005
 */


#ifndef _INGPCALCULATORS_PRDCCALCULATOR_H
#define _INGPCALCULATORS_PRDCCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gencalculator.h"

// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"


CC_USING_NS(std,pair)


/// Kernel forward declaration
class ARM_PowerReverse;
class ARM_DFBSModel;
class ARM_VolLInterpol;
class ARM_StdPortfolio;

CC_BEGIN_NAMESPACE( ARM )

/// GP forward declaration
class ARM_DateStrip;
class ARM_MarketHybridModel;
struct ARM_VanillaIrFxSwaption;


///-----------------------------------------------------------------------------
/// \class ARM_PRCSCalculator
/// \brief
///  Class that implements a Power Reverse Dual Currencies calculator
///-----------------------------------------------------------------------------
class ARM_PRCSCalculator : public ARM_GenCalculator 
{
    static const string PRDCColNamesTable [];
    static const string ExtPRDCColNamesTable [];
    static const string PRDCProfileNamesTable [];

public:

    enum PRDCColAlias
    {
        EventDate=0,
        FxResetDate,
        FxLowStrike,
        FxHighStrike,
        CpnIT,
        CpnPayDate,
        FundingStartDate,
        FundingEndDate,
        FundingLastEndDate,
        FundingPayDate,
        FundingSpread,
        FundingIT,
        FundingFlow,
        Funding,
        CpnLeverage,
        FixedCpn,
        FxLowCall,
        FxHighCall,
        FxLowCallStrip,
        FxHighCallStrip,
        RedemptionResetDate,
        RedemptionPayDate,
        RedemptionStrike,
        Redemption,
        CpnFlow,
        Cpn,
        PRDCFlow,
        PRDCSwap,
        PRDCFirstSwap,
        PRDCFirstEuropean,
        PRDCBermuda,
        PRDCOption,
    };

    enum ExtPRDCColAlias
    {
        FundingSum=0,
        CouponSum,
        FundingSum1,
        CouponSum1,
        FundingSum2,
        CouponSum2,
        ExplicitFlow
    };

    enum mdmKeysAlias
    {
        YcDomKey=0,
        YcForKey,
        ForexKey,
        YcBasisDomKey,
        YcBasisForKey,
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

        NbKeys
    };

    enum PRDCCalibType
    {
        ATMCalib=0,
        ATSFxCalib,
        ATSFxMixedCalib,
        ATSFxProfileCalib,
        ATSFxMoneynessCalib,
        ATSFxShiftedCalib,
        ATSFxMinVolCalib,
        ATMDoubleCalib,
        ATSFxEquivCalib,
		HybridBasketCalib
    };
    
    enum PRDCProfileType
    {
        FxDateStripProfile=0,
        FxLowStrikeProfile,
        FxHighStrikeProfile,
        FxLeverageProfile,
        FundingDateStripProfile,
        FundingSpreadProfile
    };
        
    /// constructor, copy constructor, assignment constructor, destructor
	ARM_PRCSCalculator(ARM_PowerReverse* powRev, ARM_Model* model,
                const ARM_ObjectVector& otherMktDatas,
                const ARM_GP_Vector& schedulerDatas = ARM_GP_Vector(0),
                const ARM_GP_Vector& truncatorDatas = ARM_GP_Vector(0),
                const ARM_StringVector& columnsToPrice=ARM_StringVector(1,PRDCColNamesTable[PRDCOption]),
                bool markovianDriftSamplerFlag = true,
                bool fxLocalModelFlag = false,
                PRDCCalibType calibType = ATMCalib, const ARM_GP_Vector& fxATSCalibDatas = ARM_GP_Vector(0),
                bool basisIRCalibFlag = true);

	ARM_PRCSCalculator( const ARM_PRCSCalculator& rhs );
	ASSIGN_OPERATOR(ARM_PRCSCalculator)
	~ARM_PRCSCalculator();

    /// Portfolio accessors
    const ARM_StdPortfolioPtr GetOSWPortfolio(mdmKeysAlias oswModelKey) const;
    const ARM_StdPortfolioPtr GetFxPortfolio() const;
    const ARM_StdPortfolioPtr GetExtraFxPortfolio() const;
    const ARM_StdPortfolio* const GetFlooredFxPortfolio() const { return itsFlooredFxOptionPF; }
    const ARM_StdPortfolio* const GetCappedFxPortfolio() const { return itsCappedFxOptionPF; }
    const ARM_StdPortfolio* const GetRedemptionFxPortfolio() const { return itsRedemptionFxOptionPF; }

	inline ARM_DFBSModel* GetAnalyticalModel()  const { return itsInputModel;}
	inline ARM_PowerReverse* GetPowerReverseSwap()  const  {  return itsInputPowRev;}

    /// CalibMethod accessors
    ARM_CalibMethod* GetFxCalib() const;
    void SetFxCalib(ARM_CalibMethod* fxCalib) const;
	const ARM_CalibMethod* GetExtraFxCalib() const;

    /// Initialisation to 0 of all columns of the deal description that could be priced
    void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	/// Core of the calculator
	virtual void CreateAndSetModel();
	virtual void CreateAndSetCalibration();
	virtual void UpdateModel();
	virtual void UpdateCalibration(bool isUpdateStrike=true);
    virtual void Calibrate() { CalibrateLocalFxModel(); }
    virtual double Price();
	virtual ARM_Vector* ComputeAll();
    virtual void CheckData();
    virtual void CheckMktData();
    virtual void ComputePricingData() const;// {};

    /// Forward FX vol computation of PRDC FX option strip
    ARM_VolLInterpol* ComputeATMFxOptionVols();


    /// PRDC deal description creation functions
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;
	virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

    /// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_PRCSCalculator(*this); }

	virtual string ExportShortName() const { return "LPRDC";}

private:

    /// Summit input
    ARM_PowerReverse* itsInputPowRev;
    ARM_DFBSModel*    itsInputModel;

    ARM_GP_Vector     itsFxResetTimeBefore;

    /// PRDC financial datas
    ARM_GP_Matrix       itsBoosterDatas;
    ARM_GP_Vector       itsRedemptionDatas;

    /// To say if an eventDate is an actual notice date
    ARM_BoolVector      itsExerciseFlag;

    /// Scheduler and truncator parameter default values may be changed
    /// from outside
    ARM_GP_Vector       itsSchedulerDatas;
    ARM_GP_Vector       itsTruncatorDatas;


    /// Flags for simpler models
    bool itsDomModelHWFlag;         /// Domestic IR model is degenerated in pure H&W version
    bool itsForModelHWFlag;         /// Foreign IR model is degenerated in pure H&W version

    /// To memorise the first event to appear in the deal description (event > asOf)
    /// Mutable because may be updated in const functions
    CC_IS_MUTABLE size_t itsFirstFutureEventIdx;


    /// Column names to get multi-prices
    ARM_StringVector itsColumnsToPrice;
	bool itsHasBeenPriced;

    /// Flag to switch between MarkovianDriftSampler and DriftedMeanRevertingSampler
    bool itsMarkovianDriftSamplerFlag;

    /// Portfolios of Fx option for floored & capped coupons
    ARM_StdPortfolio*   itsFlooredFxOptionPF;
    ARM_StdPortfolio*   itsCappedFxOptionPF;
    ARM_StdPortfolio*   itsRedemptionFxOptionPF;

    /// Flag to allow local model for FX calibration & pricing
    bool itsFxLocalModelFlag;

    /// Datas to manage calibration
    PRDCCalibType itsCalibType;
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


    //void CleanUp();
    //void CopyNoCleanUp(const ARM_PRCSCalculator& rhs);

    /// Interface with existing SUMMIT/Kernel structures
    void FillMarketDataManager(const ARM_ObjectVector& mktDatas);
    void ComputeBoosterDatas();
    ARM_CstManagerPtr ComputeCstManager();

    /// Calibration product creation and target pricing
    pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr > CreateDiagonalSwaption();
    ARM_CalibMethod* CreateIRCalibration(const ARM_StdPortfolioPtr& diagonalSwaptionPF,
            ARM_Currency* ccy, int modelIdx);

    double ComputeFxEquivStrike(size_t eventIdx);
    ARM_StdPortfolioPtr CreateFxOption();
    ARM_CalibMethod* CreateFxCalibration(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx);

	ARM_StdPortfolioPtr CreateHybridBasketOption();
	void ComputeHybridBasketPrices();
	void ComputeHybridBasketCorrelations(double evalTime,const vector< ARM_VanillaIrFxSwaption* >& hybridBaskets);
	string HybridBasketDump();

    void CreateLocalFxOption();

    void ComputeIROptionPrices(ARM_CalibMethod* calibMethod, mdmKeysAlias oswModelIdx, bool isModelHW, bool isFreezeWeights, bool isInitVolParam, bool isDiagCalib=true);
    void ComputeFxOptionPrices(ARM_CalibMethod* calibMethod, bool isFreezeWeights, bool isInitVolParam);
    void ComputeLocalFxOptionPrices();

    /// Specific calibration for local Fx model
    void CalibrateLocalFxModel();

    /// Calibration errors checking
    void CheckCalibErrors(ARM_CalibMethod* theCalib=NULL, size_t theCalibIdx=0);
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

