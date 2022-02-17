/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ccscalculator.h
 *
 *  \brief file for the Callable cross currency swaption
 *	\author  K.Belkheir & E.Ezzine
 *	\version 1.0
 *	\date January 2007
 */


#ifndef _INGPCALCULATORS_CCSCALCULATOR_H
#define _INGPCALCULATORS_CCSCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gencalculator.h"

// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"


#include "typedef.h"
CC_USING_NS(std,pair)


/// Kernel forward declaration
//class ARM_StdPortfolio;
//class ARM_Option;

CC_BEGIN_NAMESPACE( ARM )

/// GP forward declaration
//class ARM_DateStrip;

///-----------------------------------------------------------------------------
/// \class ARM_CCSCalculator
/// \brief
///  Class that implements a Power Reverse Dual Currencies calculator
///-----------------------------------------------------------------------------
class ARM_CCSCalculator : public ARM_GenCalculator 
{
protected: 
    static const string CCSColNamesTable [];
    static const string CCSProfileNamesTable [];

public:

    enum CCSColAlias
    {
        EventDate=0,		
        EndDate,
		Spot,		
		ForStartDate,
		ForEndDate,
		ForeignLeg,
		ForNominal,		
		ForExchangeNotional,
		DomStartDate,
		DomEndDate,
		DomesticLeg,
		DomNominal,		
		DomExchangeNotional,
		ForeignFlow,
		DomesticFlow,
		CCSSwapFlow,
		BasisCCSSwap,
		CCSwap,
		FirstCCSSwap,        
		CCSBermuda,
        CCSOption
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
        QDomKey,
		QForKey,
        QFxKey,

        NbKeys
    };
    
    enum CCSProfileType
    {
        FxDateStripProfile = 0,
		DomDateStripProfile,
		ForDateStripProfile,
		DomNominalProfile,
		ForNominalProfile,
		DomMarginProfile,
		ForMarginProfile,
		DomMarginTimesNominalProfile,
		ForMarginTimesNominalProfile,
    };
        
	/// constructor for deal from term sheet
	ARM_CCSCalculator(const ARM_Date& asOfDate,
		const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& domCcy,
		const ARM_Currency& forCcy,
		int domDayCount,
		int domFreq,
		const string& domResetCal,
		const string& domPayCal,
		int forFreq,
		int forDayCount,
		const string& forResetCal,
		const string& forPayCal,
		int FxResetGap,
		int stubRule,	
		const ARM_Curve& domMargin,	
		const ARM_Curve& domNominal,		
		const ARM_Curve& forMargin,
		const ARM_Curve& forNominal,				
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNoCall,
		const ARM_Curve& fees,
		const ARM_StringVector& productsToPrice	= ARM_StringVector(1,CCSColNamesTable[CCSOption]));


	///copy constructor, assignment constructor, destructor
	ARM_CCSCalculator( const ARM_CCSCalculator& rhs );
	ASSIGN_OPERATOR(ARM_CCSCalculator)
	~ARM_CCSCalculator(){};
	
	void Init(	const ARM_MarketData_ManagerRep& mktDataManager,
		const std::vector<double>& schedulerDatas			= std::vector<double>(0),
		const std::vector<double>& truncatorDatas			= std::vector<double>(0),
		bool markovianDriftSamplerFlag				= true,
		ARM_PRDCCalibType calibType					= ARM_PRCSCalibTypes::ATMCalib, 
		const std::vector<double>& fxATSCalibDatas		= std::vector<double>(0));

    /// Portfolio accessors
    const ARM_StdPortfolioPtr GetOSWPortfolio(mdmKeysAlias oswModelKey) const;
    const ARM_StdPortfolioPtr GetFxPortfolio() const;

    /// CalibMethod accessors
    ARM_CalibMethod* GetFxCalib() const;
    void SetFxCalib(ARM_CalibMethod* fxCalib) const;
	const ARM_CalibMethod* GetExtraFxCalib() const;

    /// Initialisation to 0 of all columns of the deal description that could be priced
    void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	/// Core of the calculator
	virtual void CreateAndSetModel();
	virtual void CreateAndSetCalibration();
	virtual void UpdateModel(){};
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
		" : GetStrikeAtATSFx is unimplemented for CCCS");}

	/// Underlying computing
	ARM_Vector* ComputeAll() { return NULL;};

    /// PRDC deal description creation functions
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

    /// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_CCSCalculator(*this); }

	virtual string ExportShortName() const { return "LCCCS";}

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

	ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips) const;

private:

	/// PRDC financial datas
	ARM_Date			itsStartDate;           // start date
	ARM_Date			itsFixedEndDate;        // fixed end date
	ARM_Date			itsEndDate;             // end date

						
	ARM_Currency		itsDomesticCcy;         // domestic Currency
	ARM_Currency		itsForeignCcy;			// foreign Currency
						
	int					itsDomDaycount;         // reverse dual coupon day count
    int					itsDomFreq;             // reset & pay frequency 
	int					itsFxResetGap;
    string				itsResetCalendar;         // reverse dual reset calendar
    string				itsPayCalendar;           // reverse dual payment calendar
	int					itsStubRule;            // ability to have shortstart .... 						
						
	int					itsForFreq;
	int					itsForDaycount;
	string				itsForResetCalendar;         // reverse dual Foreign reset calendar
    string				itsForPayCalendar;           // reverse dual Foreign payment calendar
						
	int					itsExerciseFreq;
	int					itsNoticeGap;
	int					itsPayRec;
						
	size_t				itsNbNoCall;
    ARM_Curve			itsFees;

	ARM_Curve			itsDomNominalCv;
	ARM_Curve			itsForNominalCv;
	ARM_Curve			itsDomSpreadCv;
	ARM_Curve			itsForSpreadCv;
						
	int					itsNbFixFlows;   // to store the number of  foxed period


	//// --------------------------------------------------
	/// ---- Vector conversions --------------------------
	/// --------------------------------------------------
	/// Convert ARM_Curves into relevant vectors
	/// this will be much easier for computations
	std::vector<double>	itsvDomSpread;
	std::vector<double>	itsvForSpread;
	std::vector<double>	itsvDomNominal;
	std::vector<double>   itsvForNominal;
	size_t			itsForSize;
	ARM_IntVector	itsvForIndex;

	std::vector<double>   itsvFixCpn;
	ARM_BoolVector	itsvCpnIsFixed;	    
	size_t			itsDomSize;
	ARM_IntVector	itsvDomIndex;
	
	ARM_BoolVector	itsvIsExerDate;
	size_t			itsExerSize;

	/// Dates Strip
	ARM_DateStripPtr itsDomDateStrip;
	ARM_DateStripPtr itsForDateStrip;
	ARM_DateStripPtr itsForexDateStrip;
	ARM_DateStripPtr itsExerciseDateStrip;

    /// Scheduler and truncator parameter default 
	/// values may be changed from outside
    std::vector<double>    itsSchedulerDatas;
    std::vector<double>    itsTruncatorDatas;


    /// Column names to get multi-prices
    ARM_StringVector itsColumnsToPrice;

    /// Flag to switch between MarkovianDriftSampler and DriftedMeanRevertingSampler
    bool itsMarkovianDriftSamplerFlag;

    /// Datas to manage calibration
    ARM_PRDCCalibType itsCalibType;
    std::vector<double> itsCalibDatas;

	/// to avoid de re-calculate the option price
	bool itsHasBeenComputed;

	/// Convert ARM_Curves into relevant vectors
	void ComputeProductVectorsFromCurves();

	// Cst Manager
	ARM_CstManagerPtr ComputeCstManager();

    /// Calibration product creation and target pricing
    pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr > CreateDiagonalSwaption();
    ARM_CalibMethod* CreateIRCalibration(const ARM_StdPortfolioPtr& diagonalSwaptionPF,
										ARM_Currency* ccy, 
										int modelIdx);

    //double ComputeFxEquivStrike(size_t eventIdx);
    virtual ARM_StdPortfolioPtr CreateFxOption();
    ARM_CalibMethod* CreateFxCalibration(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx);

	/// compute the prices for calibMethod
    void ComputeIROptionPrices(ARM_CalibMethod* calibMethod,
								mdmKeysAlias oswModelIdx, 
								bool isFreezeWeights, 
								bool isInitVolParam);
    void ComputeFxOptionPrices(ARM_CalibMethod* calibMethod, 
								bool isFreezeWeights, 
								bool isInitVolParam);

    /// Calibration errors checking
    void CheckCalibErrors();
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

