/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file maturitycapcaluclator.h
 *
 *  \brief
 *
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date September 2004
 */



#ifndef _INGPCALCULATORS_MATURITYCAPCALCULATOR_H
#define _INGPCALCULATORS_MATURITYCAPCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gencalculator.h"
#include "gpbase/countedptr.h"

class ARM_Portfolio;

CC_BEGIN_NAMESPACE( ARM )

class ARM_CurveModelParam;
class ARM_CorrelMatParam;

class ARM_MaturityCapCalculator : public ARM_GenCalculator
{
	static const string MaturityCapColNamesTable [];

public:
    enum MaturityCapColAlias
    {
        EventDate=0,
        StartDate,
        EndDate,
		FirstCouponDate,
		UnderlyingEndDate,
        PaymentDate,
		IndexStartDate,
        IT,
		Coeff,
		DF,
		Annuity,
		Yield0,
		TRINominal,
		INFINENominal,
		INFINENominalForCalib,
		RefStdCapNominal,
		Amortizing,
		Index,
		Spread,
		Yield,
		INFINEYield,
		TRIInterest,
		INFINEInterest,
		Time,
		Volatility,
		Shock,
		EstimTRIYield,
		EstimTRINominal,
		EstimINFINEYield,
		EstimINFINENominal,
		TRIMaturityCap,
		INFINEMaturityCap,
		INFINENominalFinal,
		RefStdCap,
    };

	enum mdmKeysAlias
    {
        YcKey=0,
        CfModelKey,
		MrsKey,
		BetaKey,
		CorrelKey,

        NbKeys
    };

	enum productToPriceAlias
	{
		MaturityCapPrice,
		RefStdCapPrice,
		EstimatedTRI,
		EstimatedNominal,

		NbProductsToPrice
	};

	enum ProductMode
	{
		INFINE = K_INFINE,
		TRI = K_TRI
	};

	enum CalibrationMode
	{
		EX_BOUNDARY,
		ATM,
		FLAT
	};

private:
    
public:

	/// constructor, copy constructor, assignment constructor, destructor
	ARM_MaturityCapCalculator(
		const ARM_Date& startDate,
        const ARM_Date& endDate,
		const ARM_Date& underlyingEndDate,
        int longShort,
        int capFloor,
        int resetFrequency,
		int paymentFrequency,
        const string& indexTerm,
		int	dayCount,
		int intRule,
		double spread,
		double initNominal,
		double initTRI,
		double annuity,
		ProductMode productMode, 
        double coeff,
		const ARM_Curve& amortizing,
        int resetGap,
		const string& restCal,
		const string& payCal,
		CalibrationMode calibrationMode,
		int nbIterations,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
        const ARM_MarketData_ManagerRep& mktDataManager,
        const ARM_StringVector& mdmKeys);

	// SUMMIT constructor
	ARM_MaturityCapCalculator(
		const ARM_Date& asOfDateDate,
		const ARM_Date& startDate,
        const ARM_Date& endDate,
		const ARM_Date& underlyingEndDate,
        int longShort,
        int capFloor,
        int resetFrequency,
		int paymentFrequency,
        const string& indexTerm,
		int	dayCount,
		int intRule,
		double spread,
		double initNominal,
		double initTRI,
		double annuity,
		ProductMode productMode, 
        double coeff,
		const ARM_Curve& amortizing,
        int resetGap,
		const string& restCal,
		const string& payCal,
		CalibrationMode calibrationMode,
		int nbIterations,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
		const ARM_Currency& cpnCcy
		);

	ARM_MaturityCapCalculator(const ARM_MaturityCapCalculator&);
	ARM_MaturityCapCalculator& operator=(const ARM_MaturityCapCalculator&);
	~ARM_MaturityCapCalculator();

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
	virtual void View(char* id = NULL, FILE* ficOut = NULL) const;

	int	 GetNbIterations() const { return itsNbIterations; }
	void SetNbIterations ( int nbIterations ) { itsNbIterations = nbIterations; }

    void SetMaturityCapToPrice(bool toPrice);
    bool IsMaturityCapToPrice() const;

	void SetRefStdCapToPrice(bool toPrice);
    bool IsRefStdCapToPrice() const;

	void SetEstimatedTRIToPrice(bool toPrice);
    bool IsEstimatedTRIToPrice() const;

	void SetEstimatedNominalToPrice(bool toPrice);
    bool IsEstimatedNominalToPrice() const;

	/// Get & Set MRS param for SUMMIT interface (deal dependent data)
    const ARM_CurveModelParam& GetMRS() const;
    void SetMRS(ARM_CurveModelParam* mrsParam);

	/// Get & Set Beta param for SUMMIT interface (deal dependent data)
    const ARM_CurveModelParam& GetBeta() const;
    void SetBeta(ARM_CurveModelParam* mrsParam);

	/// Get & Set Correlation param for SUMMIT interface (deal dependent data)
    const ARM_CorrelMatParam& GetCorrel() const;
    void SetCorrel(ARM_CorrelMatParam* correlationParam);

	ARM_CalibMethod* GetCFCalibMethod();
    const ARM_StdPortfolioPtr GetCFPortfolio(void);
    void SetCFPortfolio(const ARM_StdPortfolio& port);

	void InitMaturityCapForSummit(ARM_ZeroCurve* zc, 
                                  ARM_VolCurve* capVC, 
                                  ARM_VolLInterpol* capRo, 
                                  ARM_VolLInterpol* capNu, 
                                  ARM_VolLInterpol* capBeta = NULL,
                                  int SABRSigmaOrAlpha = 1);

	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;
	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	ARM_INDEX_TYPE GetIndexType();

	virtual void CreateAndSetModel();
	virtual void UpdateModel();
	virtual void CreateAndSetCalibration();
	virtual void UpdateCalibration(bool isUpdateStrike=true);
    virtual double Price();
	virtual void Calibrate();

	/// Specialised version for datas consistency
	virtual void ComputePricingData() const;
    virtual void CheckData();
	virtual void CheckMktData();
    virtual ARM_GenSecurity* GetSubGenSecurity(int indexType) const {return NULL;}

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

	/// Maturity Cap financial datas
	ARM_Date                    itsStartDate;           // start date
	ARM_Date                    itsEndDate;             // end date
	ARM_Date					itsUnderlyingEndDate;	// underlying end date

    int                         itsLongShort;			// long or short position
	int							itsCapFloor;			// cap or floor

    int							itsResetFreq;			// reset frequency
	int							itsPayFreq;				// payment frequency	

	string                      itsIndexTerm;			// index term (only MM index)
    int                         itsDayCount;			// day count
    int                         itsIntRule;             // interest rule

	double						itsSpread;
    
    double                      itsInitNominal;			// initial nominal
    double                      itsInitTRI;				// initial TRI
	double						itsAnnuity;				// annuity

	ProductMode					itsProductMode;			// product mode
	double						itsCoeff;				// coefficient

    int							itsResetGap;			// reset gap
	string                      itsResetCal;			// reset calendar
    string                      itsPayCal;				// payment calendar

	ARM_Curve					itsAmortizing;			// amortizing curve


	// MC Nb iterations
	int itsNbIterations;

	/// Flag to specify product to price
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ itsProductsToPrice;
	bool itsHasBeenPriced;

    /// Flags to enable auto-calibrations
	CalibrationMode itsCalibrationMode;

	/// Model to use
    mdmKeysAlias itsModelKey;

	ARM_StdPortfolioPtrVector itsCapFloorPortfolios;

	ARM_CalibMethodPtrVector itsVolatilityCalibMethods;

	// Pricing Data
	double itsMaturityCapPrice;
	double itsMaturityCapStdDev;
	double itsMaturityCapPricingDuration;
	double itsRefStdCapPrice;
	double itsRefStdCapStdDev;
	ARM_VectorPtr itsEstimatedTRI;
	ARM_VectorPtr itsEstimatedNominal;
	double itsBoundaryShock;

	/// Set the keys or name alias of model (usefull in multi-currency context)
	void SetModelKeys();

	ARM_StdPortfolioPtr CreateCFPortfolio();
	void ComputeCFPrices(
		bool isFreezeWeights, 
		bool isInitParam);
	void ComputeCFStrikesATMOrFLAT();
	void ComputeCFStrikesEX_BOUNDARY();

	/// calibration part
	void CreateEmptyCalibration(bool digitalCalibration,bool swaptionCalibration, bool computeExerStrikes);
};

CC_END_NAMESPACE()

#endif