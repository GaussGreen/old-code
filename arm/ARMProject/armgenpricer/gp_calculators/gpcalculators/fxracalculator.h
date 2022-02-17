/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file FXRACalculator.h
 *
 *  \brief file for the fxvanilla calculator 
 *	\author  K.Belkheir & E.Ezzine
 *	\version 1.0
 *	\date January 2007
 */


#ifndef _INGPCALCULATORS_FXRACalculator_H
#define _INGPCALCULATORS_FXRACalculator_H

#define UNIMPLEMENTED_FXRA_CALCULATOR_FUNCTION  { ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented function in FXRACalculator"); }

//gpcalculators
#include "fxvanillacalculator.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_HybridIRFX;
///-----------------------------------------------------------------------------
/// \class ARM_FXRACalculator
/// \brief
///  Class that implements a FX Vanilla calculator
///-----------------------------------------------------------------------------

class ARM_FXRACalculator : public ARM_FXVanillaCalculator
{
public:
            
	/// constructor for deal from term sheet
	ARM_FXRACalculator(const ARM_Date& asOfDate,
		const ARM_DateStrip& dateStrip,
		const ARM_Date& startDate,
		const ARM_Date& endDate,
		int expiryGap,
		int setlmentGap,
		int paymentGap,
		int frequency,
		int dayCount, 
		string resetCal,
		string payCal,
		const string fx1Name,
		const string fx2Name,
		const ARM_Currency& payCcy,
		const ARM_Curve& nominal,
		int callPut,
		const ARM_Curve& strike,
		const ARM_Curve& alpha,
		const ARM_Curve& beta,
		const ARM_Curve& strike2,
		const ARM_Curve& leverage,
		int callPut2,
		ARM_BasketType minMax,
		ARM_DigitType digitType,
		double epsilon,
		ARM_VanillaType vanillaType,
		int	intRule,
		int stubType,
		int resetTiming,
		int fixingfrequency,
		string payIdx,
		double payIdxSpread,
		string payIdxIT,
		string irIdx,
		string irIdxIT,
		const ARM_Curve& fxDownBarrierCv,
		const ARM_Curve& fxUpBarrierCv,
		const ARM_Curve& irDownBarrierCv,
		const ARM_Curve& irUpBarrierCv);

	///copy constructor, assignment constructor, destructor
	ARM_FXRACalculator( const ARM_FXRACalculator& rhs );
	ASSIGN_OPERATOR(ARM_FXRACalculator)
	~ARM_FXRACalculator()
	{
		delete itsAllFixingPrices; 
		itsAllFixingPrices = NULL;
	};
		
	void Init(	const ARM_MarketData_ManagerRep& mktDataManager,
				int Ncenter1,
				int Ncenter2);


	/// ----------------------------------------------------------
	/// ------------- pure virtual functions
	
	/// pure virtual methods to force implementation
	virtual ARM_RowInfo ColumnNames() const {return ARM_RowInfo();};
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const {return ARM_RowInfo();};
	//virtual ARM_DateStripCombiner DatesStructure() const ;
	//virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	/// pricing function
	virtual void ComputePricingData() const;
	virtual void CreateAndSetCalibration() UNIMPLEMENTED_FXVANILLA_CALCULATOR_FUNCTION;
    virtual void UpdateModel()  UNIMPLEMENTED_FXVANILLA_CALCULATOR_FUNCTION;
    virtual void UpdateCalibration(bool isUpdateStrike=true) UNIMPLEMENTED_FXVANILLA_CALCULATOR_FUNCTION;
	virtual void Calibrate() UNIMPLEMENTED_FXVANILLA_CALCULATOR_FUNCTION;
	virtual void CheckData()UNIMPLEMENTED_FXVANILLA_CALCULATOR_FUNCTION;

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

	 /// Initialisation to 0 of all columns of the deal description that could be priced
    virtual void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const UNIMPLEMENTED_FXVANILLA_CALCULATOR_FUNCTION;

	// Cst Manager
	ARM_CstManagerPtr ComputeCstManager() const {return ARM_CstManagerPtr();};

	///-----------------------------------------------------------
	///---------------not virtual functions
	/// Core of the calculator
	virtual void CreateAndSetModel();
    virtual double Price();
    virtual void CheckMktData();
	void FixingsStructure(ARM_Date startDate,
						ARM_Date endDate,
						char* ccyName,
						ARM_IRIndex irIndex,
						string irIndexTerm,
						ARM_PricingModel* model);
	ARM_Vector* ComputeAll();

    /// FXVanilla deal description creation functions
    virtual string toString(const string& indent="",const string& nextIndent="") const;

    /// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_FXRACalculator(*this); }

	virtual string ExportShortName() const { return "LCFXV";}

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetDateStrip() const  { return itsDateStrip;};

	inline const ARM_Vector* GetvAllFixingPrices() const  { return (ARM_Vector*)itsAllFixingPrices->Clone();};
	inline virtual ARM_Curve GetFXDownBarrier() const  { return itsFxDownBarrierCv;};
	inline virtual ARM_Curve GetFXUpBarrier() const  { return itsFxUpBarrierCv;};
	inline virtual ARM_Curve GetIRDownBarrier() const  { return itsIrDownBarrierCv;};
	inline virtual ARM_Curve GetIRUpBarrier() const  { return itsIrUpBarrierCv;};

	inline virtual std::vector<double> GetIRIndexStartTimes() const  { return itsIRindexStartTimes;};
	inline virtual std::vector<double> GetIRIndexEndTimes() const  { return itsIRindexEndTimes;};
	inline virtual std::vector<double> GetIRIndexTerms() const  { return itsIRindexTerms;};
	inline virtual string GetIRIdxIT() const  { return itsIrIdxIT;};
	inline virtual char* GetPayCcyName() const  { return itsPayCcy.GetCcyName();};


private:
	int	itsIntRule;
	int itsStubType;
	int itsResetTiming;
	int itsFixingfrequency;
	string itsPayIdx;
	double itsPayIdxSpread;
	string itsPayIdxIT;
	string itsIrIdx;
	string itsIrIdxIT;
	ARM_Curve itsFxDownBarrierCv;
	ARM_Curve itsFxUpBarrierCv;
	ARM_Curve itsIrDownBarrierCv;
	ARM_Curve itsIrUpBarrierCv;
	std::vector<double> itsFixingTimes;
	std::vector<double> itsIRindexResetTimes;
	std::vector<double> itsIRindexStartTimes;
	std::vector<double> itsIRindexEndTimes;
	std::vector<double> itsIRindexTerms;

	//Intermediate of Intermediate Prices: price of all fixings (used to calibrate local model for CRA double)
	ARM_Vector* itsAllFixingPrices;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

