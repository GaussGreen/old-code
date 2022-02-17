/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file fxvanillacalculator.h
 *
 *  \brief file for the fxvanilla calculator 
 *	\author  K.Belkheir & E.Ezzine
 *	\version 1.0
 *	\date January 2007
 */


#ifndef _INGPCALCULATORS_FXVanillaCalculator_H
#define _INGPCALCULATORS_FXVanillaCalculator_H

#define UNIMPLEMENTED_FXVANILLA_CALCULATOR_FUNCTION  { ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented function in FXVanillaCalculator"); }

//gpinfra
#include "gpinfra/modelparams.h"

//gpmodels
#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"

//gpcalculators
#include "gencalculator.h"
#include "forexvanilla.h"
#include "fxvanillafactory.h"

// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

// kernel
//#include <inst/stripoption.h>

#include "typedef.h"
CC_USING_NS(std,pair)


/// Kernel forward declaration
//class ARM_StdPortfolio;
//class ARM_Option;

CC_BEGIN_NAMESPACE( ARM )

/// GP forward declaration
//class ARM_DateStrip;

///-----------------------------------------------------------------------------
/// \class ARM_FXVanillaCalculator
/// \brief
///  Class that implements a FX Vanilla calculator
///-----------------------------------------------------------------------------

class ARM_FXVanillaCalculator : public ARM_GenCalculator 
{
public:

    enum mdmKeysAlias
    {
		YcBasisPayKey=0,
        Fx1ModelKey,
        Fx2ModelKey,
		IrModelKey,			// for convexity adjustment
        CorrelMatrixKey,	// correlation between Fx's
		IrFx1CorrelKey,		// correlation for convexity adjustment IR Fx1
		IrFx2CorrelKey,		// correlation for convexity adjustment IR Fx2
        NbKeys
    };
            
	/// constructor for deal from term sheet
	ARM_FXVanillaCalculator(const ARM_Date& asOfDate,
		const ARM_DateStrip& dateStrip,
		const ARM_Date& startDate,
		const ARM_Date& endDate,
		const ARM_Currency& payCcy,
		int expiryGap,
		int settlmentGap,
		int frequency,
		int dayCount, 
		const string& resetCal,
		const string& payCal,
		const string fx1Name,
		const string fx2Name,
		const ARM_Curve& leverage,
		const ARM_Curve& nominal,
		const ARM_Curve& strike,
		int callPut = K_CALL,
		ARM_VanillaType vanillaType = ARM_FXVanillaType::vanilla,
		const ARM_Curve& alpha = ARM_FlatCurve(1.0),
		const ARM_Curve& beta = ARM_FlatCurve(1.0),
		ARM_DigitType digitType = ARM_FXDigitType::centred,
		double epsilon = 1.0e-6,
		const ARM_Curve& strike2 = ARM_Curve(),
		int callPut2 = K_CALL,
		ARM_BasketType minMax = ARM_FXBasketType::max,
		const ARM_Curve& couponMin = ARM_FlatCurve(ARM_NumericConstants::ARM_INFINITY),
		const ARM_Curve& couponMax = ARM_FlatCurve(0.) );

	ARM_FXVanillaCalculator(const ARM_FxSpreadStripOption& FxStrip,
		ARM_BasketType minMax =ARM_FXBasketType::max,
		ARM_DigitType digitType = ARM_FXDigitType::centred,
		ARM_VanillaType vanillaType = ARM_FXVanillaType::vanilla);

	///copy constructor, assignment constructor, destructor
	ARM_FXVanillaCalculator( const ARM_FXVanillaCalculator& rhs );
	ASSIGN_OPERATOR(ARM_FXVanillaCalculator)
	~ARM_FXVanillaCalculator()
	{ 
		delete itsIntermediatePrices; 
		itsIntermediatePrices = NULL;
	};
	
	void Init(	const ARM_MarketData_ManagerRep& mktDataManager,
				const ARM_Curve& Kmin1,
				const ARM_Curve& Kmax1,
				const ARM_Curve& Kmin2,
				const ARM_Curve& Kmax2,
				int Nleft1,
				int Ncenter1,
				int Nright1,
				int Nleft2,
				int Ncenter2,
				int Nright2
				);


	/// ----------------------------------------------------------
	/// ------------- pure virtual functions
	
	/// pure virtual methods to force implementation
	virtual ARM_RowInfo ColumnNames() const {return ARM_RowInfo();};
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const {return ARM_RowInfo();};
	virtual ARM_DateStripCombiner DatesStructure() const ;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

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
	ARM_Vector* ComputeAll();

    /// FXVanilla deal description creation functions
    virtual string toString(const string& indent="",const string& nextIndent="") const;

    /// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_FXVanillaCalculator(*this); }

	virtual string ExportShortName() const { return "LCFXV";}

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetDateStrip() const  { return ARM_DateStripPtr(NULL);};


	///-----------------------------------------------------------
	///---------------------Accessors

protected:
	/// FXVanilla financial datas
	ARM_Date			itsStartDate;           // start date
	ARM_Date			itsEndDate;             // end date
	int					itsDayCount;			// day count
    int					itsFrequency;			// reset & pay frequency 
	
	string				itsFX1name;				//FX1 name
	string				itsFX2name;				//FX2 name
	ARM_Currency		itsPayCcy;				//payment Currency					

	int					itsExpiryGap;
	int					itsSetlmentGap;												
	string				itsResetCalendar;         // reset calendar
    string				itsPayCalendar;           // payment calendar

	int					itsCallPut;
	int					itsCallPut2;

	ARM_Curve			itsNominalCv;
	ARM_Curve			itsStrikeCv;
	ARM_Curve			itsAlphaCv;
	ARM_Curve			itsBetaCv;
	ARM_Curve			itsStrike2Cv;
	ARM_Curve			itsLeverageCv;
	ARM_Curve			itsCouponMinCv;
	ARM_Curve			itsCouponMaxCv;

	ARM_Curve			itsKmin1Cv;
	ARM_Curve			itsKmax1Cv;
	ARM_Curve			itsKmin2Cv;
	ARM_Curve			itsKmax2Cv;

	ARM_VanillaType		itsVanillaType;
	ARM_BasketType		itsMinMax;
	ARM_DigitType		itsDigitType;

	double				itsEpsilon;

	// Nb of Points for GL
	int itsNleft1;
	int itsNcenter1;
	int itsNright1;
	int itsNleft2;
	int itsNcenter2;
	int itsNright2;

	//Intermediate Prices
	ARM_Vector* itsIntermediatePrices;

	/// Dates Strip
	ARM_DateStripPtr itsDateStrip;

	/// to avoid de re-calculate the option price
	bool itsHasBeenComputed;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

