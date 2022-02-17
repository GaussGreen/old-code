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
#include <inst/stripoption.h>

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
		IrConvAdjModelKey,
        CorrelMatrixKey,
        NbKeys
    };
            
	/// constructor for deal from term sheet
	ARM_FXVanillaCalculator(const ARM_Date& asOfDate,
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
		ARM_VanillaType vanillaType);

	/// constructor from FxSpreadStripOption
	ARM_FXVanillaCalculator(const ARM_FxSpreadStripOption& FxStrip,
							ARM_BasketType minMax,
							ARM_DigitType digitType,
							ARM_VanillaType vanillaType);

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
	inline virtual ARM_GP_Vector GetvCpnNominal() const  { return ARM_GP_Vector(0);};
	inline virtual ARM_GP_Vector GetvFundNominal() const { return ARM_GP_Vector(0);};
	inline virtual ARM_GP_Vector GetvFundSpread() const { return ARM_GP_Vector(0);};

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

private:
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
	int					itsPaymentGap;
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

