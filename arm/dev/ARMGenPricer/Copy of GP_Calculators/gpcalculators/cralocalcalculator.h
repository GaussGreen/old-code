/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file CRALocalCalculator.h
 *
 *  \brief
 *
 *	\author  H. BAKHTRI & M. ABDELMOUMNI & P. LAM
 *	\version 1.0
 *	\date August 2005
 */

#ifndef _INGPCALCULATORS_CRALOCALCALCULATOR_H
#define _INGPCALCULATORS_CRALOCALCALCULATOR_H

#include "cracalculator.h"

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gencalculator.h"
#include "gpbase/countedptr.h"
#include "gpbase/gpmatrix.h"

/// kernel
#include "util/refvalue.h"
#include "util/exercise.h"
#include "ccy/currency.h"
#include "inst/optionportfolio.h"
#include "mod/bssmiled.h"
#include "mod/armfrmmodel.h"

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

class ARM_CRALocalCalculator : public ARM_CRACalculator
{
protected:

	static const string LocalFloorCRAColNamesTable [];
	static const string LocalCapFloorCRAColNamesTable [];
	static const int LocalFloorCRAProductToPriceColumns [];
	static const int LocalCapFloorCRAProductToPriceColumns [];

	static const size_t LocalCapFloorCRAExerciseProbaNb;
	static const string LocalCapFloorCRALocalExerciseProbaNamesTable [];
	static const string LocalCapFloorCRAExerciseProbaNamesTable [];

public:

	//LocalCapFloor
	enum LocalCapFloorCRAColAlias
	{
		localCfResetDate = 0,
		localCfStartDate,
		localCfMaturityDate,
		localCfFees,
		localCfFundingLeg,
		localCfCorridorLeg1,
		localCfCorridorLeg2,
		localCfCorridorLeg3,
		localCfCorridorLeg,
		localCfExoticSwap,
		localCfOption, 
		localCfBermuda,
		localCfCorridor,
		localCfFunding,
		localCfFwdCorridor,
		localCfFwdFunding,
		localCfExerSwapRate,
		localCfFrontier,
		localCfNextResetDate,
		localCfIndicExer,
		localCfLocalProba1,
		localCfProba1,
		localCfLocalProba2,
		localCfProba2,
		localCfLocalProba3,
		localCfProba3,
		localCfLocalProba4,
		localCfProba4,
		localCfLocalProba5,
		localCfProba5,
		localCfLocalProba6,
		localCfProba6,
		localCfLocalProba7,
		localCfProba7,
		localCfLocalProba8,
		localCfProba8,
		localCfLocalProba9,
		localCfProba9,
		localCfLocalProba10,
		localCfProba10,
		localCfLocalProba11,
		localCfProba11,
		localCfLocalProba12,
		localCfProba12,
		localCfLocalProba13,
		localCfProba13,
		localCfLocalProba14,
		localCfProba14,
		localCfLocalProba15,
		localCfProba15,
		localCfLocalProba16,
		localCfProba16,
		localCfLocalProba17,
		localCfProba17,
		localCfLocalProba18,
		localCfProba18,
		localCfLocalProba19,
		localCfProba19,
		localCfLocalProba20,
		localCfProba20
	};
	
	enum LocalCapFloorCRAproductToPriceAlias
	{
		localCfBermudaPrice,
		localCfCorridorPrice,
		localCfFundingPrice,
		localCfFwdCorridorPrice,
		localCfFwdFundingPrice,
		localCfCorridorLegPrices,
			
		localCfNbProductsToPrice
	};

	//LocalFloor
	enum LocalFloorCRAColAlias
	{
		localFResetDate = 0,
		localFStartDate,
		localFMaturityDate,
		localFFees,
		localFFundingLeg,
		localFCorridorDown,
		localFCorridorUp,
		localFCorridorLeg,
		localFExoticSwap,
		localFOption, 	
		localFBermuda,
		localFCorridor,
		localFFunding,
		localFFwdCorridorUp,
		localFFwdCorridorDown,
		localFFwdCorridor,
		localFFwdFunding,
	};

	enum LocalFloorCRAproductToPriceAlias
	{
		localFBermudaPrice,
		localFCorridorPrice,
		localFFundingPrice,
		localFFwdCorridorPrice,
		localFFwdFundingPrice,
		
		localFNbProductsToPrice
	};

	//Local Type
	enum LocalModelType
	{
		LocalDown	= 0,
		LocalUp		= 1,
		LocalDownUp	= 2
	};

	enum CapFloor
    {
        In = K_IN,
		Out = K_OUT,
    };

	//Constructor 
	ARM_CRALocalCalculator(const ARM_MarketData_ManagerRep& mktDataManager );

	ARM_CRALocalCalculator(	ARM_Currency& ccy,
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
							int cpnPayGap,
							int refIndexType,
							string refTerm,
							int refDayCount,
							double refCoeff,
							int localModelType,
							double meanRevMin,
							double meanRevMax,
							double betaMin,
							double betaMax,
							ARM_Vector* calibSecPFParams,
							int nbSteps,
							int flagToGenerateOSWATM,
							int localResetFreq,
							ARM_StringVector& mdmKeys,
							const ARM_MarketData_ManagerRep& mktDataManager,
							const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
							int reCalibMrs = -1,
							int reCalibBeta = -1,
							bool isStdCalib = true,
							bool isExerciseProbas=false,
							size_t exerciseProbaOffset=0);

	//CRA Spread Contructor
	ARM_CRALocalCalculator( const ARM_Currency& ccy,
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
							int cpnPayGap,
							int refIndexType1,
							const string& refTerm1,
							int refDayCount1,
							double refCoeff1,
							bool isPortfolioNoticeDays,
							int localModelType,
							const ARM_StringVector& mdmKeys,
							const ARM_MarketData_ManagerRep& mktDataManager,
							const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
							int reCalibMrs = -1,
							int reCalibBeta = -1,
							bool isStdCalib = true,
							bool isExerciseProbas=false,
							size_t exerciseProbaOffset=0);

	//Summit  CRA Constructor 
	ARM_CRALocalCalculator( ARM_OptionPortfolio* optPortfolio,
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
							int reCalibMrs,
							int reCalibBeta,
							int localModelType,
							int localResetFreq,
							bool isStdCalib = true);

	//Summit CRA Spread constructor
	ARM_CRALocalCalculator( ARM_OptionPortfolio* optPortfolio,
							ARM_StringVector& mdmKeys,
							const ARM_MarketData_ManagerRep& mktDataManager,
							const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice);

	ARM_CRALocalCalculator( const ARM_Date& asOfDate,
							ARM_Security* security);

	ARM_CRALocalCalculator(const ARM_CRALocalCalculator&);
	
	ARM_CRALocalCalculator& operator = (const ARM_CRALocalCalculator&);
	
	~ARM_CRALocalCalculator();

	virtual ARM_Object*		Clone() const;
	
	//Accessors
	int GetLocalModelType() const {return itsLocalModelType;}
	void SetLocalModelType(int localModelType){itsLocalModelType = localModelType;}

	size_t GetDescSize() const;

	//Deal Description
	virtual ARM_StringVector		PricedColumnNames() const;
	virtual ARM_RowInfo				ColumnNames() const;
	virtual ARM_RowInfo				MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	void							InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;
	virtual ARM_DateStripCombiner	DatesStructure() const;
	virtual ARM_GP_Vector*			ComputeExerciseDatesAndSetFees() const;

    //Overloaded  methods
	virtual void					UpdateCalibration(bool isUpdateStrike=true);
	virtual void					CreateAndSetModel();
	virtual void					CreateAndSetLocalModel();
	virtual void					CalibrateLocalModel();

	virtual ARM_PricingModel*		GetSFRMModel(void) const;

	void							CreateCorridorLocalPortfolio();

	ARM_StdPortfolioPtr				GetLocalPortfolio(){return itsLocalPortfolio;};
	void							SetLocalPortfolio(ARM_StdPortfolioPtr localPortfolio){itsLocalPortfolio = localPortfolio;};
	int								GetLocalResetFreq(){return itsLocalResetFreq;};

    virtual void					UpdateModel();
	virtual void					Calibrate();
	virtual double					Price();
	virtual void					ComputePricingData() const;

	//Some utilities
	virtual string					ExportShortName() const { return "LCCRA";}
	string							toString(const string& indent, const string& nextIndent) const;
	virtual string					GeneralDataToString(const string& indent = "", const string& nextIndent = "") const;
	virtual void					View(char* id = NULL, FILE* ficOut = NULL) const;

	//Utilities
    void CleanUp();
    void CopyNoCleanUp(const ARM_CRALocalCalculator& rhs);

protected:

	//Flags dedicated to local CRA Calculator
	int								itsLocalModelType;	    // LocalDown  or LocalUp	  or LocalDownUp
	int								itsLocalResetFreq; 
	ARM_StdPortfolioPtr				itsLocalPortfolio;


	/// To manage exercise probability computation
	bool			itsIsExerciseProbas;
	size_t			itsExerciseProbaOffset;
	ARM_GP_Vector	itsExerciseProbas;

};

CC_END_NAMESPACE()

#endif
