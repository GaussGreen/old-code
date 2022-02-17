/*!
 *
 * Copyright (c) IXIS CIB June 2005 Paris
 *
 *	\file vanillamepi.h
 *
 *  \brief vanilla spreadoption
 *
 *	\author  A SCHAULY
 *	\version 1.0
 *	\date June 2005
 */

#ifndef _INGPCALIB_VANILLAMEPI_H
#define _INGPCALIB_VANILLAMEPI_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "vanillaargnumeric.h"
#include "gpinfra/nummethod.h"

class ARM_PriceModel;

CC_BEGIN_NAMESPACE( ARM )
///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaMepi
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaMepi:public ARM_VanillaArgNumeric
{
	static const string VanillaMepiColNamesTable [];
	static const string VanillaMepiConstantsTable [];

public:
    enum VanillaMepiColAlias
    {
		EventDate,
		NextDate,
		MaturityDate,
		Underlying,
		CapiFactor,
		PaidCapiFactor,
		Cash,
		Portfolio,
		AvgPortfolio,
		FixedProtection,
		ActifRisque,
		CF
    };

	enum VanillaMepiCstAlias
	{
		RiskFactor,
		Strike,
		MaxBorrow,
		StartingPortfolio,
		StartingCash,
		MinInvested,
		LeverageCost,
		CashSpread,
		Fees,
		AlreadyAsianed,
		AsianDatesNb
	};

private:

public:

	ARM_VanillaMepi(
		const string& curveName,
		const string& equityModelName,
		double startDate,
		double endDate,
		long resetFreq,
		double riskFactor,
		double strike,
		double maxBorrow,
		double protectionCurveStart,
		double protectionCurveEnd,
		double startingPortfolio,
		double startingCash,
		double minInvested,
		double leverageCost,
		double cashSpread,
		double fees, 
		double alreadyAsianed,
		long AsianingPeriodNb);
		
    ARM_VanillaMepi(const ARM_VanillaMepi& arg);
    ARM_VanillaMepi& operator=(const ARM_VanillaMepi& rhs);
    virtual ~ARM_VanillaMepi();

	virtual ARM_NumMethod::GP_PricingDirection GetItsPricingDirection() const {return ARM_NumMethod::GP_FWDLOOKING;}//GP_FWDLOOKING
	
	/// Deal Description
	virtual ARM_RowInfo NumArgColumnNames() const;
	virtual ARM_RowInfo NumArgMiddleRows( size_t eventIdx, const ARM_GP_VectorPtr& eventDates ) const;
	virtual ARM_GP_VectorPtr DatesStructure() const;
	virtual const ARM_CstManagerPtr securityCstManager() const;

	/// Pricing and Hedge
	virtual double Price(ARM_PricingModel* model ) const;
	virtual double AnalyticPrice(ARM_PricingModel* model ) const;
	virtual double ComputeDelta(ARM_PricingModel* model, double bump ) const;

    virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_MEPI; }

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private:
	/// Input data
	string itsEquityModelName;
	double itsStartDate;
	double itsEndDate;
	double itsStrike;
	double itsMaxBorrow;
	double itsProtectionCurveStart;
	double itsProtectionCurveEnd;
	double itsStartingPortfolio;
	double itsStartingCash;
	double itsMinInvested;
	double itsLeverageCost;
	double itsCashSpread;
	double itsFees;
	double itsBorrowingRate;
	double itsRiskFactor;
	double itsInitialProtection;
	double itsFinalProtection;
	double itsAlreadyAsianed;
	long itsAsianingPeriodNb;

	/// Computed data
	double itsMaturity;
	double itsProtectionStep;
	ARM_GP_VectorPtr itsEventDates;


	/// Private functions
	void CopyNoCleanUp(const ARM_VanillaMepi& rhs);
	void CleanUp();
};




////////////////////////////////////////////////////////////////////////////////////////////


///                                       VanillaMepiDelta


////////////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaMepiDelta
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaMepiDelta:public ARM_VanillaArgNumeric
{
	static const string VanillaMepiColNamesTable [];
	static const string VanillaMepiConstantsTable [];

public:
    enum VanillaMepiColAlias
    {
		EventDate,
		NextDate,
		MaturityDate,
		Underlying,
		CapiFactor,
		FixedProtection,
		PaidCapiFactor1,
		Cash1,
		Portfolio1,
		AvgPortfolio1,
		ActifRisque1,
		PaidCapiFactor2,
		Cash2,
		Portfolio2,
		AvgPortfolio2,
		ActifRisque2,
		CF
    };

	enum VanillaMepiCstAlias
	{
		RiskFactor,
		Strike,
		MaxBorrow,
		StartingPortfolio,
		StartingCash,
		MinInvested,
		LeverageCost,
		CashSpread,
		Fees,
		AlreadyAsianed,
		AsianDatesNb,
		DeltaShift
	};

private:

public:

	ARM_VanillaMepiDelta(
		const string& curveName,
		const string& equityModelName,
		double startDate,
		double endDate,
		long resetFreq,
		double riskFactor,
		double strike,
		double maxBorrow,
		double protectionCurveStart,
		double protectionCurveEnd,
		double startingPortfolio,
		double startingCash,
		double minInvested,
		double leverageCost,
		double cashSpread,
		double fees, 
		double alreadyAsianed,
		long AsianingPeriodNb,
		double deltaShift);
		
    ARM_VanillaMepiDelta(const ARM_VanillaMepiDelta& arg);
    ARM_VanillaMepiDelta& operator=(const ARM_VanillaMepiDelta& rhs);
    virtual ~ARM_VanillaMepiDelta();

	virtual ARM_NumMethod::GP_PricingDirection GetItsPricingDirection() const {return ARM_NumMethod::GP_FWDLOOKING;}//GP_FWDLOOKING
	
	/// Deal Description
	virtual ARM_RowInfo NumArgColumnNames() const;
	virtual ARM_RowInfo NumArgMiddleRows( size_t eventIdx, const ARM_GP_VectorPtr& eventDates ) const;
	virtual ARM_GP_VectorPtr DatesStructure() const;
	virtual const ARM_CstManagerPtr securityCstManager() const;

	/// Pricing and Hedge
	virtual double Price(ARM_PricingModel* model ) const;
	virtual double AnalyticPrice(ARM_PricingModel* model ) const;
	virtual double ComputeDelta(ARM_PricingModel* model, double bump ) const;

    virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_MEPI; }

	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private:
	/// Input data
	string itsEquityModelName;
	double itsStartDate;
	double itsEndDate;
	double itsStrike;
	double itsMaxBorrow;
	double itsProtectionCurveStart;
	double itsProtectionCurveEnd;
	double itsStartingPortfolio;
	double itsStartingCash;
	double itsMinInvested;
	double itsLeverageCost;
	double itsCashSpread;
	double itsFees;
	double itsBorrowingRate;
	double itsRiskFactor;
	double itsInitialProtection;
	double itsFinalProtection;
	double itsAlreadyAsianed;
	double itsDeltaShift;
	long itsAsianingPeriodNb;

	/// Computed data
	double itsMaturity;
	double itsProtectionStep;
	ARM_GP_VectorPtr itsEventDates;


	/// Private functions
	void CopyNoCleanUp(const ARM_VanillaMepiDelta& rhs);
	void CleanUp();
};
CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
