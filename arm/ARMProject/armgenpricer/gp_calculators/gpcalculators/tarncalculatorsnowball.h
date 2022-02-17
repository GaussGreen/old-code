/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tarncalculatornowball.h
 *
 *  \brief
 *
 *	\author  JP Riaudel
 *	\version 1.0
 *	\date Janvier 2005
 */

#ifndef _INGPCALCULATORS_TARNCALCULATOR_SNOWBALL_H
#define _INGPCALCULATORS_TARNCALCULATOR_SNOWBALL_H

#include "tarncalculator.h"


CC_BEGIN_NAMESPACE( ARM )


class ARM_TARNCalculatorSnowBall : public ARM_TARNCalculator
{
public:

	/// constructor, copy constructor, assignment constructor, destructor
	ARM_TARNCalculatorSnowBall(const ARM_Date& startDate,
    const ARM_Date& endDate,
    const ARM_Curve& strike,
	double coupon0,
    int payRec,
    int cpnDayCount,
    int cpnFreq,
    int cpnTiming,
    const string& cpnIndexTerm,
    int cpnIndexDayCount,
    const string& cpnResetCal,
    const string& cpnPayCal,
    int cpnResetGap,
	int intRule,
	int stubRule,
	long reverse,
    const ARM_Curve& leverage,
	const ARM_Curve& couponMin,
	const ARM_Curve& couponMax,
	const ARM_Curve& levPrev,
	double lifeTimeCapTarget,
	bool globalCapFlag,
    double lifeTimeFloorTarget,
    const ARM_Curve& fundSpread,
    int fundFreq,
    int fundDayCount,
    const ARM_Curve& nominal,
	const ARM_Curve& fees,
	const ARM_Curve& fundNominal,
	const std::vector<double>& nbIterations,
	ARM_ModelType modelType,
	CalibrationMode capFloorCalibMode,
	CalibrationMode digitalCalibMode,
	bool oswCalibFlag,
	bool controlVariableFlag,
	bool digitalSmoothingFlag,
	bool smiledFRMRescallingFlag,
	const string& genType1,
	const string& genType2,
	int firstNbTimes,
	int firstNbDims,
	const string& pathScheme,
	const string& pathOrder,
	bool antithetic,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
	const ARM_Currency& cpnCcy,
    const ARM_Currency& fundingCcy,
    const ARM_Currency& basisCcy,
    const ARM_Currency& domesticCcy,
    const ARM_Currency& foreignCcy,
    const ARM_MarketData_ManagerRep& mktDataManager,
    const ARM_StringVector& mdmKeys);

	ARM_TARNCalculatorSnowBall(const ARM_TARNCalculatorSnowBall&);
	ARM_TARNCalculatorSnowBall& operator=(const ARM_TARNCalculatorSnowBall&);
	~ARM_TARNCalculatorSnowBall(){};

	virtual string CreateCpnDescription(bool isFirstExer) const;
	virtual string CreateCFCpnDdescription( bool isFirstExer, const string& cpnModelName) const;

	virtual string CreateExStrikeDescription(bool isFirstExer, const string& cpnModelName) const;
	virtual void SetCVSwapMiddleText( ARM_DealDescriptionPtr& dealDesc, size_t knockOutIndex ) const;

	double GetC0() const {return itsCoupon0;};
	double GetCoefReverse() const {return itsReverse? -1.0: 1.0;};
	virtual void UpdateCalibration(bool isUpdateStrike=true);

private:

    bool itsReverse; // true reverse, false not reverse

	double itsCoupon0;
	bool itsNeedOtherPayoff;
};

CC_END_NAMESPACE()

#endif