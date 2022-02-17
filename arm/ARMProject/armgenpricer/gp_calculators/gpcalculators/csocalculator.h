/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *      \file csocalculator.h
 *
 *  \brief
 *
 *	\author  JP Riaudel
 *	\version 1.0
 *	\date May 2005
 */

#ifndef _INGPCALCULATORS_CSOCALCULATOR_H
#define _INGPCALCULATORS_CSOCALCULATOR_H

#include "gpbase/port.h"
#include "gencsocalculator.h"

#include "gpbase/datestrip.h"

class ARM_SpreadOption;

CC_BEGIN_NAMESPACE( ARM )
struct ARM_VanillaSpreadOptionArg;

class ARM_CSOCalculator : public ARM_GenCSOCalculator
{

protected:
	static const string CSOColNamesTable [];
	static const int ProductToPriceColums[];

public:
	enum productToPriceAlias
	{
		CSOPrice,
		structPrice,
		NbProductsToPrice
	};


    enum CSOColAlias
    {
    EventDate=0,
    StartDate,
    PayDate,
    EndDateCMS1,
    EndDateCMS2,
    IT,
    CMS1,
    CMS2,
    CpnFix,
	Strike,
	Notional,
	Leverage,
    SpreadOptionLet,
	SpreadOption,
	StrikeFlow,
	PVStrikeFlow,
    FundStartDate,
    FundEndDate,
	FundSpread,
	FundFragmentNPV,
	FundNPV,
    Flow,
    TotalFlow,
    Fees,
	BermudaProfile,
	CSO/*,
	ExerciseCondition,
	NumeraireDate,
	ExerciseConditionOfIndex,
	ProbaOfExercise,
	Funding,
	Strike1,
	Strike2,
	CapCMS1,
	CapCMS2*/
    };

	enum CalibrationMode
	{
		NOCalib=0,
		VegaBestFitCalib,
		BestFitCalib,
		BootstrapCalib,
		ONEDCalib,
		TWODCalib
	};



	/// constructor, copy constructor, assignment constructor, destructor
	ARM_CSOCalculator(const ARM_Date& startDate,
					  const ARM_Date& endDate,
					  int CMSLong,
					  int CMSShort,
					  int cpnDayCount,
					  int cpnFreq,
					  int cpnResetTiming,
					  const ARM_Curve& nominal,
					  const ARM_Curve& fixCoupon,
					  const ARM_Curve& leverageLong,
					  const ARM_Curve& leverageShort,
					  const ARM_Curve& cpnMin,
					  const ARM_Curve& cpnMax,
					  const ARM_Curve& strike,
					  int fundFreq,
					  int fundDayCount,
					  const ARM_Curve& fundSpread,
					  int exerciseFreq,
					  int noticeGap,
					  int payRec,
					  const ARM_Curve& fees,
					  CalibrationMode sigmaCalibFlag,
					  CalibrationMode NDFlag,
					  bool mrsCalibFlag,
					  bool thetaCalibFlag,
					  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
					  vector<double> ModelDatas,
					  const ARM_MarketData_ManagerRep& mktDataManager,
					  const ARM_StringVector& mdmKeys);

	ARM_CSOCalculator(const ARM_Date& asOfDate,
					  const ARM_Date& startDate,
					  const ARM_Date& endDate,
					  int CMSLong,
					  int CMSShort,
					  int cpnDayCount,
					  int cpnFreq,
					  int cpnResetTiming,
					  const ARM_Curve& nominal,
					  const ARM_Curve& fixCoupon,
					  const ARM_Curve& leverageLong,
					  const ARM_Curve& leverageShort,
					  const ARM_Curve& cpnMin,
					  const ARM_Curve& cpnMax,
					  const ARM_Curve& strike,
					  int fundFreq,
					  int fundDayCount,
					  const ARM_Curve& fundSpread,
					  int exerciseFreq,
					  int noticeGap,
					  int payRec,
					  const ARM_Curve& fees);

// FIXMEFRED: mig.vc8 (23/05/2007 11:38:52):missing return type
	void InitCSOFromSummit(const ARM_MarketData_ManagerRep* mktDataManager,
					  const ARM_StringVector& mdmKeys,
					  CalibrationMode sigmaCalibFlag,
					  CalibrationMode NDFlag,
					  bool mrsCalibFlag,
					  bool thetaCalibFlag,
					  std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
					  vector<double> ModelDatas);

	ARM_CSOCalculator(const ARM_CSOCalculator& rhs);
	ARM_CSOCalculator& operator=(const ARM_CSOCalculator&);
	virtual ~ARM_CSOCalculator();

	virtual ARM_Object* Clone() const;

	void CreateEmptyCalibration();


	virtual void CreateAndSetCalibration();
	virtual void UpdateModel();
	virtual void UpdateCalibration(bool isUpdateStrike=true);
	virtual void Calibrate();
    virtual double Price();

	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure) const;

	virtual void ComputePricingData() const;
	
	const ARM_StdPortfolioPtr ARM_CSOCalculator::GetOSWPortfolio() const;
	const ARM_StdPortfolioPtr ARM_CSOCalculator::GetSOPortfolio() const;
	const ARM_StdPortfolioPtr ARM_CSOCalculator::GetCMSLONGPortfolio() const;
	const ARM_StdPortfolioPtr ARM_CSOCalculator::GetCMSSHORTPortfolio() const;
	ARM_CalibMethod* GetSWOPTCalibMethod() const;

	virtual void CreateAndSetModel();

	/// Specialised version for datas consistency
    virtual void CheckData();
	virtual void CheckMktData();

private:

	CalibrationMode			itsSigmaCalib;
	CalibrationMode			itsCalibND;
	bool		itsThetaCalib;
	bool		itsMrsCalib;
	bool		itsUpdateInitCalib;
	
	ARM_StdPortfolioPtr itsSwaptionCMS1PF;
	ARM_StdPortfolioPtr itsSwaptionCMS2PF;
	ARM_StdPortfolioPtrVector itsSwaptionCMSPF_1_2;
	ARM_StdPortfolioPtr itsSpreadOptionPF;
	ARM_StdPortfolioPtr itsSwaptionPF;
	ARM_StdPortfolioPtr itsSOAndSWOPtPF;
	
	ARM_ModelParamPtr itsTWODInitSigmaParam;
	ARM_ModelParamPtr itsTWODInitThetaParam;

    vector<double> itsModelDatas;
/*		ModelDatas[0] = SIGMA_CSO_ARM_MAX_ITER;
		ModelDatas[1] = SIGMA_DEFAULT_PRECISION;
		ModelDatas[2] = CORREL_CSO_ARM_MAX_ITER;
		ModelDatas[3] = CORREL_DEFAULT_PRECISION;
		ModelDatas[4] = CORREL_DEFAULT_PONDERATION;
*/

	virtual void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	ARM_CalibMethod* GetSigmaCalibMethod() const;

	void ComputeEquivalentSwoptStrikes(std::vector<double>& equivStrikes1, std::vector<double>& equivStrikes2);
	void GenerateEquivalentDiagSwoptStrikes(std::vector<double>& equivStrikes);
	void CreateCSOSwaptionCap();
	void CreateDiagonalSwaption();
	void ComputeDiagonalSwaptionPrice(bool isFreezeWeights, bool isArrearCpn);
	void CreateSOAndSWOPTPF();
	ARM_CalibMethodPtr CreateSigmaCalibration();

	void computeSigmaPortfolioPrices();
	void computeSOAndSWOPTPortfolioPrices();
	void computeSOPortfolioPrices();
	void computeDIAGSwaptionPortfolioPrices();
};

CC_END_NAMESPACE()

#endif //_INGPCALCULATORS_CSOCALCULATOE_H