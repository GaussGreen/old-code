/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *      \file locallculator.h
 *
 *  \brief
 *
 *	\author  JP Riaudel, Y. Khlif, A. Chaix
 *	\version 1.0
 *	\date July05
 */

#ifndef _INGPCALCULATORS_LOCALCSOCALCULATOR_H
#define _INGPCALCULATORS_LOCALCSOCALCULATOR_H

//////////////////////////////////////////////////////////////////////////
/// class : ARM_LocalCSOCalculator for CSO pricing & risk management
/// o methodology : 
///		+ Hull-White 
///		+ basket calibration 
///		+ local volatility
//////////////////////////////////////////////////////////////////////////

#include "gpbase/port.h"
#include "gencsocalculator.h"

#include "gpbase/datestrip.h"
#include "gpmodels/typedef.h"
#include "gpmodels/marketirmodel.h"
#include "gpinfra/curvemodelparam.h"

#include "crv/volcube.h"
#include "crv/correlmanager.h"
#include "mod/convadjustmanager.h"
#include "inst/forex.h"

using ARM::ARM_CurveModelParam;

class ARM_SpreadOption;
class ARM_Swaption;

CC_BEGIN_NAMESPACE( ARM )

struct ARM_VanillaSpreadOptionArg;

class ARM_LocalCSOCalculator : public ARM_GenCSOCalculator
{

protected:
	static const string LocalCSOColNamesTable [];
	static const int	LocalCSOProductToPriceColums[];
	static const string LocalCSOProbaColNamesTable [] ;
	static const string LocalCSOProbaPricedColNamesTable [];


public:
	/// ----------------------
	/// ---  Enums -----
	/// ----------------------
	enum LocalCSOproductToPriceAlias
	{
		CSOPrice,
		FundingPrice,
		FixPrice,
		FloorPrice,
		CapPrice,
		
		NbProductsToPrice
	};


    enum LocalCSOColAlias
    {
		EventDate = 0,
		StartDate,
		EndDate,
		FixFundLeg,
		FloorFundLeg,	
		CapFundLeg,
		FundingLeg,
		FixLeg,
		FloorLeg,	
		CapLeg,
		SOLeg,
		Fees,
		Option,
		CSO,
		Funding,
		Fix,
		Floor,
		Cap,
		ExerSwapRate,
		Frontier
    };
	
	enum modelsAlias
    {
        myRefModel=0,			 /// the stochastic IR model
		myLocalFloorModel,		 /// the local IR Model to price floor
		myLocalCapModel,		 /// the local Ir Model to price cap	
		myBasisMarginModel,		 /// the forward margin model for basis swap

        NbModels
    };

	enum CalibrationType
	{
		DIAG_CALIBRATION = 0,
		BASKET_CALIBRATION,
		DIAG_BASKET_CALIBRATION,
		DIAG_SPREAD_LONG,
		DIAG_SPREAD_SHORT,
		SHORT_LONG_SPREAD,
		LONG_CALIBRATION,
		DIAG_LONG_CALIBRATION
	};

	enum CalibStrikeType
	{
		ATM = 0,
		EQUIVALENT,
		ZERO,
		CAP,
		FLOOR,
		FRONTIER
	};

	/// ------------------------------------------------------------------------
	/// ---- constructor, copy constructor, destructor
	/// ------------------------------------------------------------------------
	/// contextual
	ARM_LocalCSOCalculator(const ARM_Date& startDate,
						   const ARM_Date& endDate,
						   int CMSLong,
						   int CMSShort,
						   int cpnDayCount,
						   int cpnFreq,
						   int cpnResetTiming,
						   const ARM_Curve& cpnnominal,
						   const ARM_Curve& fixCoupon,
						   const ARM_Curve& leverageLong,
						   const ARM_Curve& leverageShort,
						   const ARM_Curve& cpnMin,
						   const ARM_Curve& cpnMax,
						   const ARM_Curve& strike,
						   int fundFreq,
						   int fundDayCount,
						   const ARM_Curve& fundSpread,
						   const ARM_Curve& fundLeverage,
						   int exerciseFreq,
						   int noticeGap,
						   int payRec,
						   const ARM_Curve& fees,
						   bool switchFlag,
						   int fundingType,
						   std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
						   vector<double> ModelDatas,
						   const ARM_MarketData_ManagerRep& mktDataManager,
						   const ARM_StringVector& mdmKeys,
						   ARM_ModelType		modelType,
						   CalibrationType	calibrationType,
						   vector<CalibStrikeType>	calibStrikeType,
						   ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod = ARM_MarketIRModel::MONEYNESS);

	/// summit construction with Init
	ARM_LocalCSOCalculator(const ARM_Date& asOfDate,
						   const ARM_Date& startDate,
                           const ARM_Date& fixEndDate,
						   const ARM_Date& endDate,
						   const ARM_Currency& CpnCcy,
						   const ARM_Currency& FundCcy,
						   int CMSLong,
						   int CMSShort,
						   int cpnDayCount,
						   int cpnFreq,
						   int cpnResetTiming,
						   const ARM_Curve& cpnnominal,
						   const ARM_Curve& fixCoupon,
						   const ARM_Curve& leverageLong,
						   const ARM_Curve& leverageShort,
						   const ARM_Curve& cpnMin,
						   const ARM_Curve& cpnMax,
						   const ARM_Curve& strike,
						   int fundFreq,
						   int fundDaycount,
						   const ARM_Curve& fundNominal,
						   const ARM_Curve& fundSpread,
						   const ARM_Curve& fundLeverage,
						   int exerciseFreq,
						   int noticeGap,
						   int payRec,
						   const ARM_Curve& fees,
						   bool switchFlag = false,
						   int fundingType = 0,
						   int nbNoCall = 0);

	void InitCSOFromSummit( ARM_ZeroCurve* zc,
							ARM_VolCurve* capVol,
							ARM_VolCurve* swoptVol,
							int SABRSigmaOrAlpha,
							ARM_VolCurve* rhoCap,
							ARM_VolCurve* nuCap,
							ARM_VolCurve* betaCap,
							ARM_VolCurve* rhoSwopt,
							ARM_VolCurve* nuSwopt,
							ARM_VolCurve* betaSwopt,
							ARM_VolCurve* flatVol,
							ARM_VolCurve* convAdjustVol,
							ARM_ConvAdjustManager* convAdjustManager,
							ARM_VolCube* correlDiagCap,
							ARM_VolCube* correlDiagSwopt,
							ARM_CorrelManager* correlCorr,
							ARM_CurveModelParam* mrs,
							ARM_CurveModelParam* correl,
							ARM_CurveModelParam* volRatio,
							ARM_CurveModelParam* mrsSpread,
							std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
							vector<double> modelDatas,
							CalibrationType	calibrationType,
							vector<CalibStrikeType>	calibStrikeType,
							ARM_ModelType modelType,
							ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod,
							ARM_Forex* forex = NULL,
							ARM_ZeroCurve* fundZc = NULL,
							ARM_ZeroCurve* domBasisZc  = NULL,
							ARM_ZeroCurve* fundBasisZc = NULL);

	void InitCSOFromSummit( const ARM_MarketData_ManagerRep* mktDataManager,
							std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
							vector<double> ModelDatas,
							CalibrationType	calibrationType,
							vector<CalibStrikeType>	calibStrikeType,
							ARM_ModelType modelType,
							ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod);

	/// For CSO bi currencies, it is ugly.....
	ARM_LocalCSOCalculator(const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& CpnCcy,
		const ARM_Currency& FundCcy,
		int CMSLong,
		int CMSShort,
		int cpnDayCount,
		int cpnFreq,
		int cpnResetTiming,
		int stubRule,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& leverageLong,
		const ARM_Curve& leverageShort,
		const ARM_Curve& cpnMin,
		const ARM_Curve& cpnMax,
		const ARM_Curve& strikes,
		int fundFreq,
		int fundDayCount,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,
		int exerciseFreq,
		int noticeGap,
		int payRec,
		bool switchFlag,
		int fundingType,
		size_t nbNCall,
		const ARM_Curve& fees,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
		vector<double> ModelDatas,
		const ARM_MarketData_ManagerRep& mktDataManager,
		ARM_ModelType		modelType,
		CalibrationType	calibrationType,
		vector<CalibStrikeType>	calibStrikeType,
		ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod = ARM_MarketIRModel::MONEYNESS,
		double moyenessLevel = 1.0);
	
	/// copy
	ARM_LocalCSOCalculator (const ARM_LocalCSOCalculator& rhs);
	/// destructor
    virtual ~ARM_LocalCSOCalculator () {};
	
	/// --------------------------------------------------
	/// ---- utilities	----------------------------------
	/// --------------------------------------------------
	ARM_LocalCSOCalculator& operator=(const ARM_LocalCSOCalculator&);
	virtual ARM_Object* Clone() const {return new ARM_LocalCSOCalculator (*this);};
	
	
	/// --------------------------------------------------
	/// ---- standard calculator methods -----------------
	/// --------------------------------------------------
	void					CreateEmptyCalibration();
	virtual void			CreateAndSetCalibration();
	virtual void			CreateAndSetModel();
	virtual void			UpdateModel();
	virtual void			UpdateCalibration(bool isUpdateStrike=true) ;
	virtual void			Calibrate();
    virtual double			Price();
	virtual ARM_RowInfo		ColumnNames() const;
	virtual ARM_RowInfo		MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure) const;
	virtual void			InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const ;
	virtual void			ComputePricingData() const;
	virtual void			CheckData()    {};
	virtual void			CheckMktData() {};

	void SetDiagStrikeToFrontier(ARM_StdPortfolioPtr& swaptionPF);

	/// --------------------------------------------------
	/// ---- Cst manager creation ------------------------
	/// --------------------------------------------------
	ARM_CstManagerPtr CreateCstManager ();

	
	/// --------------------------------------------------
	/// ---- Vector conversions --------------------------
	/// --------------------------------------------------
	/// computes attributes vectors from object ARM_Curves;
	void ComputeProductVectorsFromCurves() ;
	/// tells if cpn is fix rate
	bool IsFixCoupon (size_t cpnIdx, double& fixRate) const;


	/// --------------------------------------------------
	/// ---- Portfolios ----------------------------------
	/// --------------------------------------------------
	ARM_StdPortfolioPtr CreateSwaptionPortfolio(CalibrationType type, CalibStrikeType calibStrikeType,bool isUpdateStrike = true);
	void CreateSOPortfolio();
	void ComputeSOPortfolioPrices();
	void ComputeSwaptionPortfolioPrices();

	void CreateEmptyCalibrationFor2F();
	ARM_StdPortfolioPtr CreateSwaptionPortfolioFor2F(CalibStrikeType cst);
	ARM_StdPortfolioPtr CreateSOPortfolioFor2F(bool isSOCond,bool isLong=true);

	void ComputeCalibPortfolioPricesFor2F();
	void ComputeSoPricesFor2F(ARM_StdPortfolioPtr& pf,CalibStrikeType cst);
	void ComputeSwPricesFor2F(ARM_StdPortfolioPtr& pf);

	double ComputeResidualUnderlying(int exerIdx, double& floatLegPv, double& exoLegPv);

	/// accessors: calibration PF's
	ARM_StdPortfolioPtr GetSwaptionPF()				const ;
	ARM_StdPortfolioPtr GetFwdStartSwaptionPF()		const { return itsFwdStartSwaptionPF;};
	ARM_StdPortfolioPtr GetSpreadOptionFloorPF()	const { return itsSpreadOptionFloorPF;}
	ARM_StdPortfolioPtr GetSpreadOptionCapPF()		const { return itsSpreadOptionCapPF;}
	ARM_StdPortfolioPtr GetSOPortfolio()			const ;
	ARM_StdPortfolioPtr GetCMSLONGPortfolio()		const ;
	ARM_StdPortfolioPtr GetCMSSHORTPortfolio()		const ;

	void SetOSWPortfolio(const ARM_StdPortfolio& port);


	/// --------------------------------------------------
	/// ---- variable notional swaption calibration 
	/// --------------------------------------------------
	void ComputeSOSensitivities();
	void ComputeSONormalSensitivities();
	ARM_Swaption* CreateVarNotionalSwaptionAtExer(int exerIdx,  ARM_Swaption*& fwdStartSwaption, CalibStrikeType calibStrikeType, bool isUpdateStrike = true);
	ARM_VanillaSpreadOptionArg* CreateVanillaArgSO(ARM_SpreadOption* spreadOption);

		/// Vector
	inline virtual ARM_GP_Vector GetvCpnNominal() const  { return itsvFundNominal;};
	inline virtual ARM_GP_Vector GetvFundNominal() const { return itsvInitialFundNominal;};
	inline virtual ARM_GP_Vector GetvFundSpread() const { return itsvInitialFundSpread;};	

	/// switchable
	inline bool IsSwitchable() const	{ return itsSwitchFlag;		};
	inline int GetCMSFunding() const	{ return itsCMSFunding;		};

	void ComputeExtraSONormalSensitivitiesForSwitch();

private:
	/// --------------------------------------------------
	/// ---- Method Params -------------------------------
	/// --------------------------------------------------
	
	/// calibration: diagonal swaptions or variable notional swaptions
	CalibrationType itsCalibrationType;
	vector<CalibStrikeType> itsCalibStrikeType;
	
	/// HW : 1 or 2 factors
	ARM_ModelType itsModelType;

	/// ATM or MONEYNESS
	ARM_MarketIRModel::VnsPricingMethod itsVnsPricingMethod;

	/// Level of Moyness to calculate the smiled VolBasket 
	double itsMoyenessLevel;

	/// --------------------------------------------------
	/// ---- Vector conversions --------------------------
	/// --------------------------------------------------
	/// Convert ARM_Curves into relevant vectors
	/// this will be much easier for computations
	ARM_GP_Vector	itsvFundSpread;
	ARM_GP_Vector   itsvInitialFundSpread;
	ARM_GP_Vector   itsvInitialFundNominal;
	ARM_GP_Vector	itsvFundNominal;
	ARM_GP_Vector	itsvFundLeverage;
	size_t			itsFundSize;
	
	ARM_GP_Vector	itsvCpnNominal;
	ARM_GP_Vector	itsvCpnFloor;
	ARM_GP_Vector	itsvCpnCap;
	ARM_GP_Vector	itsvCpnLeverageShort;
    ARM_GP_Vector	itsvCpnLeverageLong;
	ARM_GP_Vector	itsvCpnStrike;
	ARM_BoolVector	itsvCpnIsCapped;	// computed, but could be interfaced
	ARM_BoolVector	itsvCpnIsFloored;	// computed, but could be interfaced
	size_t			itsCpnSize;
	
	ARM_GP_Vector	itsvExerFees;		// may not be used...
	ARM_BoolVector	itsvIsExerDate;
	ARM_IntVector	itsvFundIndex;
	ARM_IntVector	itsvCpnIndex;
	size_t			itsExerSize;
			

	/// --------------------------------------------------
	/// ---- Portfolios ----------------------------------
	/// --------------------------------------------------
	/// SO portfolios used for local model calibration
	//  and variable swaptions set computation
	//
	ARM_StdPortfolioPtr itsSpreadOptionFloorPF;
	ARM_StdPortfolioPtr itsSpreadOptionCapPF;
	ARM_StdPortfolioPtr itsFwdStartSwaptionPF;
	bool itsSoPricesComputed;
		

	/// --------------------------------------------------
	/// ---- Computed prices -----------------------------
	/// --------------------------------------------------
	/// double itsCSOPrice; /// already in base class
	double itsFundingPrice;
	double itsFixPrice;
	double itsFloorPrice;
	double itsCapPrice;
		


	/// --------------------------------------------------
	/// ---- variable notional swaption calibration 
	/// --------------------------------------------------
	// vanilla arg used to store CMS date strips
	ARM_VanillaArgPtrVector itsVanillaArgSOVect;

	// sensi of cpns w.r.t long & short cms rate
	ARM_GP_Vector itsCpnLongSensi;
	ARM_GP_Vector itsCpnShortSensi;
	
	// cpn value (no discount)
	ARM_GP_Vector itsCpnValue;

	/// --------------------------------------------------
	/// ---- Exercise probability -----------------------------
	/// --------------------------------------------------
	size_t			itsNbComputedProba;
	size_t			itsFirstComputedProba;
	ARM_VectorPtr	itsProba;
	bool			itsComputeProba;


	/// --------------------------------------------------
	/// ---- Switchable feature --------------------------
	/// --------------------------------------------------
	
	// Defining flags
	// --------------------------------------------------
	bool					itsSwitchFlag;
	int						itsCMSFunding;

	// For basket calibration
	// --------------------------------------------------
	ARM_StdPortfolioPtr		itsSpreadOptionFloorPFforSwitch;
	ARM_StdPortfolioPtr		itsSpreadOptionCapPFforSwitch;

	ARM_VanillaArgPtrVector itsVanillaArgSOVectForSwitch;

	ARM_GP_Vector			itsCpnLongSensiForSwitch;
	ARM_GP_Vector			itsCpnShortSensiForSwitch;
	ARM_GP_Vector			itsCpnValueForSwitch;
};

CC_END_NAMESPACE()

#endif //_INGPCALCULATORS_LOCALCSOCALCULATOR_H
