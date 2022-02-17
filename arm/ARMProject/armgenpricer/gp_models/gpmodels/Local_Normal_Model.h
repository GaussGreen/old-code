/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Local_Normal_Model.h
 *
 *  \brief base class for local analytical normal model
 *	\author  A. Chaix
 *	\version 1.0
 *	\date June 2005
 */


#ifndef _INGPMODELS_LOCAL_NORMAL_MODEL_H
#define _INGPMODELS_LOCAL_NORMAL_MODEL_H

/// gp bases
#include "gpbase/env.h"
#include "gpbase/port.h"

/// gp models
#include "Local_Model.h"
#include "typedef.h"

class ARM_SpreadOption;
class ARM_CapFloor;
class ARM_CorridorDblCondition;


CC_BEGIN_NAMESPACE( ARM )

class ARM_Local_Normal_ModelParams;

class ARM_Local_Normal_Model : public ARM_Local_Model
{

// ---- attributes
private:
	
	///----------------------------------------------------------------------
	///-- nb stev above which we do not try to calibrate forward vol
	///----------------------------------------------------------------------
	static const double NSTDEV_NO_CALIB;

	///-----------------------------------------------
	///--- Indexes in surface lists where
	///--- are stored local vols & adjustments
	///-----------------------------------------------
	///--------------------------------------------
	/// Libor & VanillaCaplet
	static const size_t CAPLET_ADJ_SIZE;
	static const size_t CAPLET_ADJ;
	static const size_t CAPLET_VOL_SIZE;
	static const size_t CAPLET_VOL;
	///--------------------------------------------
	/// VanillaSpreadOptionlet
	static const size_t SO_ADJ_SIZE;
	static const size_t SO_ADJ_LONG;
	static const size_t SO_ADJ_SHORT;
	static const size_t SO_VOL_SIZE;
	static const size_t SO_VOL;
	///--------------------------------------------
	/// VanillaCMSCorridorlet (single condition)
public:
	static const size_t SOCOR_DOWN;
	static const size_t SOCOR_UP;

	static const size_t SOCOR_ADJ_SIZE;
private:
	static const size_t SOCOR_ADJ_LONG;
	static const size_t SOCOR_ADJ_SHORT;
	static const size_t SOCOR_ADJ_LONG_UP_FLT;
	static const size_t SOCOR_ADJ_SHORT_UP_FLT;
	static const size_t SOCOR_ADJ_LONG_DOWN_FLT;
	static const size_t SOCOR_ADJ_SHORT_DOWN_FLT;
	static const size_t SOCOR_ADJ_PAY_FLT;
public:
	static const size_t SOCOR_VOL_SIZE;
private:
	static const size_t SOCOR_VOL_UPLEFT;
	static const size_t SOCOR_VOL_UPRIGHT;
	static const size_t SOCOR_VOL_DOWNLEFT;
	static const size_t SOCOR_VOL_DOWNRIGHT;
	static const size_t SOCOR_VOL_UPLEFT_FLT;
	static const size_t SOCOR_VOL_UPRIGHT_FLT;
	static const size_t SOCOR_VOL_DOWNLEFT_FLT;
	static const size_t SOCOR_VOL_DOWNRIGHT_FLT;

	static const double SOCOR_STRIKE_SPREAD;

public:
	static const size_t SOCOR_KADJ_SIZE;
private:
	static const size_t SOCOR_KADJ_UPLEFT;
	static const size_t SOCOR_KADJ_UPRIGHT;
	static const size_t SOCOR_KADJ_DOWNLEFT;
	static const size_t SOCOR_KADJ_DOWNRIGHT;
	static const size_t SOCOR_KADJ_UPLEFT_FLT;
	static const size_t SOCOR_KADJ_UPRIGHT_FLT;
	static const size_t SOCOR_KADJ_DOWNLEFT_FLT;
	static const size_t SOCOR_KADJ_DOWNRIGHT_FLT;

	static const size_t SOCOR_MKT_RESIDUAL_RSO_PRICE;
	static const size_t SOCOR_RESIDUAL_RSO_PV_CORRECTION;
	///--------------------------------------------

	///--------------------------------------------
	/// VanillaCMSCorridorlet (double condition)
	/// Corridor types
public:
	static const size_t DBLECOR_SDOWN_RDOWN;
	static const size_t DBLECOR_SDOWN_RUP;
	static const size_t DBLECOR_SUP_RDOWN;
	static const size_t DBLECOR_SUP_RUP;

	static const size_t DBLECOR_NB_BARRIERS;

	/// Adj rates
public:
	static const size_t DBLECOR_ADJ_SIZE;

private:
	static const size_t DBLECOR_ADJ_SPREAD;
	static const size_t DBLECOR_ADJ_RATE;

	static const size_t DBLECOR_SDOWN_RDOWN_ADJ_SPREAD_FLT;
	static const size_t DBLECOR_SDOWN_RUP_ADJ_SPREAD_FLT;
	static const size_t DBLECOR_SUP_RDOWN_ADJ_SPREAD_FLT;
	static const size_t DBLECOR_SUP_RUP_ADJ_SPREAD_FLT;

	static const size_t DBLECOR_SDOWN_RDOWN_ADJ_RATE_FLT;
	static const size_t DBLECOR_SDOWN_RUP_ADJ_RATE_FLT;
	static const size_t DBLECOR_SUP_RDOWN_ADJ_RATE_FLT;
	static const size_t DBLECOR_SUP_RUP_ADJ_RATE_FLT;

	static const size_t DBLECOR_ADJ_PAY_FLT;

	/// Adj strikes
public:
	static const size_t DBLECOR_KADJ_SIZE;

private:
	static const size_t DBLECOR_SDOWN_RDOWN_KADJ_SPREAD;
	static const size_t DBLECOR_SDOWN_RUP_KADJ_SPREAD;
	static const size_t DBLECOR_SUP_RDOWN_KADJ_SPREAD;
	static const size_t DBLECOR_SUP_RUP_KADJ_SPREAD;

	static const size_t DBLECOR_SDOWN_RDOWN_KADJ_RATE;
	static const size_t DBLECOR_SDOWN_RUP_KADJ_RATE;
	static const size_t DBLECOR_SUP_RDOWN_KADJ_RATE;
	static const size_t DBLECOR_SUP_RUP_KADJ_RATE;

	static const size_t DBLECOR_SDOWN_RDOWN_KADJ_SPREAD_FLT;
	static const size_t DBLECOR_SDOWN_RUP_KADJ_SPREAD_FLT;
	static const size_t DBLECOR_SUP_RDOWN_KADJ_SPREAD_FLT;
	static const size_t DBLECOR_SUP_RUP_KADJ_SPREAD_FLT;

	static const size_t DBLECOR_SDOWN_RDOWN_KADJ_RATE_FLT;
	static const size_t DBLECOR_SDOWN_RUP_KADJ_RATE_FLT;
	static const size_t DBLECOR_SUP_RDOWN_KADJ_RATE_FLT;
	static const size_t DBLECOR_SUP_RUP_KADJ_RATE_FLT;

	/// Adj vols
public:
	static const size_t DBLECOR_VOL_SIZE;

private:
	static const size_t DBLECOR_SDOWN_RDOWN_VOL_SPREAD;
	static const size_t DBLECOR_SDOWN_RUP_VOL_SPREAD;
	static const size_t DBLECOR_SUP_RDOWN_VOL_SPREAD;
	static const size_t DBLECOR_SUP_RUP_VOL_SPREAD;

	static const size_t DBLECOR_SDOWN_RDOWN_VOL_RATE;
	static const size_t DBLECOR_SDOWN_RUP_VOL_RATE;
	static const size_t DBLECOR_SUP_RDOWN_VOL_RATE;
	static const size_t DBLECOR_SUP_RUP_VOL_RATE;

	static const size_t DBLECOR_SDOWN_RDOWN_VOL_SPREAD_FLT;
	static const size_t DBLECOR_SDOWN_RUP_VOL_SPREAD_FLT;
	static const size_t DBLECOR_SUP_RDOWN_VOL_SPREAD_FLT;
	static const size_t DBLECOR_SUP_RUP_VOL_SPREAD_FLT;

	static const size_t DBLECOR_SDOWN_RDOWN_VOL_RATE_FLT;
	static const size_t DBLECOR_SDOWN_RUP_VOL_RATE_FLT;
	static const size_t DBLECOR_SUP_RDOWN_VOL_RATE_FLT;
	static const size_t DBLECOR_SUP_RUP_VOL_RATE_FLT;

	/// Adj correls
public:
	static const size_t DBLECOR_CORREL_SIZE;

private:
	static const size_t DBLECOR_SDOWN_RDOWN_SPREAD_RATE_CORREL;
	static const size_t DBLECOR_SDOWN_RUP_SPREAD_RATE_CORREL;
	static const size_t DBLECOR_SUP_RDOWN_SPREAD_RATE_CORREL;
	static const size_t DBLECOR_SUP_RUP_SPREAD_RATE_CORREL;

	static const size_t DBLECOR_SDOWN_RDOWN_SPREAD_RATE_CORREL_FLT;
	static const size_t DBLECOR_SDOWN_RUP_SPREAD_RATE_CORREL_FLT;
	static const size_t DBLECOR_SUP_RDOWN_SPREAD_RATE_CORREL_FLT;
	static const size_t DBLECOR_SUP_RUP_SPREAD_RATE_CORREL_FLT;

	static const size_t DBLECOR_MKT_RESIDUAL_RA2_PRICE;
	static const size_t DBLECOR_RESIDUAL_RA2_PV_CORRECTION;

	static const size_t ARM_Local_Normal_Model::DBLECOR_MARKET_CORREL;

	///--------------------------------------------

	/// Gaussian correlations for copula based calibration
	static const size_t COPULA_CORREL_SIZE;
	static const size_t COPULA_MARKET_CORREL;
	static const size_t COPULA_MODEL_CORREL;


	bool itsLiborAsShiftedLognormal;

	ARM_IntVector itsParamTypes;
	ARM_IntVector itsParamNbSurfaces;

	virtual void CalibrateSpreadOptionLet (ARM_SpreadOption* spreadOption, size_t periodIdx, double targetPrice, const std::vector<double>& evalTimes);
	virtual void CalibrateSpreadOption (ARM_SpreadOption* spreadOption, double targetPrice, const std::vector<double>& evalTimes);
	virtual void CalibrateCorridorSpreadOption (ARM_SpreadOption* spreadOption, double targetPrice, const std::vector<double>& evalTimes,size_t barrierType);
	virtual void CalibrateCapFloorLet (ARM_CapFloor* capFloor, size_t periodIdx, double targetPrice, const std::vector<double>& evalTimes);
	virtual void CalibrateCapFloor (ARM_CapFloor* capFloor,  double targetPrice, const std::vector<double>& evalTimes);
	virtual void CalibrateDoubleCorridorOption(ARM_CorridorDblCondition* dblCorridor, double targetPrice, const std::vector<double>& evalTimes, size_t barrierType);

	void ComputeDoubleCorridorDigital(
					ARM_VectorPtr& spreads,
					ARM_VectorPtr& rates,
					std::vector<double>& spreadAdj,
					std::vector<double>& rateAdj,
					ARM_VectorPtr& payedRates,
					double spreadDownBarrier,
					double spreadUpBarrier,
					double rateDownBarrier,
					double rateUpBarrier,
					double evalTime,
					double resetTime,
					double sqrTtoE,
					double coefCoupon,
					size_t KAdjSpreadSdRd,
					size_t KAdjSpreadSdRu,
					size_t KAdjSpreadSuRd,
					size_t KAdjSpreadSuRu,
					size_t KAdjRateSdRd,
					size_t KAdjRateSdRu,
					size_t KAdjRateSuRd,
					size_t KAdjRateSuRu,
					size_t VolSpreadSdRd,
					size_t VolSpreadSdRu,
					size_t VolSpreadSuRd,
					size_t VolSpreadSuRu,
					size_t VolRateSdRd,
					size_t VolRateSdRu,
					size_t VolRateSuRd,
					size_t VolRateSuRu,
					size_t CorrelSpreadRateSdRd,
					size_t CorrelSpreadRateSdRu,
					size_t CorrelSpreadRateSuRd,
					size_t CorrelSpreadRateSuRu,
					ARM_VectorPtr& values) const;

public:
	///constructors/destructors
	ARM_Local_Normal_Model (const ARM_ZeroCurvePtr& zc /* A GARDER ??*/, const ARM_Local_Normal_ModelParams& params);
	ARM_Local_Normal_Model(const ARM_ZeroCurvePtr& zc, const ARM_Local_Normal_ModelParams& params,
						const ARM_IntVector& paramTypes, const ARM_IntVector& paramNbSurfaces);
	ARM_Local_Normal_Model (const ARM_Local_Normal_Model& rhs);
    ARM_Local_Normal_Model& operator = (const ARM_Local_Normal_Model& rhs);
	virtual ~ARM_Local_Normal_Model();

	/// utilities
	virtual ARM_Object* Clone() const {return new ARM_Local_Normal_Model (*this);}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LLOCN";}

	virtual void ResetModelParams ();

	static ARM_Local_Normal_ModelParams* CreateDefaultModelParams ();
	static ARM_ModelParam*  CreateDefaultForwardAdjustmentModelParam ();
	static ARM_ModelParam*  CreateDefaultVolatilityModelParam ();
	static ARM_ModelParam*  CreateDefaultShiftModelParam ();

	static ARM_Local_Normal_ModelParams* CreateDefaultModelParams(const ARM_IntVector& ParamTypes, const ARM_IntVector& ParamNbSurfaces);
	static ARM_ModelParam*  CreateDefaultModelParam(ARM_ModelParamType::ParamNb paramType, int nbSurfaces);

	///accessors
	inline const bool& LiborAsShiftedLognormal() const {return itsLiborAsShiftedLognormal;}
	inline void  ComputeLiborAsShitedLognormal() {itsLiborAsShiftedLognormal = true;}
	inline void  ComputeLiborNormal()			{itsLiborAsShiftedLognormal = false;}


	///services
	struct CorrelUnSqueezerData
	{
		double itsVarModelX;
		double itsVarModelY;
		double itsCovModelXY;
		double itsCorTargetXY;

		double itsLocalCorTargetXY;

		double itsInitX;
		double itsInitY;

		CorrelUnSqueezerData(double varModelX,double varModelY,double covModelXY,double corTargetXY,
			double localCorTargetXY,double initX,double initY)
			: itsVarModelX(varModelX),itsVarModelY(varModelY),itsCovModelXY(covModelXY),
			itsCorTargetXY(corTargetXY),itsLocalCorTargetXY(localCorTargetXY),
			itsInitX(initX),itsInitY(initY) {}
		double ComputeErr(double* x)
		{
			double localCorXY = (fabs(x[0]*x[1])*itsCorTargetXY - itsCovModelXY)
								/ sqrt((x[0]*x[0]-itsVarModelX)*(x[1]*x[1]-itsVarModelY));
			double errCor = localCorXY - itsLocalCorTargetXY;
			double errX = x[0]/itsInitX - 1;
			double errY = x[1]/itsInitY - 1;
			return errCor*errCor + 0.3*(errX*errX + errY*errY);
		}
	};

	static bool CorrelUnSqueezer(double Tcall, double Tfix,
		double spreadTarget,double rateTarget,double spreadRateCorrelTarget,
		double spreadVarModel,double rateVarModel,double spreadRateCovModel,double sTfixTcall,
		double &newSpreadKAdj,double &newRateKAdj,double &newSpreadStdDev,double &newRateStdDev,
		double &newLocalSpreadStdDev,double &newLocalRateStdDev,double &newLocalSpreadRateCorrel,
		size_t& iter);

	static bool VolUnSqueezer(
		double fwd, double targetPrice,double strike,
		double mat, int callPut, double targetVol,
		double &adj);
	
	//------------------------------
	//--- local model calibration
	//------------------------------
	virtual void CalibrateLocalModel (const ARM_Security& security, double targetPrice, const std::vector<double>& evalTimes, size_t secIdx=0);
	virtual void CalibrateLocalModel (const ARM_Security& security,	ARM_DensityFunctor& density,bool rescaling);
	
	virtual ARM_VectorPtr Func(double evalTime,const ARM_GP_VectorPtr& values) const;
	//------------------------------
	//--- implemented IR functions
	//------------------------------
	
	// libor
	virtual ARM_VectorPtr Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const ;


	// caplet
    virtual ARM_VectorPtr VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
		const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const;

	// swaption
	virtual ARM_VectorPtr VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const std::vector<double>& fixNotional,
		const std::vector<double>& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;

	// spread option
	virtual ARM_VectorPtr  VanillaSpreadOptionLet(
			const string& curveName,
			double evalTime,
			int callPut,
			double startTime,
			double endTime,
			double resetTime,
			double payTime,
			double payPeriod,
			double notional,
			double coeffLong,
			double coeffShort,
			const std::vector<double>& strikes,
			double swapLongFloatStartTime,	// adjusted ...
			double swapLongFloatEndTime,	// adjusted ...
			const std::vector<double>& swapLongFixPayTimes,
			const std::vector<double>& swapLongFixPayPeriods,
			double swapShortFloatStartTime,	// adjusted ...
			double swapShortFloatEndTime,	// adjusted ...
			const std::vector<double>& swapShortFixPayTimes,
			const std::vector<double>& swapShortFixPayPeriods,
			const ARM_PricingStatesPtr& states) const;

	// Vanilla Corridorlet function
	virtual ARM_VectorPtr VanillaCMSCorridorlet(
		const string& curveName,
		double evalTime,
		double payTime,
		double resetTime,
		double startTime,
		double endTime,
		const std::vector<double>& refIdxResettimes,
		const std::vector<double>& refIndexWeights,
		const std::vector<double>& coeff1,
		const ARM_SwapRatePtrVector& firstIndex,
		const std::vector<double>& coeff2,
		const ARM_SwapRatePtrVector& secondIndex,
		int		payIndexType,			/// K_FIXED, K_LIBOR or K_CMS
		double	coupon,					/// in case of a fixed payment (K_FIXED)
		const	ARM_SwapRate& payRate,	/// rate description (K_LIBOR or K_CMS)
		double payIndexLeverage,
		const std::vector<double>& downBarriers,
        const std::vector<double>& upBarriers,
        double  payNotional,
        int     rcvPay,
		const ARM_SwapRatePtrVector& thirdIndex, // 3rd index for double condition
		const std::vector<double>& downBarriers3,
		const std::vector<double>& upBarriers3,
        const   ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr Spread(
		const string& curveName, 
        double evalTime,
		double coeff1,
		double floatStartTime1,
        double floatEndTime1, 
		const std::vector<double>& fixPayTimes1,
        const std::vector<double>& fixPayPeriods1,
		const std::vector<double>& fwdStartTimes1,
        const std::vector<double>& fwdEndTimes1,
        const std::vector<double>& fwdPayPeriods1, 
		const std::vector<double>& floatPayTimes1,
        const std::vector<double>& floatPayPeriods1,
        const std::vector<double>& margin1,
		double coeff2,
        double floatStartTime2,
        double floatEndTime2, 
		const std::vector<double>& fixPayTimes2,
        const std::vector<double>& fixPayPeriods2,
		const std::vector<double>& fwdStartTimes2,
        const std::vector<double>& fwdEndTimes2,
        const std::vector<double>& fwdPayPeriods2, 
		const std::vector<double>& floatPayTimes2,
        const std::vector<double>& floatPayPeriods2,
        const std::vector<double>& margin2,
        const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr DoubleDigital(
		const string& modelName, 
		double evalTime,
		const ARM_VectorPtr& firstRate,
        const std::vector<double>& firstStrikeDown,
        const std::vector<double>& firstStrikeUp,
		double firstStrikeSpread,
		const ARM_VectorPtr& secondRate,
        const std::vector<double>& secondStrikeDown,
        const std::vector<double>& secondStrikeUp,
		double secondStrikeSpread,
        const ARM_PricingStatesPtr& states) const;

	// model params validation
	virtual void PreProcessing(ARM_ModelFitter& modelFitter) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter) {}
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter) {};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter) {};
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;
};

///////////////////////////////////////////////////////////////////////////////
// idea : helper for the calibration of z local normal model
//
// --> computation of normal vol + cms adjustements under a numerical model
// ideally, these services should be moved on the numerical model side... 
///////////////////////////////////////////////////////////////////////////////
class Local_Normal_Model_Calibration_Helper
{
private:
	ARM_PricingModel* itsNumericalModel;
	ARM_ModelType itsModelType;

public:
	Local_Normal_Model_Calibration_Helper (ARM_PricingModel* numericalModel);
	Local_Normal_Model_Calibration_Helper (const Local_Normal_Model_Calibration_Helper& rhs);
	virtual ~Local_Normal_Model_Calibration_Helper () {};
	
	// LIBOR rate under num. model, fixed @ expiryTime & payed @ payTime ( = E^Q(Tp)[S(Texp)] )
	double LiborRateFromNumericalModel(
									double expiryTime,
									double payTime, 
									double fwdStartTime,
									double fwdEndTime,
									double fwdPeriod,
									double strike ) const;
	// normal vol of a libor on [0, expiryTime]
	double LiborNormalVolatilityFromNumericalModel(
									double expiryTime,
									//double payTime, // not used in HW case
									double fwdStartTime,
									double fwdEndTime,
									double fwdPeriod,
									double strike ) const;

	// lognormal vol of B(.,Tstart)/B(.,Tend) on [0, expiryTime]
	double FwdZcVolatilityFromNumericalModel(
									double expiryTime,
									//double payTime, // not used in HW case
									double fwdStartTime,
									double fwdEndTime,
									double fwdPeriod,
									double strike ) const;


	
	// CMS rate under num. model, fixed @ expiryTime & payed @ payTime ( = E^Q(Tp)[S(Texp)] )
	double CmsRateFromNumericalModel(
									double expiryTime,
									double payTime, 
									double swapFloatStartTime,	// adjusted ...
									double swapFloatEndTime,	// adjusted ...
									const std::vector<double>& swapFixPayTimes,
									const std::vector<double>& swapFixPayPeriods,
									double* factor1=NULL,double* factor2=NULL) const;
									
	
	// normal vol of a CMS spread on [0, exprityTime]
	double SpreadNormalVolatilityFromNumericalModel(
									double expiryTime,
									double payTime, // not used in HW case
									double coeffLong,
									double coeffShort,
									double strike,
									double swapLongFloatStartTime,	// adjusted ...
									double swapLongFloatEndTime,	// adjusted ...
									const std::vector<double>& swapLongFixPayTimes,
									const std::vector<double>& swapLongFixPayPeriods,
									double swapShortFloatStartTime,
									double swapShortFloatEndTime,
									const std::vector<double>& swapShortFixPayTimes,
									const std::vector<double>& swapShortFixPayPeriods,
									double* factor1=NULL,double* factor2=NULL) const;

	/// CMS convexified rate when it is attached to a floating payment rate
	double FLTconvexifiedRateFromNumericalModel(
									double cms, 
									double expiryTime,
									double swapFloatStartTime,	// adjusted ...
									double swapFloatEndTime,	// adjusted ...
									const std::vector<double>& swapFixPayTimes,
									const std::vector<double>& swapFixPayPeriods,
									/// -- payed rate
									double payIndexFloatStartTime,	// adjusted ...
									double payIndexFloatEndTime,	// adjusted ...
									const std::vector<double>& payIndexFixPayTimes,
									const std::vector<double>& payIndexFixPayPeriods) const;

	/// Integrated normal covariance between a CMS spread & a rate from 0 to expiryTime
	/// Extended with optional parameters to compute CMS spread vs CMS spread normal covariance
	double SpreadRateNormalCovarianceFromNumericalModel(
									double expiryTime,
									double payTime,
									double coeffLong,
									double coeffShort,
									double strike,
									double swapLongFloatStartTime,
									double swapLongFloatEndTime,
									const std::vector<double>& swapLongFixPayTimes,
									const std::vector<double>& swapLongFixPayPeriods,
									double swapShortFloatStartTime,
									double swapShortFloatEndTime,
									const std::vector<double>& swapShortFixPayTimes,
									const std::vector<double>& swapShortFixPayPeriods,
									double rateStrike,
									double rateFloatStartTime,
									double rateFloatEndTime,
									const std::vector<double>& rateFixPayTimes,
									const std::vector<double>& rateFixPayPeriods,
									double rateCoeffLong=1.0,
									double rateCoeffShort=0.0,
									double rateShortFloatStartTime=0.0,
									double rateShortFloatEndTime=0.0,
									const std::vector<double>& rateShortFixPayTimes=std::vector<double>(0),
									const std::vector<double>& rateShortFixPayPeriods=std::vector<double>(0)) const;
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
