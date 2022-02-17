/*!
 *
 * Copyright (c) IXIS CIB July 2005 Paris
 *
 *	\file Local_SLN_Model.h
 *
 *  \brief class for Shifted LogNormal local model
 *	\author  J-M Prié
 *	\version 1.0
 *	\date July 2005
 */


#ifndef _INGPMODELS_LOCAL_SLN_MODEL_H
#define _INGPMODELS_LOCAL_SLN_MODEL_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "local_model.h"

#define UNIMPLEMENTED_SLN_PRICING_FUNCTION  { ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented pricing function in Local_SLN_Model"); return ARM_VectorPtr(); }

/// Forward declaration outside ARM namespace
class ARM_CapFloor;
class ARM_Option;
class ARM_CorridorLeg;


CC_BEGIN_NAMESPACE( ARM )

/// Forward declaration inside ARM namespace
class ARM_ModelParams;
class ARM_ModelParam;
struct ARM_PricingContext;
class ARM_PricingFunctionIR;
class ARM_PricingFunctionEquity;
class ARM_Local_SLN_ModelParams;

////////////////////////////////////////////////////
///
///	ARM_Local_SLN_Model class creates a local model
/// assuming a Shifted LogNormal diffusion
///
////////////////////////////////////////////////////
class ARM_Local_SLN_Model : public ARM_Local_Model
{
private:
	virtual void CalibrateCapFloorLet (ARM_CapFloor* capFloor, size_t periodIdx, double targetPrice, const std::vector<double>& evalTimes);
	virtual void CalibrateCapFloor (ARM_CapFloor* capFloor,  double targetPrice, const std::vector<double>& evalTimes);
	virtual void CalibrateEqFxOption (ARM_Option* option, double targetPrice, const std::vector<double>& evalTimes);
	virtual void CalibrateCorridorLet (ARM_CorridorLeg* corridor, double targetPrice, const std::vector<double>& evalTimes);
	ARM_GP_MatrixPtr itsStatisticCalibResult;
	bool itsResizeStatisticResult;


public:
    // ---- constructors/destructors
	ARM_Local_SLN_Model (const ARM_Local_SLN_ModelParams& params);
	ARM_Local_SLN_Model (const ARM_Local_SLN_Model& rhs);
    ARM_Local_SLN_Model& operator = (const ARM_Local_SLN_Model& rhs);
	virtual ~ARM_Local_SLN_Model();

    // ---- utilities
	virtual ARM_Object* Clone() const {return new ARM_Local_SLN_Model (*this);}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LLSLN";}

	virtual void ResetModelParams ();

	static ARM_ModelParams* CreateDefaultModelParams();
    static ARM_ModelParam*  CreateDefaultModelParam(ARM_ModelParamType::ParamNb paramType);

	/// Model params validation
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	virtual ARM_Date GetAsOfDate() const { return GetNumericalModel()->GetAsOfDate(); }

	//assessor
	inline const ARM_GP_MatrixPtr& GetStatisticCalibResult() const {return itsStatisticCalibResult;}
    inline ARM_GP_MatrixPtr& GetStatisticCalibResult() {return itsStatisticCalibResult;}
	inline void SetStatisticCalibResult(const ARM_GP_MatrixPtr& statisticMatrix) {itsStatisticCalibResult=statisticMatrix;}

	//------------------------------
	//--- local model calibration
	//------------------------------
	virtual void CalibrateLocalModel (const ARM_Security& security, double targetPrice, const std::vector<double>& evalTimes, size_t secIdx=0);
	virtual void CalibrateLocalModel (const ARM_GenCalculator& calculator, const std::vector<double>& evalTimes);


	//------------------------------
	//--- IR functions
	//------------------------------
	
	// libor
	virtual ARM_VectorPtr Libor( 
		const string& modelName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const;


	// caplet
    virtual ARM_VectorPtr VanillaCaplet(
		const string& modelName, 
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

	//Corridorlet
	virtual ARM_VectorPtr VanillaCorridorlet(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     indexPaymentType,
        double  fwdPaymentPeriod,
        const std::vector<double>& RefIdxResettimes,
        const std::vector<double>& RefIdxStarttimes,
        const std::vector<double>& RefIdxEndtimes,
        const std::vector<double>& RefFwdPeriods,
        const std::vector<double>& RefIndexWeight,
        double  couponMargin,
        const vector<const std::vector<double>*> downBarrierPerState,
        const vector<const std::vector<double>*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const;   


	// swaption
	virtual ARM_VectorPtr VanillaSwaption(
		const string& modelName,
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
		bool isConstantStrike = true) const UNIMPLEMENTED_SLN_PRICING_FUNCTION;

	// spread option
	virtual ARM_VectorPtr  VanillaSpreadOptionLet(
			const string& modelName,
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
			const ARM_PricingStatesPtr& states) const UNIMPLEMENTED_SLN_PRICING_FUNCTION;

	//----------------------------------
	//--- EQ/FX functions
	//----------------------------------
	
	/// Forward function
	virtual ARM_VectorPtr Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const;

	/// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const std::vector<double>& strikesPerState,
        int callPut,
	    double payTime,
	    const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context) const;

	/// Convention support for equity/fx markets
    virtual string GetSettlementCalendar(const string& modelName="") const;
    virtual double GetSettlementGap(const string& modelName="") const;
};


///////////////////////////////////////////////////////////////////////////////
// idea : helper for the calibration of z local normal model
//
// --> computation of normal vol + cms adjustements under a numerical model
// ideally, these services should be moved on the numerical model side... 
///////////////////////////////////////////////////////////////////////////////
class Local_SLN_Model_Calibration_Helper
{
private:
	ARM_PricingModel* itsNumericalModel;
	int itsModelType;
	static int SFRM;
	
public:
	Local_SLN_Model_Calibration_Helper (ARM_PricingModel* numericalModel);
	Local_SLN_Model_Calibration_Helper (const Local_SLN_Model_Calibration_Helper& rhs);
	virtual ~Local_SLN_Model_Calibration_Helper () {};
	
public:
	
	// LIBOR rate under num. model, fixed @ expiryTime & payed @ payTime ( = E^Q(Tp)[S(Texp)] )
	virtual double LiborRateFromNumericalModel(
									double evalTime,									
									double payTime, 
									double fwdResetTime,
									double fwdStartTime,
									double fwdEndTime,
									double fwdPeriod,
									double strike ) const;
	// normal vol of a libor on [0, expiryTime]
	virtual double LiborSLNVolatilityFromNumericalModel(
									double expiryTime,
									double fwdResetTime , 
									double fwdStartTime,
									double fwdEndTime,
									double fwdPeriod,
									double strike,
									double shift /*not used for SFRM*/) const;

	virtual double LiborSLNShiftFromNumericalModel(
									double expiryTime,
									double fwdResetTime , 
									double fwdStartTime,
									double fwdEndTime,
									double fwdPeriod) const;

	virtual double CovarFromNumericalModel(
								double evalTime,
								double resetTime,
								double fwdResetTime1 , 
								double fwdResetTime2) const;
};

CC_END_NAMESPACE()


#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
