#ifndef _INGPMODELS_MODELSVMMSPREAD_H
#define _INGPMODELS_MODELSVMMSPREAD_H

#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamsvec.h"
#include "gpbase/curve.h"
#include "typedef.h"

#include "gpmodels/modelparamsbgmsv1f.h"
#include "gpmodels/svbgm.h"

#include "gpinfra/pricingmodelir.h"
#include "gpcalib/numerical.h"

#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/optimization1.h"
#include "gpnumlib/levmarq.h"
#include "gpnumlib/solver.h"
#include "gpclosedforms/inverse.h"
#include "gpnumlib/ran2.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_SVMMSpread : public ARM_PricingModelIR
{
protected:

	std::vector<double>					itsResetTimes;
	std::vector<double>					itsStartTimes;
	std::vector<double>					itsEndTimes1;
	std::vector<double>					itsEndTimes2;
	std::vector<double>					itsDelta;

	std::vector<double>					itsFwdRate;

	ARM_VectorVector				itsEigenValues;
	ARM_MatrixVector				itsModelLocalRealVar;

	ARM_VanillaSecDensityPtrVector	itsCalibSecDensities;
	ARM_IntVector					itsCalibWeight;

	int								itsNbEffectiveReset;
	int								itsSimpleEulerScheme;
	
	bool							itsCalibrationStatus;
	bool							itsCalibratedStatus;

public:

	ARM_SVMMSpread ( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params = NULL);
	
	ARM_SVMMSpread ( const ARM_SVMMSpread & rhs);
	
	virtual ~ARM_SVMMSpread ();

	virtual ARM_Object*				Clone() const	{ return new ARM_SVMMSpread (*this);};
	
	ARM_SVMMSpread &				operator = (const ARM_SVMMSpread & rhs);

public:
	// méthodes virtuelles pures
	virtual bool					ValidateModelParams(const ARM_ModelParams& params) const;
	virtual void					ValidateCalibMethod(ARM_CalibMethod& calibMethod);
	virtual void					SetNumeraire(const ARM_NumerairePtr& numerairePtr);

	
    virtual void					PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void					PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void					AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void					AdviseCurrentCalib(ARM_ModelFitter& modelFitter) {};
	virtual void					Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter) {};
	
    
	virtual bool					SupportBackwardInduction() const {	return false;}
    virtual bool					SupportForwardInduction()  const {	return true;}
	virtual bool					SupportAnalyticMarginal()  const {	return false;}

	virtual void					PostInit();
	
	virtual void					setCalibrationStatus(bool calibrationStatus);

	virtual void					setCalibratedStatus(bool calibratedStatus);

	virtual ARM_PricingStatesPtr	FirstPricingStates( size_t bucketSize ) const;

	virtual ARM_PricingStatesPtr	Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

	virtual ARM_PricingStatesPtr	Induct(ARM_PricingStatesPtr& states,double toTime);

	virtual std::vector<double>&			ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos);
	
	virtual ARM_VectorPtr			ComputeNumeraireTimes(const ARM_TimeInfoPtrVector& timeInfos) const;

	virtual void					AdviseBreakPointTimes(const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb = 0);

	//////////////////////////////
	//
	//		pour passer les échéanciers !
	//
	//////////////////////////////

	virtual void					setNumericalModelFitter(ARM_NumericalModelFitter*);

	//////////////////////////////
	//
	//		pour le monte carlo
	//
	//////////////////////////////

	virtual void					ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps);

	void							ModelStateLocalStdDev(const std::vector<double>& timeSteps, const ARM_MatrixVector& localVariances);

	virtual void					NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps );

	virtual void					ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const;

	virtual void					NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const;		

	virtual void					NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances ) const;		

	virtual bool					NeedsToCholeskyDecomposeFactors( ) const { return false; }					

	virtual void					MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	void							MCFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex, const ARM_GP_Matrix& x) const;

	// Pas utile !!
	virtual void					TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;
	virtual double					VarianceToTime(double var, double minTime=0.0, double maxTime=5*K_YEAR_LEN) const;


	///////////////////////
	//
	//		les prix
	//
	///////////////////////

	virtual ARM_VectorPtr			DiscountFactor(	const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr			Spread(
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

	virtual ARM_VectorPtr			VanillaCaplet(const string& curveName, 
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

	virtual ARM_VectorPtr			VanillaSwaption(const string& curveName,
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

	virtual ARM_VectorPtr			VanillaSpreadOptionLet(const string& curveName,
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
										double swapLongFloatStartTime,
										double swapLongFloatEndTime,
										const std::vector<double>& swapLongFixPayTimes,
										const std::vector<double>& swapLongFixPayPeriods,
										double swapShortFloatStartTime,
										double swapShortFloatEndTime,
										const std::vector<double>& swapShortFixPayTimes,
										const std::vector<double>& swapShortFixPayPeriods,
										const ARM_PricingStatesPtr& states) const;
	
	/// Standard ARM object support
	virtual string					toString(const string& indent="",const string& nextIndent="") const;

	void			SetNbEffectiveReset(int nb);

	int				GetNbEffectiveReset() const;

	const std::vector<double>&	GetResetTimes() const {return itsResetTimes;};
	const std::vector<double>&	GetStartTimes() const {return itsStartTimes;};
	const std::vector<double>&	GetEndTimes1() const {return itsEndTimes1;};
	const std::vector<double>&	GetEndTimes2() const {return itsEndTimes2;};
	const ARM_VectorVector&	GetEigenValues() const {return itsEigenValues;};

	virtual size_t	ModelStatesSize() const { return itsNbEffectiveReset + 1;}	
	
	void			CalcNbEffectiveReset(const std::vector<double>& timeSteps);

	std::vector<double>	GetATMGaussianVol();

	void			CorrectIthFwd(int i, double adj);

protected:
	
	double			getTerminalTime() const;

	bool			DoesResetTimeExist(double time) const;

	ARM_VectorPtr	cmsSpread(const string& curveName, double evalTime, double startTime, double endTime1, double endTime2, const ARM_PricingStatesPtr& states) const;
};

inline void ARM_SVMMSpread::PostInit()
{
}

inline void ARM_SVMMSpread::setCalibrationStatus(bool calibrationStatus)
{
	itsCalibrationStatus = calibrationStatus;
}

inline void ARM_SVMMSpread::setCalibratedStatus(bool calibratedStatus)
{
	itsCalibratedStatus = calibratedStatus;
}

inline void ARM_SVMMSpread::SetNbEffectiveReset(int nb)
{
	itsNbEffectiveReset = nb;
}

inline int ARM_SVMMSpread::GetNbEffectiveReset() const
{
	return itsNbEffectiveReset;
}

inline void ARM_SVMMSpread::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_SVMMSpread::TreeStatesToModelStates :  not implemented ARM_SVMMSpread Model!");
}

inline double ARM_SVMMSpread::VarianceToTime(double var,double minTime,double maxTime) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_SVMMSpread::VarianceToTime not implemented ARM_SVMMSpread Model!");

}

inline double ARM_SVMMSpread::getTerminalTime() const
{
	return itsEndTimes1.Elt( itsEndTimes1.size() -1 );
}

inline bool ARM_SVMMSpread::DoesResetTimeExist(double time) const
{
	return ExistsInVector( itsResetTimes, time );
}

inline void ARM_SVMMSpread::CorrectIthFwd(int i, double adj)
{
	itsFwdRate[i] += adj;
}

CC_END_NAMESPACE()

#endif
