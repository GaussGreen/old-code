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

	ARM_GP_Vector					itsResetTimes;
	ARM_GP_Vector					itsStartTimes;
	ARM_GP_Vector					itsEndTimes1;
	ARM_GP_Vector					itsEndTimes2;
	ARM_GP_Vector					itsDelta;

	ARM_GP_Vector					itsFwdRate;

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

	virtual ARM_GP_Vector*			ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos);
	
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

	virtual void					ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps);

	void							ModelStateLocalStdDev(const ARM_GP_Vector& timeSteps, const ARM_MatrixVector& localVariances);

	virtual void					NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps );

	virtual void					ModelStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances ) const;

	virtual void					NumMethodStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances ) const;		

	virtual void					NumMethodStateGlobalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& globalVariances ) const;		

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
										const ARM_GP_Vector& fixPayTimes1,
										const ARM_GP_Vector& fixPayPeriods1,
										const ARM_GP_Vector& fwdStartTimes1,
										const ARM_GP_Vector& fwdEndTimes1,
										const ARM_GP_Vector& fwdPayPeriods1, 
										const ARM_GP_Vector& floatPayTimes1,
										const ARM_GP_Vector& floatPayPeriods1,
										const ARM_GP_Vector& margin1,
										double coeff2,
										double floatStartTime2,
										double floatEndTime2, 
										const ARM_GP_Vector& fixPayTimes2,
										const ARM_GP_Vector& fixPayPeriods2,
										const ARM_GP_Vector& fwdStartTimes2,
										const ARM_GP_Vector& fwdEndTimes2,
										const ARM_GP_Vector& fwdPayPeriods2, 
										const ARM_GP_Vector& floatPayTimes2,
										const ARM_GP_Vector& floatPayPeriods2,
										const ARM_GP_Vector& margin2,
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
										const ARM_GP_Vector& strikesPerState,
										int capFloor,
										const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr			VanillaSwaption(const string& curveName,
										double evalTime,
										double swapResetTime,
										const ARM_GP_Vector& fixNotional,
										const ARM_GP_Vector& floatNotional,
										double floatStartTime,
										double floatEndTime,
										const ARM_GP_Vector& floatResetTimes,
										const ARM_GP_Vector& floatStartTimes,
										const ARM_GP_Vector& floatEndTimes,
										const ARM_GP_Vector& floatIntTerms,
										const ARM_GP_Vector& fixPayTimes,
										const ARM_GP_Vector& fixPayPeriods,
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
										const ARM_GP_Vector& strikes,
										double swapLongFloatStartTime,
										double swapLongFloatEndTime,
										const ARM_GP_Vector& swapLongFixPayTimes,
										const ARM_GP_Vector& swapLongFixPayPeriods,
										double swapShortFloatStartTime,
										double swapShortFloatEndTime,
										const ARM_GP_Vector& swapShortFixPayTimes,
										const ARM_GP_Vector& swapShortFixPayPeriods,
										const ARM_PricingStatesPtr& states) const;
	
	/// Standard ARM object support
	virtual string					toString(const string& indent="",const string& nextIndent="") const;

	void			SetNbEffectiveReset(int nb);

	int				GetNbEffectiveReset() const;

	const ARM_GP_Vector&	GetResetTimes() const {return itsResetTimes;};
	const ARM_GP_Vector&	GetStartTimes() const {return itsStartTimes;};
	const ARM_GP_Vector&	GetEndTimes1() const {return itsEndTimes1;};
	const ARM_GP_Vector&	GetEndTimes2() const {return itsEndTimes2;};
	const ARM_VectorVector&	GetEigenValues() const {return itsEigenValues;};

	virtual size_t	ModelStatesSize() const { return itsNbEffectiveReset + 1;}	
	
	void			CalcNbEffectiveReset(const ARM_GP_Vector& timeSteps);

	ARM_GP_Vector	GetATMGaussianVol();

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
