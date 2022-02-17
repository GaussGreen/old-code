#ifndef _INGPMODELS_MODELBGMSV2F_H
#define _INGPMODELS_MODELBGMSV2F_H

#include "gpmodels/modelparamsbgmsv2f.h"

#include "gpinfra/pricingmodelir.h"
#include "gpnumlib/ran2.h"


CC_BEGIN_NAMESPACE( ARM )


class ARM_BGMSV2F : public ARM_PricingModelIR 
{
protected:
	
	bool							itsCalibrationStatus;
	bool							itsCalibratedStatus;
	
	// les échéanciers
	ARM_GP_Vector					itsResetTimes;
	ARM_GP_Vector					itsStartTimes;
	ARM_GP_Vector					itsEndTimes;
	ARM_GP_Vector					itsDelta;
	ARM_GP_Vector					itsFwdRate;
	ARM_GP_Vector					itsFromShiftedToRate;
	ARM_GP_Vector					itsATMVols;

	bool							itsAllowInterpol;	// interpolation des zc
	bool							itsInterpolZC;

	ARM_GP_Matrix					itsWeight;

	bool							itsComputeDriftAdj;
	
	ARM_VectorVector				itsEigenValues;
	ARM_MatrixVector				itsModelLocalRealVar;

	ARM_VanillaSecDensityPtrVector	itsCalibSecDensities;
	ARM_IntVector					itsCalibWeight;
	bool							itsSpreadCalib;

	bool							itsProxyStatus;
	int								itsNbEffectiveReset;
	
	mutable ARM_RandUniform_NRRan2	itsRandGen;
	mutable ARM_GP_Vector			itsU;
	mutable	double					itsPrevEvalTime;

public:

	ARM_BGMSV2F( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params = NULL, 
			bool AllowInterpol = true, bool ComputeDrift = true, bool Proxy = false);
	
	ARM_BGMSV2F( const ARM_BGMSV2F& rhs);
	
	virtual ~ARM_BGMSV2F();

	virtual ARM_Object*				Clone() const	{ return new ARM_BGMSV2F(*this);};
	
	ARM_BGMSV2F&					operator = (const ARM_BGMSV2F& rhs);

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

	virtual ARM_VectorPtr			Libor(const string& curveName,
										double evalTime,
										double fwdStartTime,
										double fwdEndTime,
										double period,
										double resetTime,
										double payTime,
										const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr			SwapRate(const string& curveName, 
										double evalTime,
										double floatStartTime,
										double floatEndTime, 
										const ARM_GP_Vector& fixPayTimes,
										const ARM_GP_Vector& fixPayPeriods,
										const ARM_GP_Vector& fwdStartTimes,
										const ARM_GP_Vector& fwdEndTimes,
										const ARM_GP_Vector& fwdPayPeriods, 
										const ARM_GP_Vector& floatPayTimes,
										const ARM_GP_Vector& floatPayPeriods,
										const ARM_GP_Vector& margin,
										bool isDbleNotional,
										const ARM_PricingStatesPtr& states) const;

	/*
	virtual ARM_VectorPtr			MaxRate(const string& curveName, 
										double evalTime,
										double floatStartTime,
										double floatEndTime, 
										const ARM_GP_Vector& fixPayTimes,
										const ARM_GP_Vector& fixPayPeriods,
										const ARM_GP_Vector& fwdStartTimes,
										const ARM_GP_Vector& fwdEndTimes,
										const ARM_GP_Vector& fwdPayPeriods, 
										const ARM_GP_Vector& floatPayTimes,
										const ARM_GP_Vector& floatPayPeriods,
										const ARM_GP_Vector& margin,
										double firstReset,
										double firstStart,
										const ARM_VectorPtr& firstRate,
										int MaxOrMin,
										int ResetFreq,
										const ARM_VectorPtr& strikes,
										int CapOrFloor,
										double RhoMinMax,
										bool IsAccrued,
										double MinAccrued,
										double MaxAccrued,
										const ARM_PricingStatesPtr& states) const;
	*/

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

	/*
	virtual ARM_VectorPtr			ImpliedVol(const string& curveName,
										double evalTime,
										double payTime,
										double period,
										double payNotional,
										double fwdResetTime,	/// used for volatility computation
										double fwdStartTime,
										double fwdEndTime,
										double fwdPeriod,
										const ARM_GP_Vector& strikesPerState,
										int capFloor,
										const ARM_PricingStatesPtr& states) const;
	*/

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

	bool			GetProxyStatus() const {return itsProxyStatus;};
	void			SetProxyStatus(bool status);
	void			SetNbEffectiveReset(int nb);
	int				GetNbEffectiveReset() const {return itsNbEffectiveReset;};
	virtual size_t	ModelStatesSize() const { return itsNbEffectiveReset + 2;}	

	const ARM_VectorVector&	GetEigenValues() const {return itsEigenValues;};
	
	void			CalcNbEffectiveReset(const ARM_GP_Vector& timeSteps);

	const ARM_GP_Vector&	GetResetTimes() const {return itsResetTimes;};
	const ARM_GP_Vector&	GetStartTimes() const {return itsStartTimes;};
	const ARM_GP_Vector&	GetEndTimes() const {return itsEndTimes;};

protected:
	
	double			getTerminalTime() const;

	bool			DoesResetTimeExist(double time) const;

	void			computeWeights();

	bool			IsOnSamePath(size_t i, size_t j) const;

	ARM_VectorPtr	DiscountFactorNoInterpol(const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr	ForwardDiscountFactor(const string& curveName,
											double evalTime, 
											double startTime,
											double endTime,
											const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr	ForwardDiscountFactorFromIdx(const string& curveName,
											double evalTime,
											size_t IdxFrom,
											size_t IdxTo,
											size_t modelNb,
											const ARM_PricingStatesPtr& states) const;

	// pour le pricing d'une swaption à partir des libors
	void			GetWeightSwaptionApprox(const string& curveName, double evalTime, const ARM_GP_Vector& floatNotional, const ARM_GP_Vector& floatStartTimes, const ARM_GP_Vector& floatEndTimes, const ARM_GP_Vector& floatPayPeriods, 
											int idx, double fwdrate, double annuity, ARM_GP_Vector& weigth) const;

	double			GetEquivSwaptionShift(int idx, const ARM_IntVector& idxs, ARM_GP_Vector& weight) const;
	double			GetEquivSwaptionLevel(int idx, const ARM_IntVector& idxs, const ARM_GP_Vector& weight, const ARM_GP_Vector * optweight = NULL) const;
	double			GetEquivSwaptionRho(int idx, const ARM_IntVector& idxs, const ARM_GP_Vector& weight, double swlevel, int UnOuDeux) const;
};

inline void ARM_BGMSV2F::PostInit()
{
}

inline void ARM_BGMSV2F::setCalibrationStatus(bool calibrationStatus)
{
	itsCalibrationStatus = calibrationStatus;
}

inline void ARM_BGMSV2F::setCalibratedStatus(bool calibratedStatus)
{
	itsCalibratedStatus = calibratedStatus;
}

inline void ARM_BGMSV2F::SetProxyStatus(bool status)
{
	itsProxyStatus = status;
	if(status == true) itsComputeDriftAdj = false;
}

inline void ARM_BGMSV2F::SetNbEffectiveReset(int nb)
{
	itsNbEffectiveReset = nb;
}

inline void ARM_BGMSV2F::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_BGMSV2F::TreeStatesToModelStates :  not implemented ARM_BGMSV2F Model!");
}

inline double ARM_BGMSV2F::VarianceToTime(double var,double minTime,double maxTime) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_BGMSV2F::VarianceToTime not implemented ARM_BGMSV2F Model!");

}

inline double ARM_BGMSV2F::getTerminalTime() const
{
	return itsEndTimes.Elt( itsEndTimes.size() -1 );
}

inline bool ARM_BGMSV2F::DoesResetTimeExist(double time) const
{
	return ExistsInVector( itsResetTimes, time );
}

inline bool ARM_BGMSV2F::IsOnSamePath(size_t i, size_t j) const
{
	return itsWeight(i,j) == 1.;
}

CC_END_NAMESPACE( )

#endif