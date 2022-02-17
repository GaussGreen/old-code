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
	std::vector<double>					itsResetTimes;
	std::vector<double>					itsStartTimes;
	std::vector<double>					itsEndTimes;
	std::vector<double>					itsDelta;
	std::vector<double>					itsFwdRate;
	std::vector<double>					itsFromShiftedToRate;
	std::vector<double>					itsATMVols;

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
	mutable std::vector<double>			itsU;
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
										const std::vector<double>& fixPayTimes,
										const std::vector<double>& fixPayPeriods,
										const std::vector<double>& fwdStartTimes,
										const std::vector<double>& fwdEndTimes,
										const std::vector<double>& fwdPayPeriods, 
										const std::vector<double>& floatPayTimes,
										const std::vector<double>& floatPayPeriods,
										const std::vector<double>& margin,
										bool isDbleNotional,
										const ARM_PricingStatesPtr& states) const;

	/*
	virtual ARM_VectorPtr			MaxRate(const string& curveName, 
										double evalTime,
										double floatStartTime,
										double floatEndTime, 
										const std::vector<double>& fixPayTimes,
										const std::vector<double>& fixPayPeriods,
										const std::vector<double>& fwdStartTimes,
										const std::vector<double>& fwdEndTimes,
										const std::vector<double>& fwdPayPeriods, 
										const std::vector<double>& floatPayTimes,
										const std::vector<double>& floatPayPeriods,
										const std::vector<double>& margin,
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
										const std::vector<double>& strikesPerState,
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

	bool			GetProxyStatus() const {return itsProxyStatus;};
	void			SetProxyStatus(bool status);
	void			SetNbEffectiveReset(int nb);
	int				GetNbEffectiveReset() const {return itsNbEffectiveReset;};
	virtual size_t	ModelStatesSize() const { return itsNbEffectiveReset + 2;}	

	const ARM_VectorVector&	GetEigenValues() const {return itsEigenValues;};
	
	void			CalcNbEffectiveReset(const std::vector<double>& timeSteps);

	const std::vector<double>&	GetResetTimes() const {return itsResetTimes;};
	const std::vector<double>&	GetStartTimes() const {return itsStartTimes;};
	const std::vector<double>&	GetEndTimes() const {return itsEndTimes;};

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
	void			GetWeightSwaptionApprox(const string& curveName, double evalTime, const std::vector<double>& floatNotional, const std::vector<double>& floatStartTimes, const std::vector<double>& floatEndTimes, const std::vector<double>& floatPayPeriods, 
											int idx, double fwdrate, double annuity, std::vector<double>& weigth) const;

	double			GetEquivSwaptionShift(int idx, const ARM_IntVector& idxs, std::vector<double>& weight) const;
	double			GetEquivSwaptionLevel(int idx, const ARM_IntVector& idxs, const std::vector<double>& weight, const std::vector<double> * optweight = NULL) const;
	double			GetEquivSwaptionRho(int idx, const ARM_IntVector& idxs, const std::vector<double>& weight, double swlevel, int UnOuDeux) const;
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