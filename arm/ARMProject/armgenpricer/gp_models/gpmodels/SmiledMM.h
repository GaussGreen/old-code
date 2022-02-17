/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
  */



#ifndef _INGPMODELS_SMM_H
#define _INGPMODELS_SMM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"

#include "gpinfra/pricingmodelir.h"
#include "gpcalib/numerical.h"

#include "gpmodels/ProcessBuilderSmiledFRM.h"
#include "gpmodels/ModelParamsSmiledFRM.h"
#include "gpmodels/argconvdefault.h"



CC_BEGIN_NAMESPACE( ARM )

class ARM_VanillaSwaptionArgSmiledFRM;
class ARM_VanillaSpreadOptionArgSmiledFRM;


class ARM_SmiledMM : public ARM_PricingModelIR
{

protected:

	inline bool		DoesResetDateExist( double date ) const				{ return ExistsInVector( itsResetDates, date );	}

	inline void		setResetTimes( const std::vector<double>& resetTimes )	{ itsResetDates = resetTimes;	}

	inline void		setEndTimes( const std::vector<double>& endTimes )		{ itsEndDates   = endTimes;		}

	inline void		setStartTimes( const std::vector<double>& startTimes )	{ itsStartDates = startTimes;	}

	inline double	getTerminalTime() const								{ return itsEndDates.Elt( itsEndDates.size() -1 ); }

	void			ComputeNumeraireTimeIndex();

	double			GetEquivalentShift(const std::vector<double>& weight) const;

	ARM_VectorPtr	GetDiffusionCoeff() const;

	ARM_VectorPtr	ComputeDSwapDRate(	const std::vector<double>& startDates,const std::vector<double>& endDates,
										const std::vector<double>& fixPayTimes,const std::vector<double>& fixPayPeriods) const;

	ARM_VectorPtr	ComputeDLogShiftSwapDRate(	double shift, 
										const std::vector<double>& startDates,const std::vector<double>& endDates,
										const std::vector<double>& fixPayTimes,const std::vector<double>& fixPayPeriods) const;

public:
	ARM_SmiledMM( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL, size_t timeStepsNb=500,size_t gridSize=501,double stdDevNb=6,bool skipPDE=false, bool allowInterpol=false, ARM_ModelParamsSmiled::CalibProxy calibProxy=ARM_ModelParamsSmiled::LocalVolatility,bool cache=false);
	ARM_SmiledMM( const ARM_SmiledMM& rhs);
	virtual ~ARM_SmiledMM();

	virtual bool	ValidateModelParams(const ARM_ModelParams& params) const;
	virtual void	ValidateCalibMethod(ARM_CalibMethod& calibMethod);
	virtual void	SetNumeraire(const ARM_NumerairePtr& numerairePtr);

	
    virtual void	PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void	PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void	AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void	AdviseCurrentCalib(ARM_ModelFitter& modelFitter){};
	virtual void	Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
	
    
	virtual bool	SupportBackwardInduction() const {	return false;}
    virtual bool	SupportForwardInduction()  const {	return true;}
	virtual bool	SupportAnalyticMarginal()  const {	return false;}

	virtual void	PostInit();
	
	virtual void	setCalibrationStatus( bool calibrationStatus )		{ itsCalibrationStatus = calibrationStatus;};

	virtual void	setCalibratedStatus( bool calibratedStatus )		{ itsCalibrated = calibratedStatus;};

	virtual ARM_PricingStatesPtr	FirstPricingStates( size_t bucketSize ) const;

	virtual ARM_PricingStatesPtr	Init(	const string& payModelName, 
											const ARM_TimeInfoPtrVector& timeInfos);

	virtual ARM_PricingStatesPtr	Induct(ARM_PricingStatesPtr& states,double toTime);

	virtual std::vector<double>&			ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos );
	
	virtual ARM_VectorPtr			ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;

	virtual void					AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb = 0 );


////////////////////////////////////////////////////////////// MC METHOD

	virtual void	ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps);
	virtual void	NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps );

	virtual void	ModelStateLocalVariances(	const std::vector<double>& timeSteps,
											ARM_MatrixVector& localVariances ) const;			

	virtual void	NumMethodStateLocalVariances(	const std::vector<double>& timeSteps,
												ARM_MatrixVector& localVariances ) const;		

	virtual void	NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
												ARM_MatrixVector& globalVariances ) const;		

	virtual bool	NeedsToCholeskyDecomposeFactors( ) const { return false; }					

////////////////////////////////////////////////////////////// END DISCRETIZATION FOR MC METHOD


    
/// Standard ARM object support
	virtual string			ExportShortName() const { return "LSBGM";}
	virtual	string			toString(const string& indent="",const string& nextIndent="") const;

/////////// NOT AVAILABLE
	virtual void	TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

	virtual double	VarianceToTime(	double var,
									double minTime=0.0,
									double maxTime=5*K_YEAR_LEN) const;

	virtual ARM_VectorPtr	VanillaCaplet(	const string& curveName, 
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

	virtual ARM_VectorPtr	VanillaSwaption(	const string& curveName,
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

	ARM_VanillaSwaptionArgSmiledFRM* GetVanillaSwaptionArg( 
		const string& curveName,
		double evalTime,
		double swapResetTime,
		double swapNotional,
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
		bool isConstantNotional = true		) const;

	ARM_VanillaSpreadOptionArgSmiledFRM* GetVanillaSpreadOptionArg( 
		const string& curveName,
		double evalTime,
		int	callPut,
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
		const std::vector<double>& swapShortFixPayPeriods) const;

	virtual double UnderlyingCorrelation(  string	underlyingType,
										   double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   double	startTime3,
										   double   endTime3,
										   double	startTime4,
										   double   endTime4) const;

	virtual double UnderlyingCovariance (  string	underlyingType,
										   double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   double	startTime3,
										   double   endTime3,
										   double	startTime4,
										   double   endTime4) const;

	virtual void ComputeUnderlyingVarCovar( double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   std::vector<double>& VarCovar,
										   std::vector<double>& swapfwd) const;

    virtual ARM_VectorPtr	VanillaDigital(	const string& curveName, 
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

	virtual ARM_VectorPtr	VanillaSpreadOptionLet(const string& curveName,
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

	virtual void	IntegratedLocalDrifts(	const	std::vector<double>& timeSteps,
													ARM_GP_MatrixPtr& relativeDrifts,
													ARM_GP_MatrixPtr& absoluteDrifts) const;

///////////////////////////////////////////////////////////////////////////////END NOT AVAILABLE

	virtual void	computeDWeights(){};
	
	virtual double	DZCDRate( double maturity, size_t i ) const;

	virtual void	MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	
	virtual ARM_VectorPtr	DiscountFactor(	const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr	Libor(	const string& curveName,
									double evalTime,
									double fwdStartTime,
									double fwdEndTime,
									double period,
									double resetTime,
									double payTime,
									const ARM_PricingStatesPtr& states) const;


protected:

	ARM_NumericalModelFitter*				itsNumericalModelFitter;
	ARM_VanillaSecDensityPtrVector			itsCalibSecDensities;

	std::vector<double>							itsResetDates;

	std::vector<double>							itsStartDates;
	std::vector<double>							itsEndDates;

	vector< ARM_ProcessBuilder* >			itsProcess;

	ARM_VectorVector						itsEigenValues;
	size_t									itsNumeraireTimeIndex;

	bool									itsCalibrationStatus;
	bool									itsCalibrated;
	bool									itsSkipPDE;
	bool									itsAllowInterpol;

	size_t									itsTimeStepsNb;
	size_t									itsGridSize;
	size_t									itsStdDevNb;
	
	ARM_ModelParamsSmiled::CalibProxy	itsCalibProxy;

	ARM_GP_VectorPtr						itsCachedTimeSteps;
	ARM_MatrixVector						itsCachedModelStatesLocalVars;

	/// --------------- fast calibration caching arguments -----------------
    CC_IS_MUTABLE ARM_VanillaSwaptionArgSmiledFRM*			itsCurrentArg;
	CC_IS_MUTABLE bool										itsNewArgNeeded;

	bool									itsCache;

	CC_IS_MUTABLE ARM_VectorPtr							itsRateWeight;

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/

