/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
  */



#ifndef _INGPMODELS_SMM_H
#define _INGPMODELS_SMM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

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

	inline void		setResetTimes( const ARM_GP_Vector& resetTimes )	{ itsResetDates = resetTimes;	}

	inline void		setEndTimes( const ARM_GP_Vector& endTimes )		{ itsEndDates   = endTimes;		}

	inline void		setStartTimes( const ARM_GP_Vector& startTimes )	{ itsStartDates = startTimes;	}

	inline double	getTerminalTime() const								{ return itsEndDates.Elt( itsEndDates.size() -1 ); }

	void			ComputeNumeraireTimeIndex();

	double			GetEquivalentShift(const ARM_GP_Vector& weight) const;

	ARM_VectorPtr	GetDiffusionCoeff() const;

	ARM_VectorPtr	ComputeDSwapDRate(	const ARM_GP_Vector& startDates,const ARM_GP_Vector& endDates,
										const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayPeriods) const;

	ARM_VectorPtr	ComputeDLogShiftSwapDRate(	double shift, 
										const ARM_GP_Vector& startDates,const ARM_GP_Vector& endDates,
										const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayPeriods) const;

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

	virtual ARM_GP_Vector*			ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos );
	
	virtual ARM_VectorPtr			ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;

	virtual void					AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb = 0 );


////////////////////////////////////////////////////////////// MC METHOD

	virtual void	ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps);
	virtual void	NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps );

	virtual void	ModelStateLocalVariances(	const ARM_GP_Vector& timeSteps,
											ARM_MatrixVector& localVariances ) const;			

	virtual void	NumMethodStateLocalVariances(	const ARM_GP_Vector& timeSteps,
												ARM_MatrixVector& localVariances ) const;		

	virtual void	NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
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
											const ARM_GP_Vector& strikesPerState,
											int capFloor,
											const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr	VanillaSwaption(	const string& curveName,
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

	ARM_VanillaSwaptionArgSmiledFRM* GetVanillaSwaptionArg( 
		const string& curveName,
		double evalTime,
		double swapResetTime,
		double swapNotional,
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
		const ARM_GP_Vector& strikes,
		double swapLongFloatStartTime,
		double swapLongFloatEndTime,
		const ARM_GP_Vector& swapLongFixPayTimes,
		const ARM_GP_Vector& swapLongFixPayPeriods,
		double swapShortFloatStartTime,
		double swapShortFloatEndTime,
		const ARM_GP_Vector& swapShortFixPayTimes,
		const ARM_GP_Vector& swapShortFixPayPeriods) const;

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
										   ARM_GP_Vector& VarCovar,
										   ARM_GP_Vector& swapfwd) const;

    virtual ARM_VectorPtr	VanillaDigital(	const string& curveName, 
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

	virtual void	IntegratedLocalDrifts(	const	ARM_GP_Vector& timeSteps,
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

	ARM_GP_Vector							itsResetDates;

	ARM_GP_Vector							itsStartDates;
	ARM_GP_Vector							itsEndDates;

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

