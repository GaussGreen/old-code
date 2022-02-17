/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
  */



#ifndef _INGPMODELS_SMFRM_H
#define _INGPMODELS_SMFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/vectormanip.h"

#include "gpinfra/pricingmodelir.h"
#include "gpcalib/numerical.h"

#include "gpmodels/ProcessBuilderSmiledFRM.h"
#include "gpmodels/ModelParamsSmiledFRM.h"
#include "gpmodels/argconvdefault.h"



CC_BEGIN_NAMESPACE( ARM )

class ARM_VanillaSwaptionArgSmiledFRM;


class ARM_SmiledFRM : public ARM_PricingModelIR
{
private:

	inline bool		DoesResetDateExist( double date ) const				{ return ExistsInVector( itsResetDates, date );	}

	inline double	getTerminalTime() const								{ return itsEndDates.Elt( itsEndDates.size() -1 ); }

	inline double	getLastResetTime() const							{ return itsResetDates.Elt( itsResetDates.size() -1 ); }

	inline void		setResetTimes( const ARM_GP_Vector& resetTimes )	{ itsResetDates = resetTimes;	}

	inline void		setEndTimes( const ARM_GP_Vector& endTimes )		{ itsEndDates   = endTimes;		}

	inline void		setStartTimes( const ARM_GP_Vector& startTimes )	{ itsStartDates = startTimes;	}


	void			computeWeights();

	void			ComputeNumeraireTimeIndex();

	double			GetEquivalentShift(const ARM_GP_Vector& weight) const;

	ARM_VectorPtr	GetDiffusionCoeff() const;

	
public:
	ARM_SmiledFRM( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL,
		size_t timeStepsNb=500,
		size_t gridSize=501,
		double stdDevNb=6,
		bool skipPDE=false, 
		bool allowInterpol=false, 
		ARM_ModelParamsSmiled::CalibProxy swaptionApprox=ARM_ModelParamsSmiled::LocalVolatilityWithRescaling,
		bool withRescaling=false);

	ARM_SmiledFRM( const ARM_SmiledFRM& rhs);
	virtual ~ARM_SmiledFRM();
	ASSIGN_OPERATOR(ARM_SmiledFRM)

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

	virtual void	setNumericalModelFitter(ARM_NumericalModelFitter*);
	virtual void	MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	
	virtual void	PostInit();
	
	virtual void	setCalibrationStatus( bool calibrationStatus )		{ itsCalibrationStatus = calibrationStatus;};

	virtual void	setCalibratedStatus( bool calibratedStatus )		{ itsCalibrated = calibratedStatus;};

	virtual ARM_PricingStatesPtr	FirstPricingStates( size_t bucketSize ) const;

	virtual size_t ModelStatesSize() const { return itsResetDates.size();}

	virtual ARM_PricingStatesPtr	Init(	const string& payModelName, 
											const ARM_TimeInfoPtrVector& timeInfos);

	virtual ARM_PricingStatesPtr	Induct(ARM_PricingStatesPtr& states,double toTime);

	virtual ARM_GP_Vector*          PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos);

	virtual ARM_GP_Vector*			ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos );
	
	virtual ARM_VectorPtr			ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;

	virtual void					AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb = 0 );

	virtual ARM_VectorPtr	DiscountFactor(	const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr	DiscountFactorNoInterpol(	const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr	ForwardDiscountFactor(const string& curveName,
											double evalTime, 
											double startTime,
											double endTime,
											const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr	Libor(	const string& curveName,
									double evalTime,
									double fwdStartTime,
									double fwdEndTime,
									double period,
									double resetTime,
									double payTime,
									const ARM_PricingStatesPtr& states) const;


	bool			FastCompute() const						{ return true;					};

	bool			IsOnSamePath(size_t i, size_t j) const	{ return ( itsWeight(i,j)==1 ); };


	ARM_VectorPtr	ForwardDiscountFactorFromIdx(	const string& curveName,
											double evalTime,
											size_t IdxFrom,
											size_t IdxTo,
											size_t modelNb,
											const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr	ComputeDSwapDLib(	const ARM_GP_Vector& startDates,const ARM_GP_Vector& endDates,
										const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayPeriods) const;

	ARM_VectorPtr	ComputeDLogShiftSwapDLib(	double shift, 
										const ARM_GP_Vector& startDates,const ARM_GP_Vector& endDates,
										const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayPeriods) const;

	double			DZCDLib( double maturity, size_t i ) const;
	
////////////////////////////////////////////////////////////// MC METHOD

	virtual void	ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps);

	void			ModelStateLocalStdDev(const ARM_GP_Vector& timeSteps, const ARM_MatrixVector& localVariances);

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
	virtual ARM_Object*		Clone() const { return new ARM_SmiledFRM(*this); };
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


	const ARM_GP_Vector&					GetResetDates() const {return itsResetDates;};
	const ARM_GP_Vector&					GetEndDates() const {return itsEndDates;};
	const ARM_VectorVector&					GetEigenValues() const {return itsEigenValues;};
	const size_t&							GetNumeraireTimeIndex() const {return itsNumeraireTimeIndex;};
	bool									GetAllowInterpol() const {return itsAllowInterpol;};

private:

	ARM_NumericalModelFitter*				itsNumericalModelFitter;	//assoc

	ARM_GP_Vector							itsResetDates;
	ARM_GP_Vector							itsStartDates;
	ARM_GP_Vector							itsEndDates;

	vector< ARM_ProcessBuilder* >			itsProcess;

	CC_IS_MUTABLE ARM_VectorVector			itsEigenValues;
	size_t									itsNumeraireTimeIndex;

	// BIDOUILLES !!
	bool									itsCalibrationStatus;
	bool									itsCalibrated;
	/////////////////////////////////////////// FIN BIDOUILLES
	bool									itsSkipPDE;
	bool									itsAllowInterpol;

	size_t									itsTimeStepsNb;
	size_t									itsGridSize;
	size_t									itsStdDevNb;
	
	/// --------------- fast calibration caching arguments -----------------
    CC_IS_MUTABLE ARM_VanillaSwaptionArgSmiledFRM*			itsCurrentArg;
	CC_IS_MUTABLE bool										itsNewArgNeeded;

	ARM_ModelParamsSmiled::CalibProxy						itsSwaptionApprox;

	CC_IS_MUTABLE ARM_GP_VectorPtr							itsCachedTimeSteps;
	CC_IS_MUTABLE ARM_MatrixVector							itsCachedModelStatesLocalVars;

	ARM_GP_Matrix											itsWeight;
	bool													itsWithRescalling;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/

