/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
  */



#ifndef _INGPMODELS_SMFRM_H
#define _INGPMODELS_SMFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

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

	inline void		setResetTimes( const std::vector<double>& resetTimes )	{ itsResetDates = resetTimes;	}

	inline void		setEndTimes( const std::vector<double>& endTimes )		{ itsEndDates   = endTimes;		}

	inline void		setStartTimes( const std::vector<double>& startTimes )	{ itsStartDates = startTimes;	}


	void			computeWeights();

	void			ComputeNumeraireTimeIndex();

	double			GetEquivalentShift(const std::vector<double>& weight) const;

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

	virtual std::vector<double>*          PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos);

	virtual std::vector<double>&			ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos );
	
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

	ARM_VectorPtr	ComputeDSwapDLib(	const std::vector<double>& startDates,const std::vector<double>& endDates,
										const std::vector<double>& fixPayTimes,const std::vector<double>& fixPayPeriods) const;

	ARM_VectorPtr	ComputeDLogShiftSwapDLib(	double shift, 
										const std::vector<double>& startDates,const std::vector<double>& endDates,
										const std::vector<double>& fixPayTimes,const std::vector<double>& fixPayPeriods) const;

	double			DZCDLib( double maturity, size_t i ) const;
	
////////////////////////////////////////////////////////////// MC METHOD

	virtual void	ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps);

	void			ModelStateLocalStdDev(const std::vector<double>& timeSteps, const ARM_MatrixVector& localVariances);

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


	const std::vector<double>&					GetResetDates() const {return itsResetDates;};
	const std::vector<double>&					GetEndDates() const {return itsEndDates;};
	const ARM_VectorVector&					GetEigenValues() const {return itsEigenValues;};
	const size_t&							GetNumeraireTimeIndex() const {return itsNumeraireTimeIndex;};
	bool									GetAllowInterpol() const {return itsAllowInterpol;};

private:

	ARM_NumericalModelFitter*				itsNumericalModelFitter;	//assoc

	std::vector<double>							itsResetDates;
	std::vector<double>							itsStartDates;
	std::vector<double>							itsEndDates;

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

