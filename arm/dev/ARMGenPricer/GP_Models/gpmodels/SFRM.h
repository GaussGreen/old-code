/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SFRM.h
 *
 *  \brief SFRM model multi-factor version!
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPMODELS_SFRM_H
#define _INGPMODELS_SFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpinfra/pricingmodelir.h"
#include "gpcalib/typedef.h"
#include "gpbase/gpmatrix.h" /// for accessor
#include <deque>
CC_USING_NS(std,deque)
/// STL
#include <map>
CC_USING_NS(std,map)

/// forward declaration
class ARM_IRIndex;
class ARM_Swaption;
class ARM_CapFloor;
class ARM_Digital;
class ARM_CorridorLeg;
class ARM_Security;
class ARM_Swaption;
class ARM_SpreadOption;


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ModelParamsSFRM;
struct ARM_VanillaSwaptionArgSFRM;
struct ARM_VanillaCapArgSFRM;
struct ARM_VanillaDigitalArgSFRM;
struct ARM_VanillaCorridorLegArgSFRM;
struct ARM_VanillaSpreadOptionletArgSFRM;
class ARM_DateStrip;
class GaussLegendre_Coefficients;


///#define __GP_SFRM_STOCH_TERMS 

////////////////////////////////////////////////////////////
/// \class ARM_SFRM
/// \brief
///  Interface class for Hull & White pricing models
////////////////////////////////////////////////////////////

class ARM_SFRM : public ARM_PricingModelIR
{
public :
	
	/// map to save DF vectors in SFRM at each maturity
	/// The map is resetted when we change the pricing states
	typedef map<double, ARM_VectorPtr> ARM_DFMap;	
    typedef ARM_DFMap::iterator ARM_DFMapIter;
	/// Mutable because GenerateFunction() will update this map
    CC_IS_MUTABLE ARM_DFMap itsDFMap;

    /// Precomputed data
	ARM_VectorPtr itsFwdResetTimes;
	ARM_VectorPtr itsFwdStartTimes;
	ARM_VectorPtr itsFwdEndTimes;
	ARM_VectorPtr itsFwdValues;
	ARM_VectorPtr itsFwdIntTerms;
	ARM_VectorPtr itsFwdShifts;
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ itsFwdStatus;
	bool itsIsFixStartDate;
	ARM_Date itsFixStartDate;
	bool itsIsFixEndDate;
	ARM_Date itsFixEndDate;


	ARM_TriangularMatrixVector itsFwdsCorrelVarCoeffs;
	ARM_VectorPtr itsFwdDrifts;
	ARM_VectorPtr itsFwdsVarCoeffs;

	/// be used for faster computation
	/// for fast computation to avoid creation and deletion of objects!
	ARM_VectorPtr itsFactors;

	/// use for location of the min fwd index!
	CC_IS_MUTABLE int itsCurrentFwdMaxIndex;
	CC_IS_MUTABLE int itsNumeraireTimeIndex;

	ARM_DateStrip* GetFloatDateStrip( const ARM_Date& startDate, const ARM_Date& endDate ) const;
	friend class ARM_ModelParamsSFRM;

	/// --------------- fast calibration caching arguments -----------------
    CC_IS_MUTABLE ARM_VanillaArg* itsCurrentArg;
	bool itsUpdateVolSwapvolFRA;

	mutable ARM_CountedPtr<GaussLegendre_Coefficients> itsLegendreCoeffs;

#ifdef __GP_SFRM_STOCH_TERMS
    // Log of the stochastic part of the forward value, used in the interpolation of the discount
    // factor
	CC_IS_MUTABLE ARM_GP_Matrix*		itsStochasticTerms;
	CC_IS_MUTABLE ARM_GP_Matrix						itsProbaChangesPoly;
	CC_IS_MUTABLE ARM_GP_Vector						itsProbaChangesVar;
#endif


	

	// Member function pointer to switch the method TreeStatesToModelStates
	typedef void (ARM_SFRM::*TreeStatesToModelStatesFunc)(ARM_PricingStatesPtr&,int) const;

	// This flag allow the non param drift;
	TreeStatesToModelStatesFunc itsTreeStatesToModelStatesFunc;

	
	ARM_VanillaSwaptionArgSFRM* ComputeCachedSwaptionData       (    ARM_Swaption* swaption         );
	ARM_VanillaCapArgSFRM*		ComputeCachedCapData			(    ARM_CapFloor* security			);
	ARM_VanillaDigitalArgSFRM*	ComputeCachedDigitalData		(    ARM_Digital* security			);
    ARM_VanillaCorridorLegArgSFRM* ComputeCachedCorridorLegData (    ARM_CorridorLeg* corridorLeg   );
	ARM_VanillaSpreadOptionletArgSFRM* ComputeCachedSpreadOptionletData(    ARM_SpreadOption* spreadoption );

	ARM_VanillaSwaptionArgSFRM* GetVanillaSwaptionArg( 
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

	ARM_VanillaSpreadOptionletArgSFRM* GetVanillaSpreadOptionletArg( 
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

	static const bool ARM_SFRM_STD_FWD;
	static const bool ARM_SFRM_AUX_FWD;

	//// function to Reset the Map
	void ResetDFMap () const;


public:
	ARM_SFRM( const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsSFRM& params, ARM_Portfolio* shiftConvPort = NULL, bool nonParamDrift = false);
	ARM_SFRM(const ARM_SFRM& rhs);
	virtual ~ARM_SFRM();
    ARM_SFRM& operator = (const ARM_SFRM& rhs);

    virtual void SetZeroCurve(const ARM_ZeroCurvePtr& zc);

    /// accessors
    inline ARM_VanillaArg* GetVanillaArg() const { return itsCurrentArg; }
    inline ARM_VectorPtr GetFwdValues() const { return itsFwdValues; }
    inline void SetFwdValues(ARM_VectorPtr& FwdValues) {itsFwdValues = FwdValues; }

	

    /// Closed form formulas
	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const;

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
        const ARM_GP_Vector& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const;

    virtual ARM_VectorPtr VanillaDigital(
		const string& curveName, 
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
    
    virtual ARM_VectorPtr VanillaCorridorlet(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     indexPaymentType,
        double  fwdPaymentPeriod,
        const ARM_GP_Vector& RefIdxResettimes,
        const ARM_GP_Vector& RefIdxStarttimes,
        const ARM_GP_Vector& RefIdxEndtimes,
        const ARM_GP_Vector& RefFwdPeriods,
        const ARM_GP_Vector& RefIndexWeight,
        double  couponMargin,
        const vector<const ARM_GP_Vector*> downBarrierPerState,
        const vector<const ARM_GP_Vector*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const;   
	
	/// Pricing Function for Restrikeable
	virtual ARM_VectorPtr VanillaCorridorletRestrikeable(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     fwdPaymentType, 
        double  fwdPaymentPeriod,
        const ARM_GP_Vector& RefIdxResettimes,
        const ARM_GP_Vector& RefIdxStarttimes,
        const ARM_GP_Vector& RefIdxEndtimes,
        const ARM_GP_Vector& RefFwdPeriods,
        const ARM_GP_Vector& RefIndexWeight,
        double  couponMargin,
        const vector<const ARM_GP_Vector*> downBarrierPerState,
        const vector<const ARM_GP_Vector*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const ;

	/// Variable Notional swaptions are priced in SFRM
	virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return isConstantSpread&&isConstantstrike;}

	virtual ARM_VectorPtr VanillaSwaption(
		const string& curveName,
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

    /// Default initialisation of the model
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
	virtual ARM_PricingStatesPtr ReInit();
    virtual bool SupportBackwardInduction() const {	return true; }
    virtual bool SupportForwardInduction()  const { return true; }
	virtual bool SupportAnalyticMarginal()  const {	return false;}

    virtual bool NeedArrowDebreuPrices() const { return true; }
    virtual bool NeedLocalDiscount() const { return false; }
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	/// non const to allow caching!
	virtual ARM_GP_Vector* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos );
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;
	//virtual ARM_GP_Vector*          PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos);

    /// Give local drifts and variances w.r.t. a given schedule
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	void SigneItsModelStateLocalStdDevs() const;

	// Compute the numeraire time index based on the numeraire
	// current date.
	void ComputeNumeraireTimeIndex() const;

	/// function to compute variances with full matrices
	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ComputeFwdsCorrelVarCoeffs(void) const;
	virtual void ComputeFwdDrifts(void);

	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

	void ARM_SFRM::CreateProbaChange(
		ARM_PricingStatesPtr& states,
		const ARM_GP_VectorPtr& proba,
		int fwdIdx,
		double var) const;

	double IntegrateProbaChange(
		int n,
		double stdDev,
		const ARM_GP_Vector& states,
		const ARM_GP_Vector& probaChanges,
		bool fromMinMax,
		double MinMax,
		double x0) const;

	void TreeStatesToModelStatesDriftMarkovian(ARM_PricingStatesPtr& states, int timeIndex) const;
	void TreeStatesToModelStatesNonParamDrift(ARM_PricingStatesPtr& states, int timeIndex) const;

	// Do we need to add the model time in discretisation scheme of the numerical
	// method ?
	virtual bool NeedModelTime() const;
	// Do we need to evaluate the states at his timeIndex? 
	virtual bool NeedStatesEval(int timeIndex) const;
	
	virtual int ModelFixTimeStep( int fixTimeStep ) const;

	/// Volatilities and correlation time steps for PDEs
	virtual ARM_GP_VectorPtr VolatilitiesAndCorrelationTimesSteps() const;

	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
	/// non implemented method as it makes no sense for SFRM right now!

	virtual void NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	void PrepareForMarkovDriftCalculation(int timeIndex, 
										 double var_At_timeIndex,
										 int itsFwdMinIndex, 
										 int NbFwds, 
										 int nbIntegrationSteps, 
										 ARM_GP_Vector*& VarCoeff,
										 ARM_GP_Vector*& SpotVarCoeff,
                                         ARM_GP_Matrix*& FwdVolCoeff,
										 double*& x, 
										 double*& w) const;

	/// non implemented method as it makes no sense for SFRM right now!
    virtual double VarianceToTime(double var,double minTime,double maxTime) const;

	virtual void ComputeUnderlyingVarCovar( double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   ARM_GP_Vector& VarCovar,
										   ARM_GP_Vector& swapfwd) const;

	virtual void ComputeSwaptionVarCovar( double	fromTime,
										   double	toTime,
										   const ARM_VanillaSwaptionArgSFRM& swaptionArg_Long,
										   const ARM_VanillaSwaptionArgSFRM& swaptionArg_Short,
										   ARM_GP_Vector& VarCovar,
										   ARM_GP_Vector& swapfwd) const;

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

	virtual ARM_VectorPtr  VanillaSpreadOptionLet(const string& curveName,
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

    virtual ARM_VectorPtr VanillaSumOption(
		const string& curveName,
		double evalTime,
		int capFloor,
		const ARM_GP_Vector& coeffs,
		const ARM_GP_Vector& fwdResetTimes,
		const ARM_GP_Vector& fwdStartTimes,
		const ARM_GP_Vector& fwdEndTimes,
		double payTime,
		const ARM_GP_Vector& fwdPeriods,
		const ARM_GP_Vector& strikesPerState,
		double volatilityRatio,
		double* sumFwd,
		double* sumVol,
		const ARM_PricingStatesPtr& states) const;

	/// function for the generic tree (2nd generation)
	virtual void VolatilitiesAndCorrelations( const ARM_GP_Vector& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	/// function for the generic tree (2nd generation)
    virtual void EulerLocalDrifts(const ARM_GP_Vector& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;
	
	/// method to convert the lognormProbability
    /// to shift curve or vice versa
    void ConvertToShiftorBetaParam(const ARM_Portfolio& ConvPort);    
	virtual void PostInit();

	/// Sizeof ModelStates
	virtual size_t ModelStatesSize() const { return itsFwdValues->size(); }

    /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter);
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalib (ARM_ModelFitter& modelFitter);
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);
    /// method to advise the break point times
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );

#ifdef __GP_SFRM_STOCH_TERMS
    inline double GetStochasticTerm(size_t stateIdx, size_t modelStateIdx) const { return (*itsStochasticTerms)(stateIdx,modelStateIdx); };
    inline void SetStochasticTerm(size_t stateIdx, size_t modelStateIdx, double value) { (*itsStochasticTerms)(stateIdx,modelStateIdx) = value; };
#endif

	void SetFixStartDate(
		const ARM_Date& fixStartDate);
	void SetFixEndDate(
		const ARM_Date& fixStartDate);

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LSFRM";}	
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
